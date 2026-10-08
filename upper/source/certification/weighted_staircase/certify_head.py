"""Arb head and entire unresolved energy for a fixed Boolean staircase.

This certifies one rounding component, not a Grothendieck endpoint.
The frozen staircase is integrated exactly on every rectangle, including
its infinite outside strips. All decimals are read as rational numbers.
"""
import argparse
import hashlib
import json
import math
import time
from fractions import Fraction
from pathlib import Path
from flint import arb, ctx

if not __debug__:
    raise RuntimeError('Run without Python -O')
ap = argparse.ArgumentParser()
ap.add_argument('--prec', type=int, default=256)
ap.add_argument('--degree', type=int, default=301)
args = ap.parse_args()
assert args.degree > 0 and args.degree % 2 == 1
ctx.prec = args.prec
HERE = Path(__file__).resolve().parent
BASE = HERE.parent.parent
source = BASE / 'competitors/staircase_3_5.json'
raw = source.read_bytes()
assert hashlib.sha256(raw).hexdigest() == '2d9a68c3ea888090e3dd2633dd89d9559aa3849a3e4d2c488971f4934d319528'
data = json.loads(raw)
assert data['weights'] == [3, 5]
D = args.degree
imax, jmax = D // 3, D // 5
indices = [(i, j) for i in range(imax + 1) for j in range(jmax + 1)
           if 0 < 3 * i + 5 * j <= D and (i + j) % 2]
den, end, rows = data['denominator'], data['end'], data['rows']
assert type(den) is int and den > 0 and type(end) is int and end > 0
assert len(rows) == den * end

def A(x):
    x = Fraction(x)
    return arb(x.numerator) / x.denominator

sqrt2, sqrt2pi = arb(2).sqrt(), (2 * arb.pi()).sqrt()

def primitive(q, n):
    x = A(q)
    phi = (-x * x / 2).exp() / sqrt2pi
    out = [(x / sqrt2).erfc() / 2]
    if n:
        hs = [arb(1)]
        if n > 1:
            hs.append(x)
        for k in range(1, n - 1):
            hs.append(x * hs[-1] - k * hs[-2])
        out.extend(phi * h for h in hs)
    return out

tic = time.time()
u = [arb(0) for _ in indices]
previous = primitive(0, jmax)
for k, row in enumerate(rows):
    roots, signs = [Fraction(v) for v in row['roots']], row['signs']
    assert all(a < b for a, b in zip(roots, roots[1:]))
    assert len(signs) == len(roots) + 1
    assert all(type(s) is int and s in (-1, 1) for s in signs)
    # Tail CDF representation avoids cancellation near large positive roots.
    mx = [arb(signs[0])] + [arb(0) for _ in range(imax)]
    for r, left, right in zip(roots, signs, signs[1:]):
        jump = right - left
        p = primitive(r, imax)
        mx[0] += jump * p[0]
        for i in range(1, imax + 1):
            mx[i] += jump * p[i]
    current = primitive(Fraction(k + 1, den), jmax)
    iy = [a - b for a, b in zip(previous, current)]
    for h, (i, j) in enumerate(indices):
        u[h] += 2 * mx[i] * iy[j]
    previous = current
    if (k + 1) % 512 == 0:
        print(f'Integrated {k + 1}/{len(rows)} strips; {time.time() - tic:.1f}s', flush=True)

# y >= end: F(x,y) = -sign(x). Its odd extension gives the negative side.
zero = primitive(0, imax)
for h, (i, j) in enumerate(indices):
    if i % 2:
        u[h] += -4 * zero[i] * previous[j]

chaos = {m: arb(0) for m in range(1, D + 1, 2)}
energy = arb(0)
factorials = [math.factorial(i) for i in range(max(imax, jmax) + 1)]
for moment, (i, j) in zip(u, indices):
    atom_energy = moment * moment / (factorials[i] * factorials[j])
    energy += atom_energy
    chaos[3 * i + 5 * j] += (-1) ** i * atom_energy
head_mass = sum((abs(c) for c in chaos.values()), arb(0))
# Parseval gives total energy 1; reflection only changes Hermite signs.
# Sum of absolute coefficients in all omitted weighted degrees is at most
# the sum of their energies, which is exactly 1 - the included energy.
unresolved = 1 - energy
assert energy > 0 and energy < 1
assert unresolved > 0 and unresolved < 1

recipe_raw = (BASE / 'upper/refined_explicit_mixture.json').read_bytes()
recipe = json.loads(recipe_raw, parse_float=Fraction)
atoms = recipe['active']
total = sum((Fraction(str(a['weight'])) for a in atoms), Fraction(0))
weighted = [a for a in atoms if a['type'] == 'weighted_3_5']
assert len(weighted) == 1
weight = Fraction(str(weighted[0]['weight'])) / total
weighted_tail = A(weight) * unresolved
out = dict(status='FIXED_STAIRCASE_HEAD_AND_FULL_TAIL_CERTIFIED',
           precision_bits=args.prec, degree=D, source_sha256=hashlib.sha256(raw).hexdigest(),
           source=str(source), verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           original_recipe_sha256=hashlib.sha256(recipe_raw).hexdigest(),
           indices=len(indices), included_energy=str(energy), head_l1=str(head_mass),
           full_tail_upper_ball=str(unresolved), exact_normalized_weight=str(weight),
           weighted_full_tail_upper_ball=str(weighted_tail),
           chaos_coefficients={str(m): str(c) for m, c in chaos.items()},
           seconds=time.time() - tic,
           proof_scope='Only the separately specified staircase replacement component. No complete rounding endpoint or universal game norm is certified.')
target = HERE / f'certificate_d{D}_p{args.prec}.json'
target.write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps({k: v for k, v in out.items() if k != 'chaos_coefficients'}, indent=2))
print('FIXED_STAIRCASE_HEAD_AND_FULL_TAIL_ACCEPTED=1')
