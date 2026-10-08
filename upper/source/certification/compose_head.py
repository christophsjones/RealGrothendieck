"""Directed polynomial composition of the certified native coefficient plane.

The final bound here is a finite-head margin, NOT a rounding-rate bound:
all additive tails must still be supplied.
"""
import argparse
import hashlib
import json
import time
from fractions import Fraction
from pathlib import Path
import numpy as np
from flint import arb, arb_poly, ctx

if not __debug__:
    raise RuntimeError('Run without Python -O')
ctx.prec = 256
HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument('--degree', type=int, choices=[301, 501, 681], default=301)
D = parser.parse_args().degree
output = HERE if D == 301 else HERE / f'degree{D}'
output.mkdir(exist_ok=True)
plane_dir = HERE / ('upper_head' if D == 301 else f'upper_head_{D}')
source = HERE / 'weighted_staircase/explicit_staircase_mixture.json'
raw = source.read_bytes()
recipe = json.loads(raw)
plane_path = plane_dir / 'certified_plane.npy'
plane_raw = plane_path.read_bytes()
plane = np.load(plane_path)
meta_path = plane_dir / 'certified_plane.json'
meta = json.loads(meta_path.read_text())
nonadd_path = HERE / 'weighted_staircase' / ('nonadditive_certificate.json' if D == 301 else f'nonadditive_certificate_d{D}.json')
nonadd = json.loads(nonadd_path.read_text())
assert meta['status'] == 'DIRECTED_FINITE_NATIVE_COEFFICIENT_PLANE'
assert meta['candidate_sha256'] == recipe['original_recipe_sha256']
assert meta['plane_sha256'] == hashlib.sha256(plane_raw).hexdigest()
assert nonadd['source_sha256'] == hashlib.sha256(raw).hexdigest()
assert nonadd['degree'] == meta['D'] == D
assert plane.shape == (D + 1, D + 1)
assert np.isfinite(plane).all() and (plane >= 0).all()

def A(x):
    x = Fraction(x)
    return arb(x.numerator) / x.denominator

weights = [Fraction(a['weight']) for a in recipe['active']]
total = sum(weights)
additive_weight = sum(w for a, w in zip(recipe['active'], weights)
                      if a['type'] == 'additive') / total
def truncate(poly):
    return arb_poly(poly.coeffs()[:D + 1])

native_rows = []
for a in range(D + 1):
    # A binary64 datum has an exact terminating binary-rational value.
    # Arb receives that value at precision256, before all compositions.
    row = [(-1) ** b * arb(float(plane[a, b])) for b in range(D + 1 - a)]
    native_rows.append(arb_poly(row))

aggregate = arb_poly([])
tic = time.time()
for index, (atom, w) in enumerate(zip(recipe['active'], weights)):
    if atom['type'] != 'additive':
        continue
    cv = atom.get('powers', {'1': atom.get('r1', '0'), '3': atom.get('r3', '0')})
    cv = {int(k): Fraction(v) for k, v in cv.items()}
    budget = max(Fraction(1), sum(abs(v) for v in cv.values()))
    rc = [arb(0)] * (max(cv) + 1)
    for degree, value in cv.items():
        rc[degree] = A(value / budget)
    R = arb_poly(rc)
    power = arb_poly([1])
    atom_head = arb_poly([])
    for a in range(D + 1):
        atom_head += truncate(native_rows[a] * power)
        power = truncate(power * R)
    aggregate += atom_head * A(w / total)
    print(f'Composed additive atom {index}; {time.time() - tic:.1f}s', flush=True)

head = {}
for m in range(1, D + 1, 2):
    head[m] = aggregate[m] + arb(nonadd['head_coefficients'][str(m)])
head_mass = sum((abs(c) for m, c in head.items() if m > 1), arb(0))
representative_margin = head[1] - head_mass
global_error = A(additive_weight) * (arb(meta['plane_l1_error']) +
                                    arb(meta['whole_kernel_cutoff_error']))
# c1 - sum_{m>=3}|cm| is 1-Lipschitz in the coefficient l1 norm.
certified_head_margin = representative_margin - global_error
after_nonadd_tail = certified_head_margin - arb(nonadd['full_tail_upper_ball'])
out = dict(status='FINITE_ROUNDING_HEAD_CERTIFIED_ADDITIVE_TAILS_MISSING',
           source_sha256=hashlib.sha256(raw).hexdigest(), degree=D,
           verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           native_plane_sha256=hashlib.sha256(plane_raw).hexdigest(),
           native_plane_metadata_sha256=hashlib.sha256(meta_path.read_bytes()).hexdigest(),
           nonadditive_certificate_sha256=hashlib.sha256(nonadd_path.read_bytes()).hexdigest(),
           arithmetic_precision_bits=ctx.prec,
           exact_normalized_additive_weight=str(additive_weight),
           representative_linear_coefficient=str(head[1]),
           representative_nonlinear_head_l1=str(head_mass),
           representative_head_margin=str(representative_margin),
           whole_head_l1_error=str(global_error),
           certified_head_margin_lower_ball=str(certified_head_margin),
           after_nonadditive_tail_lower_ball=str(after_nonadd_tail),
           allowable_additive_tail_for_KG_below_1_78=str(after_nonadd_tail - A(Fraction(50, 89))),
           head_representative_coefficients={str(m): str(c) for m, c in head.items()},
           seconds=time.time() - tic,
           scope='The true full rate is at least after_nonadditive_tail_lower_ball minus the entire omitted additive tail. That tail is not supplied here. No K_G endpoint is certified.')
(output / 'head_certificate.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps({k: v for k, v in out.items() if k != 'head_representative_coefficients'}, indent=2))
print('FINITE_ROUNDING_HEAD_ACCEPTED=1')
