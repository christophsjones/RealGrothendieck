"""Enclose every non-additive head coefficient and the full joint tail."""
import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from flint import arb, arb_series, ctx

if not __debug__:
    raise RuntimeError('Run without Python -O')
ctx.prec = 384
parser = argparse.ArgumentParser()
parser.add_argument('--degree', type=int, choices=[301, 501, 681], default=301)
D = parser.parse_args().degree
ctx.cap = D + 2
HERE = Path(__file__).resolve().parent
BASE = HERE.parent.parent
source = HERE / 'explicit_staircase_mixture.json'
raw = source.read_bytes()
data = json.loads(raw)
def A(x):
    x = Fraction(x)
    return arb(x.numerator) / x.denominator
weights = [Fraction(a['weight']) for a in data['active']]
total = sum(weights)
head = {m: arb(0) for m in range(1, D + 1, 2)}
staircase = json.loads((HERE / f'certificate_d{D}_p384.json').read_text())
staircase_hash = hashlib.sha256((BASE / 'competitors/staircase_3_5.json').read_bytes()).hexdigest()
assert staircase['source_sha256'] == staircase_hash
tail_dir = 'upper_tail' if D == 301 else f'degree{D}'
linear_tail = json.loads((HERE.parent / tail_dir / 'linear_tail_certificate.json').read_text())
assert linear_tail['cutoff'] == D
original_hash = hashlib.sha256((BASE / 'upper/refined_explicit_mixture.json').read_bytes()).hexdigest()
assert linear_tail['source_sha256'] == original_hash == data['original_recipe_sha256']
count = 0
for atom, w in zip(data['active'], weights):
    if atom['type'] == 'additive':
        continue
    w /= total
    count += 1
    if atom['type'] == 'linear':
        cv = {int(k): Fraction(v) for k, v in atom['powers'].items()}
        budget = max(Fraction(1), sum(abs(v) for v in cv.values()))
        p = [arb(0) for _ in range(D + 1)]
        for m, v in cv.items():
            p[m] = A(v / budget)
        series = arb_series(p, D + 1).asin() * (2 / arb.pi())
        for m in head:
            head[m] += A(w) * series[m]
    else:
        assert atom['type'] == 'weighted_3_5'
        assert atom['realization'] == 'exact_frozen_staircase'
        assert atom['source_sha256'] == staircase_hash
        assert w == Fraction(staircase['exact_normalized_weight'])
        for m in head:
            head[m] += A(w) * arb(staircase['chaos_coefficients'][str(m)])
tail = arb(linear_tail['weighted_linear_tail_ball']) + arb(staircase['weighted_full_tail_upper_ball'])
assert tail < A('0.000005834544')
out = dict(status='NONADDITIVE_HEAD_AND_FULL_TAIL_COMPONENT_CERTIFIED',
           source_sha256=hashlib.sha256(raw).hexdigest(), degree=D,
           precision_bits=ctx.prec, atom_count=count,
           verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           head_coefficients={str(m): str(c) for m, c in head.items()},
           full_tail_upper_ball=str(tail),
           full_tail_strict_upper='729318/125000000000',
           scope='All fifteen linear atoms and the staircase replacement atom. No additive coefficients or additive tails are included. No K_G endpoint is certified.')
filename = 'nonadditive_certificate.json' if D == 301 else f'nonadditive_certificate_d{D}.json'
(HERE / filename).write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps({k: v for k, v in out.items() if k != 'head_coefficients'}, indent=2))
print('NONADDITIVE_HEAD_AND_FULL_TAIL_ACCEPTED=1')
