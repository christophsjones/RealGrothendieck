"""Entire rare-additive tail via positive composition and native Parseval mass."""
import argparse
import hashlib
import json
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
plane = np.load(plane_path)
meta_path = plane_dir / 'certified_plane.json'
meta = json.loads(meta_path.read_text())
assert meta['candidate_sha256'] == recipe['original_recipe_sha256']
assert meta['status'] == 'DIRECTED_FINITE_NATIVE_COEFFICIENT_PLANE'
assert meta['D'] == D
assert meta['plane_sha256'] == hashlib.sha256(plane_path.read_bytes()).hexdigest()
def A(x):
    x = Fraction(x)
    return arb(x.numerator) / x.denominator
def truncate(p):
    return arb_poly(p.coeffs()[:D + 1])
weights = [Fraction(a['weight']) for a in recipe['active']]
total = sum(weights)
rows = [arb_poly([arb(float(plane[a, b])) for b in range(D + 1 - a)])
        for a in range(D + 1)]
error = arb(meta['plane_l1_error']) + arb(meta['whole_kernel_cutoff_error'])
tail = arb(0)
records = []
for index, (atom, w) in enumerate(zip(recipe['active'], weights)):
    if atom['type'] != 'additive' or index in (0, 15, 16):
        continue
    w /= total
    cv = atom.get('powers', {'1': atom.get('r1', '0'), '3': atom.get('r3', '0')})
    cv = {int(k): Fraction(v) for k, v in cv.items()}
    budget = max(Fraction(1), sum(abs(v) for v in cv.values()))
    rc = [arb(0)] * (max(cv) + 1)
    for d, v in cv.items():
        rc[d] = A(abs(v) / budget)
    R = arb_poly(rc)
    power = arb_poly([1])
    captured = arb(0)
    for a in range(D + 1):
        captured += sum(truncate(rows[a] * power).coeffs(), arb(0))
        power = truncate(power * R)
    # Exact native coefficients for F_cutoff are nonnegative and sum<=1.
    # The positive covariance majorant has total mass<=1 and dominates
    # every absolute coefficient of the actual composed cutoff kernel.
    # The original atom differs by the separately certified cutoff error.
    bound = 1 - captured + error
    assert bound > 0 and bound < 1
    contribution = A(w) * bound
    tail += contribution
    records.append(dict(index=index, exact_weight=str(w),
                        positive_head_representative_mass=str(captured),
                        entire_atom_tail_upper_ball=str(bound),
                        weighted_tail_upper_ball=str(contribution)))
assert len(records) == 13
out = dict(status='ALL_THIRTEEN_RARE_ADDITIVE_TAILS_CERTIFIED',
           source_sha256=hashlib.sha256(raw).hexdigest(), degree=D,
           verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           native_plane_sha256=hashlib.sha256(plane_path.read_bytes()).hexdigest(),
           native_plane_metadata_sha256=hashlib.sha256(meta_path.read_bytes()).hexdigest(),
           whole_weighted_rare_additive_tail_upper_ball=str(tail), atoms=records,
           scope='All additive atoms except indices0,15,16. The three heavy additive tails and universal lower-instance norm remain unproved.')
(output / 'rare_additive_tail_certificate.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps({k: v for k, v in out.items() if k != 'atoms'}, indent=2))
print('RARE_ADDITIVE_FULL_TAIL_ACCEPTED=1')
