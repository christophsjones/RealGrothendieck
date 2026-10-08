"""Assemble the full upper bound, checking all component links.

This checks the final directed arithmetic and frozen component links.
Fresh analytic generation uses the separate pipeline commands in README.md.
"""
import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from flint import arb, ctx

if not __debug__:
    raise RuntimeError('Run without Python -O')
ctx.prec = 256
HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument('--degree', type=int, choices=[501, 681], default=501)
parser.add_argument('--check-only', action='store_true', help='Validate without replacing the recorded certificate.')
args = parser.parse_args()
D = args.degree
OUT = HERE / f'degree{D}'
def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    return json.loads(p.read_text())
def A(v):
    v = Fraction(v)
    return arb(v.numerator) / v.denominator

source = HERE / 'weighted_staircase/explicit_staircase_mixture.json'
recipe = read(source)
head_path = OUT / 'head_certificate.json'
rare_path = OUT / 'rare_additive_tail_certificate.json'
heavy_name = 'heavy_tail_certificate.json' if D == 501 else 'heavy_tail_fine_certificate.json'
heavy_generator = 'certify_heavy_tail.py' if D == 501 else 'certify_heavy_tail_fine.py'
heavy_path = HERE / 'upper_tail' / heavy_name
nonadd_path = HERE / 'weighted_staircase' / f'nonadditive_certificate_d{D}.json'
head, rare, heavy, nonadd = map(read, [head_path, rare_path, heavy_path, nonadd_path])
original = HERE.parent / 'upper/refined_explicit_mixture.json'
assert recipe['original_recipe_sha256'] == heavy['source_sha256'] == sha(original)
assert head['source_sha256'] == rare['source_sha256'] == nonadd['source_sha256'] == sha(source)
assert head['degree'] == rare['degree'] == nonadd['degree'] == D
assert head['status'] == 'FINITE_ROUNDING_HEAD_CERTIFIED_ADDITIVE_TAILS_MISSING'
assert rare['status'] == 'ALL_THIRTEEN_RARE_ADDITIVE_TAILS_CERTIFIED'
assert heavy['status'] == 'RIGOROUS_HEAVY_ADDITIVE_FULL_TAIL_CERTIFIED'
assert nonadd['status'] == 'NONADDITIVE_HEAD_AND_FULL_TAIL_COMPONENT_CERTIFIED'
assert sha(HERE / f'upper_head_{D}/certified_plane.npy') == head['native_plane_sha256'] == rare['native_plane_sha256']
assert sha(HERE / f'upper_head_{D}/certified_plane.json') == head['native_plane_metadata_sha256'] == rare['native_plane_metadata_sha256']
assert sha(nonadd_path) == head['nonadditive_certificate_sha256']
assert sha(HERE / 'upper_tail/factorized_energies_certificate.json') == heavy['energy_certificate_sha256']
assert sha(HERE / 'upper_tail/early_capped_profile_certificate.json') == heavy['cap_certificate_sha256']
assert sha(HERE / 'upper_head/arc/certify_complex_arc.py') == heavy['complex_arc_verifier_sha256']
for p, generator in [(head_path, HERE / 'compose_head.py'),
                     (rare_path, HERE / 'certify_rare_additive_tails.py'),
                     (heavy_path, HERE / 'upper_tail' / heavy_generator),
                     (nonadd_path, HERE / 'weighted_staircase/certify_nonadditive.py')]:
    expected = read(p)['verifier_sha256']
    if D == 501 and expected != sha(generator):
        # Keep the completed501 checkpoint bound to its original sources.
        generator = HERE / 'degree501/checkpoint_generators' / generator.name
    assert expected == sha(generator)

additive_indices = {i for i, a in enumerate(recipe['active']) if a['type'] == 'additive'}
rare_indices = {a['index'] for a in rare['atoms']}
heavy_indices = set(heavy['heavy_indices'])
assert rare_indices.isdisjoint(heavy_indices)
assert rare_indices | heavy_indices == additive_indices
assert len(recipe['active']) - len(additive_indices) == nonadd['atom_count'] == 16
assert len(recipe['active']) == 32

# All quantities below are outward interval enclosures of the respective
# error budgets. The margin functional c1 - sum_{m>1}|cm| is 1-Lipschitz.
finite_margin = arb(head['certified_head_margin_lower_ball'])
tn = arb(nonadd['full_tail_upper_ball'])
tr = arb(rare['whole_weighted_rare_additive_tail_upper_ball'])
th = (arb(heavy['weighted_D4_L2_norm_bound']) / (14 * arb(D) ** 7).sqrt() +
      arb(heavy['weighted_cap_kernel_l1_perturbation']))
complete_margin = finite_margin - tn - tr - th
upper = Fraction(44499, 25000) if D == 501 else Fraction(1779893, 1000000)
strict_margin = Fraction('0.5618104458749009') if D == 501 else 1 / upper
assert complete_margin > A(strict_margin)
assert complete_margin > 1 / A(upper)

out = dict(status='FULL_REAL_GROTHENDIECK_UPPER_BOUND_ACCEPTED', degree=D,
           theorem=f'K_G^R < {upper} = {float(upper):.6f}',
           exact_strict_upper=str(upper),
           complete_rounding_margin_lower_ball=str(complete_margin),
           reciprocal_margin_upper_ball=str(1 / complete_margin),
           strict_margin_lower=str(strict_margin),
           head_margin_lower_ball=str(finite_margin),
           nonadditive_full_tail_upper_ball=str(tn),
           rare_additive_full_tail_upper_ball=str(tr),
           heavy_additive_full_tail_upper_ball=str(th),
           atoms_covered=32, uncovered_atoms=[],
           source_sha256=sha(source), verifier_sha256=sha(Path(__file__)),
           component_sha256={str(p.relative_to(HERE)): sha(p) for p in
                             [head_path, rare_path, heavy_path, nonadd_path]},
           proof='All full tails and the entire finite head are enclosed. Apply the odd-Dirichlet Neumann completion to H=c1*t+E: a complete absolute margin B>0 realizes exact linear rounding with rate at least B, hence K_G^R<=1/B.',
           numerical_model='Directed Arb arithmetic plus explicit IEEE754 binary64 error bounds for the native coefficient contraction.',
           lower_endpoint_status='NOT_CERTIFIED: no dimension-uniform Boolean norm upper bound for the repaired degree21 game has been obtained.',
           scope='Upper bound only. This does not certify the proposed lower endpoint or a0.001 interval.')
if not args.check_only:
    (OUT / 'UPPER_BOUND_CERTIFICATE.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
print('FULL_REAL_GROTHENDIECK_UPPER_BOUND_ACCEPTED=1')
