"""Freeze an explicit alternative recipe without modifying the original."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = HERE.parent.parent
source = BASE / 'upper/refined_explicit_mixture.json'
raw = source.read_bytes()
original = json.loads(raw, parse_float=str)
staircase = BASE / 'competitors/staircase_3_5.json'
out = {k: original[k] for k in ('active', 'interpretation', 'threshold_definition')}
out['status'] = 'EXPLICIT_MODIFIED_CANDIDATE_COMPONENTS_ONLY_CERTIFIED'
out['original_recipe_sha256'] = hashlib.sha256(raw).hexdigest()
out['modification'] = ('Only the weighted_3_5 atom is replaced by the frozen '
    'staircase pair F(x,y), G(x,y)=F(-x,y), using coordinate correlations t^3 '
    'and t^5. All probability weights and other functions are unchanged. '
    'No old numerical endpoint or nonlinear aggregate is asserted for this recipe.')
atoms = [a for a in out['active'] if a['type'] == 'weighted_3_5']
assert len(atoms) == 1 and atoms[0]['sign'] == 1
atoms[0]['source'] = '../../competitors/staircase_3_5.json'
atoms[0]['source_sha256'] = hashlib.sha256(staircase.read_bytes()).hexdigest()
atoms[0]['realization'] = 'exact_frozen_staircase'
out['weighted_staircase_definition'] = {
    'source': '../../competitors/staircase_3_5.json',
    'source_sha256': hashlib.sha256(staircase.read_bytes()).hexdigest(),
    'gaussian_coordinate_correlations': {'x': 't^3', 'y': 't^5'},
    'G': 'F(-x,y)',
    'positive_y': 'On y in [k/256,(k+1)/256), use row k x breakpoints and signs.',
    'outside': 'On y>=8, F(x,y)=-sign(x).',
    'negative_y': 'F(x,y)=-F(-x,-y).',
    'boundary_values': 'Arbitrary; all boundaries have Gaussian measure zero.'
}
target = HERE / 'explicit_staircase_mixture.json'
target.write_text(json.dumps(out, indent=2) + '\n')
print('ALTERNATIVE_RECIPE_FROZEN', hashlib.sha256(target.read_bytes()).hexdigest())
