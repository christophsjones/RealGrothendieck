#!/usr/bin/env python3
"""Exact linkage of the manuscript data to the browsable certificate sources.

This module checks finite input identities, not the analytic numerical claims.
It requires only the Python standard library. JSON decimals are never converted
to binary floating point.
"""
from __future__ import annotations
import hashlib
import json
import runpy
from fractions import Fraction as Q
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
ROWS_PATH = 'baseline/lower_push_20260906/multiple_caps/best_candidate/frozen_rows.json'
Q_PATH = 'baseline/lower_push_20260906/multiple_caps/best_candidate/frozen_q.json'
MIX_PATH = 'source/upper/refined_explicit_mixture.json'
STAIR_PATH = 'source/competitors/staircase_3_5.json'

def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)

def load(raw: bytes | str) -> Any:
    return json.loads(raw, parse_float=str)

def frac(value: Any) -> str:
    return str(Q(str(value)))

def vector(values: list[Any]) -> list[str]:
    return [frac(x) for x in values]

def check_manifest() -> dict[str, str]:
    """Check every executable/input file and reject unexpected package files."""
    entries = {}
    for line in (ROOT/'MANIFEST.sha256').read_text().splitlines():
        digest, name = line.split(None, 1)
        rel = Path(name)
        require(not rel.is_absolute() and '..' not in rel.parts, 'Unsafe manifest path')
        require(name not in entries, 'Duplicate manifest entry: '+name)
        require(len(digest)==64 and all(c in '0123456789abcdef' for c in digest),
                'Invalid SHA-256 digest')
        path = ROOT/rel
        require(path.is_file() and not any(p.is_symlink() for p in [path,*path.parents]
                                         if p != ROOT.parent), 'Missing/invalid file: '+name)
        require(hashlib.sha256(path.read_bytes()).hexdigest()==digest,
                'File hash mismatch: '+name)
        entries[name] = digest
    actual = {'verify.py','data_link.py','requirements.txt'}
    for folder in ['data','lower','upper']:
        for path in (ROOT/folder).rglob('*'):
            require(not path.is_symlink(), 'Symlink in frozen package: '+str(path))
            if path.is_file():
                actual.add(str(path.relative_to(ROOT)))
    require(set(entries)==actual, 'Source/input inventory differs from MANIFEST.sha256')
    return entries

def regenerate_definitions() -> dict[str, Any]:
    """Return exactly the mathematical input fields used by the manuscript."""
    check_manifest()
    read = lambda p: load((ROOT/'lower'/p).read_bytes())
    source_rows = read(ROWS_PATH)['rows']
    multiplier = read(Q_PATH)['coefficients_exact_decimal']
    allowances = read('inputs/target_allowances.json')['rows']
    rows = []
    for i, row in enumerate(source_rows):
        out = {
            'index': i, 'interval': vector(row['interval']),
            'ell_h': vector(row['lh']), 'ell_k': vector(row['lk']),
            't': frac(row['t']), 'c_h': frac(row['ch']), 'c_k': frac(row['ck']),
            'C': frac(allowances[i]['allowed_fibre_upper']),
            'angular': [{'N': x['N'], 'weights': vector(x['weights'])} for x in row['caps']],
            'counts': read(f'inputs/row{i}_analytic.json')['counts'],
        }
        for key, matrix in zip(('H0','H1','K0','K1'),row['matrices']):
            out[key] = [vector(v) for v in matrix]
        rows.append(out)
    lower = {
        'degree': 21, 'epsilon': '1/10000000',
        'scalar_bound': '999999998227/1773000000000',
        'q': {str(2*i+1): frac(x) for i,x in enumerate(multiplier)},
        'rows': rows,
    }
    def read_upper(path: str) -> Any:
        return load((ROOT/'upper'/path).read_bytes())
    mix = read_upper(MIX_PATH)
    stair_source = read_upper(STAIR_PATH)
    components = []
    for i, item in enumerate(mix['active']):
        typ = {'linear':'hyperplane','weighted_3_5':'staircase'}.get(item['type'],item['type'])
        row = {'index':i, 'type':typ, 'raw_weight':frac(item['weight'])}
        if typ != 'staircase':
            powers = item.get('powers',{'1':item.get('r1','0'),'3':item.get('r3','0')})
            row['raw_covariance'] = {str(d):frac(v) for d,v in powers.items() if Q(str(v))}
        components.append(row)
    threshold = mix['threshold_definition']
    upper = {
        'degree':681,
        'weight_sum':frac(sum(Q(x['raw_weight']) for x in components)),
        'threshold':dict(zip(map(str,threshold['r_normalized_probabilists_Hermite_degrees']),
                             vector(threshold['r_coefficients_exact_decimal']))),
        'components':components,
    }
    staircase = {
        'denominator':stair_source['denominator'],
        'end':stair_source['end'],
        'rows':[{'breakpoints':vector(row['roots']),'signs':row['signs']}
                for row in stair_source['rows']],
    }
    return {'lower_rows.json':lower,'upper_rounding.json':upper,'staircase.json':staircase}

def check_definitions(definitions: dict[str, Any]) -> dict[str, Any]:
    """Check exact dimensions, normalization, coverage, and rational endpoints."""
    lower = definitions['lower_rows.json']
    rows = lower['rows']
    require(len(rows)==35, 'Expected 35 lower rows')
    previous = Q(0)
    for i,row in enumerate(rows):
        require(row['index']==i, 'Wrong lower row index')
        a,b = map(Q,row['interval'])
        require(a==previous and a<b<=1, 'Intervals do not partition [0,1]')
        previous=b
        for key in ('H0','H1','K0','K1'):
            matrix = [[Q(x) for x in v] for v in row[key]]
            require(len(matrix)==11 and all(len(v)==11 for v in matrix), 'Wrong matrix dimensions')
            require(all(matrix[j][k]==matrix[k][j] for j in range(11) for k in range(11)),
                    'Non-symmetric rational matrix')
        require(len(row['ell_h'])==len(row['ell_k'])==11, 'Wrong linear dimensions')
        require([x['N'] for x in row['angular']]==[7,9], 'Wrong angular orders')
        require(all(len(x['weights'])==(x['N']-1)//2 and Q(x['weights'][0])>0
                    for x in row['angular']), 'Invalid angular inputs')
        counts=row['counts']
        require(len(counts)==36 and all(type(n) is int and n>0 for n in counts), 'Invalid panel counts')
    require(previous==1, 'Incomplete lower interval cover')
    require(Q(lower['q']['1'])==1 and all(abs(Q(q))<=1 for q in lower['q'].values()),'Invalid multiplier')
    upper=definitions['upper_rounding.json']
    components=upper['components']
    require(len(components)==32, 'Expected 32 components')
    require(sum(Q(c['raw_weight']) for c in components)==Q(upper['weight_sum'])==
            Q(25000000000000000247431,25000000000000000000000), 'Weight sum mismatch')
    for i,c in enumerate(components):
        require(c['index']==i and Q(c['raw_weight'])>0, 'Invalid component/weight')
        if i==14:
            require(c['type']=='staircase', 'Wrong staircase index')
            continue
        require(c['type']==('additive' if i<17 else 'hyperplane'), 'Wrong rounding component')
        poly={int(k):Q(v) for k,v in c['raw_covariance'].items()}
        require(all(k>0 and k%2 for k in poly), 'Covariance is not an odd polynomial')
        norm=sum(abs(v) for v in poly.values()); scale=max(Q(1),norm)
        require(sum(abs(v/scale) for v in poly.values())<=1, 'Invalid covariance normalization')
    stair=definitions['staircase.json']
    require(stair['denominator']==256 and Q(str(stair['end']))==8 and len(stair['rows'])==2048,
            'Invalid staircase strip cover')
    for row in stair['rows']:
        roots=list(map(Q,row['breakpoints']))
        require(all(a<b for a,b in zip(roots,roots[1:])), 'Unordered staircase breakpoints')
        require(len(row['signs'])==len(roots)+1 and all(x in (-1,1) for x in row['signs']),
                'Invalid staircase signs')
    v=Q(lower['scalar_bound'])
    require(v==Q(1000,1773)-Q(1,10**9) and 1/v>Q(1773,1000), 'Lower endpoint mismatch')
    margin=Q('0.5618773799453352')-Q('0.000003878150168295')-Q('0.000034097111321496')-Q('0.000007875733380418')
    require(margin==Q('0.561831528950464991') and margin>Q(1000000,1779893), 'Upper rational comparison failed')
    require(1/v>Q('1.773') and Q('1.779893')<Q('1.781'), 'Title interval failed')
    return {'lower_rows':35,'upper_components':32,'staircase_strips':2048,
            'exact_lower_endpoint':str(1/v),'upper_comparison_margin':str(margin),
            'analytic_numerical_certificates_replayed':False}

def check_data() -> dict[str, Any]:
    definitions=regenerate_definitions()
    for filename,expected in definitions.items():
        path=ROOT/'data'/filename
        require(path.is_file(), 'Missing manuscript data '+filename)
        require(load(path.read_text())==expected, 'Data-to-certificate mismatch '+filename)
    result = check_definitions(definitions)
    # Use the exact checker distributed with lower-bound-shorter-cert.zip.
    tail = runpy.run_path(str(ROOT/'lower'/'check-tail.py'))['check'](definitions['lower_rows.json'])
    require(tail['accepted'] is True and tail['checks']==70, 'Incomplete rotation tail check')
    result['rotation_tail_checks'] = tail['checks']
    result['rotation_tail_minimum_gap'] = tail['minimum_gap']
    return result

if __name__=='__main__':
    print(json.dumps(check_data(),indent=2))
    print('DATA_LINKS_CHECKED=1')
    print('NO_ANALYTIC_NUMERICAL_REPLAY=1')
