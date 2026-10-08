#!/usr/bin/env python3
"""Regenerate the complete degree-681 real Grothendieck upper certificate."""
if not __debug__:
    raise RuntimeError('Run ordinary Python without -O: proof checks must remain enabled.')

import argparse
import ctypes
import hashlib
import importlib.metadata
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import time
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BOUND = Fraction(1779893, 1000000)
RECIPE_SHA = '7e042c037c80f416e799f345badfaacaf496f04c6bf4748bfc7ae73ff06176fe'


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_manifest():
    entries = {}
    for line in (ROOT/'MANIFEST.sha256').read_text().splitlines():
        digest, name = line.split(None, 1)
        relative = Path(name)
        require(not relative.is_absolute() and '..' not in relative.parts,
                'Invalid manifest path')
        require(name not in entries, 'Duplicate manifest path')
        path = ROOT/relative
        require(path.is_file() and not path.is_symlink(), 'Missing/invalid file: '+name)
        require(sha(path) == digest, 'File hash mismatch: '+name)
        entries[name] = digest
    expected_sources = {name for name in entries if name.startswith('source/')}
    actual_sources = {str(p.relative_to(ROOT)) for p in (ROOT/'source').rglob('*') if p.is_file()}
    require(actual_sources == expected_sources and len(expected_sources) == 18,
            'The source tree must contain exactly the 18 frozen source/input files')
    require('verify.py' in entries and 'requirements.txt' in entries, 'Incomplete manifest')
    return entries


def check_work_sources(work, manifest):
    for name, digest in manifest.items():
        if name.startswith('source/'):
            copied = work/Path(name).relative_to('source')
            require(copied.is_file() and sha(copied) == digest,
                    'Copied source/input changed: '+name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'verification_runs',
                        help='Parent for a new, uniquely named run directory')
    args = parser.parse_args()
    manifest = check_manifest()
    import numpy as np
    import flint
    from flint import arb, ctx
    require(np.__version__ == '2.3.5' and flint.__version__ == '0.9.0',
            'Install the exact public dependencies in requirements.txt')
    require(np.finfo(np.float64).eps == 2**-52 and np.dtype(np.float64).itemsize == 8,
            'IEEE-754 binary64 is required')
    try:
        get_round = ctypes.CDLL(None).fegetround
    except AttributeError:
        get_round = None
    if get_round is not None:
        require(get_round() == 0, 'Floating-point rounding must be round-to-nearest')

    output = args.output.resolve()
    require(not output.is_relative_to((ROOT/'source').resolve()),
            'Output must be outside the frozen source directory')
    output.mkdir(parents=True, exist_ok=True)
    run = Path(tempfile.mkdtemp(prefix='run-', dir=output))
    work = run/'work'
    for name in manifest:
        if name.startswith('source/'):
            dest = work/Path(name).relative_to('source')
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT/name, dest)
    check_work_sources(work, manifest)
    env = os.environ.copy()
    env.pop('PYTHONPATH', None)
    env.pop('PYTHONHOME', None)
    env.update(PYTHONDONTWRITEBYTECODE='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    started = time.time()
    print('Fresh output: '+str(run), flush=True)
    (run/'RUNNING.json').write_text(json.dumps({'status':'RUNNING', 'source_manifest':manifest}, indent=2)+'\n')

    def execute(relative, arguments, log_name):
        command = [sys.executable, str(work/relative), *arguments]
        with (run/log_name).open('w') as log:
            process = subprocess.Popen(command, cwd=work, env=env,
                                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            for line in process.stdout:
                log.write(line)
                log.flush()
                print(line, end='', flush=True)
            require(process.wait() == 0, 'Regeneration failed; inspect '+str(run/log_name))

    execute('certification/weighted_staircase/freeze_recipe.py', [], 'recipe.log')
    require(sha(work/'certification/weighted_staircase/explicit_staircase_mixture.json') == RECIPE_SHA,
            'Regenerated rounding recipe differs from the certified recipe')
    execute('certification/reproduce.py', ['--regenerate', '--degree', '681'], 'regeneration.log')
    results = list((work/'certification/degree681/reproduction_logs').glob('*/RESULT.json'))
    require(len(results) == 1, 'Expected one fresh pipeline result')
    result = json.loads(results[0].read_text())
    require(result['status'] == 'FRESH_FULL_UPPER_REGENERATION_ACCEPTED', 'Full regeneration did not accept')
    require(len(result['stages']) == 12 and all(s['exit_code'] == 0 for s in result['stages']),
            'Incomplete regeneration stages')
    cert_path = work/'certification/degree681/UPPER_BOUND_CERTIFICATE.json'
    certificate = json.loads(cert_path.read_text())
    require(certificate['status'] == 'FULL_REAL_GROTHENDIECK_UPPER_BOUND_ACCEPTED', 'No accepted upper certificate')
    require(Fraction(certificate['exact_strict_upper']) == BOUND, 'Incorrect endpoint')
    require(certificate['atoms_covered'] == 32 and certificate['uncovered_atoms'] == [], 'Incomplete rounding coverage')
    ctx.prec = 256
    margin = arb(certificate['complete_rounding_margin_lower_ball'])
    require(margin > arb(BOUND.denominator)/BOUND.numerator, 'Strict final margin failed')
    check_work_sources(work, manifest)
    check_manifest()
    summary = {
        'status':'FRESH_COMPLETE_UPPER_CERTIFICATE_VERIFIED',
        'theorem':'K_G^R < 1779893/1000000 = 1.779893',
        'exact_strict_upper':str(BOUND),
        'complete_margin_lower_ball':certificate['complete_rounding_margin_lower_ball'],
        'reciprocal_margin_upper_ball':certificate['reciprocal_margin_upper_ball'],
        'atoms_covered':32, 'uncovered_atoms':[],
        'fresh_pipeline_stages':12, 'rounding_recipe_regenerated':True,
        'seconds':time.time()-started,
        'python':sys.version, 'platform':platform.platform(),
        'packages':{name:importlib.metadata.version(name) for name in ['numpy','python-flint']},
        'numerical_model':certificate['numerical_model'],
        'source_manifest':manifest,
        'complete_certificate':'work/certification/degree681/UPPER_BOUND_CERTIFICATE.json',
        'pipeline_result':str(results[0].relative_to(run)),
    }
    # Commit acceptance only after all fresh stages and final checks succeed.
    temporary = run/'VERIFIED.json.tmp'
    temporary.write_text(json.dumps(summary, indent=2)+'\n')
    temporary.replace(run/'VERIFIED.json')
    (run/'RUNNING.json').unlink()
    print('FRESH_COMPLETE_UPPER_CERTIFICATE_VERIFIED=1', flush=True)
    print('CERTIFIED_KG_UPPER=1779893/1000000', flush=True)
    print('RESULT='+str(run/'VERIFIED.json'), flush=True)


if __name__ == '__main__':
    main()
