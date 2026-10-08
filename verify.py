#!/usr/bin/env python3
"""Verify manuscript inputs and regenerate the requested numerical certificates.

No archive is executed: the browsable, hash-checked package files are copied to
fresh work directories. Only a successful complete run writes ACCEPTED.json.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
import platform
import runpy
import shutil
import subprocess
import sys
import tempfile
import time
from fractions import Fraction
from pathlib import Path
sys.dont_write_bytecode = True
from data_link import ROOT, check_data, check_manifest, require


def copy_package(kind: str, destination: Path, manifest: dict[str,str]) -> Path:
    destination.mkdir(parents=True, exist_ok=False)
    prefix=kind+'/'
    for name,digest in manifest.items():
        if name.startswith(prefix):
            target=destination/name[len(prefix):]
            target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(ROOT/name,target)
            require(hashlib.sha256(target.read_bytes()).hexdigest()==digest,
                    'Copied file differs: '+name)
    return destination


def execute(command: list[str], log_path: Path, working: Path) -> None:
    env=os.environ.copy()
    env.pop('PYTHONPATH',None)
    env.pop('PYTHONHOME',None)
    env.update(PYTHONDONTWRITEBYTECODE='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    with log_path.open('w') as log:
        process=subprocess.Popen(command,cwd=working,env=env,stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,text=True)
        for line in process.stdout:
            print(line,end='',flush=True)
            log.write(line)
            log.flush()
        require(process.wait()==0,'Verification failed; inspect '+str(log_path))


def check_lower(output: Path) -> None:
    report=json.loads((output/'VERIFICATION_REPORT.json').read_text())
    require(report['KG_theorem_accepted'] is True, 'Lower verifier did not accept')
    rows=json.loads((ROOT/'data'/'lower_rows.json').read_text())['rows']
    v=Fraction(1000,1773)-Fraction(1,10**9)
    require(len(report['rows'])==35, 'Incomplete lower report')
    for row, definition in zip(report['rows'], rows):
        require(row['index']==definition['index'], 'Lower row mismatch')
        require(Fraction(row['allowed'])==Fraction(definition['C']), 'Wrong one-dimensional allowance')
        require(Fraction(row['certified_fibre_upper'])<=Fraction(definition['C']), 'One-dimensional bound failed')
        # This also verifies the slightly stronger endpoint condition as printed.
        require(Fraction(definition['C'])+Fraction(row['non_fibre_upper'])<=v,
                'Printed lower endpoint condition failed')

def check_upper(output: Path) -> None:
    from flint import arb, ctx
    ctx.prec=256
    accepted=list(output.glob('run-*/VERIFIED.json'))
    require(len(accepted)==1, 'Expected exactly one fresh upper acceptance')
    certificate=accepted[0].parent/'work/certification/degree681/UPPER_BOUND_CERTIFICATE.json'
    data=json.loads(certificate.read_text())
    def rational(text: str):
        q=Fraction(text)
        return arb(q.numerator)/q.denominator
    tests=[('head_margin_lower_ball','0.5618773799453352','>'),
           ('nonadditive_full_tail_upper_ball','0.000003878150168295','<'),
           ('rare_additive_full_tail_upper_ball','0.000034097111321496','<'),
           ('heavy_additive_full_tail_upper_ball','0.000007875733380418','<')]
    for key, bound, relation in tests:
        value=arb(data[key]); endpoint=rational(bound)
        require(value>endpoint if relation=='>' else value<endpoint,
                'Printed numerical enclosure failed: '+key)
    require(arb(data['complete_rounding_margin_lower_ball'])>rational('1000000/1779893'),
            'Upper strict comparison failed')

def check_upper_paper_inequality() -> dict[str, object]:
    """The explicit scalar inequality displayed in the manuscript."""
    from flint import arb, ctx
    ctx.prec=256
    value = (64 * (-arb(512)).exp() * arb(32)**681
             / (345 * arb.pi() * (343 * arb.fac_ui(681)).sqrt()))
    require(value < arb(3) / 10**18, 'Manuscript Fourier-tail inequality failed')
    return {'accepted':True, 'precision_bits':256, 'upper_bound':'3/1000000000000000000',
            'expression':'64*exp(-512)*32^681/(345*pi*sqrt(343*681!))',
            'value_ball':str(value)}

def main() -> None:
    require(__debug__, 'Run ordinary Python without -O; assertions must be enabled')
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--lower',action='store_true')
    group.add_argument('--upper',action='store_true')
    group.add_argument('--all',action='store_true')
    group.add_argument('--data-only',action='store_true')
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--output',type=Path,default=ROOT/'runs')
    args=parser.parse_args()
    require(1<=args.workers<=32,'workers must be between 1 and 32')
    manifest=check_manifest()
    manifest_sha=hashlib.sha256((ROOT/'MANIFEST.sha256').read_bytes()).hexdigest()
    structural=check_data()
    print('DATA_LINKS_CHECKED=1',flush=True)
    print('EXACT_ROTATION_TAIL_CHECKS=70',flush=True)
    if args.data_only:
        print(json.dumps(structural,indent=2))
        print('NO_ANALYTIC_NUMERICAL_REPLAY=1')
        return
    chosen=['lower'] if args.lower else ['upper'] if args.upper else ['lower','upper']
    parent=args.output.resolve()
    for name in ['data','lower','upper']:
        require(not parent.is_relative_to((ROOT/name).resolve()), 'Output must be outside source/input directories')
    parent.mkdir(parents=True,exist_ok=True)
    run=Path(tempfile.mkdtemp(prefix='complete-',dir=parent))
    print('RUN_DIRECTORY='+str(run),flush=True)
    started=time.time()
    running={'status':'RUNNING','certificates_requested':chosen,'source_manifest_sha256':manifest_sha}
    (run/'RUNNING.json').write_text(json.dumps(running,indent=2)+'\n')
    (run/'DATA_CHECKS.json').write_text(json.dumps(structural,indent=2)+'\n')
    tail=runpy.run_path(str(ROOT/'lower'/'check-tail.py'))['check'](
        json.loads((ROOT/'data/lower_rows.json').read_text()))
    tail['input_sha256']=hashlib.sha256((ROOT/'data/lower_rows.json').read_bytes()).hexdigest()
    (run/'TAIL_CHECK.json').write_text(json.dumps(tail,indent=2)+'\n')
    for kind in chosen:
        package=copy_package(kind,run/(kind+'_source'),manifest)
        output=run/(kind+'_results')
        cmd=[sys.executable,str(package/'verify.py'),'--output',str(output)]
        if kind=='lower': cmd+=['--workers',str(args.workers)]
        execute(cmd,run/(kind+'.log'),package)
        (check_lower if kind=='lower' else check_upper)(output)
        if kind=='upper':
            scalar=check_upper_paper_inequality()
            (run/'UPPER_PAPER_CHECK.json').write_text(json.dumps(scalar,indent=2)+'\n')
    check_data()
    require(check_manifest()==manifest, 'Source manifest changed during verification')
    require(hashlib.sha256((ROOT/'MANIFEST.sha256').read_bytes()).hexdigest()==manifest_sha,
            'Manifest changed during verification')
    result={'status':'REQUESTED_NUMERICAL_CERTIFICATES_ACCEPTED',
            'certificates_replayed':chosen,
            'lower_bound':'1773000000000/999999998227' if 'lower' in chosen else None,
            'strict_upper_bound':'1779893/1000000' if 'upper' in chosen else None,
            'both_certificates_and_endpoint_checks_passed':len(chosen)==2,
            'exact_rotation_tail_checks':70,
            'upper_paper_scalar_checked':'upper' in chosen,
            'source_manifest_sha256':manifest_sha,
            'source_manifest':manifest,
            'seconds':time.time()-started,
            'python':sys.version,'platform':platform.platform()}
    target=run/'ACCEPTED.json'
    temporary=run/'ACCEPTED.json.tmp'
    temporary.write_text(json.dumps(result,indent=2)+'\n')
    temporary.replace(target)
    (run/'RUNNING.json').unlink()
    print('NUMERICAL_CERTIFICATES_ACCEPTED='+','.join(chosen))
    print('ACCEPTANCE_FILE='+str(target))


if __name__=='__main__':
    main()
