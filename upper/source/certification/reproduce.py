"""Explicitly choose recorded arithmetic checks or fresh full regeneration."""
import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

if not __debug__:
    raise RuntimeError('Run without Python -O')
import flint

HERE = Path(__file__).resolve().parent
ap = argparse.ArgumentParser()
mode = ap.add_mutually_exclusive_group(required=True)
mode.add_argument('--check-recorded', action='store_true')
mode.add_argument('--regenerate', action='store_true')
ap.add_argument('--degree', type=int, choices=[501, 681], default=681)
args = ap.parse_args()
D = args.degree
if args.check_recorded:
    # This recomputes endpoint arithmetic and verifies component links.
    # It deliberately does not claim to regenerate integrals or native arrays.
    stages = [('certify_upper.py', '--degree', str(D), '--check-only')]
else:
    stages = [
        (f'upper_head_{D}/certified_b_matrix.py',),
        (f'upper_head_{D}/contract_certified_head.py',),
        (f'upper_head_{D}/verify_plane_budget.py',),
        ('weighted_staircase/certify_head.py', '--degree', str(D), '--prec', '384'),
        (f'degree{D}/certify_linear_tails.py',),
        ('weighted_staircase/certify_nonadditive.py', '--degree', str(D)),
        ('upper_tail/certify_quintic_continuation.py',),
        ('upper_tail/certify_factorized_energies.py',),
        ('upper_tail/certify_heavy_tail.py' if D == 501 else
         'upper_tail/certify_heavy_tail_fine.py',),
        ('compose_head.py', '--degree', str(D)),
        ('certify_rare_additive_tails.py', '--degree', str(D)),
        ('certify_upper.py', '--degree', str(D)),
    ]

run_name = time.strftime('%Y%m%dT%H%M%SZ', time.gmtime())
logs = HERE / f'degree{D}' / 'reproduction_logs' / run_name
logs.mkdir(parents=True, exist_ok=True)
records = []
start = time.time()
for number, stage in enumerate(stages, 1):
    path = HERE / stage[0]
    if not path.is_file():
        raise FileNotFoundError(path)
    log = logs / f'{number:02d}_{path.stem}.log'
    tic = time.time()
    print(f'Running {number}/{len(stages)}: {stage[0]}', flush=True)
    with log.open('wb') as out:
        run = subprocess.run([sys.executable, str(path), *stage[1:]],
                             cwd=HERE, stdout=out, stderr=subprocess.STDOUT)
    records.append(dict(stage=list(stage), exit_code=run.returncode,
                        seconds=time.time() - tic, log=str(log)))
    if run.returncode:
        print(f'FAILED: inspect {log}', flush=True)
        raise SystemExit(run.returncode)

assert 'FULL_REAL_GROTHENDIECK_UPPER_BOUND_ACCEPTED=1' in log.read_text()
status = ('FRESH_FULL_UPPER_REGENERATION_ACCEPTED' if args.regenerate else
          'RECORDED_UPPER_COMPONENT_ARITHMETIC_ACCEPTED')
result = dict(status=status, degree=D, seconds=time.time() - start,
              python=sys.executable, python_flint_version=flint.__version__,
              stages=records,
              lower_endpoint='Not certified. This run covers only the upper endpoint.')
(logs / 'RESULT.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
print(status + '=1')
