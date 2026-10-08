#!/usr/bin/env python3
"""Build and regenerate the complete real Grothendieck certificate.

This command never accepts stored success flags in lieu of the arithmetic.
Run with ordinary Python, not python -O. It requires NumPy, Cython,
setuptools, a C++17 compiler, GMP development files and the MPFR library.
"""
if not __debug__:raise RuntimeError('Certificate verification requires Python without -O')
import os,sys,argparse,subprocess,shutil,json,hashlib,time,platform,concurrent.futures
import ctypes,ctypes.util,importlib.metadata
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def run(command,log,env,cwd=None,expect_success=True):
 with Path(log).open('w') as f:rc=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,env=env,cwd=cwd).returncode
 if expect_success and rc:raise RuntimeError(f'Command failed ({rc}); inspect {log}')
 return rc

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=4);ap.add_argument('--output',default=str(ROOT/'verification_run'));args=ap.parse_args()
 if not 1<=args.workers<=32:raise ValueError('workers must be between 1 and 32')
 out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=True)
 # The report is committed only after every component has regenerated.
 for name in ['VERIFICATION_REPORT.json','STATUS.json']:
  if (out/name).exists():(out/name).unlink()
 logs=out/'logs';logs.mkdir(exist_ok=True);build=out/'build';build.mkdir(exist_ok=True)
 env=os.environ.copy();env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
 compiler=shutil.which(os.environ.get('CXX','g++')) or shutil.which('clang++')
 if compiler is None:raise RuntimeError('A C++17 compiler is required')
 library=build/'exact_graph.so'
 run([compiler,'-O3','-std=c++17','-fPIC','-shared',str(ROOT/'tools/exact_graph.cpp'),'-o',str(library)],logs/'integer_build.log',env)
 mpfr=build/'mpfr';mpfr.mkdir(exist_ok=True)
 for name in ['balls.hpp','flint.pyx','setup.py']:shutil.copy2(ROOT/'backend'/name,mpfr/name)
 run([sys.executable,'setup.py','build_ext','--inplace','--force'],logs/'mpfr_build.log',env,cwd=mpfr)
 env['PYTHONPATH']=str(mpfr)+os.pathsep+str(ROOT/'tools')
 env['KG_EXACT_GRAPH_LIBRARY']=str(library)
 run([sys.executable,str(ROOT/'backend/test_backend.py')],logs/'mpfr_tests.log',env)
 run([sys.executable,str(ROOT/'tools/test_graph.py')],logs/'integer_tests.log',env)
 run([sys.executable,str(ROOT/'tools/test_analytic.py')],logs/'analytic_tests.log',env)
 print('ARITHMETIC_SELF_TESTS_PASSED',flush=True)
 fibres=out/'fibres';fibres.mkdir(exist_ok=True)
 def row_task(row):
  p=ROOT/'certificates'/f'row{row}.json'
  if not p.exists():raise RuntimeError(f'Missing universal fibre certificate {row}')
  run([sys.executable,str(ROOT/'tools/replay.py'),'--root',str(ROOT),'--row',str(row),'--out',str(fibres)],logs/f'row{row}.log',env)
  print(f'UNIVERSAL_FIBRE_{row}_REGENERATED',flush=True)
 with concurrent.futures.ThreadPoolExecutor(args.workers) as pool:list(pool.map(row_task,range(35)))
 # Preserve the unchanged original non-fibre checker, and give it a scratch
 # workspace with exactly the source layout it expects.
 ar=out/'arithmetic_workspace';src=ar/'lower_push_20260906/multiple_caps/best_candidate';src.mkdir(parents=True,exist_ok=True)
 checker=ar/'certify_above_1773/arithmetic';checker.mkdir(parents=True,exist_ok=True)
 for name in ['frozen_rows.json','frozen_q.json']:shutil.copy2(ROOT/'baseline/lower_push_20260906/multiple_caps/best_candidate'/name,src/name)
 shutil.copy2(ROOT/'baseline/certify_above_1773/arithmetic/certify_arithmetic.py',checker/'certify_arithmetic.py')
 run([sys.executable,str(checker/'certify_arithmetic.py'),'--precision','384'],logs/'non_fibre.log',env)
 shutil.copy2(checker/'arithmetic_certificate.json',out/'arithmetic.json')
 # This is an actual failed matrix computation, not merely a malformed JSON.
 rc=run([sys.executable,str(checker/'certify_arithmetic.py'),'--precision','384','--epsilon','0.000000000001'],logs/'insufficient_repair.log',env,expect_success=False)
 if rc==0:raise RuntimeError('Insufficient-repair negative control unexpectedly passed')
 text=(logs/'insufficient_repair.log').read_text()
 if not any(x in text for x in ['AssertionError','ArithmeticError']) or 'PD' not in text or 'row0 M2 transverse3' not in text:raise RuntimeError('Negative arithmetic failed for an unexpected reason')
 run([sys.executable,str(ROOT/'tools/test_evidence.py'),'--root',str(ROOT),'--output',str(out/'negative_controls.json'),'--arithmetic-negative-log',str(logs/'insufficient_repair.log')],logs/'negative_controls.log',env)
 run([sys.executable,str(ROOT/'tools/frozen_row_limit.py'),'--root',str(ROOT),'--output',str(out/'FROZEN_ROW_LIMIT.json')],logs/'frozen_row_limit.log',env)
 run([sys.executable,str(ROOT/'tools/assemble_final.py'),'--root',str(ROOT),'--output',str(out)],logs/'final_assembly.log',env)
 report=json.loads((out/'VERIFICATION_REPORT.json').read_text())
 if report['status']!='KG_1773_THEOREM_CERTIFIED' or report['KG_theorem_accepted'] is not True:raise RuntimeError('Final assembly did not accept')
 environment=dict(python=sys.version,platform=platform.platform(),compiler=subprocess.check_output([compiler,'--version'],text=True).splitlines()[0],workers=args.workers,analytic_backend=report['analytic_backend'],precision_bits=384,integer_source_sha256=hashlib.sha256((ROOT/'tools/exact_graph.cpp').read_bytes()).hexdigest(),mpfr_source_sha256={name:hashlib.sha256((ROOT/'backend'/name).read_bytes()).hexdigest() for name in ['balls.hpp','flint.pyx']})
 environment['python_packages']={name:importlib.metadata.version(name) for name in ['numpy','Cython','setuptools']}
 mpfrlib=ctypes.CDLL(ctypes.util.find_library('mpfr'));mpfrlib.mpfr_get_version.restype=ctypes.c_char_p
 environment['MPFR_version']=mpfrlib.mpfr_get_version().decode()
 gmplib=ctypes.CDLL(ctypes.util.find_library('gmp'))
 environment['GMP_version']=ctypes.c_char_p.in_dll(gmplib,'__gmp_version').value.decode()
 environment['verifier_sources_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((ROOT/'tools').glob('*')) if p.suffix in ['.py','.cpp']}
 environment['driver_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
 (out/'ENVIRONMENT.json').write_text(json.dumps(environment,indent=2)+'\n')
 (out/'STATUS.json').write_text(json.dumps({k:report[k] for k in ['status','KG_theorem_accepted','all_35_universal_fibres_regenerated','all_non_fibre_arithmetic_regenerated','norm_bound','certified_KG_lower','matrix_PD_checks','negative_controls_rejected']},indent=2)+'\n')
 print((logs/'final_assembly.log').read_text(),end='')
 print('REPORT='+str(out/'VERIFICATION_REPORT.json'))
if __name__=='__main__':main()
