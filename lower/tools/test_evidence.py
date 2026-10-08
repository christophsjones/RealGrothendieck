"""Fail-closed tests against malformed evidence and real arithmetic failures."""
if not __debug__:raise RuntimeError('Run without -O')
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
from fractions import Fraction as F
from copy import deepcopy
import argparse,json,tempfile,shutil
import numpy as np
from replay import load_json,replay_boxes,replay_binary,regenerated_caps,field_coordinate,regenerate_analytic,require
from graph import build,roof,cap,validate


def main(root,output,arithmetic_negative_log):
 root=Path(root);results=[]
 def reject(name,fn):
  try:fn()
  except (ValueError,ArithmeticError,AssertionError) as e:
   results.append(dict(test=name,rejected=True,reason=str(e)[:400]));return
  raise RuntimeError('Malformed/invalid evidence was accepted: '+name)

 meta=load_json(root/'inputs/row3_analytic.json');B,a,c,qe,mean=build(root/'inputs/row3.bin')
 allow=F(load_json(root/'inputs/target_allowances.json')['rows'][3]['allowed_fibre_upper'])
 error=F(meta['interpolation_and_tail_error'])+F(meta['fixed_point_input_error'])+qe
 zz=(allow-error)*(1<<44);target=zz.numerator//zz.denominator
 original=load_json(root/'certificates/row3.json')['proof']
 def check(record):return replay_boxes(B,a,c,mean,target,record)
 def trial(name,mutate):
  r=deepcopy(original);mutate(r);reject(name,lambda:check(r))
 trial('unfinished_status',lambda r:r.update(status='CERTIFICATION_IN_PROGRESS'))
 trial('nonempty_pending_set',lambda r:r.update(pending=[dict(id=999,lo=[-65536,-65536],hi=[65536,65536])]))
 trial('changed_integer_target',lambda r:r.update(target_integer=r['target_integer']+1))
 trial('unjustified_symmetry_prefix',lambda r:r.update(prefix=[r['prefix'][0]+1,1]))
 trial('changed_coordinate_entry',lambda r:r['coordinates'][0].__setitem__(0,r['coordinates'][0][0]+1))
 trial('changed_coordinate_normalizer',lambda r:r['coordinate_normalizers'].__setitem__(0,r['coordinate_normalizers'][0]+1))
 trial('duplicate_root_node',lambda r:r['nodes'].append(deepcopy(next(n for n in r['nodes'] if n['id']==0))))
 trial('missing_proof_node',lambda r:r['nodes'].pop())
 trial('rectangle_root_endpoint_gap',lambda r:next(n for n in r['nodes'] if n['id']==0)['lo'].__setitem__(0,-65535))
 split=next((n for n in original['nodes'] if n['kind']=='split'),None)
 if split:
  sid=split['id'];child=split['children'][0]
  trial('duplicate_rectangle_children',lambda r:next(n for n in r['nodes'] if n['id']==sid)['children'].__setitem__(1,child))
  trial('closed_child_coverage_gap',lambda r:next(n for n in r['nodes'] if n['id']==child)['hi'].__setitem__(split['coordinate'],split['at']-1))
 leaf=next(n for n in original['nodes'] if n['kind']=='flow');lid=leaf['id']
 def leafcert(r):return next(n for n in r['nodes'] if n['id']==lid)['certificate']
 trial('negative_cap_price',lambda r:leafcert(r)['prices'].__setitem__(0,-1))
 trial('understated_flow_bound',lambda r:leafcert(r).__setitem__('bound',leafcert(r)['bound']-1))
 trial('understated_cap_rounding_error',lambda r:leafcert(r).__setitem__('rounding_error_units',leafcert(r)['rounding_error_units']-1))
 trial('unsupported_exclusion_of_feasible_rectangle',lambda r:next(n for n in r['nodes'] if n['id']==lid).update(kind='empty',facet=[0,0,0]))
 trial('unresolved_point_mislabeled_as_completion',lambda r:next(n for n in r['nodes'] if n['id']==lid).update(kind='unresolved'))
 trial('understated_global_integer_bound',lambda r:r.update(max_upper_integer=r['max_upper_integer']-1))

 # Test a real binary subtree using exactly its parent rectangle's capped graph.
 data=load_json(root/'certificates/row2.json')['proof'];B2,a2,c2,qe2,v2=build(root/'inputs/row2.bin')
 node=next(n for n in data['nodes'] if n['kind']=='branch' and any(k['kind']=='split' and k.get('persistent') for k in n['certificate']['tree']['nodes']))
 vs=[field_coordinate(a2),v2];bb,aa,cc,err=regenerated_caps(B2,a2,c2,vs,(node['lo'],node['hi']),node['certificate']['prices'],data['prefix'][0])
 tree=node['certificate']['tree'];tt=tree['target_integer']
 def btrial(name,mutate):
  r=deepcopy(tree);mutate(r);reject(name,lambda:replay_binary(bb,aa,cc,tt,r))
 btrial('binary_missing_child',lambda r:(r['nodes'].pop(),r.update(node_count=len(r['nodes']))))
 btrial('binary_active_frontier',lambda r:r.update(active_count=1))
 btrial('binary_unjustified_reflection_prefix',lambda r:r.update(reflection_prefix=[[0,1]]))
 btrial('binary_wrong_persistent_sign',lambda r:next(n for n in r['nodes'] if n['kind']=='split' and n.get('persistent'))['persistent'][0].__setitem__(1,-next(n for n in r['nodes'] if n['kind']=='split' and n.get('persistent'))['persistent'][0][1]))
 btrial('binary_missing_persistent_assignment',lambda r:next(n for n in r['nodes'] if n['kind']=='split' and n.get('persistent'))['persistent'].pop())
 btrial('binary_understated_flow',lambda r:r['nodes'][0].__setitem__('flow',r['nodes'][0]['flow']-1))

 # Directed analytic checks are actually rerun, not simulated by a JSON flag.
 with tempfile.TemporaryDirectory(prefix='kg_negative_') as temp:
  temp=Path(temp);(temp/'baseline').symlink_to(root/'baseline',target_is_directory=True);(temp/'inputs').mkdir()
  shutil.copy2(root/'inputs/row3.bin',temp/'inputs/row3.bin')
  for field in ['tail_error','interpolation_and_tail_error','fixed_point_input_error']:
   bad=deepcopy(meta);bad[field]='0';(temp/'inputs/row3_analytic.json').write_text(json.dumps(bad))
   reject('omitted_'+field,lambda:regenerate_analytic(temp,3,temp/'replay'))
  bad=deepcopy(meta);bad['source_sha256']='0'*64;(temp/'inputs/row3_analytic.json').write_text(json.dumps(bad))
  reject('wrong_frozen_source_hash',lambda:regenerate_analytic(temp,3,temp/'replay'))

 small=np.zeros((2,2),dtype=np.int64);unary=np.zeros(2,dtype=np.int64)
 bad=small.copy();bad[0,1]=1
 reject('asymmetric_integer_graph',lambda:validate(bad,unary,0))
 reject('integer_overflow_budget',lambda:validate(small,unary,1<<61))
 reject('zero_scalar_coordinate',lambda:cap(small,unary,0,np.zeros(2,dtype=np.int64),1,-65536,65536))
 reject('reversed_cap_interval',lambda:cap(small,unary,0,np.array([1,-1],dtype=np.int64),1,1,0))

 log=Path(arithmetic_negative_log).read_text()
 require('PD check failed: row0 M2 transverse3' in log,'Missing actual insufficient-repair failure')
 results.append(dict(test='actual_non_fibre_repair_1e_minus_12',rejected=True,reason='Directed LDL rejects row0 M2 transverse3'))
 report=dict(status='ALL_NEGATIVE_CONTROLS_REJECTED',rejected_count=len(results),tests=results)
 Path(output).write_text(json.dumps(report,indent=2)+'\n')
 print('NEGATIVE_CONTROLS_REJECTED='+str(len(results)))
 return report
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--root',default=str(Path(__file__).resolve().parents[1]));ap.add_argument('--output',required=True);ap.add_argument('--arithmetic-negative-log',required=True);a=ap.parse_args();main(a.root,a.output,a.arithmetic_negative_log)
