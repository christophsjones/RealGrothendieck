"""Exact final comparison. Accepts only fresh, complete replay outputs."""
if not __debug__:raise RuntimeError('Run without -O')
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
from replay import load_json,require,sha
V=F(1000,1773)-F(1,10**9)

def assemble(root,output):
 root=Path(root);output=Path(output)
 source=root/'baseline/lower_push_20260906/multiple_caps/best_candidate/frozen_rows.json'
 qsource=source.with_name('frozen_q.json');gen=root/'baseline/certify_above_1773/arithmetic/certify_arithmetic.py'
 a=load_json(output/'arithmetic.json')
 require(a['status']=='ALL_NON_FIBRE_CONDITIONS_CERTIFIED_UNIVERSAL_FIBRES_STILL_REQUIRED','Wrong arithmetic scope')
 require(a['source_sha256']==sha(source) and a['q_source_sha256']==sha(qsource),'Arithmetic source mismatch')
 require(a['generator_sha256']==sha(gen),'Arithmetic generator changed')
 require(F(a['target_game_norm'])==V,'Non-fibre target mismatch')
 require(F(load_json(root/'inputs/target_allowances.json')['norm_target'])==V,'Target file mismatch')
 require(a['precision_bits']==384 and a['epsilon']=='1/10000000','Unexpected arithmetic parameters')
 require(a['matrix_PD_checks']==1820 and a['intervals_cover_exact_unit_interval'] is True,'Incomplete non-fibre checks')
 require(F(a['minimum_LDL_pivot_lower'])>0,'Nonpositive matrix margin')
 require(len(a['rows'])==35 and [r['index'] for r in a['rows']]==list(range(35)),'Wrong arithmetic rows')
 q=json.loads(qsource.read_text(),parse_float=str)['coefficients_exact_decimal'];require(len(q)==11 and F(str(q[0]))==1,'Multiplier normalization')
 for j in range(11):require(F(a['q'][j])==F(str(q[j])),'Multiplier coefficient mismatch')
 # The original JSON decimals are also checked through the unchanged generator's source hash.
 rows=[];stats={};last=F(0);maximum=None
 for i in range(35):
  r=load_json(output/'fibres'/f'row{i}.json');ar=a['rows'][i]
  require(r['status']=='UNIVERSAL_FIBRE_REPLAY_PASSED' and r['analytic_regenerated'] is True,'A universal fibre was not regenerated')
  require(r['row']==i and r['source_sha256']==sha(source),'Fibre source mismatch')
  require(r['certificate_sha256']==sha(root/'certificates'/f'row{i}.json'),'Fibre evidence changed after replay')
  require(r['input_sha256']==sha(root/'inputs'/f'row{i}.bin'),'Fibre input changed after replay')
  lo,hi=map(F,ar['interval']);require(lo==last and lo<hi<=1,'Scalar coverage gap');last=hi
  require(ar['psd_checks_in_row']==52 and F(ar['minimum_LDL_pivot_lower'])>0,'Incomplete row matrix checks')
  require(len(ar['caps'])==2,'Wrong number of angular caps')
  for cap in ar['caps']:
   require(cap['all_finite_omitted_coefficients_positive'] is True,'Unproved angular tail')
   require(F(cap['finite_omitted_min_lower'])>0 and F(cap['tail_cone_margin_lower'])>0,'Angular positivity failure')
  cf=F(r['certified_fibre_upper']);other=F(ar['non_fibre_cost_upper']);total=cf+other
  require(cf<=F(r['allowed']) and total<=V,'Final row exceeds norm target')
  maximum=total if maximum is None else max(maximum,total)
  rows.append(dict(index=i,interval=[str(lo),str(hi)],grid_nodes_per_channel=r['grid_nodes_per_channel'],certified_fibre_upper=str(cf),allowed=str(F(r['allowed'])),fibre_strict_slack=str(F(r['allowed'])-cf),non_fibre_upper=str(other),combined_norm_upper=str(total),strict_slack_below_v=str(V-total),statistics=r['statistics']))
  for key,val in r['statistics'].items():stats[key]=stats.get(key,0)+val
 require(last==1 and maximum is not None,'Incomplete scalar domain')
 require(V<F(1000,1773),'Requested strict rational comparison failed')
 controls=load_json(output/'negative_controls.json');require(controls['status']=='ALL_NEGATIVE_CONTROLS_REJECTED','Negative controls not complete')
 limit=load_json(output/'FROZEN_ROW_LIMIT.json')
 require(limit['status']=='FROZEN_ROW_FAMILY_LIMIT_CERTIFIED' and limit['source_sha256']==sha(source),'Frozen-family limit not verified')
 result=dict(status='KG_1773_THEOREM_CERTIFIED',KG_theorem_accepted=True,all_35_universal_fibres_regenerated=True,all_non_fibre_arithmetic_regenerated=True,root_forests_required=False,analytic_reduction='positive quadratic Bernstein interpolation',frozen_row_limit_certified=True,target_KG='1773/1000',norm_bound=str(V),certified_KG_lower=str(1/V),strict_KG_improvement=str(1/V-F(1773,1000)),maximum_assembled_norm_upper=str(maximum),source_sha256=sha(source),q_source_sha256=sha(qsource),matrix_PD_checks=1820,original_fibre_M_PD_checks=35,precision_bits=384,analytic_backend=a['python_flint_version'],finite_arithmetic='exact checked int64 / int128 and rational numbers',negative_controls_rejected=controls['rejected_count'],statistics=stats,rows=rows)
 (output/'VERIFICATION_REPORT.json').write_text(json.dumps(result,indent=2)+'\n')
 print('CERTIFIED_FIBRES=35');print('ALL_NON_FIBRE_CONDITIONS_CERTIFIED=1');print('KG_1773_THEOREM_ACCEPTED=1');print('CERTIFIED_KG_LOWER='+str(1/V))
 return result
if __name__=='__main__':
 import argparse
 ap=argparse.ArgumentParser();ap.add_argument('--root',default=str(Path(__file__).resolve().parents[1]));ap.add_argument('--output',required=True);args=ap.parse_args();assemble(args.root,args.output)
