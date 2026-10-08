"""Certified limitation of tightening only the 35 frozen universal constants.

This is not an upper bound on the Grothendieck constant or on what another
proof could establish for the multiplier. At x=381/400 the existing rows,
even with best possible universal fibre constants and zero repair cost,
cannot yield a norm bound below 5639/10000.
"""
if not __debug__:raise RuntimeError('Run without -O')
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,argparse
from analytic import A,hermites,read_row,endpoint_fraction
from flint import arb

def circle(N):
 return [tuple(sum(s[j]*(s[j+d] if j+d<N else -s[j+d-N]) for j in range(N)) for d in range(1,(N+1)//2)) for s in ([1]+[1-2*((code>>(j-1))&1) for j in range(1,N)] for code in range(1<<(N-1)))]

def run(root,output):
 root=Path(root);src=root/'baseline/lower_push_20260906/multiple_caps/best_candidate/frozen_rows.json'
 rows=json.loads(src.read_text(),parse_float=F)['rows'];assert len(rows)==35
 x=F(381,400);floor=F(5639,10000);nu=(2/arb.pi()).sqrt();h=hermites(A(0),20)
 s=[A(0)]+[nu*h[j-1]/A(j).sqrt() for j in range(1,22)];m=s+s;res=[]
 for i,r in enumerate(rows):
  M,ell=read_row(src,i)
  witness=sum((A(ell[j])*m[j] for j in range(44)),A(0))
  witness+=sum((A(M[j][k])*m[j]*m[k] for j in range(44) for k in range(44)),A(0))
  gamma=A(0);cap=F(0)
  for q in r['caps']:
   N=q['N'];w=list(map(F,q['weights']));cap+=max(sum((v*z for v,z in zip(w,cut)),F(0))/N for cut in circle(N))
   gamma+=sum((A(v)*(arb.pi()*d/N).cos() for d,v in enumerate(w,1)),A(0))
  value=witness+A(cap+F(r['ch'])+F(r['ck']))+(1-gamma)*nu*nu*A(x*x)-A(F(r['t'])*x)*nu
  assert value>A(floor),(i,str(value))
  res.append(dict(row=i,sign_fibre_value_lower=str(endpoint_fraction(witness,False)),scalar_bound_floor=str(endpoint_fraction(value,False))))
 result=dict(status='FROZEN_ROW_FAMILY_LIMIT_CERTIFIED',source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),precision_bits=384,scalar_coordinate=str(x),witness='h(s)=sign(s), k(s)=0',repair_budget_included='0',all_35_row_bounds_strictly_above=str(floor),best_possible_KG_conclusion_with_only_these_row_constants_strictly_below=str(1/floor),not_an_upper_bound_on_KG=True,rows=res)
 Path(output).write_text(json.dumps(result,indent=2)+'\n');print('FROZEN_ROW_FAMILY_LIMIT_CERTIFIED');print('NORM_FLOOR='+str(floor));print('ROW_CONSTANT_ONLY_KG_CEILING='+str(1/floor));return result
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--root',default=str(Path(__file__).resolve().parents[1]));ap.add_argument('--output',required=True);a=ap.parse_args();run(a.root,a.output)
