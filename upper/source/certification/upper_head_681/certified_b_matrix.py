"""Directed inputs and an explicit IEEE-754 error budget for the B matrix."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import json,time,math
import numpy as np
from flint import arb,ctx
from full_head_plan import HERE,rp,D
ctx.prec=512
plan=json.loads((HERE/'full_head_quadrature_plan.json').read_text());tic=time.time()
inner=[];rs=[];vrows=[];v_error=arb(0);vnorm=arb(0)
order=plan['inner_gauss_order'];roots=[arb.legendre_p_root(order,k,weight=True) for k in range(order)]
sqrt2pi=(2*arb.pi()).sqrt()
for aa,bb in plan['inner_cells']:
 a,b=arb(aa),arb(bb);center=(a+b)/2;half=(b-a)/2
 for z,weight in roots:
  x=center+half*z;weight*=half;rs.append(rp(x));h0=arb(1);h1=x
  gaussian=2*weight*(-x*x/2).exp()/sqrt2pi
  row=[];esq=arb(0);normsq=arb(0)
  for j in range(D+1):
   if j==0:h=h0
   elif j==1:h=h1
   else:
    h=(x*h1-arb(j-1).sqrt()*h0)/arb(j).sqrt();h0,h1=h1,h
   value=gaussian*h;f=float(value.mid());error=abs(value-arb(f)).upper()
   esq+=error*error;normsq+=abs(value).upper()**2;row.append(f)
  vrows.append(row);v_error+=esq.sqrt();vnorm+=normsq.sqrt()
print('inner data prepared',len(rs),'seconds',time.time()-tic,flush=True)
v=np.array(vrows);outer=[];weights=[]
order=plan['outer_gauss_order'];roots=[arb.legendre_p_root(order,k,weight=True) for k in range(order)]
for panel in range(32):
 for z,w in roots:outer.append(arb(panel)+arb('0.5')*(1+z));weights.append(w/2)
cosine=np.empty((len(outer),len(rs)));sine=np.empty_like(cosine);trig_error=arb(2)**-48
for i,s in enumerate(outer):
 for j,r in enumerate(rs):
  phase=s*r;c=phase.cos();q=phase.sin();cf=float(c.mid());sf=float(q.mid())
  if not (abs(c-arb(cf))<trig_error and abs(q-arb(sf))<trig_error):raise RuntimeError('trigonometric rounding budget failed')
  cosine[i,j]=cf;sine[i,j]=sf
 if i%64==0:print('outer nodes',i,'seconds',time.time()-tic,flush=True)
# optimize=False performs the indicated ordinary multiply-and-add contraction.
b=np.empty((len(outer),D+1));b[:,::2]=np.einsum('sn,nj->sj',cosine,v[:,::2],optimize=False);b[:,1::2]=np.einsum('sn,nj->sj',sine,v[:,1::2],optimize=False)
b*=(-1.)**(np.arange(D+1)//2)
unit=arb(2)**-53;gamm=(2*len(rs)*unit)/(1-2*len(rs)*unit)
input_error=v_error+trig_error*(vnorm+v_error)
dot_error=gamm*(1+trig_error)*(vnorm+v_error)
berror=arb(plan['inner_vector_quadrature_error'])+input_error+dot_error
u=np.empty((D+1,len(outer)));usums=[];uerror=arb(0)
for i,(s,w) in enumerate(zip(outer,weights)):
 value=w*(-s*s/2).exp()/s
 for a in range(D+1):
  if a:value*=s/arb(a).sqrt()
  f=float(value.mid());u[a,i]=f
  uerror+=abs(value-arb(f)).upper()
usuml1=[sum(abs(arb(float(x))) for x in row) for row in u]
usquare=sum(x*x for x in usuml1);usum=sum(usuml1)
assert usquare>1  # Makes the following conservative factor-error square valid.
positive_w_over_s=sum(w/s for w,s in zip(weights,outer))
a_error_inner=(2/arb.pi())*positive_w_over_s*3*berror*(1+berror)**2
a_error=arb(plan['outer_Hilbert_quadrature_error'])+arb(plan['W_projection_Fourier_tail_Hilbert'])+a_error_inner
plane_error_quadrature=2*a_error+a_error*a_error
factor=4/arb.pi()**2;ff=float(factor.mid());factor_error=abs(factor-arb(ff)).upper()
plane_error_u=factor*(2*usum*uerror+uerror*uerror)+factor_error*(usquare+uerror)**2
bound=1+berror
observed=float(np.max(np.linalg.norm(b,axis=1)))
if observed>float(bound.upper())+1e-12:raise RuntimeError('B norm check failed')
np.savez(HERE/'certified_b_inputs.npz',B=b,U=u,factor=ff)
result={'status':'CERTIFIED_B_MATRIX_INPUTS_FULL_CONTRACTION_PENDING','D':D,'nodes':len(outer),'inner_nodes_positive_half':len(rs),'candidate_sha256':plan['candidate_sha256'],'B_vector_error':berror.str(40),'B_norm_bound':bound.str(40),'observed_float_B_norm':observed,'V_rounding_vector_error':v_error.str(40),'V_sum_vector_norms':vnorm.str(40),'U_total_l1_error':uerror.str(40),'U_l1_squares_sum':usquare.str(40),'U_l1_sum':usum.str(40),'A_Hilbert_error':a_error.str(40),'plane_l1_quadrature_error':plane_error_quadrature.str(40),'plane_l1_U_input_error':plane_error_u.str(40),'whole_kernel_cutoff_error':plan['whole_kernel_cutoff_error'],'seconds':time.time()-tic,'arithmetic_assumption':'Binary64 round-to-nearest; each einsum output is the displayed N-term ordinary sum, with at most N multiplications and N additions; FMA reduces the same bound. All transcendental inputs evaluated by Arb.'}
(HERE/'certified_b_inputs.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
