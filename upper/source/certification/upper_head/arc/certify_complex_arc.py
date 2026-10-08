"""Directed universal D^4 Mehler-kernel bound on closed unit-circle arcs.

The three auxiliary coordinate pairs are grouped into chi-square_3 variables.
Gamma orthogonal polynomials diagonalize the score-polynomial moment exactly.
"""
import argparse,json,math,hashlib,time
from pathlib import Path
from fractions import Fraction
from flint import arb,acb,ctx

HERE=Path(__file__).resolve().parent;ORDER=4;ZERO=(0,0,0,0)
def indices(n,k):
 if n==1:return [(j,) for j in range(k+1)]
 return [(j,)+a for j in range(k+1) for a in indices(n-1,k-j)]
KEYS=indices(4,4)
def pochhammer(a,k):
 out=Fraction(1)
 for j in range(k):out*=a+j
 return out
alphas=[Fraction(1,2),Fraction(1,2),Fraction(3,2),Fraction(3,2)]
TRANS={}
NORMS={}
for n in KEYS:
 TRANS[n]=[]
 for k in KEYS:
  if all(ki<=ni for ki,ni in zip(k,n)):
   v=Fraction(1)
   for ni,ki,a in zip(n,k,alphas):v*=2**(ni-ki)*math.comb(ni,ki)*pochhammer(a+ki,ni-ki)
   assert v.denominator==1;TRANS[n].append((k,v.numerator))
 v=Fraction(1)
 for ki,a in zip(n,alphas):v*=4**ki*math.factorial(ki)*pochhammer(a,ki)
 assert v.denominator==1;NORMS[n]=v.numerator

def inverse(x):
 y=[1/x[0]]
 for k in range(1,5):y.append(-sum((x[j]*y[k-j] for j in range(1,k+1)),acb(0))/x[0])
 return y
def log_nonconstant(x):
 y=[acb(0)]
 for k in range(1,5):y.append((k*x[k]-sum((j*y[j]*x[k-j] for j in range(1,k)),acb(0)))/(k*x[0]))
 return y
def mul(p,q):
 out={}
 for a,v in p.items():
  for b,w in q.items():
   key=tuple(x+y for x,y in zip(a,b))
   if sum(key)<=4:out[key]=out.get(key,acb(0))+v*w
 return out

def bound(theta,r1,r3):
 t=acb(0,theta).exp();t3=t**3
 polys=[];mass=arb(1);ell=[{ZERO:acb(0)} for _ in range(4)]
 for j in range(2):
  if j==0:z=[(r1*t+r3*t3*3**k)/math.factorial(k) for k in range(5)]
  else:z=[-t/math.factorial(k) for k in range(5)]
  xp=z.copy();xp[0]+=1;xm=[-a for a in z];xm[0]+=1
  ap=inverse(xp);am=inverse(xm);lp=log_nonconstant(xp);lm=log_nonconstant(xm)
  if j==0:
   reals=[ap[0].real,am[0].real]
   if not all(a>0 for a in reals):return arb.pos_inf()
   mass/=(abs(1-z[0]**2)*reals[0]*reals[1]).sqrt()
  else:
   reals=[arb(1)/2,arb(1)/2]
   mass*=((2/(2*theta.sin()).sqrt())**3)
  for k in range(1,5):
   ell[k-1][ZERO]+=-arb(math.factorial(k))/2*(lp[k]+lm[k])*(1 if j==0 else 3)
   for h,a in enumerate([ap,am]):
    key=tuple(int(i==2*j+h) for i in range(4));ell[k-1][key]=-arb(math.factorial(k))/2*a[k]/reals[h]
 bells=[{ZERO:acb(1)}]
 for n in range(1,5):
  p={}
  for k in range(1,n+1):
   for ix,v in mul(bells[n-k],ell[k-1]).items():p[ix]=p.get(ix,acb(0))+math.comb(n-1,k-1)*v
  bells.append(p)
 orth={key:acb(0) for key in KEYS}
 for ix,v in bells[4].items():
  for key,m in TRANS[ix]:orth[key]+=m*v
 norm2=sum((NORMS[k]*(v.real*v.real+v.imag*v.imag) for k,v in orth.items()),arb(0))
 if not norm2.is_finite():return arb.pos_inf()
 return mass*norm2.sqrt()

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--start',default='.6');ap.add_argument('--step',default='.002');ap.add_argument('--bits',type=int,default=160);ap.add_argument('--r1',default='.90292');ap.add_argument('--r3',default='-.09708');args=ap.parse_args()
 ctx.prec=args.bits;r1,r3=arb(args.r1),arb(args.r3);start,step=arb(args.start),arb(args.step);end=arb.pi()/2
 tic=time.time();cells=[];a=start;total=arb(0)
 while a<end:
  b=a+step
  if b>end:b=end
  theta=arb((a+b)/2,(b-a)/2);v=bound(theta,r1,r3)
  if not v.is_finite():raise RuntimeError('nonfinite arc enclosure')
  contribution=(b-a)*v.upper()**2;total+=contribution
  cells.append({'a':a.str(50),'b':b.str(50),'D4_bound':v.upper().str(40),'squared_integral_bound':contribution.str(40)})
  a=b
  if len(cells)%50==0:print('cells',len(cells),'theta',a.str(8),'seconds',time.time()-tic,flush=True)
 cumulative=[]
 for cutoff in ['.6','.65','.7','.75','.8','.85','.9','.95','1.0','1.05','1.1','1.2']:
  c=arb(cutoff)
  if c<start:continue
  subtotal=sum((arb(row['squared_integral_bound']) for row in cells if arb(row['b'])>c),arb(0))
  # Including a possible cell that straddles the requested start is conservative.
  energy=2/arb.pi()*subtotal
  cumulative.append({'start':cutoff,'quarter_arc_L2_squared_bound':energy.str(40),'quarter_arc_L2_norm_bound':energy.sqrt().str(40),'odd_degree_501_tail_Cauchy_bound':(energy/(14*arb(501)**7)).sqrt().str(40)})
 out={'status':'DIRECTED_COMPLEX_ARC_D4_COMPONENT_CERTIFIED','r1':args.r1,'r3':args.r3,'D_operator':'t*d/dt; theta fourth derivative has the same modulus','start':args.start,'end':'pi/2','step':args.step,'precision_bits':ctx.prec,'quarter_arc_L2_squared_bound':(2/arb.pi()*total).str(40),'cumulative_bounds':cumulative,'cells':cells,'seconds':time.time()-tic,'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Closed complex arc only; combine with an independent bound on the complementary real-endpoint arc before claiming a whole-circle tail.'}
 file=HERE/f'complex_arc_r1_{args.r1.replace(".","p")}_step_{args.step.replace(".","p")}.json';file.write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k!='cells'},indent=2))
if __name__=='__main__':main()
