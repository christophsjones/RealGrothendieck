"""Directed exact Gaussian polynomial integrals for the early C6 capped profile.

Certifies T[j,l]=sum_|alpha|=j j!/alpha! E[B[alpha,l](X)^2].
These are the endpoint chain-rule energies, before the Gaussian Gram factor.
"""
import json,math,hashlib
from pathlib import Path
from fractions import Fraction
from flint import arb,arb_poly,ctx
assert __debug__;ctx.prec=1536
HERE=Path(__file__).parent;SOURCE=HERE.parents[1]/'upper/refined_explicit_mixture.json'
raw=SOURCE.read_bytes();data=json.loads(raw)
def A(x):
 q=Fraction(x);return arb(q.numerator)/q.denominator
X=arb_poly([0,1]);hs=[arb_poly([1]),X]
for k in range(1,15):hs.append((X*hs[-1]-arb(k).sqrt()*hs[-2])*(1/arb(k+1).sqrt()))
rp=sum((A(c)*hs[d] for c,d in zip(data['threshold_definition']['r_coefficients_exact_decimal'],range(1,16,2))),arb_poly())
start=A('6.75');width=A('.25');der=rp;t=[]
for j in range(16):t.append(der(start)/math.factorial(j));der=der.derivative()
chi=[1,0,0,0,0,0,-462,1980,-3465,3080,-1386,252]
cp=[t[0]]+[arb(0)]*26
for j in range(1,16):
 for k,c in enumerate(chi):cp[j+k]+=j*t[j]*c/width**k/(j+k)
cap=arb_poly(cp)
def moments(a,L,n):
 sq=arb(2).sqrt();norm=(2*arb.pi()).sqrt();fa=(-a*a/2).exp()/norm;fb=(-(a+L)**2/2).exp()/norm
 J=[((a/sq).erfc()-((a+L)/sq).erfc())/2]
 if n:J.append(fa-fb-a*J[0])
 for k in range(1,n):J.append(k*J[k-1]-a*J[k]-L**k*fb)
 assert all(x>0 for x in J)
 return J
def bell(profile):
 r=[profile]
 for j in range(4):r.append(r[-1].derivative())
 B=[[arb_poly() for k in range(5)] for n in range(5)];B[0][0]=arb_poly([1])
 for n in range(1,5):
  for k in range(1,n+1):B[n][k]=sum((math.comb(n-1,j-1)*r[j]*B[n-j][k-1] for j in range(1,n-k+2)),arb_poly())
 return B
Bc=bell(rp);Bt=bell(cap);Jc=moments(arb(0),start,200);Jt=moments(start,width,200)
def integrate(poly,J):return sum((c*J[k] for k,c in enumerate(poly.coeffs())),arb(0))
U={(0,0,0):arb(1)}
for n in range(1,5):
 for k in range(1,n+1):
  for l in range(1,n+1):
   U[n,k,l]=(2*(integrate(Bc[n][k]*Bc[n][l],Jc)+integrate(Bt[n][k]*Bt[n][l],Jt))/math.factorial(n)) if (k+l)%2==0 else arb(0)
powers={(0,0,0):arb(1)}
for _ in range(3):
 nxt={}
 for a,va in powers.items():
  for b,vb in U.items():
   c=tuple(x+y for x,y in zip(a,b))
   if c[0]<=4:nxt[c]=nxt.get(c,arb(0))+va*vb
 powers=nxt
T={(j,l):math.factorial(j)*powers.get((j,l,l),arb(0)) for j in range(1,5) for l in range(1,j+1)}
uppers={(1,1):'.099',(2,1):'.221',(2,2):'.335',(3,1):'.36',(3,2):'3.38',(3,3):'4.88',(4,1):'1.11',(4,2):'24.3',(4,3):'124',(4,4):'110'}
for k,v in T.items():assert v>0 and v<A(uppers[k]),(k,v,uppers[k])
out=dict(status='RIGOROUS_FACTORIZED_CAPPED_ENERGIES_ACCEPTED',precision_bits=ctx.prec,source_sha256=hashlib.sha256(raw).hexdigest(),verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),profile='early C6 cap, start27/4,width1/4; see early_capped_profile_certificate.json',T={str(k):str(v) for k,v in T.items()},proved_T_upper={str(k):v for k,v in uppers.items()},scope='Chain-rule energies only; combine with positive Gaussian Gram factors and a circle bound for a tail certificate.')
(HERE/'factorized_energies_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2));print('FACTORIZED_CAPPED_ENERGIES_ACCEPTED=1')
