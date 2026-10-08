from flint import arb,arb_mat,fmpq_mat,ctx
from fractions import Fraction as F
import random,json,math,ctypes
ctx.prec=384

def val(x):
 m,e=x.man_exp();return F(m)*F(2)**e

def endpoints(x):return val(x.lower()),val(x.upper())
def enc(x,lo,hi):
 l,h=endpoints(x)
 assert l<=lo<=hi<=h,(x,lo,hi)

rng=random.Random(20260907)
for _ in range(10000):
 a,b,c,d=[F(rng.randrange(-10**12,10**12),rng.randrange(1,10**12)) for k in range(4)]
 a,b=sorted([a,b]);c,d=sorted([c,d]);x=arb(a).union(arb(b));y=arb(c).union(arb(d))
 enc(x+y,a+c,b+d);enc(x-y,a-d,b-c)
 vv=[a*c,a*d,b*c,b*d];enc(x*y,min(vv),max(vv));enc(-x,-b,-a)
 enc(abs(x),0 if a<=0<=b else min(abs(a),abs(b)),max(abs(a),abs(b)))
 enc(x**2,0 if a<=0<=b else min(a*a,b*b),max(a*a,b*b))
 if not c<=0<=d:
  vv=[a/c,a/d,b/c,b/d];enc(x/y,min(vv),max(vv))
 if a>=0:
  l,h=endpoints(x.sqrt());assert l*l<=a and h*h>=b
for n in [1,3,7,12]:
 A=[[F(rng.randrange(-100,100),rng.randrange(1,100)) for j in range(n)] for i in range(n)]
 Q=fmpq_mat(A);I=Q.inv().tolist()
 assert all(sum(A[i][k]*I[k][j] for k in range(n))==int(i==j) for i in range(n) for j in range(n))
 B=arb_mat(A)*arb_mat(I)
 assert all(B[i,j].contains(int(i==j)) for i in range(n) for j in range(n))
# Independently enclose every interval product before summing, including all
# endpoint-sign configurations used by the native FMA matrix kernel.
for n in [1,2,5,12]:
 for trial in range(10):
  intervals=[];mats=[]
  for channel in range(2):
   vals=[];aa=arb_mat(n,n)
   for i in range(n):
    row=[]
    for j in range(n):
     lo,hi=sorted(F(rng.randrange(-1000,1000),rng.randrange(1,1000)) for k in range(2))
     if (i+j+trial)%11==0:lo=hi=F(0)
     row.append((lo,hi));aa[i,j]=arb(lo).union(arb(hi))
    vals.append(row)
   intervals.append(vals);mats.append(aa)
  cc=mats[0]*mats[1]
  for i in range(n):
   for j in range(n):
    low=high=F(0)
    for k in range(n):
     vv=[x*y for x in intervals[0][i][k] for y in intervals[1][k][j]]
     low+=min(vv);high+=max(vv)
    enc(cc[i,j],low,high)

# Rational Taylor enclosures of exp(1), exp(-1), and pi via Machin's identity.
s=sum((F(1,math.factorial(k)) for k in range(121)),F(0));tail=F(1,math.factorial(121))*F(122,121)
e=arb(1).exp();l,h=endpoints(e)
assert l<=s+tail and s<=h and h-l<F(1,10**110)
# Strict overlap with independently proven narrow intervals; enclosures themselves use MPFR-directed functions.
def atan_rec(x,N):
 s=sum(((-1)**j*x**(2*j+1)/F(2*j+1) for j in range(N)),F(0));t=x**(2*N+1)/F(2*N+1)
 return (s,s+t) if N%2==0 else (s-t,s)
a,b=atan_rec(F(1,5),100);c,d=atan_rec(F(1,239),100);pilo,pihi=16*a-4*d,16*b-4*c
l,h=endpoints(arb.pi());assert l<=pihi and h>=pilo
for x in [F(-7),F(-1),F(0),F(1),F(7)]:
 assert (arb(x).erfc()+arb(x).erf()).contains(1)
 assert (arb(x).exp()*arb(-x).exp()).contains(1)
 assert (arb(x).sin()**2+arb(x).cos()**2).contains(1)
print('MPFR_ENDPOINT_BACKEND_SELF_TESTS_PASSED')
