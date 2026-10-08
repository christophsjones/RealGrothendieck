"""Independent exact polynomial identities and directed Gaussian checks."""
if not __debug__:raise RuntimeError('Run without -O')
from fractions import Fraction as F
from math import factorial,comb
from analytic import Common,A,grid_nodes,hat_inputs,upper
from flint import arb_mat

def add(x,y):
 z=[0]*max(len(x),len(y))
 for j,a in enumerate(x):z[j]+=a
 for j,a in enumerate(y):z[j]+=a
 return z

def scale(x,a):return [a*t for t in x]
def mul(x,y):
 z=[0]*(len(x)+len(y)-1)
 for j,a in enumerate(x):
  for k,b in enumerate(y):z[j+k]+=a*b
 return z
he=[[1],[0,1]]
for j in range(1,42):he.append(add([0]+he[j],scale(he[j-1],-j)))
for j in range(22):
 for k in range(22):
  rhs=[0]*(j+k+1)
  for l in range(min(j,k)+1):rhs=add(rhs,scale(he[j+k-2*l],factorial(l)*comb(j,l)*comb(k,l)))
  assert mul(he[j],he[k])==rhs
C=Common();mom=C.gauss.interval(None,None)
assert mom[0].contains(1) and all(mom[j].contains(0) for j in range(1,43))
G=C.gauss.gram(None,None)
for i in range(22):
 for j in range(22):assert G[i,j].contains(int(i==j))
assert sum(C.mass,A(0)).contains(1)
Gsum=arb_mat(22,22)
for T in C.gram:Gsum=Gsum+T
for i in range(22):
 for j in range(22):assert Gsum[i,j].contains(int(i==j))
# Exact third-derivative identities and the Peano constant 1/96.
def derivative(p):return [(j+1)*p[j+1] for j in range(len(p)-1)]
for j in range(3,22):
 assert derivative(derivative(derivative(he[j])))==scale(he[j-3],j*(j-1)*(j-2))
# For 0<=t<=1/2, integral |K(s,t)| ds = t^2/8-t^3/6.
left=add(scale(mul(mul([F(1,2),-1],[F(1,2),-1]),[F(1,2),-1]),F(1,3)),scale(mul(mul([1,-1],[1,-1]),[1,-1]),F(-1,6)))
left=add(left,scale(mul([1,-1],[1,-1]),F(1,8)))
assert left==[0,0,F(1,8),F(-1,6)]
# 1/96 - integral |K| = (1-2t)^2(1+4t)/96 >= 0.
assert add([F(1,96)],scale(left,-1))==scale(mul(mul([1,-2],[1,-2]),[1,4]),F(1,96))
assert sum(v*F(1,2)**j for j,v in enumerate(left))==F(1,96)
for G,mu in zip(C.derivative,C.mass[1:-1]):
 assert (G[3,3]-6*mu).contains_zero()
 for j in range(3):
  for k in range(22):assert G[j,k].contains_zero()
# Directed Bernstein weights reproduce every quadratic exactly.
nodes=grid_nodes([2]*36)
for j in range(0,len(nodes)-2,2):
 a,m,b=nodes[j:j+3];assert m==(a+b)/2
 mu=C.gauss.interval(a,b);h=A(b-a)
 I0,I1,I2=mu[0],mu[1],A(2).sqrt()*mu[2]+mu[0]
 w0=(A(b*b)*I0-2*A(b)*I1+I2)/(h*h)
 w1=2*(-A(a*b)*I0+A(a+b)*I1-I2)/(h*h)
 w2=(A(a*a)*I0-2*A(a)*I1+I2)/(h*h)
 assert w0>0 and w1>0 and w2>0
 assert (w0+w1+w2-I0).contains_zero()
 assert (w0*A(a)+w1*A(m)+w2*A(b)-I1).contains_zero()
 assert (w0*A(a*a)+w1*A(a*b)+w2*A(b*b)-I2).contains_zero()
 cubic=w0*A(a**3)+w1*A(2*m**3-(a**3+b**3)/2)+w2*A(b**3)
 exact_cubic=A(6).sqrt()*mu[3]+3*mu[1]
 bound=A(max(abs(a),abs(b))*(b-a)).exp()*A((b-a)**3)*I0/16
 assert bound-abs(cubic-exact_cubic)>0
print('ANALYTIC_IDENTITY_TESTS_PASSED: 484 Hermite products; third derivatives; quadratic Bernstein weights; Peano constant 1/96; Gaussian Grams')
