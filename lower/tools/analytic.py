if not __debug__:
 raise RuntimeError("Certificate checks require Python without -O")
"""Directed Gaussian part of the finite-sign certificate.

All finite decimal source values and all grid endpoints are exact rationals.
The imported flint interface may be Arb or the bundled, identified MPFR backend.
"""
from fractions import Fraction as F
from pathlib import Path
import math,json,hashlib,struct,time
from flint import arb,arb_mat,ctx
ctx.prec=384

D=21;NFEATURE=22;BSHIFT=42;MSHIFT=40;QSHIFT=44

def A(x=0):return arb(F(x))
def endpoint_fraction(x,upper=True):
 y=x.upper() if upper else x.lower();m,e=y.man_exp()
 return F(m)*F(2)**e

def upper(x,denpower=90):
 q=endpoint_fraction(x);n=-((-q.numerator*(1<<denpower))//q.denominator)
 return F(n,1<<denpower)

def rounded(q,shift):
 q=F(q)*(1<<shift)
 return (2*q.numerator+q.denominator)//(2*q.denominator) if q>=0 else -rounded(-q/(1<<shift),shift)

def ceil_fraction(q,shift=80):
 q=F(q);return F(-((-q.numerator*(1<<shift))//q.denominator),1<<shift)

def hermites(x,n):
 out=[A(1)]
 if n:out.append(x)
 for j in range(1,n):out.append((x*out[-1]-A(j).sqrt()*out[-2])/A(j+1).sqrt())
 return out

class Gaussian:
 def __init__(self):
  self.norm=(2*arb.pi()).sqrt();self.rt2=A(2).sqrt();self.ends={};self.mom={};self.grams={};self.product=[]
  for j in range(22):
   self.product.append([])
   for k in range(j+1):
    self.product[j].append([(j+k-2*l,A(math.factorial(l)*math.comb(j,l)*math.comb(k,l))*A(F(math.factorial(j+k-2*l),math.factorial(j)*math.factorial(k))).sqrt()) for l in range(min(j,k)+1)])
 def endpoint(self,x):
  if x is None:return None
  x=F(x)
  if x not in self.ends:
   xx=A(x);ph=(-xx*xx/2).exp()/self.norm;h=hermites(xx,41)
   self.ends[x]=((-xx/self.rt2).erfc()/2,[ph*z for z in h],ph)
  return self.ends[x]
 def interval(self,lo,hi):
  key=(lo,hi)
  if key not in self.mom:
   l,r=self.endpoint(lo),self.endpoint(hi)
   out=[(A(1) if r is None else r[0])-(A(0) if l is None else l[0])]
   for j in range(1,43):out.append(((A(0) if l is None else l[1][j-1])-(A(0) if r is None else r[1][j-1]))/A(j).sqrt())
   assert out[0]>0
   self.mom[key]=out
  return self.mom[key]
 def gram(self,lo,hi):
  key=(lo,hi)
  if key not in self.grams:
   mu=self.interval(lo,hi);v=[[A(0) for k in range(22)] for j in range(22)]
   for j in range(22):
    for k in range(j+1):
     z=sum((c*mu[l] for l,c in self.product[j][k]),A(0));v[j][k]=z;v[k][j]=z
   self.grams[key]=arb_mat(v)
  return self.grams[key]

class Common:
 def __init__(self):
  self.gauss=Gaussian();self.ends=[F(j,2) for j in range(-18,19)]
  bounds=[None]+self.ends+[None];self.intervals=list(zip(bounds[:-1],bounds[1:]));g=self.gauss
  self.mass=[g.interval(lo,hi)[0] for lo,hi in self.intervals]
  self.gram=[g.gram(lo,hi) for lo,hi in self.intervals]
  self.flat=arb_mat([[G[j,k] for j in range(22) for k in range(22)] for G in self.gram])
  self.derivative=[]
  for G in self.gram[1:-1]:
   self.derivative.append(arb_mat([[A(j*(j-1)*(j-2)*k*(k-1)*(k-2)).sqrt()*G[j-3,k-3] if min(j,k)>=3 else A(0) for k in range(22)] for j in range(22)]))
  self.tailgram=self.gram[0]+self.gram[-1];self.tailmass=self.mass[0]+self.mass[-1]


def read_row(path,index):
 data=json.loads(Path(path).read_text(),parse_float=str);r=data['rows'][index]
 H=[[F(0) for _ in range(22)] for _ in range(22)];K=[[F(0) for _ in range(22)] for _ in range(22)]
 for X,par,B in [(H,0,r['matrices'][0]),(H,1,r['matrices'][1]),(K,0,r['matrices'][2]),(K,1,r['matrices'][3])]:
  for j in range(11):
   for k in range(11):X[2*j+par][2*k+par]=F(B[j][k])
 AA=[[(H[j][k]+K[j][k])/4 for k in range(22)] for j in range(22)]
 BB=[[(H[j][k]-K[j][k])/4 for k in range(22)] for j in range(22)]
 M=[AA[j]+BB[j] for j in range(22)]+[BB[j]+AA[j] for j in range(22)]
 Lh=[F(0) for _ in range(22)];Lk=Lh.copy();Lh[1]=F(r['t'])
 for j in range(11):Lh[2*j+1]-=2*F(r['lh'][j]);Lk[2*j+1]-=2*F(r['lk'][j])
 ell=[(Lh[j]+Lk[j])/2 for j in range(22)]+[(Lh[j]-Lk[j])/2 for j in range(22)]
 assert all(M[i][j]==M[j][i] for i in range(44) for j in range(44))
 return M,ell

def check_pd(M):
 n=len(M);L=[[A(0) for j in range(n)] for i in range(n)];d=[]
 for i in range(n):
  for j in range(i):L[i][j]=(A(M[i][j])-sum((L[i][k]*d[k]*L[j][k] for k in range(j)),A(0)))/d[j]
  di=A(M[i][i])-sum((L[i][k]*L[i][k]*d[k] for k in range(i)),A(0))
  assert di>0,('M_not_certified_positive_definite',i,str(di))
  d.append(di);L[i][i]=A(1)
 return str(min(endpoint_fraction(x,False) for x in d))

def trace(mat):return sum((mat[i,i] for i in range(mat.nrows())),A(0))
def sqrt_checked(v):
 assert v>=0,('negative_square_norm',str(v))
 return v.sqrt()

def curvature(common,M,ell):
 R=common.derivative;mass=common.mass;vals=[A(0) for _ in R]
 for c in range(2):
  v=arb_mat([[ell[c*22+j]] for j in range(22)])
  for q,r in enumerate(R):vals[q]+=sqrt_checked(mass[q+1]*(v.transpose()*r*v)[0,0])
 for bl in [0,22]:
  mat=arb_mat([row[bl:bl+22] for row in M[:22]])
  YY=[mat*r*mat.transpose() for r in R]
  flat=arb_mat([[Y[i,j] for Y in YY] for i in range(22) for j in range(22)])
  products=common.flat*flat
  for q in range(len(R)):
   vals[q]+=4*sum((sqrt_checked(mass[q+1]*mass[k]*products[k,q]) for k in range(len(mass))),A(0))
 tail=A(0)
 for c in range(2):
  v=arb_mat([[ell[c*22+j]] for j in range(22)])
  tail+=sqrt_checked(common.tailmass*(v.transpose()*common.tailgram*v)[0,0])
 for bl in [0,22]:
  mat=arb_mat([row[bl:bl+22] for row in M[:22]])
  tail+=4*sqrt_checked(common.tailmass*trace(mat*common.tailgram*mat.transpose()))
 return vals,tail


def grid_nodes(counts):
 assert len(counts)==36 and all(type(n) is int and 1<=n<=100000 for n in counts)
 assert counts==counts[::-1], 'The frozen partition must be reflection symmetric'
 out=[]
 for j,n in enumerate(counts):
  a=F(j-18,2);out.extend(a+F(k,4*n) for k in range(2*n))
 out.append(F(9)); assert out[0]==-9 and all(a<b for a,b in zip(out[:-1],out[1:]))
 return out

def integral_error(common,counts,vals,tail):
 err=tail
 for j,(n,v) in enumerate(zip(counts,vals)):
  delta=F(1,2*n);r=max(abs(common.ends[j]),abs(common.ends[j+1]))
  err+=A(r*delta).exp()*A(delta**3/96)*v
 return upper(err)

def hat_inputs(common,nodes,M,ell):
 """Features are integrated positive Bernstein bases times control functionals.

 Each pair of cells [a,mid],[mid,b] is one quadratic interpolation panel.
 Identical endpoint functionals of adjacent panels are combined.
 Reflection is exact, and the rounding error is propagated uniformly.
 """
 assert len(nodes)%4==1 and nodes==[-x for x in nodes[::-1]]
 g=common.gauss;N=len(nodes);mid=N//2
 features=[[A(0) for k in range(22)] for x in nodes]
 total=A(0)
 for j in range(mid,N-2,2):
  a,x,b=nodes[j:j+3];assert x==(a+b)/2
  h=A(b-a);mu=g.interval(a,b);I0,I1,I2=mu[0],mu[1],A(2).sqrt()*mu[2]+mu[0]
  w0=(A(b*b)*I0-2*A(b)*I1+I2)/(h*h)
  w1=2*(-A(a*b)*I0+A(a+b)*I1-I2)/(h*h)
  w2=(A(a*a)*I0-2*A(a)*I1+I2)/(h*h)
  assert w0>0 and w1>0 and w2>0
  assert (w0+w1+w2-I0).contains_zero()
  total+=2*(w0+w1+w2)
  va,vm,vb=hermites(A(a),21),hermites(A(x),21),hermites(A(b),21)
  for k in range(22):
   features[j][k]+=w0*va[k]
   features[j+1][k]+=w1*(2*vm[k]-(va[k]+vb[k])/2)
   features[j+2][k]+=w2*vb[k]
 # Reflect the half-line. The shared central endpoint has two equal weights.
 for k in range(22):
  features[mid][k]=features[mid][k]*2
  for j in range(mid+1,N):features[N-1-j][k]=features[j][k]*((-1)**k)
 assert (total-g.interval(-9,9)[0]).contains_zero()
 Bint=[[0 for x in nodes] for k in range(22)];er=[A(0) for k in range(22)];bs=[0]*22
 for j in range(N):
  for k in range(22):
   bb=features[j][k];q=rounded(endpoint_fraction(bb.mid()),BSHIFT);Bint[k][j]=q;bs[k]+=abs(q);er[k]+=abs(bb-A(F(q,1<<BSHIFT)))
 for k in range(22):assert Bint[k]==[(-1)**k*x for x in Bint[k][::-1]]
 b=[F(z,1<<BSHIFT) for z in bs]*2;e=[upper(z,80) for z in er]*2
 mint=[[rounded(v,MSHIFT) for v in row] for row in M];eint=[rounded(v,MSHIFT) for v in ell]
 ma=[[F(v,1<<MSHIFT) for v in row] for row in mint];ea=[F(v,1<<MSHIFT) for v in eint]
 error=F(0)
 for i in range(44):
  error+=abs(ea[i])*e[i]+abs(ell[i]-ea[i])*(b[i]+e[i])
  for j in range(44):
   error+=(2*b[i]*e[j]+e[i]*e[j])*abs(ma[i][j])+(b[i]+e[i])*(b[j]+e[j])*abs(M[i][j]-ma[i][j])
 return Bint,mint,eint,ceil_fraction(error),[str(z) for z in e[:22]],str(upper(total))


def save_input(path,Bint,mint,eint):
 n=len(Bint[0]);data=struct.pack('<8sIIII',b'KGHAT001',n,BSHIFT,MSHIFT,QSHIFT)
 data+=struct.pack('<'+str(44*44)+'q',*[x for row in mint for x in row]);data+=struct.pack('<44q',*eint)
 data+=struct.pack('<'+str(22*n)+'q',*[x for row in Bint for x in row])
 Path(path).write_bytes(data);return hashlib.sha256(data).hexdigest()

