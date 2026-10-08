if not __debug__:
 raise RuntimeError("Certificate checks require Python without -O")
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np,itertools,struct,tempfile
from fractions import Fraction as F
from pathlib import Path
from graph import roof,cap,build
rng=np.random.default_rng(450733)
def value(B,a,c,s):return c+int(a@s)+int(s@B@s)//2
count=0
for n in range(1,8):
 for repeat in range(14):
  B=rng.integers(-50,51,size=(n,n),dtype=np.int64);B=np.ascontiguousarray(np.triu(B,1)+np.triu(B,1).T);a=rng.integers(-50,51,size=n,dtype=np.int64);c=int(rng.integers(-100,101));bd,labels,flow=roof(B,a,c)
  vals=[]
  for signs in itertools.product([-1,1],repeat=n):
   s=np.array(signs,dtype=np.int64);vv=value(B,a,c,s);assert vv<=bd
   t=np.where(labels!=0,labels,s).astype(np.int64);assert value(B,a,c,t)>=vv,('persistency',B,a,c,s,labels)
   vals.append(vv)
  if np.all(labels!=0):assert bd==max(vals)
  weights=rng.integers(1,100,size=n,dtype=np.int64)*rng.choice([-1,1],size=n);v=np.ascontiguousarray(weights,dtype=np.int64);V=sum(abs(int(t)) for t in v)
  p=int(rng.integers(1,1<<21));l,h=-21001,38007;CB,ca,cc,err=cap(B,a,c,v,p,l,h)
  for signs in itertools.product([-1,1],repeat=n):
   s=np.array(signs,dtype=np.int64);u=F(int(v@s),V);aux=F(p,1<<20)*(u-F(l,65536))*(F(h,65536)-u)*(1<<44)
   delta=F(value(CB,ca,cc,s)-value(B,a,c,s));assert abs(delta-aux)<=err
   if F(l,65536)<=u<=F(h,65536):assert value(CB,ca,cc,s)+err>=value(B,a,c,s)
  count+=1
# Independent reconstruction of the fixed-point polynomial at every sign vector.
N=3;H=rng.integers(-(1<<35),1<<35,size=(22,N),dtype=np.int64);M=rng.integers(-(1<<36),1<<36,size=(44,44),dtype=np.int64);M=np.ascontiguousarray(np.triu(M)+np.triu(M,1).T);ell=rng.integers(-(1<<35),1<<35,size=44,dtype=np.int64)
# Positive constant features are required only to select the additional cap coordinate.
H[0]=[1<<36,1<<35,1<<34]
raw=struct.pack('<8sIIII',b'KGHAT001',N,42,40,44)+M.astype('<i8').tobytes()+ell.astype('<i8').tobytes()+H.astype('<i8').tobytes()
p=Path(tempfile.gettempdir())/'kg_graph_unit.bin';p.write_bytes(raw);B,a,c,qe,v=build(p);H0=H.astype(object);M0=M.astype(object);e0=ell.astype(object)
for signs in itertools.product([-1,1],repeat=2*N):
 s=np.array(signs,dtype=np.int64);m=np.r_[H0@s[:N],H0@s[N:]];num=int(m@M0@m)+(int(e0@m)<<42);vv=F(num,1<<124);vq=F(value(B,a,c,s),1<<44);assert abs(vv-vq)<=qe
p.unlink();print('EXACT_GRAPH_TESTS_PASSED',count,'graphs, all assignments, cap bounds, persistency, and fixed-point construction')
