"""Explicit radix-two contraction with an outward forward-error budget."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import json,time,math
from pathlib import Path
import numpy as np
from flint import arb,ctx
ctx.prec=256
HERE=Path(__file__).resolve().parent;meta=json.loads((HERE/'certified_b_inputs.json').read_text())
data=np.load(HERE/'certified_b_inputs.npz');B=data['B'];U=data['U'].copy();factor=float(data['factor']);D=meta['D'];S=len(B);N=2048
assert 3*D<N and N&(N-1)==0
unit=arb(2)**-53;twiddle=np.empty(N//2,dtype=np.complex128)
for k in range(N//2):
 angle=-2*arb.pi()*k/N;re=angle.cos();im=angle.sin();rf=float(re.mid());imf=float(im.mid())
 assert abs(re-arb(rf))<arb(2)**-48 and abs(im-arb(imf))<arb(2)**-48
 twiddle[k]=complex(rf,imf)
rev=np.array([int(f'{k:011b}'[::-1],2) for k in range(N)])

def fft(a,inverse=False):
 a=a[...,rev].copy()
 for width in [2,4,8,16,32,64,128,256,512,1024,2048]:
  shape=a.shape;v=a.reshape(*shape[:-1],N//width,width);half=width//2
  even=v[...,:half].copy();odd=v[...,half:].copy();w=twiddle[::N//width]
  if inverse:w=w.conj()
  odd*=w;v[...,:half]=even+odd;v[...,half:]=even-odd
 if inverse:a/=N
 return a

tic=time.time();G=np.empty((S,S,D+1))
for start in range(0,S,8):
 stop=min(start+8,S);g=np.zeros((stop-start,S,N),dtype=np.complex128)
 g[...,:D+1]=B[start:stop,None,:]*B[None,:,:]
 z=fft(g);z=z*z*z;z=fft(z,True)
 G[start:stop]=z[...,:D+1].real
 if start%128==0:print('FFT rows',start,'seconds',time.time()-tic,flush=True)

# Each radix-two stage has O(u) normwise error and exact operator norm sqrt2.
# With twiddle component errors <=2^-48=32u, two FFTs and a cube have
# coefficient-l1 error at most 10000*N*u*M^3 for M=||g||1, N=2048.
# The deliberately loose constant includes input products and subnormal loss.
bnorm=arb(meta['B_norm_bound']);conv_error=10000*N*unit*bnorm**6+arb('1e-290')
usq=arb(meta['U_l1_squares_sum']);float_factor=arb(factor)
fft_plane_error=float_factor*usq*conv_error
pruned=arb(0)
for a in range(D+1):
 small=U[a]<1e-18
 pruned+=sum(abs(arb(float(x))) for x in U[a,small]);U[a,small]=0
pruning_error=float_factor*(2*arb(meta['U_l1_sum'])*pruned+pruned*pruned)*bnorm**6
C=np.zeros((D+1,D+1))
for a in range(D+1):
 ix=np.flatnonzero(U[a]);uu=U[a,ix]
 sub=G[ix][:,ix,:D+1-a]
 tmp=np.einsum('s,stb->tb',uu,sub,optimize=False)
 C[a,:D+1-a]=factor*np.einsum('t,tb->b',uu,tmp,optimize=False)
 C[a,(np.arange(D+1)+a)%2==0]=0
 if a%50==0:print('C row',a,'seconds',time.time()-tic,flush=True)
assert np.isfinite(C).all()
# The two ordinary contractions have at most 2S^2+8S rounded operations
# on any contributing arithmetic path; sum absolute exact coefficients <=M^3.
k=2*S*S+8*S;gamma=k*unit/(1-k*unit)
contraction_error=gamma*float_factor*usq*(bnorm**6+conv_error)
total=(arb(meta['plane_l1_quadrature_error'])+arb(meta['plane_l1_U_input_error'])*(bnorm**6+conv_error)+fft_plane_error+pruning_error+contraction_error)
# Orthogonal parity deletion and projection onto the nonnegative orthant
# cannot increase distance from the exact nonnegative coefficient plane.
C=np.maximum(C,0)
np.save(HERE/'certified_plane.npy',C)
out={'status':'DIRECTED_FINITE_NATIVE_COEFFICIENT_PLANE','D':D,'candidate_sha256':meta['candidate_sha256'],'sum_computed_C':float(C.sum()),'plane_l1_error':total.str(40),'fft_plane_l1_error':fft_plane_error.str(40),'contraction_l1_error':contraction_error.str(40),'pruning_l1_error':pruning_error.str(40),'whole_kernel_cutoff_error':meta['whole_kernel_cutoff_error'],'seconds':time.time()-tic,'arithmetic_model':'IEEE754 binary64 round-to-nearest; explicit radix-two butterflies, Arb-certified twiddles, ordinary optimize=False einsum sums. See source for forward-error constants.'}
(HERE/'certified_plane.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2),flush=True)
