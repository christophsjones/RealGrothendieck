if not __debug__:
 raise RuntimeError("Certificate checks require Python without -O")
"""Typed wrapper for exact integer graph operations."""
import ctypes as c,os
from pathlib import Path
import numpy as np
import struct,json,hashlib
from fractions import Fraction as F
LIB=c.CDLL(os.environ.get('KG_EXACT_GRAPH_LIBRARY',str(Path(__file__).with_name('exact_graph.so'))))
I64=np.ctypeslib.ndpointer(dtype=np.int64,flags='C_CONTIGUOUS');I32=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
LIB.build_graph.argtypes=[c.c_int,I64,I64,I64,I64,I64,c.POINTER(c.c_int64),c.c_char_p,c.c_int,c.c_char_p,c.c_int];LIB.build_graph.restype=c.c_int
LIB.exact_roof.argtypes=[c.c_int,I64,I64,c.c_int64,c.POINTER(c.c_int64),c.POINTER(c.c_int64),I32,c.c_char_p,c.c_int];LIB.exact_roof.restype=c.c_int
LIB.cap_graph.argtypes=[c.c_int,I64,I64,c.c_int64,I64,c.c_int64,c.c_int64,c.c_int64,I64,I64,c.POINTER(c.c_int64),c.POINTER(c.c_int64),c.c_char_p,c.c_int];LIB.cap_graph.restype=c.c_int
LIB.validate_graph.argtypes=[c.c_int,I64,I64,c.c_int64,c.POINTER(c.c_int64),c.c_char_p,c.c_int];LIB.validate_graph.restype=c.c_int

def validate(B,a,constant):
 if B.shape!=(len(a),len(a)) or B.dtype!=np.int64 or a.dtype!=np.int64:raise ValueError('graph shape or dtype')
 budget=c.c_int64();err=c.create_string_buffer(1024)
 if LIB.validate_graph(len(a),B,a,constant,c.byref(budget),err,len(err)):raise ValueError(err.value.decode())
 return budget.value

def build(path,out=None):
 raw=Path(path).read_bytes();magic,N,bs,ms,qs=struct.unpack_from('<8sIIII',raw)
 if (magic,bs,ms,qs)!=(b'KGHAT001',42,40,44):raise ValueError('input header')
 offset=24;M=np.frombuffer(raw,dtype='<i8',count=1936,offset=offset).copy().reshape(44,44);offset+=1936*8
 ell=np.frombuffer(raw,dtype='<i8',count=44,offset=offset).copy();offset+=44*8
 H=np.frombuffer(raw,dtype='<i8',count=22*N,offset=offset).copy().reshape(22,N)
 if len(raw)!=offset+22*N*8:raise ValueError('trailing input')
 n=2*N;B=np.empty((n,n),dtype=np.int64);a=np.empty(n,dtype=np.int64);cc=c.c_int64();qe=c.create_string_buffer(256);err=c.create_string_buffer(1024)
 if LIB.build_graph(N,H,M,ell,B,a,c.byref(cc),qe,len(qe),err,len(err)):raise ValueError(err.value.decode())
 # A rational scalar coordinate, normalized to have range [-1,1].
 weights=((H[0]+(1<<17))>>18).astype(np.int64)
 if np.any(weights<0) or not int(weights.sum()):raise ValueError('bad weights')
 v=np.ascontiguousarray(np.r_[weights,-weights])
 validate(B,a,int(cc.value))
 result=(B,a,int(cc.value),F(int(qe.value),1<<124),v)
 if out is not None:np.savez(out,B=B,a=a,constant=cc.value,error_numerator=str(int(qe.value)),v=v,input_sha256=hashlib.sha256(raw).hexdigest())
 return result

def load(path):
 z=np.load(path);return (np.ascontiguousarray(z['B']),np.ascontiguousarray(z['a']),int(z['constant']),F(int(z['error_numerator']),1<<124),np.ascontiguousarray(z['v']))

def roof(B,a,constant=0):
 n=len(a);bound=c.c_int64();flow=c.c_int64();labels=np.empty(n,dtype=np.int32);err=c.create_string_buffer(1024)
 if LIB.exact_roof(n,B,a,constant,c.byref(bound),c.byref(flow),labels,err,len(err)):raise ValueError(err.value.decode())
 return bound.value,labels,flow.value

def cap(B,a,constant,v,p,lo,hi):
 n=len(a);CB=np.empty_like(B);ca=np.empty_like(a);cc=c.c_int64();errunits=c.c_int64();err=c.create_string_buffer(1024)
 if LIB.cap_graph(n,B,a,constant,v,p,lo,hi,CB,ca,c.byref(cc),c.byref(errunits),err,len(err)):raise ValueError(err.value.decode())
 validate(CB,ca,cc.value)
 return CB,ca,cc.value,errunits.value

def symmetry(B,a,v):
 n=len(a);N=n//2;p=np.r_[np.arange(N-1,-1,-1),np.arange(n-1,N-1,-1)]
 return bool(np.array_equal(a,-a[p]) and np.array_equal(v,v[p]) and np.array_equal(B,B[np.ix_(p,p)]))

if __name__=='__main__':
 import argparse,time
 ap=argparse.ArgumentParser();ap.add_argument('input');ap.add_argument('--out');args=ap.parse_args();t=time.time();B,a,cc,e,v=build(args.input,args.out);print('BUILD',len(a),float(e),time.time()-t,flush=True);print('SYMMETRY',symmetry(B,a,v),flush=True)
 bd,r,fl=roof(B,a,cc);print('ROOF',float(F(bd,1<<44)),int(np.sum(r==0)),time.time()-t,flush=True)
