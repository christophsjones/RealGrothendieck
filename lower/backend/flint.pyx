# distutils: language = c++
# cython: language_level=3
from libcpp.string cimport string
from libc.stdlib cimport malloc,free
from fractions import Fraction
import math

cdef extern from "balls.hpp" namespace "verified":
    void set_precision(long) except +
    cdef cppclass Ball:
        Ball() except +
        Ball(const Ball&) except +
        void setstr(string) except +
        void setdouble(double) except +
        @staticmethod
        Ball pi() except +
        bint iszero() const
        bint finite() const
        Ball operator+(const Ball&) except +
        Ball operator-(const Ball&) except +
        Ball operator*(const Ball&) except +
        Ball operator/(const Ball&) except +
        Ball operator-() except +
        Ball absval() except +
        Ball sqrtval() except +
        Ball power(long) except +
        Ball expval() except +
        Ball erfval() except +
        Ball erfcval() except +
        Ball logval() except +
        Ball cosval() except +
        Ball sinval() except +
        Ball lower() except +
        Ball upper() except +
        Ball midval() except +
        Ball radval() except +
        Ball hull(const Ball&) except +
        bint lt(const Ball&) const
        bint le(const Ball&) const
        bint gt(const Ball&) const
        bint ge(const Ball&) const
        bint eq(const Ball&) const
        bint contains(const Ball&) const
        double asdouble() const
        string mantissa() except +
        long exponent() except +
        string text() except +
    cdef cppclass Matrix:
        int r,c
        Matrix() except +
        Matrix(int,int) except +
        Matrix(const Matrix&) except +
        Ball get(int,int) except +
        void set(int,int,const Ball&) except +
        Matrix transpose() except +
        Matrix add(const Matrix&) except +
        Matrix sub(const Matrix&) except +
        Matrix neg() except +
        Matrix scale(const Ball&) except +
        Matrix mul(const Matrix&) except +
    cdef cppclass QMatrix:
        int r,c
        QMatrix() except +
        QMatrix(int,int) except +
        QMatrix(const QMatrix&) except +
        string get(int,int) except +
        void set(int,int,const string&) except +
        QMatrix inverse() except +
        string determinant() except +

class Context:
    _prec=384
    @property
    def prec(self): return self._prec
    @prec.setter
    def prec(self,p):
        self._prec=int(p)
        set_precision(self._prec)
    @property
    def dps(self): return int(self._prec*math.log10(2))
    @dps.setter
    def dps(self,p): self.prec=max(64,math.ceil(float(p)/math.log10(2))+5)
ctx=Context()
__version__='independent-mpfr-endpoint-0.3-fma'
fmpq=Fraction
fmpz=int

cdef arb wrap(Ball x):
    cdef arb out=arb.__new__(arb)
    out.b=x
    return out

cdef class arb:
    cdef Ball b
    def __init__(self,x=0,rad=None):
        cdef arb a,r
        if isinstance(x,arb): self.b=(<arb>x).b
        elif isinstance(x,float): self.b.setdouble(x)
        else: self.b.setstr(str(Fraction(x)).encode())
        if rad is not None:
            r=arb(rad)
            self.b=(self.b-r.b).hull(self.b+r.b)
    def __add__(self,other):
        if isinstance(other,arb_mat):return NotImplemented
        return wrap(self.b+arb(other).b)
    def __radd__(self,other):return self+other
    def __sub__(self,other):return wrap(self.b-arb(other).b)
    def __rsub__(self,other):return wrap(arb(other).b-self.b)
    def __neg__(self):return wrap(-self.b)
    def __pos__(self):return self
    def __mul__(self,other):
        if isinstance(other,arb_mat):return NotImplemented
        return wrap(self.b*arb(other).b)
    def __rmul__(self,other):return self*other
    def __truediv__(self,other):return wrap(self.b/arb(other).b)
    def __rtruediv__(self,other):return wrap(arb(other).b/self.b)
    def __pow__(self,n,mod=None):
        if mod is not None or int(n)!=n:raise ValueError('only integer powers')
        return wrap(self.b.power(int(n)))
    def __abs__(self):return wrap(self.b.absval())
    def __float__(self):return self.b.asdouble()
    def __bool__(self):return not self.b.iszero()
    def __repr__(self):return self.b.text().decode()
    def __str__(self):return self.__repr__()
    def str(self,*args,**kwargs):return self.__repr__()
    def __richcmp__(self,other,int op):
        cdef arb z
        try:z=arb(other)
        except Exception:
            if op==2:return False
            if op==3:return True
            return NotImplemented
        if op==0:return self.b.lt(z.b)
        if op==1:return self.b.le(z.b)
        if op==2:return self.b.eq(z.b)
        if op==3:return self.b.lt(z.b) or self.b.gt(z.b)
        if op==4:return self.b.gt(z.b)
        if op==5:return self.b.ge(z.b)
    def sqrt(self):return wrap(self.b.sqrtval())
    def exp(self):return wrap(self.b.expval())
    def erf(self):return wrap(self.b.erfval())
    def erfc(self):return wrap(self.b.erfcval())
    def log(self):return wrap(self.b.logval())
    def cos(self):return wrap(self.b.cosval())
    def sin(self):return wrap(self.b.sinval())
    def lower(self):return wrap(self.b.lower())
    def upper(self):return wrap(self.b.upper())
    def mid(self):return wrap(self.b.midval())
    def rad(self):return wrap(self.b.radval())
    def union(self,other):return wrap(self.b.hull(arb(other).b))
    def contains(self,other):return self.b.contains(arb(other).b)
    def contains_zero(self):return self.contains(0)
    def is_finite(self):return self.b.finite()
    def is_exact(self):return self.b.lower().eq(self.b.upper())
    def is_zero(self):return self.b.iszero()
    def man_exp(self):return int(self.b.mantissa().decode()),int(self.b.exponent())
    @staticmethod
    def pi():return wrap(Ball.pi())

cdef arb_mat wrapmat(Matrix x):
    cdef arb_mat out=arb_mat.__new__(arb_mat)
    out.m=x
    return out

cdef class arb_mat:
    cdef Matrix m
    def __init__(self,x=None,y=None,z=None):
        cdef int i,j
        cdef arb v
        if x is None:return
        if isinstance(x,arb_mat):self.m=(<arb_mat>x).m;return
        if isinstance(x,fmpq_mat):x=x.tolist()
        if isinstance(x,int):
            if y is None:y=x
            self.m=Matrix(int(x),int(y))
            if z is not None:
                if len(z)!=x*y:raise ValueError('data shape')
                for i in range(x):
                    for j in range(y):v=arb(z[i*y+j]);self.m.set(i,j,v.b)
        else:
            r=len(x);c=len(x[0]) if r else 0
            self.m=Matrix(r,c)
            for i in range(r):
                if len(x[i])!=c:raise ValueError('data shape')
                for j in range(c):v=arb(x[i][j]);self.m.set(i,j,v.b)
    def nrows(self):return self.m.r
    def ncols(self):return self.m.c
    def __getitem__(self,key):
        i,j=key;return wrap(self.m.get(i,j))
    def __setitem__(self,key,val):
        cdef arb v=arb(val)
        i,j=key;self.m.set(i,j,v.b)
    def __add__(self,other):
        if not isinstance(other,arb_mat):return NotImplemented
        return wrapmat(self.m.add((<arb_mat>other).m))
    def __radd__(self,other):
        if other==0:return self
        return self+other
    def __sub__(self,other):return wrapmat(self.m.sub((<arb_mat>other).m))
    def __neg__(self):return wrapmat(self.m.neg())
    def __mul__(self,other):
        if isinstance(other,arb_mat):return wrapmat(self.m.mul((<arb_mat>other).m))
        return wrapmat(self.m.scale(arb(other).b))
    def __rmul__(self,other):return self*other
    def __truediv__(self,other):return self*(arb(1)/other)
    def transpose(self):return wrapmat(self.m.transpose())
    def tolist(self):return [[self[i,j] for j in range(self.m.c)] for i in range(self.m.r)]
    def __repr__(self):return repr(self.tolist())
    def __str__(self):return repr(self)

cdef fmpq_mat wrapqmat(QMatrix x):
    cdef fmpq_mat out=fmpq_mat.__new__(fmpq_mat)
    out.m=x
    return out

cdef class fmpq_mat:
    cdef QMatrix m
    def __init__(self,x=None,y=None):
        cdef int i,j
        if x is None:return
        if isinstance(x,int):self.m=QMatrix(int(x),int(y if y is not None else x));return
        r=len(x);c=len(x[0]) if r else 0;self.m=QMatrix(r,c)
        for i in range(r):
            if len(x[i])!=c:raise ValueError('shape')
            for j in range(c):self.m.set(i,j,str(Fraction(x[i][j])).encode())
    def nrows(self):return self.m.r
    def ncols(self):return self.m.c
    def __getitem__(self,key):
        i,j=key;return Fraction(self.m.get(i,j).decode())
    def __setitem__(self,key,val):
        i,j=key;self.m.set(i,j,str(Fraction(val)).encode())
    def inv(self):return wrapqmat(self.m.inverse())
    def det(self):return Fraction(self.m.determinant().decode())
    def tolist(self):return [[self[i,j] for j in range(self.m.c)] for i in range(self.m.r)]
