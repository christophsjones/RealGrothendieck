#pragma once
#include <gmpxx.h>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <sstream>
extern "C" {
typedef long mpfr_prec_t; typedef long mpfr_exp_t; typedef int mpfr_sign_t;
typedef struct {mpfr_prec_t _mpfr_prec;mpfr_sign_t _mpfr_sign;mpfr_exp_t _mpfr_exp;mp_limb_t *_mpfr_d;} __mpfr_struct;
typedef __mpfr_struct mpfr_t[1]; typedef __mpfr_struct* mpfr_ptr; typedef const __mpfr_struct* mpfr_srcptr;
typedef enum {MPFR_RNDN=0,MPFR_RNDZ=1,MPFR_RNDU=2,MPFR_RNDD=3,MPFR_RNDA=4,MPFR_RNDF=5} mpfr_rnd_t;
int mpfr_fma(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);
int mpfr_dot(mpfr_ptr,const mpfr_ptr[],const mpfr_ptr[],unsigned long,mpfr_rnd_t);
void mpfr_init2(mpfr_ptr,mpfr_prec_t);void mpfr_clear(mpfr_ptr);int mpfr_set(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_set_si(mpfr_ptr,long,mpfr_rnd_t);int mpfr_set_d(mpfr_ptr,double,mpfr_rnd_t);int mpfr_set_q(mpfr_ptr,mpq_srcptr,mpfr_rnd_t);int mpfr_set_z(mpfr_ptr,mpz_srcptr,mpfr_rnd_t);int mpfr_set_str(mpfr_ptr,const char*,int,mpfr_rnd_t);
int mpfr_add(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_sub(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_mul(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_div(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_neg(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_abs(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_sqr(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_sqrt(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_exp(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_erf(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_erfc(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_log(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_cos(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_sin(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_pow_si(mpfr_ptr,mpfr_srcptr,long,mpfr_rnd_t);int mpfr_const_pi(mpfr_ptr,mpfr_rnd_t);int mpfr_cmp(mpfr_srcptr,mpfr_srcptr);int mpfr_cmp_si(mpfr_srcptr,long);int mpfr_number_p(mpfr_srcptr);int mpfr_zero_p(mpfr_srcptr);double mpfr_get_d(mpfr_srcptr,mpfr_rnd_t);mpfr_exp_t mpfr_get_z_2exp(mpz_ptr,mpfr_srcptr);int mpfr_div_2ui(mpfr_ptr,mpfr_srcptr,unsigned long,mpfr_rnd_t);const char*mpfr_get_version();
}
namespace verified {
inline long PREC=384;inline void set_precision(long p){if(p<64)throw std::runtime_error("precision below 64"); PREC=p;}
class Ball{
public: mpfr_t lo,hi;
 Ball(){mpfr_init2(lo,PREC);mpfr_init2(hi,PREC);mpfr_set_si(lo,0,MPFR_RNDN);mpfr_set_si(hi,0,MPFR_RNDN);}
 Ball(const Ball&a):Ball(){mpfr_set(lo,a.lo,MPFR_RNDD);mpfr_set(hi,a.hi,MPFR_RNDU);}
 Ball&operator=(const Ball&a){if(this!=&a){mpfr_set(lo,a.lo,MPFR_RNDD);mpfr_set(hi,a.hi,MPFR_RNDU);}return *this;}
 ~Ball(){mpfr_clear(lo);mpfr_clear(hi);}
 void setstr(std::string s){mpq_class q(s);q.canonicalize();mpfr_set_q(lo,q.get_mpq_t(),MPFR_RNDD);mpfr_set_q(hi,q.get_mpq_t(),MPFR_RNDU);}
 void setdouble(double d){if(!std::isfinite(d))throw std::runtime_error("nonfinite");mpfr_set_d(lo,d,MPFR_RNDD);mpfr_set_d(hi,d,MPFR_RNDU);}
 static Ball integer(long d){Ball b;mpfr_set_si(b.lo,d,MPFR_RNDD);mpfr_set_si(b.hi,d,MPFR_RNDU);return b;}
 static Ball pi(){Ball b;mpfr_const_pi(b.lo,MPFR_RNDD);mpfr_const_pi(b.hi,MPFR_RNDU);return b;}
 bool iszero()const{return mpfr_zero_p(lo)&&mpfr_zero_p(hi);}
 bool finite()const{return mpfr_number_p(lo)&&mpfr_number_p(hi);}
 Ball operator+(const Ball&b)const{Ball c;mpfr_add(c.lo,lo,b.lo,MPFR_RNDD);mpfr_add(c.hi,hi,b.hi,MPFR_RNDU);return c;}
 Ball operator-(const Ball&b)const{Ball c;mpfr_sub(c.lo,lo,b.hi,MPFR_RNDD);mpfr_sub(c.hi,hi,b.lo,MPFR_RNDU);return c;}
 Ball operator-()const{Ball c;mpfr_neg(c.lo,hi,MPFR_RNDD);mpfr_neg(c.hi,lo,MPFR_RNDU);return c;}
 Ball operator*(const Ball&b)const{
  Ball c;if(iszero()||b.iszero())return c;
  if(mpfr_cmp_si(lo,0)>=0&&mpfr_cmp_si(b.lo,0)>=0){mpfr_mul(c.lo,lo,b.lo,MPFR_RNDD);mpfr_mul(c.hi,hi,b.hi,MPFR_RNDU);return c;}
  if(mpfr_cmp_si(hi,0)<=0&&mpfr_cmp_si(b.hi,0)<=0){mpfr_mul(c.lo,hi,b.hi,MPFR_RNDD);mpfr_mul(c.hi,lo,b.lo,MPFR_RNDU);return c;}
  if(mpfr_cmp_si(lo,0)>=0&&mpfr_cmp_si(b.hi,0)<=0){mpfr_mul(c.lo,hi,b.lo,MPFR_RNDD);mpfr_mul(c.hi,lo,b.hi,MPFR_RNDU);return c;}
  if(mpfr_cmp_si(hi,0)<=0&&mpfr_cmp_si(b.lo,0)>=0){mpfr_mul(c.lo,lo,b.hi,MPFR_RNDD);mpfr_mul(c.hi,hi,b.lo,MPFR_RNDU);return c;}
  Ball t;mpfr_mul(c.lo,lo,b.lo,MPFR_RNDD);mpfr_mul(c.hi,lo,b.lo,MPFR_RNDU);
  mpfr_srcptr av[2]={lo,hi};mpfr_srcptr bv[2]={b.lo,b.hi};
  for(int i=0;i<2;i++)for(int j=0;j<2;j++){mpfr_mul(t.lo,av[i],bv[j],MPFR_RNDD);mpfr_mul(t.hi,av[i],bv[j],MPFR_RNDU);if(mpfr_cmp(t.lo,c.lo)<0)mpfr_set(c.lo,t.lo,MPFR_RNDD);if(mpfr_cmp(t.hi,c.hi)>0)mpfr_set(c.hi,t.hi,MPFR_RNDU);}return c;
 }
 Ball operator/(const Ball&b)const{if(mpfr_cmp_si(b.lo,0)<=0&&mpfr_cmp_si(b.hi,0)>=0)throw std::runtime_error("interval division by zero");Ball one=integer(1),v;mpfr_div(v.lo,one.lo,b.hi,MPFR_RNDD);mpfr_div(v.hi,one.hi,b.lo,MPFR_RNDU);return (*this)*v;}
 Ball absval()const{if(mpfr_cmp_si(lo,0)>=0)return *this;if(mpfr_cmp_si(hi,0)<=0)return -(*this);Ball c;mpfr_abs(c.hi,lo,MPFR_RNDU);if(mpfr_cmp(hi,c.hi)>0)mpfr_set(c.hi,hi,MPFR_RNDU);return c;}
 Ball square()const{Ball c;if(mpfr_cmp_si(lo,0)>=0){mpfr_sqr(c.lo,lo,MPFR_RNDD);mpfr_sqr(c.hi,hi,MPFR_RNDU);}else if(mpfr_cmp_si(hi,0)<=0){mpfr_sqr(c.lo,hi,MPFR_RNDD);mpfr_sqr(c.hi,lo,MPFR_RNDU);}else{Ball a=absval();mpfr_sqr(c.hi,a.hi,MPFR_RNDU);}return c;}
 Ball power(long n)const{if(n==0)return integer(1);if(n<0)return integer(1)/power(-n);if(n%2==0)return square().power(n/2);if(n==1)return *this;return (*this)*power(n-1);}
 Ball sqrtval()const{if(mpfr_cmp_si(lo,0)<0)throw std::runtime_error("sqrt negative interval");Ball c;mpfr_sqrt(c.lo,lo,MPFR_RNDD);mpfr_sqrt(c.hi,hi,MPFR_RNDU);return c;}
 Ball expval()const{Ball c;mpfr_exp(c.lo,lo,MPFR_RNDD);mpfr_exp(c.hi,hi,MPFR_RNDU);return c;}
 Ball erfval()const{Ball c;mpfr_erf(c.lo,lo,MPFR_RNDD);mpfr_erf(c.hi,hi,MPFR_RNDU);return c;}
 Ball erfcval()const{Ball c;mpfr_erfc(c.lo,hi,MPFR_RNDD);mpfr_erfc(c.hi,lo,MPFR_RNDU);return c;}
 Ball logval()const{if(mpfr_cmp_si(lo,0)<=0)throw std::runtime_error("log nonpositive interval");Ball c;mpfr_log(c.lo,lo,MPFR_RNDD);mpfr_log(c.hi,hi,MPFR_RNDU);return c;}
 Ball midval()const{Ball c;mpfr_add(c.lo,lo,hi,MPFR_RNDN);mpfr_div_2ui(c.lo,c.lo,1,MPFR_RNDN);mpfr_set(c.hi,c.lo,MPFR_RNDN);return c;}
 Ball radval()const{Ball mid=midval(),c;mpfr_sub(c.hi,mid.hi,lo,MPFR_RNDU);Ball d;mpfr_sub(d.hi,hi,mid.lo,MPFR_RNDU);if(mpfr_cmp(d.hi,c.hi)>0)mpfr_set(c.hi,d.hi,MPFR_RNDU);mpfr_set(c.lo,c.hi,MPFR_RNDN);return c;}
 Ball lower()const{Ball c;mpfr_set(c.lo,lo,MPFR_RNDN);mpfr_set(c.hi,lo,MPFR_RNDN);return c;}
 Ball upper()const{Ball c;mpfr_set(c.lo,hi,MPFR_RNDN);mpfr_set(c.hi,hi,MPFR_RNDN);return c;}
 Ball hull(const Ball&b)const{Ball c;mpfr_set(c.lo,mpfr_cmp(lo,b.lo)<0?lo:b.lo,MPFR_RNDD);mpfr_set(c.hi,mpfr_cmp(hi,b.hi)>0?hi:b.hi,MPFR_RNDU);return c;}
 Ball cosval()const{Ball m=midval(),r=radval(),v;mpfr_cos(v.lo,m.lo,MPFR_RNDD);mpfr_cos(v.hi,m.hi,MPFR_RNDU);mpfr_sub(v.lo,v.lo,r.hi,MPFR_RNDD);mpfr_add(v.hi,v.hi,r.hi,MPFR_RNDU);return v;}
 Ball sinval()const{Ball m=midval(),r=radval(),v;mpfr_sin(v.lo,m.lo,MPFR_RNDD);mpfr_sin(v.hi,m.hi,MPFR_RNDU);mpfr_sub(v.lo,v.lo,r.hi,MPFR_RNDD);mpfr_add(v.hi,v.hi,r.hi,MPFR_RNDU);return v;}
 bool lt(const Ball&b)const{return mpfr_cmp(hi,b.lo)<0;}bool le(const Ball&b)const{return mpfr_cmp(hi,b.lo)<=0;}bool gt(const Ball&b)const{return mpfr_cmp(lo,b.hi)>0;}bool ge(const Ball&b)const{return mpfr_cmp(lo,b.hi)>=0;}bool eq(const Ball&b)const{return mpfr_cmp(lo,hi)==0&&mpfr_cmp(b.lo,b.hi)==0&&mpfr_cmp(lo,b.lo)==0;}
 bool contains(const Ball&b)const{return mpfr_cmp(lo,b.lo)<=0&&mpfr_cmp(hi,b.hi)>=0;}
 double asdouble()const{Ball m=midval();return mpfr_get_d(m.lo,MPFR_RNDN);}
 std::string mantissa()const{if(mpfr_cmp(lo,hi))throw std::runtime_error("man_exp on interval");mpz_class z;mpfr_get_z_2exp(z.get_mpz_t(),lo);return z.get_str();}
 long exponent()const{if(mpfr_zero_p(lo))return 0;mpz_class z;return mpfr_get_z_2exp(z.get_mpz_t(),lo);}
 std::string text()const{std::ostringstream ss;ss.precision(17);ss<<'['<<mpfr_get_d(lo,MPFR_RNDD)<<","<<mpfr_get_d(hi,MPFR_RNDU)<<']';return ss.str();}
};
class Matrix{public:int r,c;std::vector<Ball>a;Matrix():r(0),c(0){}Matrix(int rr,int cc):r(rr),c(cc),a(rr*cc){}Ball get(int i,int j)const{return a.at(i*c+j);}void set(int i,int j,const Ball&b){a.at(i*c+j)=b;}
 Matrix transpose()const{Matrix b(c,r);for(int i=0;i<r;i++)for(int j=0;j<c;j++)b.a[j*r+i]=a[i*c+j];return b;}
 Matrix add(const Matrix&b)const{if(r!=b.r||c!=b.c)throw std::runtime_error("matrix add shape");Matrix out(r,c);for(int k=0;k<r*c;k++)out.a[k]=a[k]+b.a[k];return out;}
 Matrix sub(const Matrix&b)const{if(r!=b.r||c!=b.c)throw std::runtime_error("matrix sub shape");Matrix out(r,c);for(int k=0;k<r*c;k++)out.a[k]=a[k]-b.a[k];return out;}
 Matrix neg()const{Matrix out(r,c);for(int k=0;k<r*c;k++)out.a[k]=-a[k];return out;}
 Matrix scale(const Ball&b)const{Matrix out(r,c);for(int k=0;k<r*c;k++)out.a[k]=a[k]*b;return out;}
 Matrix mul(const Matrix&b)const{
  if(c!=b.r)throw std::runtime_error("matrix mul shape");
  Matrix out(r,b.c);std::vector<int>sa(r*c),sb(b.r*b.c);
  for(int k=0;k<r*c;k++)sa[k]=a[k].iszero()?2:(mpfr_cmp_si(a[k].lo,0)>=0?1:(mpfr_cmp_si(a[k].hi,0)<=0?-1:0));
  for(int k=0;k<b.r*b.c;k++)sb[k]=b.a[k].iszero()?2:(mpfr_cmp_si(b.a[k].lo,0)>=0?1:(mpfr_cmp_si(b.a[k].hi,0)<=0?-1:0));
  for(int i=0;i<r;i++)for(int k=0;k<c;k++){
   const Ball&x=a[i*c+k];int xs=sa[i*c+k];if(xs==2)continue;
   for(int j=0;j<b.c;j++){
    const Ball&y=b.a[k*b.c+j];int ys=sb[k*b.c+j];if(ys==2)continue;Ball&v=out.a[i*b.c+j];
    mpfr_srcptr xl,xu,yl,yu;
    if(xs==1&&ys==1){xl=x.lo;yl=y.lo;xu=x.hi;yu=y.hi;}
    else if(xs==-1&&ys==-1){xl=x.hi;yl=y.hi;xu=x.lo;yu=y.lo;}
    else if(xs==1&&ys==-1){xl=x.hi;yl=y.lo;xu=x.lo;yu=y.hi;}
    else if(xs==-1&&ys==1){xl=x.lo;yl=y.hi;xu=x.hi;yu=y.lo;}
    else if(xs==0&&ys==1){xl=x.lo;yl=y.hi;xu=x.hi;yu=y.hi;}
    else if(xs==0&&ys==-1){xl=x.hi;yl=y.lo;xu=x.lo;yu=y.lo;}
    else if(xs==1&&ys==0){xl=x.hi;yl=y.lo;xu=x.hi;yu=y.hi;}
    else if(xs==-1&&ys==0){xl=x.lo;yl=y.hi;xu=x.lo;yu=y.lo;}
    else {v=v+x*y;continue;}
    mpfr_fma(v.lo,xl,yl,v.lo,MPFR_RNDD);mpfr_fma(v.hi,xu,yu,v.hi,MPFR_RNDU);
   }
  }return out;
 }
 static void check_dot_range(const Ball&v){
  if(!v.finite())throw std::runtime_error("nonfinite dot input");
  for(mpfr_srcptr p:{v.lo,v.hi})if(!mpfr_zero_p(p)){
   // The MPFR ABI exponent is the exponent of the normalised significand.
   if(p->_mpfr_exp < -100000 || p->_mpfr_exp > 100000)throw std::runtime_error("dot exponent guard");
  }
 }
};
class QMatrix{public:int r,c;std::vector<mpq_class>a;QMatrix():r(0),c(0){}QMatrix(int rr,int cc):r(rr),c(cc),a(rr*cc){}std::string get(int i,int j)const{return a.at(i*c+j).get_str();}void set(int i,int j,const std::string&s){a.at(i*c+j)=mpq_class(s);a.at(i*c+j).canonicalize();}
 QMatrix inverse()const{if(r!=c)throw std::runtime_error("non-square");int n=r;QMatrix tmp(n,2*n);for(int i=0;i<n;i++){for(int j=0;j<n;j++)tmp.a[i*2*n+j]=a[i*n+j];tmp.a[i*2*n+n+i]=1;}
  for(int j=0;j<n;j++){int p=j;while(p<n&&tmp.a[p*2*n+j]==0)p++;if(p==n)throw std::runtime_error("singular");if(p!=j)for(int k=0;k<2*n;k++)std::swap(tmp.a[j*2*n+k],tmp.a[p*2*n+k]);mpq_class pivot=tmp.a[j*2*n+j];for(int k=j;k<2*n;k++)tmp.a[j*2*n+k]/=pivot;for(int i=0;i<n;i++)if(i!=j){mpq_class v=tmp.a[i*2*n+j];if(v==0)continue;for(int k=j;k<2*n;k++)tmp.a[i*2*n+k]-=v*tmp.a[j*2*n+k];}}
  QMatrix out(n,n);for(int i=0;i<n;i++)for(int j=0;j<n;j++)out.a[i*n+j]=tmp.a[i*2*n+n+j];return out;}
 std::string determinant()const{if(r!=c)throw std::runtime_error("non-square");int n=r;std::vector<mpq_class>b=a;mpq_class det=1;for(int j=0;j<n;j++){int p=j;while(p<n&&b[p*n+j]==0)p++;if(p==n)return "0";if(p!=j){det=-det;for(int k=j;k<n;k++)std::swap(b[j*n+k],b[p*n+k]);}mpq_class v=b[j*n+j];det*=v;for(int i=j+1;i<n;i++){mpq_class fac=b[i*n+j]/v;for(int k=j+1;k<n;k++)b[i*n+k]-=fac*b[j*n+k];}}return det.get_str();}
};
}
