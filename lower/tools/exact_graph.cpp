// Exact fixed-point construction and integer maximum-flow certificates.
// No floating-point operation occurs in this file.
#include <vector>
#include <queue>
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <cstring>
#include <limits>
using I=__int128_t;using U=__uint128_t;using L=int64_t;
static const I LIM=(I(1)<<126)-1;
static I ab(I x){return x<0?-x:x;}
static L narrow(I x){if(x<-(I(1)<<62)||x>(I(1)<<62))throw std::runtime_error("integer range");return (L)x;}
static std::string text(I x){bool neg=x<0;U y=neg?U(-x):U(x);std::string z;do{z.push_back('0'+y%10);y/=10;}while(y);if(neg)z.push_back('-');std::reverse(z.begin(),z.end());return z;}
static void message(char* buf,int cap,std::string s){if(cap>0){std::strncpy(buf,s.c_str(),cap-1);buf[cap-1]=0;}}
static L rndshift(I x,int shift){if(shift<=0||shift>124)throw std::runtime_error("shift");I z=(ab(x)+(I(1)<<(shift-1)))>>shift;return narrow(x<0?-z:z);}
static L rndratio(I x,I d){if(d<=0)throw std::runtime_error("denominator");I v=ab(x),q=v/d,r=v%d;if(2*r>=d)q++;return narrow(x<0?-q:q);}
extern "C" int build_graph(int N,const L* H,const L* M,const L* ell,L* B,L* a,L* constant,char* qerror,int qcap,char* error,int ecap){
 try{
  if(N<2||N>10000)throw std::runtime_error("node count");int n=2*N;
  I hmax=0,mrow=0;
  for(int j=0;j<22;j++)for(int k=0;k<N;k++)hmax=std::max(hmax,ab(H[j*N+k]));
  for(int j=0;j<44;j++){I s=0;for(int k=0;k<44;k++){if(M[j*44+k]!=M[k*44+j])throw std::runtime_error("asymmetric M");s+=ab(M[j*44+k]);}mrow=std::max(mrow,s);}
  if(hmax>(I(1)<<48)||mrow>(I(1)<<52))throw std::runtime_error("input too large");
  I ubound=mrow*hmax;if(hmax && ubound>LIM/(88*hmax))throw std::runtime_error("dot product overflow bound");
  std::vector<I> C(44*n);
  for(int r=0;r<44;r++)for(int j=0;j<n;j++){
   int ch=j/N,k=j%N;I v=0;for(int t=0;t<22;t++)v+=I(M[r*44+ch*22+t])*H[t*N+k];C[r*n+j]=v;
  }
  I qe=0,diag=0;
  for(int i=0;i<n;i++){
   int ch=i/N,k=i%N;I lin=0;for(int t=0;t<22;t++)lin+=I(ell[ch*22+t])*H[t*N+k];
   a[i]=rndshift(lin,38);qe+=ab(lin-I(a[i])*(I(1)<<38))*(I(1)<<42);
   for(int j=i;j<n;j++){
    I val=0;for(int t=0;t<22;t++)val+=I(H[t*N+k])*C[(ch*22+t)*n+j];
    if(i==j){if(__builtin_add_overflow(diag,val,&diag)||diag>LIM||diag<-LIM)throw std::runtime_error("diagonal overflow");B[i*n+j]=0;}
    else{L z=rndshift(val,79);B[i*n+j]=z;B[j*n+i]=z;qe+=ab(2*val-I(z)*(I(1)<<80));}
   }
  }
  *constant=rndshift(diag,80);qe+=ab(diag-I(*constant)*(I(1)<<80));message(qerror,qcap,text(qe));return 0;
 }catch(const std::exception& e){message(error,ecap,e.what());return 1;}
}
struct Edge{int to,rev;L cap,initial;};
class Flow {
 public:int n;std::vector<std::vector<Edge>>g;std::vector<int>level,it;
 Flow(int n,int reserve):n(n),g(n),level(n),it(n){for(int i=0;i<n-2;i++)g[i].reserve(reserve);g[n-2].reserve(n);g[n-1].reserve(n);}
 void add(int u,int v,L cap,bool undirected=false){if(cap<0)throw std::runtime_error("negative capacity");if(!cap)return;int a=g[u].size(),b=g[v].size();if(u==v)throw std::runtime_error("self arc");g[u].push_back({v,b,cap,cap});g[v].push_back({u,a,undirected?cap:0,undirected?cap:0});}
 bool bfs(int S,int T){std::fill(level.begin(),level.end(),-1);std::queue<int>q;q.push(S);level[S]=0;while(!q.empty()){int v=q.front();q.pop();for(const auto&e:g[v])if(e.cap>0&&level[e.to]<0){level[e.to]=level[v]+1;q.push(e.to);}}return level[T]>=0;}
 L dfs(int v,int T,L f){if(v==T)return f;for(int&j=it[v];j<(int)g[v].size();j++){auto&e=g[v][j];if(e.cap>0&&level[e.to]==level[v]+1){L z=dfs(e.to,T,std::min(f,e.cap));if(z>0){e.cap-=z;g[e.to][e.rev].cap+=z;return z;}}}return 0;}
 L solve(int S,int T){L total=0,z;while(bfs(S,T)){std::fill(it.begin(),it.end(),0);while((z=dfs(S,T,L(1)<<61))>0)total+=z;}check(S,T,total);return total;}
 void check(int S,int T,L value){
  if(level[S]<0||level[T]>=0)throw std::runtime_error("bad cut");I cut=0;
  for(int u=0;u<n;u++){
   I balance=0;
   for(int j=0;j<(int)g[u].size();j++){
    const auto&e=g[u][j];if(e.to<0||e.to>=n||e.rev<0||e.rev>=(int)g[e.to].size())throw std::runtime_error("bad reverse link");
    const auto&r=g[e.to][e.rev];if(r.to!=u||r.rev!=j||e.cap<0||I(e.cap)+r.cap!=I(e.initial)+r.initial)throw std::runtime_error("infeasible flow");
    balance+=I(e.initial)-e.cap;
    if(level[u]>=0&&level[e.to]<0)cut+=e.initial;
   }
   I expected=u==S?I(value):(u==T?-I(value):I(0));if(balance!=expected)throw std::runtime_error("flow conservation");
  }
  if(cut!=value)throw std::runtime_error("flow cut inequality not tight");
 }
};
extern "C" int exact_roof(int n,const L* B,const L* a,L constant,L* bound,L* flowvalue,int* labels,char* error,int ecap){
 try{
  if(n<0||n>20000)throw std::runtime_error("graph size");I base=constant,absolute=ab(constant);
  for(int i=0;i<n;i++){absolute+=ab(a[i]);base+=ab(a[i]);if(B[i*n+i])throw std::runtime_error("nonzero graph diagonal");for(int j=i+1;j<n;j++){if(B[i*n+j]!=B[j*n+i])throw std::runtime_error("asymmetric graph");base+=ab(B[i*n+j]);absolute+=ab(B[i*n+j]);}}
  if(absolute>=(I(1)<<60))throw std::runtime_error("flow overflow budget");
  int S=2*n,T=S+1;Flow f(2*n+2,n+2);
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
   L z=B[i*n+j];if(z>0){f.add(i,j,z,true);f.add(n+i,n+j,z,true);}else if(z<0){f.add(i,n+j,-z,true);f.add(n+i,j,-z,true);}
  }
  for(int i=0;i<n;i++){
   if(a[i]>0){f.add(S,i,a[i]);f.add(n+i,T,a[i]);}else if(a[i]<0){f.add(i,T,-a[i]);f.add(S,n+i,-a[i]);}
  }
  *flowvalue=f.solve(S,T);*bound=narrow(base-*flowvalue);
  for(int i=0;i<n;i++)labels[i]=(f.level[i]>=0&&f.level[n+i]<0)?1:((f.level[i]<0&&f.level[n+i]>=0)?-1:0);
  return 0;
 }catch(const std::exception&e){message(error,ecap,e.what());return 1;}
}
// Adds a rounded quadratic cap mu*(u-lo)*(hi-u), u=sum(v_i*s_i)/V.
// Inputs: mu=p/2^20; lo=l/2^16; hi=h/2^16.  Output error is in graph units.
extern "C" int cap_graph(int n,const L* B,const L* a,L constant,const L*v,L p,L l,L h,L* CB,L* ca,L* cc,L* error_units,char* error,int ecap){
 try{
  if(p<0||p>(L(1)<<30)||l<-65536||h>65536||l>h)throw std::runtime_error("cap parameters");
  I V=0,squares=0,vmax=0;for(int i=0;i<n;i++){V+=ab(v[i]);squares+=I(v[i])*v[i];vmax=std::max(vmax,ab(v[i]));}
  if(V==0||V>(I(1)<<30)||vmax>(I(1)<<24))throw std::runtime_error("cap coordinate");
  I den1=(I(1)<<20)*V*V,den2=(I(1)<<36)*V,den3=I(1)<<52;
  I scale=I(1)<<44;L nonexact=0;
  for(int i=0;i<n;i++){
   I num=I(p)*(l+h)*v[i]*scale;L z=rndratio(num,den2);ca[i]=narrow(I(a[i])+z);if(num%den2)nonexact++;
   CB[i*n+i]=0;
   for(int j=i+1;j<n;j++){
    I bnum=-2*I(p)*v[i]*v[j]*scale;L zz=rndratio(bnum,den1);L newval=narrow(I(B[i*n+j])+zz);CB[i*n+j]=newval;CB[j*n+i]=newval;if(bnum%den1)nonexact++;
   }
  }
  if(I(p)*squares>LIM/scale)throw std::runtime_error("cap constant overflow bound");
  I num1=-I(p)*squares*scale,num2=-I(p)*l*h*scale;
  L c1=rndratio(num1,den1),c2=rndratio(num2,den3);if(num1%den1)nonexact++;if(num2%den3)nonexact++;
  *cc=narrow(I(constant)+c1+c2);*error_units=(nonexact+1)/2;return 0;
 }catch(const std::exception&e){message(error,ecap,e.what());return 1;}
}

extern "C" int validate_graph(int n,const L* B,const L* a,L constant,L* budget,char* error,int ecap){
 try {
  if(n<0||n>20000)throw std::runtime_error("graph size");
  I total=ab(constant);
  for(int i=0;i<n;i++){
   total+=ab(a[i]);if(B[i*n+i])throw std::runtime_error("nonzero graph diagonal");
   for(int j=i+1;j<n;j++){
    if(B[i*n+j]!=B[j*n+i])throw std::runtime_error("asymmetric graph");
    total+=ab(B[i*n+j]);
   }
  }
  if(total>=(I(1)<<60))throw std::runtime_error("graph integer overflow budget");
  *budget=(L)total;return 0;
 }catch(const std::exception&e){message(error,ecap,e.what());return 1;}
}
