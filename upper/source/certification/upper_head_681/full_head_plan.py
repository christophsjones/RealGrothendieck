"""Directed Gauss quadrature plan for a simultaneous Hermite-head enclosure.

All quadrature error bounds are in Hilbert norm, not per coefficient.
"""
import json,math,hashlib
from pathlib import Path
from flint import arb,arb_poly,ctx
ctx.prec=512
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parents[1]/'upper/refined_explicit_mixture.json'
raw=json.loads(SOURCE.read_text(),parse_float=str)
cs=raw['threshold_definition']['r_coefficients_exact_decimal']
hs=[arb_poly([1]),arb_poly([0,1])]
for k in range(1,15):hs.append((arb_poly([0,1])*hs[-1]-arb(k).sqrt()*hs[-2])*(1/arb(k+1).sqrt()))
rp=sum((arb(c)*hs[d] for c,d in zip(cs,range(1,16,2))),arb_poly([]))
r=rp.coeffs();D=681;L=arb('7.3');T=arb(32);ninner=32;nouter=24
rho=arb('0.99');eps_inner=arb('1e-10');sqrt2pi=(2*arb.pi()).sqrt()

def shifted(c):
 return [sum((r[j]*math.comb(j,k)*c**(j-k) for j in range(k,len(r))),arb(0)) for k in range(len(r))]

def inner_error(a,b):
 c=(a+b)/2;h=(b-a)/2;shift=shifted(c);radius=arb('1.25')*h;y=arb('0.75')*h
 dr=sum((k*abs(shift[k])*radius**(k-1) for k in range(1,len(shift))),arb(0))
 imag_r=y*dr;xmin=max(arb(0),a-h/4)
 mv=rho**(-arb(D)/2)*(1-rho*rho)**arb('-.25')/sqrt2pi
 mv*=(-xmin*xmin/(2*(1+rho))+y*y/(2*(1-rho))+T*imag_r).exp()
 # Bernstein ellipse parameter two. Positivity of GL weights and
 # Chebyshev truncation give error <=8*length*M*2^(-2n).
 return 8*(b-a)*mv*arb(2)**(-2*ninner)

pending=[(arb(0),L)];cells=[];total=arb(0)
while pending:
 a,b=pending.pop();error=inner_error(a,b)
 if error < eps_inner*(b-a)/(2*L):
  cells.append((a,b,error));total+=2*error
 else:
  c=(a+b)/2;pending.append((c,b));pending.append((a,c))
  if len(pending)+len(cells)>100000:raise RuntimeError('quadrature plan exceeded cell budget')
cells.sort(key=lambda v:float(v[0].mid()))

# Average tilted moments give a much sharper outer analytic bound than max|r|.
# Every outer GL panel has length 1, so its parameter-two ellipse has
# maximum imaginary part 3/8.
y=arb(3)/8;j0=arb(0);j2=arb(0)
for k in range(7300):
 a=arb(k)/1000;b=arb(k+1)/1000;c=(a+b)/2;h=(b-a)/2
 v=rp(arb(c,h));upper=abs(v).upper();wt=(b-a)*(-a*a/2).exp()/sqrt2pi
 tilt=(2*y*upper).exp();j0+=2*wt*tilt;j2+=2*wt*upper*upper*tilt
outer_M=2*arb(2).sqrt()/arb.pi()*(y*y).exp()*((1+4*y*y).sqrt()*j0**arb('1.5')+3*j2.sqrt()*j0)
outer_error=8*T*outer_M*arb(2)**(-2*nouter)

# Fourier tail for projection onto W chaos a<=D: Minkowski and incomplete
# gamma integrals. |B_j| vector has Hilbert norm <=1.
tail=arb(0)
for a in range(D+1):
 if a==0:
  term=(-T*T/2).exp()/(T*T)
 else:
  term=arb(2)**(arb(a)/2-1)*(T*T/2).gamma_upper(arb(a)/2)/arb.fac_ui(a).sqrt()
 tail+=term
tail*=2/arb.pi()
cutoff_error=2*(3*(L/arb(2).sqrt()).erfc()).sqrt()

out={'status':'DIRECTED_QUADRATURE_PLAN_ONLY_NO_COEFFICIENTS_YET','candidate_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'precision_bits':ctx.prec,'D':D,'L':'7.3','T':'32','inner_gauss_order':ninner,'outer_gauss_order':nouter,'inner_cells':[[a.str(170),b.str(170)] for a,b,e in cells],'inner_nodes_positive_half':len(cells)*ninner,'outer_nodes':32*nouter,'inner_vector_quadrature_error':total.str(40),'outer_Hilbert_quadrature_error':outer_error.str(40),'outer_analytic_bound':outer_M.str(40),'tilted_J0_bound':j0.str(40),'tilted_J2_bound':j2.str(40),'W_projection_Fourier_tail_Hilbert':tail.str(40),'whole_kernel_cutoff_error':cutoff_error.str(40)}
(HERE/'full_head_quadrature_plan.json').write_text(json.dumps(out,indent=2))
print({k:v for k,v in out.items() if k!='inner_cells'})
