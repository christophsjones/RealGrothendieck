"""Directed L1-kernel perturbation certificate for a C5 quintic tail continuation.

This defines a new, explicitly separate threshold; frozen source is read only.
"""
import json,math,hashlib
from fractions import Fraction
from pathlib import Path
from flint import arb,ctx
assert __debug__;ctx.prec=256
HERE=Path(__file__).parent;SOURCE=HERE.parents[1]/'upper/refined_explicit_mixture.json'
raw=SOURCE.read_bytes();data=json.loads(raw)
def A(x):
 q=Fraction(x);return arb(q.numerator)/q.denominator
coeff=[A(x) for x in data['threshold_definition']['r_coefficients_exact_decimal']]
def hermites(x,n):
 h=[arb(1),x]
 for k in range(1,n):h.append((x*h[-1]-arb(k).sqrt()*h[-2])/arb(k+1).sqrt())
 return h[:n+1]
def taylor(center):
 h=hermites(A(center),15)
 return [sum((c*arb(math.factorial(d)//math.factorial(d-j)).sqrt()*h[d-j] for c,d in zip(coeff,range(1,16,2)) if d>=j),arb(0))/math.factorial(j) for j in range(16)]
t7=taylor(7)
assert all(c>0 for c in t7)
assert t7[0]>A('13.98')
# |r(x)|<=3.4 on [-4,4], using exact panel centers and Taylor remainders (polynomial is finite).
for k in range(128):
 center=Fraction(-4)+Fraction(2*k+1,32)
 # 128 panels each of width 1/16; halfwidth 1/32.
 tc=taylor(center)
 enclosure=tc[0]+sum((abs(tc[j])*A(Fraction(1,32))**j for j in range(1,16)),arb(0))*arb(0,1)
 assert abs(enclosure)<A('3.4'),(k,enclosure)
sqrt2=arb(2).sqrt()
pa=(A(7)/sqrt2).erfc();pb=(A(4)/sqrt2).erfc();pw=(A('7.18')/sqrt2).erfc()
disagreement=6*pa*pb+3*pa*pw
kernel_error=4*disagreement.sqrt()
assert kernel_error<A('0.000000125')
additive_weight=sum(Fraction(str(a['weight'])) for a in data['active'] if a['type']=='additive')/sum(Fraction(str(a['weight'])) for a in data['active'])
mixture_error=A(additive_weight)*kernel_error
out=dict(status='RIGOROUS_QUINTIC_CONTINUATION_PERTURBATION_ACCEPTED',precision_bits=ctx.prec,source_sha256=hashlib.sha256(raw).hexdigest(),verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),definition='Keep r(x) for |x|<=7. For x>7 use sum_{j=0}^5 r^(j)(7)(x-7)^j/j!. Extend oddly to x<-7.',smoothness='C5 on R, piecewise polynomial degree15 inside and degree5 outside',taylor_at7=[str(c) for c in t7],certified_inside_bound='|r(x)|<17/5 for |x|<=4',certified_outside_bound='r(x)>=699/50 and r_tilde(x)>=699/50 for x>=7',disagreement_probability_upper_ball=str(disagreement),each_additive_kernel_l1_error_ball=str(kernel_error),proved_each_kernel_error_upper='1/8000000',normalized_additive_weight=str(additive_weight),whole_mixture_l1_error_ball=str(mixture_error),scope='Perturbation only; no full-tail or complete rounding-rate certificate is asserted.')
(HERE/'quintic_continuation_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2));print('QUINTIC_CONTINUATION_PERTURBATION_ACCEPTED=1')
# A bounded C6 cap has the same disagreement event and the same certificate.
chi=[1,0,0,0,0,0,-462,1980,-3465,3080,-1386,252]
assert sum(chi)==0
assert [k*chi[k] for k in range(1,12)]==[0,0,0,0,0]+[-2772*(-1)**j*math.comb(5,j) for j in range(6)]
assert all(sum(math.factorial(k)//math.factorial(k-j)*chi[k] for k in range(j,12))==0 for j in range(6))
cap=[t7[0]]+[arb(0)]*26
for j in range(1,16):
 for k,v in enumerate(chi):cap[j+k]+=j*t7[j]*v/(j+k)
plateau=sum(cap,arb(0));assert plateau>A('13.98')
capout=dict(out,status='RIGOROUS_BOUNDED_C6_CAP_PERTURBATION_ACCEPTED',definition='Keep r on |x|<=7. On x=7+s with 0<=s<=1 set rcap=r(7)+integral_0^s rprime(7+u)chi(u)du. For x>=8 keep rcap(8) constant. Extend oddly.',smoothness='C6 on R; bounded, with compactly supported derivative',chi_coefficients=chi,cap_power_coefficients=[str(v) for v in cap],plateau_ball=str(plateau))
(HERE/'capped_profile_certificate.json').write_text(json.dumps(capout,indent=2)+'\n')
print('CAPPED_PROFILE_PERTURBATION_ACCEPTED=1',plateau)
# An earlier, shorter cap controls rare high derivatives at still-small full-kernel cost.
early=Fraction(27,4);width=Fraction(1,4);te=taylor(early)
assert all(c>0 for c in te);assert te[0]>A(7)
for k in range(112):
 center=Fraction(-7,2)+Fraction(2*k+1,32);tc=taylor(center)
 enclosure=tc[0]+sum((abs(tc[j])*A(Fraction(1,32))**j for j in range(1,16)),arb(0))*arb(0,1)
 assert abs(enclosure)<A('1.65')
pe=(A(early)/sqrt2).erfc();pbe=(A('3.5')/sqrt2).erfc();pwe=(A('3.7')/sqrt2).erfc()
de=6*pe*pbe+3*pe*pwe;ee=4*de.sqrt();assert ee<A('0.000001')
ce=[te[0]]+[arb(0)]*26
for j in range(1,16):
 for k,v in enumerate(chi):ce[j+k]+=j*te[j]*v/A(width)**k/(j+k)
plateaue=sum((c*A(width)**j for j,c in enumerate(ce)),arb(0))
eo=dict(capout,status='RIGOROUS_EARLY_C6_CAP_PERTURBATION_ACCEPTED',definition='Keep r on |x|<=27/4. On x=27/4+s for 0<=s<=1/4 use r(27/4)+integral_0^s rprime(27/4+u)chi(4u)du. Keep constant for x>=7 and extend oddly.',cap_start=str(early),cap_width=str(width),cap_power_coefficients=[str(c) for c in ce],plateau_ball=str(plateaue),certified_inside_bound='|r(x)|<33/20 for |x|<=7/2',certified_outside_bound='Both profiles>=7 for x>=27/4',disagreement_probability_upper_ball=str(de),each_additive_kernel_l1_error_ball=str(ee),proved_each_kernel_error_upper='1/1000000',whole_mixture_l1_error_ball=str(A(additive_weight)*ee))
(HERE/'early_capped_profile_certificate.json').write_text(json.dumps(eo,indent=2)+'\n')
print('EARLY_CAPPED_PROFILE_PERTURBATION_ACCEPTED=1',ee,plateaue)
