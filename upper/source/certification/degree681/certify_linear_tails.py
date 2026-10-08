"""Directed full absolute-tail certificate for the frozen Gaussian-linear atoms.

Run with python-flint, without -O. No quadrature is used. The additive and
weighted-chaos atoms are explicitly outside this component certificate.
"""
import hashlib,json,math
from pathlib import Path
from fractions import Fraction
from flint import arb,arb_series,ctx
assert __debug__
ctx.prec=256;ctx.cap=683
HERE=Path(__file__).parent
SOURCE=HERE.parents[1]/'upper/refined_explicit_mixture.json'
raw=SOURCE.read_bytes();data=json.loads(raw);D=681
def A(x):
 q=Fraction(x);return arb(q.numerator)/q.denominator
weights=[Fraction(str(a['weight'])) for a in data['active']];total=sum(weights)
result=[];tail=arb(0)
for i,(atom,w) in enumerate(zip(data['active'],weights)):
 if atom['type']!='linear':continue
 coeff={int(k):Fraction(str(v)) for k,v in atom['powers'].items()}
 budget=max(Fraction(1),sum(abs(v) for v in coeff.values()));coeff={k:v/budget for k,v in coeff.items()}
 if set(coeff)=={1,3} and coeff[1]>Fraction(4,5) and coeff[3]<0:
  # Prove |Re P(z)|<=eta<1 on |z|<=R by its cubic Chebyshev boundary polynomial.
  R=Fraction(6,5) if coeff[1]<Fraction(9,10) else Fraction(21,20)
  eta=Fraction(19,20) if R==Fraction(6,5) else Fraction(49,50)
  aa=coeff[1]*R;bb=-coeff[3]*R**3;C=aa+3*bb
  assert abs(aa-bb)<eta
  if C<12*bb:assert C**3<27*bb*eta**2
  # Otherwise the unique positive interior critical point is outside [0,1].
  M=(2/arb.pi())*A(aa+bb)/(1-A(eta)**2).sqrt()
  bound=M*A(R)**(-D-2)/(1-A(R)**-2)
  method='COMPLEX_DISK_REAL_STRIP_CAUCHY'
  details=dict(radius=str(R),real_part_bound=str(eta),boundary_modulus_bound=str(M))
 else:
  # Nonnegative coefficient majorant M(t)=(2/pi)asin(sum |p_d|t^d).
  # Its total mass is <=1, and every absolute coefficient of the atom is bounded by M.
  p=[arb(0)]*(D+1)
  for d,v in coeff.items():p[d]=A(abs(v))
  series=arb_series(p,D+1).asin()*(2/arb.pi())
  head=sum((series[m] for m in range(1,D+1,2)),arb(0))
  bound=1-head
  assert bound>0
  method='POSITIVE_ARCSINE_COEFFICIENT_MAJORANT'
  details=dict(majorant_head_mass=str(head))
 contribution=A(w/total)*bound;tail+=contribution
 result.append(dict(atom_index=i,normalized_weight=str(w/total),method=method,details=details,unweighted_tail_ball=str(bound),weighted_tail_ball=str(contribution)))
assert tail<A('0.000004')
out=dict(status='RIGOROUS_LINEAR_ATOM_TAIL_COMPONENT_ACCEPTED',cutoff=D,precision_bits=ctx.prec,source_sha256=hashlib.sha256(raw).hexdigest(),verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),linear_atom_count=len(result),weighted_linear_tail_ball=str(tail),proved_weighted_linear_tail_upper='1/250000',atoms=result,scope='Only Gaussian-linear atoms; additive and weighted(3,5) tails are not covered.')
(HERE/'linear_tail_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2));print('LINEAR_ATOM_FULL_TAIL_ACCEPTED=1')
