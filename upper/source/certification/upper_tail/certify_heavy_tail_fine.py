"""Full heavy-additive absolute tail from certified endpoint energies and complex arcs."""
import json,math,hashlib,sys
from pathlib import Path
from fractions import Fraction
from flint import arb,acb,ctx
assert __debug__;ctx.prec=192
HERE=Path(__file__).parent;BASE=HERE.parents[1];SOURCE=BASE/'upper/refined_explicit_mixture.json'
sys.path.insert(0,str(HERE.parent/'upper_head/arc'))
from certify_complex_arc import bound as complex_bound
raw=SOURCE.read_bytes();data=json.loads(raw);energies=json.loads((HERE/'factorized_energies_certificate.json').read_text());cap=json.loads((HERE/'early_capped_profile_certificate.json').read_text())
assert energies['status']=='RIGOROUS_FACTORIZED_CAPPED_ENERGIES_ACCEPTED'
assert cap['status']=='RIGOROUS_EARLY_C6_CAP_PERTURBATION_ACCEPTED'
assert energies['source_sha256']==cap['source_sha256']==hashlib.sha256(raw).hexdigest()
def A(x):
 q=Fraction(x);return arb(q.numerator)/q.denominator
T={tuple(map(int,k.strip('()').split(','))):A(v) for k,v in energies['proved_T_upper'].items()}
S=[[0]*5 for _ in range(5)];S[0][0]=1
for n in range(1,5):
 for j in range(1,n+1):S[n][j]=S[n-1][j-1]+j*S[n-1][j]
def convolution(a,b):return [sum((a[j]*b[n-j] for j in range(n+1)),arb(0)) for n in range(5)]
def endpoint(theta,r1,r3):
 t=acb(0,theta).exp();z=r1*t+r3*t**3;rho=abs(z).upper()
 if not rho<A('.99'):return arb.pos_inf()
 d=1-rho*rho
 G=[2/arb.pi()*math.prod(range(1,2*n,2))/(d**n*d.sqrt()) for n in range(4)]
 M={(p,0):G[p-1] for p in range(1,5)}
 for j in range(1,5):
  for p in range(5-j):M[p,j]=sum(((T[j,l]*G[p+l-1]).sqrt() for l in range(1,j+1)),arb(0))**2
 v=[arb(0)]+[abs(r1*t+r3*3**k*t**3).upper()/math.factorial(k) for k in range(1,5)]
 powers=[[arb(1)]+[arb(0)]*4]
 for p in range(4):powers.append(convolution(powers[-1],v))
 out=arb(0)
 for p in range(5):
  for k in range(5):
   coeff=A(Fraction(24,math.factorial(p)*math.factorial(k)))*powers[p][4-k]
   for j in range(k+1):out+=coeff*S[k][j]*M.get((p,j),arb(0))
 return out.upper()
weights=[Fraction(str(a['weight'])) for a in data['active']];total=sum(weights);heavy=[]
for i in [0,15,16]:
 atom=data['active'][i];assert atom['type']=='additive'
 a=Fraction(str(atom['r1']));b=Fraction(str(atom['r3']));budget=max(Fraction(1),abs(a)+abs(b))
 heavy.append((i,A(weights[i]/total),A(a/budget),A(b/budget)))
step=Fraction(1,4000);end=arb.pi()/2;cells=[];energy=arb(0);k=0
while A(k*step)<end:
 a=A(k*step);b=A((k+1)*step)
 if b>end:b=end
 theta=arb((a+b)/2,(b-a)/2);combined=arb(0);chosen=[]
 for i,w,r1,r3 in heavy:
  e=endpoint(theta,r1,r3)
  g=complex_bound(theta,r1,r3).upper() if a>=A('.6') else arb.pos_inf()
  use=e if e<g else g
  assert use.is_finite()
  combined+=w*use;chosen.append(dict(atom=i,method='endpoint' if use==e else 'complex',bound=use.str(18)))
 contribution=(b-a)*combined.upper()**2;energy+=contribution
 if k%25==0:print('cell',k,'theta',a.str(5),'bound',combined.str(8),flush=True)
 cells.append(dict(a=a.str(30),b=b.str(30),weighted_D4_bound=combined.upper().str(30),chosen=chosen))
 k+=1
norm2=2/arb.pi()*energy;norm=norm2.sqrt();perturb=sum((w for _,w,_,_ in heavy),arb(0))*arb(cap['each_additive_kernel_l1_error_ball']).upper()
tails={}
for D in [301,401,501,601,701,1001]:
 value=norm/(14*arb(D)**7).sqrt();tails[str(D)]=dict(capped_tail_bound=value.str(30),original_tail_bound=(value+perturb).str(30))
out=dict(status='RIGOROUS_HEAVY_ADDITIVE_FULL_TAIL_CERTIFIED',precision_bits=ctx.prec,source_sha256=hashlib.sha256(raw).hexdigest(),verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),energy_certificate_sha256=hashlib.sha256((HERE/'factorized_energies_certificate.json').read_bytes()).hexdigest(),cap_certificate_sha256=hashlib.sha256((HERE/'early_capped_profile_certificate.json').read_bytes()).hexdigest(),complex_arc_verifier_sha256=hashlib.sha256((HERE.parent/'upper_head/arc/certify_complex_arc.py').read_bytes()).hexdigest(),heavy_indices=[0,15,16],circle_symmetry_factor='2/pi times integral from0topi/2',cells=len(cells),weighted_D4_L2_norm_bound=norm.str(40),weighted_cap_kernel_l1_perturbation=perturb.str(30),tail_bounds=tails,scope='Only three heavy additive atoms; combine with separately certified full head, rare additive tails, and nonadditive tails.')
(HERE/'heavy_tail_fine_certificate.json').write_text(json.dumps(out,indent=2)+'\n');(HERE/'heavy_tail_fine_cells.json').write_text(json.dumps(cells,indent=2)+'\n')
print(json.dumps(out,indent=2));print('HEAVY_ADDITIVE_FULL_TAIL_ACCEPTED=1')
