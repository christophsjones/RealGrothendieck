"""Re-evaluate every accumulated plane error and bind the output hashes."""
import json,hashlib
from pathlib import Path
import numpy as np
from flint import arb,ctx
ctx.prec=256
HERE=Path(__file__).resolve().parent
meta=json.loads((HERE/'certified_b_inputs.json').read_text());old=json.loads((HERE/'certified_plane.json').read_text());u=np.load(HERE/'certified_b_inputs.npz')['U'];C=np.load(HERE/'certified_plane.npy')
S=meta['nodes'];N=2048;unit=arb(2)**-53;M=arb(meta['B_norm_bound'])**2
eta=1100*unit/(1-1100*unit);eta20=20*unit/(1-20*unit)
d=arb(N).sqrt()*(eta*(1+unit)+unit)
derived=arb(N).sqrt()*(3*(1+d)**2*(eta*(1+unit)+unit)+eta20*(1+d)**3+eta*(1+eta20)*(1+d)**3)*M**3
allocated=10000*N*unit*M**3+arb('1e-290')
assert derived<allocated
factor=arb(float(np.load(HERE/'certified_b_inputs.npz')['factor']))
usq=arb(meta['U_l1_squares_sum']);usum=arb(meta['U_l1_sum'])
fft_error=factor*usq*allocated
drop=sum(abs(arb(float(x))) for x in u[u<1e-18]);prune=factor*(2*usum*drop+drop*drop)*M**3
k=2*S*S+8*S;gamma=k*unit/(1-k*unit)
contract=gamma*factor*usq*(M**3+allocated)
input_error=arb(meta['plane_l1_U_input_error'])*(M**3+allocated)
total=arb(meta['plane_l1_quadrature_error'])+input_error+fft_error+prune+contract
assert total<arb('3e-7')
old.update({'status':'DIRECTED_FINITE_NATIVE_COEFFICIENT_PLANE','precision_bits_error_budget':ctx.prec,'plane_l1_error':total.str(50),'fft_plane_l1_error':fft_error.str(50),'contraction_l1_error':contract.str(50),'pruning_l1_error':prune.str(50),'FFT_derived_per_pair_l1_error':derived.str(50),'FFT_allocated_per_pair_l1_error':allocated.str(50),'FFT_allocation_passed':True,'sum_computed_C':float(C.sum()),'plane_sha256':hashlib.sha256((HERE/'certified_plane.npy').read_bytes()).hexdigest(),'input_array_sha256':hashlib.sha256((HERE/'certified_b_inputs.npz').read_bytes()).hexdigest(),'input_metadata_sha256':hashlib.sha256((HERE/'certified_b_inputs.json').read_bytes()).hexdigest(),'error_verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
(HERE/'certified_plane.json').write_text(json.dumps(old,indent=2));print(json.dumps(old,indent=2))
