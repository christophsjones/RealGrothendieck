"""Directed non-fibre certificate for the frozen degree-21 lower candidate.

No universal fibre assertion is made. All finite input decimals are exact
rationals. Every PSD/trigonometric check uses Arb balls and strict inequalities.
"""
from pathlib import Path
from decimal import Decimal
from fractions import Fraction as F
import json, hashlib, argparse, time, sys
import flint
if not __debug__:
    raise RuntimeError("Run without -O: certificate checks must not be disabled")
from flint import arb,ctx

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SOURCE=ROOT/'lower_push_20260906/multiple_caps/best_candidate/frozen_rows.json'
QSOURCE=ROOT/'lower_push_20260906/multiple_caps/best_candidate/frozen_q.json'
ap=argparse.ArgumentParser();ap.add_argument('--precision',type=int,default=384);ap.add_argument('--epsilon',default='0.0000001');args=ap.parse_args()
ctx.prec=args.precision
assert ctx.prec>=128
raw=json.loads(SOURCE.read_text(),parse_float=Decimal)
qraw=json.loads(QSOURCE.read_text(),parse_float=Decimal)
D=21;n=11;even=list(range(0,22,2));odd=list(range(1,22,2))
q=[F(x) for x in qraw['coefficients_exact_decimal']]
assert len(q)==11 and q[0]==1
EPS=F(args.epsilon);assert EPS>0
TARGET=F(1000,1773)-F(1,10**9)

def a(x):
    x=F(x)
    return arb(x.numerator)/arb(x.denominator)

def fracstr(x):
    x=F(x)
    return str(x.numerator) if x.denominator==1 else f'{x.numerator}/{x.denominator}'

def exact_lower(x):
    m,e=x.lower().man_exp();return F(int(m))*(F(2)**int(e))

def exact_upper(x):
    m,e=x.upper().man_exp();return F(int(m))*(F(2)**int(e))

def floor_decimal(x,digits=12):
    scale=10**digits;return F((x.numerator*scale)//x.denominator,scale)

def decimalstr(x,digits=12):
    assert (x*10**digits).denominator==1
    val=(x*10**digits).numerator;sign='-' if val<0 else '';val=abs(val)
    return sign+str(val//10**digits)+'.'+str(val%10**digits).zfill(digits)

def rows_for_circle(N):
    out=set();k=(N-1)//2
    for code in range(1<<(N-1)):
        s=[1]+[1-2*((code>>(j-1))&1) for j in range(1,N)]
        out.add(tuple(sum(s[j]*(s[j+d] if j+d<N else -s[j+d-N]) for j in range(N)) for d in range(1,k+1)))
    return sorted(out)

CUTS={N:rows_for_circle(N) for N in [7,9]}
COS={N:[(arb.pi()*j/N).cos() for j in range(1,(N+1)//2)] for N in CUTS}
for N,cc in COS.items():
    assert all(x>0 and x<1 for x in cc)
    assert all(cc[j]>cc[j+1] for j in range(len(cc)-1))

checks=0
min_pivot=None

def prove_pd(M,label):
    """Unpivoted interval LDL^T: strict positive pivots prove positive definite."""
    global checks,min_pivot
    m=len(M);L=[[arb(0) for _ in range(m)] for _ in range(m)];dd=[]
    for j in range(m):
        d=M[j][j]-sum((L[j][k]**2*dd[k] for k in range(j)),arb(0))
        if not d>0:raise ArithmeticError(f'PD check failed: {label}, pivot {j}: {d}')
        dd.append(d);L[j][j]=arb(1)
        lo=exact_lower(d)
        if min_pivot is None or lo<min_pivot:min_pivot=lo
        for i in range(j+1,m):
            L[i][j]=(M[i][j]-sum((L[i][k]*L[j][k]*dd[k] for k in range(j)),arb(0)))/d
    checks+=1
    return min(exact_lower(v) for v in dd)

nu=(arb(2)/arb.pi()).sqrt()
assert len(raw['rows'])==35
intervals=[[F(x) for x in row['interval']] for row in raw['rows']]
assert intervals[0][0]==0 and intervals[-1][1]==1
assert all(lo<hi for lo,hi in intervals)
assert all(intervals[i][1]==intervals[i+1][0] for i in range(len(intervals)-1))
start=time.time();result=[];rational_rows=[]
for idx,row in enumerate(raw['rows']):
    gamma={m:arb(0) for m in odd};cap=F(0);cap_records=[]
    for cp in row['caps']:
        N=cp['N'];assert N in CUTS and not cp['joint'] and not cp['rho']
        ww=[F(x) for x in cp['weights']];cc=COS[N]
        assert len(ww)==len(cc) and ww[0]>0
        cb=max(sum((w*z for w,z in zip(ww,cut)),F(0))/N for cut in CUTS[N])
        cap+=cb
        vals=[]
        for m in range(D+2,102,2):
            gm=sum((a(w)*c**m for w,c in zip(ww,cc)),arb(0))
            assert gm>0,(idx,N,m,gm)
            vals.append(exact_lower(gm))
        rem=a(ww[0])-sum((a(max(F(0),-ww[j]))*(cc[j]/cc[0])**101 for j in range(1,len(cc))),arb(0))
        assert rem>0,(idx,N,rem)
        for m in odd:gamma[m]+=sum((a(w)*c**m for w,c in zip(ww,cc)),arb(0))
        cap_records.append(dict(N=N,weights=[fracstr(w) for w in ww],exact_cut_bound=fracstr(cb),cut_row_count=len(CUTS[N]),all_finite_omitted_coefficients_positive=True,finite_omitted_min_lower=str(min(vals)),tail_cone_margin_lower=str(exact_lower(rem))))
    MM=[[[F(x) for x in line] for line in M] for M in row['matrices']]
    assert all(M[i][j]==M[j][i] for M in MM for i in range(n) for j in range(n))
    ellh=[F(x) for x in row['lh']];ellk=[F(x) for x in row['lk']]
    ch=F(row['ch']);ck=F(row['ck']);tt=F(row['t']);CC=F(row['C'])
    modified=[[[a(M[i][j]+(EPS if i==j else 0)) for j in range(n)] for i in range(n)] for M in MM]
    hs={m:arb(0) if m==1 else a(q[(m-1)//2])-gamma[m] for m in odd}
    ks={m:-a(q[(m-1)//2])-gamma[m] for m in odd}
    row_pivots=[]
    for k,M in enumerate(modified):
        row_pivots.append(prove_pd(M,f'row{idx} M{k}'))
        row_pivots.append(prove_pd([[(arb(10) if i==j else arb(0))-M[i][j] for j in range(n)] for i in range(n)],f'row{idx} 10I-M{k}'))
    for k,ix,rs,co in [(0,even,range(1,22,2),hs),(1,odd,range(2,22,2),hs),(2,even,range(1,22,2),ks),(3,odd,range(2,22,2),ks)]:
        for r in rs:
            S=[[modified[k][i][j]-(co.get(ix[i]+r,arb(0)) if i==j else arb(0)) for j in range(n)] for i in range(n)]
            row_pivots.append(prove_pd(S,f'row{idx} M{k} transverse{r}'))
    for k,co,ell,cost in [(1,hs,ellh,ch),(3,ks,ellk,ck)]:
        S=[[modified[k][i][j]-(co[odd[i]] if i==j else arb(0)) for j in range(n)]+[a(ell[i])] for i in range(n)]
        S.append([a(x) for x in ell]+[a(cost+EPS)])
        row_pivots.append(prove_pd(S,f'row{idx} mean M{k}'))
    lo,hi=intervals[idx]
    A=(1-gamma[1])*nu**2;B=-a(tt)*nu
    assert A>0,(idx,A)  # all frozen rows are convex: endpoints suffice exactly
    pvals=[A*a(x)**2+B*a(x) for x in [lo,hi]]
    scalar_upper=max(exact_upper(v) for v in pvals)
    # Parent proves ORIGINAL fibre objective <= C_allowed. Joint Parseval then
    # raises its modified-matrix counterpart by at most EPS.
    other_upper=scalar_upper+cap+ch+ck+3*EPS
    allowed=floor_decimal(TARGET-other_upper)
    assert allowed+other_upper<=TARGET
    ref_excess=max(F(0),F(raw['history'][-1]['violations'][idx]))
    slack=allowed-CC-ref_excess
    result.append(dict(index=idx,interval=[fracstr(lo),fracstr(hi)],epsilon=fracstr(EPS),original_fibre_C_allowed=decimalstr(allowed),original_stored_C=fracstr(CC),located_excess_reference=fracstr(ref_excess),allowance_above_stored_C_and_located_excess=fracstr(slack),non_fibre_cost_upper=fracstr(other_upper),directed_scalar_quadratic_upper=fracstr(scalar_upper),exact_cap_sum=fracstr(cap),gamma1=str(gamma[1]),minimum_LDL_pivot_lower=fracstr(min(row_pivots)),psd_checks_in_row=52,caps=cap_records))
    rational_rows.append(dict(index=idx,interval=[fracstr(lo),fracstr(hi)],matrices=[[[fracstr(x) for x in line] for line in M] for M in MM],lh=[fracstr(x) for x in ellh],lk=[fracstr(x) for x in ellk],t=fracstr(tt),stored_C=fracstr(CC),required_original_fibre_C=decimalstr(allowed),ch_original=fracstr(ch),ck_original=fracstr(ck),diagonal_shift=fracstr(EPS),caps=cap_records))
    print(f'row {idx:02d} arithmetic passed; original fibre allowance {decimalstr(allowed)}; slack over located max {float(slack):.9g}',flush=True)
assert checks==35*52,(checks,35*52)
summary=dict(status='ALL_NON_FIBRE_CONDITIONS_CERTIFIED_UNIVERSAL_FIBRES_STILL_REQUIRED',precision_bits=ctx.prec,python_version=sys.version,python_flint_version=flint.__version__,q=[fracstr(x) for x in q],target_game_norm=fracstr(TARGET),strict_target_below_1000_over_1773=fracstr(F(1000,1773)-TARGET),epsilon=fracstr(EPS),matrix_PD_checks=checks,intervals_cover_exact_unit_interval=True,minimum_LDL_pivot_lower=fracstr(min_pivot),minimum_allowance_above_stored_C_and_located_excess=fracstr(min(F(r['allowance_above_stored_C_and_located_excess']) for r in result)),source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),q_source_sha256=hashlib.sha256(QSOURCE.read_bytes()).hexdigest(),generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seconds=time.time()-start,rows=result)
(HERE/'arithmetic_certificate.json').write_text(json.dumps(summary,indent=2)+'\n')
(HERE/'rational_fibre_inputs.json').write_text(json.dumps(dict(status='FROZEN_RATIONAL_INPUTS_WITH_CERTIFIED_ARITHMETIC_THRESHOLDS',q=[fracstr(x) for x in q],source_sha256=summary['source_sha256'],rows=rational_rows),indent=2)+'\n')
print('ALL_NON_FIBRE_CONDITIONS_CERTIFIED=1',flush=True)
