# Certified arithmetic reduction for K_G > 1.773

This packet certifies every non-fibre condition for the frozen degree-21 candidate. The global fibre inequalities are a separate, still-required hypothesis. No Grothendieck bound is asserted by this packet alone.

The exact rational target is

    v = 1000/1773 − 1/10^9.

The exact rational multiplier q has q_1=1 and the terminating decimal coefficients frozen in the supplied input. Every JSON decimal, including the matrix entries and slab endpoints, is parsed as a terminating rational. The generator binds both input hashes.

## Required universal fibre statements

Let h_j=He_j/sqrt(j!) be the normalized probabilists' Hermite polynomials, and let a measurable pair (u(S),v(S)) take values in {(±1,0),(0,±1)}. For j=0,...,21 set

    m_j = E[u(S)h_j(S)],  n_j = E[v(S)h_j(S)],  S~N(0,1).

For each row i, the original frozen fibre objective is

    F_i(m,n) = m_even^T HE_i m_even + m_odd^T HO_i m_odd
             + n_even^T KE_i n_even + n_odd^T KO_i n_odd
             + t_i m_1 − 2 lh_i^T m_odd − 2 lk_i^T n_odd.

The separate proof must establish F_i≤C_i for every such measurable pair, where `required_original_fibre_C` in `rational_fibre_inputs.json` gives an exact terminating rational C_i. Those thresholds already reserve all arithmetic and matrix perturbation costs. The matrices in that file are the ORIGINAL frozen matrices; the fibre prover should not add the diagonal shift itself.

## Exact circle caps and omitted tails

For an odd N, enumerate every antipodal sign pattern of length 2N with first sign +1, and calculate the integer numerators

    A_d(s) = sum_{j=0}^{N−1} s_j s_{j+d},  d=1,...,(N−1)/2.

For the exact rational angular weights w, the generator replaces the solver's stored cap with the exact rational maximum B=max_s sum_d w_d A_d(s)/N. The Gaussian rotation argument gives

    sum_{m odd} gamma_m [W_m(h)+W_m(k)] ≤ B,
    gamma_m = sum_d w_d cos(pi d/N)^m.

This follows by applying the finite cut inequality to f and g and using h=(f+g)/2, k=(f−g)/2. The cap sums for N=7 and N=9 are used directly. All trigonometric values are enclosed by Arb.

For every individual cap, strict positivity is checked for odd m=23,25,...,101. The generator also proves

    w_1 − sum_{d≥2} max(0,−w_d)(c_d/c_1)^101 > 0,
    1>c_1>c_2>...>0.

Dividing gamma_m by c_1^m proves gamma_m>0 for every integer m≥101. Thus every omitted odd Hermite degree is controlled, with no numerical truncation hypothesis.

## Directed matrix repair

Set epsilon=1/10^7. Each of the four matrices in each row is replaced by M'=M+epsilon I, and ch,ck are replaced by ch+epsilon,ck+epsilon. The generator checks, by directed interval LDL^T with all pivots strictly positive:

1. M'>0 and 10I−M'>0 for all four matrices;
2. every transverse-chaos domination matrix from the finite model;
3. both augmented mean matrices with their shifted bottom-right entries.

There are 52 strict positive-definiteness checks per row, 1,820 in total. The upper bounds 10I−M'>0 are additional numerical checks and are not needed by the reduction.

The transverse diagonals are retained exactly through total degree 21. At larger total degrees the actual residual is −Gamma_m≤0, so replacing it by zero is safe. For transverse degrees beyond 21, positive semidefiniteness of M' provides the same domination. Terms with distinguished-coordinate Hermite degree above 21 have nonpositive residual and may be dropped.

Joint Parseval gives

    sum_{j=0}^{21}(m_j²+n_j²) ≤ E(u²+v²)=1.

Consequently the matrix repair increases the fibre objective by at most epsilon. The two mean costs increase by 2epsilon, so the total arithmetic repair costs at most 3epsilon.

For clarity about the mean-block sign: the positive semidefinite block [[M'−diag(c),ell],[ell^T,cost']] is evaluated at (y,−1), giving y^T(diag(c)−M')y≤−2ell^T y+cost'. This matches the minus signs in the displayed fibre objective.

## Exact scalar coverage

After alignment of the distinguished Gaussian direction, ||Pi_1 h||=nu x with nu=sqrt(2/pi) and x∈[0,1]. The exact rational slab endpoints form a gap-free partition of [0,1]. The generator proves (1−Gamma_1)nu²>0 in every row, hence the scalar quadratic

    (1−Gamma_1)nu²x² − t_i nu x

is convex and its maximum is attained at a slab endpoint. Both endpoints are enclosed directly with Arb. Their rigorous rational upper bound, plus the exact cap sums, ch+ck and 3epsilon, is denoted A_i.

The allowed original fibre threshold is the downward rounding to twelve decimal places of v−A_i. Therefore, exactly,

    C_i+A_i ≤ v < 1000/1773.

If all 35 required universal fibre inequalities are established, the standard joint-fibre reduction then gives ||T_q||_(infinity→1)≤v, and consequently K_G≥1/v>1.773. The present certificate establishes the arithmetic implication and its exact thresholds, not those universal fibre inequalities.

## Reproduction and negative control

Run `certify_arithmetic.py` with the Arb-capable Python runtime used in this workspace. At 384-bit precision it prints `ALL_NON_FIBRE_CONDITIONS_CERTIFIED=1` and writes the certificate and rational fibre inputs. Replacing epsilon by 1e−12 is a negative control: it fails on row 0, KE, transverse degree 3. The unmodified floating matrices are therefore not silently treated as positive semidefinite.
