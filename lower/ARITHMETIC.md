# Arithmetic and acceptance boundary

## Analytic calculations

The recorded analytic replay uses 384-bit intervals with MPFR endpoints.
`backend/balls.hpp` represents an interval by a lower and an upper endpoint.
Rational inputs are rounded outward. Addition and subtraction use the
corresponding endpoint formulas; multiplication checks signs or the four
endpoint products; division rejects a denominator interval containing zero.
Square root, exponential, error functions, and pi use MPFR's directed rounding.
The matrix product uses sign-aware directed fused multiply-add operations.
Exact rational linear algebra uses GMP rational numbers where required by the
unchanged original checker.

`backend/flint.pyx` is the narrow compatibility interface consumed by the
original programs. The name of that interface does not mean that the arithmetic
is Arb: its version string identifies the independent MPFR implementation.
There is no claim that this delivery executed an alternative Arb backend.

The Gaussian proof needs no numerical quadrature: polynomial integrals are
reduced to the Gaussian CDF/density, normalized-Hermite antiderivatives, and
exact polynomial product identities. Both infinite tails are included. A
positive quadratic Bernstein rule gives a finite sign problem; its middle
control functional is 2p(midpoint)-(p(left)+p(right))/2. The integrated Peano
remainder uses exp(Rh) h^3/96 and the third Gaussian derivative. The explicit
interpolation and tail enclosure is added, not ignored. Input rounding and
integer graph quantization have separate uniform error allowances.

## Integer calculations

`tools/exact_graph.cpp` performs no floating-point arithmetic. It constructs
integer quadratics from the fixed-point input, computes rounded rational
quadratic caps with uniform error allowances, and produces exact graph-cut
upper bounds. Signed 128-bit integers are used for intermediate products;
stored graph coefficients and capacities use signed 64-bit integers. Explicit
input and aggregate bounds reject computations outside the safe range.

The flow algorithm's result is checked independently by verifying capacities,
reverse-edge consistency, conservation at every vertex, and equality between
the flow value and the returned cut capacity. This gives a weak-duality upper
bound with equality, rather than reliance on the optimizer's return status.
The persistency rule used before binary branching is justified in the accompanying manuscript.

`tools/replay.py` is independent of the numerical search stopping logic. It
reconstructs the original graph, exact scalar coordinates and any symmetry
prefix; checks complete closed rectangle covers; regenerates every quadratic
cap; and traverses every binary subtree, checking both children, persistent
assignments, and all terminal inequalities. Active, missing, repeated,
extraneous, or unsupported nodes prevent acceptance.

The final rational assembler requires all 35 newly regenerated universal
fibres, the freshly regenerated 1,820 non-fibre matrix checks, angular-cap and
omitted-degree checks, the exact scalar-interval cover, and the negative
controls. A historical success flag cannot replace one of these computations. The additional frozen-row limitation check evaluates an explicit valid sign witness at x=381/400 using directed lower bounds. This supplementary limitation is not an upper bound on the Grothendieck constant.

## Tests and limits of the certificate claim

The replay includes exact-rational endpoint arithmetic tests, interval matrix
products, exact inverse identities, 484 Hermite product identities, Gaussian
Gram, third-derivative, Peano-kernel, and positive Bernstein identities, and exhaustive small integer graph tests.
Malformed-evidence controls test failure on incomplete covers and branches,
incorrect integer bounds, unjustified persistent assignments, insufficient
analytic errors, and broken source linkage. A separate insufficient-repair
control executes the original matrix checker with repair 10^-12 and must fail
at the designated matrix, while the positive run uses 10^-7.

These tests supplement the mathematical argument and implementation review.
They are not a machine-checked proof of the compiler, MPFR, GMP, NumPy, Python,
the operating system, or the processor. The certificate is a reproducible
computer-assisted proof with explicit source and data, not a Lean/Coq formal
proof. No assertion is made that unexecuted arithmetic backends or platforms
have passed replay. The exact versions actually used are recorded separately.
