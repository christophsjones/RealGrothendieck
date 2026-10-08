# Rounding functions calculated by the primal/dual method

This directory describes the rounding functions obtained by the example run of the primal/dual method. This is different from the mixture of rounding functions used in the upper bound certificate.

The final rounding mixture is a distribution on eight **pairs** $(f_r,g_r)$, with

$$
\mu=\sum_{r\in\{42,107,199,265,266,269,270,271\}}w_r\,\delta_{(f_r,g_r)}.
$$

The row numbers identify the saved pool entries. The formulas below describe the frozen mixture used in the subsequent tail calculations. They do not incorporate a new optimization or inverse-correlation preprocessing.

| Pool row $r$ | Weight $w_r$ | Form |
|---:|---:|---|
| 42 | 0.07659708908805163 | Graded linear signs |
| 107 | 0.38070941881309583 | Graded linear signs |
| 199 | 0.3422214571899304 | Graph thresholds, $n=8$ |
| 265 | 0.027934013874205684 | Two-variable polynomial signs |
| 266 | 0.04989420667246783 | Two-variable polynomial signs |
| 269 | 0.05838336532409958 | Symmetric polynomial signs, $n=32$ |
| 270 | 0.017645044351596674 | Symmetric polynomial signs, $n=32$ |
| 271 | 0.04661540468655233 | Symmetric polynomial signs, $n=8$ |

Throughout, use orthonormal probabilists' Hermite polynomials

$$
h_k(x)=\frac{\text{He}_k(x)}{\sqrt{k!}},\qquad
h_0=1,\quad h_1=x,\quad h_2=\frac{x^2-1}{\sqrt2},\quad
h_3=\frac{x^3-3x}{\sqrt6}.
$$

All displayed coordinates are independent standard Gaussians at a single input. A coordinate of grade $d$ has correlation $t^d$ with its corresponding coordinate at another input of correlation $t$. Thus grade-3 Gaussian coordinates below are separate Gaussian variables, not $h_3$ of a grade-1 coordinate. A polynomial term has weighted degree given by the Hermite orders multiplied by their coordinate grades. The convention at a zero of a defining polynomial is immaterial up to Gaussian null sets.

**Rows 42 and 107: graded linear signs.** For each row, let $r_d$ be the constants in the following table, put $s=\sum_d|r_d|$, and use independent grade-$d$ Gaussians $X_d$ and independent private grade-1 Gaussians $U,V$. Set

$$
\begin{aligned}
f&=\text{sign}\!\left(\sum_d\sqrt{|r_d|}\,X_d+\sqrt{1-s}\,U\right),\\
g&=\text{sign}\!\left(\sum_d\text{sign}(r_d)\sqrt{|r_d|}\,X_d+\sqrt{1-s}\,V\right).
\end{aligned}
$$

| $d$ | $r_d$ for row 42 | $r_d$ for row 107 |
|---:|---:|---:|
| 1 | 0.8921904501321847 | 0.8912326131182807 |
| 3 | -0.10780954986781527 | -0.1087673868817175 |
| 5 | 0.0 | 6.070620111124092e-16 |
| 7 | 0.0 | -5.085788243086013e-16 |
| 11 | 0.0 | 2.1631621080173375e-17 |
| 13 | 0.0 | 9.398467226683276e-17 |

All unlisted $r_d$ vanish. Row 42 has $s=1$, so it needs no private noise. Row 107 retains the tiny saved coefficients and residual private noise; discarding them gives the simpler approximation

$$
f_{107}\approx\text{sign}\bigl(\sqrt{0.8912326131182807}\,X_1+
\sqrt{0.1087673868817175}\,X_3\bigr),\qquad
g_{107}\approx\text{sign}\bigl(\sqrt{0.8912326131182807}\,X_1-
\sqrt{0.1087673868817175}\,X_3\bigr).
$$

The exact correlation-kernel formula for these linear pairs is

$$H_r(t)=\frac2\pi\arcsin\!\left(\sum_d r_dt^d\right).$$

**Row 199: an eight-coordinate graph threshold.** Let $X_1,\ldots,X_8,Z_1$ have grade 1 and let $Z_3$ have grade 3. With

$$a=0.9174075833182247,$$

$$
\begin{aligned}
P(x)={}&0.0032620292873087736\,h_1(x)
+0.09737613326063567\,h_3(x)\\
&+0.05450952157793496\,h_5(x)
+0.07453953245559392\,h_7(x)\\
&+0.055787167650815435\,h_9(x)
+0.01730672745975997\,h_{11}(x),
\end{aligned}
$$

set

$$
\begin{aligned}
f_{199}(X,Z_1,Z_3)&=\text{sign}\!\left(\sqrt a\,Z_1+\sqrt{1-a}\,Z_3-\sum_{i=1}^8P(X_i)\right),\\
g_{199}(X,Z_1,Z_3)&=\text{sign}\!\left(\sqrt a\,Z_1-\sqrt{1-a}\,Z_3+\sum_{i=1}^8P(X_i)\right).
\end{aligned}
$$

The sum of $P(X_i)$ is not divided by $\sqrt8$.

**Rows 265 and 266: polynomial signs on two graded Gaussians.** Let $X$ have grade 1 and $Y$ have grade 3. Write

$$
\mathcal I=\{(a,b)\in\mathbb Z_{\ge0}^2: a+3b\le15,\ a+3b\text{ odd}\}.
$$

There are 27 indices. For $r=265,266$, set

$$
f_r(X,Y)=\text{sign}\!\left(\sum_{(a,b)\in\mathcal I}A^{(r)}_{ab}h_a(X)h_b(Y)\right),\qquad
g_r(X,Y)=\text{sign}\!\left(\sum_{(a,b)\in\mathcal I}B^{(r)}_{ab}h_a(X)h_b(Y)\right).
$$

The full saved coefficient arrays are:

| $(a,b)$ | $A^{(265)}_{ab}$ | $B^{(265)}_{ab}$ | $A^{(266)}_{ab}$ | $B^{(266)}_{ab}$ |
|:---:|---:|---:|---:|---:|
| (0,1) | -0.2773497878438673 | 0.27734222003970815 | 0.2784474431131192 | -0.2787580289150735 |
| (0,3) | -2.8600042030220553e-06 | 2.859798489170806e-06 | 0.0011282455340979808 | -0.0011323411070702945 |
| (0,5) | 0.00048538302335466644 | -0.00048529470876766955 | 0.00011895542567642515 | -0.00011573545874141544 |
| (1,0) | -0.898829097704409 | -0.898830795076349 | -0.8981116429722209 | -0.8980412951252327 |
| (1,2) | -0.021819778248902806 | -0.021818728271059875 | -0.02405736283831636 | -0.024104270633874915 |
| (1,4) | 0.0035019383033385865 | 0.0035016571788729716 | 0.0035405668033340245 | 0.003553904107990601 |
| (2,1) | -0.10877059845128823 | 0.10876826065553306 | 0.10745759485813684 | -0.10755086851135208 |
| (2,3) | 0.00627674296811334 | -0.006276313143583789 | -0.006038959132896776 | 0.00605721344966644 |
| (3,0) | -0.28800993687613297 | -0.2880123484110201 | -0.2875178659151116 | -0.2874190395369313 |
| (3,2) | -1.0067942063166342e-05 | -1.0067513783797783e-05 | -0.005185102187583879 | -0.005193977578162221 |
| (3,4) | 0.003250601278895458 | 0.0032503296594498117 | 0.003930910764383742 | 0.003945657647445324 |
| (4,1) | -0.05151474059315139 | 0.051513936007660036 | 0.05387719210873339 | -0.05391054066299014 |
| (4,3) | 0.01169688402457609 | -0.011696150372846837 | -0.012699504766242114 | 0.01273327803831609 |
| (5,0) | -0.10625555263521388 | -0.10625715126173382 | -0.10496118975717113 | -0.10489673513411275 |
| (5,2) | 0.016428418932375503 | 0.016427805669016228 | 0.017163957163590073 | 0.017189573084029575 |
| (6,1) | -1.759025462493703e-05 | 1.759008622503123e-05 | 0.008605539465432072 | -0.00860864088347251 |
| (6,3) | 0.00816481888276673 | -0.008164351162673607 | -0.010249069636176503 | 0.010273368357773456 |
| (7,0) | -0.04447926302027252 | -0.044480232730481524 | -0.04680071847488868 | -0.046759438239799264 |
| (7,2) | 0.02483959183248632 | 0.024838808560419923 | 0.027867676027806876 | 0.02790235515995418 |
| (8,1) | 0.023737441181601135 | -0.0237373599570856 | -0.023598154368231643 | 0.02360061427116618 |
| (9,0) | -1.35896910653415e-05 | -1.3590083725806834e-05 | -0.006737755942893505 | -0.006729931077153171 |
| (9,2) | 0.014967753291271488 | 0.014967359327352034 | 0.018934407175534535 | 0.018953847719348677 |
| (10,1) | 0.03140765766741424 | -0.031407748399418664 | -0.033644593584140674 | 0.033639157139247755 |
| (11,0) | 0.01671578667198192 | 0.01671638923518355 | 0.01694859256023539 | 0.01692422674304487 |
| (12,1) | 0.016822327876304815 | -0.01682248733589806 | -0.020329027425059324 | 0.02032008344780426 |
| (13,0) | 0.02031056995629511 | 0.020311450707410784 | 0.02238798006956791 | 0.022349663330041836 |
| (15,0) | 0.010116888301381502 | 0.010117403632011399 | 0.012707187232646296 | 0.01268193246433082 |

Their two polynomial coefficient arrays are distinct. In particular, a reflected copy of the left polynomial would not reproduce the saved pair exactly.

**Rows 269, 270 and 271: symmetric polynomial signs.** Here $n=32,32,8$, respectively. The coordinates $X_1,\ldots,X_n,Z_1$ have grade 1 and $Z_3$ has grade 3. For an integer partition $\lambda$ into positive parts, let $s=|\lambda|$, $\ell=\ell(\lambda)$, and let $m_d(\lambda)$ be the multiplicity of the part $d$. Define

$$
M_{n,\lambda}=\frac{n!}{(n-\ell)!\prod_{d\ge1}m_d(\lambda)!},\qquad
\Phi_{n,\lambda}(X)=\frac1{\sqrt{M_{n,\lambda}}}
\sum_{\alpha\in\text{Orb}_n(\lambda)}\prod_{i=1}^nh_{\alpha_i}(X_i),
$$

where $\text{Orb}_n(\lambda)$ contains each distinct permutation of $\lambda$ padded with zeros exactly once. For the empty partition, $M=\Phi=1$. Set

$$
\mathcal J_n=\{(\lambda,k,j):\ell(\lambda)\le n,\quad
1\le |\lambda|+k+3j\le15\text{ is odd}\}.
$$

Each of these three pairs has the formula

$$
\begin{aligned}
p_r(X,Z_1,Z_3)&=\sum_{(\lambda,k,j)\in\mathcal J_n}
A^{(r)}_{\lambda,k,j}\,\Phi_{n,\lambda}(X)h_k(Z_1)h_j(Z_3),\\
f_r&=\text{sign}p_r,\qquad
g_r(X,Z_1,Z_3)=f_r(-X_1,\ldots,-X_n,Z_1,-Z_3).
\end{aligned}
$$

Equivalently, the polynomial defining $g_r$ multiplies each displayed coefficient by $(-1)^{|\lambda|+j}$.

For rows **270** and **271**, the $A^{(r)}$ are exactly the saved `left_coefficients`, indexed by `projection_indices`. The saved `right_coefficients` are identical to the corresponding left arrays. There are **2174** coefficients for row 270 and **2047** for row 271. These complete decimal arrays are in the accompanying `grothendieck_rounding_functions.json`; no coefficients are discarded.

For row **269**, the coefficient formula in the same normalized basis is

$$
A^{(269)}\_{\lambda,k,j}=
\widetilde q\_{s+k+3j}(-1)^{s+j}
\sqrt{M_{32,\lambda}}\ C_{\lambda,k+j}
\sqrt{\binom{k+j}{j}}\ a\_\*^{k/2}(1-a\_\*)^{j/2},\qquad s=|\lambda|,
$$

where $a_*=0.9999981420252191$. The source-game vector, in degree order $1,3,\ldots,15$, is

$$
\begin{aligned}
\widetilde q=\bigl(&1,\ -0.8810891383815884,\ 0.5383995613288928,\\
&-0.3399978923773829,\ 0.041522890683156295,\\
&0.24816851838209933,\ -0.3144302166125606,\ 0.12448997441746398\bigr).
\end{aligned}
$$

The 1422 frozen constants $C_{\lambda,k}$ are stored as triples `[lambda,k,value]` in row 269's `nonzero_raw_C`; unlisted constants are zero. These saved decimal constants define the function. The source-game vector above is distinct from the final degree-15 game; changing it to the final game would change the rounding function.

The constants originated as numerical estimates of

$$
\mathbb E\!\left[\text{sign}\!\left(W-\sum_{i=1}^{32}P_*(X_i)\right)
\left(\prod_{i=1}^{\ell(\lambda)}h_{\lambda_i}(X_i)\right)h_k(W)\right],
$$

with independent standard Gaussians and

$$
\begin{aligned}
P_*(x)={}&0.002602096290313642\,h_1(x)+0.09052612817694039\,h_3(x)\\
&+0.047181878574677365\,h_5(x)+0.06955009070152741\,h_7(x)\\
&+0.05677369820796234\,h_9(x)+0.019723270575535756\,h_{11}(x).
\end{aligned}
$$

This expectation explains their origin; the actual frozen function uses the saved constants rather than recomputing those expectations.

Finally, the mixture's correlation kernel is

$$
H_\mu(t)=\sum_r w_r\,\mathbb E[f_r(\mathbf G)g_r(\mathbf G^{(t)})],
$$

using the graded correlation convention above. All formulas and weights are extracted from `final_audit_04/mixture/mixture.json`. The JSON includes both floating-point weights and their exact rational reconstructions for the saved finite coefficient table.
