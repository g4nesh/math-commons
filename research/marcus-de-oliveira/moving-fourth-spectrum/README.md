# Moving the fourth eigenvalue: positive and excluded regions

Task `mdo-moving-fourth-eigenvalue`, agent `ganesh-search`, owner `g4nesh`. The coordinator supplied the separate complex-square, near-zero obstruction and Claude-inequality replays, which were inspected and cross-checked here. These are written proofs with exact executable certificates under one GitHub owner, **not independent-owner review or Lean verification**. No novelty claim or general MDO resolution is made.

Fix the archived four-cycle $U$ and second spectrum $b$. We establish:

| Fourth eigenvalue $z=a_3$ | Result, uniform over arbitrary $(a_0,a_1,a_2)\in\mathbb C^3$ |
|---|---|
| $z=1-i+t$, $-3/4\le t\le2$ | Twelve explicit rational weights are strictly positive and reproduce the determinant. Both endpoints are included. |
| $z=1-i+x+iy$, $|x|,|y|\le1/1024$ | The same support admits strictly positive weights; an exact rank identity and inverse-perturbation bound prove the whole closed square. |
| $z=x+iy$, $|x|,|y|<26085/58260907$ | No probability distribution on **all 24 permutations** can work uniformly over the other three eigenvalues. |

The last statement is about a shared distribution, not pointwise MDO. Every probe used in its separating functional has a separate convex certificate valid for every complex $z$. Regions not listed are not classified.

## 1. Data, source and quantifiers

The fixed data are

$$U=\frac1{625}\begin{pmatrix}
-31&492&-240&300\\
-300&225&500&0\\
240&-180&225&500\\
-492&-256&-180&225
\end{pmatrix},\qquad b=(3-2i,2,1+2i,-i).$$

Set $B=UD_bU^*$, $a=(a_0,a_1,a_2,z)$ and $d(a)=\det(D_a+B)$. Exact unitarity is established by the archived source checks; $D_a$ and $B$ are normal for all complex spectral parameters. A permutation $p$ maps input columns to output rows, with pairing vertex

$$v_p(a)=\prod_{j=0}^3(a_{p(j)}+b_j).$$

The [parent certificate](https://github.com/g4nesh/math-commons/blob/fa8b4ea56e894f2287ce1a8f1c6a04059679de2b/research/marcus-de-oliveira/four-cycle-spectral/README.md#4-coordinator-supplement-three-complex-variables-and-an-exact-limit) gives twelve positive masses at $z=1-i$, uniform over the other three complex eigenvalues. This task continues those masses with $z$. It does not seek one constant distribution for all four variables: the parent's seven-diagonal separator already excludes that stronger requirement.

We use the polynomial-coefficient organization discussed in Claude's [PMSH proposition H2](https://github.com/g4nesh/math-commons/blob/8c1b1fbf189e1b369aae3d034da1fbfc67ce7d64/research/marcus-de-oliveira/minors/pmsh.md). Its Hermitian/collinear conclusions are **not premises** here. Our $b$ is noncollinear: the real planar cross product of $b_0-b_1=1-2i$ and $b_3-b_1=-2-i$ is $-5\ne0$. All identities below are verified directly for this fixed complex spectrum.

## 2. The exact coefficient equations

For $S\subseteq\{0,1,2\}$, the coefficient of $\prod_{i\in S}a_i$ in $d(a)$ is

$$c_S(z)=\det\big(B+\operatorname{diag}(0,0,0,z)\big)[S^c,S^c].$$

The corresponding coefficient in a pairing vertex is

$$c_{S,p}(z)=\prod_{j:\,p(j)\notin S}\big(b_j+z\,\mathbf1_{p(j)=3}\big).$$

Each is affine in $z$: there is only one varying diagonal entry or one varying pairing factor. Matching these eight complex coefficients, together with total mass one, gives 17 real equations when $z=1-i+t$ with $t\in\mathbb R$:

$$M(t)w(t)=r(t),\qquad M(t)=M_0+tM_1,\quad r(t)=r_0+tr_1. \tag{1}$$

Both the twelve-column support matrix and the full 24-column matrix have rank 12 at the base. The eight complex coefficient identities are equivalent to the determinant identity for all complex triples, by multi-affinity. Matching only a finite collection of arbitrary triples would not suffice without this degree argument.

## 3. Explicit nonconstant weights on a closed segment

Use the following support in this order. The row contains the coefficients of $N_j(t)$ in ascending powers $1,t,t^2,t^3,t^4$; an omitted fourth-degree coefficient is zero.

| $p$ | Coefficients of $N_p(t)$ |
|---|---|
| $(0,3,2,1)$ | $24025,38440,16337,1922$ |
| $(1,0,2,3)$ | $291600,550560,259888,34528$ |
| $(1,0,3,2)$ | $694025,1222440,516737,55522$ |
| $(1,3,0,2)$ | $630000,1008000,428400,50400$ |
| $(1,3,2,0)$ | $634375,1029000,462175,70350,2800$ |
| $(2,0,1,3)$ | $974025,1474440,600737,66722$ |
| $(2,0,3,1)$ | $115975,73560,34063,9278$ |
| $(2,3,1,0)$ | $350000,546000,207200,8400,-2800$ |
| $(3,0,1,2)$ | $3975975,6361560,2703663,318078$ |
| $(3,1,0,2)$ | $810000,1296000,550800,64800$ |
| $(3,1,2,0)$ | $175625,225000,24225,-30750,-5600$ |
| $(3,2,1,0)$ | $1090000,1800000,836400,132000,5600$ |

Define

$$D(t)=390625(t+1)(2t+5)(t+5),\qquad w_p(t)=\frac{N_p(t)}{D(t)},$$

and assign zero mass to the twelve omitted permutations. The [exact coefficients](derived.json) include the chosen elimination rows and unreduced Cramer determinant polynomial. At $t=0$, these weights recover every parent mass. They are nonconstant; for example

$$w_{(1,0,2,3)}(0)=\frac{11664}{390625},\qquad
w_{(1,0,2,3)}(1)=\frac{40592}{1171875}.$$

### Polynomial proof of reconstruction

The executable verifies, coefficient by coefficient,

$$\sum_p N_p(t)=D(t),\qquad
\sum_p N_p(t)c_{S,p}(1-i+t)=D(t)c_S(1-i+t)
\quad(S\subseteq\{0,1,2\}). \tag{2}$$

These are ordinary exact polynomial identities; the second set has degree at most five after reduction. Dividing by a nonzero $D$ and expanding in the three free complex eigenvalues proves

$$d(a_0,a_1,a_2,1-i+t)=\sum_p w_p(t)v_p(a_0,a_1,a_2,1-i+t). \tag{3}$$

The derivation used exact elimination, Cramer's rule and polynomial interpolation. Each unreduced determinant has degree at most 12, so 13 rational evaluations determine it exactly. A common polynomial factor was cancelled. **Replay does not trust the derivation:** it checks (2) directly, reconstructs the recorded minor polynomial using its degree bound, and verifies the signs below. No numerical or optimization solver is required.

### Positivity for every $t\in[-3/4,2]$, including endpoints

All three denominator factors are positive on this interval. For each numerator $N(t)=\sum_{k=0}^n c_kt^k$, put $t=L+(R-L)x$ with $L=-3/4$, $R=2$, $0\le x\le1$. The exact conversion to Bernstein form is

$$N(L+(R-L)x)=\sum_{i=0}^n\beta_i\binom ni x^i(1-x)^{n-i},$$

$$\alpha_j=\sum_{k=j}^n c_k\binom kj L^{k-j}(R-L)^j,
\qquad \beta_i=\sum_{j=0}^i\alpha_j\frac{\binom ij}{\binom nj}.$$

Every computed $\beta_i$ is strictly positive; the complete rational lists are in [verification.json](verification.json). The smallest across all twelve lists is $114359/32$. Bernstein basis functions are nonnegative and sum to one throughout the **closed** unit interval. Therefore every numerator is at least $114359/32>0$, including at $x=0,1$.

Also $D(t)\le D(2)=73828125$ because its positive factors increase on this interval. Consequently every supported mass obeys the explicit uniform bound

$$w_p(t)\ge\frac{114359}{2362500000}>0.$$

Together with (2), this proves the requested nonnegative convex representation on the entire closed segment $a_3\in[1/4,3]-i$, uniformly over all three other complex eigenvalues. Direct endpoint, Boolean-corner and off-grid complex evaluations are extra regressions, not substitutes for the polynomial and Bernstein proof.

### Singular formulas and support limits

The denominator vanishes at $t=-1,-5/2,-5$, all outside the certified interval. Exact elimination on the specified support gives:

| $t$ | $a_3$ | Coefficient rank | Augmented rank | What follows |
|---|---|---:|---:|---|
| $-1$ | $-i$ | 11 | 12 | No signed solution on this twelve-permutation support. |
| $-5/2$ | $-3/2-i$ | 11 | 11 | The equations are linearly consistent; the displayed rational formula is singular. No nonnegative solution is asserted or excluded here. |
| $-5$ | $-4-i$ | 11 | 12 | No signed solution on this twelve-permutation support. |

A vanishing elimination determinant or formula denominator is not itself infeasibility. Failure of this support leaves other permutations available. The segment above is a sufficient interval, not a claimed maximal feasible interval.

## 4. A closed two-dimensional positive neighborhood

The coordinator's [positive_box.py](positive_box.py) proves a complementary region, using the same twelve supports. Let $z=1-i+x+iy$. Write the normalization plus eight Boolean-corner equations as a real $17\times13$ augmented matrix, affine in $x,y$.

**Global compatibility identity.** Every $13\times13$ minor has degree at most 13 in each real variable. Exact augmented ranks are at most 12 on the $14\times14$ grid $x,y\in\{0,\ldots,13\}$. Tensor-product polynomial interpolation therefore forces every such minor to vanish identically. Thus the augmented system has rank at most 12 for every $x,y$. This is an identity check with an explicit degree bound, not sampled evidence for an unbounded claim.

At the center, select twelve independent coefficient rows, giving an invertible matrix $A_0$ and positive base masses $w^0$. Put $W=\operatorname{diag}(w^0)$ and write $w=w^0+Wv$. After multiplying by $W^{-1}A_0^{-1}$, the selected equations become

$$(I+xM_x+yM_y)v=xr_x+yr_y.$$

The checker obtains exact norm bounds

$$H=\max_i\sum_j(|(M_x)_{ij}|+|(M_y)_{ij}|)=\frac{875739203}{1229375},$$

$$K=\max_i(|(r_x)_i|+|(r_y)_i|)=\frac{6720}{4639}.$$

For the closed square $|x|,|y|\le\rho=1/1024$, we have $\rho H<1$. The Neumann-series inverse bound proves selected-system invertibility and

$$\|v\|_\infty\le\frac{\rho K}{1-\rho H}
=\frac{2753800000}{592463385761}<1.$$

Hence every $w_j=w_j^0(1+v_j)$ is strictly positive, even on the square boundary. The selected coefficient matrix has rank 12, while the full augmented rank is at most 12; consequently the selected solution satisfies all rows. The eight corner identities then imply the identity for every complex triple by multi-affinity. The same masses are used at every triple; interpolation factors are not probability weights.

The optional `--offset X Y` command returns exact masses at rational displacements in the certified square. Tests verify all four square corners again using the complementary-principal-minor coefficient basis rather than the corner basis.

## 5. An open region where no shared distribution exists

The coordinator's [zero_fiber.json](zero_fiber.json) supplies six binary triples $\varepsilon_r\in\{0,1\}^3$ and rational real-linear functionals $\ell_r$. For a fixed permutation define

$$g_p(z)=\sum_r\ell_r\big(d(\varepsilon_r,z)-v_p(\varepsilon_r,z)\big).$$

Because the determinant and vertices are affine in complex $z$, each gap is a real-affine function $g_p(x+iy)=c_p+u_px+v_py$. [zero_fiber.py](zero_fiber.py) checks every permutation and obtains

$$\min_p c_p=\frac{1210344}{15625}>0,\qquad
R=\min_{p:\,|u_p|+|v_p|>0}\frac{c_p}{|u_p|+|v_p|}
=\frac{26085}{58260907}.$$

For $|x|,|y|<R$, each $g_p(x+iy)>0$. If a distribution $w(z)$ represented the determinant for every complex triple at that fixed $z$, it would represent these six triples. Applying the same sum of real functionals would yield $0=\sum_p w_p(z)g_p(z)>0$, a contradiction. This excludes **every** distribution on all 24 permutations throughout the stated open square around zero. No boundary claim is made at $|x|=R$ or $|y|=R$.

Each probe also carries its **own** positive distribution that reconstructs its determinant for every complex $z$. Both affine coefficients are verified exactly, including by a separate direct cofactor calculation in the tests. Thus even these witness inputs individually satisfy MDO. The obstruction concerns simultaneously fitting several inputs with one shared distribution. It does not imply a determinant outside its pointwise spectral hull.

## 6. Concurrent PMSH work is context, not a premise

Claude's [November 25 artifact](https://github.com/g4nesh/math-commons/blob/38eb4ff58c9b3d1caefb3ab372f9dd37cca15ab2/research/marcus-de-oliveira/pmsh/equivalence.md) supplies a proposed computer-assisted equivalence between PMSH(4) and its inequality (★). That artifact explicitly leaves the universal inequality open. We do not independently certify its equivalence argument here.

The coordinator's [claude_star.py](claude_star.py) only replays that precise inequality on this fixed **real** four-cycle and its row/column/transpose relabelings. All 1,152 labeled variants pass with exact rational arithmetic: original slack $5779/15625$, minimum slack $15136/390625$. This finite cross-check neither proves (★) for arbitrary unitaries nor supplies a Hermitian/collinear hypothesis for our noncollinear $b$. None of the region proofs above uses it.

## 7. Reproduction and selected handoff

Run from the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 research/marcus-de-oliveira/moving-fourth-spectrum/moving.py
PYTHONDONTWRITEBYTECODE=1 python3 research/marcus-de-oliveira/moving-fourth-spectrum/positive_box.py
PYTHONDONTWRITEBYTECODE=1 python3 research/marcus-de-oliveira/moving-fourth-spectrum/positive_box.py --offset 1/1024 1/1024
PYTHONDONTWRITEBYTECODE=1 python3 research/marcus-de-oliveira/moving-fourth-spectrum/zero_fiber.py
PYTHONDONTWRITEBYTECODE=1 python3 research/marcus-de-oliveira/moving-fourth-spectrum/claude_star.py
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s research/marcus-de-oliveira/moving-fourth-spectrum -p 'test_moving.py' -v
```

`moving.py --derive --output derived.json` re-derives the coefficients; use the full owned output path when running from the repository root. `moving.py --report verification.json` regenerates verification details similarly. Default replay uses the stored certificate next to the script and does not perform an optimization. The code uses standard-library exact rationals and the archived Gaussian-rational arithmetic; no floating tolerance, external model or paid compute is used. The main executable enforces a 60-second CPU cap, and all other replays have explicit finite loops. The research budget was one 20-active-minute pass. Python 3.9.6 was tested; checked artifacts contain no interpreter-version or timing metadata.

All **12 tests passed** locally. They cover exact re-derivation and saved-report replay, coefficient identities, Bernstein conversion/signs, closed endpoints and off-grid complex triples, source masses, corruption/inexact/out-of-range rejection, singular-support distinctions, both region certificates through separate coefficient constructions, all four positive-square corners, and the finite Claude-inequality replay. These checks do not replace independent-owner review.

**Selected next question:** follow the diagonal path $a_3=s(1-i)$ from the positive neighborhood at $s=1$ toward the excluded neighborhood at $s=0$. Construct piecewise rational certificates and exact separating functionals while allowing support changes. Do not assume feasibility is monotone or that failure of the current twelve supports excludes the full 24-permutation problem. This would narrow the boundary of the shared-distribution property; pointwise MDO remains a different question.

The coordinator owns all Git, CI, shared research navigation and actual nonoverlapping usage. No separate usage import is performed by this agent.
