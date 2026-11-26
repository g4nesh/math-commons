# Exact fourth-eigenvalue boundary

For the archived four-cycle unitary $U$ and noncollinear spectrum $b$ below, put
$a_3=s(1-i)$, where $0\le s\le1$. There exists **one probability distribution on
permutations that gives the determinant identity for every complex triple
$(a_0,a_1,a_2)$ if and only if $s\in[\alpha,1]$**, where $\alpha$ is the unique
positive root of

$$f(s)=1808s^3+9419s^2+8398s-1400.$$

The exact isolating interval is $143/1000<\alpha<144/1000$; its decimal value is
approximately $0.143106204978043899$. The decimal is descriptive, not a premise.
The boundary itself is feasible, with eleven positive masses and one zero mass
in the first twelve-term formula.

This classifies a **uniform coefficient-matching problem for this particular
matrix and second spectrum**. Failure below $\alpha$ means a single mixture
cannot work for all three free eigenvalues. It does not exhibit any particular
triple whose determinant lies outside its own permutation hull. It is not a
counterexample to the Marcus–de Oliveira conjecture or a general solution.

## Fixed data and permutation convention

$$
U=\frac1{625}\begin{pmatrix}
-31&492&-240&300\\
-300&225&500&0\\
240&-180&225&500\\
-492&-256&-180&225
\end{pmatrix},\qquad b=(3-2i,2,1+2i,-i).
$$

Set $B=U\operatorname{diag}(b)U^*$ and $A=\operatorname{diag}(a_0,a_1,a_2,s(1-i))$.
The script checks exact unitarity and the archived source data. Both matrices
are normal. This $b$ is not collinear, so guarantees requiring Hermitian or
collinear spectra are not used.

A permutation $p$ maps **column to row**. Its vertex is

$$v_p(a)=\prod_{j=0}^3(a_{p(j)}+b_j).$$

The question is whether $w_p(s)\ge0$, $\sum_pw_p(s)=1$ can satisfy
$\det(A+B)=\sum_pw_p(s)v_p(a)$ for all $(a_0,a_1,a_2)\in\mathbb C^3$.
Weights may depend on $s$, but not on that triple.

## Reduction to exact coefficient equations

For each $S\subseteq\{0,1,2\}$, the coefficient of $\prod_{i\in S}a_i$ in the
determinant is the principal minor on the complementary coordinates of
$B+\operatorname{diag}(0,0,0,s(1-i))$. In the vertex polynomial it is

$$\prod_{j:p(j)\notin S}\left(b_j+\mathbf1_{p(j)=3}s(1-i)\right).$$

Normalization followed by the real and imaginary parts of these eight
coefficients gives a real system $M(s)w=d(s)$ with 17 rows and 24 columns.
Every entry is affine in $s$. The ordering is normalization first, then subsets
by increasing size and lexicographic order, with real before imaginary.
Equality of these coefficients proves equality for every complex triple;
finitely many sampled triples alone would not prove the interval result.

## Positive side: three exact supports cover the entire interval

[positive.json](positive.json) stores three supports, each with twelve
permutations, an integer-coefficient denominator $D(s)$, and twelve
integer-coefficient numerators $N_j(s)$. Coefficients are in ascending powers.
Their weights are $N_j(s)/D(s)$.

Let $\beta$ be the unique positive root of

$$g(s)=41s^3+212s^2+183s-36.$$

It satisfies $164/1000<\beta<165/1000$, approximately $0.1644110316555803695$.
Both $f$ and $g$ are strictly increasing on $[0,1]$, since every nonconstant
coefficient is positive. Exact endpoint signs prove these isolations and
$\alpha<\beta<97/100<1$.

The three certified supports cover:

| Formula | Closed interval | Sign argument |
|---|---|---|
| First | $[\alpha,\beta]$ | Two numerators are $-729g(s)$ and $16f(s)$; every other numerator and the denominator is strictly positive on $[143/1000,165/1000]$. |
| Second | $[\beta,97/100]$ | Its first numerator is $g(s)q(s)$; $q$, all other numerators, and the denominator are strictly positive on $[164/1000,97/100]$. |
| Third | $[97/100,1]$ | Every numerator and the denominator is strictly positive throughout the interval. |

Here

$$q(s)=18954+14580s-16038s^2-14580s^3-2916s^4.$$

For each positivity assertion, `boundary.py` converts the polynomial exactly to
Bernstein form on the stated rational interval and verifies that **all** control
coefficients are positive. Bernstein basis functions are nonnegative and sum to
one, so this proves positivity everywhere on the interval, not only on a grid.
For the first formula, $D(s)=390625(4+s)(1-s^2)$, consistently positive on its
certified interval. The replay checks the general denominator sign directly.

For each of the three formulas, all 17 identities

$$\sum_j M_{ij}(s)N_j(s)=d_i(s)D(s)$$

are checked as exact polynomial coefficient identities. Thus normalization and
all eight complex determinant coefficients hold identically. This gives 51
real polynomial identities.

At $\alpha$, the $16f$ mass is exactly zero and the other eleven are positive.
At $\beta$, each neighboring formula has one zero and eleven positive masses.
Either formula works there; they need not select the same distribution.
At $97/100$, either neighboring formula works and all its masses are positive.
Consequently there is no unproved endpoint or gap between formulas. No
continuity of the particular selected weights is asserted. The switch at
$97/100$ is a convenient certificate choice, not a second feasibility boundary.

## Negative side: one functional excludes every permutation

The coordinator independently derived a cofactor functional from the first
positive support. This contribution was regenerated and checked here using
exact rational arithmetic, including all 24 columns.

Select rows

$$R=(0,5,2,3,4,1,6,9,8,10,11,12)$$

of $M(s)$, and the twelve columns of the first support, to form a square matrix
$C(s)$. Let $L_s(x)$ be the determinant obtained by replacing its column with
zero-based index 1 by $x_R$, divided by

$$h(s)=\frac{12544}{5}s^3(s+4).$$

For $s>0$ this denominator is positive, and $L_s$ is a real linear functional of
$x$. [dual.json](dual.json) records the exact polynomial values of $L_s$ on
all 24 columns of $M$ and on $d$. The numerator determinants have degree at most
12 because all their entries are affine. Thirteen exact interpolation values
therefore establish every recorded identity everywhere. The verifier computes
those 325 determinants directly; no optimizer is needed for replay.

Every column polynomial is nonnegative throughout $[0,144/1000]$, certified by
nonnegative Bernstein coefficients. The target polynomial is exactly

$$L_s(d(s))=32(1+3s+s^2)f(s).$$

For $0<s<\alpha$ this is negative, whereas every convex combination of column
values is nonnegative. Hence $M(s)w=d(s)$ has no nonnegative solution there,
even allowing all 24 permutations. This is a separating-functional proof for
an entire interval, rather than an inference of monotonicity from sampled LPs.

At $s=0$ the displayed cofactor normalization vanishes, so it is not used.
A separate exact 17-row functional is stored as `zeroDual`. The verifier checks
$y^TM(0)\le0$ on every column and

$$y^Td(0)=3472/33125>0.$$

This excludes the endpoint and completes the exact classification $[\alpha,1]$.
At $s=\alpha$ the cofactor target value is zero, consistently with the explicit
convex certificate; the strict separation claim deliberately excludes it.

## Replay and use

Python 3.9 or newer, standard library only. From the repository root:

```sh
python3 -B research/marcus-de-oliveira/fourth-boundary/boundary.py
python3 -B -m unittest discover -s research/marcus-de-oliveira/fourth-boundary -p 'test_*.py' -v
python3 -B research/marcus-de-oliveira/fourth-boundary/boundary.py --weights 3/20
```

The first command verifies source data, the positive formulas, both separating
certificates, and 32 direct Boolean-corner reconstructions. Its output matches
[verification.json](verification.json) exactly. Ten tests also exercise
non-Boolean complex triples, all three formulas, rejected inexact inputs,
tampered masses/duals, missing columns, interpolation degree bounds, and the singular zero endpoint.
Local validation used Python 3.9.6; the report excludes machine/version metadata
so it can replay identically elsewhere.

The importable `weights_at(s)` accepts an integer or `Fraction`, rejects
parameters outside the feasible segment, and returns `(column_to_row, mass)`
pairs. For the algebraic endpoint $\alpha$, the first JSON formula itself is
an exact certificate in $\mathbb Q(\alpha)$; the rational-input API does not
pretend to encode an irrational number as a float.

Optional discovery, separate from the solver-free proof replay:

```sh
python3 -B research/marcus-de-oliveira/fourth-boundary/discover.py 0 7/50 3/20 1/4 1
```

`discover.py` is a bounded exact phase-I simplex with Bland pivot selection.
It verifies each returned dual or feasible distribution against all 17 original
rows. Discovery initially considered the 21 rational parameters $k/20$ and
then $1/8,13/100,7/50$ near the transition. The three formulas were obtained by
exact Cramer interpolation and common-polynomial cancellation from the supports
at $3/20,1/4,1$. Those samples chose useful supports; they are not used as
interval evidence. The independent coefficient and sign checks above supply
that evidence. The task's twenty-active-minute limit covers this bounded pass;
no paid compute was used.

## Attribution and remaining scope

The archived fixture and coefficient system come from the repository's
[single-cycle construction](../single-cycle/README.md),
[three-variable fiber certificate](../four-cycle-spectral/README.md), and
[previous moving-fourth-eigenvalue work](../moving-fourth-spectrum/README.md).
The latter is pinned at
[ac81a919](https://github.com/g4nesh/math-commons/blob/ac81a919c68fc685e1e104b350bac9ca2499c5c5/research/marcus-de-oliveira/moving-fourth-spectrum/README.md).
That work gave a positive neighborhood near $1-i$ and an excluded neighborhood
near zero. Here the positive basis continuation and the coordinator's cofactor
dual close the entire intervening path. The general polynomial-matching
viewpoint also appears in [PMSH notes](../pmsh/equivalence.md); their spectral
hypotheses are not used to justify this noncollinear fixture.

The theorem is the precisely stated fixed-data path classification. No claim
of overall literature novelty, formal proof-assistant verification, general
MDO resolution, or independent-owner approval is made. The computational
cross-checks were produced under the same GitHub owner; review by a different
owner remains pending. A useful next question is the shape of the feasible
region for unrestricted complex $a_3$, rather than assuming a radial or convex
geometry from this one segment.
