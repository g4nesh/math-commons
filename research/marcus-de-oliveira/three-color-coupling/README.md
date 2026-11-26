# A three-color coupling beyond common-minor certificates

This gives an explicit **nonempty open subclass of the missing pattern
$(2,2,2)$** in Conjecture F of the [minors report](../minors/report.md#4-the-three-eigenvalue-conjecture).
For every unitary in a stated ball about an exact dense $6\times6$ matrix,
there is a strictly positive distribution on all 90 colorings with two sites
of each color. The distribution matches every grouped squared-minor statistic.
At the same time, **every unitary in that ball has no distribution matching all
individual squared minors**.

Thus the grouped certificate method succeeds on a class that the full
common-minor method cannot cover. The result is a precise subclass, not a proof
of Conjecture F for all unitaries. Padding and a proved pendant-gate operation
also give structured classes in arbitrary dimension. No literature novelty
claim or independent-owner verification is asserted.

## Exact statement

Let $\mathcal K$ be the 90 maps $\kappa:\{0,\ldots,5\}\to\{0,1,2\}$ with
exactly two occurrences of each color. Column blocks are
$J_0=\{0,1\}$, $J_1=\{2,3\}$, $J_2=\{4,5\}$.
The exact matrix $U_0$, all 90 positive masses $w^0_\kappa$, all 90 canonical
signed masses, and all 90 correction coefficients are in
[interior.json](interior.json). Their deterministic construction is below.

Write

$$\epsilon=10^{-24},\qquad \delta_0=\frac1{6\cdot10^{25}},$$

and let $\Gamma$ be the exact positive rational number

$$\Gamma=\frac{64066719632048727152801657518817785285266217556555697629472}
{1564688922443906437535942723796925003281411937642189687515625}.$$

It is approximately $0.04094534$. For **every complex unitary** $U$ satisfying

$$\max_{i,j}|U_{ij}-(U_0)_{ij}|\le\delta_0,$$

an explicit distribution $w(U)$ satisfies:

1. Every mass is at least $\epsilon/2=1/(2\cdot10^{24})$.
2. For every subset $S$ and every three-component count vector $\nu$,

   $$\sum_{\kappa:\,|S\cap\kappa^{-1}(j)|=\nu_j\ \forall j}w_\kappa(U)
   =\sum_{T:\,|T\cap J_j|=\nu_j\ \forall j}|\det U[S,T]|^2.$$

3. A fixed 21-term common-minor inequality is violated by at least $\Gamma/2$.

The radius is deliberately conservative and very small. It certifies an actual
neighborhood; it is not an estimate of the full feasible region or its volume.
The closed ball contains a nonempty relatively open subset of $U(6)$.
The result permits arbitrary complex values for the three eigenvalues and
arbitrary complex eigenvalues of the other matrix. They need not be collinear,
real, or distinct.

## 1. Canonical signed masses and a fixed correction

For an arbitrary unitary $U\in U(n)$ with column blocks $J_0,\ldots,J_{m-1}$ of
sizes $r_0,\ldots,r_{m-1}$, fix a coloring $\kappa$ with these counts and put
$X_j=\kappa^{-1}(j)$. Let $\eta(X)$ be the sign of the row permutation that
lists the $X_j$ consecutively, each in increasing order. Define

$$c_\kappa(U)=\eta(X)\overline{\det U}\prod_{j=0}^{m-1}\det U[X_j,J_j],
\qquad r_\kappa(U)=\operatorname{Re}c_\kappa(U).$$

These are the canonical signed coloring weights discussed in the report;
they need not be nonnegative. They sum to one and match all grouped statistics.
Here is a direct proof of the latter assertion.

Let $Q_j=U\Pi_{J_j}U^*$ and introduce a separate variable $x_{ij}$ for every
coordinate and color. Expanding the determinant by the column blocks gives

$$\det\left(\sum_j D_{x_j}Q_j\right)
 =\sum_\kappa c_\kappa(U)\prod_i x_{i,\kappa(i)}.$$

Set $x_{ij}=c_j$ when $i\in S$ and $x_{ij}=1$ otherwise. The matrix has identity
rows outside $S$, so its determinant is
$\det(\sum_j c_jQ_j)[S,S]$. Cauchy–Binet expands this as the polynomial whose
coefficients are the grouped squared minors. They are real. Taking real parts
therefore proves the required statistics for $r(U)$ and proves their sum is one.
This argument is an identity of polynomials, not a positivity assertion.

At our exact base $U_0$, the verifier supplies a positive law $w^0$ with the same
statistics. Define the fixed vector

$$h_\kappa=w^0_\kappa-r_\kappa(U_0).$$

Every grouped sum of $h$ is zero. Consequently the explicit formula

$$\boxed{w_\kappa(U)=r_\kappa(U)+h_\kappa}$$

has the correct statistics and total mass one for **every** unitary $U$;
nonnegativity is the only remaining issue. This fixed correction is essential:
21 canonical real masses at $U_0$ are negative. The exact coefficients $h$ are
stored in `interior.json`, rather than depending on a new optimization at each
matrix.

## 2. A quantitative stability proof

Assume $U,V$ are unitary and $\max|U_{ij}-V_{ij}|\le\delta$. Every minor of a
unitary has modulus at most one. Telescope a determinant one column at a time
and use Hadamard's inequality. For an $r\times r$ submatrix,

$$|\det U[S,T]-\det V[S,T]|\le r\sqrt r\,\delta.$$

All unchanged truncated columns have Euclidean norm at most one; the changed
column has norm at most $\sqrt r\,\delta$. For the full determinant in dimension
six this is at most $18\delta$, using $\sqrt6<3$. Each two-minor changes by at
most $4\delta$. Telescoping the product of the full determinant and three
two-minors therefore gives

$$|r_\kappa(U)-r_\kappa(V)|\le30\delta.$$

Since all base masses are at least $\epsilon$, the corrected weights satisfy
$w_\kappa(U)\ge\epsilon-30\delta_0=\epsilon/2$ throughout the stated ball.
This is a proof for all complex unitaries in the neighborhood. A finite replay
of nearby matrices is only an additional regression check.

More generally, this argument proves an explicit stability lemma for any
strict pattern certificate in any dimension, with the bound
$n\lceil\sqrt n\rceil+\sum_j r_j\lceil\sqrt{r_j}\rceil$ in place of 30.
It requires a certificate with every allowed coloring mass strictly positive;
it does not infer interior feasibility from a certificate with zero masses.

## 3. Exact construction and verification of the base

Start with $\operatorname{diag}(U_4,I_2)$, where $U_4$ is the archived
[four-cycle obstruction](../single-cycle/README.md#5-a-four-cycle-with-no-common-minor-certificate):

$$U_4=\frac1{625}\begin{pmatrix}
-31&492&-240&300\\-300&225&500&0\\
240&-180&225&500\\-492&-256&-180&225
\end{pmatrix}.$$

Left-multiply, in this order, by gates on coordinate pairs

$$(0,4),(1,5),(2,4),(3,5),(0,5),(1,4),(2,5),(3,4),(4,5).$$

Each gate has block

$$\begin{pmatrix}c&s\zeta\\-s\overline\zeta&c\end{pmatrix},
\quad c=\frac{9999}{10001},\quad s=\frac{200}{10001},$$

with $\zeta=1,i,1,i,1,i,1,i,1$, respectively. This defines an exact
Gaussian-rational unitary. All 36 entries of $U_0$ are nonzero.

An exact phase-I calculation found $w^0_\kappa\ge10^{-24}$ for all 90 colorings,
with equality for some. The integer grouped-statistic matrix has 358 rows and
rank 43. Discovery required 124 rational pivots; it is not trusted for the
proof. The portable verifier independently computes **all 924 individual
squared minors**, groups them, and verifies all 358 identities directly from
the stored positive masses. It also recomputes the canonical masses and checks
that every grouped sum of the stored correction is exactly zero.

## 4. A common-minor obstruction throughout the same ball

Use the four singleton terms and nine two-minor terms in the archived
four-cycle cover, on indices $0,1,2,3$:

| Rows | Columns |
|---|---|
| $\{0\}$ | $\{0\}$ |
| $\{1\}$ | $\{3\}$ |
| $\{2\}$ | $\{1\}$ |
| $\{3\}$ | $\{2\}$ |
| $\{0,1\}$ | $\{0,2\},\{0,3\},\{1,3\}$ |
| $\{0,2\}$ | $\{0,1\},\{0,3\},\{1,2\}$ |
| $\{0,3\}$ | $\{0,1\},\{0,2\},\{2,3\}$ |

Add the eight singleton terms with row in $\{4,5\}$ and column in
$\{0,1,2,3\}$. Call the sum of these 21 squared minors $H(U)$.

For any permutation of six coordinates, either one of the first four columns
escapes to row 4 or 5, giving an added term equal to one, or those four columns
stay in the first four rows and the original cover is at least one. Thus every
permutation satisfies $H(P)\ge1$. The checker also enumerates all 720
permutations, with minimum exactly one. Any full common-minor distribution
would force $H(U)\ge1$.

Direct exact evaluation gives $H(U_0)=1-\Gamma<1$. For two unitary matrices in
entrywise distance at most $\delta$, a singleton squared modulus changes by at
most $2\delta$. A two-minor squared modulus changes by at most $8\delta$,
using the previous determinant estimate and both moduli at most one.
There are twelve singleton terms and nine two-minor terms, so

$$|H(U)-H(U_0)|\le96\delta.$$

The verifier checks $\delta_0\le\Gamma/192$. Hence
$H(U)\le1-\Gamma/2<1$ throughout the same closed ball. The distinction between
pattern certificates and full common-minor certificates is therefore stable
on a genuine neighborhood, not confined to the single archived matrix.

## 5. Consequences in arbitrary dimension

For any pattern $(r_0,r_1,r_2)$ with all $r_j\ge2$, take a unitary in the proved
six-dimensional ball and append $r_j-2$ fixed coordinates of color $j$.
Extend each coloring by these fixed colors. Principal-minor identities factor
across the direct sum, so this is a pattern certificate in dimension
$r_0+r_1+r_2$. It still cannot have a full common-minor certificate: zero
cross-block singleton marginals would force a hypothetical permutation law to
preserve the six-dimensional block almost surely, contradicting the core
obstruction. This padded class is structured; it is not claimed to be open
in the full higher-dimensional unitary group.

There is also a useful non-block-diagonal closure operation. Given any pattern
certificate for $U$, append a new identity coordinate of a chosen color, and
**left-multiply** by an arbitrary complex two-coordinate unitary acting on
that new coordinate and one old coordinate. Extend the coloring law by the
fixed new color, then independently swap those two coordinates with probability
equal to the squared modulus of the gate's off-diagonal entry.

To prove this rule, use Cauchy–Binet. The gate support and the old matrix support
meet in at most one coordinate, so for each pair of endpoint subsets there is
at most one nonzero intermediate-subset term, exactly as in the
[one-coordinate gluing proof](../single-cycle/README.md#3-one-coordinate-gluing-theorem).
After taking squared moduli and grouping the input columns by color counts,

$$G_{\text{gate}(U\oplus1)}(S,\nu)
 =\sum_R |\det(\text{gate})[S,R]|^2G_{U\oplus1}(R,\nu).$$

The random coordinate swap has precisely that subset transition law, so it
preserves every required grouped statistic. Iteration gives connected
structured examples in arbitrarily large dimension with all multiplicities
at least two. The operation is stated only for a **new coordinate and left
multiplication**; arbitrary gate interleavings or right multiplication across
color blocks are not established. The full-common-minor obstruction was proved
uniformly for the six-dimensional ball and its direct-sum padding; it is not
asserted after unrestricted pendant mixing.

Finally, Cauchy–Binet converts any of these pattern certificates to

$$\det(D_a+U\operatorname{diag}(\beta_0,\beta_0,\beta_1,\beta_1,\beta_2,\beta_2)U^*)
 =\sum_\kappa w_\kappa(U)\prod_i(a_i+\beta_{\kappa(i)}).$$

The analogous count-adjusted formula holds in larger dimension. These are
original permutation hull points with repeated spectral values, so this proves
MDO on the stated classes for every complex $a$ and $\beta$.

## Replay, provenance and limits

Python 3.9+ and the standard library only; local validation used Python 3.9.6.
From the repository root:

```sh
python3 -B research/marcus-de-oliveira/three-color-coupling/coupling.py
python3 -B -m unittest discover -s research/marcus-de-oliveira/three-color-coupling -p 'test_*.py' -v
```

The solver-free output matches [verification.json](verification.json).
Nine tests include changed/omitted masses, malformed values, incorrect
corrections, nonunitary/out-of-ball inputs, a distinct exact complex neighbor
with changed weights, direct grouped-minor reconstruction there, and an
original determinant check with three noncollinear complex spectral values.
The pendant rule is also independently replayed on a genuinely complex
seven-dimensional example with pattern $(3,2,2)$, checking all 3,432 squared
minors against the constructed coloring law.
`corrected_weights(U)` evaluates the explicit corrected law at an exact
Gaussian-rational unitary in the ball. The mathematical theorem also includes
irrational complex unitaries.

Optional discovery reconstruction, unnecessary for verification:

```sh
python3 -B research/marcus-de-oliveira/three-color-coupling/discover.py --output /tmp/three-color-rederived.json
```

This bounded exact phase-I program reproduces all stored masses and corrections.
Its output was checked against `interior.json`. Discovery used explicit fixed
gates and exact linear algebra, not another random LP campaign. The initial
floor $10^{-30}$ was increased to the final $10^{-24}$ and the certificate was
recomputed. The exact solver-free checker, not an optimizer success flag, is
the evidence for feasibility. No paid compute or fabricated usage is involved.

Attribution: the pattern formulation, canonical signed candidate, and
Conjecture F come from the existing [minors report, §§2–4](../minors/report.md).
The obstruction starts from the existing four-cycle cover, and the pendant
argument uses the existing unique-intermediate-subset lemma. The additions here
are the explicit strict interior certificate, its polynomial correction law
with quantitative positivity, and the simultaneous open-neighborhood
separation from every full common-minor certificate. All tests are author-side
checks under the same GitHub owner; independent-owner review remains pending.
The unrestricted $(2,2,2)$ pattern and Conjecture F in general remain open here.
