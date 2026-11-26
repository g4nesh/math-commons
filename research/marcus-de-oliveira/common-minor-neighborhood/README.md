# A quantitative neighborhood of the common-minor locus in every dimension

Fix $n\ge3$ and two labeled complex spectra $a,b$. For a permutation taking columns to rows put $v_\pi=\prod_j(a_{\pi(j)}+b_j)$, and let $K$ be their planar convex hull. Assume $K$ has positive area and each edge contains exactly two permutation **labels**, its endpoints. Interior coincidences do not affect the proof.

**Theorem.** There is an explicit $\xi>0$ depending only on the spectra such that a unitary $U$ satisfies the Marcus–de Oliveira inclusion whenever a distribution on permutations approximates every squared minor of $U$ to absolute accuracy $\xi$. It follows that the inclusion holds in a uniform open neighborhood of the entire common-minor locus in $U(n)$. The edge hypothesis holds on an open dense set of spectral pairs. The neighborhood includes unitaries that admit no exact common-minor distribution.

Write $A=\operatorname{diag}(a)$, $B=\operatorname{diag}(b)$, $D(U)=\det(A+UBU^*)$. Define
$$
E(U)=\min_{w\ge0,\;\sum w=1}\max_{|S|=|T|}
\left|\,|\det U[S,T]|^2-\sum_{\pi(T)=S}w_\pi\right|.
$$
This minimum over a finite simplex is attained. For each edge choose its outward **unit** support functional $\ell:\mathbb C\to\mathbb R$ and maximum $h$. Let $g$ be the minimum of $h-\ell(v_\nu)$ over all edges and all labels other than their two endpoints. Any positive lower bound for this minimum can be used below. Take $R\ge\max(1,\max|a|+\max|b|)$ and set
$$
\rho=\min\left(1,\frac{3g}{16n^3R^n}\right),\quad
\eta=\min\left(\frac1{16},\frac{\rho^2}{512n}\right),\quad
\xi=\min\left(\frac{\eta}{n},\frac{g\eta}{2\binom{2n}{n}R^n}\right).
$$
Then $E(U)\le\xi$ implies $D(U)\in K$. In particular every violation for these spectra requires $E(U)>\xi$. These deliberately conservative constants give a test using only the finite permutation products and minors. No unknown global extremum of the determinantal range enters them.

## Which spectra qualify

Each permutation product is a holomorphic polynomial. For three distinct permutations the ratio $(v_3-v_1)/(v_2-v_1)$ cannot be constant. Indeed the coefficients that select one $b_j$ and all $a$'s except $a_i$ recover the entries of the corresponding permutation matrix. A constant ratio would make three distinct permutation matrices collinear. With real entries that would make one a convex combination of the other two, contrary to extremality. The same coefficient test makes the products pairwise distinct as polynomials. Hence for generic complex spectra no two values coincide and no three are collinear: outside the coincidences each nonconstant holomorphic ratio cannot take only real values on an open set. There are finitely many triples. The edge condition itself is open by finiteness and strict edge gaps.

## Support near an endpoint block

Fix an edge, relabel the columns to put one endpoint at the identity, and call the other permutation $\tau$. Let $F(U)=h-\ell(D(U))$. If $\tau$ is a transposition, let $H$ be the block subgroup with one two-dimensional block on its moved coordinates and scalar blocks elsewhere. Otherwise let $H$ be the diagonal unitary torus. In the latter case we use the analogous torus at the other permutation as well, after relabeling. For every $V\in H$, $F(V)=0$. In the two-site case the determinant is exactly the convex combination of the two endpoint products with the squared stay and swap amplitudes as weights.

Use the Frobenius metric on skew-Hermitian tangent vectors. In $U=\exp(X)V$, normal directions to $H$ have zero diagonal blocks. Write $X_{ij}$ for the rectangular block connecting two different blocks, with its negative adjoint in the opposite position. Conjugating by arbitrary block-scalar phases fixes $V$ and preserves $F$ everywhere. The normal derivative therefore vanishes there. Mixed Hessian terms for different unordered block pairs vanish: phase averaging multiplies those pairs by distinct nonconjugate characters, up to orientation. Equivalently average over independent signs on the blocks and flip the sign on a block incident to exactly one of the pairs.

For a single off-diagonal block write $s=\|X_{ij}\|_F$. Along the one-block rotation $\exp(tX)V$ the following estimates hold.

- For two singleton blocks the swap probability is $\sin^2(ts)$. Any disjoint original two-dimensional block has its own stay/swap distribution. The product distribution assigns at least $\sin^2(ts)$ total mass to nonendpoint permutations, since swapping the two singletons is never either endpoint.
- For a singleton $k$ coupled to the two-dimensional block, the varying part is a three-dimensional unitary. Every three-dimensional unitary admits a common-minor law: decompose the doubly stochastic singleton array by Birkhoff, and use complementary minors for the order-two entries. The determinant expansion then gives the corresponding probability law on products. The probability of moving $k$ is $1-|U_{kk}|^2=\sin^2(ts)$ since $X_{ij}$ is a single complex vector. Both endpoint permutations fix $k$.
- At a diagonal torus used for a nontransposition edge, a single swap also always gives a nonendpoint.

In every case $F(\exp(tX)V)\ge g\sin^2(ts)$. Summing the second variations with the mixed terms absent gives $F''(0)\ge g\|X\|_F^2$. This estimate holds uniformly over the full subgroup, including the entire one-block rotation between two tied permutations.

The third derivative of $D(\exp(tX)V)$ has absolute value at most $(2n\|X\|_{\rm op})^3R^n$. In the determinant differentiated by its columns there are $n^3$ ordered choices of differentiated column. An undifferentiated column has norm at most $R$ and an order-$k$ derivative has norm at most $R(2\|X\|_{\rm op})^k$ by the iterated commutator. Hadamard's inequality bounds each term. Taylor's theorem now gives
$$
F(\exp(X)V)\;\ge\;\tfrac14 g\|X\|_F^2\qquad (X\perp\mathfrak h,\ \|X\|_F\le\rho).
$$
The geodesic distance to $H$ equals the length of a shortest normal geodesic: $H$ is compact, the metric is bi-invariant, and the first variation at a minimizing point is perpendicular to $H$. For two unitaries the geodesic distance is at most $(\pi/2)$ times their Frobenius chord distance, by diagonalizing their relative unitary and applying the scalar chord inequality on $[-\pi,\pi]$.

## Turn an approximate law into a support bound

Take a witnessing distribution with residual at most $\xi$. Write $p,q$ for its masses on the two labels of the chosen edge and $\varepsilon=1-p-q$. The diagonal expansion and Cauchy–Binet express $D(U)$ as a sum over equal-cardinality subsets $S,T$, with coefficient $|\det U[S,T]|^2$ times $\prod_{i\notin S}a_i\prod_{j\in T}b_j$. Expanding each pairing product gives the same monomials with permutation probabilities in place of these coefficients. Thus
$$
|D(U)-\sum w_\nu v_\nu|\le\binom{2n}{n}R^n\xi .
$$
If $\varepsilon\ge\eta$, the weighted product lies at least $g\eta$ below this support, and the error is at most $g\eta/2$.

Suppose instead $\varepsilon<\eta$. Put $e=\varepsilon+n\xi\le2\eta\le1/8$. Each column $j$ has squared mass outside $\{j,\tau(j)\}$ at most $e$, by summing its singleton errors. Its two designated entries have squared moduli at least $p-\xi$ and $q-\xi$, respectively.

1. If $\tau$ is a transposition, the total squared norm in the off-diagonal blocks for $H$ is at most $ne$. Polar decomposition on each diagonal block produces a block unitary within Frobenius distance $\sqrt{2ne}$ of $U$. Indeed the block singular values $d$ are at most one and $(1-d)^2\le1-d^2$; the sum of the latter quantities equals the total squared off-block norm.
2. If $\tau$ has a cycle of length at least three, take $k=\tau(j)$ in this cycle. The permitted row sets for columns $j,k$ intersect just at $k$. In their orthogonality identity every other row is outside at least one permitted set. Cauchy–Schwarz on the two resulting sums bounds $|\bar U_{kj}U_{kk}|\le2\sqrt e$.
3. Otherwise a nontransposition has at least two disjoint two-cycles. Take $j,k$ in distinct cycles. Neither endpoint maps columns $\{j,k\}$ to rows $\{j,\tau(k)\}$. Thus this two-minor has modulus at most $\sqrt{\varepsilon+\xi}\le\sqrt e$. The off-diagonal cross term in its determinant has modulus at most $e$ by singleton leakage. This bounds $|U_{jj}U_{\tau(k),k}|\le\sqrt e+e\le2\sqrt e$.

In cases 2 and 3 these inequalities give $(p-\xi)_+(q-\xi)_+\le4e$. The larger factor is at least $(1-\varepsilon)/2-\xi\ge1/4$. Hence the smaller endpoint mass is at most $16e+\xi$. Every column then has an entry of squared modulus at least $1-19e$ in the larger permutation's position. Choosing its phase columnwise gives a monomial unitary within distance $\sqrt{38ne}$: the squared chord to the phase is $2-2|U_{ij}|\le2(1-|U_{ij}|^2)$. The bound remains valid when $1-19e$ is negative.

In all cases the geodesic distance to the corresponding endpoint subgroup is at most $2\sqrt{38ne}\le2\sqrt{76n\eta}\le\rho$. The normal support estimate proves $F(U)\ge0$. Taking all edges proves the theorem. In the limiting exact equality case the same classification shows a common-minor unitary on the supporting line must be monomial or lie in the single active two-site block, up to the label changes.

For $V$ with an exact common law, the same law approximates $U$ with residual at most $2n\|U-V\|_{\rm op}$. The determinant telescope on each size-$k$ compression bounds the amplitude change by $k\|U-V\|_{\rm op}$, and both amplitudes have modulus at most one. This proves the asserted neighborhood with operator radius $\xi/(2n)$. The bounds persist locally in the spectra because edge identities are stable and $g,R$ have positive and finite uniform bounds on a sufficiently small closed spectral neighborhood. In particular the same radius applies around every acyclic two-site gate circuit with a common-minor law, even for perturbations that add cycles. The threshold may tend to zero as spectra approach a degenerate hull.

## An exact irreducible example outside the common-minor locus

Use $a=(-4-4i,3i,2-2i,-3-2i,4-3i)$ and $b=(-1+i,2+5i,4,5-5i,2-i)$. The 120 permutation products have seven hull edges. An exact unit-normal gap lower bound is $37438/379$. The [standard-library verifier](check.py) recomputes the threshold, checks a normal Hessian on an active two-site edge with a complex phase, and directly verifies the final determinant against all seven edges.

Let $R_{ij}(t)$ be the rotation with cosine $(1-t^2)/(1+t^2)$ and sine $2t/(1+t^2)$. On four sites set $W=R_{03}R_{23}R_{12}R_{01}$ with this same parameter. Any permutation law must satisfy
$$
m_{02,23}-m_{0,2}-m_{3,1}-m_{01,13}\le0,
$$
where $m_{S,T}$ is its subset probability and digits denote increasing subsets. If columns $23$ map to $02$, either column 2 maps to 0 or column 3 maps to 0. In the second case column 1 maps to 3 or to 1, accounting for the other two terms. For a fifth site subtract also its probability of moving. This gives a valid inequality for all 120 permutations. Direct two-minor expansion on $W$ gives the left side
$$
\frac{16t^4(t^4-6t^2+1)}{(1+t^2)^6}.
$$
The verifier also proves this rational identity by evaluation at 33 distinct rational parameters. After clearing $(1+t^2)^{16}$, both sides have degree at most 32 because each matrix entry comes from four gates and the largest squared minor has order two. The expression is positive for arbitrarily small positive $t$. The fourth-order violation illustrates why exact minor feasibility is too restrictive even very close to a commuting basis.

Take $t=\min(1/128,\xi/(24n))$ for the computed $n=5$ threshold. Append an identity site and then a complex gate between sites 3 and 4 with parameter $t^5$. The resulting matrix has connected bipartite support, so for the two simple spectra it has no common reducing block. Each gate has operator distance at most twice its parameter from the identity, giving residual at most $2n(8t+2t^5)\le\xi$. Exact arithmetic checks that the extended separating inequality is still violated, although the original determinant lies strictly inside the hull. Run `python3 -I -B research/marcus-de-oliveira/common-minor-neighborhood/check.py`.

## Relation to earlier results

Kovačec, Bebiano and da Providência studied the local shape of determinantal ranges for real orthogonal matrices at permutation points under a fiber hypothesis. Mulero-Martínez studied support functionals on the unitary group. Our statement concerns an explicit LP-residual tolerance on the whole complex unitary group for fixed regular spectra. The classical local geometry, complementary-minor identities and single-transposition segment are used in this setting. No priority claim is made. [On the corners of certain determinantal ranges](https://www.sciencedirect.com/science/article/pii/S0024379507001644) · [A variational framework for determinantal inequalities of normal matrices](https://www.sciencedirect.com/science/article/pii/S0024379526001266)

This is a written proof with exact finite replays. Independent mathematical review is requested. Spectra with larger sets of supporting labels and unitaries with larger residuals remain outside this theorem.
