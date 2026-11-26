# Exact pole compression and interlacing lifts

Shitov proved MDO for Hermitian $A$ when at most three eigenvalues of $B$ are nonreal, by Möbius inversion, interlacing and dimension lifting. We give a direct finite determinant factorization at the pole, then apply the interlacing argument using a computable smaller spectrum. This allows the four-dimensional cross-ratio result in the repository to be used in higher dimensions. The algebra below is independent of that computational result, and the four-exception corollary carries its review status. [Shitov, Theorem 19 and Lemma 17](https://www.mathnet.ru/eng/im9292) · [Merged four-dimensional result](../star-halfline/README.md)

## The exact deletion identity

Work in an eigenbasis of the normal matrix $B$. Let $b$ be an eigenvalue on a coordinate set $E$, take $E$ to be its entire multiplicity, and let $R=E^c$. Suppose $A+bI$ is invertible. Put
$$G=(A+bI)^{-1},\quad C=G[R,R],\quad
\gamma_i=(a_i+b)^{-1},\quad \beta_j=(b_j-b)^{-1}\ (j\in R),\qquad
S=\det(A+bI)\prod_{j\in R}(b_j-b).$$
Then, identically,
$$\det(A+B)=S\det(C+D_\beta). \tag{1}$$
Indeed factor $A+bI$ on the left. The columns on $E$ of $I+G(B-bI)$ are identity columns. Expanding them leaves $I_R+C\operatorname{diag}(b_j-b)$, which gives (1). The pairing vertices are **exactly** the numbers
$$S\prod_{j\in R}(\gamma_{\iota(j)}+\beta_j),\qquad \iota:R\hookrightarrow[n]. \tag{2}$$
Every injection extends by pairing the leftover eigenvalues to the repeated $b$; conversely every full pairing restricts to an injection. Thus global hull membership is equivalent to membership of $\det(C+D_\beta)$ in the hull of the unscaled injective products. Equation (1) by itself assumes no normality of $C$. It even holds for nonnormal $A$ when the indicated inverse exists; the MDO conclusions concern normal pairs.

## Lifting an ordinary lower-dimensional certificate

If $A$ is Hermitian and $b$ is real, $G$ and its principal compression $C$ are Hermitian. Let $\alpha$ be the spectrum of $C$. Each permutation vertex of the ordinary hull $\Delta(\alpha,\beta)$ is a convex combination of the injective products in (2), by interlacing. Here is an explicit common interpolation construction.

For deletion of one coordinate, order the eigenvalues $g_1\le\cdots\le g_n$ and $g_i\le h_i\le g_{i+1}$. Construct a random injection for $(h_i)$ from the largest index downward. Assign $h_{n-1}$ either $g_{n-1}$ or $g_n$ with the ordinary interpolation probabilities giving mean $h_{n-1}$. Delete the selected $g$. The remaining $h_1,\dots,h_{n-2}$ still interlace the remaining list in either branch. Repeat this rule. Equal interpolation endpoints give the same value; select either index deterministically. Induction, conditional on the largest assignment, gives
$$\mathbb E\prod_i(g_{\iota(i)}+t_i)=\prod_i(h_i+t_i)$$
for all complex $t_i$, with the same injection law. For several deletions apply Cauchy interlacing along a chain of consecutive principal compressions. Sample successive injections independently and compose them. A subset of the intervening factors can be omitted, since its product identity is a coefficient of the full multi-affine identity. This justifies composition even when an intermediate injection uses fewer than all coordinates.

Consequently any ordinary lower-dimensional certificate for $(C,D_\beta)$ supplies a certificate upstairs. All inversions and permutations in this argument are finite. Removing distinct real eigenvalues individually preserves the Hermitian/real hypotheses on the remaining real coordinates. In generic nonsingular data one can therefore delete real coordinates down to dimension three, recovering Shitov's known conclusion; the full known theorem also covers degenerate inputs by closedness. For the local four-dimensional criterion below we state the invertibility assumptions directly. A singular datum is also covered if it is obtained as a limit of covered nonsingular data.

## A four-exception family in every dimension

Fix $A=A^*$, a real $b$ with $A+bI$ invertible, and a four-dimensional subspace $V$ with orthogonal projection $P$. Let $B$ equal $bI$ on $V^\perp$ and take **any normal operator on $V$**. Compute the four eigenvalues $h_1<h_2<h_3<h_4$ of $P(A+bI)^{-1}|_V$. If
$$\frac{(h_2-h_1)(h_4-h_3)}{(h_3-h_1)(h_4-h_2)}\in[1/6,5/6],$$
then MDO for $(A,B)$ follows from Theorem OS in the [four-dimensional report](../star-halfline/README.md). It applies to every orientation and every spectrum of the four exceptional eigenvalues. If some equal $b$, perturb them within their eigenbasis and take a limit. Orientation of $V$ inside the original space is arbitrary subject to the displayed compression test. No condition on the cross-ratio of the **original** $A$ spectrum is substituted. The higher-dimensional implication is exact even if Theorem OS is treated as a conditional premise pending independent mathematical review.

**Circular variant.** The single-pole result needs only collinearity of the eigenvalues of $G$. Then $G=\eta I+\omega H$ with $H$ Hermitian and $|\omega|=1$. Compression preserves this form; use interlacing in the real coordinate of the line. Every conditional interpolation identity transports under $x\mapsto\eta+\omega x$, and the four-point cross-ratio is affine invariant. In particular this covers any normal $A$ whose spectrum lies on a generalized circle (circle or line) through $-b$ and avoids the pole, since $a\mapsto1/(a+b)$ sends that generalized circle to a line. The outside eigenvalue $b$ may now be complex. The four arbitrary exceptional eigenvalues of $B$ need not be on a generalized circle.

For distinct real exterior eigenvalues one can first delete them successively by (1), checking the inverses at each step, and apply the same cross-ratio test to the final Hermitian matrix. The reduced eigenbasis can vary. Elimination depends on the eliminated real eigenvalues and their eigenspaces. Their uneliminated complex eigenvalues are free as long as the indicated scalar reciprocals exist. This construction does not cover arbitrary compressed spectra whose ratio is outside the certified interval.

## Checking the cross-ratio with interval bounds

Suppose $C$ has one real eigenvalue in each of four ordered, disjoint intervals $[u_i,v_i]$. Then the following deliberately conservative inequalities ensure the cross-ratio condition:
$$\frac{(u_2-v_1)(u_4-v_3)}{(v_3-u_1)(v_4-u_2)}\ge1/6,\qquad
\frac{(v_2-u_1)(v_4-u_3)}{(u_3-v_1)(u_4-v_2)}\le5/6.$$
The occupancy hypothesis on $C$ must be established for the chosen subspace, for example by perturbation from a known compression with separated simple eigenvalues. Interlacing alone generally permits collisions between two compressed eigenvalues. The bounds follow from monotonicity of each positive distance in its endpoints.

## Two exact order-six examples

The two fixtures have a Hermitian or circular original spectrum and four distinct nonreal exceptional eigenvalues. Their compressed spectra are respectively $(2,4,6,8)$ and $(2+i,4+i,6+i,8+i)$. Both cross-ratios equal $1/4$. The reduced eigenbasis is a flat order-four Householder, whose eight oriented three-cycles give a direct certificate. For each such pairing two independent binary interlacings lift it to four original pairings of mass $1/32$ apiece. The checker evaluates all 32 on the original spectra and independently computes both sides of (1). This instance-level test does not need the four-dimensional cross-ratio theorem. Run `python3 -B research/marcus-de-oliveira/pole-compression/check.py`. [Exact checker and full spectra](check.py)
