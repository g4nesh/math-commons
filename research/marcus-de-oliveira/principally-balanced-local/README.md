# A local principal-minor law around every simple balanced Hermitian center

**Theorem.** Let $H_0$ be an $n\times n$ Hermitian matrix with distinct eigenvalues, and suppose that all its principal minors of each fixed order are equal. There is an open neighborhood $\mathcal N$ of $H_0$ in the real vector space of **all Hermitian matrices** and continuous strictly positive probabilities $w_p(H)$, $p\in S_n$, such that

$$\det H[S,S]=\sum_{p\in S_n}w_p(H)\prod_{i\in S}\lambda_{p(i)}(H)\qquad(S\subseteq[n],\ H\in\mathcal N).$$

Here the real eigenvalues are ordered continuously near their distinct values at $H_0$. The same law serves all subset orders simultaneously. Consequently, for **every** $a\in\mathbb C^n$,

$$\det(\operatorname{diag}(a)+H)=\sum_p w_p(H)\prod_i(a_i+\lambda_{p(i)}(H)).$$

The weights depend on $H$ and its spectrum, but not on the arbitrary opposing spectrum $a$. Both eigenvectors and eigenvalues of $H$ may change. This does not assert a law independent of both spectra at a fixed relative eigenbasis, nor a global theorem for every Hermitian matrix.

A concrete center exists in every dimension:

$$ (H_n)_{jk}=\begin{cases}i&j<k,\\0&j=k,\\-i&j>k.\end{cases} $$

The argument below proves the theorem for any simple principally balanced Hermitian center. The term *principally balanced* and this skew-sign family are established literature; the local argument is presented with a proof and no priority claim. The evidence is a written proof with an exact finite replay, not a Lean theorem or an independent-owner review.

## 1. The only affine relations for a distinct real spectrum

For fixed distinct real $\lambda_1,\ldots,\lambda_n$, define

$$v_p(\lambda)_S=\prod_{i\in S}\lambda_{p(i)},\qquad S\subseteq[n].$$

In this note, **$p$ sends a coordinate position to an eigenvalue label**. This is the inverse of the column-to-row convention used in some other reports; the entire permutation set and convex hull are unchanged. Include the empty subset coordinate, with value one. For every $k$,

$$\sum_{|S|=k}v_p(\lambda)_S=e_k(\lambda).$$

We prove these are all the affine relations. The resulting affine dimension is

$$d=2^n-n-1.$$

**Injective-grid lemma.** A multiaffine polynomial in $m$ variables that vanishes on every injective $m$-tuple from a set of at least $m+1$ distinct real numbers is identically zero.

*Proof.* Induct on $m$, with the constant case $m=0$ immediate. Write the polynomial as $f_0+x_m f_1$. Fix any injective tuple in the first $m-1$ slots. At least two distinct choices for the last slot remain. Vanishing at both makes both coefficients zero at the fixed tuple. Apply induction to $f_0,f_1$, using the same set of available numbers. This proves every coefficient is zero. No assertion about arbitrary-degree polynomials is being made.

**Affine-span lemma.** If a multiaffine polynomial $F(x_1,\ldots,x_n)$ is constant on the permutation orbit of a distinct $\lambda$, then $F$ is symmetric.

*Proof.* For any two coordinates $i,j$, multiaffinity gives

$$F(x)-F((ij)x)=(x_i-x_j)G(x_{[n]\setminus\{i,j\}}),$$

where $G$ is multiaffine and independent of $x_i,x_j$. For each injective assignment of the other $n-2$ variables from the $n$ distinct spectral values, use the two remaining values for $x_i,x_j$. Constancy of $F$ on the orbit and their inequality force $G$ to vanish. The injective-grid lemma proves $G=0$. Thus $F$ is invariant under every transposition. The cases $n=1$ and constant $G$ are included separately or directly.

Symmetric multiaffine polynomials have one coefficient per subset size, so they form the $(n+1)$-dimensional span of $e_0,\ldots,e_n$. Linear functionals constant on the vectors $v_p$ are exactly such polynomials. The annihilator of their difference span therefore has dimension $n+1$, proving the dimension claim. In particular that difference span is the fixed space

$$W=\{z\in\mathbb R^{2^n}:\ \sum_{|S|=k}z_S=0\text{ for }k=0,\ldots,n\}.$$

This establishes full affine span at **every** distinct spectrum, including spectra containing zero or having vanishing elementary symmetric functions. It is not merely a generic-rank calculation.

## 2. An explicit continuous correction of the uniform law

Let $N=n!$ and

$$b(\lambda)=\frac1N\sum_p v_p(\lambda),\qquad b(\lambda)_S=\frac{e_{|S|}(\lambda)}{\binom n{|S|}}.$$

For each degree choose one reference subset; delete those $n+1$ coordinates. The resulting projection $\pi:W\to\mathbb R^d$ is an isomorphism, since the deleted coordinate is the negative sum of the retained coordinates of that degree.

Fix a reference permutation $p_0$. At the center spectrum $\lambda^0$, the span lemma supplies $d$ other permutations $p_1,\ldots,p_d$ for which the square matrix

$$T(\lambda)=\big[\pi(v_{p_j}(\lambda)-v_{p_0}(\lambda))\big]_{j=1}^d$$

is invertible at $\lambda^0$. These permutations stay fixed. Its entries are polynomials in the spectral values, so $T(\lambda)$ remains invertible nearby, with continuous inverse.

Let $m(H)_S=\det H[S,S]$ and

$$e(H)=m(H)-b(\lambda(H)).$$

All principal minors of a Hermitian matrix are real. The characteristic-polynomial coefficient formula gives $\sum_{|S|=k}m(H)_S=e_k(\lambda(H))$ for every $k$, hence $e(H)\in W$ exactly, even when the spectrum moves. Define

$$c(H)=T(\lambda(H))^{-1}\pi e(H).$$

Set

$$w_{p_j}=\frac1N+c_j\quad(1\le j\le d),\qquad
w_{p_0}=\frac1N-\sum_{j=1}^d c_j,$$

and set every other weight to $1/N$. They sum to one, and their coefficient vector is $b+e=m(H)$ because $\pi$ is injective on $W$.

At a principally balanced center, the characteristic-polynomial formula says precisely $m(H_0)=b(\lambda^0)$. Thus $c(H_0)=0$. Ordered eigenvalues of Hermitian matrices are continuous, and simplicity keeps their labels distinct locally. Principal minors and the matrix inverse above are continuous too. Therefore the weights remain strictly positive on a full open Hermitian neighborhood. For example it suffices to impose

$$\|c(H)\|_1<\frac1{n!}.$$

One may shrink the neighborhood so that every weight exceeds $1/(2n!)$. For $n=1$, $d=0$, the sole weight is one and the result holds everywhere; no nonempty inverse is required.

**A computable certificate.** Given a candidate $H$ with distinct real eigenvalues, any chosen invertible $T(\lambda)$ and the exact inequality above certify the law, regardless of proximity to a balanced center. Alternatively, with $T_0=T(\lambda^0)$, the two bounds

$$\|T_0^{-1}(T(\lambda)-T_0)\|_1\le\tfrac12,
\qquad \|T_0^{-1}\pi e(H)\|_1<\frac1{2n!}$$

suffice by the Neumann-series inverse bound. The matrix norm is the induced column-sum norm, and the vector norm is the sum of absolute values. These are finite tests suitable for exact algebraic arithmetic or rigorous interval bounds. This report does not claim a uniform numerical radius independent of the center or dimension.

Finally, expanding $\det(\operatorname{diag}(a)+H)$ by complementary principal minors and substituting the coefficient identities gives the determinant formula in the theorem. Since these are polynomial identities in all complex $a_i$, the conclusion is uniform over unbounded opposing spectra, not just a compact set of them.

## 3. A principally balanced simple center in every dimension

Every increasing-index principal submatrix of $H_n=iK_n$ is exactly $H_k$, where $K_n$ has ones above and minus ones below the diagonal. Thus $H_n$ is principally balanced and Hermitian. This is the skew-symmetric family already appearing in Huang–Oeding, Theorem 1.1(2), multiplied by $i$.

To prove simplicity directly, let $H_nv=\lambda v$ with real $\lambda$. Subtract the equations in consecutive rows. The resulting relation is

$$(\lambda+i)v_{j+1}=(\lambda-i)v_j.$$

Since $\lambda+i\ne0$, the entire eigenvector is determined by its first coordinate; that coordinate cannot vanish for a nonzero eigenvector. Each eigenspace is therefore at most one-dimensional. Hermitian diagonalizability then proves every eigenvalue is simple. This works also at $\lambda=0$ and for $n=1$.

For an optional closed form, write $r=(\lambda-i)/(\lambda+i)$, so $v_j=r^{j-1}v_1$. The first-row equation together with the geometric sum gives $r^n=-1$. Conversely each root gives the same recurrence and first-row equation. The eigenvalues are therefore

$$\lambda_k=\cot\frac{(2k-1)\pi}{2n},\qquad k=1,\ldots,n.$$

The proof of simplicity does not require this formula or a numerical eigenvalue computation.

## 4. Exact off-center example and replay

The [exact checker](check.py) constructs a non-balanced Hermitian perturbation of $H_4$ by a rational orthogonal rotation of its first two coordinates, with cosine $999999/1000001$ and sine $2000/1000001$. It preserves the distinct spectrum $\{1+\sqrt2,\sqrt2-1,1-\sqrt2,-1-\sqrt2\}$ but changes its order-two principal minors. It constructs the fixed-basis correction in $\mathbb Q(\sqrt2)$, verifies every subset identity and positivity, and exports all 24 probabilities in [witness.json](witness.json), with the coefficient convention stated above.

The correction has exact norm

$$\|c\|_1=\frac{2999997000}{1000002000001}\sqrt2<\frac1{48},$$

and all 24 weights exceed $1/48$. The checker proves these inequalities by rational comparisons in the ordered quadratic field; no floating-point value is used. It also checks the ranks $0,1,4,11,26$ at distinct integer spectra of sizes one through five. Those finite checks support the implementation; the unbounded rank result is the proof in Section 1.

Run from the repository root with Python's standard library:

```sh
python3 -B research/marcus-de-oliveira/principally-balanced-local/check.py
python3 -B -m unittest discover -s research/marcus-de-oliveira/principally-balanced-local -p 'test_*.py'
```

The first command regenerates the exact witness in memory and requires equality with the archived JSON. The optional `--write` flag intentionally rewrites that JSON only after every exact check passes. This source-only contribution does not modify shared CI. An exact off-center witness is not a numerical lower bound on the radius of the whole Hermitian neighborhood.

## 5. Literature and remaining scope

- Huajun Huang and Luke Oeding, [*Symmetrization of Principal Minors and Cycle-Sums*, arXiv:1510.02515v2](https://arxiv.org/abs/1510.02515v2), Theorem 1.1(2), classify the skew-symmetric balanced family containing $K_n$. The family itself is credited to that literature; it is not introduced here.
- Octavio Arizmendi, Franz Lehner and Amnon Rosenmann, [*On matrices in finite free position*, Linear Algebra and its Applications 724 (2025), 137–170](https://doi.org/10.1016/j.laa.2025.06.016), Definition 4.2, Lemma 4.5 and Theorem 4.7, give the principally balanced terminology, the equal-minor characteristic-polynomial coefficients, and the diagonal/principally-balanced finite-free pairing. The exact center identity is thus established prior work.
- Yuki Ueda, [*Finite free position of maximal abelian star-subalgebras of the matrix algebra*, arXiv:2609.39677v1](https://arxiv.org/abs/2609.39677v1), studies the stronger requirement that one fixed relative eigenbasis realize uniform finite-free averaging for **every** pair of diagonal inputs. Its dimension restriction concerns that different quantifier structure; here the positive law is allowed to vary with the Hermitian matrix and its spectrum.

A bounded primary-source search did not locate this precise full-neighborhood formulation in the cited papers. That is not a novelty determination. The neighborhood theorem above is proved directly, with the prior center result credited. No unrestricted principal-minor Schur–Horn theorem or Marcus–de Oliveira conjecture is claimed. Repeated-spectrum centers require a separate argument because the full-span lemma can fail there.
