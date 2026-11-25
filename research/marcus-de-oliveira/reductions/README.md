# Exact reductions and a rational-witness route

Task: `mdo-reductions` · Author: `ganesh-reductions` (GitHub owner `g4nesh`).

**Evidence level:** mathematical derivations for review, not a proof of the general Marcus–de Oliveira conjecture (MOC), not a formal verification, and no novelty claim. Written November 24, 2026 (America/Phoenix). The [objective brief](../README.md) fixes the original conjecture. These reductions are intended to make subsequent proof attempts and counterexample searches precise.

## 1. The relative unitary is essential

Let normal matrices have spectral decompositions

\[
A=V D_a V^*,\qquad B=W D_b W^*,\qquad U=V^*W.
\]

Then

\[
\det(A+B)=\det(D_a+U D_b U^*)=:d_{a,b}(U).
\]

Conversely, every diagonal pair and every unitary \(U\) produces a normal pair \(A=D_a\), \(B=U D_b U^*\) with exactly those spectra, with multiplicity. Thus the diagonal-orbit formulation and the original normal-pair formulation are equivalent. This argument never asserts that the original pair commutes. It uses no inverse of \(A\), \(B\), or \(A+B\), so singular sums and zero or repeated eigenvalues are included.

For fixed spectra define

\[
z_\sigma(a,b)=\prod_{i=1}^n(a_i+b_{\sigma(i)}),\quad
H(a,b)=\operatorname{conv}\{z_\sigma:\sigma\in S_n\}.
\]

Each \(z_\sigma\) is attained by a permutation unitary. The image of the compact connected space \(U(n)\) under \(d\) is compact and connected. Neither property implies that this image is convex. The conjecture asks whether this image lies in \(H\); it does not identify the image with the hull.

## 2. Support functions: both directions, all boundary cases

For a complex number \(w\), put \(L_w(z)=\operatorname{Re}(\overline w z)\). For any finite nonempty set \(Z\subset\mathbb C\),

\[
x\in\operatorname{conv}Z
\quad\Longleftrightarrow\quad
L_w(x)\leq\max_{z\in Z}L_w(z)\ \text{for every }w\in\mathbb C.
\tag{1}
\]

The forward implication follows by applying the real linear map \(L_w\) to a convex combination. For the reverse, suppose \(x\notin H=\operatorname{conv}Z\). Compactness gives a nearest point \(y\in H\). Let \(w=x-y\ne0\). For any \(z\in H\), the segment \(y+t(z-y)\) belongs to \(H\). The right derivative at \(t=0\) of its squared distance to \(x\) is nonnegative, giving

\[
L_w(z-y)\leq0.
\]

But \(L_w(x-y)=|w|^2>0\). Therefore \(L_w(x)>\max_{z\in H}L_w(z)\), contradicting the right side of (1). The argument works for point and segment hulls as well as two-dimensional polygons.

After dividing a nonzero \(w\) by its modulus, (1) yields exactly

\[
\operatorname{Re}(e^{-i\theta}d_{a,b}(U))
\leq\max_{\sigma}\operatorname{Re}(e^{-i\theta}z_\sigma(a,b))
\quad\text{for all }\theta\in\mathbb R.
\tag{2}
\]

Finite angle samples generally do not imply (2). For one explicitly computed hull, checking all its facet inequalities is sufficient if the hull is two-dimensional; point and segment hulls need their own affine-span and endpoint checks. A modulus comparison against the farthest vertex loses directional information.

## 3. Symmetries and their exact scope

The following transformations preserve membership, and every stated equivalence has an inverse.

| Transformation | Effect on determinant and permutation points | Scope |
| --- | --- | --- |
| \(a_i\mapsto a_i+t,\ b_j\mapsto b_j-t\) | Both are unchanged. | Every complex \(t\); opposite shifts only. |
| \(a_i\mapsto c a_i,\ b_j\mapsto c b_j\) | Both multiply by \(c^n\). | Equivalence for \(c\ne0\); \(c=0\) collapses to the trivial zero pair. |
| Conjugate \(a,b,U\) entrywise | Both conjugate. | Reflects the hull in the real axis. |
| Swap \(a,b\), replace \(U\) by \(U^*\) | The determinant is unchanged; permutation labels invert. | Conjugate \(D_a+UD_bU^*\) by \(U^*\). |
| Relabel either spectrum and permute the corresponding rows or columns of \(U\) | Only the labeling changes. | Multiplicities are retained. |
| Replace \(U\) by \(LUR\) for diagonal unitary \(L,R\) | \(d\) is unchanged. | \(R\) commutes with \(D_b\), and \(L\) commutes with \(D_a\). |

In the last row, \(D_a+LURD_bR^*U^*L^*=L(D_a+UD_bU^*)L^*\). Choosing \(R=\operatorname{diag}((\det U)^{-1},1,\ldots,1)\) shows that ranging over \(SU(n)\) is sufficient. This phase freedom does not let us replace a general complex unitary by a real orthogonal matrix.

A useful numerical normalization first uses \(t=b_1\), giving \(b'_1=0\). Set

\[
R=\max\{\max_i|a_i+b_1|,\ \max_j|b_j-b_1|\}.
\]

If \(R=0\), all pairwise sums are zero and MOC holds for this pair. Otherwise scale both spectra by \(1/R\). This preserves membership, but rescales absolute numerical gaps by \(R^{-n}\). It does **not** preserve Gaussian-rational entries in general because \(R\) can be irrational; use a positive rational upper bound instead when rationality matters.

For the universal problem over **all spectra**, proving only the direction \(\theta=0\) in (2) is sufficient: for another angle use \(c=e^{-i\theta/n}\), which multiplies every determinant and vertex by \(e^{-i\theta}\). This does not reduce a fixed-spectra numerical membership check to one direction; its spectra would change.

## 4. Exact rational data suffice for finding a counterexample

Write \(\mathbb Q(i)=\{p+iq:p,q\in\mathbb Q\}\). A skew-Hermitian matrix satisfies \(K^*=-K\). Define

\[
C(K)=(I-K)(I+K)^{-1}.
\tag{3}
\]

**Lemma 1 (rational Cayley construction).** For every skew-Hermitian \(K\), \(I+K\) is invertible and \(C(K)\) is unitary. If \(K\) has entries in \(\mathbb Q(i)\), so does \(C(K)\).

**Proof.** If \((I+K)v=0\), then the real part of \(v^*(I+K)v\) is \(\|v\|^2=0\), since \(v^*Kv\) is purely imaginary. Thus \(v=0\). Taking adjoints in (3) gives

\[
C(K)^*=(I-K)^{-1}(I+K).
\]

All these matrices and their inverses commute because they are rational functions of \(K\). Multiplication gives \(C(K)^*C(K)=I\). Inverting a nonsingular matrix over the field \(\mathbb Q(i)\) stays in that field. ∎

**Lemma 2 (density).** Rational Cayley unitaries \(C(K)\), with skew-Hermitian \(K\in\mathbb Q(i)^{n\times n}\), are dense in \(U(n)\).

**Proof.** If a unitary \(U\) has no eigenvalue \(-1\), then

\[
K=(I-U)(I+U)^{-1}
\]

is skew-Hermitian and \(C(K)=U\). For example, the spectral theorem reduces the first assertion to the fact that \((1-\lambda)/(1+\lambda)\) is purely imaginary for \(|\lambda|=1\), \(\lambda\ne-1\). Approximate the imaginary diagonal entries of \(K\) by rational imaginary numbers and the upper triangular entries by Gaussian rationals, setting the lower entries to the negative conjugates. The resulting skew-Hermitian \(K_m\) converge to \(K\); continuity of (3) gives \(C(K_m)\to U\).

An arbitrary unitary is a limit of unitaries without eigenvalue \(-1\): multiply it by phases tending to one while avoiding the finitely many forbidden phases determined by its spectrum. Approximating each of these unitaries by the previous construction proves density. ∎

**Proposition (rational witness reduction).** If MOC fails in a fixed dimension \(n\), it fails for Gaussian-rational spectra \(a,b\), a rational Cayley unitary, and a separating functional

\[
L_{p+iq}(z)=p\operatorname{Re}z+q\operatorname{Im}z,
\qquad p,q\in\mathbb Q,
\]

with a strictly positive rational gap against every one of the \(n!\) permutation points.

**Proof.** A failure yields a strict separator \(w\) by (1). Its gap

\[
g(a,b,U,w)=L_w(d_{a,b}(U))-
\max_{\sigma\in S_n}L_w(z_\sigma(a,b))
\tag{4}
\]

is positive. Determinants, products, and a finite maximum of continuous functions are continuous, so positivity persists on a neighborhood of these data. Gaussian-rational spectra and functionals are dense, and Lemma 2 supplies a rational Cayley unitary in that neighborhood. All determinants and vertex coordinates are then Gaussian rational. Hence all evaluated functional values and their positive gap are rational. ∎

The contrapositive is also useful: a proof covering all Gaussian-rational spectra and rational Cayley unitaries would prove the full conjecture, by this continuity argument. Checking **finitely many** such inputs cannot do so.

This gives a theoretical enumeration that halts if a counterexample exists: enumerate dimension and rational entries fairly; build skew-Hermitian \(K\); compute (3) and all vertices; then test rational functionals. Every individual operation is exact and finite, and a genuine violation eventually appears in such an enumeration. There is no termination guarantee if MOC is true, no effective denominator bound here, and no practical complexity claim. The implemented verifier currently caps dimension at six, so it alone does not implement this unrestricted enumeration.

## 5. Boundary cases and the remaining gap

These reductions accommodate repeated eigenvalues, singular sums, repeated vertices, and point or segment hulls. Lemma 2 handles the missing Cayley chart at eigenvalue \(-1\) by a limit; it does not discard those matrices. A strict violation is stable under sufficiently small perturbations, while exact boundary membership usually is not stable under arbitrary floating error.

What remains is a universal bound on (4), or one certified positive instance. None of the reductions supplies that bound. In particular, diagonal phase symmetry does not eliminate all complex phases, compactness provides a maximizer but does not make it a permutation matrix, and a path through unitary rotations does not guarantee convex interpolation of determinant values at its endpoints.

## 6. Reproduction and next handoff

This file is a proof draft to check line by line. The algebraic statements have no computational prerequisites. The adjacent [exact verifier](../certificates/README.md) supplies executable checks of individual witnesses; its test suite does not certify these universal derivations.

Suggested next task for another owner:

1. Check both directions of (1), including lower-dimensional hulls.
2. Check the Cayley signs, inverse formula, missing-eigenvalue limit, and continuity argument in (4).
3. Implement a bounded rational-Cayley sampler, using small Gaussian-integer skew-Hermitian matrices and rational spectra, and pass each resulting witness to the exact verifier. Preserve bit-size and CPU limits; denominator growth can dominate cost.
4. Investigate a bound on (4) for a substantive family. Treat a successful small family as a lemma, not a proof for all unitary matrices.

Background: the standard orbit formulation and established special-case context are described in [Bebiano and da Providência (2025)](https://www.mdpi.com/2227-7390/13/5/711). The proofs above are written out here so the research does not rely on an unexamined reformulation or a source's displayed formula. Independent review under a different GitHub owner is still pending.
