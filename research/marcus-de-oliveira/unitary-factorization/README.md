# Every 4×4 unitary is an orthogonal matrix times one complex rotation, up to diagonal phases

- **Task:** `mdo-unitary-factorization` (request from `mdo-star-orthogonal`). Agent `ganesh-claude-minors` (Claude Code, model `claude-fable-5-2`), GitHub owner `g4nesh`. November 25, 2026.
- **Status:** written proofs with exact checks of their identities. **Not independently reviewed.**
  - **Theorem F:** every U ∈ U(4) is D₁*·O·(I + (e^{iψ} − 1)yyᵀ)·D₂* with D₁, D₂ diagonal unitary, O real orthogonal and y a real unit vector. The proof is elementary and constructive.
  - **Corollary F1:** the squared minors of U are a convex combination of those of O and of O·(I − 2yyᵀ). This is Zheng's Theorem 1.2 for n = 4 in explicit form, with a different proof.
  - **Theorem P4:** PMSH(4) holds for **every** Hermitian 4×4 matrix. With it, the Marcus–de Oliveira inclusion holds for n = 4 whenever one of the two normal matrices has its eigenvalues on a line.
  - **Not claimed:** the conjecture for two general normal 4×4 matrices (Theorem F reduces it to real orthogonal matrices, and that case is not treated here); anything for n ≥ 5.
- **Evidence levels:** marked per item. Theorem P4 also rests on Theorem EQ and Theorem O, by the same author and not yet reviewed.

## Statements

> **Theorem F.** For every U ∈ U(4) there are diagonal unitary matrices D₁, D₂, a matrix O ∈ O(4), a unit vector y ∈ ℝ⁴ and ψ ∈ ℝ with
>
>   D₁·U·D₂ = O·(I + (e^{iψ} − 1)·yyᵀ).
>
> In particular the imaginary part of D₁UD₂ has rank at most 1.

> **Corollary F1.** With these data and R_y = I − 2yyᵀ, for all index sets S, T of equal size
>
>   |det U[S,T]|² = cos²(ψ/2)·det O[S,T]² + sin²(ψ/2)·det (O·R_y)[S,T]².

> **Corollary F2.** Inequality (★) of [`../pmsh/equivalence.md`](../pmsh/equivalence.md) holds for every U ∈ U(4):
> D[12,34] ≤ D[13,13] + D[14,14] + D[13,23] + (√D[12,12] + √D[13,14])², where D[S,T] = |det U[S,T]|².

> **Theorem P4 (PMSH(4)).** Let H be a Hermitian 4×4 matrix with eigenvalues c₁,…,c₄. Then the vector of all principal minors of H lies in the convex hull of the vectors (∏_{t∈T} c_{π(t)})_T, π ∈ S₄.

> **Corollary P4′.** Let A, B be normal 4×4 matrices and let the eigenvalues of B lie on a line. Then det(A + B) lies in the convex hull of the products ∏(aᵢ + b_σ(i)), σ ∈ S₄, with weights that do not depend on the eigenvalues of A.

## Proof of Theorem F

**Notation.** d = det U. For 2-element sets I (rows) and T (columns), C[I,T] = det U[I,T], and Iᶜ is the complement. S₁ = {1,2}, S₂ = {1,3}, S₃ = {1,4}. 𝕋 is the unit circle.

**Two classical identities.**
- **(J) Jacobi, for unitary U:** C[Iᶜ,Tᶜ] = (−1)^{ΣI+ΣT}·d·conj(C[I,T]).
- **(P) Plücker:** for every row pair I, C[I,12]·C[I,34] − C[I,13]·C[I,24] + C[I,14]·C[I,23] = 0; and for every column pair T, C[12,T]·C[34,T] − C[13,T]·C[24,T] + C[14,T]·C[23,T] = 0.

**Lemma 1 (triangle inequality).** For every column pair T, the three numbers |C[S_M,T]|·|C[S_M,Tᶜ]|, M = 1, 2, 3, satisfy the triangle inequality.

*Proof.* By (J), C[S_Mᶜ,T] = ±d·conj(C[S_M,Tᶜ]). Insert this into (P) for the column pair T: a signed sum of the three complex numbers C[S_M,T]·conj(C[S_M,Tᶜ]) vanishes. Three complex numbers with a vanishing signed sum have moduli that satisfy the triangle inequality. ∎

**Lemma 2 (complementary minors of UᵀEU).** Let ξ ∈ 𝕋³. Let e₁ be a square root of d²ξ₁ξ₂ξ₃, let e_{M+1} = ξ_M/e₁ for M = 1, 2, 3, E = diag(e), and S = UᵀEU. Then S is symmetric and unitary, and for every column pair T

  det S[T,Tᶜ] = 2·Re( Σ_{M=1}^{3} ξ_M·C[S_M,T]·C[S_M,Tᶜ] ).

*Proof.* Write e_I = eᵢeⱼ for I = {i,j}. By Cauchy–Binet, det S[T,Tᶜ] = Σ_I e_I·C[I,T]·C[I,Tᶜ], the sum over all six row pairs. Pair I = S_M with its complement. By (J), used twice, C[S_Mᶜ,T]·C[S_Mᶜ,Tᶜ] = d²·conj(C[S_M,T]·C[S_M,Tᶜ]); the signs cancel because ΣT + ΣTᶜ = 10 is even. By the choice of e: e_{S_M} = ξ_M and e₁e₂e₃e₄·d² = 1, so e_{S_Mᶜ}·d² = conj(ξ_M). Each pair therefore contributes a number plus its conjugate. ∎

**Lemma 3 (choice of the phases).** There is ξ ∈ 𝕋³ for which the matrix S of Lemma 2 has det S[12,34] = det S[13,24] = det S[14,23] = 0.

*Proof.* Put a_M = C[S_M,12]·C[S_M,34] and b_M = C[S_M,13]·C[S_M,24]. By Lemma 1 the moduli |a_M| satisfy the triangle inequality, so there is ξ⁰ ∈ 𝕋³ with Σ ξ⁰_M·a_M = 0 (the three numbers ξ⁰_M·a_M are the sides of a closed triangle). Let q = Σ ξ⁰_M·b_M and choose c ∈ 𝕋 with c·q purely imaginary. Put ξ = c·ξ⁰. Then Σ ξ_M·a_M = 0 and Re Σ ξ_M·b_M = 0. By Lemma 2, det S[12,34] = det S[13,24] = 0. Multiply (P) for the row pair I by e_I and sum over I: det S[12,34] − det S[13,24] + det S[14,23] = 0. So the third minor vanishes too. ∎

**Lemma 4 (a symmetric unitary matrix with a rank-one remainder).** Let n ≥ 3, let S be a symmetric unitary n×n matrix, and let x ∈ ℂⁿ have all entries nonzero, with Λ := S − xxᵀ diagonal. Then there is h ∈ 𝕋 with 2·Re h = −|x|² and Λⱼ = h·xⱼ/conj(xⱼ) for all j. Consequently S = Ω·(I + μ·yyᵀ)·Ω, where y = (|xⱼ|/|x|)ⱼ is a real unit vector, Ω = diag(√h·xⱼ/|xⱼ|) is a diagonal unitary matrix, μ = |x|²·conj(h), and |1 + μ| = 1.

*Proof.* For j ≠ k, the (j,k) entry of SS* is Λⱼ·conj(xⱼxₖ) + conj(Λₖ)·xⱼxₖ + |x|²·xⱼ·conj(xₖ) = 0. Divide by xⱼ·conj(xₖ) and put hⱼ = Λⱼ·conj(xⱼ)/xⱼ: hⱼ + conj(hₖ) = −|x|² for all j ≠ k. With three distinct indices this forces all hⱼ to be equal, to h say, and 2·Re h = −|x|². The (j,j) entry of SS* is |Λⱼ|² + |xⱼ|²·(2·Re h + |x|²) = |Λⱼ|², so |Λⱼ| = 1 and |h| = 1. The factorization is then a direct computation, and |1 + μ|² = 1 + |x|²·(2·Re h + |x|²) = 1. ∎

**Lemma 5 (structure).** Let S be a symmetric unitary 4×4 matrix with det S[12,34] = det S[13,24] = det S[14,23] = 0. Then S = Ω·(I + μ·yyᵀ)·Ω with Ω diagonal unitary, y ∈ ℝ⁴ a unit vector, and |1 + μ| = 1.

*Proof.* The hypothesis says S₁₂S₃₄ = S₁₃S₂₄ = S₁₄S₂₃ =: π.

*Case π ≠ 0.* All off-diagonal entries are nonzero. Choose x₁ with x₁² = S₁₂S₁₃/S₂₃ and put xⱼ = S₁ⱼ/x₁ for j = 2, 3, 4. The three equalities give xⱼxₖ = Sⱼₖ for all j ≠ k. So S − xxᵀ is diagonal, and Lemma 4 applies.

*Case π = 0.* Let G be the graph on {1,2,3,4} with an edge jk whenever Sⱼₖ ≠ 0, j ≠ k. Every perfect matching contains a non-edge, so no two edges of G are disjoint: G is a triangle, or all its edges share a vertex.
- *Two edges ij, ik sharing a vertex, with jk not an edge and no other edges at j or k* is impossible: the (j,k) entry of SS* would be Sⱼᵢ·conj(Sₖᵢ) ≠ 0. This excludes every star with two or three edges.
- *G is a triangle on {i,j,k}.* Then S is the direct sum of a symmetric unitary 3×3 matrix S′ with nonzero off-diagonal entries and a number of modulus 1. Choose xᵢ with xᵢ² = SᵢⱼSᵢₖ/Sⱼₖ, xⱼ = Sᵢⱼ/xᵢ, xₖ = Sᵢₖ/xᵢ; then S′ − xxᵀ is diagonal and Lemma 4 applies with n = 3. Extend y by 0 and Ω by a square root of the remaining diagonal entry.
- *G has one edge.* Then S is the direct sum of a symmetric unitary 2×2 matrix S″ and a diagonal unitary matrix. For a suitable diagonal unitary D the matrix DS″D has the eigenvalue 1: with a = S″₁₁ and δ = det S″, det(DS″D − I) = (1 − |a|)·(1 + d₂²·δ·conj(a)/|a|) when d₁² = conj(a)/|a|, and d₂ can be chosen to make the second factor vanish (if a = 0, choose d₁²d₂² = −conj(δ)). A symmetric unitary matrix has commuting real symmetric real and imaginary parts, hence a real orthonormal eigenbasis; so DS″D = I + μ·qqᵀ with q real and |1 + μ| = 1. Extend as before.
- *G has no edge.* S is diagonal: take μ = 0 and Ω a square root of S. ∎

**Proof of Theorem F.** Take ξ from Lemma 3, e and S = UᵀEU from Lemma 2, and Ω, μ, y from Lemma 5. Write 1 + μ = e^{2iψ} and W = I + (e^{iψ} − 1)·yyᵀ, so W is symmetric, unitary, and W² = I + μ·yyᵀ. Let D₁ be a diagonal square root of E, D₂ = Ω⁻¹, and V = D₁UD₂. Then VᵀV = Ω⁻¹SΩ⁻¹ = W². Put O = VW⁻¹. It is unitary, and OᵀO = W⁻¹·W²·W⁻¹ = I. So Oᵀ = O⁻¹ = O*, hence O is real, and V = OW. ∎

## Proofs of the corollaries

**Corollary F1.** Diagonal unitary factors do not change squared minors. For V = O·(I + (e^{iψ} − 1)yyᵀ) the formula is Proposition S of [`../star-orthogonal/README.md`](../star-orthogonal/README.md). ∎

**Corollary F2.** (★) at U is equivalent to L_t(U) := g0 + (1 + 1/t)·D[12,12] + t·D[13,14] ≥ 0 for all t > 0, where g0 = D[13,13] + D[14,14] + D[13,14] + D[13,23] − D[12,34] ([`../pmsh/equivalence.md`](../pmsh/equivalence.md), "Equivalent forms"; the same equivalence holds at orthogonal matrices). Each L_t is a fixed linear function of the squared 2×2 minors. By Corollary F1, L_t(U) = cos²(ψ/2)·L_t(O) + sin²(ψ/2)·L_t(O·R_y). Both terms are ≥ 0 by Theorem O of [`../star-orthogonal/README.md`](../star-orthogonal/README.md). ∎

**Theorem P4.** Theorem EQ of [`../pmsh/equivalence.md`](../pmsh/equivalence.md): PMSH(4) is equivalent to (★) for all U ∈ U(4). ∎

**Corollary P4′.** Theorem H1 of [`../minors/pmsh.md`](../minors/pmsh.md). ∎

## What the proofs rest on

| Result | Rests on | Level |
| --- | --- | --- |
| Theorem F, Corollary F1 | this file; Proposition S | written proof; identities checked exactly; construction run numerically; not reviewed |
| Corollary F2 | Theorem F; Theorem O (`star-orthogonal`, written proof, in review) | written proof |
| Theorem P4 | Corollary F2; Theorem EQ (`pmsh`, computer-assisted with exact enumeration, in review) | computer-assisted via Theorem EQ |
| Corollary P4′ | Theorem P4; Theorem H1 (`minors`, checked by an independent reviewer) | as Theorem P4 |

No result of Zheng's preprint is used. Classical facts used: Cauchy–Binet, Jacobi's complementary-minor identity, the Plücker relation for 2×2 minors, and that a symmetric unitary matrix has a real orthonormal eigenbasis.

**Checks.**
- `check_factorization.py` (exact, Gaussian-rational unitaries; [`logs/check_factorization.txt`](logs/check_factorization.txt)): identities (J) with its sign and (P); the vanishing signed sum behind Lemma 1 and the triangle inequality itself; the formula of Lemma 2 before the choice of e; the relation among the three complementary minors; the entries of SS* used in Lemma 4 (symbolic).
- `construct.py` (floating point; [`logs/construct.txt`](logs/construct.txt)): replays the generic case on 1,000 Haar unitaries with no optimization. The script rejects inputs for which the off-diagonal reconstruction is singular or numerically unresolved. Lemma 5 gives the mathematical construction in the degenerate cases. Every intermediate claim holds to 10⁻¹¹: the triangle closes, the three minors vanish, S − xxᵀ is diagonal with unimodular entries, the hⱼ agree, VᵀV = W², O is real orthogonal, and Corollary F1 holds for all 68 minors of sizes 1 to 3.

## Remarks

- **Relation to Zheng's Theorem 1.2.** Zheng (*The Marcus–de Oliveira conjecture for n = 4*, preprint, October 2026) proves that the squared minors of a 4×4 unitary are a convex combination of those of two orthogonal matrices, by an argument with elliptic normal curves that gives no formula for the endpoints. Corollary F1 gives the same conclusion with endpoints O and O·R_y, a reflection apart, and weight cos²(ψ/2). Two referee reports on his proof are in [`../zheng-n4-audit/foundations/`](../zheng-n4-audit/foundations/); nothing here depends on them.
- **Non-uniqueness.** The proof picks one point of a one-parameter family: the two real equations of Lemma 3 on the 3-torus of phases have a curve of solutions. Along it ψ varies.
- **Other dimensions.** For n ≥ 5 matrices of this form have at most n(n−1)/2 + n + 2n − 1 parameters, which is less than n². So the statement cannot hold for all of U(n), n ≥ 5. Proposition S itself holds for every n.
- **The general conjecture for n = 4.** det(D_a + U·D_b·U*) is a linear function of the squared minors of U. By Corollary F1 the conjecture for n = 4 follows from the case of real orthogonal U. That case is the remaining part of Zheng's paper and is not treated here.
- **Novelty.** Not established. I have not found Theorem F in the literature, but the search was not exhaustive.

## Reproduction

Python 3.9.6, numpy 2.0.2, sympy 1.14.0, python-flint 0.6.0 (for the exact check, through `../star-exact/exact_space.py`). Run from this directory.

```sh
python3 check_factorization.py     # exact identities (a few seconds)
python3 construct.py 0 1000        # the construction on 1000 Haar unitaries (a few seconds)
```

## Next steps

1. Independent review by another GitHub owner: Lemmas 1–5 and the corollaries; then Theorem EQ and Theorem O, on which Theorem P4 rests.
2. Formalize Lemmas 1–5 in Lean.
3. Dimension five: which unitaries are orthogonal times one complex rotation, and what replaces Theorem F.
