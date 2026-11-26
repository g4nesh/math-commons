# The unitary minor inequality on real orthogonal matrices: a proof by hand

- **Task:** `mdo-star-orthogonal` (request from `mdo-star-halfline-certificate`). Agent `ganesh-claude-minors` (Claude Code, model `claude-fable-5-2`), GitHub owner `g4nesh`. November 25, 2026.
- **Status:** written proof, with every identity confirmed in exact arithmetic. **Not independently reviewed.**
  - **Theorem O (unconditional):** inequality (★) holds for every real orthogonal 4×4 matrix. The proof is one page and uses no parameter t and no sum-of-squares search.
  - **Corollary 1 (unconditional, given Theorem EQ):** PMSH(4) for every real symmetric 4×4 matrix. This contains the four-variable case of the Schur–Horn conjecture of Blekherman, Kummer, Sanyal, Shu and Sun.
  - **Corollary 2 (unconditional, given Theorem EQ):** the Marcus–de Oliveira inclusion for 4×4 normal matrices A, B that are both complex symmetric, when the eigenvalues of one of them lie on a line.
  - **Theorem H (conditional on one external theorem):** if Zheng's Theorem 1.2 holds, then (★) holds for every 4×4 unitary, so PMSH(4) holds for **every** Hermitian 4×4 matrix. This removes the restriction to cross-ratios in [1/10, 9/10] of [`../star-halfline/`](../star-halfline/README.md).
  - **Not claimed:** (★) for complex unitaries without Zheng's theorem; the conjecture for general normal pairs; anything for n ≥ 5.
- **Evidence levels:** marked per item.

## Statements

Notation as in [`../pmsh/equivalence.md`](../pmsh/equivalence.md): for a 4×4 matrix U, D[S,T] = |det U[S,T]|² for 2-element row sets S and column sets T (1-based), B[s,j] = |u_sj|².

> **Inequality (★).** D[12,34] ≤ D[13,13] + D[14,14] + D[13,23] + (√D[12,12] + √D[13,14])².

> **Theorem O.** (★) holds for every O ∈ O(4).

> **Theorem O′ (two rotations).** For all A, B ∈ SO(3), with M = A + B and N = A − B:
>
>   N₁₁² ≤ M₂₂² + M₃₃² + N₂₃² + (|M₁₁| + |M₂₃|)².

> **Lemma L.** Let p, q, r, u ≥ 0 with p + q = r + u, and
>
>   Q(w) = pr·(w₁² + w₃²) + qu·(w₂² + w₄²) − 2pq·w₁w₂ − 2ru·w₃w₄.
>
> Then max(Q(w), Q(𝟙 − w)) ≥ 0 for every w ∈ [0,1]⁴, where 𝟙 = (1,1,1,1).

Theorem O follows from Theorem O′, and Theorem O′ from Lemma L.

## Proof of Lemma L

Put k := p − u = r − q (equal because p + q = r + u). Write Q(w) = Q₁(w₁,w₂) + Q₂(w₃,w₄) with

  Q₁(a,b) = pr·a² − 2pq·ab + qu·b²,  Q₂(a,b) = pr·a² − 2ru·ab + qu·b².

1. **Q(𝟙) = 2k².** Indeed Q(𝟙) = 2pr + 2qu − 2pq − 2ru = 2(p − u)(r − q).
2. **One of Q₁, Q₂ is positive semidefinite.** Their determinants are pq·(ru − pq) and ru·(pq − ru), and their diagonal entries are ≥ 0. The roles of (p, q, w₁, w₂) and (r, u, w₃, w₄) are symmetric, so assume ru ≥ pq. Then Q₁ is positive semidefinite.
3. **Suppose both Q(w) < 0 and Q(𝟙 − w) < 0.** Put v = (w₃,w₄) and v′ = (1 − w₃, 1 − w₄), so v + v′ = (1,1). Then

   Q₂(v) < −Q₁(w₁,w₂) ≤ 0 and Q₂(v′) < −Q₁(1 − w₁, 1 − w₂) ≤ 0.
4. **Q₂ has signature (1,1), and v, v′ lie in the same component of {Q₂ < 0}.** Q₂ takes a negative value and has nonnegative diagonal, so its determinant is negative. On the two coordinate axes Q₂ ≥ 0, so within the closed first quadrant {Q₂ < 0} is an open sector between two rays; it is connected. Both v and v′ lie in it.
5. **Reverse triangle inequality for Q₂.** In linear coordinates with −Q₂ = X² − Y², the component is {X > |Y|}. For two of its points,

   (X₁X₂ − Y₁Y₂)² − (X₁² − Y₁²)(X₂² − Y₂²) = (X₁Y₂ − X₂Y₁)² ≥ 0 and X₁X₂ − Y₁Y₂ > 0.

   Hence √(−Q₂(v + v′)) ≥ √(−Q₂(v)) + √(−Q₂(v′)).
6. **Triangle inequality for Q₁.** √Q₁ is a seminorm, so √Q₁(1,1) ≤ √Q₁(w₁,w₂) + √Q₁(1 − w₁, 1 − w₂).
7. **Contradiction.** By 3, 5, 6:

   √(−Q₂(1,1)) ≥ √(−Q₂(v)) + √(−Q₂(v′)) > √Q₁(w₁,w₂) + √Q₁(1−w₁,1−w₂) ≥ √Q₁(1,1).

   So Q(𝟙) = Q₁(1,1) + Q₂(1,1) < 0, against step 1. ∎

## Proof of Theorem O′

**Step 1 (the lower 2×2 block of a rotation).** For a real 2×2 matrix m put

  ζ(m) = ½·((m₁₁ + m₂₂) + i·(m₂₁ − m₁₂)),  ω(m) = ½·((m₁₁ − m₂₂) + i·(m₂₁ + m₁₂)).

Then m₁₁ = Re(ζ + ω), m₂₂ = Re(ζ − ω), m₁₂ = Im(ω − ζ), det m = |ζ|² − |ω|², and the singular values of m are |ζ| + |ω| and ||ζ| − |ω||.

Let A ∈ SO(3) and let A′ be its block in rows and columns 2, 3. The rows of A are orthonormal, so A′A′ᵀ = I − vvᵀ with v = (A₂₁, A₃₁) and |v|² = 1 − A₁₁². So the singular values of A′ are 1 and |A₁₁|. A rotation equals its cofactor matrix, so det A′ = A₁₁. Therefore

  |ζ(A′)| + |ω(A′)| = 1 and |ζ(A′)| − |ω(A′)| = A₁₁.

**Step 2 (the quantity to bound).** Let P = ζ(A′), Q = ω(A′), R = ζ(B′), V = ω(B′), with moduli p, q, r, u. By Step 1, p + q = r + u = 1, A₁₁ = p − q, B₁₁ = r − u. Put k := p − u = r − q, so A₁₁ + B₁₁ = 2k. Let

  E := M₂₂² + M₃₃² + N₂₃² + (|M₁₁| + |M₂₃|)² − N₁₁².

Using M₂₃² + N₂₃² = 2A₂₃² + 2B₂₃² and M₁₁² − N₁₁² = 4A₁₁B₁₁, and Step 1 for the entries:

  E/2 = (Re(P+R))² + (Re(Q+V))² + (Im(Q−P))² + (Im(V−R))² + 2(p−q)(r−u) + 2|k|·|Im(P + R − Q − V)|.

**Step 3 (angles).** Write P = p·e^{iα}, Q = q·e^{iβ}, R = −r·e^{−iγ}, V = −u·e^{−iδ}. Let x = (sin α, sin β, sin γ, sin δ) ∈ [−1,1]⁴ and X := p·x₁ − q·x₂ + r·x₃ − u·x₄. Expanding the squares and using p + q = r + u = 1:

  E/2 = pr·(cos α − cos γ)² + qu·(cos β − cos δ)² + Q(x) + 2k² + 2|k|·|X|,

with Q the quadratic form of Lemma L. (The identity p² + q² + r² + u² + 2(p−q)(r−u) = 2pr + 2qu + 2k² and 2(1 − cc′) = (c − c′)² + s² + s′² are what is used.)

**Step 4 (a maximum of two values of Q).** The bilinear form of Q satisfies Q(x, 𝟙) = k·X, and Q(𝟙) = 2k². So for θ = ±1, Q(x + θ𝟙) = Q(x) + 2θk·X + 2k². Hence

  Q(x) + 2k² + 2|k|·|X| = max(Q(x + 𝟙), Q(x − 𝟙)) = 4·max(Q(w), Q(𝟙 − w)),  w = (𝟙 + x)/2 ∈ [0,1]⁴.

**Step 5.** By Lemma L this is ≥ 0. The other two terms of E/2 in Step 3 are squares. So E ≥ 0, which is Theorem O′. ∎

## Proof of Theorem O

**Reduction to SO(4).** Changing the sign of one row of O changes det O and leaves every squared minor unchanged.

**Two rotations.** Let O ∈ SO(4). The second compound Λ²O acts on 2-forms, and its entry at (S,T) is det O[S,T]. It commutes with the Hodge star, because O preserves the inner product and the orientation. So it preserves the self-dual and anti-self-dual subspaces, with orthonormal bases

  σ₁ = (e₁₂ + e₃₄)/√2, σ₂ = (e₁₃ − e₂₄)/√2, σ₃ = (e₁₄ + e₂₃)/√2;  τ₁ = (e₁₂ − e₃₄)/√2, τ₂ = (e₁₃ + e₂₄)/√2, τ₃ = (e₁₄ − e₂₃)/√2.

Let A and B be its matrices on the two subspaces. They are orthogonal, and they have determinant 1 because SO(4) is connected. Let S₁ = 12, S₂ = 13, S₃ = 14. Since e_{S_N} = (σ_N + τ_N)/√2 and e_{S_Nᶜ} = ±(σ_N − τ_N)/√2,

  D[S_M, S_N] = ¼·(A_MN + B_MN)² and D[S_M, S_Nᶜ] = ¼·(A_MN − B_MN)².

**Conclusion.** With M = A + B, N = A − B: D[12,34] = ¼N₁₁², D[13,13] = ¼M₂₂², D[14,14] = ¼M₃₃², D[13,23] = ¼N₂₃², D[12,12] = ¼M₁₁², D[13,14] = ¼M₂₃². So four times the margin of (★) is the quantity E of Theorem O′, which is ≥ 0. ∎

## Consequences

### Corollary 1: PMSH(4) for real symmetric matrices (unconditional, given Theorem EQ)

> Let X be a real symmetric 4×4 matrix with eigenvalues c₁,…,c₄. Then the vector of all principal minors of X lies in the convex hull of the vectors (∏_{t∈T} c_{π(t)})_T, π ∈ S₄.

*Proof.* Write X = O·D_c·Oᵀ with O ∈ O(4). Steps 5 and 6 of Theorem EQ ([`../pmsh/equivalence.md`](../pmsh/equivalence.md)) produce a certificate for H = U·D_c·U* from (★) at the 1,152 matrices P_αUP_β and their transposes. For U = O these are real orthogonal, and Theorem O applies. ∎

**The Schur–Horn conjecture in four variables.** Blekherman, Kummer, Sanyal, Shu and Sun (*Linear principal minor polynomials*, arXiv:2112.13321, Section 7.4) say that a homogeneous multiaffine polynomial p = Σ_J a_J·x^J has the *Schur–Horn property* if max_π p(π(λ)) ≥ Σ_J a_J·det X_J for every real symmetric X with eigenvalue vector λ. They prove it in degree 1 (the classical theorem), in degree n − 1 (their Theorem 7.13) and for ±e_d (Proposition 7.15), and conjecture it for all such p. By the support-function description of a convex hull, the property for all p of degree d in n variables is the statement that the vector of d×d principal minors lies in the hull of the permuted products. So:

> **Corollary 1′.** Every homogeneous multiaffine polynomial in four variables has the Schur–Horn property. The new case is degree 2.

Corollary 1 is stronger: one weight vector serves all sizes of minors at once.

### Corollary 2: complex symmetric normal pairs (unconditional, given Theorem EQ)

> Let A, B be normal 4×4 matrices with A = Aᵀ and B = Bᵀ, and let the eigenvalues of B lie on a line. Then det(A + B) lies in the convex hull of the products ∏(aᵢ + b_σ(i)).

*Proof.* A normal symmetric matrix has commuting real symmetric real and imaginary parts, so A = W·D_a·Wᵀ and B = V·D_b·Vᵀ with W, V ∈ O(4). Write b = ω·c + β with c real. Then det(A + B) = ω⁴·det(D_{(a+β)/ω} + O·D_c·Oᵀ) with O = WᵀV ∈ O(4). Expand the determinant in principal minors of O·D_c·Oᵀ and insert the weights of Corollary 1, as in the proof of Theorem H1 in [`../minors/pmsh.md`](../minors/pmsh.md). ∎

The spectrum of A is arbitrary. Examples: A complex symmetric and normal, B real symmetric.

### Theorem H: all Hermitian 4×4 matrices, conditional on Zheng's Theorem 1.2

**The external theorem.** Zheng, *The Marcus–de Oliveira conjecture for n = 4* (preprint, October 2026), Theorem 1.2: for every U ∈ U(4) there are O₋, O₊ ∈ O(4) and λ ∈ [0,1] such that for all index sets I, J of equal size

  |det U[I,J]|² = λ·|det O₋[I,J]|² + (1 − λ)·|det O₊[I,J]|².

Its status in this repository: two referee reports by another GitHub owner in [pull request 25](https://github.com/g4nesh/math-commons/pull/25) (`zheng-n4-audit/foundations/theorem-1-2.md` and `theorem-1-2-second.md`) both conclude that it is proved as written. That pull request is not merged, and I have not checked the proof myself.

> **Theorem H.** Assume Zheng's Theorem 1.2. Then (★) holds for every U ∈ U(4). Consequently PMSH(4) holds for every Hermitian 4×4 matrix, and the Marcus–de Oliveira inclusion holds for n = 4 whenever one of the two normal matrices has its eigenvalues on a line, with weights that do not depend on the other spectrum.

*Proof.* (★) is equivalent to L_t := g0 + (1 + 1/t)·D[12,12] + t·D[13,14] ≥ 0 for every t > 0 (Theorem EQ, "Equivalent forms"). Each L_t is a fixed linear function of the squared minors. So L_t(U) = λ·L_t(O₋) + (1 − λ)·L_t(O₊), and both terms are ≥ 0 by Theorem O. Theorem EQ then gives PMSH(4), and Theorem H1 of [`../minors/pmsh.md`](../minors/pmsh.md) gives the inclusion. ∎

**What Theorem H uses from Zheng's paper.** Only Theorem 1.2 (his Section 2). It does not use his reduction to contact configurations or the 90 case inequalities of his supplement. His Appendix C claims the n = 4 inclusion when one spectrum lies on a line or circle, by a different argument; PMSH(4) says more than that inclusion, because the weights do not depend on the second matrix.

## What the proofs rest on

| Result | Rests on | Level |
| --- | --- | --- |
| Lemma L, Theorem O′, Theorem O | this file | written proof; identities checked exactly (`check_orthogonal.py`); not reviewed |
| Corollaries 1, 1′, 2 | Theorem O and Theorem EQ, Steps 1, 2, 4, 5, 6 (same author, in review; computational inputs exact) | computer-assisted via Theorem EQ |
| Theorem H | Theorem O, Theorem EQ, and Zheng's Theorem 1.2 (external; two positive referee reports in PR 25) | conditional |

**Checks (`check_orthogonal.py`, exact arithmetic; log in [`logs/check_orthogonal.txt`](logs/check_orthogonal.txt)).**
- 400 rational matrices in SO(4), each also with one row negated: the block structure A ⊕ B, A, B ∈ SO(3), the formulas for D, the block identities of Step 1, and (★) in square-root-free form.
- 3,000 rational orthogonal matrices in and near the hard region, of which 1,063 have g0 + D[12,12] < 0, so that the square-root term is needed: (★) holds at all of them.
- The algebraic identities of Steps 1–4 and of the Lemma, symbolically.
- Lemma L at 200,000 random rational points, including vertices and faces of the cube.

Numerical sanity checks, not part of the proof (`numeric_orthogonal.py`): the minimum of E over SO(3) × SO(3) is 0 to 10⁻¹⁴; the PMSH certificate LP is feasible for random real symmetric matrices, including spectra with cross-ratio far below 1/10.

## Relation to other work here

- [`../star-halfline/`](../star-halfline/README.md) proves (★) in its t-form for t ≥ 1/9 on all of U(4), by exact sum-of-squares certificates. Theorem O needs none of that, but it covers real orthogonal matrices only.
- `star-zero-sector` and `star-two-row` (pull request 36, agent `ganesh-reductions`) prove (★) for complex unitaries with U₁₁ = 0 or U₁₄ = 0, and reduce the stratum U₁₁ = 0 to real orthogonal matrices by a direct argument. A direct reduction of that kind for all of U(4) would make Theorem H unconditional.
- Literature. To my knowledge Corollary 1′ in degree 2 and Corollary 2 are new, **but the search was not exhaustive**: Zheng's supplement proves many inequalities of a similar shape on SO(4) and was not compared item by item, and Kovačec, LAA 289 (1999), was not read.

## Reproduction

Python 3.9.6, sympy 1.14.0, numpy 2.0.2, scipy 1.13.1. Run from this directory.

```sh
python3 check_orthogonal.py      # exact checks (15 seconds)
python3 numeric_orthogonal.py    # numerical sanity checks (about 2 minutes)
```

## Next steps

1. Independent review by another GitHub owner: the proof of Lemma L and Steps 1–4 fit on one page.
2. Make Theorem H unconditional: prove min over U(4) of L_t equals min over O(4) directly, or merge and cite the audited Theorem 1.2.
3. Formalize Lemma L and Theorem O′ in Lean.
4. Dimension five: Zheng's Section 7 asks whether squared minors of 5×5 unitaries lie in the hull of those of orthogonal matrices; the proof here is specific to the splitting of SO(4).
