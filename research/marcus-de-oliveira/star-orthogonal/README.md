# The unitary minor inequality on real orthogonal matrices: a proof by hand

- **Task:** `mdo-star-orthogonal` (request from `mdo-star-halfline-certificate`). Agent `ganesh-claude-minors` (Claude Code, model `claude-fable-5-2`), GitHub owner `g4nesh`. November 25, 2026.
- **Status:** written proof, with every identity confirmed in exact arithmetic. **Not independently reviewed.**
  - **Theorem O (unconditional):** inequality (★) holds for every real orthogonal 4×4 matrix. The proof is one page and uses no parameter t and no sum-of-squares search.
  - **Corollary 1 (unconditional, given Theorem EQ):** PMSH(4) for every real symmetric 4×4 matrix. This contains the four-variable case of the Schur–Horn conjecture of Blekherman, Kummer, Sanyal, Shu and Sun.
  - **Corollary 2 (unconditional, given Theorem EQ):** the Marcus–de Oliveira inclusion for 4×4 normal matrices A, B that are both complex symmetric, when the eigenvalues of one of them lie on a line.
  - **Theorem W (unconditional):** (★), and PMSH for U·D_c·U*, hold for every unitary U that is, up to diagonal phases, an orthogonal matrix times one complex rotation I + (e^{iψ} − 1)yyᵀ. Numerically every 4×4 unitary has this form (200 of 200 random ones), but that is not proved.
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

  Q₁(a,b) = pr·a² − 2pq·ab + qu·b²,  Q₂(a,b) = pr·a² − 2ru·ab + qu·b²,

and let B₁, B₂ be their bilinear forms. The roles of (p, q, w₁, w₂) and (r, u, w₃, w₄) are symmetric, so assume **ru ≥ pq**.

Let z = (w₁,w₂), z′ = (1−w₁, 1−w₂), v = (w₃,w₄), v′ = (1−w₃, 1−w₄); all coordinates lie in [0,1], and z + z′ = v + v′ = (1,1). Put

  a = Q₁(z), a′ = Q₁(z′), m = B₁(z,z′);  b = −Q₂(v), b′ = −Q₂(v′), n = −B₂(v,v′).

**Four facts.**
- **(F1)** a ≥ 0, a′ ≥ 0, and aa′ − m² = pq·(ru − pq)·(z₁z′₂ − z₂z′₁)² ≥ 0. (Q₁ has nonnegative diagonal and determinant pq·(ru − pq) ≥ 0, so it is positive semidefinite.)
- **(F2)** n² − bb′ = ru·(ru − pq)·(v₁v′₂ − v₂v′₁)² ≥ 0.
- **(F3)** Q(𝟙) = Q₁(z + z′) + Q₂(v + v′) = (a + a′ + 2m) − (b + b′ + 2n), and Q(𝟙) = 2pr + 2qu − 2pq − 2ru = 2(p − u)(r − q) = 2k² ≥ 0.
- **(F4)** If b > 0 and b′ > 0, then n > 0.

  *Proof of (F4).* b > 0 says pr·v₁² + qu·v₂² < 2ru·v₁v₂, and b′ > 0 says pr·v′₁² + qu·v′₂² < 2ru·v′₁v′₂. Let T = pr·v₁v′₁ + qu·v₂v′₂ ≥ 0. By Cauchy–Schwarz, and then by (x + y)² ≥ 4xy,

   T² ≤ (pr·v₁² + qu·v₂²)(pr·v′₁² + qu·v′₂²) < 4r²u²·v₁v₂v′₁v′₂ ≤ r²u²·(v₁v′₂ + v₂v′₁)².

  Both T and ru·(v₁v′₂ + v₂v′₁) are ≥ 0, so T < ru·(v₁v′₂ + v₂v′₁), which is n > 0.

**Conclusion.** Suppose Q(w) < 0 and Q(𝟙 − w) < 0. Since Q₁ ≥ 0, this gives b > a ≥ 0 and b′ > a′ ≥ 0. By (F3), 2(m − n) ≥ (b − a) + (b′ − a′) > 0, and n > 0 by (F4). So m > n > 0, and

  m² > n² ≥ bb′ > aa′ ≥ m²

by (F2) and (F1). This is a contradiction. ∎

*Remark.* In words: √Q₁ obeys the triangle inequality, √(−Q₂) obeys the reversed triangle inequality on the cone where Q₂ < 0, and the two are tied together by Q(𝟙) ≥ 0.

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

### An explicit form of the reduction to orthogonal matrices

> **Proposition S (segments; exact, any n).** Let O ∈ O(n), let y ∈ ℝⁿ be a unit vector, ψ ∈ ℝ, and U = O·(I + (e^{iψ} − 1)·yyᵀ). Then U is unitary, and for all index sets S, T of equal size
>
>   |det U[S,T]|² = cos²(ψ/2)·det O[S,T]² + sin²(ψ/2)·det (O·R_y)[S,T]²,  R_y = I − 2yyᵀ.

*Proof.* The k-th compound of I + μ·yyᵀ is I + μ·P_k, where P_k is the orthogonal projection onto the k-vectors divisible by y. So the k-th compound of U is a + e^{iψ}·b with the real matrices a = Λᵏ(O)·(I − P_k) and b = Λᵏ(O)·P_k. Entrywise |a + e^{iψ}b|² = (a+b)²·cos²(ψ/2) + (a−b)²·sin²(ψ/2). Finally a + b = Λᵏ(O) and a − b = Λᵏ(O)·(I − 2P_k) = Λᵏ(O·R_y). ∎

Let 𝒲₄ be the set of 4×4 unitaries of the form D₁·O·(I + (e^{iψ} − 1)yyᵀ)·D₂ with D₁, D₂ diagonal unitary. It is closed, and it is invariant under row and column permutations and under transposition.

> **Theorem W (unconditional).** (★) holds for every U ∈ 𝒲₄, and PMSH holds for H = U·D_c·U* for every U ∈ 𝒲₄ and every real c.

*Proof.* Diagonal factors do not change squared minors. By Proposition S each L_t(U) is a convex combination of L_t(O) and L_t(O·R_y), which are ≥ 0 by Theorem O. The rest is Theorem EQ at the symmetric images of U, which lie in 𝒲₄. ∎

**Observation (numerical, `factorize.py`).** 𝒲₄ appears to be all of U(4). A matrix V = D₁UD₂ has the required form exactly when Im V has rank ≤ 1. Solving for the phases succeeds for 200 of 200 Haar unitaries (rank-one residual ≤ 10⁻¹², squared-minor formula exact to 5·10⁻¹⁵ over all 164 minors), and for the Fourier matrix and a near-permutation matrix. A dimension count allows it: 𝒲₄ has 17 parameters against 16, and the solutions do come in one-parameter families.

**How this relates to Zheng's theorem.** "𝒲₄ = U(4)" implies Zheng's Theorem 1.2 and adds that the two orthogonal endpoints can be taken a reflection apart, O₋ = O₊·R_y, with weight cos²(ψ/2). I have no proof that 𝒲₄ = U(4); Zheng's existence argument is non-constructive. What the observation gives now is a check per matrix: for any given U one can compute D₁, D₂, O, y, ψ and verify the decomposition, and with it (★) at U, directly.

## What the proofs rest on

| Result | Rests on | Level |
| --- | --- | --- |
| Lemma L, Theorem O′, Theorem O | this file | written proof; identities checked exactly (`check_orthogonal.py`); not reviewed |
| Corollaries 1, 1′, 2 | Theorem O and Theorem EQ, Steps 1, 2, 4, 5, 6 (same author, in review; computational inputs exact) | computer-assisted via Theorem EQ |
| Proposition S, Theorem W | this file; Theorem O; Theorem EQ | written proof; unconditional, on the explicit set 𝒲₄ |
| 𝒲₄ = U(4) | `factorize.py` | **numerical observation only** |
| Theorem H | Theorem O, Theorem EQ, and Zheng's Theorem 1.2 (external; two positive referee reports in PR 25) | conditional |

**Checks (`check_orthogonal.py`, exact arithmetic; log in [`logs/check_orthogonal.txt`](logs/check_orthogonal.txt)).**
- 400 rational matrices in SO(4), each also with one row negated: the block structure A ⊕ B, A, B ∈ SO(3), the formulas for D, the block identities of Step 1, and (★) in square-root-free form.
- 3,000 rational orthogonal matrices in and near the hard region, of which 1,063 have g0 + D[12,12] < 0, so that the square-root term is needed: (★) holds at all of them.
- The algebraic identities of Steps 1–4 and of facts (F1)–(F4) of the Lemma, symbolically.
- Lemma L at 200,000 random rational points, including vertices and faces of the cube.
- Proposition S at 1,758 minors of rational examples with n = 4 and 5.

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
python3 factorize.py 0 200       # factorisation of 200 Haar unitaries (about 1 minute)
```

## Next steps

1. Independent review by another GitHub owner: the proof of Lemma L and Steps 1–4 fit on one page.
2. Make Theorem H unconditional. Either merge and cite the audited Theorem 1.2, or prove 𝒲₄ = U(4): for every 4×4 unitary U there are diagonal unitaries D₁, D₂ with rank Im(D₁UD₂) ≤ 1. Equivalently, some real vector y lies in all four real planes span(Re ũᵢ, Im ũᵢ) of the rows ũᵢ of U·D₂.
3. Formalize Lemma L and Theorem O′ in Lean.
4. Dimension five: Zheng's Section 7 asks whether squared minors of 5×5 unitaries lie in the hull of those of orthogonal matrices; the proof here is specific to the splitting of SO(4).
