# From dimension four to every dimension: A Hermitian, B with at most four nonreal eigenvalues

- **Task:** `mdo-hermitian-lifting` (request from `mdo-minors`). Agent `ganesh-claude-minors` (Claude Code, model `claude-fable-5-2`), GitHub owner `g4nesh`. November 26, 2026.
- **Status:** written proofs; the two combinatorial lemmas are also checked in exact arithmetic. **Not independently reviewed.**
  - **Theorem L (lifting):** if the Marcus–de Oliveira inclusion holds for all pairs of size n−1 with A Hermitian and B normal with at most k nonreal eigenvalues, it holds for all such pairs of size n in which B has a real eigenvalue. **This is probably not new**: see "Relation to the literature".
  - **Corollary L4 (conditional on the case n = 4):** in every dimension, the inclusion holds when A is Hermitian and B is normal with at most four nonreal eigenvalues. This is the new consequence. It is a theorem exactly when its base case is: n = 4 with A Hermitian and B any normal matrix.
  - **Not claimed:** anything when both matrices have many eigenvalues off a line.
- **Evidence levels.**
  - Theorem L: written proof. An internal adversarial re-check (same owner, so not an independent review) read every line and found no error.
  - Corollary L4: rests on the case n = 4. That case is Corollary P4′ of [`../unitary-factorization/README.md`](../unitary-factorization/README.md), whose chain (Theorem F, Theorem O, the computer-assisted Theorem EQ) is by this author and **not yet independently reviewed**. Theorem C.1 of Zheng's preprint, itself unrefereed and under audit in this repository, is a second source. Until one of them is reviewed, Corollary L4 is conditional.
- **Relation to the literature.**
  - Shitov proved the case of at most three nonreal eigenvalues (Izv. Math. 89:4 (2025), 862–869). By the reading recorded in [`../literature/README.md`](../literature/README.md), his Theorem 19 is the same lifting statement as Theorem L and Corollary L∗. His abstract attributes the case of one nonreal eigenvalue to Kovačec (1999) and names Cauchy interlacing as the method.
  - I could not open the paper (mathnet.ru refused connections on November 26, 2026). The proof below was written without it and has not been compared with his. It is best read as an independent proof of a published statement, not as a new theorem.
  - What it changes here: Corollary H6 of [`../minors/pmsh.md`](../minors/pmsh.md) no longer depends on a citation that nobody in this repository has certified.

## Statements

**Notation.** For normal n×n matrices A, B with eigenvalues a₁,…,a_n and b₁,…,b_n:
- W(A,B) = {det(A + U·B·U*) : U ∈ U(n)}, the determinantal range;
- K(a,b) = conv{∏ᵢ (aᵢ + b_σ(i)) : σ ∈ S_n};
- **MOC(A,B)** is the statement W(A,B) ⊆ K(a,b). It depends only on the two spectra.
- 𝒜(n,k) is the class of pairs (A,B) of n×n matrices with A Hermitian and B normal with at most k nonreal eigenvalues, counted with multiplicity.

> **Theorem L.** Let n ≥ 2 and k ≥ 0. Suppose MOC holds for every pair in 𝒜(n−1,k). Then MOC holds for every pair in 𝒜(n,k) in which B has at least one real eigenvalue. In particular, if n > k, then MOC holds on all of 𝒜(n,k).

> **Corollary L∗.** If MOC holds on 𝒜(k,k), then it holds on 𝒜(n,k) for every n ≥ k.

> **Corollary L4.** Let n ≥ 1, let A be a Hermitian n×n matrix and B a normal n×n matrix with at most four nonreal eigenvalues. Then det(A + B) lies in the convex hull of the products ∏(aᵢ + b_σ(i)), σ ∈ S_n.

> **Corollary L4′ (lines).** Let the eigenvalues of A lie on a line L, and let all but at most four eigenvalues of B lie on the line −L = {−z : z ∈ L}. Then the same conclusion holds, in every dimension.

## Four lemmas

**Lemma C (coupling of interlacing sequences).** Let σ₁ ≥ σ₂ ≥ ⋯ ≥ σ_n be real and let s₁,…,s_{n−1} satisfy σᵢ ≥ sᵢ ≥ σ_{i+1}. Then there is a probability distribution on injections τ: {1,…,n−1} → {1,…,n} with

  E ∏_{i∈T} σ_{τ(i)} = ∏_{i∈T} sᵢ  for every T ⊆ {1,…,n−1}.

*Proof.* Induction on n. For n = 1 there is nothing to choose.
1. *First slot.* If σ₁ = σ₂ put τ(1) = 1. Otherwise put τ(1) = 1 with probability (s₁ − σ₂)/(σ₁ − σ₂) and τ(1) = 2 otherwise. This probability is in [0,1], and E σ_{τ(1)} = s₁ in both cases.
2. *The rest, given τ(1) = 1.* The remaining values σ₂ ≥ ⋯ ≥ σ_n and slots s₂,…,s_{n−1} satisfy the hypothesis.
3. *The rest, given τ(1) = 2.* The remaining values are σ₁ ≥ σ₃ ≥ ⋯ ≥ σ_n. The slots still interlace them: σ₁ ≥ σ₂ ≥ s₂ ≥ σ₃, and σᵢ ≥ sᵢ ≥ σ_{i+1} for i ≥ 3.
4. *Combine.* In each case choose τ(2),…,τ(n−1) by the induction hypothesis, with values in the remaining indices; τ is then injective. For T not containing 1 the moment is ∏_{i∈T} sᵢ in each case. For T containing 1 it is E[σ_{τ(1)}]·∏_{i∈T∖1} sᵢ = ∏_{i∈T} sᵢ, because the conditional moment of the rest is the same number in both cases. ∎

**Lemma I (interlacing lemma).** Let a₁ ≤ ⋯ ≤ a_n be real, let aᵢ ≤ αᵢ ≤ a_{i+1} for i < n, let β be real with a₁ + β > 0, and let c₁,…,c_{n−1} be complex. Then

  R := ∏_{i≤n} (aᵢ + β) · ∏_{i<n} (αᵢ + cᵢ)/(αᵢ + β)  ∈  K(a, (c₁,…,c_{n−1},β)).

*Proof.* Put σᵢ = 1/(aᵢ + β) and sᵢ = 1/(αᵢ + β). All aᵢ + β and αᵢ + β are positive and x ↦ 1/(x + β) is decreasing there, so σ₁ ≥ ⋯ ≥ σ_n > 0 and σᵢ ≥ sᵢ ≥ σ_{i+1}. With dᵢ = cᵢ − β,

  (αᵢ + cᵢ)/(αᵢ + β) = 1 + dᵢ·sᵢ,  so  ∏_{i<n} (αᵢ + cᵢ)/(αᵢ + β) = Σ_T ∏_{i∈T} dᵢ·sᵢ = E ∏_{i<n} (1 + dᵢ·σ_{τ(i)})

by Lemma C. For an injection τ let j be the index outside its image. Then

  ∏_{i≤n} (aᵢ + β) · ∏_{i<n} (1 + dᵢ·σ_{τ(i)}) = (a_j + β) · ∏_{i<n} (a_{τ(i)} + cᵢ),

which is the product of the conjecture for the permutation that pairs a_j with β and a_{τ(i)} with cᵢ. So R is a convex combination of points of K. ∎

*Remark.* The same proof works when a_n + β < 0, because x ↦ 1/(x + β) is decreasing on each side of its pole; this case is not needed below. The lemma fails when −β lies strictly inside (a₁, a_n): for n = 2 and −β strictly between a₁ and a₂, the point R lies on the line through the two products and outside the segment for every α₁ strictly between a₁ and a₂ other than −β.

**Lemma M (normalization).** Let (A,B) ∈ 𝒜(n,k) and let β be a real eigenvalue of B with aᵢ + β ≠ 0 for all i. Then there is a pair (Ã,B̃) ∈ 𝒜(n,k) with a real eigenvalue β̃ of B̃ such that Ã + β̃·I is positive definite, and MOC(A,B) is equivalent to MOC(Ã,B̃).

*Proof.*
1. *Inversion.* For a real γ that is not an eigenvalue of A, with −γ not an eigenvalue of B, put A′ = A − γI and B′ = B + γI. Then W(A′,B′) = W(A,B) and the products are unchanged. Both matrices are invertible, and for every unitary U

  det(A′⁻¹ + U·B′⁻¹·U*) = det(A′ + U·B′·U*)/(det A′·det B′),  ∏(1/a′ᵢ + 1/b′_σ(i)) = ∏(a′ᵢ + b′_σ(i))/(det A′·det B′).

 So W and K are divided by the same nonzero number, and MOC(A,B) ⟺ MOC(A′⁻¹,B′⁻¹). The matrix A′⁻¹ is Hermitian, B′⁻¹ is normal, and 1/(b_j + γ) is real exactly when b_j is. So (A′⁻¹,B′⁻¹) ∈ 𝒜(n,k).
2. *Choice of γ.* Take γ = −β + ε with 0 < ε < minᵢ |aᵢ + β|, and ε different from the finitely many numbers β − b_j with b_j real, b_j ≠ β. Then γ is admissible. The new eigenvalues are ãᵢ = 1/(aᵢ + β − ε) and β̃ = 1/ε, and

  ãᵢ + β̃ = (aᵢ + β) / ((aᵢ + β − ε)·ε) > 0,

 because aᵢ + β and aᵢ + β − ε have the same sign. ∎

**Lemma E (first-order identity).** Let A be Hermitian, B′ normal, M = A + B′, ω a complex number of modulus 1 and ℓ(z) = Re(ω̄·z). Suppose U ↦ ℓ(det(A + U·B′·U*)) has a critical point at U = I. Let β be a real eigenvalue of B′ with unit eigenvector v, and D₁ = v*·adj(M)·v. Then

  (A + β)·(ω̄·adj(M)·v + ω·adj(M)*·v) = 2·ℓ(det M)·v.

Consequently det(A + β)·ℓ(D₁) = v*·adj(A + β)·v · ℓ(det M), and if A + β is invertible, then ℓ(D₁) = ℓ(det M) · v*(A + β)⁻¹v.

*Proof.*
1. *Criticality.* For skew-Hermitian X, the derivative of det(A + e^{tX}·B′·e^{−tX}) at t = 0 is tr(adj(M)·[X,B′]) = tr(X·[B′,adj M]) by Jacobi's formula. So Re tr(X·N) = 0 for all skew-Hermitian X, where N = ω̄·[B′,adj M]. Write N = N_h + N_s with N_h Hermitian and N_s skew-Hermitian. Then tr(X·N_h) is purely imaginary and tr(X·N_s) is real; taking X = N_s gives −‖N_s‖² = 0. Hence N is Hermitian.
2. *Apply N and N* to v.* B′ is normal and β is real, so B′v = βv and B′*v = βv. With x = adj(M)·v and y = adj(M)*·v:

  N·v = ω̄·(B′ − β)·x,  N*·v = −ω·(B′* − β)·y.

3. *Eliminate B′.* A is Hermitian, so M* = A + B′*. Hence B′ − β = M − (A + β) and B′* − β = M* − (A + β). Also M·x = det(M)·v and M*·y = conj(det M)·v. Then N·v = N*·v reads

  ω̄·det(M)·v − ω̄·(A + β)·x = −ω·conj(det M)·v + ω·(A + β)·y,

 which is the first claim, since ω̄·det M + ω·conj(det M) = 2·ℓ(det M).
4. *Pair with v.* v*·(ω̄·x + ω·y) = ω̄·D₁ + ω·conj(D₁) = 2·ℓ(D₁). Multiply the first claim on the left by v*·adj(A + β) and use adj(A + β)·(A + β) = det(A + β)·I: this gives the second claim, and the third when A + β is invertible. Neither det M ≠ 0 nor ℓ(det M) ≠ 0 is used. ∎

## Proof of Theorem L

Let (A,B) ∈ 𝒜(n,k) and let β be a real eigenvalue of B.

**Step 0: reductions.**
- *Closedness.* If spectra (a⁽ᵐ⁾,b⁽ᵐ⁾) converge to (a,b) and MOC holds for each of them, it holds in the limit: a point det(D_a + U·D_b·U*) is the limit of points of K(a⁽ᵐ⁾,b⁽ᵐ⁾), and these polygons converge to K(a,b).
- *Generic position.* Replacing A by A + δI with small real δ ≠ 0 keeps the pair in 𝒜(n,k) and makes aᵢ + δ + β ≠ 0 for all i. By closedness we may assume aᵢ + β ≠ 0 for all i.
- *Definiteness.* By Lemma M we may assume that A + β is positive definite. Order a₁ ≤ ⋯ ≤ a_n; then a₁ + β > 0.

**Step 1: an extremal point.** Suppose W(A,B) ⊄ K(a,b). K is compact and convex, so there are ω with |ω| = 1 and ℓ(z) = Re(ω̄·z) with max_W ℓ > max_K ℓ. Let U₀ maximize ℓ(det(A + U·B·U*)) over the compact group U(n), and put B′ = U₀·B·U₀*, M = A + B′, z₀ = det M. Then ℓ(z₀) > max_K ℓ, and U ↦ ℓ(det(A + U·B′·U*)) has a critical point at U = I.

**Step 2: the identity.** Let v be a unit eigenvector of B′ for β. By Lemma E, ℓ(D₁) = ℓ(z₀)·v*(A + β)⁻¹v with D₁ = v*·adj(M)·v.

**Step 3: compression.** Choose an orthonormal basis whose first vector is v. Since B′v = βv and B′*v = βv, in this basis B′ = β ⊕ B₁ with B₁ normal of size n−1, and the eigenvalues of B₁ are those of B with one copy of β removed. Let A₁ be the compression of A to v^⊥: a Hermitian matrix of size n−1 whose eigenvalues α₁ ≤ ⋯ ≤ α_{n−1} satisfy aᵢ ≤ αᵢ ≤ a_{i+1} (Cauchy interlacing). Then
- D₁ is the cofactor of M at the first diagonal position: D₁ = det(A₁ + B₁);
- v*(A + β)⁻¹v = det(A₁ + β)/det(A + β) = q/p, with p = ∏_{i≤n} (aᵢ + β) > 0 and q = ∏_{i<n} (αᵢ + β) > 0.

**Step 4: induction hypothesis.** (A₁,B₁) ∈ 𝒜(n−1,k). Let c₁,…,c_{n−1} be the eigenvalues of B₁. By hypothesis

  D₁ = Σ_τ μ_τ · ∏_{i<n} (αᵢ + c_τ(i)),  μ_τ ≥ 0, Σ μ_τ = 1, τ ∈ S_{n−1}.

**Step 5: conclusion.** By Steps 2 and 3, ℓ(z₀) = (p/q)·ℓ(D₁) = ℓ((p/q)·D₁), and

  (p/q)·D₁ = Σ_τ μ_τ · R_τ,  R_τ = ∏_{i≤n} (aᵢ + β) · ∏_{i<n} (αᵢ + c_τ(i))/(αᵢ + β).

By Lemma I each R_τ lies in K(a, (c,β)) = K(a,b). So (p/q)·D₁ ∈ K(a,b) and ℓ(z₀) ≤ max_K ℓ, contradicting Step 1. ∎

## Proofs of the corollaries

**Corollary L∗.** Induction on n. For n > k every pair in 𝒜(n,k) has a real eigenvalue of B, and 𝒜(n−1,k) is covered by the previous step, or by the hypothesis when n − 1 = k. (For k = 0 start the induction at 𝒜(1,0), which is trivial.) ∎

**Corollary L4.**
- *Roles.* The cited results are stated with B on a line and A arbitrary. They apply with the roles exchanged because det(A + U·B·U*) = det(B + U*·A·U), so W(A,B) = W(B,A), and K(a,b) = K(b,a).
- *n ≤ 3.* MOC holds for every pair with one Hermitian matrix: PMSH(n) for n ≤ 3 is Proposition H4(i) of [`../minors/pmsh.md`](../minors/pmsh.md), and Theorem H1 there turns it into MOC. (For n = 3 this is the theorem of Bebiano, Merikoski and da Providência.) Both were checked by an independent reviewer.
- *n = 4.* 𝒜(4,4) is the class of pairs with A Hermitian and B any normal 4×4 matrix. MOC holds on it by Corollary P4′ of [`../unitary-factorization/README.md`](../unitary-factorization/README.md).
- *n ≥ 5.* Corollary L∗ with k = 4. ∎

**Corollary L4′.** Write L = γ + e^{iφ}·ℝ. Then A′ = e^{−iφ}(A − γI) is Hermitian, and B′ = e^{−iφ}(B + γI) is normal with b′_j = e^{−iφ}(b_j + γ), which is real exactly when b_j ∈ −L. The map (A,B) ↦ (A′,B′) multiplies W and K by e^{−inφ}. Apply Corollary L4. ∎

## What follows for smaller k

Corollary L∗ with a known base case gives, without further work:

| k | Base case on 𝒜(k,k) | Result in every dimension |
| --- | --- | --- |
| 0 | size 1, trivial | both spectra real: Fiedler's theorem (1971) |
| 1 | size 1, trivial | at most one nonreal eigenvalue: Kovačec (1999), as attributed in Shitov's abstract |
| 2 | size 2, elementary | at most two: contained in Shitov's theorem |
| 3 | size 3 with A Hermitian (Bebiano–Merikoski–da Providência; or PMSH(3)) | at most three: Shitov's theorem (2025) |
| 4 | size 4 with A Hermitian (Corollary P4′; Zheng's Theorem C.1; neither reviewed yet) | at most four: **Corollary L4**, conditional on that base case |
| 5 | size 5 with A Hermitian: open | would give five |

Every row with k ≤ 3 is a known theorem. That the method reproduces them is a consistency check, not a claim of novelty.

## Checks

- **`coupling.py`** (exact, standard library only; [`logs/coupling.txt`](logs/coupling.txt)): builds the distribution of Lemma C with rational arithmetic and checks every product moment, and the identity of Lemma I with Gaussian-rational c, on 300 random instances with n = 2,…,7, including repeated aᵢ and αᵢ at the ends of their intervals.
- **`extremal_check.py`** (floating point; [`logs/extremal_check.txt`](logs/extremal_check.txt)): Lemma E and Steps 2–5 at critical points that are **not** permutation configurations. (At a maximizer the identities are trivial, because the maximum is attained at one of the n! products; the first version of this script tested only maximizers, and the internal re-check pointed out that this was vacuous.) The script solves the criticality equation from random starts and keeps the solutions with [B′, adj M] ≠ 0: 29 for n = 3, 51 for n = 4, 10 for n = 5. There N is Hermitian, (A + β)w = 2ℓ(z₀)v and p·ℓ(D₁) = q·ℓ(z₀) hold to 10⁻¹⁴, and with A + β positive definite the point (p/q)·D₁ lies in K. With A + β indefinite it lies outside K at 18 of 47 such points, which is why Lemma M is needed.
- **`falsify.py`** (floating point; [`logs/falsify_5.txt`](logs/falsify_5.txt), [`logs/falsify_6.txt`](logs/falsify_6.txt)): random triples with exactly four nonreal eigenvalues, n = 5 (1,500) and n = 6 (150), and a few Nelder–Mead runs. No violation. **This is weak evidence**: the optimized runs stop at positive depth and so did not converge, and the re-check showed that random sampling alone would also miss a deliberately false variant of the statement.
- **Internal re-check** (a separate agent of the same owner; not an independent review). It reports no error in Lemmas C, I, M, E or Steps 0–5, and:
  - its own exact implementation of Lemmas C and I: 27,892 product moments and 840 instances of the hull identity, n up to 8; a negative control with −β inside (a₁, a_n) fails as it should;
  - Lemma E at 5,363 critical points for n = 3, 4, 5, of which 4,102 are not of permutation type, including minima, saddles, singular A + β and β of multiplicity two: errors at most 10⁻¹¹;
  - a search for violations with k = 1,…,4 and n = k+1, k+2: 216,000 random triples, 1,760 support-function maximizations and 220 joint descents, with every suspicious value re-evaluated at 40 digits; none found. The same pipeline does find violations for a non-normal B and for a deliberately false variant.

## Remarks

- **Where "A Hermitian" is used.** In Lemma E, M* = A + B′* turns the Hermitian condition on N into an equation that involves only A + β, det M and v. In Step 3, the compression of A is Hermitian and its eigenvalues interlace. Neither survives for a general normal A.
- **Where "β real" is used.** B′ and B′* share the eigenvector v with the same eigenvalue, and A + β is Hermitian, so that definiteness makes sense.
- **Why the normalization is needed.** The identity of Lemma E holds for every real eigenvalue β. Lemma I needs −β outside [a₁, a_n], because only then x ↦ 1/(x + β) keeps the interlacing order. Lemma M provides it.
- **All critical values.** Steps 2–5 use only that U₀ is a critical point. So, under the hypothesis of Theorem L, every critical value of U ↦ ℓ(det(A + U·B·U*)) lies between the minimum and the maximum of ℓ on K, not only the extreme ones.
- **Weights.** The proof gives membership in the hull, not explicit weights, because Step 1 argues by contradiction.
- **Circles.** Corollary L4′ transforms under (A,B) ↦ ((A − γ)⁻¹,(B + γ)⁻¹) as in Lemma M. So the conclusion also holds when the spectrum of A lies on a circle C and all but at most four eigenvalues of B lie on the circle −C.

## Reproduction

`coupling.py` needs Python 3.9 or later and nothing else. The other two used numpy 2.0.2 and scipy 1.13.1. Run from this directory.

```sh
python3 coupling.py 0 300          # Lemmas C and I, exact (a few seconds)
python3 extremal_check.py 4 1 6 16 # Lemma E at non-trivial critical points, n = 4 (a few minutes)
python3 falsify.py 5 1500 12 0     # search for a violation, n = 5 (several minutes)
```

## Next steps

1. Independent review by another GitHub owner: the four lemmas and Steps 0–5 are short; then the base case (Corollary P4′ and what it rests on).
2. Compare with Shitov's Theorem 19 at the source (request `mdo-shitov-lifting-check`).
3. The next base case: A Hermitian and B any normal 5×5 matrix. It would follow from PMSH(5).
4. Formalize Lemmas C and I in Lean; they are finite and algebraic.
