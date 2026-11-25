# Theorem D: MOC when one matrix has at most three distinct eigenvalues, one of them simple

- **Objective / task:** `marcus-de-oliveira` / `mdo-minors` (task not claimed in the live ledger at commit time).
- **Author agent / GitHub owner:** Claude Code session (model `claude-opus-5-5`), operated by `g4nesh`. Not a registered ledger agent.
- **Evidence level:** complete written proof. It has not been independently reviewed.
- **Exact artifact commit:** the commit that adds this file on branch `claude/mdo-minors-compound-certificates`.
- **Status and known gaps:**
  - Not found in the literature screened, but **Drury, "On a theorem of Wielandt and the compounds of unitary matrices", LAA 199 (1994) 391–412, was not read**. It proves a Birkhoff-type theorem for a restricted set of compound entries and may overlap.
  - Kovačec, LAA 201 (1994) is the special case where γ has multiplicity 1.

## Statement

Let n ≥ 1 and let A, B ∈ ℂ^{n×n} be normal with eigenvalues a₁,…,aₙ and b₁,…,bₙ (with multiplicity). Suppose B has **at most three** distinct eigenvalues, and if it has exactly three, **at least one of them is simple**. Then

  det(A+B) ∈ conv{ z_σ = ∏ᵢ (aᵢ + b_{σ(i)}) : σ ∈ Sₙ }.

More precisely, det(A+B) = Σ_σ λ_σ z_σ for a probability vector λ that depends on the relative unitary U but **not** on a, β, γ, δ.

The same holds with the roles of A and B exchanged, because both det(A+B) and the set {z_σ} are symmetric under A ↔ B.

## Definitions

- **Normal form.** Write A = V D_a V*, B = W D_b W*, U = V*W. Then det(A+B) = d(U) := det(D_a + U D_b U*); see [`../reductions/`](../reductions/README.md) §1.
- **Index sets.** The hypothesis is equivalent to the following: for some index m, the other n−1 eigenvalues take at most two values. Take δ := b_m to be a simple eigenvalue if B has three distinct eigenvalues, and any eigenvalue otherwise. Order the eigenvectors so that b_j = β for j ∈ J₀, b_j = γ for j ∈ J, and b_m = δ, with J₀ ⊔ J ⊔ {m} = [n] and k := |J| ≥ 0.
  - **Distinctness is never used.** β, γ, δ may coincide; the identities below are polynomial in them.
  - If the n−1 remaining eigenvalues take one value, any split works.
  - If J₀ = ∅, β is an arbitrary auxiliary complex number, used only as a translation.
- **Notation.** For X ⊆ [n], a^X := ∏_{i∈X} aᵢ (empty product = 1). M[S,T] is the submatrix on rows S and columns T.
- **Column blocks.** V := U[:,J] ∈ ℂ^{n×k} and V′ := U[:,J∪{m}] ∈ ℂ^{n×(k+1)}. Both have orthonormal columns.
- **Projection determinantal process.** For Y ∈ ℂ^{n×r} with orthonormal columns, P^Y is the probability measure on r-subsets of [n] with P^Y(X) = |det Y[X,:]|². It sums to 1 by Cauchy–Binet: Σ_X |det Y[X,:]|² = det(Y*Y) = 1. For r = 0 it is the point mass at ∅.

## Dependencies

| Lemma or external theorem | Exact statement / source | Required hypotheses | Why they hold here |
| --- | --- | --- | --- |
| Cauchy–Binet | det(PQ) = Σ_{\|K\|=r} det P[:,K] det Q[K,:] for P ∈ ℂ^{r×N}, Q ∈ ℂ^{N×r} | none | standard |
| Prop. 1 ([report](report.md) §2) | d(U) = Σ_{\|S\|=\|T\|} \|det U[S,T]\|² a^{[n]∖S} b^T | U unitary | proved in the report via Cauchy–Binet |
| Lemma 3 (below) | P^Y(R ⊆ 𝒳) = det(Y[R,:] Y[R,:]*) | Y has orthonormal columns | V, V′ are column blocks of a unitary |
| Lyons, Publ. Math. IHÉS 98 (2003), **Theorem 6.2** | E finite, H₁ ⊂ H₂ subspaces of ℓ²(E) ⇒ P^{H₁} ≼ P^{H₂} (stochastic domination for increasing events) | E finite; nested subspaces | E = [n]; range V ⊂ range V′. Lyons' P^H is determined by P^H(R ⊆ 𝒳) = det(P_H[R,R]), which by Lemma 3 equals our P^V for H = range V |
| Strassen's theorem, finite form | μ ≼ ν on the finite poset (2^E, ⊆) ⇒ there is a coupling (𝒳, 𝒴) with 𝒳 ~ μ, 𝒴 ~ ν, 𝒳 ⊆ 𝒴 a.s. | finite poset | 2^{[n]} is finite |

**Lemma 3.** If Y ∈ ℂ^{n×r} has orthonormal columns and R ⊆ [n], then Σ_{X ⊇ R, |X| = r} |det Y[X,:]|² = det(Y[R,:] Y[R,:]*).

*Proof.*
1. **Degenerate cases.** If |R| > r, or Y[R,:] has rank < |R|, both sides are 0. On the left, every such det Y[X,:] has linearly dependent rows.
2. **Reduce to a triangular block.** Otherwise take an LQ factorization Y[R,:] = [L 0] Q* with Q ∈ U(r) and L lower-triangular and invertible. Replace Y by YQ. This keeps orthonormal columns, multiplies every maximal minor by det Q (modulus 1), and leaves Y[R,:]Y[R,:]* unchanged.
3. **Factor each minor.** For X = R ⊔ Z, put the rows of R first. Then (YQ)[X,:] is block lower-triangular with diagonal blocks L and W[Z,:], where W := (YQ)[[n]∖R, |R|+1..r]. So |det Y[X,:]|² = |det L|²·|det W[Z,:]|².
4. **Sum over Z.** By Cauchy–Binet, Σ_Z |det W[Z,:]|² = det(W*W). Now W*W = I − ((YQ)[R, |R|+1..r])*((YQ)[R, |R|+1..r]) = I, because that block of YQ is 0.
5. **Conclude.** The sum is |det L|² = det(Y[R,:]Y[R,:]*). ∎

## Proof

1. **Translate.**
   - Put a′ = a + β·𝟙 and b′ = b − β·𝟙, so every product aᵢ + b_{σ(i)} = a′ᵢ + b′_{σ(i)} is unchanged.
   - Then b′ = c on J, e at m and 0 on J₀, where c = γ − β and e = δ − β.
   - d(U) is also unchanged, because D_a + U D_b U* = D_{a′} + U D_{b′} U*.
2. **Expand d(U).**
   - In Prop. 1, b′^T = 0 unless T ⊆ J ∪ {m}.
   - For a fixed S, split the T ⊆ J∪{m} with |T| = |S| by whether m ∈ T. Cauchy–Binet gives Σ_{T⊆J} |det U[S,T]|² = det(V[S,:]V[S,:]*) =: π(S).
   - Likewise the sum over all T ⊆ J∪{m} is π′(S) := det(V′[S,:]V′[S,:]*). So the T ∋ m part is π′(S) − π(S). Hence
   d(U) = Σ_{S⊆[n]} a′^{[n]∖S} ( c^{|S|} π(S) + c^{|S|−1} e (π′(S) − π(S)) ).
   - **Edge cases.** For S = ∅ the bracket is 1, since π(∅) = π′(∅) = 1. The second term is present only when m can lie in T, which forces |S| ≥ 1, so c^{|S|−1} is a polynomial.
3. **Hull points.** A permutation σ is determined, up to the irrelevant order inside blocks, by X = σ⁻¹(J) and j = σ⁻¹(m). For such a pair
   z_{X,j} = ∏_{i∈X}(a′ᵢ + c) · (a′ⱼ + e) · ∏_{i∉X∪{j}} a′ᵢ = Σ_{S ⊆ X∪{j}} a′^{[n]∖S} c^{|S∩X|} e^{|S∩{j}|}.
4. **Coupling.**
   - range V ⊆ range V′, so Lyons' Theorem 6.2 gives P^V ≼ P^{V′}. If k = n−1, P^{V′} is the point mass at [n], and the domination is trivial.
   - Strassen gives (𝒳, 𝒴) with 𝒳 ~ P^V, 𝒴 ~ P^{V′} and 𝒳 ⊆ 𝒴 a.s.
   - Since |𝒳| = k and |𝒴| = k+1 a.s., 𝒴 = 𝒳 ⊔ {𝒿} for a unique random index 𝒿.
5. **Take expectations.** By step 3,
   E z_{𝒳,𝒿} = Σ_S a′^{[n]∖S} ( c^{|S|} P(S ⊆ 𝒳) + c^{|S|−1} e P(S ⊆ 𝒴, 𝒿 ∈ S) ).
   - Lemma 3 gives P(S ⊆ 𝒳) = π(S).
   - "S ⊆ 𝒴 and 𝒿 ∉ S" is the same event as "S ⊆ 𝒳". So P(S ⊆ 𝒴, 𝒿 ∈ S) = P(S ⊆ 𝒴) − P(S ⊆ 𝒳) = π′(S) − π(S), again by Lemma 3.
   - Comparing with step 2, E z_{𝒳,𝒿} = d(U).
6. **Conclude.** Define λ_σ := P(𝒳 = σ⁻¹(J), 𝒿 = σ⁻¹(m)) / (|J|!·|J₀|!), spread uniformly over the permutations giving the same pair. This is a probability vector on Sₙ with Σ_σ λ_σ z_σ = d(U). It depends only on U, through the coupling, so d(U) lies in the closed convex hull. ∎

No inverse of A, B or A+B is used, so singular sums are included. Repeated or coincident values (β = γ, γ = δ, …) only make some z_σ coincide; the identity in steps 2–5 is polynomial and holds for all values.

## Computational assistance

None is needed for the proof. **Independent numerical check:** [`moc/verify_thmE.py`](moc/verify_thmE.py) (seed 42) solves for a coupling with exactly the marginals P^V and P^{V′} by LP. It then compares Σ P(X,Y) z with d(U) for random complex a, β, γ, δ and Haar U.

- Observed maximum relative error: ≤ 2·10⁻¹⁴ over 160 instances, for (n, k) ∈ {(4,1), (4,2), (5,2), (5,3), (6,2), (6,3), (6,4), (7,3)}.
- Environment: Python 3.9.6, numpy 2.0.2, scipy 1.13.1.

## Adversarial checks

- **Pattern LPs.** [`moc/pattern_lp.py`](moc/pattern_lp.py) found the a-independent certificate feasible in every trial for the patterns this theorem covers: (2,1,1), (3,1,1), (2,2,1), (3,2,1), (4,1,1), 120–300 Haar trials each.
- **Hypothesis is sharp in one direction.** Without a simple eigenvalue, but with three values, the statement is open (Conjecture F, [report](report.md) §4). With four values, a-independent certificates provably fail ([proposition-e.md](proposition-e.md)). So the "one simple eigenvalue" hypothesis cannot be replaced by "any three values" using this proof.
- **Edge cases checked in the argument:** k = 0 (the rank-one case, de Oliveira / Bebiano–Queiró), J₀ = ∅, and coincident values.

## Conclusion and limits

The theorem establishes MOC for every normal A paired with any normal B having at most three distinct eigenvalues, at least one of them simple when there are three, in every dimension n. (With at most two distinct eigenvalues this recovers Merikoski–Virtanen 1993.) It says nothing about B with three or more values besides the simple one, nor about multiplicity patterns (p,q,r) with all parts ≥ 2.

## Independent review

None yet. A reviewer under a different GitHub owner should check:

- Lemma 3;
- that Lyons' P^H coincides with P^V as defined here;
- the event identity in step 5;
- the bookkeeping in step 6.
