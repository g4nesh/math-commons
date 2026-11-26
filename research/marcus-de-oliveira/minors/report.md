# The Marcus–de Oliveira Conjecture: a progress report

Task `mdo-minors` (claimed in the live ledger on 25 Nov 2026, after the first-pass commit 78908ca; see [README](README.md)) · Author: `ganesh-claude-minors`, a Claude Code session (model `claude-opus-5-5`) operated by GitHub owner `g4nesh` · 24 Nov 2026 (America/Phoenix) · **revised 25 Nov 2026** after the independent review [`../reviews/mdo-minors.md`](../reviews/mdo-minors.md): bottom-line item 3, §5.2, §5.5, §5.6.

**Evidence level:** mixed and marked per item. An independent-owner review found the theorems and exact certificates correct and asked for changes to the n = 4 search summary; this revision makes them and has not itself been reviewed. All code and raw results are in [`moc/`](moc/).

**Bottom line.** The conjecture is **not solved** here, and no accepted resolution exists for any n ≥ 4. A claimed n = 4 proof (Zheng, October 2026) is recorded by [`../literature/`](../literature/README.md) as *claimed, audit pending*. Nothing below depends on that claim or audits it. This report contains:

1. **One new theorem**, to my knowledge. MOC holds whenever one of the two matrices has three distinct eigenvalues, one of which is simple (Theorem D). It contains Kovačec's 1994 "n−2 equal eigenvalues" theorem as a special case.
2. **A sharp dichotomy.** The "probabilistic-certificate" method is the mechanism behind almost all known exact cases. It always fails once one matrix has ≥ 4 distinct eigenvalues (Proposition E). I conjecture it always succeeds with ≤ 3 (Conjecture F), and give an equivalent reformulation as a statement about thinnings of determinantal point processes (Proposition G).
3. **A reduction of n = 4 to a small region, with numerical evidence at few points.** A counterexample for n = 4 can only occur for unitaries whose compound vector lies in one of 48 disjoint, symmetry-equivalent caps. Together these cover about 1.9% of U(4). The caps sit beyond the exotic facets of an explicit 13-dimensional polytope. Fixed-U searches (270,000 optimizations of the eigenvalues) at 271 distinct points of that region found no counterexample; 270 of the 540 unitaries used are symmetric copies of the deepest point. The values closest to zero, about −7·10⁻⁵ after re-optimization, occur at shallow points next to a repeated eigenvalue. No unitary was varied inside a cap, and no uniform margin is claimed.
4. **Targeted counterexample searches for n = 5, 6, 7.** None found.

Key computational claims were cross-checked by independent reruns, by 50-digit arithmetic, or against published values; each item says which. Each statement below is tagged:

- **[Known]**: published; I give the citation and re-derive it in a common framework.
- **[New]**: not found in the literature I could access (caveats given).
- **[Conjecture]**
- **[Numerical]**: evidence, not proof.

---

## 1. Statement and status

Let A, B be n×n normal matrices with eigenvalues a₁,…,aₙ and b₁,…,bₙ. For σ ∈ Sₙ put
  z_σ = ∏ᵢ (aᵢ + b_{σ(i)}).

> **MOC** (Marcus 1973, de Oliveira 1982): det(A+B) ∈ conv{ z_σ : σ ∈ Sₙ }.

WLOG A = D_a = diag(a) and B = U D_b U* with U unitary. Write Δ(U) = det(D_a + U D_b U*).

Status, from a literature survey done in this session (citations in §8):

- n = 2: Marcus.
- n = 3: Bebiano–Merikoski–da Providência, LMA 1987.
- n = 4: **no accepted proof**. It is stated to be open in Drury's ELA 2004 paper, in zbMATH reviews of the 2019 and 2025 papers, and in Mulero-Martínez, LAA 2026. Zheng's October 2026 preprint claims the full n = 4 case; the repository's literature audit records it as "claimed n = 4 proof; audit pending". My own survey missed this preprint.
- Known special cases:
  - Hermitian A, B (Fiedler 1971).
  - Spectra collinear or concyclic (Bebiano–Kovačec–da Providência 1994; Drury 2006/07).
  - One matrix with ≤ 2 distinct eigenvalues (Merikoski–Virtanen 1993).
  - n−2 equal eigenvalues (Kovačec 1994).
  - All eigenvalues real except three eigenvalues of B (Shitov 2025).
- No counterexample to the original conjecture has been published. Goyal's example, recorded in Kovačec's DMUC p2611 v2, refutes only a stronger polygon-enclosure conjecture; see [`../literature/`](../literature/README.md).

---

## 2. Framework

**Proposition 1 (Cauchy–Binet expansion) [Known; de Oliveira 1982].**
  Δ(U) = Σ_{k=0}^{n} Σ_{|S|=|T|=k} |det U[S,T]|² · a^{[n]∖S} · b^{T},
where a^{X} = ∏_{i∈X} aᵢ. With U = P_σ (so (P_σ)_{i,σ(i)} = 1) this gives z_σ, because |det P_σ[S,T]|² = [σ(S) = T].

*Proof.* Write D_a + U D_b U* = M N with M = [I | U] (n×2n) and N = diag(D_a, D_b)·[I; U*] (2n×n). Apply Cauchy–Binet to det(MN). A column set J = S′ ⊔ (n+T) gives the product ± det U[S′ᶜ,T] · a^{S′} b^T · (± conj det U[S′ᶜ,T]). Both signs come from the same Laplace expansion, so they agree. Set S = S′ᶜ. ∎
(Verified numerically to 5·10⁻¹⁵ for n ≤ 5: `moc/test_basic.py`.)

**Definition (certificate).** A probability vector λ on Sₙ is a *certificate for U* if
  Σ_σ λ_σ [σ(S) = T] = |det U[S,T]|² for all S, T.
By Proposition 1, a certificate gives Δ(U) = Σ λ_σ z_σ **for all a, b simultaneously**, so MOC holds for this U. Whether every U has a certificate is the Merikoski–Virtanen question (1989). Drury disproved it for n = 4 in 1992.

**Pattern certificates.** Suppose b is constant on the blocks J₀,…,J_{m−1} of a partition of [n] (a "multiplicity pattern"). Then Δ depends on U only through the grouped sums
  G_U(S, κ) = Σ_{T : |T∩J_j| = κ_j ∀j} |det U[S,T]|².
A pattern certificate is a probability λ on colourings κ : [n] → {0,…,m−1} with |κ⁻¹(j)| = |J_j| whose grouped statistics match. Equivalently:
  E_λ[∏_{i∈S} c_{κ(i)}] = det(Σ_j c_j Q_j)[S,S] for all S, as polynomials in c,
where Q_j = U Π_{J_j} U* is the eigenprojection onto the j-th eigenspace. A pattern certificate gives MOC for every A whenever B has that pattern.

**Proposition 2 (Haar average) [elementary].** E_U Δ(U) = (1/n!) Σ_σ z_σ, the centroid of the hull.
*Proof.* By symmetry, E|det U[S,T]|² = 1/C(n,k), because Σ_T |det U[S,T]|² = 1. ∎
This is consistent with the Marcus–Spielman–Srivastava quadrature formula for finite free convolution.

---

## 3. Proven cases, by one probabilistic method

Throughout, V ∈ ℂ^{n×k} has orthonormal columns. P^V denotes the *projection determinantal point process* with P^V(𝒳 = X) = |det V[X,:]|² for |X| = k.

**Lemma 3.** For R ⊆ [n]: Σ_{X ⊇ R, |X| = k} |det V[X,:]|² = det(V[R,:] V[R,:]*). That is, P^V(R ⊆ 𝒳) = det K[R,R] with K = VV*.

*Proof.*
1. If |R| > k, or V[R,:] is rank-deficient, both sides are 0.
2. Otherwise take an LQ factorization V[R,:] = [L 0]Q* with Q unitary, and replace V by VQ. This changes no modulus of a maximal minor and does not change V[R,:]V[R,:]*.
3. For X = R ⊔ Y the matrix (VQ)[X,:] is block lower-triangular, so |det| = |det L|·|det W[Y,:]| with W = (VQ)[Rᶜ, r+1..k].
4. Sum over Y. By Cauchy–Binet the sum is det(W*W). Now W*W = I − (VQ)[R, r+1..k]*(VQ)[R, r+1..k] = I, since that block of VQ is zero.
5. Hence the total is |det L|² = det(V_R V_R*). ∎

**Theorem A (n ≤ 3) [Known: Bebiano–Merikoski–da Providência 1987].**
For n = 3, Jacobi's complementary-minor identity gives |det U[S,T]| = |u_{s,t}| for |S| = |T| = 2, where s = Sᶜ and t = Tᶜ. Take any Birkhoff decomposition of the unistochastic matrix (|u_ij|²). It matches every layer of the compound vector, so it is a certificate. ∎
(LP check: 200/200 random U feasible.)

**Theorem B (two distinct eigenvalues) [Known: Merikoski–Virtanen 1993].**
Let b = β on Jᶜ and γ on J, with |J| = k. Then
  Δ = Σ_{|X|=k} |det U[X,J]|² · ∏_{i∈X}(aᵢ+γ) ∏_{i∉X}(aᵢ+β),
an explicit convex combination with weights given by the squared Plücker coordinates of the eigenspace.

*Proof.*
1. Shift a′ = a + β and c = γ − β.
2. By Proposition 1 and Cauchy–Binet, Δ = Σ_S a′^{Sᶜ} c^{|S|} det(V_S V_S*) with V = U[:,J].
3. Expanding the right-hand side gives Σ_R a′^{Rᶜ} c^{|R|} Σ_{X⊇R} |det V[X,:]|².
4. Apply Lemma 3. ∎

(Verified to 6.5·10⁻¹⁵ for all n ≤ 7 and all k: `moc/verify_misc.py`.)

**Theorem C (an eigenvalue of multiplicity ≥ n−2) [Known: Kovačec 1994].**
Shift so that b = (c₁, c₂, 0, …, 0). Proposition 1 then reduces the problem to the following. Find λ_{ij} ≥ 0 (i ≠ j) with row sums rᵢ = |u_{i1}|², column sums |u_{i2}|², and λ_{ij} + λ_{ji} = w_{ij} := |det U[{i,j},{1,2}]|².

This is a fractional orientation problem on Kₙ with edge weights w. By max-flow/min-cut (Gale–Hakimi), it is feasible iff r(Z) ≥ w(E(Z)) for every Z ⊆ [n]. By Cauchy–Binet,
  w(E(Z)) = det(V_Z* V_Z) = r(Z)s(Z) − |Σ_Z ū_{i1}u_{i2}|² ≤ r(Z)s(Z) ≤ r(Z).
The column sums follow from Lemma 3 with |R| = 1. ∎
*(A short self-contained flow proof. I have not seen Kovačec's text, which reportedly goes through restricted compounds, so the two may coincide.)*

**Theorem D (three distinct eigenvalues, one of them simple) [New, to my knowledge].**
Let B have eigenvalue β with multiplicity n−k−1, γ with multiplicity k ≥ 1, and a simple eigenvalue δ. Then MOC holds for (A, B) for every normal A. In fact Δ is a convex combination of the z_σ with weights independent of a.

*Proof.*
1. **Setup.** Let J be the γ-columns and m the δ-column. Put V = U[:,J], V′ = U[:,J∪{m}], c = γ−β, e = δ−β, a′ = a+β.
2. **Expand Δ.** Proposition 1 and Cauchy–Binet give
  Δ = Σ_S a′^{Sᶜ} [ c^{|S|} π(S) + c^{|S|−1} e (π′(S) − π(S)) ],
  with π(S) = det(V_S V_S*) and π′(S) = det(V′_S V′_S*). For S = ∅ the bracket is 1, since π′(∅) = π(∅) = 1.
3. **The hull points.** The hull point for positions X (receiving γ) and j (receiving δ) is
  z_{X,j} = Σ_S a′^{Sᶜ} c^{|S∩X|} e^{|S∩{j}|} [S ⊆ X∪{j}].
4. **The coupling.** range V ⊂ range V′, so by **Lyons, Thm 6.2** (Publ. IHÉS 98, 2003) P^V is stochastically dominated by P^{V′}. By Strassen's theorem there is a coupling (𝒳, 𝒴) with 𝒳 ~ P^V, 𝒴 ~ P^{V′} and 𝒳 ⊆ 𝒴. Cardinalities are k and k+1, so 𝒴 = 𝒳 ∪ {j} for a unique random j.
5. **Match the coefficients.** For each S:
   - P(S ⊆ 𝒳) = π(S), by Lemma 3.
   - P(S ⊆ 𝒴, j ∈ S) = P(S ⊆ 𝒴) − P(S ⊆ 𝒳) = π′(S) − π(S), because "S ⊆ 𝒴 and j ∉ S" is the same event as "S ⊆ 𝒳".
   Hence E[z_{𝒳,j}] = Δ. ∎

- *Scope.* It contains Kovačec's theorem (the case k = 1). It also gives new cases such as multiplicities (2,2,1) for n = 5 and (3,2,1) or (2,3,1) for n = 6. MOC is symmetric in A ↔ B, so either matrix may carry the pattern.
- *Verification.* I solved the coupling LP with exactly Lyons' marginals and checked Δ = Σ P(X,Y) z: relative error ≤ 2·10⁻¹⁴ over 160 random instances, n = 4…7 (`moc/verify_thmE.py`).
- *Novelty caveat.* The survey found no such statement. Zheng's 2026 preprint (Theorem 1.2) also works with squared minors for n = 4; I have not read it. However, Drury, "On a theorem of Wielandt and the compounds of unitary matrices" (LAA 199, 1994), proves a Birkhoff-type theorem for "a certain restricted subset" of compound entries. I could not access its full text, so an overlap cannot be ruled out.

**Proposition E (≥ 4 distinct eigenvalues ⇒ certificates fail) [New argument, building on Drury 1992].**
Let a multiplicity pattern have ≥ 4 parts. Then some U ∈ U(n) admits no pattern certificate.

*Proof.*
1. Take Drury's U₄ ∈ U(4), which has no certificate. Embed it on rows 1–4 and on one column from each of four distinct blocks; put a permutation elsewhere.
2. For S ⊆ {1,…,4}, every nonzero minor has T inside those four columns.
3. The statistic for S = {1,…,4} forces the colouring restricted to those rows to be a bijection onto the four blocks, almost surely.
4. Its lower statistics then form a certificate for U₄, a contradiction. ∎
Numerically, the LP gap is the same 1.686·10⁻² for the patterns (1,1,1,1), (2,1,1,1), (2,2,1,1), (3,1,1,1) and (1⁵) (`moc/verify_misc.py`).

---

## 4. The three-eigenvalue conjecture

> **Conjecture F.** Every pattern with at most three parts admits a pattern certificate for every U. Consequently, MOC holds whenever A or B has at most three distinct eigenvalues.

With Theorems B and D and Proposition E, this would say: *the certificate method works exactly up to three distinct eigenvalues.* The first open pattern is (2,2,2), n = 6.

**Proposition G (thinning reformulation) [New].**
A pattern certificate for (Q₀, Q₁, Q₂) exists iff there is a random pair of disjoint sets (𝒳₁, 𝒳₂) with the following property. Thin 𝒳₁ independently with retention s and 𝒳₂ with retention t. The union is then the determinantal process DPP(sQ₁ + tQ₂), for **every** s, t ∈ [0,1].

*Proof.*
1. Shift c₀ to 0. Both sides are homogeneous of degree |S| in (c₁, c₂).
2. The left side evaluated at (s,t) is P(S ⊆ 𝒳₁^{(s)} ∪ 𝒳₂^{(t)}). The right side is det(sQ₁+tQ₂)[S,S].
3. A polynomial identity on [0,1]² is a polynomial identity everywhere.
4. Inclusion probabilities determine a point process on a finite set. ∎

For two types, independent thinning of a single DPP is classical: thinning P^K gives P^{tK}. The conjecture asks for the two-type analogue.

**Remark (why a canonical formula cannot work) [New observation].**
The natural determinantal candidate is the generating polynomial det(I + (D_y − I)Q₁ + (D_w − I)Q₂). It has the correct statistics. Its coefficients are, up to sign, ∏_j det U[X_j, J_j] / det U. These are complex in general, and their real parts reduce to the Merikoski–Virtanen signed weights for (1,…,1). This suggests that any proof of Conjecture F will be non-canonical, in the style of Birkhoff or Strassen.

**Evidence [Numerical].** Pattern certificates were tested by LP. "Adversarial" means a ray-shooting depth s*(U) minimized over U; s* < 1 means no certificate.

| pattern | parts | random LP trials | adversarial min s* | status |
|---|---|---|---|---|
| (2,2), (3,1), (3,2) | 2 | 300 each, all feasible | — | Theorem B |
| (2,1,1), (3,1,1), (2,2,1) | 3 | 300 each, all feasible | — | Theorem D |
| (3,2,1), (4,1,1) | 3 | 120 each, all feasible | — | Theorem D |
| **(2,2,2)** | 3 | 120, all feasible | **1.000000** (80 restarts) | Conjecture F |
| **(3,2,2)** | 3 | — | **1.000000** (16 restarts) | Conjecture F |
| (1,1,1,1) | 4 | 4/300 infeasible | 0.881 | fails (Drury) |
| (2,1,1,1) | 4 | 10/300 infeasible | 0.856 | fails (Prop. E) |
| (2,2,1,1) | 4 | 1/120 infeasible | 0.798 | fails (Prop. E) |
| (1,1,1,1,1) | 5 | 44/300 infeasible | — | fails |

The same adversarial code finds failures for each 4-part pattern tested within 8 restarts. It never found one for a 3-part pattern (96 restarts in total).

A strictly stronger statement also survived adversarial testing. Replace the eigenprojections by an arbitrary 3-outcome **POVM** (E₀+E₁+E₂ = I, E_j ⪰ 0). Then min s* = 1.000000 for n = 3, 4, 5 with various ranks. The 4-outcome projection case (n = 4) gives 0.800. This version is closed under compression, which suggests an inductive proof might exist (`moc/povm_*.py`).

---

## 5. The case n = 4

**5.1 The compound polytope [Numerical; facet count matches Kovačec–Ribeiro, arXiv:1604.06112].**
The vectors p_σ = ([σ(S)=T])_{S,T} span a 13-dimensional affine space. Their convex hull Q₄ has 24 vertices and 82 facets, in three types:

- **16 facets:** |u_ij|² ≥ 0.
- **18 facets:** |det U[S,T]|² ≥ 0 for 2×2 minors. Jacobi's identity pairs the 36 minors into 18.
- **48 exotic, simplicial facets** (13 vertices each). They form a **single orbit** under left/right S₄ × transposition, and are the only facets a unitary can violate. One representative, in 1-based indices:

  **g(U) = |u₁₁|² + |u₃₃|² + |det U[{1,4},{3,4}]|² − |det U[{1,2},{1,4}]|² ≥ 0.**

  On vertices, g takes the values 0 (13 times), 1 (10 times) and 2 (once).

- **Cross-check with Zheng's Example 1.3.** Zheng's functional ψ, as transcribed in [`../literature/`](../literature/README.md), coincides with one of these 48 facets on all 24 vertices, hence on the whole affine span. Our facet code independently reproduces ψ(O*) = −158/2401 ≈ −0.065806 for his rational orthogonal O* (`moc/crosscheck_zheng.py`). The complex-unitary optimum in §5.3 goes slightly deeper (−0.067442).

All facet-normal (direction) computations in §5 are done in the 13 affine coordinates (M1[:3,:3], N[:2,:2]) used by `moc/reduced4.py`.

**5.2 Reduction lemma [Rigorous, given the facet list of Q₄].**
If MOC fails at (U, a, b) with n = 4, then g(P_α U P_β) < 0, or the same holds with U replaced by U*, for some permutations α, β.
*Proof.* If p_U ∈ Q₄ there is a certificate. The trivial facets always hold. The symmetries (U,a,b) ↦ (P_αUP_β, permuted a, permuted b) and (U,a,b) ↦ (U*, b, a) preserve Δ and the hull, and they act transitively on the exotic facets. ∎
*Dependence.* The proof needs the facet list of §5.1 to be complete. `moc/facets4.py` obtains it from a floating-point convex hull. It has since been confirmed in exact integer arithmetic: 82 facets from the orientation determinants of all 1,961,256 vertex 14-subsets ([`../star-halfline/q_facets_exact.py`](../star-halfline/q_facets_exact.py)).

**5.3 Geometry of the bad region [Numerical].**
- **Size and shape (400,000 Haar samples):**
  - The bad region has measure ≈ 1.91%.
  - No sample violated two exotic facets at once, so the caps are disjoint.
  - In the sample, every bad p_U lies at distance ≥ 0.371 from all vertices (diameter of Q₄ = 3.0). So, empirically, there is no near-permutation subtlety.
  - The projections of bad points onto the violated facet lie inside the facet simplex for 95.6% of samples, and within 0.007 of it for all; every barycentric coordinate is ≤ 0.743.
- **Deepest violation:** min g = −0.067442. Equivalently, Drury's functional 1 − g reaches **1.067442**, reproducing Drury's 1992 value. The extremal U has a striking structure: |U|² has three entries 0.32390 and one entry 0.02830 in every row, and the four small entries form a permutation pattern.

**5.4 Why compound failure does not give an MOC failure [Numerical].**
- A MOC functional Re(e^{−iθ}·Δ) has vertex values v_σ = Re(e^{−iθ}z_σ(a,b)). Over all (a,b,θ), the largest correlation of these centred vertex values with the exotic-facet normal is **0.95197164**. This was found identically by 1,000 restarts with two optimizers and log-uniform scales.
- So no MOC functional points along the facet normal. A functional aligned with it would need 13 of the 24 points z_σ on one supporting line.

**5.5 Targeted search [Numerical].** *Revised after the independent review ([`../reviews/mdo-minors.md`](../reviews/mdo-minors.md), finding 1); re-analysis in `moc/reanalyse4.py`.*
- **What was run:** 540 unitaries × 500 eigenvalue starts each = 270,000 optimizations of (a, b, θ) **with U held fixed**. The objective is
  J = (Re e^{−iθ}Δ − max_σ Re e^{−iθ}z_σ)/width, with analytic gradients verified to 10⁻¹⁰.
- **What the 540 unitaries are:** only **271 distinct cap points** up to symmetry.
  - The 270 "deepened" unitaries all converged to the global extremal orbit: their symmetry-invariant depth (minimum of g over P_αUP_β and transposes) is −0.067442 for every one, and their sorted |u_ij|² profiles agree to 1.5·10⁻⁵.
  - The other 270 are raw Haar violators, with invariant depth between −0.0436 and −0.0001.
  - No search varies U inside a cap. The docstring of `moc/bad4.py` used to announce a joint refinement of U; it was never implemented.
- **Result: no counterexample in these runs. Every J is negative.**
  - At the extremal point: J between −0.0287 and −0.0280 over its 270 copies.
  - Over the raw points: the largest stored value is J = −2.59·10⁻⁴. **It is not a converged maximum.** Re-optimizing the three best records at the same fixed U (restarted Nelder–Mead) gives:

    | stored J | re-optimized J | invariant depth of U | note |
    | --- | --- | --- | --- |
    | −2.59·10⁻⁴ | −7.7·10⁻⁵ | −0.0020 | two aᵢ coalesce (relative gap 2·10⁻⁷) |
    | −2.79·10⁻⁴ | −6.6·10⁻⁵ | −0.0011 | two aᵢ nearly coalesce (9·10⁻⁶) |
    | −2.88·10⁻⁴ | −1.0·10⁻⁴ | −0.0062 | two bᵢ coalesce (1·10⁻⁶) |

    All three stay negative, and all three sit next to a proven case (a repeated eigenvalue).
- **Withdrawn:** the earlier sentence "corr(J, violation depth) = −0.93; deeper violations give more negative J". The correlation is produced by the two clusters (copies of the extremal point against shallow points). Within the 270 raw points it is −0.17 (−0.18 with the invariant depth). The stored violation is a Euclidean distance in non-symmetric coordinates and is not symmetry-invariant: it ranges from 0.0294 to 0.0389 over the 270 copies of the one extremal point.
- **One line experiment (`moc/margin4.py`), with its limits:**
  - Push an artificial point outward from the facet along the facet normal through one copy of the deepest violator (x_t = foot + t·normal, in the reduced coordinates) and replace Δ by Σμ_σ(x_t)z_σ.
  - The optimizer detects a "fake counterexample" once t ≳ 0.065. That copy of the true point sits at t = 0.0389, a ratio of about 1.66 on this line.
  - This shows that the search finds positive J readily when positive J exists. It is **not** an invariant safety margin: the depth 0.0389 depends on which symmetric copy and which coordinates are used.

**5.6 A route to a computer-assisted proof of n = 4 (not carried out).**
By 5.2 it suffices to show, for the single facet g:
  max_{U ∈ cap_g} Re e^{−iθ}Δ ≤ max_σ Re e^{−iθ}z_σ for all (a,b,θ).
What the data do and do not show:
- The cap is a compact semialgebraic set of dimension ≤ 9, away from the vertices.
- At the deepest point of the cap the shortfall is clear (J ≈ −0.028). At shallow points it is small (down to about −7·10⁻⁵), next to configurations with a repeated eigenvalue, which are proven cases.
- The earlier claim of a uniform shortfall with a 1.66× margin is **withdrawn**: fixed-U searches at 271 points cannot support it.

Making this rigorous would need certified enclosures of the cap and an interval branch-and-bound over (a,b,θ) (13 real dimensions after symmetries), together with a separate local argument near the degenerate limits. For the case where one matrix is Hermitian, a different and exact route has since been carried out for part of the spectra: see [`../pmsh/`](../pmsh/README.md), [`../star-degree4/`](../star-degree4/README.md) and [`../star-halfline/`](../star-halfline/README.md).

---

## 6. Searches for n = 5, 6, 7 [Numerical]

- **Bad-region fractions:** compound certificates fail for ≈ 16% of Haar U(5) (200 of 1,249 samples) and ≈ **80%** of Haar U(6) (36 of 45 samples, so a rough estimate). For large n the certificate method is the exception, not the rule.
- **Targeted searches** (optimize (a,b,θ) at bad U):

  | n | unitaries | starts per U | max J |
  |---|---|---|---|
  | 5 | 200 | 400 | −3.06·10⁻⁴ |
  | 6 | 36 | 150 | −6.80·10⁻⁴ |
  | 7 | 18 (Haar) | 60 | −1.33·10⁻³ |

- **Unconstrained joint searches** (60k / 15k / 3k / 150 starts for n = 4/5/6/7) all converged to U within 10⁻¹³ of a permutation matrix, i.e. Δ at a hull vertex. At 50-digit precision J ≤ 0 for all of them; the apparent J ≈ +10⁻¹² values were floating-point roundoff (`moc/verify_best.py`).

---

## 7. What would move this further

1. **Prove Conjecture F**, perhaps via Proposition G. A two-type "Strassen for DPP thinnings" theorem would give MOC whenever one matrix has ≤ 3 distinct eigenvalues.
2. **Carry out §5.6.** Of everything here, this is the most plausible route to n = 4.
3. **Check Drury 1994** for overlap with Theorem D.
4. **Relation to Zheng's claim.** If Zheng's n = 4 proof survives audit, §5 is consistent with it but adds nothing to its rigor. Conversely, §5.2–5.5 may help an auditor see where the hard region lies: the 48 caps and the deepest violator.

---

## 8. Key references

These were checked by the survey agent at abstract, zbMATH-review or full-text level, as indicated in its report.

- M. Fiedler, *Bounds for the determinant of the sum of Hermitian matrices*, Proc. AMS 30 (1971).
- G. N. de Oliveira, *Normal matrices (research problem)*, LMA 12 (1982).
- N. Bebiano, J. K. Merikoski, J. da Providência, *On a conjecture of G. N. de Oliveira on determinants*, LMA 20 (1987) — n = 3.
- J. K. Merikoski, A. Virtanen, LAA 121 (1989); LAA 187 (1993) — two distinct eigenvalues.
- S. W. Drury, *A counterexample to a question of Merikoski and Virtanen on the compounds of unitary matrices*, LAA 168 (1992).
- A. Kovačec, *On a conjecture of Marcus and de Oliveira*, LAA 201 (1994) — n−2 equal eigenvalues.
- S. W. Drury, *On a theorem of Wielandt and the compounds of unitary matrices*, LAA 199 (1994).
- N. Bebiano, A. Kovačec, J. da Providência, LAA 197/198 (1994); S. W. Drury, ELA 15 (2006) — essentially Hermitian / circle–line cases.
- S. W. Drury, *The external vertices conjecture in case n = 4*, ELA 11 (2004).
- A. Kovačec, P. Ribeiro, arXiv:1604.06112 — 82 facets of the n = 4 compound polytope.
- Y. Shitov, Izv. Math. 89 (2025) — all real except three eigenvalues of B.
- N. Bebiano, J. da Providência, *Revisiting the Marcus–de Oliveira Conjecture*, Mathematics 13 (2025) — survey.
- R. Lyons, *Determinantal probability measures*, Publ. Math. IHÉS 98 (2003), Thm 6.2 — https://arxiv.org/abs/math/0204325
- D. Zheng, *The Marcus–de Oliveira conjecture for n = 4*, preprint, October 2026, doi:10.13140/RG.2.2.36790.25923 — **claimed** proof; not read by me; audit status in [`../literature/`](../literature/README.md).
- A. Kovačec, DMUC preprint p2611 v2 (4 Nov 2026) — stronger polygon conjecture and Goyal's counterexample to it; see [`../literature/`](../literature/README.md).

---

## 9. Reproducing

Requires `numpy`, `scipy` and `mpmath`. Run from `research/marcus-de-oliveira/minors/`, e.g. `python3 moc/test_basic.py`. Tested with Python 3.9.6, numpy 2.0.2, scipy 1.13.1 (HiGHS LP solver) and mpmath 1.3.0 on macOS (Darwin 24.6.0, 11 cores). Exact commands, seeds and compute are listed in [README.md](README.md).

| file | purpose |
|---|---|
| `moclib.py` | core: Haar sampling, Δ, hull points, Cauchy–Binet, certificate LP |
| `test_basic.py` | Prop. 1 check; certificate LP rates for n = 3, 4, 5 |
| `reduced4.py`, `facets4*.py`, `exotic4.py` | n = 4 polytope, 82 facets, exotic facet formula |
| `orbit_align.py`, `align4.py` | single-orbit check; 0.95197 alignment bound |
| `crosscheck_zheng.py` | ψ from Zheng's Example 1.3 is one of the 48 exotic facets; reproduces −158/2401 |
| `geom4.py`, `bad4.py`, `analyse4.py`, `margin4*.py` | bad-region geometry, 270k targeted search, depth correlation, safety margin |
| `bsearch.py`, `gradcheck.py`, `verify_best.py` | batched analytic-gradient search; gradient check; 50-digit re-verification |
| `badn.py`, `fixedU.py` | targeted searches, n = 5, 6, 7 |
| `pattern_lp.py`, `pattern_adv.py` | pattern certificates (random and adversarial) |
| `povm_lp.py`, `povm_adv.py` | POVM generalization |
| `verify_thmE.py`, `verify_misc.py` | Theorem D construction; Theorem B weights; Prop. E embedding |
| `runs/` | raw JSON / NPY results and logs |
