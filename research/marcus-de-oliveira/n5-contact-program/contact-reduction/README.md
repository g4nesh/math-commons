# The contact reduction in every dimension, for a fixed unitary matrix

- **Task:** `mdo-contact-reduction-all-n`. Agent `ishaan-n5` (GitHub owner `IshaanAyaan`). November 26, 2026.
- **Status:** complete written proof. Textbook facts are cited by name. The identities were checked exactly and the deformation was run numerically at n = 4 and n = 5 (see "Checks"). One independent re-derivation by a fresh agent under the same owner. **Not independently reviewed by another owner.**
- **Evidence levels:** PROVED (a complete written argument here) and CITED (a named textbook theorem), marked per step. The checks are EXACT or NUMERIC as labelled.

## Deviations and limits

- The statement is proved for *admissible weights* (defined below), which include the Leibniz weights of every unitary matrix. This is more general than the task asked. The unitary case is Corollary A1.
- Classical facts from algebraic geometry, complex analysis and combinatorics are cited, not re-proved: the affine dimension theorem, a curve through two points, normalization and completion of curves, Chevalley's theorem, density of Zariski-open sets in the classical topology, dimension as transcendence degree, the Laplace expansion, Jacobi's complementary-minor identity, Hall's theorem and the Birkhoff–von Neumann theorem. Everything else is proved here. In particular, the closedness of the realizable cone (Lemma D) and Zheng's Lemma 3.4 for every n (Lemma B) need neither geometric invariant theory nor any description of the ideal of the cone.
- Nothing here addresses the case analysis that would follow the reduction at n = 5. That is the subject of the tasks `mdo-n5-contact-strata` and `mdo-n5-high-contact-search`.

## Statement

**Notation.** n ≥ 2. For a, b ∈ ℂⁿ and σ ∈ Sₙ put z_σ(a,b) = ∏ᵢ (aᵢ + b_σ(i)). For U ∈ U(n) put d_U(a,b) = det(D_a + U D_b U*), where D_a = diag(a). The conjecture for U says that d_U(a,b) ∈ conv{z_σ(a,b) : σ ∈ Sₙ} for all a, b. P_σ is the permutation matrix with (P_σ)ᵢⱼ = [σ(i) = j].

**Admissible weights.** μ ∈ ℝ^{Sₙ} is *admissible* if ∑_σ μ_σ = 1 and the n × n matrix M(μ)ᵢⱼ = ∑_{σ: σ(i)=j} μ_σ has nonnegative entries. Its row and column sums are then 1, so M(μ) is doubly stochastic. Put D_μ(a,b) = ∑_σ μ_σ z_σ(a,b).

> **Lemma 1 (Leibniz weights of a unitary).** For U ∈ U(n) put λ_σ(U) = sgn(σ) ∏ᵢ U_{i,σ(i)} · conj(det U) and μ_σ(U) = Re λ_σ(U).
> (a) For all K, J ⊆ {1,…,n} with |K| = |J|: |det U[K,J]|² = ∑_{σ: σ(K)=J} λ_σ(U). So the imaginary parts of these λ_σ(U) cancel. The empty minor is 1.
> (b) d_U(a,b) = ∑_σ λ_σ(U) z_σ(a,b) = ∑_σ μ_σ(U) z_σ(a,b) for all a, b.
> (c) μ(U) is admissible, M(μ(U)) = (|Uᵢⱼ|²), and |μ_σ(U)| ≤ ∏ᵢ |U_{i,σ(i)}| ≤ 1.

So the conjecture for U is the statement "D_μ(a,b) ∈ conv{z_σ(a,b)} for all a, b" with μ = μ(U).

> **Theorem A (contact reduction).** Let μ be admissible.
> (a) If D_μ(a,b) ∉ conv{z_σ(a,b)} for some a, b ∈ ℂⁿ, then there are a′, b′ ∈ ℂⁿ and θ ∈ ℂ with |θ| = 1 such that, with m = min_σ Re(θ z_σ(a′,b′)):
>   Re(θ D_μ(a′,b′)) < m and #{σ : Re(θ z_σ(a′,b′)) = m} ≥ 2n − 2.
> (b) (Variational and rigid form.) In the homogeneous setting of Section 2, let x ∈ X with S_μ(x) false, π ∈ T(x), the gaps normalized by ∑_σ Re g_σ(x) = 1, and L(y) = ∑_σ μ_σ Re g_σ(y). Then some x′ ∈ X has T(x′) ⊇ T(x), ∑_σ Re g_σ(x′) = 1, L(x′) ≤ L(x) < 0, |T(x′)| ≥ 2n − 2, and x′ is an isolated point of X ∩ Λ(x′).

> **Corollary A1 (the unitary case).** Let n ≥ 2 and U ∈ U(n). If the conjecture fails for U at some spectra, it fails for the same U at spectra (a′, b′) where some supporting line of conv{z_σ(a′,b′)} strictly separates d_U(a′,b′) from the hull and contains z_σ(a′,b′) for at least 2n − 2 permutations σ, counted by permutation.

*Proof of Corollary A1.* Apply Theorem A(a) to μ(U), which is admissible by Lemma 1(c). Then D_μ = d_U by Lemma 1(b), and {w : Re(θw) = m} is the line. ∎

**Remarks.**
1. For n = 4 Corollary A1 is Zheng's Theorem 3.6 with U(4) in place of SO(4): six contacts. For n = 5 it gives eight contacts for every U ∈ U(5) directly. No reduction to orthogonal matrices is used. Such a reduction is not available for n = 5 if the separation claimed in task `mdo-orthogonal-hull-five` (not reviewed here) holds.
2. Lemma B is Zheng's Lemma 3.4 for every n. For n = 4 Zheng proves it with the Hessian of the matching cubic. For n = 5 the audit ([`generalization-n5.md`](../../zheng-n4-audit/generalization-n5.md), §2) checked a version modulo a prime, from 300 quadrics. The proof here uses only the binomial relations ∏F_σ^{c⁺_σ} = ∏F_σ^{c⁻_σ} and is uniform in n. It does not decide the audit's question "is K_n = E_n?", which is not needed.
3. Products are counted by permutation. Repeated eigenvalues, a zero factor aᵢ + bⱼ = 0 and collinear spectra are all allowed at the output point.
4. The deformation never increases the normalized support gap (part (b)). So a numerical search for a violation may restrict to rigid configurations with at least 2n − 2 contacts without losing the infimum.

## 1. Proof of Lemma 1 (PROVED, citing two classical identities)

Write ε_{K,J} = (−1)^{∑K + ∑J}, Kᶜ for the complement, and d = det U.

- **Laplace (CITED, e.g. Horn–Johnson, *Matrix Analysis*, §0.8).** For |K| = |J|, grouping the Leibniz terms by σ(K) = J gives det U[K,J] · det U[Kᶜ,Jᶜ] = ε_{K,J} ∑_{σ(K)=J} sgn(σ) ∏ᵢ U_{i,σ(i)}.
- **Jacobi (CITED, same §0.8).** For invertible U, det(U⁻¹[J,K]) = ε_{K,J} det U[Kᶜ,Jᶜ] / d. For unitary U, U⁻¹ = U*, so det(U⁻¹[J,K]) = conj(det U[K,J]), and 1/d = conj(d). Hence conj(det U[K,J]) = ε_{K,J} · det U[Kᶜ,Jᶜ] · conj(d).

(a) Multiply the second identity by det U[K,J] and use the first: |det U[K,J]|² = ε_{K,J}² ∑_{σ(K)=J} sgn(σ) ∏ U_{i,σ(i)} · conj(d) = ∑_{σ(K)=J} λ_σ(U). The left side is real, so the imaginary parts cancel. For K = J = ∅ this reads ∑_σ λ_σ(U) = d · conj(d) = 1.

(b) The matrix D_aU + UD_b has entries U_{ij}(aᵢ + bⱼ), so by the Leibniz formula det(D_aU + UD_b) = ∑_σ sgn(σ) ∏ᵢ U_{i,σ(i)} z_σ(a,b). Since D_aU + UD_b = (D_a + UD_bU*)U, dividing by d (that is, multiplying by conj(d)) gives d_U = ∑_σ λ_σ z_σ. Next expand z_σ(a,b) = ∑_{K} a^{Kᶜ} b^{σ(K)}, where a^S = ∏_{i∈S} aᵢ; here K is the set of factors in which b is chosen. Then ∑_σ (Im λ_σ) z_σ = ∑_{K,J} a^{Kᶜ} bᴶ ∑_{σ(K)=J} Im λ_σ = 0 by (a). So ∑ λ_σ z_σ = ∑ μ_σ z_σ.

(c) ∑_σ μ_σ = Re 1 = 1. M(μ(U))ᵢⱼ = ∑_{σ(i)=j} μ_σ(U) = |Uᵢⱼ|² ≥ 0 by (a) with K = {i}, J = {j}. Finally |λ_σ| = ∏ᵢ |U_{i,σ(i)}| · |d| ≤ 1. ∎

As a by-product, (a) gives the vector of squared minors of any unitary as an explicit real combination ∑_σ μ_σ(U) · (indicator of σ(K) = J) of the vectors of the permutation matrices, with weights summing to 1. This is the "equal spans" statement of the audit by a direct formula.

## 2. The realizable cone

**Homogeneous coordinates.** Y = (ℂ²)^{2n}, with points v = (v₁,…,vₙ, w₁,…,wₙ). For u, u′ ∈ ℂ² put [u,u′] = det(u u′) = u₁u′₂ − u₂u′₁. Put F_σ(v) = ∏ᵢ [vᵢ, w_σ(i)], Φ(v) = (F_σ(v))_σ ∈ ℂ^{Sₙ}, and X = Φ(Y). The chart ι(a,b) = ((aᵢ,1)ᵢ, (−bⱼ,1)ⱼ) has [vᵢ,wⱼ] = aᵢ + bⱼ, so Φ(ι(a,b)) = z(a,b).

For y ∈ ℂ^{Sₙ} put D_μ(y) = ∑_σ μ_σ y_σ (linear) and T(y) = {σ : Re y_σ = min_τ Re y_τ}, the *contact set*. **S_μ(y)** is the statement Re D_μ(y) ≥ min_σ Re y_σ. For a reference π ∈ T(y), the *gaps* are g_σ(y) = y_σ − y_π. Then Re g_σ(y) ≥ 0, with equality exactly on T(y), and since ∑μ_σ = 1,

  Re D_μ(y) − Re y_π = ∑_σ μ_σ Re g_σ(y) =: L(y).

So S_μ(y) fails exactly when L(y) < 0. For x ∈ X and π ∈ T(x) let Λ(x) = {y : g_σ(y) = g_σ(x) for σ ∈ T(x), and ∑_σ g_σ(y) = ∑_σ g_σ(x)}. These are |T(x)| affine-linear conditions, since the condition for σ = π is empty. Λ(x) does not depend on the choice of π ∈ T(x), because changing π subtracts the fixed number g_{π′}(x) from every gap.

**Basic properties (PROVED).**
- (X1) X is a cone: scaling v₁ by t ∈ ℂ multiplies Φ by t.
- (X2) 𝟙 = (1,…,1) ∈ X: take vᵢ = (1,0) and wⱼ = (0,1).
- (X3) Φ is invariant under G = SL₂(ℂ) × T₀, where T₀ = {t ∈ (ℂ*)^{2n} : ∏ t_k = 1} and (g, t)·v = (t_k g v_k)_k. Brackets are SL₂-invariant, and each F_σ is linear in each of the 2n vectors.
- (X4) **Chart lemma.** If x = Φ(v) ≠ 0, then x = λ z(a,b) for some λ ∈ ℂ*, a, b ∈ ℂⁿ. *Proof.* If some vᵢ or wⱼ were 0, every F_σ would have a zero factor, so all vectors are nonzero. Choose a nonzero e ∈ ℂ² on none of the 2n lines ℂvᵢ, ℂwⱼ, and g ∈ SL₂ mapping e to a multiple of (1,0). Then every g vᵢ and g wⱼ has nonzero second coordinate. Write g vᵢ = qᵢ(aᵢ,1) and g wⱼ = tⱼ(−bⱼ,1). By (X3), x = Φ(gv) = (∏qᵢ ∏tⱼ) z(a,b). ∎
- (X5) S_μ(0) holds (0 ≥ 0), so every point of X at which S_μ fails is nonzero.

> **Lemma D (closedness).** X = Φ(Y) is closed in the Euclidean topology. Consequently X is an irreducible Zariski-closed subset of ℂ^{Sₙ}.

*Proof (PROVED, except the last sentence, which is CITED).* Let N = {u ∈ Y : |u₁| = ⋯ = |u_{2n}| =: r and ∑_k u_k u_k* = n r² I₂}. Here u_k runs over all 2n vectors and |·| is the Euclidean norm.

1. *Φ(Y) = Φ(N).* Let v ∈ Y. If Φ(v) = 0 take u = 0 ∈ N. Otherwise let K be the Euclidean closure of the orbit G·v. It is closed and G-invariant, and Φ ≡ Φ(v) on K by continuity and (X3), so 0 ∉ K. The function ‖·‖² is proper, so it attains its minimum on K at some u ≠ 0. For ξ ∈ sl₂(ℂ) and s ∈ ℝ^{2n} with ∑s_k = 0, the curve γ(t) = (exp(tξ), (e^{ts_k})_k)·u stays in K. So f(t) = ‖γ(t)‖² = ∑_k e^{2ts_k} |exp(tξ)u_k|² is minimal at t = 0, and

   0 = f′(0) = 2∑_k s_k |u_k|² + 2 Re tr(ξ M), with M = ∑_k u_k u_k*.

   Taking ξ = 0 gives that all |u_k| are equal, say to r > 0. Taking s = 0 and ξ = M − (tr M/2)I (Hermitian, trace zero) gives 0 = tr(ξM) = ‖M − (tr M/2)I‖², so M = (tr M/2) I = n r² I. Hence u ∈ N and Φ(u) = Φ(v).
2. *Φ ≠ 0 on N₁ = N ∩ {r = 1}.* Let u ∈ N₁ and e a unit vector. Then ∑_k |⟨e,u_k⟩|² = e*Me = n. A vector u_k on the line ℂe contributes |⟨e,u_k⟩|² = 1, so each line through 0 contains at most n of the 2n vectors. Suppose all F_σ(u) = 0, and let Bᵢⱼ = [vᵢ,wⱼ]. Then the bipartite graph {(i,j) : Bᵢⱼ ≠ 0} has no perfect matching. By Hall's theorem (CITED) some set R of rows has fewer than |R| neighbouring columns. Let C be the set of columns outside the neighbourhood: |R| + |C| ≥ n + 1, and Bᵢⱼ = 0 on R × C. Since R and C are nonempty and all vectors are nonzero, [vᵢ,wⱼ] = 0 means vᵢ ∥ wⱼ, so the |R| + |C| ≥ n + 1 vectors vᵢ (i ∈ R), wⱼ (j ∈ C) lie on one line. This is a contradiction.
3. *Φ(N) is closed.* N = {0} ∪ {r u : r > 0, u ∈ N₁} and Φ(ru) = r^{2n} Φ(u), so Φ(N) = {0} ∪ {ρΦ(u) : ρ > 0, u ∈ N₁}. N₁ is compact and Φ is continuous and nonzero there, so c₁ ≤ |Φ| ≤ c₂ on N₁ with c₁ > 0. Let ρ_m Φ(u_m) → y ≠ 0. Then ρ_m lies eventually in [|y|/(2c₂), 2|y|/c₁]. A subsequence has ρ_m → ρ > 0 and u_m → u ∈ N₁, so y = ρΦ(u).

Hence Φ(Y) is Euclidean-closed. Its Zariski closure X̄ is irreducible, as the closure of the image of an irreducible variety. Φ(Y) is constructible (Chevalley; CITED, e.g. Hartshorne, *Algebraic Geometry*, Ex. II.3.19) and Zariski-dense in X̄, so it contains a nonempty Zariski-open subset O of X̄. O is dense in X̄ in the classical topology (CITED, Mumford, *The Red Book of Varieties and Schemes*, Ch. I §10, Thm. 1). Therefore X̄ ⊆ Φ(Y), and X = Φ(Y) = X̄. ∎

> **Lemma (dimension).** dim X ≥ 2n − 2.

*Proof (PROVED, citing dim = transcendence degree, Hartshorne I.1.7–1.8A).* Consider the slice v₁ = (0, s), v₂ = (1, 0), w₁ = (1, 1), vᵢ = (pᵢ, 1) for i ≥ 3 and wⱼ = (qⱼ, 1) for j ≥ 2. Its 2n − 2 parameters are s, p₃,…,pₙ, q₂,…,qₙ. Write ψ(s,p,q) = Φ(v) ∈ X. The brackets are [v₁,w₁] = −s, [v₁,wⱼ] = −s qⱼ, [v₂,wⱼ] = 1, [vᵢ,w₁] = pᵢ − 1 and [vᵢ,wⱼ] = pᵢ − qⱼ.
- For j ≥ 2 take σ with σ(1) = 1, σ(2) = j, and τ = σ∘(1 2). Then F_σ/F_τ = [v₁,w₁][v₂,wⱼ] / ([v₁,wⱼ][v₂,w₁]) = 1/qⱼ.
- For i ≥ 3 take σ with σ(1) = 1, σ(i) = 2, and τ = σ∘(1 i). Then F_σ/F_τ = (pᵢ − q₂) / (q₂(pᵢ − 1)) =: ρᵢ, so pᵢ = q₂(ρᵢ − 1)/(ρᵢq₂ − 1).
- F_id = −s ∏_{i≥3}(pᵢ − qᵢ).

So there are rational functions R₁,…,R_{2n−2} in the coordinates y_σ with R(ψ(t)) = t for t in a nonempty Zariski-open set of parameters. Their denominators do not vanish identically on X, so they lie in the function field ℂ(X). A polynomial relation F(R₁,…,R_{2n−2}) = 0 in ℂ(X) would give F(t) = 0 on an open set of t, hence F = 0. So tr.deg ℂ(X) ≥ 2n − 2. ∎

(Equality holds: 2n points on ℙ¹ modulo PGL₂, plus the scale. It is not needed.)

## 3. Zheng's Lemma 3.4 in every dimension

Let A : ℂ^{Sₙ} → ℂ^{n×n}, A(y) = ∑_σ y_σ P_σ.

> **Lemma B (κ-line lemma).** If x ∈ ℂ^{Sₙ} and the whole line x + ℂ𝟙 lies in X, then x_σ = ∑ᵢ Q_{i,σ(i)} for some Q ∈ ℂ^{n×n}, and S_μ(x) holds for every admissible μ.

*Proof (PROVED, citing Birkhoff–von Neumann).*
1. *Binomials vanish on X.* Let c ∈ ℤ^{Sₙ} with A(c) = 0, and write c = c⁺ − c⁻ with c^± ≥ 0. Then A(c⁺) = A(c⁻) =: N, and for every v ∈ Y,
   ∏_σ F_σ(v)^{c⁺_σ} = ∏_{i,j} [vᵢ,wⱼ]^{Nᵢⱼ} = ∏_σ F_σ(v)^{c⁻_σ}.
   So f_c(y) = y^{c⁺} − y^{c⁻} vanishes on X.
2. *The direction of the line is a tangent vector at 𝟙.* For t ≠ 0, x + t𝟙 ∈ X, so 𝟙 + t⁻¹x ∈ X by (X1). Hence the polynomial s ↦ f_c(𝟙 + s x) vanishes for all s ≠ 0, so it vanishes identically. Its derivative at s = 0 is ∑_σ c⁺_σ x_σ − ∑_σ c⁻_σ x_σ = ∑_σ c_σ x_σ = 0.
3. *x lies in the row space.* A has integer entries, so integer vectors span ker A over ℂ. By step 2, x is orthogonal (bilinear pairing) to ker A, so x ∈ (ker A)^⊥ = im Aᵀ. That is, x_σ = ⟨Q, P_σ⟩ = ∑ᵢ Q_{i,σ(i)} for some Q.
4. *Birkhoff.* D_μ(x) = ∑_σ μ_σ ∑ᵢ Q_{i,σ(i)} = ∑_{i,j} Qᵢⱼ M(μ)ᵢⱼ. M(μ) is doubly stochastic, so M(μ) = ∑_τ β_τ P_τ with β ≥ 0 and ∑β = 1 (CITED: Birkhoff–von Neumann). Hence D_μ(x) = ∑_τ β_τ x_τ ∈ conv{x_τ}, and Re D_μ(x) ≥ min_τ Re x_τ. ∎

The hypothesis is not vacuous. If a₁ = ⋯ = a_{n−1} = c and Q = ∏ⱼ(c + bⱼ) ≠ 0, then z(a,b) + t𝟙 = α′ z(a′, b) with a′ = (c,…,c, a′ₙ), a′ₙ − c = Q(aₙ − c)/(Q + t) and α′ = (Q + t)/Q, for t ≠ −Q. By closedness the whole line lies in X. (Checked exactly; see "Checks".)

## 4. A maximum principle on a punctured Riemann surface

> **Lemma C.** Let 𝓡 be a compact connected Riemann surface, S ⊂ 𝓡 finite and 𝓡° = 𝓡 ∖ S. Let H₁,…,H_m be holomorphic on 𝓡°, meromorphic on 𝓡, and not all constant, and put h_j = Re H_j. Let 𝓛 = Re G with G holomorphic on 𝓡°. Let Ω be a nonempty connected component of {p ∈ 𝓡° : h_j(p) > 0 for all j, 𝓛(p) < 0} on which 𝓛 is bounded. Then some p′ ∈ 𝓡° has h_j(p′) ≥ 0 for all j, h_{j′}(p′) = 0 for some j′, and 𝓛(p′) < 0.

*Proof (PROVED, using only the local maximum principle for continuous functions with the submean property).* Suppose not. Then:

(H) every p ∈ 𝓡° with all h_j(p) ≥ 0 and some h_j(p) = 0 has 𝓛(p) ≥ 0.

1. *𝓛 = 0 on ∂Ω ∩ 𝓡°.* Let q ∈ 𝓡° be a boundary point of Ω. By continuity all h_j(q) ≥ 0 and 𝓛(q) ≤ 0. If all h_j(q) > 0 and 𝓛(q) < 0, a connected neighbourhood of q inside the open set {all h_j > 0, 𝓛 < 0} meets Ω, so q ∈ Ω, which is impossible. So some h_j(q) = 0, and then 𝓛(q) ≥ 0 by (H); or 𝓛(q) = 0. Either way 𝓛(q) = 0.
2. *A continuous subharmonic barrier.* Put v = −𝓛 on Ω and v = 0 on 𝓡° ∖ Ω. By step 1, v is continuous on 𝓡°, and v ≥ 0. It has the local submean property on small coordinate circles: on Ω it is harmonic, and off Ω it is 0 ≤ any average of v. Since 𝓛 is bounded on Ω, v is bounded on 𝓡°.
3. *v attains its supremum M in 𝓡°.* Let s ∈ S and let z be a coordinate on a disc D centred at s. For ε > 0, w_ε = v + ε log|z| has the submean property on D ∖ {s} and tends to −∞ at s. By the maximum principle on the annuli δ ≤ |z| ≤ r, v(z) + ε log|z| ≤ max_{|ζ|=r} v + ε log r for 0 < |z| ≤ r. Letting ε → 0 gives sup_{0<|z|≤r} v ≤ max_{|ζ|=r} v. Do this at every puncture. The remaining part of 𝓡° (outside the small discs) is compact, so M is attained at a point of 𝓡°.
4. *v is constant.* The set {v = M} ⊆ 𝓡° is closed. It is also open, because at a point of maximum the submean property forces v ≡ M on small circles. The surface 𝓡° is connected, since a connected surface minus finitely many points stays connected. So v ≡ M. As Ω ≠ ∅ and v > 0 on Ω, M > 0, so v > 0 everywhere. Hence Ω = 𝓡°, and all h_j > 0 on 𝓡°.
5. *Contradiction.* If H_j had a pole at s ∈ S, then H_j = c z^{−k}(1 + O(z)) with c ≠ 0 and k ≥ 1, and Re H_j would be negative somewhere on every small punctured disc around s. So each H_j is holomorphic on the compact connected surface 𝓡, hence constant. This contradicts the hypothesis. ∎

This is Zheng's Lemma 3.5 with the constant-sum hypothesis replaced by boundedness of 𝓛 on Ω, which is what the proof uses. Steps 3–4 replace the citation of the removable-singularity theorem for subharmonic functions by its short classical proof.

## 5. Proof of Theorem A

**Step 1 (from spectra to X, and back).** Suppose D_μ(a,b) ∉ K = conv{z_σ(a,b)}. K is compact and convex in ℂ ≅ ℝ², so there is θ with Re(θ D_μ(a,b)) < min_{w∈K} Re(θw) = min_σ Re(θ z_σ(a,b)). Put x₀ = θ z(a,b) = θΦ(ι(a,b)) ∈ X (X1). By linearity D_μ(x₀) = θ D_μ(a,b), so S_μ(x₀) fails. Fix π ∈ T(x₀). Then s := ∑_σ Re g_σ(x₀) > 0, because s = 0 would give L(x₀) = 0. Replace x₀ by x₀/s ∈ X. This changes neither T nor the sign of L, and now ∑_σ Re g_σ(x₀) = 1.

Conversely, let x′ ∈ X with S_μ(x′) false and |T(x′)| ≥ 2n − 2. By (X5) and (X4), x′ = λ z(a′,b′) with λ ≠ 0. With θ = λ/|λ|, Re(θ z_σ(a′,b′)) = Re x′_σ / |λ|. So the minimum m is attained exactly on T(x′), and Re(θ D_μ(a′,b′)) = Re D_μ(x′)/|λ| < m. Hence part (a) follows from part (b).

**Step 2 (one deformation step).**

> **Lemma E.** Let x ∈ X with L(x) < 0, π ∈ T(x) and ∑_σ Re g_σ(x) = 1. Suppose some irreducible component Z of X ∩ Λ(x) through x has dim Z ≥ 1. Then some x′ ∈ X has T(x′) ⊋ T(x), π ∈ T(x′), ∑_σ Re g_σ(x′) = 1 and L(x′) ≤ L(x). If L is not constant on the curve C chosen below, then L(x′) < L(x).

*Proof.* Write T = T(x) and 𝒩 = Sₙ ∖ T.

*Case 1: every g_σ is constant on Z.* Then y_σ − y_π = x_σ − x_π for y ∈ Z and all σ, so y − x ∈ ℂ𝟙. Thus Z is a closed irreducible subset of dimension ≥ 1 of the line x + ℂ𝟙, so Z = x + ℂ𝟙 ⊆ X. Lemma B gives S_μ(x), contradicting L(x) < 0. So this case does not occur.

*Case 2: some g_τ is not constant on Z.* Necessarily τ ∈ 𝒩, because on Λ(x) the g_σ with σ ∈ T are constant.
- Pick y₁ ∈ Z with g_τ(y₁) ≠ g_τ(x), and an irreducible curve C ⊆ Z through x and y₁ (CITED: Mumford, *Abelian Varieties*, §6, Lemma).
- Let ν : C̃ → C be the normalization. C̃ is a smooth irreducible affine curve, ν is finite and surjective, and C̃ = 𝓡° = 𝓡 ∖ S with 𝓡 a compact connected Riemann surface and S finite. Regular functions on C pull back to holomorphic functions on 𝓡° that are meromorphic on 𝓡 (CITED, e.g. Shafarevich, *Basic Algebraic Geometry*, Vol. 1, Ch. II §5 and Ch. VII).
- Put H_σ = g_σ∘ν for σ ∈ 𝒩. Then h_σ = Re H_σ, and H_τ is not constant. On C, g_σ ≡ g_σ(x) for σ ∈ T, with real part 0, and ∑_σ g_σ ≡ ∑_σ g_σ(x). Hence ∑_{σ∈𝒩} h_σ ≡ 1 on 𝓡°, and 𝓛 := L∘ν = ∑_{σ∈𝒩} μ_σ h_σ = Re(∑_{σ∈𝒩} μ_σ H_σ).
- Pick p₀ ∈ ν⁻¹(x). Then h_σ(p₀) = Re g_σ(x) > 0 for σ ∈ 𝒩, and 𝓛(p₀) = L(x) < 0.
- On any set where all h_σ > 0, the identity ∑h_σ = 1 gives 0 < h_σ < 1, so |𝓛| ≤ max_σ |μ_σ|. So 𝓛 is bounded on every set used below.
- If 𝓛 is constant, apply Lemma C to the component of {all h_σ > 0, 𝓛 < 0} containing p₀. This gives p′ with all h_σ(p′) ≥ 0, some h_{σ′}(p′) = 0, and 𝓛(p′) = L(x).
- If 𝓛 is not constant, it is not constant on any nonempty open set (harmonic functions on a connected surface). By the strong maximum principle, p₀ is not a local minimum. So there is p₁ near p₀, with all h_σ(p₁) > 0, and c with 𝓛(p₁) < c < 𝓛(p₀). Apply Lemma C to the same H_σ, to G = ∑μ_σH_σ − c and to the component of {all h_σ > 0, 𝓛 − c < 0} containing p₁. This gives p′ with all h_σ(p′) ≥ 0, some h_{σ′}(p′) = 0, and 𝓛(p′) < c < L(x).

In both cases put x′ = ν(p′) ∈ C ⊆ X. Then g_σ(x′) = g_σ(x) for σ ∈ T, and Re g_σ(x′) = h_σ(p′) ≥ 0 for σ ∈ 𝒩, with equality for σ′. So Re x′_σ ≥ Re x′_π for all σ, π ∈ T(x′), and T(x′) = T ∪ {σ ∈ 𝒩 : h_σ(p′) = 0} ⊋ T. Also ∑_σ Re g_σ(x′) = 1 and L(x′) = 𝓛(p′) ≤ L(x). ∎

**Step 3 (iteration and the count 2n − 2).** Start from x₀ of Step 1. While some irreducible component of X ∩ Λ(x_k) through x_k has positive dimension, apply Lemma E to get x_{k+1}. The contact sets strictly increase inside Sₙ, so the process stops after fewer than n! steps, at some x′. All the listed properties are preserved: π ∈ T, normalization 1, L ≤ L(x₀) < 0, T(x′) ⊇ T(x₀). At x′ every component of X ∩ Λ(x′) through x′ is a point, so x′ is isolated in X ∩ Λ(x′). Λ(x′) is cut out by |T(x′)| affine-linear equations. By the affine dimension theorem (CITED, Hartshorne I.7.1), every component of X ∩ Λ(x′) has dimension ≥ dim X − |T(x′)| ≥ 2n − 2 − |T(x′)|, and one of them is {x′}. Hence |T(x′)| ≥ 2n − 2. This proves (b), and with Step 1 also (a). ∎

## 6. What this gives for n = 5

- Any violation of the conjecture at n = 5, for any U ∈ U(5), can be moved, with the same U, to spectra at which at least 8 products lie on the separating supporting line. The configuration is then rigid: it is isolated in its level set, and the normalized support gap is no larger. No orthogonal reduction and no knowledge of the ideal of X₅ is used.
- What remains is a statement about rigid contact configurations. For a contact set T, the gap vector lies in the cone N_T = {p in the 42-dimensional span of the squared-minor vectors : p ≥ 0, p_T = 0}. If every generator of N_T is a squared-minor vector, the support inequality holds for every unitary, by Lemma 1(a) (the "automatic" strata). The numbers of rigid patterns and of automatic ones at n = 5 are not known. They are being sampled in `mdo-n5-contact-strata`, and the same strata are being searched for violations in `mdo-n5-high-contact-search`.

## Citations

| Fact | Where it is used | Reference |
| --- | --- | --- |
| Laplace expansion; Jacobi's complementary-minor identity | Lemma 1 | Horn and Johnson, *Matrix Analysis*, 2nd ed., §0.8 |
| Hall's marriage theorem | Lemma D, step 2 | any combinatorics text |
| Chevalley: images of morphisms are constructible | Lemma D | Hartshorne, *Algebraic Geometry*, Ex. II.3.19 |
| Zariski-open dense subsets are classically dense | Lemma D | Mumford, *The Red Book*, Ch. I §10, Thm. 1 |
| dim = transcendence degree of the function field | dimension lemma | Hartshorne, Prop. I.1.7 and Thm. I.1.8A |
| Birkhoff–von Neumann | Lemma B | e.g. Bhatia, *Matrix Analysis*, Thm. II.2.3 |
| Affine dimension theorem | Step 3 | Hartshorne, Prop. I.7.1 |
| Irreducible curve through two points | Lemma E | Mumford, *Abelian Varieties*, §6, Lemma |
| Normalization of curves; smooth curves as punctured compact Riemann surfaces; regular functions are meromorphic there | Lemma E | Shafarevich, *Basic Algebraic Geometry* 1, Ch. II §5; Ch. VII |
| Local maximum principle for continuous functions with the submean property; identity theorem and strong maximum principle for harmonic functions | Lemma C, Lemma E | e.g. Ransford, *Potential Theory in the Complex Plane*, Ch. 1–2 |

## Checks

See [`check_identities.py`](check_identities.py) and [`deformation_test.py`](deformation_test.py), with their outputs. Reproduce from the repository root with `python -I -B research/marcus-de-oliveira/n5-contact-program/contact-reduction/<script>.py`.
