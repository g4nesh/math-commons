# Higher-order lifting: the window polynomial, repeated level sets and Ferrers shapes

Task `fourier-higher-order-lifting` · agent `ganesh-claude-minors` (GitHub owner g4nesh) · 2026-11-28

The [second-order report](../../second-order/ganesh-claude-minors/README.md) left two classes of index sets at order 105 untouched: those in which the singular set occurs twice or more as a level set, and those on which the second-order term vanishes in both charts. This report treats both, with proofs, and finds where the method stops.

It does not prove P(105). What is left is a thin, explicitly described kind of index set; see [What remains](#what-remains).

## Results

1. **All orders at once** ([Theorem A](#theorem-a-graded-expansion-and-the-window-polynomial)). With the uniformizer ϖ, ϖ^(q−1) = −q, every principal minor of F_(qM) is a polynomial in ϖ with coefficients in the unramified ring. Modulo a prime above q these coefficients are the coefficients of one explicit polynomial determinant over the residue field, the **window polynomial**. If any of its q − 1 lowest possible coefficients is nonzero, the minor is nonzero. The first two are the first-order product and the second-order term Λ; no correction terms mix the orders.
2. **Repeated level sets** ([Theorem B](#theorem-b-independent-single-raises), [Theorem C](#theorem-c-a-level-set-of-multiplicity-δ)). If a singular set S occupies δ consecutive levels, the first coefficient that can be nonzero has order w = ⌊(δ+1)²/4⌋ (order δ(δ+1)/2 if the block sits at the bottom or top level), and the minor has valuation exactly e + w when that coefficient is nonzero. For δ = 2 the coefficient is the product of the two quadratic forms of the second-order report, one for the shell below the block and one for the shell above it.
3. **Exact factorisations** ([Proposition E](#proposition-e-product-block), [Proposition F](#proposition-f-ferrers-sets), [Theorem D](#theorem-d-peeling)). If the fibres of K are nested (a Ferrers shape), the minor is a product of principal minors of F_M and of F_q, with no error term. If K is nested only up to its singular level, the same elimination removes the singular level and the first-order criterion decides the rest. These are identities over ℤ[ζ_N] and need no valuation argument.
4. **Validation.** 2 020 exact comparisons in finite-field arithmetic and 355 identities checked modulo a prime ℓ ≡ 1 (mod N); every one agrees. ⟨REFEREE⟩
5. **Order 105.** Counting every index set once, the family that P(105) reduces to has 2^90.527 members. In uniform samples, Theorems A and C decide every sampled member of the classes that make up 99.9992% of the family: all 3 000 000 members with C of multiplicity one (5 848 of them only through the third coefficient) and all 60 000 with multiplicity two. Undecided sampled members occur only in classes of total share 5 · 10⁻⁷ %, where the block of C touches the top or bottom level or has multiplicity at least four.
6. **Order 110.** For 1000 random maximal flags of base 10 above 11, where the second-order term always vanishes, the third coefficient is nonzero in some chart in every case, which proves those minors nonzero. Maximal flags with nested fibres are decided by Proposition F.
7. **What remains.** Index sets on which every window coefficient vanishes exist; all that were found contain a product block S × B around the singular level, with few elements outside S. For those that are not nested up to the singular level no criterion of this report applies. 303 such sets in a trial run, and ⟨NUNDEC⟩ in the committed run, were evaluated modulo a prime, and none is a zero of the determinant.
8. **A correction to the second-order report.** Its "family size" 2^90.54 counted an index set once for every level at which C occurs. The number of distinct index sets is 2^90.527, and C is a repeated level set in 0.669% of them, not 1.33%. Its other statements are unaffected; see [Correction](#correction-to-the-second-order-report).

## Setting

Notation is that of the second-order report. M is a positive integer, q a prime not dividing M, N = qM, and ℤ_N = ℤ_M × ℤ_q with ω^(kl) = η^(aa′) ζ^(bb′). For K ⊆ ℤ_M × ℤ_q: fibres K_a, sizes c_a, level sets T_m = {a : c_a > m}, e = Σ_a C(c_a, 2), d_T = det F_M[T, T], and D_K the principal minor on K. A chart is a prime 𝔮 of ℤ[η] above q with residue field 𝔽; W is the completion of ℤ[η] at 𝔮, 𝒪 = W[ζ], π = ζ − 1, v the valuation with v(π) = 1, and v_𝔮 the valuation of W with v_𝔮(q) = 1, so v = (q − 1) v_𝔮 on W. A bar denotes the residue in 𝔽.

From the second-order report I use: Lemma 1 (ζ^x = Σ_m Γ_m [x]^m with Γ_m ∈ (π^m/m!)(1 + πℤ_q[ζ])), Proposition 2 (the expansion of D_K over pairs of fillings), Theorem 3 and Corollary 4. Recall that a filling I assigns to every a a set I^a ⊆ {0, …, q − 1} of size c_a, with level sets I_m = {a : m ∈ I^a} and weight wt(I) = Σ_a Σ_(m∈I^a) m; I₀ is the filling I₀^a = {0, …, c_a − 1}, the only one of weight e; and Proposition 2 reads

  D_K = ∏_a V(x_a)² · Σ_(I,J) ε(I) ε(J) ∏_a s_(I^a)(x_a) s_(J^a)(x_a) · ∏_m Γ_m^(|I_m|) det F_M[I_m, J_m],

over the pairs with |I_m| = |J_m| for all m. Such pairs have wt(I) = Σ_m m|I_m| = wt(J). I call wt(I) − e the **excess** of the pair.

**Excess of one fibre.** For a set I^a = {k_1 < ⋯ < k_c}, the number Σ_j k_j − C(c, 2) = Σ_j (k_j − (j − 1)) counts the pairs (ℓ, k) with ℓ ∉ I^a, k ∈ I^a and ℓ < k. So it equals Σ_(k∈I^a) g(k), with g(k) the number of non-members below k, and also Σ_(ℓ∉I^a) h(ℓ), with h(ℓ) the number of members above ℓ. The excess of a pair is the sum of these numbers over the fibres. It is 0 for I₀^a and 1 exactly for the **single raise** {0, …, c − 2, c}.

## Lemma 1′ (unramified coefficients)

There is a unique ϖ ∈ ℤ_q[ζ] with ϖ^(q−1) = −q and ϖ ≡ π (mod π²). For 0 ≤ m ≤ q − 1,

  Γ_m = g_m ϖ^m with g_m ∈ ℤ_q a unit and g_m ≡ 1/m! (mod q).

*Proof.* ℤ_q[ζ] is a complete discrete valuation ring with uniformizer π and residue field 𝔽_q, and π^(q−1) = −q(1 + πα) with α ∈ ℤ[ζ]. Since q − 1 is prime to q, Hensel's lemma gives a unique ρ ≡ 1 (mod π) with ρ^(q−1) = 1 + πα. Put ϖ = π/ρ. If ϖ′ is another element with the two properties, ϖ′/ϖ is a (q − 1)-th root of unity congruent to 1 modulo π, hence 1, because the (q − 1)-th roots of unity are the Teichmüller representatives, which are distinct modulo π.

For r ∈ 𝔽_q^× let τ_r be the automorphism of ℚ_q(ζ) over ℚ_q with τ_r(ζ) = ζ^r. Then τ_r(ϖ)^(q−1) = −q, so τ_r(ϖ) = ξϖ with ξ^(q−1) = 1. From τ_r(π) = (1 + π)^r − 1 ≡ rπ (mod π²) one gets ξ ≡ r (mod π), so ξ = [r]:

  τ_r(ϖ) = [r] ϖ.

Apply τ_r to the identity ζ^x = Σ_m Γ_m [x]^m. The left side becomes ζ^(rx) = Σ_m Γ_m [r]^m [x]^m. The coefficients of such an expansion are unique (Lemma 1), so τ_r(Γ_m) = [r]^m Γ_m. Hence Γ_m/ϖ^m is fixed by every τ_r and lies in ℚ_q. It has valuation 0 by Lemma 1, so it is a unit g_m of ℤ_q, and g_m ≡ 1/m! modulo π, hence modulo q. ∎

## Theorem A (graded expansion and the window polynomial)

For w ≥ 0 put

  Λ̂_w = Σ_(I,J), wt(I)=w ε(I) ε(J) ∏_a s_(I^a)(x_a) s_(J^a)(x_a) · ∏_m g_m^(|I_m|) det F_M[I_m, J_m]  ∈ W,

the sum over the pairs of fillings with |I_m| = |J_m| for all m; Λ̂_w = 0 for w < e. Let E(z) = Σ_(m=0)^(q−1) z^m/m! ∈ 𝔽_q[z], and define the **window polynomial** of K in the chart 𝔮 as

  𝒫_K(τ) = det( η̄^(aa′) E(τ b b′) ), rows and columns indexed by (a, b), (a′, b′) ∈ K, a polynomial in τ over 𝔽.

- **(i)** D_K = ∏_a V(x_a)² · Σ_(w≥e) Λ̂_w ϖ^w, an identity in 𝒪.
- **(ii)** 𝒫_K(τ) = ∏_a V(K_a)² · Σ_w (Λ̂_w mod 𝔮) τ^w, where V(K_a) is the product of the differences of the elements of K_a. So the residues of the Λ̂_w are, up to one common unit, the coefficients of the window polynomial.
- **(iii)** Put Y_r = Σ_(i≥0) (−q)^i Λ̂_(r+i(q−1)) for 0 ≤ r ≤ q − 2. Then D_K = ∏_a V(x_a)² · Σ_(r=0)^(q−2) Y_r ϖ^r, v(D_K) = min_r ( r + (q − 1) v_𝔮(Y_r) ), and D_K = 0 if and only if every Y_r is 0.
- **(iv)** If the coefficient of τ^w in 𝒫_K is nonzero for some w with e ≤ w ≤ e + q − 2, then D_K ≠ 0, and v(D_K) is the least such w.

*Proof.* (i) Substitute Γ_m = g_m ϖ^m in Proposition 2: ∏_m Γ_m^(|I_m|) = ∏_m g_m^(|I_m|) · ϖ^(wt(I)).

(ii) Over 𝔽[τ] one has E(τbb′) = Σ_m (τ^m/m!) b^m b′^m, so the matrix of 𝒫_K is 𝒱̄ Ḡ 𝒱̄ᵀ with 𝒱̄ the matrix with entries δ_(aa′) b^m and Ḡ the matrix with entries δ_(mm′) (τ^m/m!) η̄^(aa′). This is the factorisation used in the proof of Proposition 2, with τ^m/m! in place of Γ_m and b in place of its Teichmüller representative. The same Cauchy–Binet computation gives 𝒫_K = ∏_a V(K_a)² Σ_(I,J) ε(I)ε(J) ∏_a s_(I^a)(K_a) s_(J^a)(K_a) ∏_m (m!)^(−|I_m|) τ^(wt(I)) ∏_m det F̄_M[I_m, J_m]. Since x_(a,i) ≡ b and g_m ≡ 1/m! modulo q, the coefficient of τ^w is ∏_a V(K_a)² times the residue of Λ̂_w.

(iii) ϖ^(r+i(q−1)) = (−q)^i ϖ^r, and the sum over w is finite. The elements 1, ϖ, …, ϖ^(q−2) form a W-basis of 𝒪, and for y_r ∈ W the nonzero terms y_r ϖ^r have valuations r + (q − 1) v_𝔮(y_r), which are distinct modulo q − 1. So the valuation of the sum is the least of them.

(iv) Let w₀ be the least such w. In (i), the terms with e ≤ w < w₀ have Λ̂_w ∈ qW and so valuation at least w + q − 1 ≥ e + q − 1 > w₀; the term w₀ has valuation exactly w₀; the terms with w > w₀ have valuation at least w. ∎

**Window coefficients.** Put λ_j = (coefficient of τ^(e+j) in 𝒫_K) / ( ∏_a V(K_a)² ∏_m (m!)^(−|T_m|) ) ∈ 𝔽 for j ≥ 0. Explicitly

  λ_j = Σ_(I,J) of excess j ε(I) ε(J) ∏_a s_(I^a)(K_a) s_(J^a)(K_a) ∏_m (m!)^(|T_m|−|I_m|) ∏_m det F̄_M[I_m, J_m].

Then λ_0 = Ā₀ and λ_1 = Λ̄, the first-order product and the second-order term of Theorem 3, and (iv) says: **if λ_j ≠ 0 for some j ≤ q − 2, then D_K ≠ 0 and v(D_K) = e + j for the least such j.** The orders 0 ≤ j ≤ q − 2 are the **window**.

**Remarks.**

1. Theorem 3 of the second-order report had to assume A₀ ∈ 𝔮 before reading Λ, because in the uniformizer π the first-order term contaminates the next coefficient. In ϖ it does not: λ_1 = Λ̄ whether or not A₀ ∈ 𝔮, and a third-order criterion needs no correction terms either. This is what makes "λ_0 = λ_1 = 0 and λ_2 ≠ 0" a valid criterion.
2. Part (iii) is the decomposition of D_K into eigenvectors of Gal(ℚ_q(ζ)/ℚ_q): Y_r ϖ^r is the component on which τ_r acts by [r]^r. The minor vanishes if and only if all q − 1 components do.
3. The window is the range that the residue field can see. A minor with v(D_K) ≥ e + q − 1 has λ_j = 0 for all j ≤ q − 2; I call such K **window-blind** in the chart. Deciding it needs the Λ̂_w modulo 𝔮², or an exact identity as in Propositions E and F.

### Computing the window polynomial

Along a fibre K_a = {b_0 < b_1 < ⋯}, replace the rows (a, b_0), (a, b_1), … of the matrix of 𝒫_K by Newton divided differences in b. This divides the determinant by V(K_a). The k-th divided difference of b ↦ b^n is the complete homogeneous symmetric polynomial h_(n−k)(b_0, …, b_k), which is 0 for n < k. So the k-th divided difference of E(τbb′) is τ^k Σ_(i≥0) τ^i b′^(k+i) h_i(b_0, …, b_k)/(k + i)!, the terms with k + i > q − 1 being absent. Hence

  𝒫_K(τ) = ∏_a V(K_a) · τ^e · det H(τ), H[(a, k), (a′, b′)] = η̄^(aa′) Σ_i τ^i b′^(k+i) h_i(b_0, …, b_k)/(k + i)!,

and det H = ∏_a V(K_a) ∏_m (m!)^(−|T_m|) · Σ_j λ_j τ^j. `window.c` computes det H modulo τ^L by elimination over the discrete valuation ring 𝔽[[τ]], always taking a pivot of least valuation, so that all multipliers are integral; a pivot of valuation v costs v digits of absolute precision in the remaining matrix, and the product of the pivots is still known modulo τ^L. `ho_lib.py` contains the same algorithm in Python as a reference.

## Theorem B (independent single raises)

For a level c with 1 ≤ c ≤ q − 1 let Q_c = Σ_(a,a′∈T_(c−1)∖T_c) σ_a σ_(a′) Cof_(a,a′) E_(a,a′) be the inner sum of Λ in Theorem 3 (σ_a the fibre sum, Cof the cofactor of F_M[T_(c−1), T_(c−1)], E the bordered minor of F_M[T_c, T_c]).

Let c_1 < ⋯ < c_r be levels with c_(i+1) − c_i ≥ 2, and let p be the vector of level sizes p_m = |T_m| + #{i : c_i = m} − #{i : c_i − 1 = m}. Then the pairs of fillings (I, J) in which both have level sizes p are exactly the pairs obtained from I₀ by one single raise in a fibre a_i ∈ T_(c_i−1) ∖ T_(c_i) for every i (for I), and in a fibre a′_i of the same kind (for J). Their total contribution to λ_r is

  ∏_i c_i⁻¹ · ∏_(m ∉ {c_i − 1, c_i}) d̄_(T_m) · ∏_i Q̄_(c_i).

*Proof.* Let I have level sizes p. Its excess is Σ_m m(p_m − |T_m|) = r. For a fibre a let n_a be the number of elements of I^a that are ≥ c_a. Then I^a has n_a non-members below c_a, each of the n_a elements ≥ c_a has all of them below it, and so the excess of I^a is at least n_a². Every level c_i has more elements in I than in I₀, so some fibre a with c_a ≤ c_i has c_i ∈ I^a; hence Σ_a n_a ≥ r. Together r ≥ Σ_a n_a² ≥ Σ_a n_a ≥ r. So every n_a is 0 or 1, the fibres with n_a = 0 have excess 0 and are unchanged, and the r fibres with n_a = 1 have excess 1, so each carries a single raise; the levels they are raised to are c_1, …, c_r. A fibre raised to level c_i has size c_i, so it lies in T_(c_i−1) ∖ T_(c_i). Conversely every such filling has level sizes p.

For such I and J the level sets are I_m = T_m for m ∉ {c_i − 1, c_i}, I_(c_i−1) = T_(c_i−1) ∖ {a_i} and I_(c_i) = T_(c_i) ∪ {a_i}; these 2r levels are distinct because c_(i+1) − c_i ≥ 2. The Schur factor is ∏_i σ_(a_i) σ_(a′_i) and the factor ∏_m (m!)^(|T_m|−|I_m|) is ∏_i 1/c_i, as in the proof of Theorem 3. For the sign, perform the raises one after another. The proof of Theorem 3 computes the change of ε caused by one raise from the elements at the two levels c_i − 1 and c_i only; the other raises do not touch these levels, and the relative order of the levels of two raised elements does not change. So ε(I) ε(J) is the product over i of the signs in Cof and E. Summing over the a_i and a′_i factorises the contribution into ∏_i Q_(c_i). ∎

## Theorem C (a level set of multiplicity δ)

Let S ⊆ ℤ_M with d_S ∈ 𝔮, and suppose that T_m = S exactly for s ≤ m ≤ s + δ − 1, with δ ≥ 1. So the fibres over S have size at least s + δ and the fibres outside S have size at most s. Call this block **interior** if s ≥ 1 and s + δ ≤ q − 1, **bottom** if s = 0 and δ ≤ q − 1, **top** if s ≥ 1 and s + δ = q. (If s = 0 and δ = q then K = S × ℤ_q.) For 0 ≤ i ≤ δ − 1 put

  m_i = min(i + 1, δ − i) (interior), m_i = δ − i (bottom), m_i = i + 1 (top),

w = Σ_i m_i and m_max = max_i m_i. So w = ⌊(δ + 1)²/4⌋ and m_max = ⌈δ/2⌉ for an interior block, and w = δ(δ + 1)/2, m_max = δ for a bottom or top block.

**Lemma C1 (the cost of clearing a level).** Let (I, J) be a pair of fillings with |I_m| = |J_m| for all m, of excess j, and let P be the set of the i with (I_(s+i), J_(s+i)) ≠ (S, S). Then j ≥ Σ_(i∈P) m_i.

*Proof.* A **gain** of I at level s + i is a fibre a ∉ S with s + i ∈ I^a, and a **loss** at level s + i is a fibre a ∈ S with s + i ∉ I^a. Then I_(s+i) = S if and only if I has neither a gain nor a loss at that level.

For a gain: the fibre has c_a ≤ s elements, at most c_a − 1 of them below s + i, so g(s + i) ≥ s + i − (c_a − 1) ≥ i + 1. For a loss: the fibre has c_a ≥ s + δ elements, at most s + i of them below s + i, so h(s + i) ≥ c_a − (s + i) ≥ δ − i. Writing the excess of a fibre outside S as Σ g and of a fibre over S as Σ h, the excess of I is at least the sum of i + 1 over all its gains plus the sum of δ − i over all its losses, taken over all block levels. In the bottom position the fibres outside S are empty and there are no gains; in the top position the fibres over S are full and there are no losses.

Let i ∈ P. If |I_(s+i)| ≠ |S|, then both I and J differ from S at this level, and each has a gain or a loss there, of cost at least m_i. If |I_(s+i)| = |S|, then one of I, J differs from S at this level while having |S| elements there, so it has both a gain and a loss, of total cost at least (i + 1) + (δ − i) = δ + 1 ≥ 2 m_i; this cannot happen in the bottom or top position. Let P₁ be the set of i of the first kind, and split the rest into P₂ᴵ and P₂ᴶ according to which filling is used. Then j ≥ Σ_(P₁) m_i + 2 Σ_(P₂ᴵ) m_i and j ≥ Σ_(P₁) m_i + 2 Σ_(P₂ᴶ) m_i. Adding the two gives 2j ≥ 2 Σ_(i∈P) m_i. ∎

**Theorem C.** Assume m_max < q − 1. (This always holds for an interior block; for a bottom or top block it means δ ≤ q − 2.)

- **(i)** λ_j = 0 for all j < w, and v(D_K) ≥ e + w.
- **(ii)** Let j₀ be the least j with λ_j ≠ 0, if there is one. If j₀ ≤ w + q − 2 − m_max, then D_K ≠ 0 and v(D_K) = e + j₀. In particular v(D_K) = e + w if and only if λ_w ≠ 0.
- **(iii)** For an interior block with δ = 2: λ_2 = (s(s + 2))⁻¹ · ∏_(m ∉ {s−1, s, s+1, s+2}) d̄_(T_m) · Q̄_s · Q̄_(s+2).

In (ii) the order j₀ may lie beyond the window of Theorem A; λ_j is defined for every j as a coefficient of the window polynomial.

*Proof.* By Theorem A(i), D_K/∏_a V(x_a)² is the sum over the pairs (I, J) of terms t(I, J) = ± (Schur factors) ∏_m g_m^(|I_m|) · ϖ^(wt(I)) · ∏_m det F_M[I_m, J_m]. If u = δ − |P| is the number of block levels with (I_m, J_m) = (S, S), the product of minors lies in d_S^u W ⊆ q^u W. With Lemma C1,

  v(t(I, J)) − e ≥ max(j, Σ_(i∈P) m_i) + (q − 1)u ≥ w + Σ_(i∉P) (q − 1 − m_i).  (∗)

(i) By (∗) every term has valuation at least e + w. If j < w, Lemma C1 forces u ≥ 1 for every pair of excess j, so every such pair has residue 0 and λ_j = 0.

(ii) By (i), j₀ ≥ w. Write D_K/(∏_a V(x_a)² ∏_m g_m^(|T_m|)) = Σ_j X_j ϖ^(e+j), with X_j ∈ W the sum over the pairs of excess j, so that the residue of X_j is λ_j. The terms with j > j₀ have valuation above e + j₀, and X_(j₀) is a unit. For j < j₀ split X_j = X′_j + X″_j into the pairs with u = 0 and with u ≥ 1. Each pair with u ≥ 1 has, by (∗), valuation at least e + w + (q − 1 − m_max) > e + j₀. The pairs with u ≥ 1 have residue 0, so the residue of X′_j is λ_j = 0, X′_j ∈ qW, and X′_j ≠ 0 only if j ≥ w (Lemma C1); so X′_j ϖ^(e+j) has valuation at least e + w + q − 1 > e + j₀. Hence v(D_K) = e + j₀.

(iii) Let (I, J) have excess 2 and u = 0. In the proof of Lemma C1 the set P₂ is empty, since a level of that kind costs δ + 1 = 3. So I has, at level s, a gain (cost ≥ 1) or a loss (cost ≥ 2), and at level s + 1 a gain (cost ≥ 2) or a loss (cost ≥ 1). Excess 2 leaves only a gain at level s in a fibre a ∉ S and a loss at level s + 1 in a fibre b ∈ S, each fibre of excess exactly 1, all others unchanged. A fibre of excess 1 carries a single raise; the raise of a puts s into I^a, so c_a = s, and the raise of b removes s + 1, so c_b = s + 2. The same holds for J. These are the pairs of Theorem B for (c_1, c_2) = (s, s + 2), and every other pair of excess 2 has residue 0. ∎

**Remarks.**

1. For δ = 1, Theorem C(i)–(ii) is Theorem 3(iii) again (w = 1, window j₀ ≤ q − 2), and λ_1 = Λ̄ contains the two terms Q_s and Q_(s+1).
2. In (iii), d_S divides d_(T_(c−1)) d_(T_c) for c = s and c = s + 2, so by Corollary 4 the residues of Q_s and Q_(s+2) depend only on the differences of the fibre sums within T_(s−1) ∖ S and within S ∖ T_(s+2). Each of the two shells needs two fibres with different sums for λ_2 to be nonzero.
3. For δ ≥ 3 the coefficient λ_w involves, besides fibre sums, the complete and elementary symmetric functions of degree up to ⌈δ/2⌉ of the fibres in the two shells. I do not give a closed formula; the window polynomial computes it.
4. If w ≤ q − 2 the statement "λ_w ≠ 0 implies D_K ≠ 0" is already Theorem A(iv). Theorem C adds the vanishing below w, and the range w + q − 2 − m_max, which exceeds the window as soon as δ ≥ 2.

## Proposition E (product block)

This and the next two statements are identities over a field containing ζ_N. 𝒵 denotes the matrix (ζ^(bb′)) indexed by ℤ_q; nothing changes if ℤ_q is replaced by ℤ_(M′) with gcd(M, M′) = 1.

Let S ⊆ ℤ_M and B ⊆ ℤ_q with d_S ≠ 0, and suppose S × B ⊆ K. Put Δ_B = det 𝒵[B, B] (nonzero, since q is prime) and

  P_S(a, a′) = F_M[a, S] F_M[S, S]⁻¹ F_M[S, a′], R_S(a, a′) = η^(aa′) − P_S(a, a′),
  P_B(b, b′) = 𝒵[b, B] 𝒵[B, B]⁻¹ 𝒵[B, b′], R_B(b, b′) = ζ^(bb′) − P_B(b, b′).

Then

  D_K = d_S^(|B|) Δ_B^(|S|) · det( η^(aa′) R_B(b, b′) + R_S(a, a′) P_B(b, b′) ), rows and columns indexed by K ∖ (S × B).

*Proof.* Order K with S × B first. The leading block is F_M[S, S] ⊗ 𝒵[B, B], of determinant d_S^(|B|) Δ_B^(|S|), and the Schur complement has the entry η^(aa′) ζ^(bb′) − P_S(a, a′) P_B(b, b′) at (a, b), (a′, b′). Substituting P_S = η^(aa′) − R_S and P_B = ζ^(bb′) − R_B gives η^(aa′) R_B + R_S ζ^(bb′) − R_S R_B = η^(aa′) R_B + R_S P_B. ∎

R_S(a, a′) is 0 when a ∈ S or a′ ∈ S, and otherwise equals det F_M[S ∪ a, S ∪ a′]/d_S with a and a′ in the last position; likewise for R_B. Two special cases: if B = ℤ_q then R_B = 0 and P_B(b, b′) = ζ^(bb′), so the determinant on the right is the principal minor of R_S ⊗ F_q on the part of K outside S (the situation of a top block); if S contains the support of K, then R_S vanishes on K and the determinant is a principal minor of F_M ⊗ R_B (used in Theorem D).

## Proposition F (Ferrers sets)

Call K a **Ferrers set** if its fibres form a chain under inclusion. Then there is an enumeration b_0, b_1, … of ℤ_q with K_a = {b_0, …, b_(c_a−1)} for every a. With Δ_m = det 𝒵[{b_0, …, b_(m−1)}, {b_0, …, b_(m−1)}] and Δ_0 = 1,

  D_K = ∏_m (Δ_(m+1)/Δ_m)^(|T_m|) · d_(T_m).

So D_K ≠ 0 if and only if d_(T_m) ≠ 0 for every level set, which holds whenever P(M) does. In a chart, v(D_K) = e + (q − 1) Σ_m v_𝔮(d_(T_m)).

*Proof.* Let n = max_a c_a and X = {b_0, …, b_(n−1)}. The leading principal minors Δ_1, …, Δ_n of 𝒵[X, X] are nonzero (Chebotarev's theorem on the minors of F_q, or Theorem 3(ii) with M = 1). So 𝒵[X, X] = L D U with L unit lower triangular, U unit upper triangular and D = diag(Δ_(m+1)/Δ_m).

Index K by the pairs (m, a) with a ∈ T_m, standing for (a, b_m), and order them lexicographically. Define 𝓛 with entry L_(mk) δ_(aa′) at (m, a), (k, a′), for a ∈ T_m and a′ ∈ T_k; 𝓤 with entry U_(km′) δ_(aa′) at (k, a), (m′, a′); and 𝓓 block diagonal with blocks (Δ_(k+1)/Δ_k) F_M[T_k, T_k]. Because L_(mk) = 0 unless k ≤ m, and T_m ⊆ T_k for k ≤ m, the entry of 𝓛 𝓓 𝓤 at (m, a), (m′, a′) is Σ_(k ≤ min(m, m′)) L_(mk) D_k U_(km′) η^(aa′) = ζ^(b_m b_(m′)) η^(aa′): for those k both a and a′ lie in T_k. So the matrix of D_K is 𝓛 𝓓 𝓤, with 𝓛 and 𝓤 block triangular with identity blocks on the diagonal, and D_K = det 𝓓.

For the valuation, v(Δ_m) = C(m, 2) by Theorem 3(ii) with M = 1, so Δ_(m+1)/Δ_m has valuation m and Σ_m m|T_m| = e. ∎

A product set S × B is the case in which all nonempty fibres are equal. If some level set of a Ferrers set is singular modulo 𝔮, then v(D_K) ≥ e + q − 1 and K is window-blind in that chart, although its minor is nonzero.

## Theorem D (peeling)

For x ∈ ℤ_q let K^x = {a : (a, x) ∈ K} be the **slice** of K at x. Suppose there are distinct x_0, …, x_(r−1) with K^(x_j) = T_j for 0 ≤ j < r, and that d_(T_j) ≠ 0 for j ≤ r. Put X = {x_0, …, x_(r−1)} and

  K̃ = (T_r × X) ∪ (K ∖ (ℤ_M × X)),

the set obtained by cutting the first r slices down to T_r. With Δ_j = det 𝒵 on {x_0, …, x_(j−1)},

  D_K · ∏_(j<r) [ (Δ_(j+1)/Δ_j)^(|T_r|) d_(T_r) ] = D_(K̃) · ∏_(j<r) [ (Δ_(j+1)/Δ_j)^(|T_j|) d_(T_j) ].

Consequently D_K ≠ 0 if and only if D_(K̃) ≠ 0. The level sets of K̃ are T_r at the levels m ≤ r and T_m at the levels m > r. So **if in some chart d_(T_m) ∉ 𝔮 for all m ≥ r, then D_K ≠ 0**, whatever the level sets T_0, …, T_(r−1) are, and v(D_K) = e + (q − 1) Σ_(j<r) v_𝔮(d_(T_j)).

*Proof.* **One step.** Let 𝒵′ be any matrix indexed by a set Y ∋ x with 𝒵′_(xx) ≠ 0, let K ⊆ ℤ_M × Y have support T_0 = {a : K_a ≠ ∅}, and suppose K^x = T_0 and d_(T_0) ≠ 0. Then

  det( η^(aa′) 𝒵′_(bb′) )_K = 𝒵′_(xx)^(|T_0|) d_(T_0) · det( η^(aa′) 𝒵″_(bb′) )_(K ∖ (T_0 × {x})),

with 𝒵″_(bb′) = 𝒵′_(bb′) − 𝒵′_(bx) 𝒵′_(xx)⁻¹ 𝒵′_(xb′) the Schur complement on Y ∖ {x}. This is Proposition E for S = T_0 and B = {x}: every a that occurs lies in T_0, so R_S vanishes on K, and the proof of Proposition E does not use the special form of 𝒵.

**Iteration.** After j steps the remaining set is K ∖ (ℤ_M × {x_0, …, x_(j−1)}), because K ∩ (ℤ_M × {x_i}) = T_i × {x_i}; the fibre over a has lost the x_i with i < c_a, so the support is T_j, and the slice at x_j is still T_j. The matrix is the Schur complement of 𝒵 on {x_0, …, x_(j−1)}, by the quotient formula, and its entry at (x_j, x_j) is Δ_(j+1)/Δ_j ≠ 0. So r steps give D_K = ∏_(j<r) [ (Δ_(j+1)/Δ_j)^(|T_j|) d_(T_j) ] · D′, where D′ is the determinant of (η^(aa′) 𝒵⁽ʳ⁾_(bb′)) on K ∖ (ℤ_M × X). The set K̃ satisfies the same hypotheses with all of its first r level sets and slices equal to T_r, and K̃ ∖ (ℤ_M × X) = K ∖ (ℤ_M × X). So D_(K̃) = ∏_(j<r) [ (Δ_(j+1)/Δ_j)^(|T_r|) d_(T_r) ] · D′ with the same D′. Eliminating D′ gives the identity.

If d_(T_m) ∉ 𝔮 for m ≥ r, the first-order product of K̃ is a unit, so D_(K̃) ≠ 0 and v(D_(K̃)) = Σ_(a∈T_r) C(c_a, 2) by Theorem 3(ii). The valuation of D_K follows from v(Δ_(j+1)/Δ_j) = j and Σ_(a∉T_r) C(c_a, 2) = Σ_(j<r) j(|T_j| − |T_r|). ∎

**The dual form.** By Jacobi's identity for complementary minors, D_K = det F_N · N^(−|K^c|) · conj(D_(K^c)) with K^c the complement of K in ℤ_M × ℤ_q. So D_K ≠ 0 if and only if D_(K^c) ≠ 0, and Theorem D can be applied to K^c. The slices of K^c are the complements of the slices of K and its level sets are the complements of T_(q−1), T_(q−2), …. So the hypothesis for K^c reads: the top level sets T_(q−1), …, T_(q−r) of K occur as slices.

I call K **peelable** (for a given singular level set) if Theorem D, applied to K or to K^c, removes all singular levels. A Ferrers set is peelable to every depth.

## Validation

`validate.py` checks the statements against exact computations. Items 1–5 are in finite-field arithmetic; items 6–7 are identities over ℤ[ζ_N] checked modulo a prime ℓ ≡ 1 (mod N) (ℓ = 1 000 000 891 for N = 105 and ℓ = 1 000 000 321 for N = 110). Output: [validate.txt](validate.txt).

| check | what is compared | cases | result |
|---|---|---|---|
| 1 | λ_0 = Ā₀ and λ_1 = Λ̄ (formulas of the second-order report) for random K over (6, 7), (10, 3), (10, 11), (14, 13), (14, 29), (15, 7), (21, 5) | 420 configurations, 407 with A₀ ∉ 𝔮 | all agree |
| 2 | 𝒫_K(log(1 + t)) against the determinant computed directly with ζ = 1 + t in 𝔽[t]/(t^(q−1)), all q − 1 coefficients | 180 configurations over six pairs | all agree |
| 3 | `window.c` against the Python reference | 105 configurations | all agree |
| 4 | Theorem C(iii): λ_0 = λ_1 = 0 and λ_2 equals the product formula, for a singular carrier at two consecutive levels | 540 (configuration, chart) pairs over (15, 7), (10, 11), (14, 13), (14, 29); the product is nonzero in 440 | all agree |
| 5 | Theorem C(i): λ_j = 0 for j < w, for blocks of every multiplicity and position | 880 (configuration, chart) pairs on (15, 7) and (10, 11) | all vanish; λ_w ≠ 0 in most (per block in the output) |
| 6 | Proposition F on random Ferrers sets in ℤ_15 × ℤ_7, ℤ_21 × ℤ_5, ℤ_35 × ℤ_3, ℤ_10 × ℤ_11, ℤ_7 × ℤ_15; the identity of Theorem D on random sets with r nested levels; nonvanishing for peelable members of the order-105 family | 150 + 90 + 40 sets | all hold |
| 7 | Proposition E on random sets containing S × B | 75 sets | all hold |

Check 2 is the test of Lemma 1′ and Theorem A(ii): it compares all coefficients below t^(q−1), not only the leading one, on configurations with e ≤ q − 2. Check 1 tests in particular that λ_1 = Λ̄ when A₀ ∉ 𝔮. The valuation statements of Theorem A(iv) and Theorem C(ii) for orders beyond what 𝔽[t]/(t^(q−1)) shows are proved above but not tested by `validate.py`. ⟨REFEREE2⟩

## Order 105

ℤ₁₀₅ = ℤ₁₅ × ℤ₇ and C = {3, 5, 6, 9, 10, 12}. As in the earlier reports, P(105) is equivalent to D_K ≠ 0 for the K with some level set equal to C; in both charts above 7 the singular subsets of ℤ₁₅ are the translates of C and of its complement, and no level set other than C is singular for such K (no translate of C is disjoint from C). For a member of the family, the levels with T_m = C form a block s, …, s + δ − 1. `apply105.py` (output: [apply105.txt](apply105.txt)):

0. **Census.** The family has 1 784 380 694 596 407 357 786 806 016 members (2^90.527). By the multiplicity δ of C as a level set: δ = 1 in 99.331%, δ = 2 in 0.668%, δ = 3 in 0.00075%, δ ≥ 4 in 2.2 · 10⁻⁷ %. The largest classes (s, δ) are (3, 1) with 68.7% and (4, 1) with 29.8%.
1. **Multiplicity one** (uniform sample of 3 000 000 members). The second-order term is nonzero in some chart for 2 994 152 of them (99.805%). For each of the other 5 848 the third coefficient λ_2 is nonzero in at least one chart. So Theorem A(iv) proves the minor nonzero for all 3 000 000.
2. **Multiplicity two**, classes (3, 2), (2, 2) and (4, 2) (20 000 uniform members each); the other classes with δ = 2 have total share 7 · 10⁻⁶ % of the family and are in item 4. Here w = 2 and Theorem C covers the orders 2 ≤ j ≤ 6. The product λ_2 of Theorem C(iii) is nonzero in some chart for 96.6%, 90.3% and 86.0% of the members of the three classes; all others have λ_3 ≠ 0, except 30 members with least order 4, 5 or 6. All 60 000 are decided.
3. **Multiplicity three**, classes (2, 3), (3, 3), (1, 3) (1000 uniform members each): w = 4, covered orders 4 ≤ j ≤ 7, and all 3000 are decided, 2850 of them at order 4.
4. **The small classes.** Where the block of C touches the bottom or the top level, or δ ≥ 4, part of each sample is undecided:

   | class (s, δ) | share of the family | w | orders covered | undecided in a sample of 1000 |
   |---|---|---|---|---|
   | (1, 2) | 6.6 · 10⁻⁶ % | 2 | 2..6 | 0 |
   | (5, 2), top | 2.4 · 10⁻⁷ % | 3 | 3..6 | 44 |
   | (2, 4) | 2.1 · 10⁻⁷ % | 6 | 6..9 | 2 |
   | (4, 3), top | 5.0 · 10⁻⁸ % | 6 | 6..8 | 163 |
   | (1, 4) | 4.5 · 10⁻⁹ % | 6 | 6..9 | 2 |
   | (3, 4), top | 1.0 · 10⁻⁹ % | 10 | 10..11 | 311 |
   | (1, 5) | 2.0 · 10⁻¹² % | 9 | 9..11 | 23 |
   | (2, 5), top | 8.1 · 10⁻¹³ % | 15 | 15 | 729 |
   | (0, 2), bottom | 1.2 · 10⁻¹³ % | 3 | 3..6 | 199 |
   | (0, 3), bottom | 4.9 · 10⁻¹⁴ % | 6 | 6..8 | 471 |
   | (0, 4), bottom | 3.8 · 10⁻¹⁵ % | 10 | 10..11 | 755 |
   | (0, 5), bottom | 3.3 · 10⁻¹⁷ % | 15 | 15 | 982 |

   The classes (1, 6) and (0, 6) are outside Theorem C (m_max = q − 1) and were not sampled; (0, 7) is the single set C × ℤ₇.
5. **Totals.** 3 075 000 members were sampled. The 3 060 000 from the classes with share above 0.001% are all decided. 3 681 members of the small classes are not. Weighted by class size, the undecided share seen in these samples is 1.9 · 10⁻⁸ % of the family, roughly 2^58 index sets.

"Decided" means: some window coefficient in the range covered by Theorem A or Theorem C is nonzero in at least one of the two charts, so the minor is proved nonzero. "Undecided" means only that these theorems do not apply; it does not mean that the minor vanishes.

The figure in item 5 counts what uniform sampling can see. It is not an upper bound for what is left: the window-blind members described in the next section have relative size far below 1/3 000 000 inside the large classes and do not show up in such samples.

## Order 110

For base 10 above 11 there are four charts, and the chain {5, 8, 9} ⊂ {5, 7, 8, 9} ⊂ {3, 5, 6, 7, 8, 9} ⊂ {2, 3, 5, 6, 7, 8, 9} has one member singular in each. A **maximal flag** has occupancies 10, 9, …, 1 along the order 5, 8, 9, 7, 3, 6, 2, 0, 1, 4; all occupancies are distinct, every shell has one fibre, and the second-order report showed Λ ∈ 𝔮 in all four charts. `apply110.py` (output: [apply110.txt](apply110.txt)):

1. **Random maximal flags** (1000, fibres uniformly random of the prescribed sizes). The least order of a nonzero window coefficient is 2 in 815 to 847 of them, depending on the chart, and at most 6 in every chart for every flag. Over the four charts it is 2 for all 1000. By Theorem A(iv) each of these minors is nonzero, with valuation e + 2 in the chart that shows it.
2. **Nested maximal flags** (100). All window coefficients of order at most 9 vanish in all four charts, as Proposition F predicts, and the minor is nonzero by Proposition F because P(10) holds.

This bears on [Question 7.3 of Wang–Zhang](https://arxiv.org/abs/2608.22119v1): on maximal flags, where first and second order give nothing, the third-order coefficient decided every random example, and the flags it cannot decide include the nested ones, which have an exact formula. This is one chain and a sample. It is not a proof of P(110), and I have not classified the maximal flags on which λ_2 vanishes in all charts.

## What remains

A proof of P(105) along these lines needs every member of the family to be reached by some criterion. The uniform samples above are reached without exception, but a uniform sample of 2^90.5 sets cannot see a subfamily of relative size 2^(−30). `blind105.py` therefore looks at structured members directly (output: [blind105.txt](blind105.txt)).

⟨BLIND105⟩

What I can and cannot say about the undecided members:

- **No zero was found.** Every undecided member met in these runs has a nonzero minor modulo a prime ℓ ≡ 1 (mod 105), so none is a counterexample to P(105).
- **They contain a product block.** In all three families, every member whose window coefficients vanish to order 13 in both charts contains C × B with |B| = s + 1, or satisfies the same condition for its complement. I have not proved that this is necessary in general. It is not sufficient: a member that contains such a block and has many fibres outside C is normally decided at order 1 or 2 (third family in the output).
- **Why the window is blind there.** By Proposition E, D_K = d_C^(|B|) · Δ_B^6 · det(F ⊗ R_B + R_C ⊗ P_B). The factor d_C^(|B|) has valuation (q − 1)|B|. The matrix R_C has the pole 1/d_C on the fibres outside C and P_B has rank |B|, so the determinant can cancel at most as many factors d_C as the fibres outside C allow. With few elements outside C a factor d_C survives and v(D_K) ≥ e + q − 1. This is an explanation, not a proved statement.
- **The next step** is a criterion at depth one in q for the determinant in Proposition E: the residue of its coefficient of d_C^(−ρ), with ρ the rank of the pole part. For Ferrers sets and peelable sets this recovers Proposition F and Theorem D. I have not carried it out.

Two further limits:

- Theorem C does not cover a bottom or top block with δ ≥ q − 1 (m_max = q − 1). At order 105 these are the classes (s, δ) = (0, 6) and (1, 6), of total share 7.5 · 10⁻¹⁸ %.
- For δ ≥ 3 and for the third-order coefficient at multiplicity one there is no closed formula here, only the window polynomial. A proof for a whole class of index sets needs the formula, or a different argument for nonvanishing.

## Correction to the second-order report

The second-order report states that P(105) reduces to "about 2^90.5" index sets and gives the family size 1 796 332 607 012 463 666 508 846 848. That number is Σ_t #{K : T_(t−1) = C}, which counts a set once for every level at which C occurs. The number of distinct index sets is 1 784 380 694 596 407 357 786 806 016 (2^90.527). The sampler of that report drew (set, level) pairs uniformly, so its percentages are shares of pairs. In terms of distinct sets:

- C is a repeated level set in 0.669% of the family (the report says 1.33%);
- the second-order criterion decides 99.14% of the family (99.805% of the members of multiplicity one, in a sample of 3 000 000), where the report says 98.5%;
- Λ vanishes in both charts on 0.195% of the members of multiplicity one, which is 0.194% of the family (the report says 0.19%).

The theorems, the validation and the qualitative conclusions of that report are not affected.

## Reproduce

```sh
cc -O2 -o window window.c
python3 validate.py ./window > validate.txt                   # about 15 seconds
python3 apply105.py ./window 3000000 20000 > apply105.txt     # about 11 minutes
python3 blind105.py ./window 1500 > blind105.txt              # about 8 minutes
python3 apply110.py ./window 1000 > apply110.txt              # about half a minute
```

Python 3.9, standard library only; a C compiler for `window.c`. `ho_lib.py` imports `so_lib.py` from the second-order folder, and `apply105.py` imports the tables of `apply105.py` there. Random seeds are fixed in the scripts. The compiled binary is not committed.
