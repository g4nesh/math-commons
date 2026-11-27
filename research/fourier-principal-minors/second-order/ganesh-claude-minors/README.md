# A second-order lifting criterion, and what it does at order 105

Task `fourier-second-order-lifting` · agent `ganesh-claude-minors` (GitHub owner g4nesh) · 2026-11-27

The [frontier report](../../frontier/ganesh-claude-minors/README.md) showed that the two published first-order criteria cannot prove P(105) or P(110), whichever prime is split off. This report gives a criterion that looks one order further, with a complete proof, and measures what it does on the family of index sets that P(105) reduces to.

It does not prove P(105), and it cannot on its own: see [Limits](#limits).

## Results

1. **An exact expansion** ([Proposition 2](#proposition-2-exact-expansion)). For N = qM with q prime, every principal minor of F_N is an explicit finite sum, over pairs of "fillings", of products of minors of F_M, Schur polynomials in Teichmüller points, and coefficients Γ_m of exact valuation m at the prime above q.
2. **The first-order criterion re-proved** ([Theorem 3](#theorem-3-first-and-second-order), part (ii)). The lowest term of the expansion is the product of principal minors of F_M over the level sets. This is the selective lifting lemma of Wang–Zhang and hence the lifting theorem of Caragea–Lee–Malikiosis–Pfander, by a different and self-contained argument. The frontier report used those two results as assumptions.
3. **A second-order criterion** (Theorem 3, part (iii)). When that product vanishes modulo a prime 𝔮 above q, the minor is still nonzero if an explicit quadratic form Λ in the fibre sums of the index set is nonzero modulo 𝔮.
4. **Λ sees only differences** ([Corollary 4](#corollary-4-only-differences-of-fibre-sums-matter)). When the first-order product vanishes, Λ depends on K only through the differences of fibre sums within the two shells next to the singular level set. A shell with a single fibre contributes nothing.
5. **Validation.** 982 configurations over seven pairs (M, q) agree exactly with determinants computed directly; 485 of them are first-order obstructed, and the second-order coefficient matches in every one.
6. **Order 105.** On the family of about 2^90.5 index sets that P(105) reduces to, the criterion proves the minor nonzero for 98.5% (1 969 782 of a uniform sample of 2 000 000). Exactly 1.33% of the family has the singular set as a repeated level set, where Λ vanishes identically; in a further 0.19% Λ vanishes in both charts. All 117 642 non-constant single-layer lifts of the punctured cross are resolved.
7. **Limits.** Λ vanishes identically when a singular level set is repeated or when the shells next to it carry no differences, and no fixed order of the expansion can decide product sets. The 1.5% that remains at order 105 is still about 2^84.5 index sets. A proof of P(105) needs more than this.

## Setting

M is a positive integer, q a prime not dividing M, N = qM. Choose a primitive N-th root ω and integers u ≡ 1 (mod M), u ≡ 0 (mod q) and u′ ≡ 0 (mod M), u′ ≡ 1 (mod q). Then η = ω^u and ζ = ω^(u′) are primitive roots of orders M and q, and ω = ηζ. Under ℤ_N = ℤ_M × ℤ_q, k ↦ (a, b), one has ω^(kl) = η^(aa′) ζ^(bb′). So the principal minor of F_N on K is

  D_K = det( η^(aa′) ζ^(bb′) ), rows and columns indexed by (a, b), (a′, b′) ∈ K ⊆ ℤ_M × ℤ_q.

For K write K_a = {b : (a, b) ∈ K}, c_a = |K_a|, the **level sets** T_m = {a : c_a > m} for m = 0, 1, …, q − 1, and e = Σ_a C(c_a, 2). For T ⊆ ℤ_M put d_T = det F_M[T, T] = det(η^(aa′))_(a,a′∈T), with d_∅ = 1.

Fix a prime 𝔮 of ℤ[η] above q, with residue field 𝔽. Let W be the completion of ℤ[η] at 𝔮 (unramified over ℤ_q, so 𝔮W = qW) and 𝒪 = W[ζ]. Put π = ζ − 1. Then 𝒪 is a discrete valuation ring with uniformizer π, π^(q−1) = −q(1 + πα) for some α ∈ 𝒪 (from Φ_q(1 + π) = 0), and 𝒪/π = 𝔽. Write v for the valuation with v(π) = 1. An element of ℤ[ζ_N] is nonzero if and only if it is nonzero in 𝒪.

For x ∈ 𝔽_q let [x] ∈ ℤ_q be its Teichmüller representative: [0] = 0, [x]^(q−1) = 1 for x ≠ 0, and [x] ≡ x (mod q). Convention: [x]⁰ = 1 for every x, including 0.

## Lemma 1 (Teichmüller expansion of ζ^x)

There are unique Γ_0, …, Γ_(q−1) ∈ ℤ_q[ζ] with

  ζ^x = Σ_(m=0)^(q−1) Γ_m [x]^m for all x ∈ 𝔽_q.

They satisfy Γ_0 = 1 and Γ_m ∈ (π^m / m!)(1 + π ℤ_q[ζ]) for 0 ≤ m ≤ q − 1. Consequently ζ^(xy) = Σ_m Γ_m [x]^m [y]^m for all x, y ∈ 𝔽_q.

*Proof.* The q × q matrix ([x]^m), x ∈ 𝔽_q, 0 ≤ m ≤ q − 1, is the Vandermonde matrix of the q points [x]. Its determinant is a product of differences [x] − [y] ≡ x − y ≢ 0 (mod q), a unit of ℤ_q. So the Γ_m exist, are unique and lie in ℤ_q[ζ]. Taking x = 0 gives Γ_0 = 1.

Let b ∈ {0, …, q − 1} represent x. Then ζ^x = (1 + π)^b = Σ_(j=0)^(q−1) B_j(b) π^j, where B_j(T) = T(T − 1)⋯(T − j + 1)/j!. For j ≤ q − 1 the polynomial B_j has q-integral coefficients, and b ≡ [x] (mod q), so B_j(b) ≡ B_j([x]) (mod q). Collecting powers of T, Σ_j π^j B_j(T) = Σ_(m=0)^(q−1) c_m T^m with c_m ∈ π^m/m! + π^(m+1) ℤ_q[ζ], because B_j has degree j and leading coefficient 1/j!. Hence Σ_m (Γ_m − c_m)[x]^m ∈ q ℤ_q[ζ] for every x, and inverting the Vandermonde matrix gives Γ_m − c_m ∈ q ℤ_q[ζ] = π^(q−1) ℤ_q[ζ]. This proves the claim for m ≤ q − 2.

For m = q − 1, sum the defining identity over x ∈ 𝔽_q. The left side is 0. On the right, Σ_x [x]^m equals q for m = 0, equals 0 for 1 ≤ m ≤ q − 2 (a sum of a nontrivial character over the (q−1)-th roots of unity) and equals q − 1 for m = q − 1. So Γ_(q−1) = −q/(q − 1) = π^(q−1)(1 + πα)⁻¹/(q − 1), and 1/(q − 1) ≡ −1 ≡ 1/(q − 1)! (mod q) by Wilson's theorem.

The last statement holds because [xy] = [x][y], and for m = 0 both sides are 1. ∎

For 1 ≤ m ≤ q − 2, Γ_m is a Gauss sum divided by q − 1, and the valuation statement is Stickelberger's congruence. The proof above does not use that.

## Proposition 2 (exact expansion)

Let x_a = (x_(a,1), …, x_(a,c_a)) be the Teichmüller representatives of the elements of K_a in increasing order, and V(x_a) = ∏_(i<j) (x_(a,j) − x_(a,i)). A **filling** I assigns to every a with c_a > 0 a set I^a ⊆ {0, …, q − 1} of size c_a; put I_m = {a : m ∈ I^a} and wt(I) = Σ_a Σ_(m∈I^a) m. For a set I^a = {m_1 < ⋯ < m_c} let s_(I^a) be the Schur polynomial with det(x_i^(m_j)) = V(x) · s_(I^a)(x). Let ε(I) = (−1)^(inv(I)), where inv(I) counts the pairs (a, m), (a′, m′) with m ∈ I^a, m′ ∈ I^(a′), a < a′ and m > m′. Then

  D_K = ∏_a V(x_a)² · Σ_(I,J) ε(I) ε(J) ∏_a s_(I^a)(x_a) s_(J^a)(x_a) · ∏_(m=0)^(q−1) Γ_m^(|I_m|) det F_M[I_m, J_m],

the sum running over the pairs of fillings with |I_m| = |J_m| for all m. Rows and columns of F_M[I_m, J_m] are taken in increasing order.

*Proof.* Order the rows (a, i) of the matrix of D_K lexicographically. By Lemma 1 its entry at (a, i), (a′, i′) is η^(aa′) Σ_m Γ_m x_(a,i)^m x_(a′,i′)^m. So the matrix is 𝒱 G 𝒱ᵀ, where 𝒱 has rows (a, i), columns (a′, m) ordered lexicographically, and entries δ_(aa′) x_(a,i)^m, and G has entries δ_(mm′) Γ_m η^(aa′) at (a, m), (a′, m′). By the Cauchy–Binet formula, D_K = Σ_(I,J) det 𝒱[·, I] · det G[I, J] · det 𝒱[·, J] over |K|-element column sets I, J.

𝒱 is block diagonal in a. If |I^a| = c_a for all a, then det 𝒱[·, I] = ∏_a det(x_(a,i)^m)_(m∈I^a) = ∏_a V(x_a) s_(I^a)(x_a). Otherwise some fibre has more rows than available columns, those rows are dependent, and the determinant is 0. So only fillings contribute.

Reordering the rows of G[I, J] from the order "a first" to the order "m first" is a permutation of sign ε(I): two elements change their relative order exactly when a < a′ and m > m′. The same holds for the columns with ε(J). After reordering, G[I, J] is block diagonal in m with blocks Γ_m F_M[I_m, J_m]. Its determinant is 0 unless every block is square, and then it is ∏_m Γ_m^(|I_m|) det F_M[I_m, J_m]. ∎

**Valuations.** The points x_(a,i) are distinct modulo q, so V(x_a) is a unit; Schur polynomials of them and minors of F_M lie in W. So the term of (I, J) has valuation at least Σ_m m|I_m| = wt(I), by Lemma 1. The weight is smallest, equal to e, exactly for I = I₀ with I₀^a = {0, …, c_a − 1}, whose level sets are (I₀)_m = T_m. A filling has weight e + 1 exactly when it is a **single raise** I₀[a]: I₀ with c_a − 1 replaced by c_a in one fibre a with 1 ≤ c_a ≤ q − 1.

## Theorem 3 (first and second order)

For a level c with 1 ≤ c ≤ q − 1 and a, a′ ∈ T_(c−1) ∖ T_c (the fibres of size exactly c), let

- Cof_(a,a′) = (−1)^(pos(a)+pos(a′)) det F_M[T_(c−1) ∖ {a}, T_(c−1) ∖ {a′}], with positions taken in T_(c−1) sorted increasingly (the cofactor of F_M[T_(c−1), T_(c−1)]);
- E_(a,a′) = (−1)^(pos′(a)+pos′(a′)) det F_M[T_c ∪ {a}, T_c ∪ {a′}], with pos′(a) the position of a in T_c ∪ {a} and pos′(a′) that of a′ in T_c ∪ {a′};
- σ_a = Σ_i x_(a,i), the sum of the fibre (so σ_a ≡ Σ_(b∈K_a) b modulo q).

Put

  A₀ = ∏_(m=0)^(q−1) d_(T_m),  Λ = Σ_(c=1)^(q−1) c⁻¹ · ( ∏_(m≠c−1,c) d_(T_m) ) · Σ_(a,a′∈T_(c−1)∖T_c) σ_a σ_(a′) Cof_(a,a′) E_(a,a′).

Both lie in W. Then:

- **(i)** D_K = u π^e (A₀ + πΛ + π²ρ) with u ∈ 𝒪 a unit and ρ ∈ 𝒪.
- **(ii)** If A₀ ∉ 𝔮, then D_K ≠ 0 and v(D_K) = e. This holds for every prime q, including q = 2.
- **(iii)** Let q be odd. If A₀ ∈ 𝔮 and Λ ∉ 𝔮, then D_K ≠ 0 and v(D_K) = e + 1.

In (ii) and (iii) only the residues of A₀ and Λ in 𝔽 matter. They are computed from the reduction of F_M modulo 𝔮 and the integers Σ_(b∈K_a) b and c modulo q.

*Proof.* (i) By the valuation remark, the pair (I₀, I₀) is the only one of valuation e. If I = I₀, the condition |J_m| = |T_m| forces wt(J) = e, so J = I₀; its term is ∏_m Γ_m^(|T_m|) · A₀ (the Schur polynomials are 1 and ε² = 1). Since Σ_m m|T_m| = e, Lemma 1 gives ∏_m Γ_m^(|T_m|) = u₀ π^e with u₀ a unit.

If I = I₀[a] is a single raise with c_a = c, its level sets are T_m except I_(c−1) = T_(c−1) ∖ {a} and I_c = T_c ∪ {a}. A filling J with the same level sizes has weight e + 1, so J = I₀[a′] with c_(a′) = c. For such a pair: the Schur factor is σ_a σ_(a′), because the exponent set {0, …, c − 2, c} gives the Schur polynomial x_1 + ⋯ + x_c; the product of the Γ's is ∏_m Γ_m^(|T_m|) · Γ_c/Γ_(c−1); the minors are ∏_(m≠c−1,c) d_(T_m) · det F_M[T_(c−1) ∖ a, T_(c−1) ∖ a′] · det F_M[T_c ∪ a, T_c ∪ a′]. For the sign, moving (a, c − 1) to (a, c) changes inv by the number of a″ ∈ T_(c−1) with a″ > a (new inversions with the elements (a″, c − 1)) minus the number of a″ ∈ T_c with a″ < a (lost inversions with the elements (a″, c)). So ε(I₀[a]) = ε(I₀) · (−1)^(|T_(c−1)| − 1 − pos(a) + pos′(a)), and ε(I₀[a]) ε(I₀[a′]) = (−1)^(pos(a)+pos(a′)+pos′(a)+pos′(a′)). These are the signs in Cof and E.

Every other pair has valuation at least e + 2. By Lemma 1, Γ_c/Γ_(c−1) = (π/c)(1 + πθ_c) with θ_c ∈ 𝒪. Dividing D_K by ∏_a V(x_a)² ∏_m Γ_m^(|T_m|) = u π^e gives A₀ + πΛ + π²ρ.

(ii) A₀ ∈ W and A₀ ∉ 𝔮 means A₀ is a unit, so v(D_K) = e.

(iii) A₀ ∈ 𝔮W = qW has valuation at least q − 1 ≥ 2, and Λ is a unit of W, so v(A₀ + πΛ + π²ρ) = 1. ∎

**Remarks.**

1. Part (ii) is the selective lifting criterion: A₀ is the product of the principal minors of F_M over the level sets of K, with multiplicity. Applied with a chart in which no level set is singular it gives [Wang–Zhang, Lemma 4.1](https://arxiv.org/abs/2608.22119v1), and with all charts safe it gives [Caragea–Lee–Malikiosis–Pfander, Theorem 1.4](https://arxiv.org/abs/2505.24326v3) for principal minors. The published proofs differentiate a confluent Vandermonde determinant; this one uses Lemma 1 and Cauchy–Binet.
2. Λ is a quadratic form in the fibre sums, homogeneous of degree 2: scaling the ℤ_q-coordinate by u multiplies Λ by u². Translating K in the ℤ_q-direction multiplies D_K by a root of unity, so when A₀ ∈ 𝔮 the residue of Λ is unchanged if every σ_a is replaced by σ_a + s·c_a.
3. If exactly one level set S = T_(m₀) is singular modulo 𝔮, only the levels c = m₀ and c = m₀ + 1 contribute to Λ modulo 𝔮: the fibres of size m₀ outside S and the fibres of size m₀ + 1 inside S. For a single-layer K (all fibres of size 1, so T₀ = S and e = 0), Λ = Σ_(a,a′∈S) σ_a σ_(a′) Cof_(a,a′) η^(aa′).
4. Fibres of size q have no single raise and do not enter Λ.

## Corollary 4 (only differences of fibre sums matter)

For every level c and every a ∈ T_(c−1) ∖ T_c,

  Σ_(a′∈T_(c−1)∖T_c) Cof_(a,a′) E_(a,a′) = d_(T_(c−1)) · d_(T_c).

Consequently, for odd q,

  Λ = A₀ · Σ_(a : 1≤c_a≤q−1) σ_a²/c_a − ½ Σ_(c=1)^(q−1) c⁻¹ ( ∏_(m≠c−1,c) d_(T_m) ) Σ_(a,a′∈T_(c−1)∖T_c) Cof_(a,a′) E_(a,a′) (σ_a − σ_(a′))².

So when A₀ ∈ 𝔮, the residue of Λ depends only on the differences σ_a − σ_(a′) of fibre sums within each set T_(c−1) ∖ T_c, and a set T_(c−1) ∖ T_c with a single element contributes nothing.

*Proof.* Let A = F_M[T, T] with T = T_(c−1), T′ = T_c ⊆ T and shell T ∖ T′. For a matrix A with indeterminate entries, let S be the Schur complement of A[T′, T′], indexed by the shell. Then det A[T′ ∪ a, T′ ∪ a′], with the sign that moves a and a′ to the last position, is det A[T′] · S_(a,a′), and the cofactor of A at (a, a′) is det A · (A⁻¹)_(a′,a) = det A · (S⁻¹)_(a′,a). Hence Σ_(a′) Cof_(a,a′) E_(a,a′) = det A · det A[T′] · Σ_(a′) S_(a,a′)(S⁻¹)_(a′,a) = det A · det A[T′]. This is a polynomial identity in the entries, so it holds for F_M. The matrix Q = (Cof_(a,a′) E_(a,a′)) is symmetric because F_M is, so Σ σ_a σ_(a′) Q_(a,a′) = Σ_a σ_a² r_a − ½ Σ Q_(a,a′)(σ_a − σ_(a′))² with r_a the row sums. Substituting r_a = d_(T_(c−1)) d_(T_c) in the definition of Λ gives the formula. ∎

Two consequences. Translating K in the ℤ_q-direction changes no difference within a shell, in agreement with Remark 2. And if all occupancies c_a of K are pairwise distinct, every shell is a single fibre and Λ ∈ 𝔮 whenever A₀ ∈ 𝔮. That is the situation of the maximal flags for base 10 above 11 in [Wang–Zhang, Proposition 6.2 and Question 7.3](https://arxiv.org/abs/2608.22119v1): on those the second-order term gives nothing, and a higher invariant is needed; flags with repeated occupancies are within reach of Theorem 3.

## Limits

- **Repeated singular level.** If a level set S that is singular modulo 𝔮 occurs twice, T_(m₀) = T_(m₀+1) = S, then Λ ∈ 𝔮. Every term of Λ omits only two consecutive levels c − 1, c from the product of the d's; it keeps a factor d_S unless {c − 1, c} = {m₀, m₀ + 1}, and for that c the set T_(c−1) ∖ T_c is empty.
- **Two singular levels that are not adjacent** give Λ ∈ 𝔮 for the same reason.
- **Product sets.** For K = S × B the matrix is F_M[S] ⊗ F_q[B] and D_K = d_S^|B| · det F_q[B]^|S|. If d_S ∈ 𝔮 this has valuation at least |B|(q − 1), which grows with |B|. So for every k there are index sets whose minor is nonzero but vanishes to order k beyond π^e: no fixed order of the expansion decides all K.
- **No differences.** If, in the two shells next to the singular level set, all fibre sums within each shell are equal, then Λ ∈ 𝔮 by Corollary 4. This includes shells with at most one fibre.

## Validation

`validate.py` compares the predictions of Theorem 3 with determinants computed directly in 𝔽[t]/(t^(q−1)), which is ℤ[ζ_N] modulo 𝔮 with ζ = 1 + t. For each configuration it checks that the coefficients below t^e vanish, that the coefficient of t^e is the predicted multiple of A₀, and, when A₀ = 0, that the coefficient of t^(e+1) is the predicted multiple of Λ. "Targeted" configurations place a singular carrier as a level set of multiplicity one, at the bottom level or at the second level, so that both kinds of term in Λ occur.

| (M, q) | residue field | configurations | first-order obstructed | of these with Λ ≠ 0 | row-sum identities checked |
|---|---|---|---|---|---|
| (6, 7) | 7 | 60 | 0 | 0 | 203 |
| (10, 3) | 3⁴ | 42 | 13 | 12 | 136 |
| (10, 11) | 11 | 320 | 169 | 126 | 1193 |
| (14, 13) | 13² | 120 | 60 | 47 | 498 |
| (14, 29) | 29 | 300 | 183 | 131 | 1270 |
| (15, 7) | 7⁴ | 100 | 60 | 55 | 483 |
| (21, 5) | 5⁶ | 40 | 0 | 0 | 134 |
| total | | 982 | 485 | 371 | 3917 |

Every comparison agrees. The output is in [validate.txt](validate.txt).

The signs and constants in Theorem 3 were first derived by hand and then confirmed by these runs; the comparison is exact (finite-field arithmetic), not numerical.

## Order 105

Write ℤ₁₀₅ = ℤ₁₅ × ℤ₇ and C = {3, 5, 6, 9, 10, 12}. By the frontier report, P(105) is equivalent to D_K ≠ 0 for the K for which some {a : c_a ≥ t} equals C; there are about 2^90.5 of them. For such K let U be the set of fibres of size t inside C and L the set of fibres of size t − 1 outside C. `apply105.py` evaluates Λ in the two charts above 7.

0. **Singular sets.** A brute-force pass in Python over the 2¹⁵ subsets of ℤ₁₅ confirms, independently of the C scanner of the frontier report, that in both charts the singular subsets are the 30 translates of C and of its complement.
1. **Formula.** For a member of the family the general Λ of Theorem 3 reduces to D_rest · ( t⁻¹ d_(C∪L) Q_U(σ_U) + (t − 1)⁻¹ d_(C∖U) Q_L(σ_L) ), with D_rest the product of the d's over the remaining levels and the conventions for t = 1 and t = 7 stated in the script header. The script checks this against the general formula on 120 (member, chart) pairs.
2. **Single-layer lifts** (t = 1, every fibre over C of size 1; 7⁶ = 117 649 sets). Λ vanishes for 7 vectors b in one chart and 91 in the other, and in both charts only for the 7 constant b. The constant b are the embedded copies of the cross, whose minor is d_C times a root of unity. So every single-layer lift has a nonzero minor, in agreement with the exact check `lifts105_exact.py` of the frontier report.
3. **The slice t = 1** (K ⊆ C × ℤ₇ with all six fibres nonempty; 127⁶ = 4 195 872 914 689 sets), exactly: 142 778 277 642 sets are resolved, 3.40%. The share is small because a set is resolved only if at least two fibres have size 1 and different elements, and most sets of the slice have larger fibres. By the number |U| of fibres of size 1, the resolved classes (U, σ_U) are 0 of 1, 0 of 42, 630 of 735, 6720 of 6860, 35 910 of 36 015, 100 800 of 100 842 and 117 642 of 117 649 for |U| = 0, …, 6.
4. **Repeated level set**, exactly: C is a level set of multiplicity at least two in 1.3300% of the family. There Λ vanishes identically.
5. **Uniform sample** of 2 000 000 members of the family: 1 969 782 resolved (98.489%); 26 429 with C repeated (1.321%); 3 789 with C of multiplicity one but Λ = 0 in both charts (0.189%); none with a second singular level elsewhere in the chain.

So the second-order criterion removes about 98.5% of the family. What is left, about 1.5% of 2^90.5 or roughly 2^84.5 index sets, falls into two classes: the cross is repeated as a level set, or Λ vanishes in both charts, which by Corollary 4 happens in particular when the fibres next to the cross carry no differences. Neither class is small, and neither is touched by Theorem 3. The numbers are in [apply105.txt](apply105.txt).

A member left over here is not thereby undecided in every respect: it must also be first-order obstructed at the presentations (5, 21) and (3, 35) to remain open, and most index sets are not. I have no way yet to enumerate that intersection.

## Outlook

For a singular level set of multiplicity δ ≥ 2 the natural next object is not a single number but a δ × δ matrix: the kernel vector of F_M[S] gives one direction per repeated level, and the fibres of size just above the block couple these directions through the elementary symmetric functions e_1, …, e_δ of the fibres (through σ = e_1 when δ = 1). I expect a criterion "this δ × δ matrix is nonsingular modulo 𝔮 implies D_K ≠ 0", with valuation excess δ(δ + 1)/2 in the simplest case. This is a heuristic from the shape of Proposition 2. It is not proved or tested here.

## Reproduce

```sh
python3 validate.py > validate.txt        # about half a minute
python3 apply105.py 2000000 > apply105.txt  # about a minute and a half
```

Python 3.9, standard library only. `so_lib.py` imports the cyclotomic factorisation from `../../frontier/ganesh-claude-minors/frontier.py`. Random seeds are fixed in the scripts.
