# In dimension five and above, squared minors of a unitary need not be a mixture of orthogonal ones: full proof and independent check

- **What this is.** A complete write-up and independent check of the claim made in ledger task `mdo-orthogonal-hull-five`.
- **Credit.** The functional Φ, the unitary V and its exact properties (their Proposition A), the claim r₅(V) ∉ conv r₅(O(5)) (their Theorem N5), the face argument for n ≥ 5 (their Corollary N), the statement "r_n(U(n)) ⊆ conv r_n(O(n)) if and only if n ≤ 4" with its small-n argument (their Corollary Z) and the description of V's spectrum are the work of agent `ganesh-claude-minors` (GitHub owner `g4nesh`). Their folder is `research/marcus-de-oliveira/orthogonal-hull-five` on branch `claude/mdo-orthogonal-hull-five`; commit `18501ba` was read here as data. None of their code was run.
- **Who checked.** Ledger agent `ishaan-air` (GitHub owner `IshaanAyaan`), November 25, 2026. This folder contains three component checks, `h5-verify` ([`check_exact.py`](check_exact.py)), `h5-bound` ([`bound/`](bound/README.md)) and `h5-structure` ([`structure.md`](structure.md)), plus this write-up. For the write-up, every exact checker was re-run and the reports were reconciled.
- **Status.** A computer-assisted proof. Every step of Theorem 1 is either exact arithmetic or a proof given here; the n = 4 case of Corollary 2 is cited. It has **not yet been reviewed by a third GitHub owner**.
- **Levels.** Every claim carries one of three labels. **EXACT** means a complete proof, or a computation in exact rational or algebraic-number arithmetic. **NUMERIC** means floating-point evidence only. **UNPROVED** means a conjecture, or a statement cited and not checked here.

## 1. The result

> **Theorem 1 (EXACT, computer-assisted).** Let Φ and V be as in §2–§3.
> 1. V is unitary and Φ(V) = 37/25.
> 2. Φ(O) ≤ B < 11371/10000 for every real orthogonal 5×5 matrix O, both components included. Here B = 113702917/10⁸ = 1.13702917 exactly.
> 3. Hence r₅(V) ∉ conv r₅(O(5)). More generally, r₅(U) ∉ conv r₅(O(5)) for every U in the open set {U ∈ U(5) : Φ(U) > 11371/10000}, which contains V. At V the linear functional Φ separates by at least 37/25 − 11371/10000 = 0.3429.
> 4. For every n ≥ 5, r_n(V ⊕ I_{n−5}) ∉ conv r_n(O(n)).

> **Corollary 2 (Corollary Z of `ganesh-claude-minors`).** The inclusion r_n(U(n)) ⊆ conv r_n(O(n)) holds if and only if n ≤ 4.
> - n ≥ 5: this is Theorem 1.4 (EXACT).
> - n ≤ 3: an elementary argument, §5.3 (EXACT).
> - n = 4: Theorem 1.2 of D. Zheng, "The Marcus–de Oliveira conjecture for n = 4", preprint, October 2026, DOI 10.13140/RG.2.2.36790.25923 (cited, **not re-proved here**). For n = 4 the corollary is therefore conditional on that theorem.

For n ≥ 5 this answers the question from Section 7 of Zheng's preprint that the other owner's README cites: no n ≥ 5 has the inclusion (the wording of that question was not checked against the preprint; see §6). It is **not** a counterexample to the Marcus–de Oliveira conjecture (§6).

**Proved along the way (all EXACT).**
- Φ(O) = 1 − 6·det(O)·E_S(O) on O(5), and Φ(U) = 1 − 6·conj(det U)·E_S(U) on U(5) (§4.1).
- The exact structure of the 60-element set S (§4.3).
- V is a strict local maximum of Φ on U(5) modulo diagonal phases (§3).
- The sharp bound Φ ≤ 1 holds on the 25 six-dimensional families of O ∈ O(5) that have an entry ±1 (§4.4).
- Φ is not, up to an affine change, a determinant functional of the kind that appears in the Marcus–de Oliveira conjecture (§6).

**NUMERIC only.**
- max Φ over O(5) is 1, the value at the 60 permutation matrices outside S; the maximisers form a 6-dimensional set.
- max Φ over U(5) is 37/25.
- At V, det(A + B) lies in the Marcus–de Oliveira hull for every spectrum tested.

**UNPROVED.**
- The sharp bound Φ ≤ 1 on O(5), equivalently E_S ≥ 0 on SO(5).
- Global maximality of 37/25 on U(5).
- The other owner's own certificate for Φ ≤ 13/10. It was not checked here. Theorem 1 does not use it.

## 2. Definitions

**Squared-minor vector.** For an n×n complex matrix M and index sets I, J ⊆ {0, …, n−1} with |I| = |J| (the empty set included), let det M[I,J] be the minor with rows I and columns J. Then

  r_n(M) = ( |det M[I,J]|² )_{|I|=|J|} ∈ ℝ^{C(2n,n)}.

For n = 5 this vector has 252 coordinates. We write O(n) for the real orthogonal and U(n) for the unitary n×n matrices, and conv for the convex hull. The hull conv r_n(O(n)) is compact, because O(n) is compact and r_n is continuous.

**The functional Φ** (due to `ganesh-claude-minors`). Indices are taken in ℤ₅. For a 5×5 matrix M let

  Φ(M) = − Σ_v |m_vv|²
     + Σ_v ( |det M[{v+1,v−1},{v+2,v−2}]|² + |det M[{v+2,v−2},{v+1,v−1}]|² )
     − Σ_v Σ_{d=1,2} ( |det M[{v,v+d},{v,v−d}]|² + |det M[{v,v−d},{v,v+d}]|² ).

It is a linear function of r₅(M) with 35 coefficients ±1 on 35 distinct coordinates (EXACT, [`check_exact.py`](check_exact.py) Part A).

Φ comes from the cyclic one-factorisation of K₆ on ℤ₅ ∪ {∞}. The pairs {v+1, v−1} and {v+2, v−2} form a perfect matching M_v of ℤ₅ ∖ {v}. Adding the edge {v, ∞} to M_v gives one of the five one-factors of that one-factorisation.

**Symmetries (EXACT).** Φ is invariant under the four maps below. The first two are checked in Part A; the last two leave every squared minor unchanged.
- M ↦ P_g M P_gᵀ for g in the affine group AGL(1,5) = F₂₀ = {x ↦ ax + b};
- transposition;
- complex conjugation;
- M ↦ D₁ M D₂ with D₁, D₂ diagonal unitary.

**Permutation matrices.** P_σ has entries (P_σ)_{i,σ(i)} = 1. The minor det P_σ[I,J] equals ±1 when σ(I) = J and 0 otherwise.

**The set S and the alternating sum E_S.** Let

  S = {σ ∈ S₅ : Φ(P_σ) = −5}  and  E_S(M) = Σ_{σ∈S} sgn(σ) Π_i M_{i,σ(i)}.

For a set T ⊆ S₅ write E_T in the same way. Also write per_H(M) = Σ_{σ∈H} Π_i M_{i,σ(i)}. D₅ denotes the dihedral subgroup {x ↦ ±x + b} of F₂₀.

## 3. The unitary V and the exact value Φ(V) = 37/25

The unitary and the value below were found by `ganesh-claude-minors` (their Proposition A). With s = √6 and t = i√15,

```
            ⎡ 4    4s     4s     4s     4s   ⎤
            ⎢ 4s  −1−t   −6−2t  −6+2t   9+t  ⎥
V = 1/20 ·  ⎢ 4s  −6−2t  −1+t    9−t   −6+2t ⎥
            ⎢ 4s  −6+2t   9−t   −1+t   −6−2t ⎥
            ⎣ 4s   9+t   −6+2t  −6−2t  −1−t  ⎦
```

The matrix was transcribed twice, independently: into [`V_exact.json`](V_exact.json) as coordinates in the basis 1, s, t, st of ℚ(√6, i√15), and into [`bound/exact_phi_at_V.py`](bound/exact_phi_at_V.py) as sympy. Both transcriptions agree with the matrix above and with the other owner's README.

> **Proposition 3 (EXACT).** V·V* = I and Φ(V) = 37/25.

*Check.* Two independent implementations give the same answer:
- [`check_exact.py`](check_exact.py) Part C does the arithmetic in ℚ(√6, i√15) with the standard library only. It also includes a sympy cross-check; that cross-check ends with `nsimplify`, so it is not exact by itself, and the exact evidence is the field arithmetic here and the script below.
- [`bound/exact_phi_at_V.py`](bound/exact_phi_at_V.py) uses sympy.

In both, every entry of VV* − I and the difference Φ(V) − 37/25 reduce to 0.

**Further exact facts about V (EXACT).** All three are stated in the other owner's Proposition A and "What V is". Item 1 is also checked by [`check_exact.py`](check_exact.py) Part C; all three are checked by the snippet in the appendix.
1. The squared moduli of the entries are 1/25 on the diagonal and 6/25 off it. So V has no entry of modulus 1, which matters in §4.4.
2. The squared 2×2 minors take only the values 3/125, 8/125, 12/125, 18/125 and 27/125.
3. V is symmetric, tr V = 0, and (V − I)(V + I)(V² + V + I) = 0.

From item 3, the eigenvalues lie in {1, −1, ω, ω²} with ω = e^{2πi/3}. The imaginary part of the trace forces ω and ω² to have equal multiplicity, and then the real part of the trace forces the spectrum (1, 1, −1, ω, ω²). Write V = X + iY with X and Y real symmetric. Unitarity gives XY = YX, so X and Y share a real orthonormal eigenbasis, and

  V = Q·diag(1, 1, −1, ω, ω²)·Qᵀ  with Q ∈ O(5).

This description is due to `ganesh-claude-minors`; the deduction is re-done here by hand. In particular **det V = −1**, so V ∉ SU(5) and −V ∈ SU(5). Also E_S(V) = 2/25 and E_S(−V) = −2/25.

> **Proposition 3′ (EXACT).** V is a critical point of Φ on U(5). The Hessian there is negative semidefinite, and its kernel is exactly the 9-dimensional tangent space of the orbit {D₁VD₂ : D₁, D₂ diagonal unitary}. Hence V is a strict local maximum of Φ on U(5) modulo that orbit.

*Check.* [`check_exact.py`](check_exact.py) Part D computes the exact first and second Taylor coefficients of ε ↦ Φ(exp(εX)V) over a basis of 𝔲(5) with entries in the field. All 25 first coefficients vanish. The Hessian annihilates the 10 orbit tangent vectors, which span 9 dimensions. On a complementary 16-dimensional subspace, an exact LDLᵀ over ℚ(√6) has 16 negative pivots. Φ is constant on the orbit, so the Morse–Bott argument gives strict local maximality transversally.

Global maximality of 37/25 on U(5) is NUMERIC only:
- 1,400 local maximisations in `h5-verify` (1,000 from Haar starts and 400 from starts near O(5)) and 40 more from Haar starts for this write-up all ended at 1.4800000000.
- The 40 maximisers sampled in `h5-verify` all had the form D₁PV′QD₂, with V′ ∈ {V, V̄} and P, Q permutation matrices.

## 4. The bound on O(5)

### 4.1 The Leibniz form of squared-minor functionals

> **Lemma 4 (EXACT).** Let M be invertible, I, J ⊆ {0,…,n−1} with |I| = |J|, and ε = (−1)^{ΣI+ΣJ}. Then
>   Σ_{σ : σ(I)=J} sgn(σ) Π_i M_{i,σ(i)} = det M · det M[I,J] · det(M⁻¹)[J,I].
> In particular the sum equals det(O)·det O[I,J]² for O ∈ O(n), and det(U)·|det U[I,J]|² for U ∈ U(n).

*Proof.* The generalised Laplace expansion along the rows I groups the Leibniz terms of det M by the set σ(I). The group with σ(I) = J sums to ε·det M[I,J]·det M[Iᶜ,Jᶜ]. Jacobi's complementary-minor identity gives det M[Iᶜ,Jᶜ] = ε·det M·det(M⁻¹)[J,I]. Multiply, and use ε² = 1. For O ∈ O(n) take M⁻¹ = Oᵀ, so that det(Oᵀ)[J,I] = det O[I,J]. For U ∈ U(n) take M⁻¹ = U*, so that det(U*)[J,I] is the complex conjugate of det U[I,J]. ∎

**Machine checks.**
- EXACT: all 252 pairs (I,J) on six rational orthogonal matrices, three of each determinant ([`check_exact.py`](check_exact.py) Part B).
- NUMERIC, orthogonal: the consequence Φ(O) = 1 − 6·det(O)·E_S(O) (Corollary 5) on 200 Haar matrices, error 7·10⁻¹⁵ ([`bound/sanity_numeric.py`](bound/sanity_numeric.py)).
- NUMERIC, unitary form: 50 Haar unitaries, error 3·10⁻¹⁵, checked for this write-up.

> **Corollary 5 (EXACT).** Let L = Σ c_{IJ}·|det M[I,J]|² and put v_σ = L(r(P_σ)) = Σ_{σ(I)=J} c_{IJ}. Then
>   L(O) = det(O)·Σ_σ sgn(σ)·v_σ·Π_i O_{i,σ(i)} on O(n),  and  L(U) = conj(det U)·Σ_σ sgn(σ)·v_σ·Π_i U_{i,σ(i)} on U(n).
> Φ takes only the values 1 and −5 on permutation matrices (§4.3), so v = 1 − 6·1_S and
>   Φ(O) = 1 − 6·det(O)·E_S(O) on O(5),  Φ(U) = 1 − 6·conj(det U)·E_S(U) on U(5).

*Proof.* Sum Lemma 4 against c_{IJ}, and use |det U| = 1. ∎

**The determinant −1 component (EXACT).** Let D = diag(−1,1,1,1,1). Squared minors are unchanged by a column sign, so Φ(OD) = Φ(O). Each Leibniz product contains exactly one entry of column 0, so E_S(OD) = −E_S(O). Hence the two components of O(5) give the same values of Φ. On SO(5):
- Φ ≤ 1 ⟺ E_S ≥ 0;
- Φ ≤ 13/10 ⟺ E_S ≥ −1/20;
- Φ ≤ B ⟺ E_S ≥ −(B − 1)/6 ≈ −0.022838.

At V: det V = −1, so Φ(V) = 1 + 6·E_S(V), with E_S(V) = 2/25 and E_S(−V) = −2/25 for −V ∈ SU(5).

### 4.2 The certificate for Φ ≤ B

> **Theorem 6 (EXACT, computer-assisted).** For every O ∈ O(5), Φ(O) ≤ B, where B = c + ‖R‖₁ = 113702917/10⁸ = 1.13702917 < 11371/10000.

The certificate is [`bound/certificate.json`](bound/certificate.json) (sha256 `90cdf5de…2d3c`), built by the `h5-bound` check. Write x_ij (0 ≤ i, j ≤ 4) for the 25 real entries. On real matrices Φ is a quartic polynomial in the x_ij. Let g_k be the 30 orthogonality polynomials:
- Σ_j x_ij·x_kj − δ_ik for i ≤ k;
- Σ_i x_ij·x_il − δ_jl for j ≤ l.

The file supplies three things:
- the rational number c = 1137021/1000000;
- 346 rational symmetric matrices G_b, each with a vector v_b of monomials of degree ≤ 3, in sizes 42, 26, 13, 6, 5 and 2;
- 30 rational multiplier polynomials h_k, of degree ≤ 4, with 8,310 terms in total.

All numbers are integers over the common denominator 10⁹. Define the polynomial R by

  c − Φ(x) = Σ_b v_b(x)ᵀ·G_b·v_b(x) + Σ_k h_k(x)·g_k(x) + R(x).

**What the checker establishes in exact arithmetic.** [`bound/check_certificate.py`](bound/check_certificate.py) uses the Python standard library only and shares no code with the generator.
1. It builds Φ from the definition by 1×1 and 2×2 determinants and confirms that the 35 minors are distinct. It builds the g_k from their definitions. From the file it takes only c, the pairs (v_b, G_b) and the coefficients of the h_k.
2. Every G_b is positive definite: an exact LDLᵀ without pivoting has all pivots > 0, which is Sylvester's criterion.
3. It expands the right-hand side and computes R exactly. R has 4,804 nonzero coefficients, the largest of absolute value 2.4·10⁻⁷. The sum of their absolute values is ‖R‖₁ = 8.17·10⁻⁶.
4. B = c + ‖R‖₁ satisfies B ≤ 11371/10000 < 13/10 < 37/25, by exact comparisons of fractions. In fact ‖R‖₁ = 817/10⁸ exactly, so B = 113702917/10⁸ and the first inequality is strict.

*Proof of Theorem 6.* Let O ∈ O(5).
- Each g_k(O) = 0.
- Each v_b(O)ᵀG_b v_b(O) ≥ 0.
- Every row of O is a unit vector, so |O_ij| ≤ 1. Every monomial therefore has absolute value ≤ 1 at O, and |R(O)| ≤ ‖R‖₁.

Hence Φ(O) = c − Σ_b v_b(O)ᵀG_b v_b(O) − R(O) ≤ c + ‖R‖₁ = B. ∎

Two remarks on soundness:
- The argument uses nothing about how the blocks were chosen or how the numbers were obtained. Whatever rational data the file holds, the checker computes R from that data, so a wrong file can only produce a larger B or a failed check. Positive *semi*definiteness would suffice; positive definiteness is what the checker tests.
- The relations g_k hold on both components of O(5), so the determinant −1 component needs no separate treatment.

**How it was built** (not needed for verification). The source is [`bound/generate_certificate.py`](bound/generate_certificate.py); it uses cvxpy with Clarabel and takes about 35 to 70 s. Re-running it in a separate copy reproduced `certificate.json` byte for byte (same sha256). Φ and the g_k are invariant under sign changes of single rows and columns, so the monomials of degree ≤ 3 split by the parities of their row and column degrees. A semidefinite program minimises c with all G_b ⪰ 10⁻⁵·I. Its optimum is 1.1370208. The data were rounded to 10⁻⁹, and the rounding error is exactly R.

**Independent spot checks run for this write-up** (scripts not committed). These read only the JSON fields and do not use the checker's code.
- Floating-point eigenvalues: the smallest eigenvalue over all 346 blocks is 9.998·10⁻⁶ > 0 (NUMERIC).
- Exact evaluation at 8 rational orthogonal matrices (Cayley transforms, 4 of each determinant): c − Φ(O) − Σ_b v_b(O)ᵀG_b v_b(O) has absolute value at most 2.39·10⁻⁷ ≤ ‖R‖₁, as the identity requires (EXACT at those points).
- The same quantity in floating point at 300 Haar matrices: at most 2.47·10⁻⁷ (NUMERIC).

**Relation to the other owner's certificate.** The branch of `ganesh-claude-minors` contains a different certificate. According to its README it proves 13/10 − Φ = SOS + (an exact element of the orthogonality ideal), with 225 blocks, a 40-element symmetry group and a 37 kB file. Its README marks it as not independently reviewed. It was **not run or audited here** (UNPROVED as far as this check goes). Theorem 1 uses only Theorem 6.

**Limit of the method (NUMERIC).** Degree-6 certificates of this shape cannot go below about 1.1365:
- the other owner reports 1.1365 without an eigenvalue margin (their log: 1.136544, solver status `optimal_inaccurate`, for a symmetry-reduced shape);
- here the optimum is 1.13702 with margin 10⁻⁵;
- an exploratory SO(5) variant that adds the relations O = cof(O) and det O = 1 reached 1.13654.

The sharp value 1 is out of reach at this degree.

### 4.3 The set S (EXACT)

All of the following is integer enumeration over S₅ ([`check_exact.py`](check_exact.py) Part A; the D₅ double-coset count was redone for this write-up).

> **Proposition 7 (EXACT).** Φ(P_σ) ∈ {1, −5}, and each value is taken on exactly 60 permutations. The set S where Φ(P_σ) = −5 is
>   S = D₅ ∪ {all 10 transpositions} ∪ {all 20 three-cycles} ∪ {the 20 four-cycles not in F₂₀}.

**Counts by cycle type.**

| cycle type | in S | in the complement Sᶜ |
|---|---|---|
| identity | 1 | 0 |
| (2) | 10 (all) | 0 |
| (2,2) | 5 (the reflections x ↦ −x + b) | 10 |
| (3) | 20 (all) | 0 |
| (3,2) | 0 | 20 (all) |
| (4) | 20 | 10 (the maps x ↦ 2x + b and x ↦ 3x + b) |
| (5) | 4 (the translations x ↦ x + b, b ≠ 0) | 20 |

**Further properties.**
- S has 30 even and 30 odd elements, so S ≠ A₅. It contains the identity and the transposition (0 1).
- S is closed under inversion, but it is not a subgroup.
- {g : gSg⁻¹ = S} = F₂₀, and S ∩ F₂₀ = D₅.
- Only g = e satisfies gS = S or Sg = S.
- F₂₀ has two double cosets, of sizes 20 and 100. S meets both and contains neither.
- S is not a union of D₅ double cosets: D₅·(0 1)·D₅ has 50 elements, of which 30 lie in S.

**Closed forms (EXACT, from [`structure.md`](structure.md) §2 and its appendix check (1)).** Let θ be the action of S₅ on the six one-factorisations of K₆. Then

  sgn(σ)·1_S(σ) = 1/6 + sgn(σ)/6 + sgn(σ)·fix(σ)/3 + 1[σ ∈ F₂₀] − fix(θσ)/3.

Consequently, on SO(5),

  E_S(O) = 1/6 + ⅓·Σ_i O_ii² + per(O)/6 + per_{F₂₀}(O) − ⅓·Σ_{a=0}^{5} per_{H_a}(O),

where H₀ = F₂₀, …, H₅ are the six conjugates of F₂₀.

**Positivity certificates of a simple kind cannot work (EXACT).** 1_S is not a nonnegative combination of the indicators of {σ : σ(I) = J}. If it were, Corollary 5 would write E_S(U) on SU(5) as a nonnegative combination of squared minors. But E_S(−V) = −2/25 < 0. So any proof of E_S ≥ 0 on SO(5) must use more than nonnegative combinations of squared minors.

### 4.4 The sharp bound on the matrices with an entry ±1 (EXACT)

These results are from [`structure.md`](structure.md) §3.1. The identities were re-run for this write-up (appendix checks (3) and (4) there).

> **Theorem 8 (EXACT).** Let O ∈ SO(5), and take all indices in ℤ₅.
> - (a) If O e_v = e_v, put A = {v±1} and B = {v±2}. Then E_S(O) = 1 − det O[A,B]² = Σ_{|J|=2, J≠B} det O[A,J]². Hence Φ(O) = −5 + 6·det O[A,B]² ∈ [−5, 1].
> - (b) If O e_v = e_w with w ≠ v, put c = (v+w)/2 and d = (c+v)/2. Then E_S(O) = O_cc² + det O[{d,v},{c,d}]². Hence Φ(O) ≤ 1.
> - (c) Consequently Φ(O) ≤ 1 for every O ∈ O(5) that has an entry ±1.

*Proof of (a), for v = 4.* Then O = diag(Q, 1) with Q ∈ SO(4). Write Q as x ↦ p·x·q̄ with unit quaternions p, q. An exact sympy expansion gives the identity of forms of bidegree (4,4)

  E_S = (|p|²|q|²)² − t²,  where t = ½(det O[A,B] + det O[B,A]).

Jacobi's identity in SO(4) gives det Q[A,B] = (−1)^{ΣA+ΣB}·det Q[B,A]. Since ΣA + ΣB = 0+1+2+3 is even, t = det O[A,B]. The rows A of O are orthonormal, so Cauchy–Binet gives Σ_J det O[A,J]² = 1.

*Proof of (b), for (w, v) = (0, 4).* The same parametrisation, applied to the block O[{1,2,3,4},{0,1,2,3}], gives the exact identity E_S = |p|²|q|²·O₂₂² + det O[{3,4},{2,3}]².

*Remaining cases.* Conjugation by F₂₀ preserves S, hence E_S. F₂₀ is sharply 2-transitive and its maps preserve midpoints, so (a) and (b) transfer from the cases above to every v and every pair (w, v). For (c), let O ∈ O(5) have O e_v = ±e_w, and let D_k be the sign change of coordinate k. Row and column sign changes leave Φ unchanged. If O e_v = −e_w, replace O by D_x·O·D_v with x ≠ w. This keeps det O and gives O e_v = e_w. If det O = −1, replace O by O·D_y with y ≠ v. This keeps O e_v and gives det O = +1. Then (a) or (b) applies. ∎

**Unitary version (EXACT, conditional on Zheng's Theorem 1.2, which is cited and not re-proved).** Let U ∈ U(5) have an entry of modulus 1. Then r₅(U) is a linear image of r₄ of the complementary 4×4 unitary block. By Zheng's theorem that r₄ lies in conv r₄(O(4)), and Theorem 8 then gives Φ(U) ≤ 1. So every unitary separated by Φ has all entries of modulus < 1, as V does: its squared moduli are 1/25 and 6/25.

### 4.5 The sharp bound Φ ≤ 1 on O(5): numerical picture (NUMERIC)

- **Maximum.** Local maximisations all end at 1. In `h5-verify`, 2,000 per component ended between 1.0000000000000 and 1.0000000000002. [`bound/sanity_numeric.py`](bound/sanity_numeric.py) ran 20, and 20 more ran with a new seed for this write-up. All of these ended between 1 − 2·10⁻¹¹ and 1. The largest value over 20,000 Haar samples is 0.992, and 0.977 with the new seed.
- **Maximiser set.** At 80 random maxima the 10×10 Hessian has exactly 6 eigenvalues below 6·10⁻⁶ in absolute value, and the others are at least 0.62. PCA of nearby maximisers agrees. None of the 80 is a signed permutation matrix.
- **Permutation maximisers.** At the 60 permutation maximisers the Hessian has rank 4 (the 20 with a fixed point) or rank 3 (the 40 without). The second-order expansion behind this is EXACT ([`structure.md`](structure.md) §3.2). At the rank-3 points E_S grows like t⁴ along the kernel.
- **Shape of the zero set.** The zero set of E_S is not a coset of a subgroup, and E_S has no continuous symmetry.
- **Minimum.** The minimum of Φ is −5, on both O(5) and U(5).

## 5. The separation theorem and the extension to n ≥ 5

### 5.1 Dimension five

> **Theorem 9 (EXACT, computer-assisted; Theorem N5 of `ganesh-claude-minors`).** r₅(V) ∉ conv r₅(O(5)). More precisely, the linear functional Φ satisfies
>   sup over conv r₅(O(5)) of Φ ≤ B < 11371/10000 < 37/25 = Φ(V).

*Proof.* Φ is linear in r₅, so its supremum over the convex hull equals its supremum over r₅(O(5)), which is at most B by Theorem 6. Proposition 3 gives Φ(V) = 37/25. ∎

Φ is continuous, so every U with Φ(U) > 11371/10000 is separated as well; this is an open set containing V. The set of all U ∈ U(5) with r₅(U) ∉ conv r₅(O(5)) is also invariant under M ↦ D₁PMQD₂, for permutation matrices P, Q and diagonal unitaries D₁, D₂: the diagonal factors leave r₅ unchanged, and M ↦ PMQ permutes the coordinates of r₅ and maps O(5) onto itself, so it maps the hull onto itself.

### 5.2 Every n ≥ 5

> **Corollary 10 (EXACT; the argument is that of `ganesh-claude-minors`, re-checked step by step).** For every n ≥ 5, W = V ⊕ I_{n−5} satisfies r_n(W) ∉ conv r_n(O(n)).

*Proof.* Let R = {5, …, n−1}, and suppose r_n(W) = Σ_k w_k·r_n(O_k) with w_k > 0, Σ_k w_k = 1 and O_k ∈ O(n).
1. The (R,R) coordinate of r_n(W) is 1. For orthogonal O, Hadamard's inequality and the fact that rows have norm 1 give det O[R,R]² ≤ Π_{i∈R} ‖row_i(O[R,R])‖² ≤ 1. So every O_k has det O_k[R,R]² = 1.
2. Equality in step 1 forces every row of O_k[R,R] to have norm 1. So the rows R of O_k vanish outside the columns R, and O_k[R,R] is orthogonal. Its columns then have norm 1, so the columns R vanish outside the rows R. Hence O_k = A_k ⊕ B_k with A_k ∈ O(5).
3. For I, J ⊆ {0,…,4}, the (I,J) coordinates then give r₅(V) = Σ_k w_k·r₅(A_k), which contradicts Theorem 9. ∎

The argument needs no bound on O(n) beyond the n = 5 separation. Here Theorem 6 provides that separation, so Corollary 10 does not depend on the other owner's unreviewed certificate.

### 5.3 Small n

> **Corollary 11.**
> - (a) n ≤ 3 (EXACT): r_n(U(n)) ⊆ conv r_n(O(n)).
> - (b) n = 4 (cited): the inclusion is Zheng's Theorem 1.2. The other owner's folder `unitary-factorization` (Corollary F1) also gives it. Neither proof was checked here.

*Proof of (a)* (the argument of `ganesh-claude-minors` in their Corollary Z, written out). For n = 3, Jacobi's identity gives |det U[I,J]| = |U_{i′j′}| for 2×2 minors, where {i′} = Iᶜ and {j′} = Jᶜ, and |det U| = 1. So r₃(U) = F(X), where X = (|U_ij|²) is doubly stochastic and F is an affine map that is the same for every unitary. Birkhoff's theorem gives X = Σ λ_k P_k. Since P_k is orthogonal with squared-entry matrix P_k, we get r₃(U) = Σ λ_k F(P_k) = Σ λ_k r₃(P_k). For n = 2, r₂(U) = r₂(R_θ) for the rotation with cos²θ = |u₁₁|². The case n = 1 is trivial. ∎

Theorem 1.4 and Corollary 11 together give Corollary 2.

## 6. What this means, and what it does not mean

**For Zheng's question.** The other owner's README presents Corollary 2 (their Corollary Z) as the answer to a question in Section 7 of Zheng's preprint (reference in §1), about the n ≥ 5 for which every r_n(U) lies in conv r_n(O(n)). I did not re-read the preprint, so the exact wording of that question is not checked here. For that question as described, the answer is: for none (Theorem 1.4).

**For the Marcus–de Oliveira conjecture (MOC).** Let A = W_A·D_a·W_A* and B = W_B·D_b·W_B* be normal, and put U = W_A*·W_B. Then det(A + B) = det(D_a + U·D_b·U*). Cauchy–Binet gives

  det(D_a + U·D_b·U*) = Σ_{|K|=|L|} a_{Kᶜ}·b_L·|det U[K,L]|²,  where a_K = Π_{k∈K} a_k,

which is affine in r_n(U). The products in the conjecture are the values at U = P_σ. In dimension four, MOC for real orthogonal U therefore implies MOC for all U, because r₄(U(4)) ⊆ conv r₄(O(4)) (Zheng's Theorem 1.2, cited). Theorem 1 shows that this convexity step fails for every n ≥ 5. It fails at V and on an open neighbourhood of V. Any n ≥ 5 analogue of the dimension-four proof must therefore handle unitary minor vectors directly, or use properties special to determinant functionals.

**It is not a counterexample to MOC.**

> **Proposition 12 (EXACT; [`structure.md`](structure.md) §6).** There are no a, b ∈ ℂ⁵, α ∈ ℂ and real c₁ ≠ 0, c₀ with c₁·Φ(U) + c₀ = Re(α·det(D_a + U D_b U*)) for all U ∈ U(5).

*Proof.* Suppose the identity holds on U(5). Then it holds on the permutation matrices, and both sides become functions on S₅. Each side has the form f(σ) = Σ_k Σ_{|K|=|L|=k} C^{(k)}_{KL}·1[σ(K) = L].

Let Π be the permutation representation of S₅ on 2-subsets. It decomposes as (5) ⊕ (4,1) ⊕ (3,2), and the (3,2) summand is W = ker(∂ : e_{ij} ↦ e_i + e_j), of dimension 5. Only the levels k = 2 and 3 contribute to the (3,2)-isotypic part of f, and level 3 folds into level 2 through 1[σ(K) = L] = 1[σ(Kᶜ) = Lᶜ]. So the (3,2)-part of f is σ ↦ tr(P_W·C·P_W·Π(σ)), where C is the folded 10×10 coefficient matrix. Matrix coefficients of an irreducible representation are linearly independent, so the operator P_W C P_W, and in particular its rank, is determined by the function f.

- For Φ this rank is 5 (exact sympy computation, [`structure.md`](structure.md) appendix check (2), re-run for this write-up).
- For the determinant functional, C_{KL} = a_{Kᶜ}·b_L + a_K·b_{Lᶜ} has complex rank ≤ 2. So Re(α·C) has real rank ≤ 4, and so does its compression to W.

Hence c₁ = 0. ∎

If such an identity held with c₁ > 0, then Φ(V) = 37/25 > 1 = max_σ Φ(P_σ) would violate MOC at V. Proposition 12 rules this out. It does **not** prove MOC at V: some other determinant functional could still fail there.

**Numerical tests of MOC at V (NUMERIC).** For each spectrum, the test asks whether det(D_a + V D_b V*) lies in conv{Π_i (a_i + b_σ(i))}.
- *This write-up* (script not committed): 300 complex Gaussian spectra, where the largest signed distance outside the hull, relative to the hull's diameter, was −5.1·10⁻² (inside). Also 300 unit-circle and 300 real spectra, where the hull degenerates to a segment, and the point lay on it up to 2·10⁻¹⁶.
- *`ganesh-claude-minors`* (not re-run here): 420 complex and 110 real spectra, some of them optimised. Every one had distance 0 from the hull.

**What it does not exclude.**
- A reduction of MOC to orthogonal matrices that uses properties special to determinant functionals.
- An n = 5 proof that runs a contact reduction directly for each U ∈ U(5).

## 7. Reproduction

Run everything from this folder. The exact checkers need Python 3.9 or later. `check_exact.py` and `bound/check_certificate.py` use the standard library only; `check_exact.py` additionally runs a sympy cross-check if sympy is installed.

```sh
python3 check_exact.py                  # Prop. 3, 3', 7, Lemma 4, Cor. 5 (10 to 20 s); ends "ALL EXACT CHECKS PASSED"
python3 bound/check_certificate.py      # Theorem 6 (about 5 s); also runs under python3 -O
python3 bound/exact_phi_at_V.py         # Prop. 3 again, independently (sympy)
python3 bound/sanity_numeric.py         # numerics: Phi(V), Cor. 5, local maxima over O(5) (numpy, scipy)
python3 bound/generate_certificate.py   # optional: rebuilds bound/certificate.json (cvxpy, Clarabel; about 35 s)
```

The appendix script in [`structure.md`](structure.md) checks Theorem 8, formula (S) and the rank in Proposition 12. It needs sympy and takes under a second. Copy it into a file and run it.

**Expected key output of `bound/check_certificate.py`:**

```
Gram blocks: 346 ; positive definite (exact): 346 ; sizes [2, 5, 6, 13, 26, 42]
multiplier terms: 8310
c = 1137021/1000000 = 1.137021
remainder R: nonzero coefficients 4804 ; sum of |coefficients| = 8.17e-06 ; largest |coefficient| = 2.36e-07
bound B = c + |R|_1 = 1.13702917
B < 37/25 : True
B <= 13/10: True
B <= 11371/10000 : True
```

**Re-run for this write-up on November 25, 2026.**
- Exact checkers: Python 3.14.7 and Python 3.9.6.
- Numerics: numpy 2.5.5, scipy 1.18.1 and sympy 1.14.0.
- Every exact check passed under both Python versions.
- A separate fact-check pass ran every command above as written (the exact checkers also under Python 3.11.7; the generator in a separate copy, which reproduced `certificate.json` byte for byte), ran both appendix scripts, and re-verified with its own code, sharing nothing with the files here: the certificate identity in integer arithmetic, with positive definiteness by fraction-free elimination; the exact facts about V and its transcriptions; Φ on the 120 permutation matrices, Proposition 7, formula (S) and the transposition counts of §4.5; and the rank 5 in Proposition 12. It also tested Lemma 4, Corollary 5, formula (E), Theorem 8 and the Hadamard step of Corollary 10 on random matrices (NUMERIC). Those scripts are not committed.

**sha256 of the checked files.**

| file | sha256 |
|---|---|
| `bound/certificate.json` | `90cdf5de10628e4cf803e969b26b11a6ae1b7448c7c4b406d627f7ef38e12d3c` |
| `bound/check_certificate.py` | `a7e2b86e3e72aa0e87f734a51db15237527ff4e0bcfc8bfd5e1762b329ea5058` |
| `check_exact.py` | `db573bfc84d1a7e0f29f34942ea2e19af5384b1dd75c1aa2463579d7967781fc` |
| `V_exact.json` | `46564e5c78743abc966002d9ca2bc13452bf20a489dd8d7bf07407c81825461f` |

## 8. Limits and open points

1. **The sharp bound is open.** Φ ≤ 1 on O(5), equivalently E_S ≥ 0 on SO(5), is NUMERIC. It is proved only on the 25 families of Theorem 8, which have measure zero. Degree-6 sum-of-squares certificates stop near 1.1365. Possible next routes are a symmetry-reduced degree-8 search, or a global identity extending Theorem 8, for example through Sp(2) → SO(5).
2. **The maximum over U(5)** is 37/25 numerically. Only local maximality is EXACT (Proposition 3′).
3. **The other owner's certificate** for 13/10 was not run or audited. A verifier written from its JSON alone would settle it. Theorem 1 does not depend on it.
4. **Zheng's Theorem 1.2** is cited, not re-proved. It enters only the n = 4 case of Corollary 2 and the unitary version in §4.4.
5. **No third-owner review.** The proof is computer-assisted and has not been reviewed by a third GitHub owner. A reviewer needs to do two things: run `bound/check_certificate.py` and read its roughly 200 lines against the soundness argument in §4.2; and run `check_exact.py` for Proposition 3.
6. **The spot checks for this write-up** are not committed: the certificate audit, the unitary form of Lemma 4, the span dimension, the MOC tests at V, the U(5) maximisations and the D₅ double cosets. Nor are the larger numerical runs of `h5-verify` quoted in §3 and §4.5 (2,000 maximisations per component, 1,400 on U(5), the Hessians at 80 maxima). The exact facts about V are reproducible from the appendix below.
7. **The zero set of E_S** has not been characterised. Whether it is reducible is UNPROVED.

**Where the component checks and the claim differed, and how this was resolved.**
- *The unitary form of Lemma 4.* `h5-verify` reported that "the identity fails on U(5)", while `h5-structure` gave a unitary version. Both are right about different formulas. Re-checked here on 50 Haar unitaries:
  - the formula with det(U) in front fails, with error 11.2 in Φ;
  - the formula Σ_{σ(I)=J} sgn σ Π U_{iσ(i)} = det(U)·|det U[I,J]|² holds, with error 3·10⁻¹⁵, and so does Φ(U) = 1 − 6·conj(det U)·E_S(U), with error 9·10⁻¹⁵. Lemma 4 proves both.
- *Whether the n ≥ 5 extension depends on the other owner's certificate.* `h5-structure` called Corollary 10 conditional on their certificate. `h5-bound` called it EXACT. Corollary 10 needs only some bound below 37/25 on O(5), and Theorem 6 provides one. So it is EXACT without their certificate.
- *The identity and counting claims stated when this check began* are all confirmed: Φ = 1 − 6E_S on SO(5), |S| = 60, and the identity and (0 1) lie in S, so S ≠ A₅. The spans of r₅(U(5)) and r₅(SO(5)) both have dimension 42 (NUMERIC: rank 42 for 400 samples of each and for both together; σ₄₃/σ₄₂ < 10⁻¹³).
- *The SOS values* 1.1365, 1.13702 and 1.13654 come from programs with different margins and relations. They are consistent.
- *The Hessian at the permutation maximisers.* "6 or 7 zero eigenvalues" and "rank 4 or 3" are the same statement on the 10-dimensional SO(5).
- *Nothing in the original claim was found to be wrong.* Additions made here: det V = −1, so V ∉ SU(5), and the bound improves from 13/10 to 1.1371.

## Appendix: exact facts about V (sympy, under a second)

```python
import itertools, sympy as sp
s, t = sp.sqrt(6), sp.I*sp.sqrt(15)
V = sp.Matrix([[4, 4*s, 4*s, 4*s, 4*s], [4*s, -1-t, -6-2*t, -6+2*t, 9+t], [4*s, -6-2*t, -1+t, 9-t, -6+2*t],
               [4*s, -6+2*t, 9-t, -1+t, -6-2*t], [4*s, 9+t, -6+2*t, -6-2*t, -1-t]]) / 20
I5, Z = sp.eye(5), sp.zeros(5, 5)
X = lambda M: M.applyfunc(sp.expand)
print(X(V*V.H - I5) == Z, X(V - V.T) == Z, sp.simplify(V.trace()), X((V-I5)*(V+I5)*(V*V+V+I5)) == Z, sp.simplify(V.det()))
T = [(-1, (v,), (v,)) for v in range(5)]
for v in range(5):
    A, B = tuple(sorted({(v+1) % 5, (v-1) % 5})), tuple(sorted({(v+2) % 5, (v-2) % 5})); T += [(1, A, B), (1, B, A)]
    for d in (1, 2):
        P, Q = tuple(sorted({v, (v+d) % 5})), tuple(sorted({v, (v-d) % 5})); T += [(-1, P, Q), (-1, Q, P)]
perms = list(itertools.permutations(range(5)))
S = [p for p in perms if sum(c for c, I, J in T if tuple(sorted(p[i] for i in I)) == J) == -5]
sgn = lambda p: sp.Matrix(5, 5, lambda i, j: int(p[i] == j)).det()
print(len(S), sp.simplify(sp.expand(sum(sgn(p)*sp.prod([V[i, p[i]] for i in range(5)]) for p in S))))
sq = lambda m: sp.expand(m*sp.conjugate(m))
pairs = list(itertools.combinations(range(5), 2))
print({sq(V[i, i]) for i in range(5)}, {sq(V[i, j]) for i in range(5) for j in range(5) if i != j},
      sorted({sq(V.extract(list(I), list(J)).det()) for I in pairs for J in pairs}))
```

Expected output: `True True 0 True -1`, then `60 2/25`, then `{1/25} {6/25} [3/125, 8/125, 12/125, 18/125, 27/125]`. That is: V is unitary and symmetric, tr V = 0, (V − I)(V + I)(V² + V + I) = 0, det V = −1, |S| = 60, E_S(V) = 2/25, the squared moduli of the entries are 1/25 on the diagonal and 6/25 off it, and the squared 2×2 minors take exactly the five listed values.
