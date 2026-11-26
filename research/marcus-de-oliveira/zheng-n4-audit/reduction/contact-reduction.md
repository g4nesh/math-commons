Evidence level: EXACT for every finite identity listed in Part 3 (own sympy code, integer/rational/Gaussian-rational arithmetic); PROOF-READ for the logic of Prop. 3.3 (direction used), Lemmas 3.4, 3.5, Theorem 3.6, Lemma 4.2, (4.7), Prop. 4.4 and Section 6; CITED for the classical facts named in the table (affine GIT quotient surjectivity, Kempe generation or the null-cone substitute, removable singularities and the maximum principle for subharmonic functions, the affine dimension theorem, normalisation and completion of curves, Birkhoff–von Neumann); out of scope here and taken as given: Corollary 2.2, Prop. 4.1, the orbit classification (214 / 124 / 90, Appendix E), Prop. 4.3.

# Audit of the contact reduction (paper Sections 3, 4.2–4.3 logic, 6)

Auditor label `sec3`. Paper: D. Zheng, "The Marcus–de Oliveira conjecture for n = 4", 24 Oct 2026. Page numbers below are the printed page numbers (same as the PDF page index).

## Deviations and limits

- The optional numeric experiment on the mechanism of Theorem 3.6 was **not done**. Everything in Part 3 is exact; nothing in this report is only numeric.
- The surjectivity of the affine GIT quotient and Kempe's degree-one generation theorem were **not re-proved**. They are recorded as CITED, with the precise statement needed. Part 2a says which pieces I could prove directly.
- The classical potential-theory facts used in Lemma 3.5 are CITED, not re-derived. These are the removable-singularity theorem for subharmonic functions bounded above near a point, and the maximum principle on a compact connected Riemann surface.
- I did not check Prop. 4.1 (facet completeness), the orbit count 214 = 124 + 90, membership of Appendix E, or any of the 90 cases. The reduction uses these as black boxes; this report checks only that the logic connecting them is sound and that the black boxes are stated in the form the logic needs.
- I did not use the shared library. My P was built independently by Plücker straightening and checked as a polynomial identity. It is in `check_section3_P.json` for comparison with the library.

## Verdict

The implication

  [Cor. 2.2] + [Prop. 4.3 for T_1..T_90] + [124 minor-only orbits, with Prop. 4.1 and Appendix E] ⟹ MOC(4)

is **correct modulo the listed classical citations**. I found no logical gap.

Weak points, none of which breaks the proof:

1. **Closedness (p. 16, the paragraph after (3.1)).** The claim that every nonzero x ∈ X is α·(matching vector of two finite spectra) is load-bearing. It is needed for the point produced by Theorem 3.6 and in the proof of Prop. 4.3. It is not proved in the paper: it is the classical fact that the image of the matching map is closed (GIT surjectivity plus Kempe, or a null-cone substitute). Mathematically it is sound and standard, but it is a citation, not a proof.
2. **(3.9) is stated as an equality X = {∇C = 0}, citing HMSV.** Only the inclusion X ⊆ {∇C = 0} is used in Section 3. I proved that inclusion EXACTLY by symbolic substitution. The equality is unnecessary for the reduction.
3. **Misprint on p. 25 (4.7).** For g·σ = ρσρ′ the printed g·O = R_{ρ′}^{-1} O R_ρ^{-1}, i.e. (g·O)_{ij} = O_{ρ′^{-1}(i), ρ^{-1}(j)}. That is covariant for σ ↦ ρσρ′^{-1}, not for ρσρ′. A correct choice is (g·O)_{ij} = O_{ρ′(i), ρ^{-1}(j)}. Because ρ′ ranges over all of S_4 this is harmless. Covariance itself was checked EXACTLY on generators.
4. **(3.8) holds only for det O = +1.** For det O = −1 the sign flips (checked EXACTLY). The paper uses it only on SO(4), so this is fine. The det −1 images in Prop. 4.4 are correctly repaired by the column flip on p. 25.

## Part 1. Self-contained derivation

Notation: v = (v_1..v_8) ∈ (C²)⁸, [ij] = det(v_i, v_j), and B_0..B_13 are the noncrossing matchings (A.1). Ψ(v) = (B_0(v), .., B_13(v)) ∈ C¹⁴ = V_C, with V_R = R¹⁴. S_4 is listed lexicographically, and p_σ = ∏_i [i, 4+σ(i)]. For x ∈ V_C, F_σ(x) = (Px)_σ, where P is the integer 24×14 matrix with p_σ = Σ_a P_{σa} B_a.

**D1 (matching basis, P).** Every perfect-matching product is an integer combination of the B_a, via repeated Plücker uncrossing [ac][bd] = [ab][cd] + [ad][bc] for a<b<c<d. The resulting P satisfies p_σ = Σ P_{σa} B_a identically in 16 homogeneous variables. Its rank is 14, and the printed 14×14 submatrix has det 1. Hence P is injective, and for x ∈ V_C the 24 numbers F_σ(x) determine x. *EXACT.*

**D2 (κ).** κ = e_13 gives Pκ = 1, so κ is the unique such vector by injectivity. κ = Ψ(configuration with v_1..v_4 = (0,1), v_5..v_8 = (1,1)), so κ ∈ X. *EXACT.*

**D3 (X).** X is the Zariski closure of Ψ((C²)⁸).

- X is irreducible, because it is the closure of the image of an irreducible variety. No GIT is needed.
- X is a C*-cone: scaling column 1 by θ multiplies Ψ by θ.
- dim X = 6. The differential dΨ has rank exactly 6 at a rational configuration, so the generic rank is ≥ 6 and dim X ≥ 6. Ten polynomial vector fields lie in ker dΨ everywhere: three from sl₂, and seven from the torus T₀ = {t ∈ (C*)⁸ : ∏ t_i = 1}, which fixes Ψ because every matching is multilinear of degree one in each column. These ten fields are independent at the same point, hence on a dense open set, so the generic rank is ≤ 6. In characteristic 0 the generic rank equals the dimension of the image closure. *EXACT (ranks) + standard fact.*
- Only dim X ≥ 6 and irreducibility are used later.

**D4 (finite spectral representation; see 2a).** Every x ∈ X∖{0} equals α·Ψ(v(a,b)), with v(a,b) = ((a_i,1)_{i≤4}, (−b_j,1)_{j≤4}), finite a, b ∈ C⁴ and α ≠ 0. Then F_σ(x) = α∏_i(a_i + b_{σ(i)}), which is (3.1) (checked EXACTLY on random data). *CITED: Ψ((C²)⁸) = X; the rest is PROOF-READ.*

**D5 (gaps).** Fix π*. Define g_σ(x) = F_σ(x) − F_{π*}(x) and c(x) = Re x − Re F_{π*}(x)·κ ∈ V_R. P is real, so (Pc(x))_σ = Re g_σ(x).

- Changing π* to another contact π′ (a permutation with gap 0) changes c by Re(F_{π*} − F_{π′})κ = 0.
- The common kernel of all g_σ is Cκ, by injectivity of P.

*EXACT.*

**D6 (Leibniz sum, (3.5)).** Φ_c(O) = Σ_σ sgn σ (Pc)_σ ∏_i O_{iσ(i)} for O ∈ SO(4). For x = αΨ(v(a,b)), since O is real:

  Φ_{c(x)}(O) = Re(α Σ_σ sgn σ ∏ O_{iσ(i)} ∏(a_i + b_{σ(i)})) − Re F_{π*}(x)·det O.

The inner sum is det(D_a O + O D_b) = det(D_a + O D_b Oᵗ)·det O. With det O = 1 this gives (3.5). *EXACT (3 random Gaussian-rational instances × all 24 references, rational Cayley O) and PROOF-READ.*

**D7 (Prop. 3.3, direction ⇐, the only direction used).** Assume (3.6) for all x ∈ X, all π*, and all O ∈ SO(4). Fix spectra a, b and O ∈ SO(4), and put x₀ = Ψ(v(a,b)) and D = det(D_a + O D_b Oᵗ).

- For θ ∈ C∖{0}, θx₀ ∈ X (D3) and F_σ(θx₀) = θz_σ.
- Choose π* minimising Re θz_σ. Then Pc(θx₀) ≥ 0, and (3.5) gives Re θD ≥ min_σ Re θz_σ.
- If D ∉ K = conv{z_σ}, which is compact and convex in C = R², strict separation gives θ with Re θD < min_K Re θz. That is a contradiction. Hence D ∈ K, which is MOC for A = D_a, B = O D_b Oᵗ, O ∈ SO(4).
- If O ∈ O(4) has det −1, replace O by O·diag(−1,1,1,1), which leaves O D_b Oᵗ unchanged.
- Cor. 2.2 then gives the unitary case.

The quantifier "for all θ ≠ 0" is exactly the support-function characterisation of a compact convex set, so nothing is lost. *PROOF-READ.*

Direction ⇒ uses D4 and is not needed for Section 6.

**D8 (entry subspace, (3.7), (3.8)).** For Q ∈ C^{4×4}, E(Q) is the solution of PE(Q) = (Σ_i q_{iσ(i)})_σ.

- It exists. The system is consistent on all 24 rows, and the solution equals the coordinates printed in App. A.2.
- dim E = 10, and E = {m₀ = m₁ = m₅ = m₆ = 0}.
- For real Q and O ∈ SO(4): Φ_{E(Q)}(O) = Σ_{i,j} q_{ij} O_{ij} cof_{ij}(O) = Σ q_{ij} O_{ij}², since cof(O) = det(O) O^{-ᵀ} = O. With det O = −1 the sign flips.
- E ∩ V_R = E(R^{4×4}), because E(Re Q) = Re E(Q).
- Birkhoff: (O_{ij}²) is doubly stochastic, so it equals Σ λ_σ P_σ with λ ≥ 0 and Σλ = 1. Hence Φ_{E(Q)}(O) = Σ λ_σ (PE(Q))_σ ≥ 0 whenever PE(Q) ≥ 0.

*EXACT (symbolic Q, rational O) + Birkhoff (classical).*

**D9 (cubic, Lemma 3.4).**

- C (A.2) is transcribed from the p. 36 image. ∇C∘Ψ ≡ 0 as polynomials in the 16 homogeneous variables, so X ⊆ {∇C = 0}.
- ∇C(κ) = 0 and H_κ = Hess C(κ) = e₀e₆ᵗ + e₆e₀ᵗ − e₁e₅ᵗ − e₅e₁ᵗ, of rank 4, with ker H_κ = E.
- ∇C(x+tκ) = ∇C(x) + tH_κx exactly. The t² term is ∇C(κ) = 0.

All *EXACT*. Lemma 3.4 follows:

1. If x + Cκ ⊆ X, then ∇C vanishes at t = 0 and t = 1, so H_κx = 0.
2. H_κ is real and H_κκ = 0, so H_κ c(x) = 0. Hence c(x) ∈ E ∩ V_R.
3. If Pc(x) ≥ 0, D8 gives Φ_{c(x)}(O) ≥ 0 on SO(4).

*PROOF-READ, uses only X ⊆ Sing C.*

**D10 (Lemma 3.5).** The argument is checked in detail in 2c. *PROOF-READ + CITED classical potential theory.*

**D11 (Theorem 3.6).** The argument is checked in detail in 2d. *PROOF-READ + CITED (affine dimension theorem, curve selection, normalisation).*

**D12 (Lemma 4.2, logic).** Identify c ∈ V_R with the affine function v ↦ v(c) on aff(Q), which has dim 13. This is a bijection onto all affine functions, since P is injective, dim V_R = 14, and κ maps to the constant 1.

- If Q = {f_j ≥ 0, j = 0..81} (Prop. 4.1), the affine Farkas lemma gives c = λ₀κ + Σ λ_j f_j with all λ ≥ 0 for every c ≥ 0 on Q.
- Evaluate at the barycentre v_T of the T-vertices: 0 = λ₀ + Σ λ_j f_j(v_T), with every term ≥ 0. So λ₀ = 0, and λ_j = 0 unless f_j vanishes at all T-vertices (each f_j is ≥ 0 at each vertex).
- The reverse inclusion is trivial.

*PROOF-READ; rests on Prop. 4.1 (not checked here).*

**D13 (symmetries, (4.7)).** Relabelling the 8 columns by τ ∈ S₈ commutes with SL₂ × T₀ and acts on V_C by an integer matrix M_τ with M_τΨ(v) = Ψ(v∘τ), so M_τX = X. Checked EXACTLY for the generators: a transposition and a 4-cycle in each block, and the block swap i ↔ 4+i.

- Each M_τ fixes κ and permutes the rows of P. The block swap sends F_σ to F_{σ⁻¹} with no sign, because each bipartite bracket flips sign and there are four of them.
- So M_τ c_{π*}(x) = c_{τ·π*}(M_τx). Gap vectors of points of X go to gap vectors of points of X, and c ∈ N_T ⟺ M_τc ∈ N_{τ·T}.
- With the extended Leibniz sum Φ_c(O) = Σ (Pc)_σ sgn σ ∏O_{iσ(i)}/det O on O(4), Φ_{Mc}(O′) = Φ_c(O) for an explicit O′ (a row or column permutation of O, or Oᵗ). Flipping one column of O′ leaves the extended sum unchanged and lands in SO(4).

*EXACT on generators + PROOF-READ*; the general element is a composite.

**D14 (Prop. 4.4).** Suppose x ∈ X, Pc(x) ≥ 0, and at least 6 zeros. Pick T ⊆ zeros with |T| = 6, so c(x) ∈ N_T. Pick g ∈ G with gT = T_rep, so M_g c(x) = c_{gπ*}(M_g x) ∈ N_{T_rep} with M_g x ∈ X.

- If T_rep is minor-only, Lemma 4.2 writes M_g c(x) as a nonnegative combination of minor vectors, and (4.2) gives Φ ≥ 0 on SO(4).
- Otherwise T_rep = T_ν and Prop. 4.3 applies. Prop. 4.3 has a non-strict hypothesis (c ∈ N_{T_ν}, contact set ⊇ T_ν), which is what is needed.
- The reference gπ* need not lie in T_ν. This is harmless: by D5 the gap vector equals c_{π′} for any π′ ∈ T_ν.
- Return by D13: Φ_{c(x)}(O) = Φ_{M_g c(x)}(O′) ≥ 0 with O′ ∈ SO(4) after a column flip.

*PROOF-READ*; rests on Prop. 4.1, the orbit classification and Appendix E (not checked here). (4.2) itself was checked EXACTLY for all 52 (I,J) at two rational O.

**D15 (Section 6).** Suppose MOC(4) fails. By D7 some x ∈ X, π* and O₀ ∈ SO(4) have Pc(x) ≥ 0 and Φ_{c(x)}(O₀) < 0. Theorem 3.6 gives x′ ∈ X with the same O₀ and π*, at least 6 contacts, and Φ_{c(x′)}(O₀) < 0. D14 gives Φ_{c(x′)}(O₀) ≥ 0, a contradiction. *PROOF-READ.*

## Part 2. Attack log

### 2a. X: irreducibility, dimension, finite spectra

- **Irreducible and C*-stable:** elementary (D3).
- **dim X = 6:** EXACT (D3); only ≥ 6 is used.
- **Nonzero ⇒ finite spectra (D4).**
  - (i) *Unstable configurations give 0.* Suppose some column is 0, or some point of P¹ carries ≥ 5 of the 8 columns. Then every perfect matching pairs two of those columns, since the cluster has ≥ 5 members and only ≤ 3 partners outside it, so every element of R₁ vanishes.
  - (ii) *Semistable configurations give nonzero vectors.* If no point carries more than 4 columns, list the columns grouped by position and pair entry k with entry k+4. Each class has size ≤ 4, so every pair joins distinct points, that matching product is nonzero, and Ψ(v) ≠ 0.
  - (iii) *Chart change.* If x = Ψ(v) ≠ 0, choose a point of P¹ off the ≤ 8 positions and an SL₂ element sending it to ∞. Rescaling each column to second coordinate 1 multiplies Ψ by a nonzero scalar. This yields α ≠ 0 and finite a, b. Repeated points (repeated eigenvalues, a_i = −b_j) are allowed, and strictly semistable configurations (exactly 4 coincident) are covered.
  - (i)–(iii) are elementary and PROOF-READ.
- **What is not elementary is that every point of X is Ψ(v) for some v, i.e. Ψ((C²)⁸) is closed.** With G = SL₂ × T₀ (reductive), R := C[(C²)⁸]^G is the diagonal ring ⊕_d R_d. Then:
  - (G1) the quotient map (C²)⁸ → Spec R is surjective (Hilbert; Mumford–Fogarty–Kirwan, GIT, Thm 1.1: CITED, source not re-read);
  - (G2) Spec R → X is surjective. This follows either from Kempe's theorem (R = C[R₁], so Spec R ≅ X), as the paper cites [15, §2.1], or from finiteness of R over C[R₁]. Finiteness holds because by (i) and the first fundamental theorem for SL₂, every v killing R₁ kills all of R₊. (In a bracket monomial of multidegree (d..d), a 5-cluster needs 5d nonvanishing brackets, but there are only 4d brackets.) This uses the FFT, which is also classical.
- **No weaker fact suffices.** Theorem 3.6 outputs a point of the closed curve C ⊂ X where the hypothesis of Lemma 3.5 fails. That condition is a system of equalities and a sign condition, so it cannot be perturbed into Ψ((C²)⁸). Prop. 4.3 is then applied to that point.
- The paper's paragraph on p. 16 is a correct sketch that presupposes (G1) and (G2).
- **Status: CITED (classical), not a gap.**

### 2b. Lemma 3.4: inclusion versus equality; (3.10)–(3.12); Birkhoff

- Only X ⊆ {∇C = 0} is used. Lemma 3.4 needs ∇C = 0 at x, at x + κ and at κ.
- The inclusion is EXACT: the 14 partials composed with Ψ are the zero polynomial in 16 homogeneous variables. Hence they vanish on X, the closure of the image.
- The equality X = Sing(C) (HMSV) is not needed anywhere in Section 3. Nor is it needed in the Group IV/V sketches on pp. 26–27, which use ∇C(x) = 0 and C(Re x) = 0 for x ∈ X, both consequences of the inclusion.
- (3.10), rank 4, (3.11) and (3.12) are EXACT. The kernel equals E, with dim E = 10 and the App. A.2 coordinates reproduced.
- (3.8) needs det O = 1. EXACT: the sign is −1 for det O = −1. Lemma 3.4 is stated on SO(4), so this is consistent.
- Birkhoff–von Neumann is applied to the doubly stochastic (O_{ij}²) with weights on the row sums Σ_i q_{iσ(i)} = (Pc)_σ ≥ 0. Correct.

### 2c. Lemma 3.5 (maximum principle)

- **Boundary points in Σ°.** Let q ∈ cl(Ω)∖Ω. By continuity all h_j(q) ≥ 0 and L(q) ≤ 0.
  - If all h_j(q) > 0 and L(q) < 0, then q ∈ U (the open set {all h_j > 0, L < 0}). A small disc around q lies in U and meets Ω, so q is in the same component. That contradicts q ∉ Ω.
  - So either some h_j(q) = 0, and the hypothesis gives L(q) ≥ 0, or L(q) = 0. Either way L(q) = 0.
  - A point of U in a different component of U is not in cl(Ω), so the gluing is unaffected.
- **Gluing.** v = −L on Ω and 0 elsewhere is continuous on Σ°.
  - On Ω, v is harmonic, so it has the mean-value property on small circles.
  - Off Ω, v(q) = 0 ≤ the mean of v ≥ 0.
  - A continuous function with the local sub-mean-value property is subharmonic. Correct.
- **Removable singularities.** v is bounded above near each p ∈ S. It equals −L on Ω near p, which is bounded by hypothesis, and 0 elsewhere. A subharmonic function bounded above on a punctured disc extends subharmonically, since points are polar (classical; CITED).
- **Compact surface.** An upper-semicontinuous function attains its maximum on compact Σ. By the strong maximum principle on connected Σ, v is constant.
  - If Ω ≠ ∅ the constant is positive, so v > 0 on Σ°, forcing Ω = Σ°. Σ° is connected because Σ is connected and S is finite.
  - Then all Re H_j > 0 on Σ°. A pole of order m at p would give Re(c z^{−m}(1 + O(z))) taking both signs on small circles. So every H_j is holomorphic on compact Σ, hence constant. This contradicts "not all constant".
  - **The case Ω = Σ° is covered by this step.** If S = ∅, all H_j are constant from the start.
- **Singular points of the curve.** After normalisation, Σ° is a smooth affine curve mapping onto C, and the functions are pullbacks. The hypothesis on Σ° is equivalent to the hypothesis on C, and the point where it fails maps to a point of C ⊂ X. Branches at singular points are separate points of Σ°, which is harmless.
- **"Not all constant" (checked in 2d).**
- **Status: PROOF-READ, correct.**

### 2d. Theorem 3.6

- **Normalisation.** If Σ_σ Re g_σ = 0 with all gaps ≥ 0, then Pc = 0, so c = 0 and Φ = 0, a contradiction. Scaling by a real r > 0 keeps x in X and scales c and Φ by r.
- **Bookkeeping.**
  - The contact set is by definition the exact zero set of the gaps. π* is always in it, so k ≥ 1, and every non-contact gap is > 0.
  - If k ≥ 6, stop.
  - A non-unique reference is irrelevant: π* stays fixed, g_{π*} ≡ 0, and by D5 the gap vector does not depend on which contact is the reference.
- **Dimension.** The conditions are g_σ(y) = g_σ(x) for the k−1 contacts σ ≠ π*, plus Σ_σ g_σ(y) = Σ_σ g_σ(x). These are k complex affine-linear equations, possibly dependent. X is irreducible with dim ≥ 6 and contains x. By the affine dimension theorem (dim of each component of X ∩ L ≥ dim X + dim L − 14, with codim L ≤ k), every component Y through x has dim ≥ 6 − k ≥ 1. Only dim X ≥ 6 is needed. CITED, standard.
- **Exceptional case.** If every g_σ is constant on Y, then P(y − x) ∈ C·1, so y − x ∈ Cκ. Hence Y ⊆ x + Cκ, and since Y is positive-dimensional and irreducible, Y = x + Cκ ⊆ X. Lemma 3.4 gives Φ_{c(x)}(O₀) ≥ 0, a contradiction. Correct.
- **Curve selection.** Let f = g_σ be nonconstant on Y.
  - The paper's general-hyperplane argument is sound. Each general hyperplane through x contains no positive-dimensional component through x of Y, nor of the level set Z = Y ∩ {f = f(x)}, which has dim Y − 1. After dim Y − 1 cuts, a curve component C through x would lie in Z, which has dimension 0 at that stage, if f were constant on C.
  - A simpler citation also works: any two points of an irreducible variety lie on an irreducible curve. Take y ∈ Y with f(y) ≠ f(x).
- **Old contacts keep gap 0.** g_σ is fixed on C to g_σ(x), which has real part 0.
- **Applying Lemma 3.5.**
  - The H_j are the g_σ with σ outside the current contact set. The nonconstant f is among them, because the old-contact g's and g_{π*} are constant. Its pullback to the normalisation is nonconstant.
  - Σ_j Re H_j ≡ 1 on C, because the complex sum is fixed and the old real parts are 0.
  - L(y) = Φ_{c(y)}(O₀) = Σ_σ sgn σ ∏_i O₀[i,σ(i)]·Re g_σ(y) is a real combination of the h_j with |λ_σ| ≤ 1. On Ω each h_j ∈ (0,1), so |L| ≤ 1. The boundedness hypothesis holds globally.
  - At x all h_j > 0 and L < 0, so Ω ∋ x is nonempty.
  - Lemma 3.5 therefore forces its hypothesis to fail. Some point x′ ∈ C has all h_j ≥ 0, some h_j = 0, and L(x′) < 0.
- **Properties of x′.**
  - It lies in C ⊂ X.
  - All gaps are ≥ 0: old ones are 0, the others are ≥ 0.
  - The contact set strictly grows.
  - Φ_{c(x′)}(O₀) < 0 at the same O₀.
  - Σ Re g_σ(x′) = 1, so x′ ∉ Cκ and in particular x′ ≠ 0.
  - It need not lie in the component Ω, which is irrelevant.
  - Its finite spectral representation is D4/2a (CITED).
- **Induction.** k strictly increases from k ≥ 1, so at most 5 rounds reach k ≥ 6. Sound.
- **Status: PROOF-READ, correct.**

### 2e. Prop. 3.3 and the fourth root

- The "for all θ ≠ 0" family of support inequalities is equivalent to membership in the compact convex hull (D7). Correct.
- Fourth root in Prop. 4.3: every bracket of v(ta, tb) is t times the bracket of v(a, b) (bipartite and non-bipartite alike), so Ψ(v(ta,tb)) = t⁴Ψ(v(a,b)) = αΨ(v(a,b)) = x. EXACT as a polynomial identity in t. Legitimate. Injectivity of P is not even needed.
- The footnote's x ↦ ix convention is correct: Im(iz) = Re z, and X is a C*-cone.

### 2f. Prop. 4.4 and symmetry

- **G preserves X.** Each generator is a column relabelling. M_τΨ(v) = Ψ(v∘τ) is EXACT.
- **G preserves κ and permutes the F_σ.** EXACT: the block swap gives σ ↦ σ⁻¹ with sign +1.
- **Gap vectors go to gap vectors.** M_τ c_{π*}(x) = c_{τπ*}(M_τx).
- **Leibniz covariance holds on O(4) with the det-normalised extension.** EXACT on generators. The det = −1 repair by a column flip is EXACT.
- **Misprint (minor) on p. 25.** The displayed g·O does not match σ ↦ ρσρ′ (it matches ρσρ′^{-1}). A correct version is (g·O)_{ij} = O_{ρ′(i), ρ^{-1}(j)}.
- **The non-strict form of Prop. 4.3 is exactly what is used.** Its hypothesis is c(x) ∈ N_{T_ν}: zeros on T_ν, ≥ 0 elsewhere, larger contact sets allowed (p. 21). For the 124 orbits, Lemma 4.2 likewise needs only c ∈ N_T.
- **The reference permutation is free (see D14).**
- **Remaining dependencies, not checked here:**
  - Prop. 4.1, needed for Lemma 4.2;
  - that the 124 minor-only representatives and T_1..T_90 together meet every one of the 214 G-orbits of 6-subsets;
  - the F_T computations.

### 2g. Other findings

- No hidden genericity assumption in Section 3. Degenerate spectra (coincident eigenvalues, a_i + b_j = 0, unstable limits giving 0, strictly semistable limits) are all handled by 2a. The burden passes to Prop. 4.3, which must hold for **every** finite spectral pair. The paper says so explicitly on p. 21 ("eigenvalues may coincide").
- No strictness is used on a boundary. The only strict inequalities are h_j > 0 and L < 0 defining Ω, and these are handled correctly at ∂Ω (2c).
- No quantifier swap. O₀ is fixed through Theorem 3.6, and Prop. 4.4 is applied for that same O₀ (indeed for all O). Theorem 3.6 needs (3.6) at no point; it is a pure existence statement.
- Section 6 uses only the ⇐ direction of Prop. 3.3, which does not need D4 for its input point. That point is Ψ(v(a,b)) scaled by θ, so it is in the image by construction. D4 is needed only for the output of Theorem 3.6.
- Not in the reduction's path, informational: the 14 partials of C are linearly independent (det 1 on the listed monomials), and C(T_i m) = −C(m) for the 7 adjacent transpositions, with T_i² = I. Both EXACT.

## Part 3. Exact replays (check_section3.py, all PASS; output in check_section3_output.txt)

1. (A.1) is exactly the set of 14 noncrossing perfect matchings of 1..8 (out of 105).
2. P was built by Plücker straightening. p_σ = Σ P_{σa}B_a holds as a polynomial identity in 16 homogeneous variables. P is an integer matrix of rank 14, and the App. A 14×14 submatrix has det = 1. P is dumped to check_section3_P.json, with rows in lexicographic one-line order.
3. Pκ = 1 with κ = e₁₃, and κ = Ψ(4×(0,1), 4×(1,1)).
4. ∇C∘Ψ ≡ 0 and C∘Ψ ≡ 0 symbolically.
5. Hess C(κ) as in (3.10), rank 4, kernel {m₀=m₁=m₅=m₆=0}; ∇C(κ) = 0; polarisation (3.12) exactly.
6. E(Q) solves (3.7) on all 24 rows, equals the App. A.2 coordinates, and dim E = 10.
7. (3.8) holds with symbolic Q at 3 rational O ∈ SO(4); the sign flips for det −1.
8. All 52 indicator vectors 1{σ(I)=J} lie in P(V_R), giving 34 distinct minor vectors. (4.2) holds for all of them at 2 rational O.
9. (3.1) and (3.5), with (Pc)_σ = Re g_σ, on 3 random Gaussian-rational (a, b, α) × 24 references. The fourth-root identity holds.
10. rank dΨ = 6 at a rational point, and 10 independent group directions lie in the kernel, so dim X = 6.
11. Symmetry generators: integer matrices fixing κ; M_τΨ = Ψ∘τ; row permutation of P (σ ↦ σ⁻¹ for the swap); (4.7) on O(4); det −1 column flip.

## Step / status / what it rests on

| Step | Status | Rests on |
|---|---|---|
| P, rank 14, det-1 minor, κ | EXACT | own straightening + polynomial identity |
| X irreducible, C*-cone | PROOF-READ | closure of image of (C²)⁸ |
| dim X = 6 (only ≥ 6 used) | EXACT + standard | Jacobian rank; generic rank = dim (char 0) |
| nonzero x ∈ X ⇒ αΨ(v(a,b)), finite a, b | CITED | GIT surjectivity of affine quotient (MFK Thm 1.1) + Kempe [15 §2.1] or FFT null-cone finiteness; chart change PROOF-READ |
| (3.1), (3.5), gaps (Pc)_σ = Re g_σ | EXACT | det O = 1 |
| Prop. 3.3 ⇐ | PROOF-READ | separation in R²; Cor. 2.2 (given); O(4)→SO(4) by a column sign |
| Prop. 3.3 ⇒ (not used in Sec. 6) | PROOF-READ | needs D4 (CITED) |
| X ⊆ Sing C | EXACT | symbolic substitution |
| X = Sing C (3.9) | CITED, **not needed** | HMSV Prop 2.2(b), Cor 4.2 |
| (3.10), (3.11), (3.12), E, (3.7) | EXACT | own code |
| (3.8) | EXACT | det O = +1 essential |
| Lemma 3.4 | PROOF-READ | inclusion only + Birkhoff–von Neumann (classical) |
| Lemma 3.5 | PROOF-READ | removable singularity for bounded-above subharmonic functions; maximum principle (classical, CITED) |
| curve normalisation/completion | CITED (standard) | normalisation of affine curves; compact Riemann surface |
| Theorem 3.6 dimension count | PROOF-READ | affine dimension theorem (CITED, standard) |
| Theorem 3.6 curve selection | PROOF-READ | general hyperplanes, or curve through two points (standard) |
| Theorem 3.6 induction, new point | PROOF-READ | D4 for finite spectra of x′ |
| Lemma 4.2 | PROOF-READ | affine Farkas; **Prop. 4.1 (not checked here)** |
| (4.2) minor vectors | EXACT | own code, 2 rational O |
| symmetry action, (4.7), det −1 fix | EXACT (generators) + PROOF-READ | p. 25 g·O formula misprinted (harmless) |
| Prop. 4.4 | PROOF-READ | orbit classification 214/124/90 and Appendix E (**not checked here**); Prop. 4.3 (given, non-strict form) |
| Section 6 | PROOF-READ | all of the above |
