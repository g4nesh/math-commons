Evidence level: EXACT for every algebraic identity in Appendix B; PROOF-READ for the existence steps (Birkhoff, Carathéodory, max-flow min-cut, divisibility by an irreducible polynomial); CITED (classical) only for those three named theorems; NUMERIC for adversarial searches.

# Referee report: Appendix B (Theorem B.1, Lemma B.2, Corollary B.3)

Paper: Zheng, "The Marcus–de Oliveira conjecture for n = 4", pages 37–41.
Parallel supplement notes compared: `final_referee_shared_zero_birkhoff_20261022` (3 pp., = Theorem B.1) and
`final_referee_single_collision_divisor_compensation_20261023` (3 pp., Lemma 1 = Corollary B.3, Theorem 2 = Lemma B.2).
Own code: `appB_checks.py` (runs in about 45 s, deterministic seeds, nothing from the supplement is executed or imported); output in `appB_checks_output.txt`.

## Deviations and limits

- No deviation from the assigned scope.
- Symbolic unitary matrices are not used directly. Instead: (i) every step that is a polynomial identity *before* unitarity is imposed was checked symbolically with generic entries (entries and conjugates as independent symbols); (ii) every step that uses unitarity was replayed exactly on Gaussian-rational unitary matrices (Cayley transforms) plus structured degenerate ones; (iii) the remaining steps (the Jacobi complementary-minor identity, Cauchy–Binet for a Gram determinant) are classical and were checked on samples only. Identity (B.1) for a given U is proved exactly, not sampled: both sides are polynomials of degree at most 1 in each of a1..a3, e1..e3, so equality on the grid {0,1}^6 proves it. The same holds for Lemma B.2 on a 320-point grid (degree ≤ 1 in a0..a3, α, β and ≤ 4 in γ).
- Birkhoff–von Neumann, Carathéodory and max-flow min-cut are classical external theorems. I checked their hypotheses. I also replaced each one by an explicit exact computation on every sample: a greedy exact Birkhoff decomposition, exact Carathéodory reduction, and an exact Edmonds–Karp flow.
- Downstream uses: I read Section 1.2 (Prop. 1.4), Sections 3.4, 4.3 and 5.5, Appendix F of the paper, and the collision passages of five supplement case notes (inventory 40, 54, 55, 60, 85). I did **not** audit any case-specific claim that "exactly one denominator factor vanishes" at a given puncture. That is for the case auditors.
- I did not audit the supplement's extra Remark 1 and Theorem 2 (failure of a uniform law for n ≥ 5). They are not part of Appendix B and are not used for n = 4.

## 1. Statements restated

Main-text convention: normal A and B have eigenvalues a_i and b_j, and det(A+B) = det(D_a + U D_b U*) with U unitary. Appendix B writes e = −b, indexes everything by 0..3 (main-text index = App. B index + 1), and permutations act on {0,1,2,3}.

**Theorem B.1.** Fix U ∈ U(4). There is a probability vector (μ_σ)_{σ∈S4} that depends on U alone (and on the convention that the coincident pair has index 0 in both spectra), such that for all d, a1, a2, a3, e1, e2, e3 ∈ ℂ

  det(diag(d,a1,a2,a3) − U diag(d,e1,e2,e3) U*) = Σ_σ μ_σ ∏_{i=0}^{3} (a_i − e_{σ(i)}), with a0 = e0 = d,

with at most 10 nonzero μ_σ and Σ_{σ(0)=0} μ_σ = |u00|². In main-text notation: if a_i + b_j = 0 for some i, j, then for every U there is one law μ (after relabelling i, j to 0) with det(D_a + U D_b U*) = Σ μ_σ ∏(a_k + b_σ(k)) for every spectrum with a_i + b_j = 0. In particular MOC(4) holds at such spectra.

**Lemma B.2.** Fix U ∈ U(4) and set e = (α, β, γ, γ), so that b_2 = b_3 in main-text notation. There is a probability law μ on S4, depending only on U (in fact only on the two columns p, q of U that carry α and β), with det(diag(a) − U diag(e) U*) = Σ μ_σ ∏(a_i − e_σ(i)) for all a0..a3, α, β, γ ∈ ℂ. By relabelling the same holds for any designated pair e_i = e_j, and, via U ↦ U* and σ ↦ σ⁻¹, for any designated pair a_i = a_j.

**Corollary B.3.** Fix O ∈ SO(4) and one bracket d among the 28 brackets [A_i,C_j], [A_i,A_j], [C_i,C_j] of homogeneous spinors A_k, C_k ∈ ℂ². Put Π_σ = ∏_k [A_k, C_σ(k)] and Δ_O = Σ_σ sgn σ Π_σ ∏_k O_{k,σ(k)}.
(a) There is a law μ that depends on O and on the chosen bracket, and a polynomial Q_O, with Δ_O − Σ μ_σ Π_σ = d · Q_O. This is (B.9).
(b) Let A(t), C(t) be holomorphic near t = 0, and let E(t) = d(A(t),C(t)) R(t) with R holomorphic, R(0) ≠ 0 and E ≢ 0. Fix η ∈ ℂ and τ ∈ S4. Put g_σ = η(Π_σ − Π_τ)/E and L = Re(η(Δ_O − Π_τ)/E). Then L is bounded below near t = 0 on {all Re g_σ ≥ 0}. The bound depends on O, η and the curve germ, and is not uniform in any of them. Nothing is asserted when a second factor of E vanishes at t = 0.

## 2. Step table

| # | Step (page/eq.) | Check | Status |
|---|---|---|---|
| B1.1 | d-shift reduction to a0 = e0 = 0 (p.37) | trivial algebra; (B.1) with d ≠ 0 replayed at 2 random Gaussian-rational points per U (224 points) | EXACT |
| B1.2 | (B.3) row marginals via 2×3 Gram/Cauchy–Binet (p.38 display) | Gram matrix entries and det = s − \|u_i0\|² on all 112 exact U | EXACT (samples) + PROOF-READ |
| B1.3 | (B.3) column marginals: Σ_i \|u0k uil − u0l uik\|² = \|u0k\|²+\|u0l\|² (p.38) | exact on all U | EXACT (samples) + PROOF-READ |
| B1.4 | (B.3) total 2s | exact on all U; also symbolic: column-3 marginal is implied by the others | EXACT |
| B1.5 | (B.4) W doubly stochastic, W00 = 0 | exact on all U | EXACT |
| B1.6 | Birkhoff decomposition, π(0) ≠ 0 on support | classical (CITED); replaced by an exact greedy decomposition on every sample (always ≤ 9 terms; reconstructs W exactly) | EXACT (samples), PROOF-READ |
| B1.7 | φ involution; m(φ(π)) = lower-right block of P_π (p.38) | enumerated over all 18 π with π(0) ≠ 0 | EXACT (complete) |
| B1.8 | (B.5) and (B.2) | exact on all U | EXACT |
| B1.9 | (B.7) principal-minor + Cauchy–Binet expansion | generic symbolic identity with U and conj(U) as independent 4×4 symbol matrices (a0 = e0 = 0) | EXACT (symbolic, general) |
| B1.10 | complementary minor \|det U[{0,k,l},[3]]\|² = \|u_i0\|² | exact on all U (classical Jacobi identity for unitary U) | EXACT (samples) + PROOF-READ |
| B1.11 | (B.8) from (B.7) and (B.3); coefficient ½ (p.39) | symbolic, with m_ij, \|u0j\|², \|ui0\|², s as free symbols subject to (B.3) | EXACT (symbolic, general) |
| B1.12 | ∏(a_i − e_σ(i)) = Σ H_ij m_ij(σ) (p.39) | all 24 σ, symbolic | EXACT (complete) |
| B1.13 | averaging gives (B.1) | follows from B1.11, B1.12 and (B.5) by linearity; (B.1) itself proved on the {0,1}^6 grid for every sample U | EXACT |
| B1.14 | s = 0 case (p.40) | PROOF-READ; exact on 4 block samples with \|u00\| = 1 and on the 6 signed permutation matrices fixing 0 | EXACT (samples) + PROOF-READ |
| B1.15 | face {W00 = 0} has affine dim 8; Carathéodory ⇒ ≤ 9 terms; φ injective ⇒ ≤ 10 weights | dimension 8 computed exactly; exact Carathéodory reduction implemented; max nonzero weights seen = 10 | EXACT + PROOF-READ |
| B1.16 | "In particular" MOC when a_i + b_j = 0 | immediate from (B.1) with b = −e | PROOF-READ |
| B2.1 | Σ w_ij = 1, Σ_{j≠i} w_ij = a'_i + b'_i | exact on 42 U; the first is the Lagrange identity | EXACT (samples) + PROOF-READ |
| B2.2 | subset bound Σ_{i<j∈S} w_ij = ‖p_S‖²‖q_S‖² − \|⟨p_S,q_S⟩\|² ≤ ‖p_S‖² | identity and inequality exact for all 15 subsets S on every sample | EXACT (samples) + PROOF-READ |
| B2.3 | cut capacity ≥ 1, max-flow min-cut ⇒ θ (p.40–41) | cut argument PROOF-READ (it is correct: a finite cut must contain both endpoints of its edge nodes). Exact Edmonds–Karp flow of value 1; θ ≥ 0, θ_ij+θ_ji = w_ij, out-mass a'_i, in-mass b'_i | EXACT (samples) + PROOF-READ |
| B2.4 | rank-two principal-minor formula, k = 1..4 (p.41) | generic symbolic (p, q, conj as independent symbols), all 15 index sets | EXACT (symbolic, general) |
| B2.5 | placement expectation equals that formula | symbolic in P(α∈I), P(β∈I), P(both) | EXACT (symbolic, general) |
| B2.6 | split θ_ij over the two γ-placements; Lemma identity | proved on the 320-point grid for 42 U | EXACT |
| B2.7 | coincident-a version via U*, σ ↦ σ⁻¹, (−1)^4 = 1 | proved on the grid for 10 U | EXACT |
| B3.1 | Δ_O on chart = det(diag(a) − O diag(e) Oᵀ) (needs det O = 1) | symbolic in a, e for rational O ∈ SO(4) | EXACT |
| B3.2 | relabelled B.1/B.2 laws vanish on {d = 0}; (B.9) divisibility | for 2 rational O ∈ SO(4) and all 28 brackets: Δ_O − Σ μ Π_σ = d·Q_O in ℚ[16 variables]. Control: Δ_O − Π_id is not divisible (28/28 for the generic O) | EXACT |
| B3.3 | irreducibility of d; density; Nullstellensatz | PROOF-READ (d is a 2×2 determinant in four distinct variables) | PROOF-READ |
| B3.4 | L = Σ μ Re g_σ + Re(ηQ_O/R) and lower bound | PROOF-READ (one line; uses Σ μ = 1). Numerically the remainder L − Σ μ Re g_σ stays bounded (≤ 9) as r → 0 on 300 curves | PROOF-READ + NUMERIC |

## 3. Gaps, errors and remarks

No mathematical gap or error was found in Appendix B. Remarks follow, none of which affects correctness.

1. (p.37, (B.2)) (B.2) is not an extra choice: **every** law satisfying (B.1) satisfies it. The nine H_ij are linearly independent polynomials (rank 9, EXACT), so (B.1) is equivalent to (B.5). Each σ with σ(0) ≠ 0 has exactly two ones in m(σ), and each σ with σ(0) = 0 has none. Summing (B.5) therefore gives 2 Σ_{σ(0)≠0} μ_σ = 2s. The LP solutions confirm this numerically (defect ≤ 2·10⁻⁷).
2. (p.40, last paragraph) "the Birkhoff decomposition … can retain at most nine permutation matrices" should read "can be chosen with at most nine". Carathéodory has to be applied inside the face (a polytope of dimension 8 whose vertices are the 18 permutation matrices with π(0) ≠ 0). The bound 10 then also needs φ to be injective, which it is (an involution), and needs the image of φ to avoid the added 0-fixing permutation, which it does. Wording only.
3. (p.37 and p.40) "depending only on U" means existence of some law for each U. The laws are not unique (Birkhoff and flow decompositions are not), and no continuity or measurability in U is claimed or proved. No downstream use needs continuity: see §4.
4. (p.41–42, Cor. B.3) The lower bound is local and depends on O, η and the germ. It needs only Re g_σ ≥ 0 for σ in the support of μ. It says nothing when two factors of E vanish together. The supplement note states this limitation explicitly; the paper's statement ("E = d·R, R(0) ≠ 0") encodes it correctly.
5. (p.41, Cor. B.3, numeric illustration) At a double coincidence (a0 = e0 and a1 = e1, E = d1·d2) the single-law remainder L − Σ μ Re g_σ grows like 1/r (5·10⁶ at r = 10⁻⁶ in the 300-curve test). So the corollary's mechanism really does not extend to double factors. This is not a counterexample to anything. On realizable curves, min L over the region was ≥ 0 at r = 10⁻⁶ in every instance where the region was nonempty (289 single, 295 double), as MOC itself predicts.
6. Main text versus appendix: Corollary B.3 is stated for O ∈ SO(4). For det O = −1 the chart identity B3.1 changes sign and (B.9) fails. That is consistent with the paper's reduction to SO(4) (Cor. 2.2), but it is worth recording.
7. Comparison with the supplement notes: the shared-zero note's Theorem 1 matches B.1 step for step, with c in place of e. The divisor note's Lemma 1 is B.3 in the "−η Im" convention, which equals Re(iη ·) and so is the same statement. It obtains Q by linear division in an affine patch instead of by irreducibility; the two are equivalent. Its Theorem 2 matches B.2 and states the flow criterion as an "if and only if"; only "if" is needed. I found no discrepancy between the paper and the supplement on these three statements.

## 4. How the statements are used downstream, and whether they suffice

- **Prop. 1.4(ii)** (p.4–5) cites B.1 for MOC when a_i + b_j = 0. B.1 gives this directly. It suffices.
- **Section 5 (T85)** does not use Appendix B. Its negative-part bound is (5.17). **Appendix F** (p.61) lists B.1, B.2 and B.3 as inputs to the Group VI curve arguments in the supplement.
- **Pointwise collision exclusion** (67 case-note extracts mention cross or shared-zero collisions; 11 mention within-spectrum repeated pairs). The notes say coincident spectral points "satisfy the support inequality" so that a counterexample "may be assumed eight-distinct" (for example inventory85 p.1, inventory54 p.2, inventory69 p.1). This needs only: for a fixed O, and for x ∈ X represented by spinors with one coincidence, Φ_{c(x)}(O) = Σ μ_σ (Pc(x))_σ ≥ 0. B.1 and B.2 give exactly that, because Φ_{c(x)}(O) = Re(α Δ_O) − Re F_{π*} and Σ μ = 1. B.1 and B.2 are stated for finite eigenvalues. A coincidence at ∞ needs an SL₂ change of coordinates first, which fixes every bracket and therefore every Π_σ and Δ_O. Appendix B does not say this explicitly, but it is immediate (PROOF-READ). The same pointwise statement serves the boundary hypothesis of Lemma 3.5 ("L ≥ 0 when some h_j = 0") at boundary points where a collision occurs. **Sufficient.**
- **Removable punctures (Cor. B.3)** are invoked explicitly in 7 case notes, two of them in two versions (inventory 40, 54, 55, 60, 62, 80, and 85 in both the final-referee and the corrected root version). In each, the note asserts that exactly one factor of a four-bracket denominator E vanishes at the puncture, and it treats simultaneous-factor corners separately (for example inventory85's "y = 1" corner). Lemma 3.5 needs −L bounded above on Ω near each puncture, and Ω ⊂ {all Re g_σ > 0}. Cor. B.3(b) gives exactly this, for a fixed O (O is fixed in Lemma 3.5), for any vanishing order, and even when other brackets outside E also vanish (inventory 40 and 55 rely on this; correct, since (B.9) is an ambient polynomial identity). **Sufficient**, provided each case verifies three things:
  - (i) The normalized gaps on its curve are exactly η(Π_σ − Π_τ)/E, with a constant η and the same holomorphic lifts in Π_σ and E. E must have multidegree one in every spinor so that rescaling a lift cancels.
  - (ii) R(0) ≠ 0 at that puncture.
  - (iii) Every puncture with two vanishing factors has its own argument.

  These are case-level claims and were not checked here. Not needed and not provided: continuity of μ in U, uniformity in O, any statement at a double coincidence.

## 5. Verdicts

- **Theorem B.1: proved as written.** That includes "at most ten" and (B.2). The latter holds for every representing law, not only the constructed one. EXACT for all identities, PROOF-READ for Birkhoff and Carathéodory (classical).
- **Lemma B.2: proved as written.** EXACT for all identities, PROOF-READ for max-flow min-cut (classical; replaced by explicit exact flows on samples).
- **Corollary B.3: proved as written.** (B.9) is EXACT on 2 rational O ∈ SO(4) for all 28 brackets, and the general argument is PROOF-READ. The scope is a single vanishing denominator factor and a fixed O. Downstream case proofs must, and say they do, handle double-factor punctures separately.

## 6. Search sizes (NUMERIC, seeds fixed in the script)

- B.1, paper's construction in floats: 20000 Haar U plus 305 degenerate U. The degenerate set covers signed permutations, 2+2 blocks, permuted blocks, |u00| = 1, u00 = 0, near s = 0 (ε down to 10⁻⁹), near-permutation U, real orthogonal U and the DFT matrix. Max relative residual on the proving grid 7.6·10⁻¹³.
- B.1, independent LP feasibility (HiGHS; the LP does not use the reduction to m_ij): 3000 Haar plus 305 degenerate U. All feasible, residual 0, LP-vertex support ≤ 10.
- B.2, LP on the 320-point grid: 600 Haar plus 102 degenerate U. All feasible.
- B.3 curves: 300 single-coincidence and 300 double-coincidence germs, radii 10⁻² to 10⁻⁶, about 6000 to 22000 angles per radius, concentrated near the boundary rays of the region.
