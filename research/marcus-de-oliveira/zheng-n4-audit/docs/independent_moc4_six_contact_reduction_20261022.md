# Audit: independent_moc4_six_contact_reduction_20261022

Evidence level: EXACT (all identities, counts and coefficient transfers), PROOF-READ (Lemma 2, Lemma 4, Theorem 5, the case logic of Props. 10 and 13), CITED (paper Prop. 1.4(i)/Thm C.1 for repeated spectra, Section 3 / HMSV for (M1), Thm 1.2 for the unitary forms), NUMERIC (falsification).

PDF: `proofs/manin_transfer_20261021/independent_moc4_six_contact_reduction_20261022.pdf`, 8 pages. Primary proof of paper cases T_1 (Prop. 10, p. 6) and T_2 (Prop. 13, p. 7). Batch C01-T1-T3.

## Deviations and limits

- None of the supplement's scripts were run (Remarks 1, 9, 11, 14 name them). Every check below is our own code in `replay/independent_moc4_six_contact_reduction_20261022__check.py`, using the shared library `lib/moc4lib.py`.
- The supplement's facet numbering (`root_permutation_polytope_complete_20261022.json`, loaded as data once during exploration) agrees index-for-index with paper Appendix D / `moc4lib.F_TABLE` (all 82 value vectors are identical, and the 24 permutations come in the same lexicographic order). So "f_j" means the same facet in the note, the paper and the library. The final replay script does not read the supplement.
- The cubic F of (M2) is taken to be the paper's C of (A.2). The note does not write F down.
- Facet completeness (82 facets, Remark 9) is taken from the library's own exact certificate, not from the note's 6846-simplex certificate.
- (M1) is not re-proved: irreducibility, the dimension upper bound, and "every nonzero point of X is λ times the matching vector of finite spectra". These are CITED to paper Section 3, GIT and HMSV.

## A. Statements

Permutations are indexed 0..23 in lexicographic one-line order, and F_σ = λ ∏_i (a_i + b_σ(i)) = λ ∏_i (A_i − C_σ(i)) with C = −b. Throughout, "cost" means the gap vector g_σ = Re(F_σ − F_τ) ≥ 0 for a minimizing τ.

- **Lemma 2 (p. 2).** C is an irreducible affine curve, f_1..f_m and h are regular on C, not all f_j are constant, and Σ Re f_j ≡ 1. Write g_j = Re f_j and L = Re h. Assume L is bounded on {all g_j ≥ 0}, and assume L ≥ 0 at every finite point with all g_j ≥ 0 and some g_j = 0. Then L ≥ 0 wherever all g_j > 0.
- **Lemma 4 (p. 3).** If x + Cκ ⊂ X, then Re Φ_O(x) ≥ min_σ Re(Px)_σ for every O ∈ O(4).
- **Theorem 5 (p. 4).** Fix O ∈ O(4). Suppose the support inequality Re Φ_O(x) ≥ min_σ Re(Px)_σ holds for every x ∈ X with at least 6 contacts. Then it holds for every x ∈ X. This is the paper's Theorem 3.6 in contrapositive form.
- **Corollary 7.** Under the same hypothesis for every O ∈ O(4), MOC(4) holds for all normal pairs. This needs paper Thm 1.2.
- **Prop. 8 (p. 5).** There are 214 G-orbits of 6-subsets, and exactly 90 have a non-minor generator. The others are nonnegative or carry only the zero cost. Also claimed: 113 orbits are exact six-vertex faces; all 90 representatives are exact faces; the non-minor generator counts are 1, 2, 3, 4 in 51, 35, 3, 1 orbits; 9 of the 90 orbits contain an image of {2,7,19,20}.
- **Prop. 10 (p. 6) = paper T_1 = {0,1,2,3,6,7}.** Take any spectra a, b ∈ C^4 (repeated eigenvalues and a_i + b_j = 0 allowed), any λ, and suppose the contact set contains T_1 (non-strict). Then the support inequality holds for every unitary U. Hence (*) holds for every O ∈ SO(4).
- **Prop. 13 (p. 7) = paper T_2 = {0,1,2,3,6,8}.** Same quantifiers. Away from the degenerations of Prop. 10, the gap vector is a nonnegative combination of the 11 minor generators of T_2. So (*) holds.
- Statement check: both contact sets are literally the paper's T_1 and T_2 (EXACT). Both statements are non-strict ("contains T"), cover every O ∈ O(4), and cover every spectrum. Repeated eigenvalues within one spectrum go to the one-circle theorem (CITED, paper Prop. 1.4(i)). A zero contact product goes to the entry face. Other cross collisions a_i + b_j = 0 are carried inside the main argument; only A_1 = C_3 can actually occur in Prop. 10, and it is treated.

## B. Dependencies

| ref | kind | used for | status |
|---|---|---|---|
| (M1) X irreducible of dim 6, P rank 14, nonzero points = finite spectra | core-doc (paper §3, HMSV, GIT) | Thm 5 dimension count and realization | rank 14 and dim ≥ 6 EXACT (lib); the rest CITED |
| (M2) X ⊂ Sing C, ∇C(κ) = 0, eq. (2) | core-doc (paper (3.10)-(3.12), A.1) | Lemma 4 | EXACT: (2) re-derived symbolically; X ⊂ Sing C from lib |
| (M3) ker H_κ = span{e_ij}, rank 4 | core-doc (paper (3.11)) | Lemma 4 | EXACT |
| (M4) Φ_O(κ) = 1, Φ_O(e_ij) = O_ij^2, Φ_O(x) = λ det(...) | core-doc (paper (3.5), (3.8)) | Lemma 4, Thm 5 | EXACT in lib; the general-O proof is in lib/README |
| Two-orthogonal theorem | core-doc (paper Thm 1.2) | the unitary forms of Cor. 7 and Props. 10, 13 | CITED; not needed for (*) |
| Birkhoff's theorem | external | Lemma 4 | CITED (standard) |
| Removable singularities and the maximum principle for bounded subharmonic functions on compact Riemann surfaces | external | Lemma 2 | standard; PROOF-READ |
| Facet completeness (Prop. 4.1) and Lemma 4.2 | core-doc | Props. 8, 13 | EXACT (lib certificate) |
| One-circle theorem: one spectrum with ≤ 3 distinct points | common-lemma (paper Prop. 1.4(i), Thm C.1) | Props. 10 and 13, repeated spectra | CITED |
| Entry face {σ: σ(i) = j}: 15 generators, all minor | internal | Props. 10 and 13, zero contact product | EXACT for all 16 (i, j) |

## C. Steps

| # | claim | kind | status | how |
|---|---|---|---|---|
| L2.1 | u = −L on D, 0 elsewhere, is continuous, bounded and subharmonic on the finite normalization | analytic | PROOF-READ | At boundary points of D either L = 0 or some g_j = 0, and then the hypothesis plus continuity force L = 0. The sub-mean-value property holds off D trivially and on ∂D because u ≥ 0. |
| L2.2 | it extends across the finitely many punctures, is constant, so D is everything | analytic | PROOF-READ | Bounded subharmonic functions have removable isolated singularities; maximum principle on a compact surface |
| L2.3 | then 0 ≤ g_j ≤ 1, so the f_j have no poles and are constant, a contradiction | analytic | PROOF-READ | Near a pole of order k, Re f is unbounded above and below |
| L4 | H_κ x = 0 ⇒ Re x ∈ span e_ij ⇒ Φ_O(Re x) = Σ t_σ (P Re x)_σ ≥ min | identity + Birkhoff | EXACT for (2), (3); PROOF-READ for the rest | `M2_*`, `M3_*` |
| L4b | h_{O,τ} is a real combination of the gap forms; L is bounded on the simplex | linear algebra | EXACT | The differences of P rows have rank 13, so their common kernel is Cκ |
| T5.1 | all 24 values equal ⇒ zero gap | identity | PROOF-READ | |
| T5.2 | a component Y through x has dim ≥ 6 − r ≥ 1 | dimension | PROOF-READ | Affine dimension theorem; r ≤ 5 equations |
| T5.3 | gap map constant on Y ⇒ Y = x + Cκ ⇒ Lemma 4 | algebra | PROOF-READ | |
| T5.4 | there is an irreducible curve C ⊂ Y through x on which some gap is nonconstant | algebraic geometry | PROOF-READ (argument supplied by auditor) | Pick y ∈ Y with a different gap value. Any two points of an irreducible variety lie on an irreducible curve in it (Mumford, Abelian Varieties §6 lemma). The note's "generic hyperplane sections" is only a sketch (presentation). |
| T5.5 | on C the contact gaps stay 0 and the sum stays 1 | identity | PROOF-READ | (5) fixes the complex values; their real parts are 0 |
| T5.6 | a new zero on C means ≥ r+1 contacts, so the induction hypothesis applies; Lemma 2 then gives L(x) ≥ 0 | logic | PROOF-READ | τ remains a minimizer on C |
| P8 | 214 / 90 / 113 / all 90 are faces / 51,35,3,1 / 9 | finite count | EXACT | `Prop8_*` |
| P8b | cost vanishing on T = nonnegative combination of the facets in A(T) | convexity | EXACT (lib) + PROOF-READ | {c : Pc ≥ 0} is pointed and its extreme rays are the 82 facets |
| P10.0 | zero contact product ⇒ the whole entry face is in the contact set, and its 15 generators are all minor | identity + count | EXACT | `entry_faces_*` |
| P10.1 | the six factorizations and the three ratios; nonvanishing comes from contact factors and within-spectrum distinctness only | identity | EXACT | `P10_D**`, `P10_ratio*` |
| P10.2 | Möbius normalization keeps the products up to a common factor | algebra | PROOF-READ | (A−C) ↦ det·(A−C)/((rA+s)(rC+s)), and the pole is chosen on the circle away from the 8 points |
| P10.3 | imaginary numerators xy(x−y)v and −tu(u−t)q | identity | EXACT | |
| P10.4 | u = 0 (A_1 = C_3) case: D01, D02 formulas ⇒ Re L = 0 ⇒ q = 0 | identity + logic | EXACT + PROOF-READ | |
| P10.5 | Re L = 0 ⇒ all σ with σ(0) = 0 are contacts ⇒ entry face | logic | PROOF-READ | |
| P13.1 | D02 = λ(A0−C0)(A3−C3)(A1−A2)(C1−C2) ≠ 0; D79 D02 = D68 D13 ⇒ g7 = g9 | identity | EXACT | `P13_D79_D02_eq_D68_D13` |
| P13.2 | J = generators(T_2); f17 + f9 = f4 + f10 + f16; ℓ values; the β-transfer is an exact identity with formal coefficients | identity | EXACT | `P13_*` |

## D. Falsification (NUMERIC)

- T_1 stratum: `random_point_on_contact_stratum(T_1, margin=1e-2)` with seeds 0-3 returned 4 points whose contact set is exactly T_1. In every one, a_2 = a_3 to 1e-15. So all of them are repeated-spectrum points, which Prop. 10 sends to the one-circle theorem. This agrees with Prop. 10's claim that the nondegenerate T_1 stratum always picks up the full entry face. (The library README's "point with contact set exactly T_1" is such a repeated-spectrum point. It is not a counterexample.)
- T_2 stratum, seeds 0-5: 6 points with contact set exactly T_2. All are nondegenerate: within-spectrum separations range from 1.0e-2 to 6.0e-1. |g7 − g9| ≤ 2e-13 at every point, which confirms (7) numerically. The minimum of Φ over SO(4) (20 BFGS restarts) is ≥ 1e-16 at every point. No violation.
- The same run also covers T_3 stratum points (seeds 0-3): the minimum of Φ is ≥ 1e-14.

## E. Verdict

**confirmed-modulo-citations.** Every identity, count and coefficient transfer was replayed exactly, and every analytic step was read line by line. The citations are (M1) realization, irreducibility and dim ≤ 6 (paper §3, HMSV); the one-circle theorem (paper Prop. 1.4(i), Thm C.1); paper Thm 1.2 (only for the unitary forms); and Birkhoff.

Presentation notes, no substantive gap:
- p. 4, existence of the curve C: the "successive generic hyperplane sections retaining x" step is only sketched. The auditor supplied the argument (see T5.4).
- p. 4, Remark 6 says strictly semistable configurations satisfy the inequality "for the same reason". This is not used by the proof, because Theorem 5's hypothesis already includes those points.
- p. 6, Prop. 10 asserts that the six differences used are nonzero. This holds using only contact factors and within-spectrum distinctness (checked).
