Evidence level: EXACT for every displayed identity, factorization, resultant, Newton-polygon and series claim, the normal cone and all boundary-tie reductions; PROOF-READ for the puncture analysis, the maximum-principle and constant-exception arguments and the Rouché perturbation; CITED for the earlier cases T_37, T_45, T_57, T_65, Appendix B (Theorem B.1, Lemma B.2, Corollary B.3), Lemma 4.2 / Prop. 4.1, (3.11) and Birkhoff; NUMERIC only for a brief falsification search.

# Audit of `root_inventory85_native_uniform_corrected_20261023` (primary proof of paper case T_86)

PDF: `proofs/manin_transfer_20261021/root_inventory85_native_uniform_corrected_20261023.pdf`, 7 pages. Title: "A two-chart contact curve with single-divisor compensation". Auditor label C43-T86.

Replay: `replay/root_inventory85_native_uniform_corrected_20261023__check.py` (own code, nothing from the supplement executed, imported or read; uses `lib/moc4lib.py`). Output: `..._check_output.txt`, 58 EXACT checks, all pass, plus the NUMERIC section.

## Deviations and limits

- No deviation from the assigned scope. I did not use the supplement's verifier `final_referee_inventory85_native_uniform_exact_20261023.py` or its JSON. Every identity was re-derived from the definitions in the native chart (A = (0,1,z,w), C = (∞,v,x,t), F_σ = ∏ [A_i, C_σ(i)]).
- The cited earlier cases (T_37, T_45, T_57, T_65 and their alternatives) and the single-collision lemma (= paper Corollary B.3) were not re-audited here. Appendix B is taken from the audit's `foundations/appendix-b.md`. I checked only that each use matches the statement.
- The main statement was searched numerically only briefly, as the batch prompt instructed. The effort went into the intermediate identities and the puncture list.
- Equation numbers below are the document's own. The document reuses (1)-(12) for two different sets of displays. For example, "(3)" is both the bound θ ≤ g23 on p. 1 and Res_y(N, A) on p. 4. I write (3θ) for the first and (3R) for the second.

## A. Statements

**Theorem 1 (p. 1).** Let T = {4, 5, 8, 9, 12, 13} (0-based lexicographic indices). For every x ∈ X with P c(x) ≥ 0 and (P c(x))_π = 0 for all π ∈ T (contacts non-strict, "tied permutations include (1)"), and for every real orthogonal O: Φ_{c(x)}(O) ≥ 0. The spectra are arbitrary complex, and repeated eigenvalues and a_i + b_j = 0 are included: the proof handles them at the outset.

- Code: `CONTACT_SETS[86] == T` (check A1). So this is literally the paper's T_86, with the identity group element.
- Strength: at least (*) for T_86. Contacts are non-strict. All O ∈ SO(4) are covered, and the O(4) form is equivalent by the paper's column-flip remark. All spectra are covered: degenerate spectra by Prop. 1.4 / Appendix B (see B.6 below), eight-distinct spectra by the curve argument.

There are no other labelled statements. The single-collision lemma is used, not proved, here.

## B. Dependencies

| # | dependency | kind | where used | check of the match |
|---|---|---|---|---|
| B1 | Lemma 4.2 (contact cone = cone over the f_j vanishing on T) and Prop. 4.1 (82 facets) | core (paper, replayed in lib) | the cone (2), the bound θ ≤ g23, the minor cone at the six other indices | `generators(T)` = [0,1,5,45,59,67,68,70,72], equal to (2); rank 8; the contact rows have rank 6 (A2-A5) |
| B2 | minor vectors have Φ = det(O[I,J])² ≥ 0 (4.2) | core | the seven minor rays | names f0 = n1,4, f1 = n12,34, f5 = n12,24, f45 = n13,14, f59 = n14,14, f67 = n2,2, f72 = n12,12 match the document's 0-based row/column table (A3) |
| B3 | earlier six-contact cases for the ties j ∈ {2,3,6,7,10,11,14,15,16,17,18,19} | paper cases | ordinary boundary (4) | for each j some g ∈ G maps an earlier-stage (Table 6, stages 1-7; T_86 is stage 8) case into T ∪ {j} (A11): j = 2, 3 → T_37 (stage 1; T_39 also fits); j = 6, 7, 16, 17 → T_65 (stage 6; T_66, T_69, T_75 also fit); j = 10, 11 → T_57 (stage 6; T_58, T_78 also fit); j = 14, 15, 18, 19 → T_45 (stage 7; T_46, T_69, T_76, T_81 also fit). Explicit g are printed in the output. CITED |
| B4 | ties at j ∈ {0, 1, 20, 21, 22, 23} | internal | ordinary boundary | f68 and f70 are both > 0 there (A8), and every generator of T ∪ {j} is a minor (A9). EXACT |
| B5 | Corollary B.3 (single-collision-divisor compensation) | common lemma (paper App. B; established in `foundations/appendix-b.md`) | punctures where exactly one factor of (16) vanishes | form matches: g_i = −η Im H_i = Re(iη(Π_i − Π_12)/E) and L_O = Re(iη(Δ_O − Π_12)/E) with τ = 12, complex scalar iη, E = Π_13 − Π_12 in the same primitive lifts (S16), Δ_O = Σ sgnσ ∏O_{iσ(i)} Π_σ (S23). The cofactor R(0) ≠ 0 is verified fibre by fibre (C15-C18). The two-factor corner y = 1 is handled separately (C19). CITED |
| B6 | Lemma B.2 (within-spectrum pair law) and Theorem B.1 (shared zero), i.e. Prop. 1.4 | common lemma (established) | spectra that are not eight-distinct | each gives μ with det = Σ μ_σ z_σ for every spectrum with the given coincidence, hence Φ ≥ 0 whenever P c ≥ 0. CITED |
| B7 | (3.11): ker Hess C(κ) = the 10-dimensional singleton span E; ∇C(κ) = 0 | core (paper, replayed by lib) | the constant exception | statement as needed. CITED (lib EXACT) |
| B8 | Birkhoff-von Neumann | external, classical | singleton cost | (O_ij²) is doubly stochastic. CITED |
| B9 | removable singularities for bounded subharmonic functions; maximum principle on a compact Riemann surface; Rouché | classical analysis | max-principle step, perturbation | hypotheses checked (C21, C23) |

## C. Steps

| # | claim (location) | kind | status | how checked |
|---|---|---|---|---|
| C1 | (2): the normal cone of T has the nine generators 0,1,5,45,59,67,68,70,72; rank 8; contact rank 6 (p. 1) | identity | EXACT | A2, A4, A5 |
| C2 | all except 68, 70 are minors, with the rows/columns printed (p. 1) | identity | EXACT | A3 |
| C3 | (3θ): θ = c68 + c70 ≤ g23 for every nonnegative representation (p. 1) | inequality | EXACT | f68(23) = f70(23) = 1 and minors are ≥ 0, so g23 = Σ c_j f_j(23) ≥ θ (A6) |
| C4 | each additional tie in (4) contains an earlier proved six-contact case (p. 1-2) | reduction | EXACT + CITED | A11 (see B3) |
| C5 | at the other six noncontact indices both non-minor generators are positive, so a tie there kills both non-minor coefficients (p. 2) | identity | EXACT | A7-A9 |
| C6 | repeated values and cross collisions satisfy the inequality (p. 2) | citation | CITED | B6 |
| C7 | Möbius normalization (5): products and response pick up one common scalar (p. 2) | identity | PROOF-READ + EXACT | each F_σ contains every spinor exactly once, brackets are SL_2-invariant, and the response is the Leibniz combination Δ_O = Σ sgnσ ∏O_{iσ(i)} F_σ (S23, exact on Gaussian-rational spectra). The scalar is absorbed into α |
| C8 | (6) E = F13 − F12 = x(t−v)(z−w); (16) E = ±[A0,C2][A1,C0][A2,A3][C1,C3] in lifts (p. 2-3) | identity | EXACT | S1, S16 |
| C9 | contacts ⇔ r, s, u, h real, with H4 = u − r, H8 = h − s, H12 = 0, H13 = 1 (7) (p. 2) | identity | PROOF-READ | immediate from (7). The gaps −η Im H_i are the paper's Re(α(F − F_12)) with α = iη/E |
| C10 | (8) universal inverse; y, y−1, N, L ≠ 0 at eight-distinct points; (9) (p. 2) | identity | EXACT + PROOF-READ | S3: r = (t−1)(x−v)/(x(t−v)), s = v(1−x)/(x(t−v)). With v = yx the two ratio equations are linear in (x,t), with determinant yL (S5). Their solution is (8) (S6). If L = 0, consistency forces N = 0 (S8, adjugate), which (9) (S7) excludes together with r ≠ 0. y ≠ 0, 1 and N ≠ 0 come from v ≠ A0, v ≠ x and x ≠ A0 |
| C11 | (10) H9 = h ⇔ C9 = W_D w − W_N = 0 (p. 3) | identity | EXACT | numer(H9 − h) = ±L·N·C9 (S9) |
| C12 | bilinear determinant (h−s)(h−1)K² (up to sign); (11); (12); h = s ⇒ z = y and second denominator L(w−t); h = 1 ⇒ w = y and first denominator L(z−t); hence the charts cover every eight-distinct source (p. 3) | identity + case split | EXACT + PROOF-READ | S10-S15. Exhaustive: chart 1 fails only if W_D = 0, and then C9 = 0 forces W_N = 0, which needs a zero bilinear determinant (K ≠ 0), i.e. h ∈ {s, 1}. At h = 1 this needs z = t (excluded); at h = s it needs z = y |
| C13 | (13)-(14): in chart 1, H5 = u is a 174-term equation of bidegree (4,2), with no discarded factor (p. 3) | identity | EXACT | numer(H5 − u) in the lifts (13) = K·L·C, with C of bidegree (4,2) and 174 terms (S17). The discarded factors K and L are nonzero at every eight-distinct point |
| C14 | (15) generic set proper (nonzero at (2,3,5,7)); Res_y(C_z², C_z) ≠ 0 excludes vertical components (p. 3) | identity + PROOF-READ | EXACT | G1. A vertical component y = y0 would need C_z²(y0) = C_z(y0) = 0 |
| C15 | generic y: only [A2,A3] can vanish in (16) (p. 4) | analytic | PROOF-READ | in the lifts (13), [A0,C2] = N, [A1,C0] = 1, [C1,C3] = −LK (S16). (W_N, W_D) is primitive for y ≠ 1, generic h (bilinear determinant). Points z = ∞ or w = ∞ are handled by primitive local lifts and are at most a zero of [A2,A3] |
| C16 | y = 0 and y = ∞: C2 = ∞, resp. C1 = ∞, and the other two indicated brackets stay nonzero (p. 4) | analytic | PROOF-READ (argument supplied by auditor) | y = 0: N(0) = −s, L(0) = 1 − r − s ≠ 0, K = 1. y = ∞: x → −1/(s−1), t → (s−1−r)/(s−1), v → ∞. Only [A2,A3] can vanish |
| C17 | roots of N: z finite, z ≠ w, w = ∞ allowed, only [A0,C2] vanishes (p. 4) | analytic + identities | EXACT + PROOF-READ | (3R) R1, (4) R2, (5) R3, (6) R4, (17) S13. N and L have no common root (Res = −r²s, R5), and N(1) = −r ≠ 0, so L, K ≠ 0 and t − v = K/L ≠ 0. C_z² ≠ 0 there (needs s + 1 − h ≠ 0, in (15)). z ≠ y and Lz ≠ M give z ≠ w by (17) |
| C18 | root of L: (7), (8), (18), Res_y(L,(r−1)y+s); finite branch: only [C1,C3] vanishes; infinite branch: rescaling (19), unit denominator (20) (p. 4-5) | identities + analytic | EXACT + PROOF-READ | R6-R10, L1-L6. At L = 0: N = −sK, M = (1−s)K, so all of C1..C3 are at ∞ (lift (N,0),(N,0),(M,0)). [A2,A3] = (z−y)(Lz−M) ≠ 0 on the finite branch. Z* = −B_L/F_L ≠ 0, M_L. L·w → Z*(1−h)K/(Z*+(s−h)K) (L3). C' = (∞, N, N/y, M) (L4). Rescaled denominator ±N Z*(Z*−M)K ≠ 0 (L6) |
| C19 | corner y = 1: Newton support deg_ε + 2deg_b ≥ 4, weight-4 part ε⁴P(b/ε²) with P as printed (21); (23); expansions (11), (12), (24); E of exact order 4; F12, F23 of order ≥ 4; hence H23, θ and L_O⁻ are bounded (p. 5) | identities + analytic | EXACT + PROOF-READ | Y1-Y7. C(1,z) = r²(s+u−h)(z−1)² (S18), so both branches have z → 1. P has nonzero leading coefficient and P(0) ≠ 0, so the lower Newton boundary is the single edge (4,0)-(0,2), and both branches are b = Yε² + O(ε³) with Y a simple root (disc P ≠ 0 in (15)). w − y = (z−y)(1 − (Lz−M)/W_D) gives W independently of the higher terms of z (Y5). W − Y = Y(rY−(s−1))/(s−h−rY) ≠ 0 (Y6). The lower bound L_O ≥ −θ·max(|Φ_f68(O)|, |Φ_f70(O)|) ≥ −const·g23 holds on the supported part |
| C20 | single-divisor lemma applies at every puncture with exactly one vanishing factor; probabilities fixed for O (p. 3-4) | citation | CITED (B5) | form and cofactor checked as in B5 |
| C21 | compact maximum principle: v = −L_O on a negative strict-support component, 0 elsewhere, is continuous and subharmonic, bounded at every puncture, hence constant (p. 5) | analytic | PROOF-READ | Boundary points with E ≠ 0 have L_O = 0: a boundary point with L_O < 0 would be a new tie or a collision point of X with P c ≥ 0, covered by C4-C6. Gluing max(−L_O, 0) is standard. Every point of the compact normalization lies over some y ∈ P¹, and C15-C19 list every fibre, so the list of punctures is complete |
| C22 | constant exception: a positive constant makes the component all of the curve, H_i have no poles (Im of a Laurent leading term changes sign), hence are constant; x/E = fκ + v0; f nonconstant ⇒ ∇C(v0) = 0 and Hess C(κ)v0 = 0 ⇒ singleton cost ⇒ Birkhoff; f constant ⇒ (26) makes the curve a point (p. 5-6) | analytic + identities | PROOF-READ + EXACT + CITED | ∇C(v0 + fκ) = ∇C(v0) + f·Hess C(κ)v0 + f²∇C(κ) (symmetric trilinear form), PROOF-READ; (3.11) CITED (B7); (26) EXACT (S22) |
| C23 | nongeneric real parameters: the germs are nonzero (C(1,z) = r²(s+u−h)(z−1)²; [y⁴z]C + [y⁴]C = −rs at h = s+u; numer(H5 − u) = −L·C5, C5 with 88 terms; [y⁰w⁰]C5(y,y,w)|_{h=s} = s²); Rouché gives a nearby generic counterexample (p. 6) | identities + analytic | EXACT + PROOF-READ | S18-S21. Chart 1 is valid at every eight-distinct source except h = s, z = y (C12), where chart 2 reduces to z ≡ y (S15). The perturbed zero satisfies all four ratio equations because the numerators are K·L·C, ±L·N·C9 and −L·C5 with nonvanishing cofactors. The contacts are exact because r', s', u', h' are real; strict gaps, eight-distinctness and L_O < 0 are open conditions. Generic real parameters are dense, since (15) is a proper algebraic set |

Presentation notes (no gap):
- p1 (p. 1-6) Equation numbers (1)-(12) are used twice, so "By (3)" on p. 5 means the θ-bound of p. 1, not the resultant (3) of p. 4. Recoverable.
- p2 (p. 2) The proof that L ≠ 0 is compressed into "The only potential ambiguity, L = N = 0, is impossible". The missing step is that the (x,t) system has determinant yL and adjugate component N (C10). Recoverable in two lines.

## D. Falsification (NUMERIC)

- Main statement: 8 points on the T_86 stratum (lib `random_point_on_contact_stratum`, π* = 12, seeds 0-7, margin 1e-3 for odd seeds, so the contact set is exactly T for 4 of them). Minimum of Φ_c over SO(4) by 15-restart BFGS per point. Worst value −4.3e-13, at round-off level and attained at contact-permutation matrices. No violation (threshold −1e-8).
- Chart round trip: 400 points (seed 11), random real (r,s,u,h) and complex y, both roots z of C(y,·) = 0, reconstructed through (8) and (11). The four ratios come back with relative error ≤ 7.3e-12. This confirms the parametrization numerically, in addition to the exact checks.
- The intermediate inequalities (θ ≤ g23, the boundedness at y = 1) are proved exactly or by explicit order counts, so there was nothing further to break numerically.

## E. Verdict

**confirmed-modulo-citations.** Every displayed identity, resultant, factorization and series claim was replayed exactly and is correct. The puncture list is complete over all of P¹_y: generic y, y = 0, ∞, the roots of N, the root of L (both branches), and the y = 1 corner, which is the only fibre with two vanishing factors. The corner is handled not by Corollary B.3 but by the exact θ ≤ g23 bound together with exact order counts. Each use of Corollary B.3 has the required form, the same lifts and a nonvanishing cofactor. The chart cover and the Rouché perturbation are sound.

Citations: T_37, T_45, T_57, T_65 (or the listed alternatives, all in stages 1-7); Appendix B (B.1, B.2, B.3; established); Lemma 4.2 / Prop. 4.1 and (3.11) (paper, lib-replayed); Birkhoff. No gap found.
