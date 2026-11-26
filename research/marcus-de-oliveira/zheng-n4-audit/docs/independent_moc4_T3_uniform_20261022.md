# Audit: independent_moc4_T3_uniform_20261022

Evidence level: EXACT for the finite linear algebra and identity (3); EXACT symbolic checks for the chart factorizations and the inverse chart (6); EXACT evaluation at 297 random rational points for the 24 gap formulas and identities (7)-(11) (probabilistic identity test, see Deviations); PROOF-READ for the sign/branch analysis; CITED for the SOC inequality (1) and the boundary theorems; NUMERIC for falsification.

PDF: `proofs/manin_transfer_20261021/independent_moc4_T3_uniform_20261022.pdf`, 4 pages. This is the primary proof of paper case T_3 = {0,1,2,3,7,8} (Theorem 1, p. 1). Batch C01-T1-T3.

## Deviations and limits

- **Chart identities.** Identities (7)-(11) and "all contacts of T have zero gap" are identities between rational functions of (x, t, u, r, μ, Re v, Im v). Full sympy simplification ran past 400 s and was abandoned. Instead, each identity was evaluated with exact Gaussian-rational arithmetic at 297 random rational points (seed 20261125). One third of the points have t, x < 0. Every fifth point is on the stratum r = −1, u = X, which is where (11) lives. This is exact arithmetic, but as a proof of an identity it is probabilistic (Schwartz-Zippel): the chance of a false pass is negligible but not zero. The ratio factorizations D13/D02, D23/D01, D28 and D17, and the inverse chart (6), were checked fully symbolically.
- **Sign convention for ρ.** The point C_3 = ∞ contributes the bracket det((A_i,1),(1,0)) = −1. With that sign we reproduce the note's D28 and D17 exactly, and (7)-(11) then hold with ρ = −μ Im v. Without the sign, (7)-(11) hold with ρ = +μ Im v, but D28 and D17 flip sign. The note is therefore internally inconsistent by one global sign. This is harmless, because ρ is a free nonzero real parameter and the proof derives its sign (ρ < 0) from g_4 > 0. It is recorded as a presentation gap.
- **SOC inequality (1).** Its proof was not audited; it is CITED. Batch document final_referee_T2_scalar_compensation_20261022 proves only the scalar inequality, and that note says explicitly that the bridge to Leibniz sums is not given there. Our own evidence for (1) is numeric only.
- The supplement's checker and JSON were not executed or used. The facet numbering equals the paper's (see the six-contact audit).

## A. Statements

**Theorem 1.** Take any spectra A, C = −b ∈ C^4 and κ ∈ C with F_σ = κ ∏(A_i − C_σ(i)) and g_σ = Re(F_σ − F_0). Suppose g_σ ≥ 0 for all σ and g_σ = 0 for σ ∈ T = {0,1,2,3,7,8} (non-strict contacts). Then Φ_O(g) ≥ 0 for every O ∈ O(4).

For eight distinct projective points, g lies in cone(squared minors) + {a f1 + b f11 + c f17 : a, b, c ≥ 0, ab ≥ c²}. Repeated spectra and A_i = C_j are handled by the cited boundary theorems.

This is (*) for paper T_3 (the contact set is literally T_3, EXACT). It covers all O and all spectra, given the citations.

## B. Dependencies

| ref | kind | used for | status |
|---|---|---|---|
| (1): Φ_O(a f1 + b f11 + c f17) ≥ 0 when ab ≥ c², a, b, c ≥ 0 ("144-orbit SOC") | common-lemma | the A < 0 case | CITED. Group check (EXACT): {f1, f11, f17} lies in the same G-orbit of triples, with h ↦ h, as {f9, f1, f22} of final_referee_T2_scalar_compensation and {f11, f39, f17}. The orbit has 144 elements (stabilizer 8). The scalar inequality there is proved, but the Leibniz bridge is not. NUMERIC: the minimum over SO(4) × τ of Φ(τ f1 + f11/τ + f17) is 2.4e-16 (80 BFGS restarts), so the inequality is tight and was not violated. |
| repeated-spectrum theorem | common-lemma (paper Prop. 1.4(i), Thm C.1) | two equal points within a spectrum | CITED |
| shared-zero theorem, A_i = C_j | common-lemma (paper Prop. 1.4(ii), Thm B.1) | cross collisions | CITED |
| facet list and Lemma 4.2 | core-doc | source facets | EXACT (lib) |
| two-orthogonal theorem | core-doc (paper Thm 1.2) | unitary transfer only | CITED, not needed for (*) |

## C. Steps

| # | claim | kind | status |
|---|---|---|---|
| 1 | The source facets are {0,1,2,3,4,5,6,11,15,16,17}; only f17 is non-minor; the span has dim 8; the coordinates (g4, g6, g9, g10, g12, g13, g16, g18) are injective on the span | finite | EXACT |
| 2 | cone(10 minors) has 12 facets. Exactly one is negative on f17, and it equals A = −g4 + g6 + g10 − g12 + g13 + g18. The other 11 are ≥ 0 on all 11 sources. So minor cone = N_T ∩ {A ≥ 0} | finite | EXACT |
| 3 | identity (3) on the whole span | identity | EXACT |
| 4 | A < 0 and A² ≤ g6·g22 ⇒ (3) + (1) ⇒ Φ ≥ 0 | logic | PROOF-READ (+ CITED (1)) |
| 5 | Normalization A = (z,0,1,x), C = (w,v,c,∞) exists for 8 distinct points, up to a common factor | algebra | PROOF-READ |
| 6 | D13/D02 = x and D23/D01 = c/v; contacts ⇒ x, t real | identity | EXACT (symbolic) |
| 7 | D28 = −tv(v−w)(z−1) and D17 = z(v−w)(x−tv), with the −1 bracket convention; r, u real and nonzero | identity | EXACT (symbolic) |
| 8 | Chart (6) inverts (5). The denominators are nonzero. x, t ≠ 1 and x, t, u, r ≠ 0 follow from distinctness | identity + logic | EXACT + PROOF-READ. In the D > 0 argument, t+u = u+x = 0 actually gives z = v (A_0 = C_1), and then w = z. The note's "z = w" is also true. |
| 9 | All 6 contacts of T have zero gap in the chart | identity | EXACT (random-point) |
| 10 | v real ⇒ all gaps vanish | identity | EXACT (random-point) |
| 11 | (7): g4, g6, g9 = −x g6, g10 = (r+1) g4, g14, A = g13 + g22/t | identity | EXACT (random-point); the global sign of ρ is a convention, see Deviations |
| 12 | t > 0 ⇒ A ≥ 0. t < 0 ⇒ x < 0, ρ < 0, r ≥ −1, ru < 0 | sign logic | PROOF-READ |
| 13 | (8) k/T = e − j, A = f − k/T, g21 = Xb + A; (9); both identities in (10) | identity | EXACT (random-point) |
| 14 | r > 0 branch ⇒ X > T ⇒ f > e ⇒ A > 0 | sign logic | PROOF-READ |
| 15 | −1 ≤ r < 0 branch: T < X is impossible (u > X ⇒ j < 0; r = −1, u = X uses (11)) | sign logic | PROOF-READ; (11) EXACT (random-point) |
| 16 | T ≥ X, A < 0 ⇒ α ≤ k/T and α ≤ Xb ≤ Tb ⇒ A² ≤ g6 g22 | inequality | PROOF-READ |
| 17 | Case split exhaustive: distinct points / repeated within a spectrum / A_i = C_j; v real or not; t > 0 or t < 0; r > 0 or −1 ≤ r < 0 | case split | PROOF-READ (r ≠ 0, t ≠ 0 by distinctness; r ≥ −1 forced) |

## D. Falsification (NUMERIC)

- Chart sampling with seed 7: 102473 random chart points, of which 400 had all gaps ≥ 0. Of those, 8 had A < 0, and the minimum of g6·g22 − A² over them was 5.1e-3 > 0, as step 16 predicts. The minimum of Φ over SO(4) at 25 feasible points (8 restarts each) was 3.9e-15. No violation.
- SOC (1): minimum 2.4e-16 over SO(4) × τ. Over 5000 random O, the minimum of Φ(f17) + 2√(Φ(f1)Φ(f11)) was 4.3e-4.

## E. Verdict

**confirmed-modulo-citations.** The citations are SOC inequality (1), the repeated-spectrum theorem (paper Thm C.1), the shared-zero theorem (paper Thm B.1), and Thm 1.2 for the unitary form.

Presentation gaps:
- p. 2-3, (7): ρ = μ Im v has a global sign that is inconsistent with the note's own D28 and D17 convention. Harmless.
- p. 2: the D > 0 justification names the wrong coincidence. The conclusion is still correct.

Caveat on strength: the chart identities (7)-(11) were checked by exact evaluation at random points, not by symbolic simplification.
