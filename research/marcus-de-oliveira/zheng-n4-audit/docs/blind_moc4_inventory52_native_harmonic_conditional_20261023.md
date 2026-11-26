Evidence level: EXACT for the chart identities, the amplitudes, the curve restrictions (6), the resultants, the infinity expansions (8) and the corner formula, rebuilt from the definitions with sympy. EXACT for the tie reductions, in the whole_assembly replay. PROOF-READ for the harmonic / minimum-principle argument, Lemma 2 and the Rouché recovery. CITED: the seven-contact hypothesis, earlier paper cases, Appendix B, Birkhoff.

# blind_moc4_inventory52_native_harmonic_conditional_20261023 (5 pp.)

Replay: `replay/blind_moc4_inventory52_native_harmonic_conditional_20261023__check.py`. Output: all PASS.

## Deviations and limits

- Convention used: F_sigma = prod_i (A_i - Cpt_sigma(i)), with the factor at the infinite point dropped. With this choice E = F0 - F6 = t(y - z)(w - 1) holds exactly as printed.
- Not replayed: the explicit homogeneous inverse (7) (E = -J W^2 (W - 1) H L), and the real-Taylor forms (10) and (11). Only F0/E = -rT/H was replayed.
- The claim that all 24 coordinates are combinations of the eight amplitudes and H9 was checked only through the generator identity of the cone (exclusive gaps). That suffices, because the gap vector lies in the 8-dimensional contact space.

## A. Statement

**Theorem 1.** Assume (*) for every x with contacts containing T7 = T6 ∪ {17}. Then (*) holds for T6 = {0,1,4,6,9,12}, which is paper T_53 literally. All spectra are covered:
- degenerate spectra (a_i + b_j = 0, repeated eigenvalues) through Appendix B (Thm B.1, Lemma B.2), which give MOC directly there;
- added ties through earlier cases and the hypothesis.

## B. Dependencies

Ties j ∈ {5, 7, 11, 14, 15, 18} on T6: EXACT. T6 ∪ {j} contains orbit members of earlier cases:
- j = 5 and 14: T_4;
- j = 7: T_18;
- j = 11 and 15: T_26;
- j = 18: T_42 and T_48.

Other dependencies:
- Tie 17: the hypothesis.
- Ties where f35 > 0: d = 0, so a sum of minors. EXACT: f35 vanishes off T6 exactly at {5,7,11,14,15,17,18}.
- Birkhoff (classical), Rouché (classical).

## C. Steps

| # | claim | status |
|---|---|---|
| 1 | E = t(y - z)(w - 1); (3) solves the r, s, u tree equations | EXACT |
| 2 | V ≠ 0 and the denominators are nonzero for distinct points | PROOF-READ |
| 3 | J inverts to z = rW/ZD | EXACT |
| 4 | d is the exclusive gap at 2; H_2 = J exactly; f9 and f0 amplitudes exact; f21 and f7 amplitudes up to the real constants -s and -u, which are allowed | EXACT |
| 5 | C(J, W) has bidegree (3,3); the five restrictions (6) | EXACT |
| 6 | Res_W(C, L) and Res_W(C, H) as printed; Q_H example with discriminant -15616 | EXACT |
| 7 | finite poles (J = 0, W = 0, W = 1, L = 0, W = ∞) all have real J, so d → 0 and lim inf ℓ >= 0 on the support region | PROOF-READ, from 5 and 6 |
| 8 | J = ∞: E(0, T) = T((r-1)T + s); corner -K^2 (r-1)(K(h+r) - 1); these three branches exhaust (cubic in W, leading coefficient ≠ 0; [J^3]C vanishes only at W = 1) | EXACT + PROOF-READ |
| 9 | Lemma 2: an amplitude -J + g(eps) with g real-Taylor admits no support region near the puncture | PROOF-READ. Near the puncture the support region is inside {d = 0}, where ℓ >= 0 |
| 10 | (9) is consistent with f9amp = -r - s/T | PROOF-READ (checked by hand) |
| 10' | (10) and (11) | UNCHECKED |
| 11 | height divisor H = 0: every H_sigma is holomorphic there, so ℓ_O extends harmonically and the puncture is removable | PROOF-READ. The Birkhoff argument is not needed |
| 12 | no vertical or horizontal curve components, so J is nonconstant | PROOF-READ (C(J,1) = -s^2 and [J^3]C) |
| 13 | Rouché perturbation to generic (r, s, u, h): the perturbed point keeps strict support and Phi < 0 | PROOF-READ |
| 14 | minimum principle on the punctured normalized component; a constant negative ℓ forces the support region to be the whole component, which contradicts Lemma 2 at a J pole | PROOF-READ |
| 15 | collision points on the curve are ordinary points of X (the matching vector is defined), not punctures | PROOF-READ (auditor remark) |

## D. Falsification

Only the brief numeric test in the whole_assembly report. The main statement is a special case of MOC, so it was not searched separately.

## E. Verdict

**confirmed-modulo-citations.** Cited:
- the seven-contact hypothesis (supplied by the low, middle, third and fourth notes);
- paper cases T_4, T_18, T_26 and T_42 / T_48;
- Appendix B.

Unreplayed but minor: (7), (10) and (11).
