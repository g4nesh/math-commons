Evidence level: EXACT for every chart identity, the boundary relations, (4), (5) and the capacity identity; PROOF-READ for the sign chain (2)-(3); CITED for the SOC lemma (f1,f9,f18) (NUMERIC support), T_4 and Appendix B.

# independent_moc4_T8_double_boundary_20261022 (3 pages)

## Deviations and limits

- Replay: `replay/independent_moc4_T8_double_boundary_20261022__check.py`, 20 EXACT checks plus 1 NUMERIC check. All pass. The companion script was not used.
- The SOC lemma for (f1, f9, f18) is cited. I tested it numerically only: min over SO(4) of Phi_f18 + 2 sqrt(Phi_f1 Phi_f9) was ≈ 7e-6 >= 0 on 20000 random SO(4) matrices (seed 7). This is the condition for the lemma to hold, and it appears sharp.
- The document is logically redundant for T_9. Theorems 2 and 3 of the uniform document each cover points that also have the other contact.

## A. Statement

**Theorem 1.** Every x in X with P c(x) >= 0 and contacts containing {0,1,2,3,11,12,18,20} has Phi_{c(x)}(O) >= 0 for all O in SO(4). Repeated points and cross-collisions are covered by citing Appendix B. Real v: g4 = 0, so the coefficient of f18 vanishes and the point is minor (the document cites an "earlier circle theorem", but the minor argument is immediate). The theorem covers all spectra modulo these citations.

## B. Dependencies

| ref | kind | status |
|---|---|---|
| SOC (f1,f9,f18): a f18 + K f1 + C f9 is safe when CK >= a^2 | common-lemma | CITED, NUMERIC support |
| T_4 via T_9∪{18,20,7} and T_9∪{18,20,23} (b = 0 or g23 = 0) | paper-case | containment EXACT (inside T_9∪{18,7} and T_9∪{18,23}), case CITED (Stage 1) |
| Appendix B | common-lemma | CITED |
| "three-plus-three circle theorem" for real v | other-doc | not needed: g4 = 0 gives a minor combination directly |

## C. Steps

| # | claim | status |
|---|---|---|
| 1 | profile identity g = g23 f1 + g21 f2 + g20 f7 + g9 f8 + g7 f9 + g6 f12 + g18 f15 + (g4-g18) f18 on N_T | EXACT |
| 2 | inverse z, w from H12 = 1-t-r and H18 = s | EXACT |
| 3 | g4 = ρx/N, g7 = ρB/O, g13 = ρr(x-1)/O, g23 = ρrtxMH/Dp; M g4 + s g7 - Q g20 - M g11 = 0; (s-P) g7 - M g23 - D g20 = 0 | EXACT |
| 4 | a/b = -s/M, g23/b = (s-P)/M; M ≠ 0 | EXACT + PROOF-READ |
| 5 | g14 = g9 + a (forms); b = g9 - x g14 (in fact g7 - g9 + x g14 ≡ 0 on the chart); g9/b = (s-1)/M, g14/b = -1/M | EXACT |
| 6 | sign chain (2): M < 0, 0 < s <= 1, s < P, x < 1, t < 1 | PROOF-READ (each inequality checked by hand) |
| 7 | g21 = g13 - a - b ⇒ g13 > 0. The signs of r, B and H follow. Cases x > 0 are excluded. g21/b = (x-1)(r-P)(s-P)/(MB) (EXACT) ⇒ r >= P; t > 0 is excluded ⇒ (3) | EXACT + PROOF-READ |
| 8 | identities (4), (5) | EXACT (polynomial identities in v, vbar) |
| 9 | J > 0 (since r+t >= 1, s-P < 0, t < 0), and J g6/b = P{C + s(P-r)} ⇒ C >= 0 | PROOF-READ |
| 10 | capacity (b g23 - a^2)/b^2 = C/M^2 >= 0 | EXACT |
| 11 | SOC closes | CITED |

## D. Falsification

The only intermediate inequality stronger than the conclusion is b g23 >= a^2. Its derivation is exact plus a short sign chain. The SOC lemma was tested numerically (above). No violation was found.

## E. Verdict

**confirmed-modulo-citations.** Citations: the SOC lemma for (f1, f9, f18), paper T_4, and Appendix B. No gaps found beyond presentation: "the checker supplies explicit row and column permutations" is replaced by the exact containment checks.
