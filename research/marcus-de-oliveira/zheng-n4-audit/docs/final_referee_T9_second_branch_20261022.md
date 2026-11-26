Evidence level: EXACT for the identities of Lemma 1, (2), (3), Lemma 2 (4), the conj(zeta) forms and the Theorem 1 algebra (20 exact checks in `replay/final_referee_T9_second_branch_20261022__check.py`). PROOF-READ for the sign arguments. NUMERIC for the conclusions.

# final_referee_T9_second_branch_20261022: region x<0, t>1 of {0,1,3,5,7,8} (paper T_10)

PDF: `proofs/manin_transfer_20261021/final_referee_T9_second_branch_20261022.pdf`, 3 pages. This is a branch document.

## Deviations and limits

- The text extract drops conjugation bars; I read the formulas from the page images.
- The chart comes from the T9 main document and is replayed there.
- The four-coefficient inequality (1) is CITED. Its source is unresolved; see the T9 main report.
- The note's closing remark says that the region x > 1, t < 0 "still requires its own capacity proof". That region is supplied by the first-branch note.

## A. Statements

Setting: T = {0,1,3,5,7,8}, eight distinct points, v nonreal, region x < 0, t > 1, r > 0.

- **Lemma 1 (p. 1).** A supporting source with d > 0 has u < 0.
- **Lemma 2 (p. 2).** With x = -X, t = T+1, u = -U, normalized g2 = 1, V = TH - t, C0 = GN + A0 V and S0 = XG + A0 H: d C0 = N(m S0 + GH).
- **Theorem 1 (p. 2).** If d > 0, g18 >= 0 and g22 >= 0, then a = g14, B = g9 and c - d = g19 - d are positive, and aB(c - d) > d^2(a + B + c + d). The cost is therefore safe by the four-coefficient inequality.

## B. Dependencies

- T9 main chart and decomposition: EXACT in the main check.
- Four-coefficient inequality on (f4, f6, f2, f28): CITED, source unresolved. NUMERIC: no violation, and the inequality is sharp.

## C. Steps

| # | claim | status |
|---|---|---|
| 1 | Lemma 1 identities: g10 N0 NY = M = -xh{...}; d\|Qv\|^2 = r^2 D2 + rtuL Nx; D2 = M/h - tu(1-x) N0; D2 + L NY = M/x; g9 = xy(r/Nx - tu/NY); g10 = g20 + d | EXACT |
| 2 | Lemma 1 signs: rho = -1 is a positive rescaling (y < 0); M > 0; h < 0 and x < 0 give D2 < 0 and D2 + L NY < 0; g9 >= 0 gives r >= r0 >= 0; the decomposition d\|Q\|^2/r = (r - r0) D2 + r0(D2 + L NY) < 0 (EXACT) contradicts d > 0 | PROOF-READ |
| 3 | (2) and (3): g14 = Tm - d, g19 = m + G, g9 = X(G+H), g18 = 1 - d, g22 = (Xm - XH - d)/A0 | EXACT |
| 4 | Lemma 2 (4), together with Q/(zeta(zeta+T)) = C0 + S0 conj(zeta) and K = -d - m conj(zeta) | EXACT |
| 5 | m >= H + d/X, so C0 > 0; derivation of (5) | EXACT / PROOF-READ |
| 6 | V <= 0 contradicts (5), so H > t/T and m > H > 1; then a > 0, c > d (since m > 1 >= d) and B > 0 | PROOF-READ |
| 7 | capacity difference = X(G+H) D0 - L0 d^2; T m(G+H) - L0 = G(Tm-1) + m(TH-t) > 0; D0 > 0; the lower bound via X >= d/(m-H) | EXACT / PROOF-READ |

## D. Falsification

- Region II with u < 0 (seed 7, 400000 draws): 71 supporting points with d > 0. The minimum capacity difference was 92, so no violation.
- Region x < 0, t > 1, u > 0 (same draws with u -> -u): 0 supporting points with d > 0, consistent with Lemma 1.

## E. Verdict

**confirmed-modulo-citations.** I found no gap. The only external input is the four-coefficient inequality, whose source is unresolved.
