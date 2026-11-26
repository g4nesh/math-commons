Evidence level: PROOF-READ for Theorem 1's sign argument, with every algebraic rearrangement re-derived EXACT by sympy; CITED for the spectral chart inputs (coordinate formulas (1)-(2), the no-division identity (3), the values (4)-(6) of the cubic P, s != 0); NUMERIC sanity of the sign bookkeeping.

# blind_moc4_T48_c0_native_interlacing_20261023 (3 pp.): interlacing on the c = 0 seven-contact face

Auditor label C26-T49. Script: `replay/blind_moc4_T48_c0_native_interlacing_20261023__check.py` (0 EXACT failures).

## Deviations and limits

- The spectral chart (parameters r, s, u, h, J, W), the coordinate formulas (1)-(2), identity (3) and the cubic P(J) are defined in the "retained original-product chart" and in `independent_true_T48_c0_native_Heron`. Neither is in my batch, and this document does not derive them. I could not replay them, so they are CITED. Everything this document derives from them I checked.
- The scope is canonical eight-role sources on the c = 0 face (T_49 u {17} = Table 5 set #23, EXACT).

## A. Statements

- Theorem 1. For a realizable source with eight distinct spectral points, canonical form, c = 0 and q, p, j, e, d, b, a > 0, the chart parameter r is negative.
- Corollary 2. Every realizable canonical c = 0 source satisfies the support inequality. The proof uses allocation III low-j for j <= d, and for j > d this theorem, then Heron, then the upper-Heron chain.

## B. Dependencies

| dependency | status |
|---|---|
| (1)-(2) coordinate formulas, (3) mN(Rv+1) = tH, J a root of the real cubic P, values (4)-(6) of P, s = 0 excluded by eight-distinctness | CITED (chart of `independent_true_T48_c0_native_Heron` / the retained chart) |
| allocation III low-j (j <= d) | CITED, not in this batch |
| coordinate Heron theorem | CITED |
| upper-Heron chain (this batch) | EXACT (see its report); its hypothesis j >= d matches the branch j > d |
| collisions: Theorem C.1 (repetition) and Theorem B.1 (cross) | C.1 CITED; B.1 established |
| zero faces | embeddings EXACT (allocation-II report), proofs CITED |

## C. Steps of Theorem 1 (assume r > 0 and derive a contradiction)

| # | claim | status |
|---|---|---|
| 1 | from (1): t > 0 and Rv + 1 < 0 (since b + d > 0, q + d > 0, y = j > 0) | PROOF-READ |
| 2 | v = 0 is impossible; v < 0 gives R > 0 and l = n + r(1-v) > 0, so p < 0 by (2), a contradiction. Hence v > 0, R < -1/v, l > 0 | PROOF-READ; EXACT for l = n - r(v-1) |
| 3 | s > 0: then n > t, and b > 0 gives rt > n, so r > 1 and m > 0. From (3), r m N q/y = rm(t-n) - (r-1) l t < 0, contradicting q >= 0 | PROOF-READ; EXACT for the rearrangement (with H = -rm + (r-1)l) |
| 4 | s = 0 excluded (CITED); s < 0 gives u > 1, then e > 0 gives v > 1 and u > t/(v-1); rt > n > 0 | PROOF-READ |
| 5 | delta > 0, t/(v-1) - delta = (l - s)/(v-1) > 0, so -u lies left of -delta and -t/v | EXACT + PROOF-READ |
| 6 | P(-u) > 0 by the four-term decomposition of its last factor, each term positive | EXACT (identity) + PROOF-READ (signs) |
| 7 | v delta - t = r - l and v delta - r - t = -l, so P(-delta) and P(-t/v) have opposite signs, ordered as claimed for l < r and l > r | EXACT + PROOF-READ |
| 8 | two distinct real roots, degree <= 3 and P nonzero, so P has no nonreal root; this contradicts the nonreal J | PROOF-READ |
| 9 | wall l = r: perturb t only. The strict inequalities persist, nearby polynomials are real-rooted, the limit is nonzero (P(-u) > 0), and Hurwitz/Rouche gives no nonreal root in the limit | PROOF-READ. The argument concerns coefficients only, so it is valid without the perturbed object being realizable |

## D. Falsification

NUMERIC: 2000 random points of the open region given by (4)-(5) (seed 5). The sign pattern P(-u) > 0, P(left) < 0 < P(right) with -u leftmost, computed from the stated factored values, held at 2000/2000. This checks only the bookkeeping, not the CITED values of P.

## E. Verdict

**confirmed-modulo-citations.** The argument is correct given its CITED chart inputs, which this document does not derive and I could not replay. Presentation: the equation labels (4) and (5) are used twice, once for the scalar inequalities and once for the values P(-delta) and P(-t/v). Corollary 2's case j <= d rests on allocation III low-j, which is outside this batch.
