Evidence level: PROOF-READ for the case logic, the chart coverage and the exhaustion. EXACT (in the whole_assembly replay) for the three four-coefficient transports. UNCHECKED for every coefficient table of this note: the 9 L_ij tables, the 12 + 12 controller tables, the 27 three-capacity coefficients with factor removal, the box tables, and the L00 sub-certificate.

# blind_moc4_inventory52_plus17_fourth_products_20261023 (3 pp.)

## Deviations and limits

- Time budget: none of this note's polynomial certificates were replayed. My own pipeline for the low and middle notes (own N, R, Q; G and octic validated) could be extended to it, but that was not done.
- Replay file `replay/blind_moc4_inventory52_plus17_fourth_products_20261023__check.py` checks only the chart-coverage inequalities, symbolically.

## A. Statement

**Theorem 1.** On the realizable strictly positive T7 face, assume a <= b <= c and d^2 < ab. Then at least one of A >= capA, B >= capB, C >= capC holds, with the caps as in (1). The support inequality then follows from the four-coefficient atom.

## B. Dependencies

- The real-branch criterion (3) from the low-products note: confirmed-mod-cit here. The use is the same, and chart independent.
- G, the octic F*, Prop. C.6.
- The four-coefficient atom: EXACT orbit check in the whole_assembly replay.

## C. Steps

| # | claim | status |
|---|---|---|
| 1 | Case 1 (a <= 1): parametrization (4) covers 0 < a <= 1, ab > 1, c >= b | PROOF-READ (ab = 1 + u) |
| 2 | Case 1: A0 < Acap, since m - abc > 0 | PROOF-READ |
| 3 | Case 1 lower interval A <= A0 with B < Bcap: nine positive L_ij tables imply G > 0 | UNCHECKED |
| 4 | Case 1 upper interval: twelve controller tables imply non-realizable by (3) | UNCHECKED |
| 5 | Case 2 charts (6) for a > 1 | PROOF-READ (coverage verified by hand) |
| 6 | Case 2: 27 three-capacity coefficients positive after removal of positive factors | UNCHECKED |
| 7 | Case 3 box (7): coverage, m - a > 0 and m - abc > 0 | coverage PROOF-READ; tables UNCHECKED |
| 8 | Case 3 tables and the L00 sub-certificate under the third failure | UNCHECKED |
| 9 | exhaustion: a <= 1, then a >= 2, then 1 < a < 2 with b >= 2, then a, b < 2 with c >= 3, then the box | PROOF-READ |
| 10 | the logic "G > 0 on an interval means no root, hence not realizable" | PROOF-READ (G = 0 is a necessary condition) |

## D. Falsification

See the whole_assembly report (numeric sample of region 4).

## E. Verdict

**could-not-check.** The logic and the exhaustion read correctly. All the polynomial sign certificates the proof rests on are unchecked here.
