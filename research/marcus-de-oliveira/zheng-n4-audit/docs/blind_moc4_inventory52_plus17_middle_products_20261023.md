Evidence level: EXACT (all twelve chart tables, term counts, signs and markers, rebuilt by own code). Real-branch input CITED from the low-products note, audited in this batch: confirmed-modulo-citations.

# blind_moc4_inventory52_plus17_middle_products_20261023 (2 pp.)

Replay: `replay/blind_moc4_inventory52_plus17_middle_products_20261023__check.py`. It imports the low-products replay. Output: `..._check_output.txt`, all PASS.

## Deviations and limits

- G and the octic F* come from the supplement data, validated as in the low-products report.
- The four-coefficient transport (7,21;9;35) was checked in the whole_assembly replay.

## A. Statement

**Theorem 1.** On the realizable strictly positive T7 face, assume a <= b <= c and ac <= d^2 < bc. Then A >= d(b+d)(c+d)/(bc - d^2). Hence Phi >= 0, through the four-coefficient atom on (f7,f21;f9;f35) plus nonnegative minors. All spectra are covered, as in the low-products note.

## B. Dependencies

- Real-branch criterion: if gamma > 0, R beta - N alpha > 0 and R F* - N alpha < 0, then the point is not realizable. This is exactly what the low-products note proves (steps 12-14 there), independent of the chart. The use matches.
- The four-coefficient atom (so4-atoms.md §3.4); the orbit was checked EXACT.

## C. Steps

| # | claim | status |
|---|---|---|
| 1 | the chart 0 < X <= Y <= 1 < Z covers the chamber (bc = Z^2, ac = Y^2, ab = X^2) | PROOF-READ |
| 2 | A = UF/(V(1+F)), F >= 0, covers 0 <= A < U/V, with V = Z^2 - 1 > 0 | PROOF-READ |
| 3 | each controller has degree <= 3 in A (gamma has degree 2) | EXACT |
| 4 | formula (3) for V^3(1+F)^3 P(UF/(V(1+F))) | EXACT (built that way and expanded) |
| 5 | term counts [4473,4824,5131,4374], [5220,5913,6435,5687], [5179,6624,7299,7514]; signs +, +, -; markers 4096 x^18 z^3, 2048 x^18 z^3, -384 x^18 z^3 C | EXACT |
| 6 | strictness on the walls W0 = 0, y = 0: the markers involve only x, z, C > 0 | PROOF-READ |
| 7 | failure of (1) is excluded, so (1) holds; the equality endpoint is covered by the atom | PROOF-READ + CITED atom |

## D. Falsification

See the whole_assembly report.

## E. Verdict

**confirmed-modulo-citations.** Cited: the real-branch criterion of the low-products note, Prop. C.6 and the four-coefficient atom.
