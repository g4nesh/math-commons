# root_moc4_T7_uniform_20261022 (primary proof of paper case T_8)

Evidence level: EXACT for the finite cone algebra, every chart and trace identity, the discriminants, the denominator divisibility and the exhaustiveness of the source curve (41 checks + 3 INFO, 0 failures). PROOF-READ for the harmonic minimum-principle argument. CITED for the circle, collision and seven-contact inputs.

PDF: proofs/manin_transfer_20261021/root_moc4_T7_uniform_20261022.pdf (4 pages). Replay: `replay/root_moc4_T7_uniform_20261022__check.py` and its output.

## Deviations and limits

- **Sign convention.** The document's products give the bracket of A0 = ∞ with every C_j the value -1 (F0/D = -x(v-1)). The paper's bracket gives +1. The replay uses the document's sign. A global sign only flips mu, and every ratio identity is unaffected. With the paper's sign, R1 - R3 = -E, and the document's mu = -1 becomes mu = +1.
- The chart and trace identities are exact rational-function identities (python-flint) in (E, conj E, r, x, s, t). I added the exhaustiveness of (4), which the document only asserts. For any v, the E- and s-definitions are solved by the stated w and z. The t-definition minus t then equals [(v-r)P0(E) + t{(r-1)E + s(x-1)}] / (a quadratic). This was checked exactly.
- The harmonic argument was read, not machine-checked. The seven-contact cases for j = 4, 10, 18, 19, 22 were resolved to paper cases but not audited. No numerical search on the main statement was run.

## A. Statements

- **Theorem 1 (p. 1).** For x in X with g >= 0 and g = 0 on T = {0,2,5,6,7,8} (literally paper T_8; checked), Phi_{c(x)}(O) >= 0 for all O. Contacts are non-strict, and every collision layer is included.
- The distinct-label chart A = (∞, v, xv, z), C = (r, 1, 0, w) has the source curve (4) in one complex variable E.
- Real v is sent to the "three-A/three-C common-circle theorem". Repeated and shared labels are sent to the repeated-spectrum and shared-value inputs.

## B. Dependencies

| ref | kind | used for | status |
|---|---|---|---|
| Paper Theorem B.1, Lemma B.2 | core-doc | shared and repeated labels | CITED (project-audited) |
| "entire three-A/three-C common-circle theorem" | external / other-doc | v real (∞, v, xv and r, 1, 0 all on R) | CITED. Not resolved to a precise paper statement: Appendix C treats one whole spectrum on a circle, and the reference meant here should be named |
| seven-contact j = 4, 10, 18, 19, 22 | paper-case | boundary of Omega | resolved: some G-image of T_2/T_4/T_7 (j=4), T_3/T_6/T_10/T_11 (j=10), T_11/T_16 (j=18), T_2/T_13 (j=19) and T_3/T_12 (j=22) lies in T ∪ {j}. All are consistent with the document's closed list (historical T0-T6, T9, T10 = paper T_1, T_2, T_6, T_3, T_4, T_5, T_7, T_10, T_11). CITED |
| root_moc4_T7_boundary16 / boundary23 | branch | j = 16, 23 | audited in this batch |
| SOC atom (F03) on (f2, f39; f46) | common-lemma | j = 16, 23 | CITED. The G-orbit relation to (f1, f9; f22) was checked EXACT |
| Paper Section 3 (each nonzero x in X is spectrally realizable; X closed) | core-doc | full-height continuation | CITED |
| Theorem 1.2 | core-doc | real-to-unitary | CITED |

## C. Steps

| # | claim | status |
|---|---|---|
| 1 | generators {0,1,2,3,4,5,6,39,45,46}, the only nonminor one is 46; d values on the generators; the two replacements; d <= 0 gives a minor decomposition | EXACT |
| 2 | (1) holds on all of W. On N_T ∩ {d >= 0}, all eight coefficients are >= 0 (checked on the generators of that cone, including f46+f0 and f46+f1) | EXACT |
| 3 | (F2-F0)/D, (F6-F0)/D, (F8-F0)/D; x, r and mu real; R5 = -sx and R7 = x(r-1) - t | EXACT (+ PROOF-READ for the reality deductions) |
| 4 | (2): R1 - R3 = E, R9, the x-combination identity, and R16 - (1-s)R12 real. Consequences g9 = -rb, a = -xd, g16 = (1-s)g12 | EXACT |
| 5 | signs (3): x < 0, r < 0, 0 < s <= 1; b > 0 | PROOF-READ |
| 6 | (4) inverts the chart and is exhaustive (auditor's check); P0(E) ≠ 0 for nonreal E | EXACT + PROOF-READ |
| 7 | trace identities (5) and (6) | EXACT |
| 8 | (7): L > 0 and t > 0 | PROOF-READ (from the exact (5)) |
| 9 | the indices j with f46(j) = 0 are exactly {4,10,16,18,19,22,23}; for j with f46(j) > 0 the nonminor coefficient is forced to 0 | EXACT |
| 10 | the reduced denominators of all R_j divide D·L3·L10. disc_t(disc_E D) = -16r^2 s(s-1)(x-1), with leading coefficient 1 in t, and (r-t)^2 at s = 1. L10 is a nonzero constant when rx = 1 | EXACT |
| 11 | v has denominator P0, with disc P0 = ((1-s)x-1)^2 + 2t(1-(1+s)x) + t^2 > 0; F0/D = -x(v-1) | EXACT |
| 12 | dR has denominator D only and is finite at ∞; R3·L3 is a polynomial; D(E3) = rs(x-1)L/(r-1)^2 ≠ 0 | EXACT |
| 13 | continuation: the vector is holomorphic on the upper half-plane, lies in X (closed) and is nonzero (R1 - R3 = E) | PROOF-READ |
| 14 | (8): L_O >= -C_O d on supported points with d > 0; d -> 0 at every real point and at ∞ (via d <= g3 at the roots of D) | PROOF-READ |
| 15 | minimum principle on each component Omega | PROOF-READ |

## D. Falsification

Every identity was replayed exactly, and none failed. No counterexample search was run on the main statement.

## E. Verdict

**confirmed-modulo-citations.** The analytic argument is sound as written. Every pole is listed, including the roots of D, where the bound d <= g3 handles tangential approach. Infinity is covered, and the curve stays in X through the chart exceptions. The citations are Appendix B, the three-A/three-C circle theorem (precise reference to be named), the five seven-contact paper cases above, the SOC atom and Theorem 1.2. Corollary B.3 is not used.
