# root_moc4_T6_uniform_20261022 (primary proof of paper case T_7)

Evidence level: EXACT for every identity, plus the forward direction of the chart (38 checks, 0 failures). PROOF-READ for the sign argument. CITED for the circle theorem and the shared-zero / repeated-spectrum theorems.

PDF: proofs/manin_transfer_20261021/root_moc4_T6_uniform_20261022.pdf (3 pages). Replay: `replay/root_moc4_T6_uniform_20261022__check.py` and its output.

## Deviations and limits

- The chart identities are exact rational-function identities (python-flint). They use the same Re/Im reduction as the T_6 replay, because kappa is i mu over a real-coefficient expression.
- I supplied the forward direction of the chart myself. Contacts 0, 1, 2, 4 give h and t real. Contact 8 is equivalent to s being real. With w = v + s(v-r)/r, (F0-F7)/(F0-F2) = -(sX + 1/(k-1)) with X = (k/r - 1/v)/(k-1) - 1, so contact 7 is equivalent to l being real. All of this was checked exactly.
- No numerical search on the main statement was run (budget).

## A. Statements

- **Theorem 1 (p. 1).** For x in X with g >= 0 and g = 0 on T = {0,1,2,4,7,8} (literally paper T_7; checked), Phi_{c(x)}(O) >= 0 for all O. Contacts are non-strict and all spectra are covered. Distinct labels with nonreal v give a minor decomposition. Repeated or shared labels go to the shared-zero and repeated-spectrum theorems (paper Theorem B.1 and Lemma B.2). Real v puts the whole spectrum A on R ∪ {∞}, which goes to the circle theorem (paper Appendix C, Theorem C.1).

## B. Dependencies

| ref | kind | used for | status |
|---|---|---|---|
| Paper Theorem B.1 and Lemma B.2 | core-doc | shared and repeated labels | CITED (project-audited) |
| Paper Appendix C, Theorem C.1 (one whole spectrum on a generalized circle) | core-doc | v real | CITED |
| Paper Lemma 4.2 / Prop. 4.1 | core-doc | generator list | lib EXACT |
| Paper Theorem 1.2 | core-doc | unitary transfer | CITED |

## C. Steps

| # | claim | status |
|---|---|---|
| 1 | generators {0,1,2,3,4,5,6,11,21,22,23,26}; the nonminor ones are 22 and 26, both positive at permutation 3 | EXACT |
| 2 | dimension 8 with the stated coordinates; identity (2) for g19 | EXACT |
| 3 | 10 minor generators of rank 8; exactly 14 dual rays; exactly 2 negative on a source generator, namely A and B of (1) | EXACT |
| 4 | contacts force h and t real; kappa = i mu/((z-c)(v-r)); contacts 7 and 8 give l and s real (auditor's identities) | EXACT |
| 5 | the stated ratio formula for (1-c)(z-r)/[r(z-c)] | EXACT |
| 6 | converse: chart (3) makes all six contact gaps vanish identically | EXACT |
| 7 | excluded parameter values, k ≠ 0, nonvanishing denominators | PROOF-READ |
| 8 | identities (4) and (5) for g3, g6, g10, g14 = e, g18 and g9 | EXACT |
| 9 | (6) and (7): both expressions for A and for B | EXACT |
| 10 | g3 = 0 gives a minor decomposition; e = 0 gives A, B >= 0 | PROOF-READ (+ EXACT sign of f22 and f26 at 3) |
| 11 | ratio identities G/g10, g6/e, g10/g6 and g18/e, which give the signs (8) | EXACT + PROOF-READ |
| 12 | interpolation identity, alpha/D = 1 - h(1 + W/D), and k - l - 1 = -D - W | EXACT |
| 13 | contradiction argument for m >= 0 (the k/D < 0 sub-cases), then A, B >= 0 | PROOF-READ |

## D. Falsification

Every identity was replayed exactly, and none failed. The inequality chain consists of strict sign consequences, each read line by line. No numerical search was run.

## E. Verdict

**confirmed-modulo-citations.** The citations are Appendix B (B.1, B.2), Appendix C (Theorem C.1) and Theorem 1.2. Corollary B.3 is not used.
