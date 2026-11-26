# Audit: root_inventory51_bounded_native_whole_20261023 (primary proof of paper case T_52)

Evidence level: EXACT for the contact set and its stated representative, the nine-ray cone, the tie list and earlier-stage pullbacks for every possible new tie (`replay/root_inventory51_bounded_native_whole_20261023__check.py`, 16 checks, all pass); PROOF-READ for the assembly; CITED for the conditional reduction (audited in this batch), the seven-element cases A, C, G, the six-element cases for ties 0, 7, 16, the bounded-negative-part lemma, Appendices B and C and Theorem 1.2.

PDF: `proofs/manin_transfer_20261021/root_inventory51_bounded_native_whole_20261023.pdf`, 1 page.

## Deviations and limits

- One-page assembly; the analytic content is in `root_inventory51_bounded_native_conditional_20261023` (separate report).
- The six-element pullbacks for ties 0, 7, 16 are not named in either document; the auditor resolved them.
- `independent_shared_coreG_uniform_20261023` (proof of G) not read; CITED. The document's claim that it does not depend on cases 49 and 51 (inventory) was not checked beyond Table 6 (G is stage 1).

## A. Statements

**Theorem 1 (p. 1).** For every x in X with P c(x) >= 0 whose contacts contain T = {1,3,5,10,11,14} (literally paper T_52; non-strict), Phi_{c(x)}(O) >= 0 for all real orthogonal O, and the unitary statement via Theorem 1.2. The stated representative {3,4,5,7,9,12} lies in the orbit of T (EXACT). "Including boundary configurations under the stated hypotheses": degenerate spectra are handled in the conditional document by "retained full results" (Appendices B, C), CITED. This is (*) for T_52.

## B. Dependencies

| dependency | kind | resolution | status |
|---|---|---|---|
| root_inventory51_bounded_native_conditional_20261023 | branch | same batch | CITED |
| G = {0,1,2,7,9,10,23} | seven-contact | T_10^(7), stage 1; tie 20 | CITED |
| A = {0,1,2,7,9,10,17} | seven-contact | T_9^(7); ties 6, 13 | CITED |
| C = {0,1,2,7,9,16,17} | seven-contact | T_11^(7); tie 12 | CITED |
| ties 0, 16 | paper cases | T_10, T_11, T_12, T_24, T_39 (stages 1-5) | EXACT existence; CITED |
| tie 7 | paper cases | T_42, T_46 (stages 1-2) | EXACT existence; CITED |
| bounded-negative-part lemma, App. B, C, Thm 1.2 | | | CITED |

T_52 is stage 6; every input is earlier.

## C. Steps

| # | claim | kind | status |
|---|---|---|---|
| 1 | contact set = T_52; representative in orbit | identity | EXACT |
| 2 | sole non-minor coefficient dominated by g4 and g15 | inequality | EXACT (f30(4) = f30(15) = 1) + PROOF-READ |
| 3 | only W = 0 is a pole of H4; fibre quadratic with finite points; resultant keeps D14 nonzero there; H15 bounded | identity / analytic | EXACT (conditional script) + PROOF-READ |
| 4 | ordinary boundaries reduce to earlier cases or A, C, G | case split | EXACT |
| 5 | inverse, nonconstant last germ, Rouche restoration, harmonic lemma | analytic | PROOF-READ (conditional report) |

## D. Falsification

See the conditional report.

## E. Verdict

**confirmed-modulo-citations.** Presentation remark: pullbacks for ties 0, 7, 16 not named (exist; earlier stages).
