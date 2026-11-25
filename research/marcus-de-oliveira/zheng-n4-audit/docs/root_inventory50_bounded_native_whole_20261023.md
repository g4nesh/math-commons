# Audit: root_inventory50_bounded_native_whole_20261023 (primary proof of paper case T_51)

Evidence level: EXACT for the contact set, the eight-generator structure, the tie list and the existence of earlier-stage pullbacks for every possible new tie (`replay/root_inventory50_bounded_native_whole_20261023__check.py`, 15 checks, all pass); PROOF-READ for the assembly; CITED for the conditional reduction (audited in this batch), the seven-element cases B, C, D, the six-element cases used for ties 0, 7, 16, the bounded-negative-part lemma, Appendices B and C and Theorem 1.2.

PDF: `proofs/manin_transfer_20261021/root_inventory50_bounded_native_whole_20261023.pdf`, 1 page.

## Deviations and limits

- One-page assembly; the analytic content is in `root_inventory50_bounded_native_conditional_20261023` (separate report, same batch).
- The document does not name the six-element cases used for new ties 0, 7 and 16 ("a strictly prior entire family"). The auditor resolved them by search over the symmetry group (below).
- `blind_moc4_coreD_uniform_20261023` (proof of D) was not read; CITED.

## A. Statements

**Theorem 1 (p. 1).** For every x in X with P c(x) >= 0 whose contacts contain T = {1,3,5,10,11,12} (literally paper T_51; non-strict, extra contacts allowed), Phi_{c(x)}(O) >= 0 for all real orthogonal O, hence for unitary U by Theorem 1.2. "All admissible parameter walls are included." Degenerate spectra are covered through the conditional document's appeal to the "global theorems" (Appendices B, C), CITED. This is (*) for T_51.

## B. Dependencies

| dependency | kind | resolution | status |
|---|---|---|---|
| root_inventory50_bounded_native_conditional_20261023 | branch | same batch: confirmed-modulo-citations | CITED |
| D = {0,1,2,7,9,16,23} | seven-contact | Table 5 T_12^(7), stage 1; tie 20 | CITED (blind_moc4_coreD_uniform) |
| B = {0,1,2,7,9,10,16} | seven-contact | T_8^(7), stage 1; ties 6, 13 | CITED |
| C = {0,1,2,7,9,16,17} | seven-contact | T_11^(7), stage 1; tie 14 | CITED |
| ties 0 and 16 | paper cases | earlier-stage sets mapping into T+{j}: T_9, T_11, T_17, T_38, T_65, T_76 (and Table 5 T_3^(7)) | EXACT existence; CITED |
| tie 7 | paper cases | T_42, T_65, T_71 | EXACT existence; CITED |
| bounded-negative-part lemma, App. B, App. C, Thm 1.2 | common / core | | CITED |

T_51 is stage 7 in Table 6; every pullback found is from stages 1-6.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | contact set is paper T_51 | identity | EXACT | |
| 2 | eight exclusive nonnegative coefficients; one non-minor, d = g4 | identity | EXACT | generators(T) = f0,f1,f2,f8,f10,f12,f30,f59; f30 is the only generator nonzero at 4 |
| 3 | only W = 0 needs a second bound; two branches there; H4 + H14 holomorphic on both | identity / analytic | EXACT (in conditional script) + PROOF-READ | |
| 4 | every ordinary new tie maps to a prior family or B, C, D | case split | EXACT | ties {0,6,7,13,14,16,20}; ties with f30 > 0 force d = 0 |
| 5 | inverse explicit, nonconstant coordinate, last germ nonconstant, Rouche restoration | analytic | PROOF-READ (conditional report) | |
| 6 | D proved elsewhere discharges the hypothesis | assembly | PROOF-READ / CITED | |

## D. Falsification

See the conditional report (stratum points and exact random instances of L_O >= d f30(O)).

## E. Verdict

**confirmed-modulo-citations.** Presentation remark: the pullbacks for ties 0, 7, 16 are not named; they exist and precede T_51 in Table 6 (auditor's search).
