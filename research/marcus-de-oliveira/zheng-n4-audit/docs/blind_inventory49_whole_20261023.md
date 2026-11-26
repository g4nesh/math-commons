# Audit: blind_inventory49_whole_20261023 (primary proof of paper case T_50)

Evidence level: EXACT for the contact set, the representative action, the ten-source decomposition, the bound (1) and the seven literal tie actions (`replay/blind_inventory49_whole_20261023__check.py`, 25 checks, all pass); PROOF-READ for the assembly logic; CITED for the conditional reduction (audited separately in this batch: confirmed-modulo-citations), the seven-element cases F, G, A, B, the six-element cases T_10, T_43, T_66, the bounded-negative-part lemma, Appendices B and C and Theorem 1.2 of the paper.

PDF: `proofs/manin_transfer_20261021/blind_inventory49_whole_20261023.pdf`, 3 pages.

## Deviations and limits

- This document is a summary plus assembly; the analytic content lives in `blind_inventory49_bounded_native_conditional_20261023`, audited in its own report (same batch). Its identities were replayed there.
- "The wrapper certificate" and "the exact certificate" referred to in the document were not used.
- Stage order was checked against Table 6 of the paper as transcribed from the text extract (p. 63).

## A. Statements

**Theorem 1 ('whole', p. 1).** For every x in X with P c(x) >= 0 and (P c(x))_pi = 0 for all pi in T = {1,3,5,6,10,11} (contacts non-strict; extra contacts, zero coefficients and all spectra allowed), Phi_{c(x)}(O) >= 0 for every O in SO(4) (stated for real orthogonal O; O(4) is equivalent by Section 4.2), and then for unitary U via the paper's Theorem 1.2.

Check against the paper: T is literally paper T_50 (Appendix E) (EXACT); the statement is exactly (*) for T_50, non-strict, all O, all spectra (repeated eigenvalues via Appendix C, a_i + b_j = 0 via Appendix B, both cited, p. 2 of the document: "Repeated spectral points and cross-spectrum common points are covered by the retained one-circle and shared-zero theorems").

Restated sub-claims: the ten-source nonnegative decomposition (EXACT, as in the conditional document); bound (1) -L_O <= M min(g4,g15) with M = max(0,-Phi_{f30}(O)); the tie table; the inverse parametrization summary (H4, H15, H0, C(W,Z), W = 0 branch); the wall restoration.

## B. Dependencies

| dependency | kind | resolution | status |
|---|---|---|---|
| blind_inventory49_bounded_native_conditional_20261023 | branch | audited in this batch: confirmed-modulo-citations | CITED (same batch) |
| F = {0,1,2,7,9,10,14} | seven-contact | Table 5 T_7^(7), stage 1 | CITED |
| G = {0,1,2,7,9,10,23} | seven-contact | Table 5 T_10^(7), stage 1 | CITED |
| A = {0,1,2,7,9,10,17}, B = {0,1,2,7,9,10,16} | seven-contact | T_9^(7), T_8^(7), stage 1 | CITED |
| "T 9 : {0,1,3,5,7,8}" | paper case | T_10 (stage 2); action maps it into T+{0} | EXACT action, CITED result |
| "T 42 : {1,2,4,6,7,8}" | paper case | T_43 (stage 1); into T+{7} | EXACT action, CITED |
| "T 65 : {2,4,5,6,9,12}" | paper case | T_66 (stage 7); into T+{16} | EXACT action, CITED |
| bounded-negative-part harmonic lemma; one-circle; shared-zero Birkhoff | common lemmas | final_referee_bounded_bad_harmonic_lemma; paper App. C, App. B | CITED |
| simultaneous two-real-response theorem | core | paper Theorem 1.2 | CITED |

Every cited contact set is earlier than T_50 (stage 8) in Table 6 (EXACT check of the transcribed table).

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | T is the literal tuple; representative {1,3,4,6,7,8} -> T by the stated action | identity | EXACT | lib group action with the document's convention |
| 2 | contact rank 6, kernel dimension 8, basis [f0,f1,f2,f3,f8,f10,f27,f59] | identity | EXACT | |
| 3 | intervals meet; ten-source nonnegative decomposition; relations f4, f30 | identity + inequality | EXACT / PROOF-READ | replayed in the conditional script |
| 4 | every source except 30 is a minor; d = 0 gives Phi >= 0 | identity | EXACT | F_IS_MINOR |
| 5 | bound (1) | inequality | PROOF-READ + 30 EXACT instances | f30(4) = f30(15) = 1 |
| 6 | extra tie with f30 > 0 forces d = 0; other ties exactly {0,7,12,13,14,16,20} | case split | EXACT | |
| 7 | the seven literal actions (table p. 2) | identity | EXACT | images computed on all 24 entries; each inside T+{j} |
| 8 | multiple extra ties are included | PROOF-READ | containment in T+{j} implies containment in any larger set |
| 9 | degenerate spectra via one-circle and shared-zero theorems | citation | CITED | |
| 10 | inverse parametrization, H4, H15, H0, C, W = 0 branch, (12) | identity | EXACT | replayed in the conditional script |
| 11 | nonconstant coordinate, harmonic exclusion, Rouche restoration | analytic | PROOF-READ (in the conditional report; one clarification of the puncture set supplied by the auditor) | |
| 12 | F and G discharge the two hypotheses of the conditional theorem | assembly | PROOF-READ | the conditional theorem's hypotheses are exactly the F and G inequalities |
| 13 | unitary assertion | citation | CITED | paper Thm 1.2 |

## D. Falsification

Bound (1) on 30 exact random instances (pass). The T_50 stratum test is reported in the conditional document's report (6 points, worst min Phi = -1.0e-15, no violation).

## E. Verdict

**confirmed-modulo-citations.** Citations: the conditional reduction (confirmed-modulo-citations in this batch), F, G, A, B (Table 5, stage 1), T_10, T_43, T_66, the bounded-negative-part lemma, Appendices B and C, Theorem 1.2.

Presentation remark: the six-element inputs are labelled by inventory index ("T 9", "T 42", "T 65" are paper T_10, T_43, T_66); the explicit tuples make them unambiguous.
