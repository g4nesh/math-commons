# Audit of blind_moc4_T38_uniform_20261022 (primary proof of paper case T_39)

Evidence level: EXACT for the source cone, the three generator identities, the C/Q formulas of both branches, the sign of C on the g70 hard region and the group map of the T24 core; PROOF-READ for the branch reduction, domain exhaustiveness, Bernstein/SOC peeling logic and the closure arguments; NUMERIC for G40 and G70 non-vanishing below the stated endpoints; UNCHECKED for the 16 + 24 + 60 coefficient tables, the 799 Sturm groups and G40, G70 themselves; CITED for the capacities, the T27/T28 allocation lemmas and Theorem 1.2.

Auditor label C15-T38-T40. Document: proofs/manin_transfer_20261021/blind_moc4_T38_uniform_20261022.pdf, 6 pages. Replay: `replay/blind_moc4_T38_uniform_20261022__check.py` (0 failed checks).

## Deviations and limits

- About 15 minutes. I did not rebuild G40 (1455 terms) or G70 (1113 terms) as polynomials, nor any cleared table, Sturm sequence, KKT chart table or the identity k0 = C G6. All are UNCHECKED.
- As a substitute, I exactly checked that -(Q^2/4) det Hess at h = -C/Q is a rational square at 4 points per branch. I also ran exact Sturm root counts of G along random lines: G40 in j on [0, J] in domains 1 and 7, and G70 in a on [0, min(J, 1/s)). For G70, J was computed by numeric 1-D minimization.
- Domain 8 (ab >= 1, ac >= 1, allocation endpoint) and domains 2-6 were not tested numerically.

## A. Statements

Theorem 1 (p. 1). For every x in X with P c(x) >= 0, contact set containing T_39 = {2,3,5,7,9,12}, and every O in SO(4) (and through Theorem 1.2 every unitary), Phi_{c(x)}(O) >= 0. Contacts are non-strict, and zero coefficients are included. The opening sentence restricts the scope to "spectrally realizable radial necessity, with their existing complex-spectral limiting scope". The actual argument uses only x in X ⊂ Sing C, so it covers all spectra, including degenerate ones. The contact set equals paper T_39 literally (EXACT).

Internal claims: three generator identities; branch representations (1), (2) on 8-element bases; Q and C formulas; the capacity list for g40 transported from the T24 core and the (0,7,12,40) four-coefficient capacity; Le, Lp bounds for k0; Bernstein endpoint certificates; the H/chain resource capacity J for g70 and its KKT charts; SOC peeling G70 >= (1-as)^2 (1-a/J)^2.

## B. Dependencies

| kind | ref | used for | status |
|---|---|---|---|
| core-doc | T24 core (paper T_25 generators {0,1,2,6,53,54,55,59,67,70}), mapped (f59,f53,f55,f67,f54) -> (f5,f2,f67,f0,f40) | g40 capacities: old-H, parallel, coherent, four-coefficient E_ab | group map EXACT (g = ((0,3,1,2),(1,2,3,0),0)); capacities CITED |
| core-doc | four-coefficient capacity with tuple (0,7,12,40), P_ac | g40 | CITED |
| core-doc | "proved T27 allocation lemma" (paper T_28 doc) | 8 rational domains, allocation endpoint | CITED |
| core-doc | "proved T28 resource lemma" (paper T_29 doc), H/chain resource allocation | g70 capacity J and KKT charts | CITED; numeric J used in tests |
| core-doc | old SOC (2,5,70), old-H (12,2,70), (39,2,70), chains | g70 | CITED |
| common-lemma | Theorem 1.2 | U(4) | CITED |

Note: T_25, T_28, T_29 and T_39 are all in stage 1 of the paper's Table 6, which says proofs within a stage are independent. This document uses pointwise capacity lemmas proved in the T_25, T_28 and T_29 documents. These are not their case propositions, so there is no cycle, but the stage table is inaccurate at this level of detail (presentation).

## C. Steps

| # | claim | kind | status |
|---|---|---|---|
| 1 | generators of T_39 = {0,1,2,5,7,12,39,40,67,70}; non-minor 40, 70 | identity | EXACT |
| 2 | f40+f70 = f0+f5+f12+f39+f67; f1+f40 = f0+f7+f39; f7+f70 = f1+f5+f12+f67 | identity | EXACT |
| 3 | the reduction procedure yields (1) or (2) with nonnegative coefficients, or a minor sum; both bases have rank 8 | logic + identity | PROOF-READ + EXACT |
| 4 | g40: Q = d^2+de+(d+e)(p+q)-ej; C has 21 positive coefficients including d^3, so Q = 0 is excluded | identity | EXACT |
| 5 | g70: C = bc(a-U) - T0, Q = (a-U)M - dU - as; -MC identity; MT0 - bc has nonnegative coefficients containing b and c; M = 0 gives C = -(1+p+r)(1+q) | identity | EXACT |
| 6 | det Hess = -4 G^2 / Q^2 (after Q^4 clearing) | identity | UNCHECKED as a polynomial; EXACT square property at 8 points |
| 7 | g40 capacity list and the core map | citation | group map EXACT; capacities CITED |
| 8 | k0 = C G6, bilinear G6; four coefficient functions nonnegative after the Le, Lp shifts (16 tables) | certificate | UNCHECKED |
| 9 | domain cover for g40: ab<1 split as {a,b<=1}, {a>=1}, {b>=1} and by ac<1/ac>1; ab>1 with ac<1; ab,ac>=1 | case split | PROOF-READ (exhaustive; I also checked that the a>=1 endpoint is <= 1+1/e) |
| 10 | seven direct domains (24 tables) and 799 Sturm groups; allocation domain (3 tables) | certificate | UNCHECKED |
| 11 | Bernstein identity G40 = (1-j/J)^2 k0 + ... gives G40 > 0 on 0 <= j < J | logic | PROOF-READ |
| 12 | closed boundaries and zero coefficients via continuity of non-strict coefficient inequalities and strict k0 | analytic | PROOF-READ (non-strict inequalities pass to limits; k0 >= C g00 (1-Le e)(1-Lp p) > 0 stays strict) |
| 13 | g70 capacity J (H/chain allocation) and its derivation | citation + logic | CITED; minimization algebra (delta_i = min(p/2, beta_1) etc.) PROOF-READ |
| 14 | SOC peeling tables L0, B1, B2 (60 tables) and the identity (1+aS)^2 G70(a, S/(1+aS)) = L0 + S B1 + S^2 B2 | certificate | UNCHECKED |
| 15 | zero resources via eps-perturbation and closedness of the SAFE set | analytic | PROOF-READ |

## D. Falsification

- NUMERIC (exact Sturm on random lines, seed 39): G40 has no root in j on [0, J] on 6 random domain-1 lines (J = 1+1/e) and 6 random domain-7 lines (J = 1/e - 1/E_ab). G70 has no root in a on [0, min(J, 1/s)) on 6 random lines.
- No stratum sampling was done for this document.

## E. Verdict

could-not-check. No gap found. The structural steps are EXACT, the logical skeleton is PROOF-READ, and the numerical tests of the decisive positivity claims pass. The certificate tables were not replayed. Presentation remarks: the Table 6 stage-independence wording, and the vague scope phrase on p. 1.
