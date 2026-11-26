# Audit: root_moc4_T21_T22_uniform_20261022 (paper cases T_22 and T_23)

Evidence level: EXACT for every identity (generators, 7 pair and 14 absorption replacements, normal forms, 6 cubic identities, rectangles, swap factorizations); PROOF-READ for the elimination procedure and the real-part lemma; CITED for collisions (paper App. B) and the unitary reduction (Thm 1.2); NUMERIC falsification only.

PDF: proofs/manin_transfer_20261021/root_moc4_T21_T22_uniform_20261022.pdf (3 pages). Replay: `replay/root_moc4_T21_T22_uniform_20261022__check.py`.

## Deviations and limits

- Supplement labels T21, T22 are the paper's T_22 = {2,3,4,8,10,12} and T_23 = {2,3,5,9,11,12} (same sets, identity group element).
- The supplement's cubic F is identified with the paper's C (A.2); all six identities hold with the stated signs sigma, so the identification is exact.
- Numeric stratum points were almost all collision points (T_22: 0/12 eight-distinct; T_23: 1/12), so the eight-distinct argument was tested numerically only at one point.

## A. Statements

Theorem 1 (p. 1). For T in {T_22, T_23}, every x in X with P c(x) >= 0 vanishing on T (extra contacts allowed) has Phi_{c(x)}(O) >= 0 for all O in SO(4) (hence O(4), U(4)). Eight-distinct layer: c(x) is a nonnegative sum of minor generators. Collision layers (repeated eigenvalue, a_i + b_j = 0): by citation of App. B. All spectra covered.

## B. Dependencies

Same as for root_moc4_T20_uniform_20261022: X subset Sing C (paper A.1, replayed in lib), paper App. B (Thm B.1, Lemma B.2) for collisions, Thm 1.2 for the unitary step. No other contact case.

## C. Steps

| # | claim | kind | status |
|---|---|---|---|
| 1 | S21 = generators(T_22) = {0,1,2,6,10,53,54,55,56,67,69,70,72}, non-minor {54,56,69,70}; S22 = generators(T_23) = {0,1,2,7,9,12,39,40,41,67,68,72}, non-minor {40,68} | finite | EXACT |
| 2 | F2 F10 = F4 F8, F3 F11 = F5 F9 | identity | EXACT |
| 3 | F2-F4, F2-F8, F3-F5, F3-F9 are products of two within-spectrum and two cross label differences | identity | EXACT |
| 4 | common real height H = 0 for distinct labels, so c = Re x' | algebra | EXACT + PROOF-READ (as in T_21 audit) |
| 5 | C(Re x) = 0 on X | lemma | PROOF-READ (polarization) |
| 6 | 7 pair replacements (every pair of non-minor types, right sides minor) | identity | EXACT |
| 7 | 14 absorption replacements, right sides minor and free of the same type's resources | identity | EXACT |
| 8 | elimination terminates with at most one non-minor type and, if it survives, with all its resources exhausted | combinatorial | PROOF-READ: each pair step removes a type; each absorption removes the non-minor type or a resource, resources never regenerate |
| 9 | the surviving support is exactly the listed 8-tuple lambda for each of the 6 non-minor types | finite | EXACT |
| 10 | sigma C(sum r_i f_{lambda_i}) = P24 or P33, six identities | identity | EXACT |
| 11 | P24, P33 have all coefficients positive and contain d^3 with coefficient 1, so P >= d^3 > 0 | identity | EXACT |
| 12 | contradiction with step 5, so eight-distinct costs are minor combinations | logic | PROOF-READ |
| 13 | collisions; O(4), U(4) | citation | CITED |

## D. Falsification

Library stratum points, 12 seeds per case (pi* = 2): worst min Phi over SO(4) is -5.9e-15 (T_22) and -2.7e-10 (T_23), 8 restarts each. Only one point (T_23, seed 9, label separation 0.068) had distinct labels; it had min Phi 3.0e-13 >= 0. No violation.

## E. Verdict

confirmed-modulo-citations (App. B, Thm 1.2, A.1). No gap. Presentation notes: the text's sentence "The separately nonnegative for the support functional collision layers ..." is garbled; the proof gives no pi* choice explicitly ("gaps anchored at that contact") but any contact in T works since all contacts share the height H.
