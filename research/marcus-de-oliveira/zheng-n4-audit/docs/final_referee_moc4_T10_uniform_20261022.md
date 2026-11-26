Evidence level: EXACT for the contact set and group element, generators, A and B values, the seven minor replacements, both branch decompositions and the capacity-orbit maps (22 exact checks in `replay/final_referee_moc4_T10_uniform_20261022__check.py`). PROOF-READ for the case split and the minor-case elimination. CITED for both branch proofs and the global inputs.

# final_referee_moc4_T10_uniform_20261022: six contacts {0,1,2,4,7,9} (paper T_11)

PDF: `proofs/manin_transfer_20261021/final_referee_moc4_T10_uniform_20261022.pdf`, 2 pages. This is the primary proof of paper case T_11.

## Deviations and limits

- The B < 0 branch is proved in `root_moc4_T10_B_branch_20261022`, which is not in this batch. It is CITED and was not audited.
- The A < 0 branch is audited in this batch (`final_referee_T10_A_branch_20261022`).
- The SOC and four-coefficient inequalities are CITED, source unresolved. Numeric evidence is in the T9 main report.

## A. Statements

**Theorem 1 (p. 1).** For every x in X (all spectra) whose gap vector satisfies g >= 0 and g_j = 0 for j in T = {0,1,2,4,7,9} (extra contacts allowed), Phi_{c(x)}(O) >= 0 for every O. T is not literally T_11. The element g = (rho, rho', eps) = (id, (0,1,3,2), 0) of the group maps T onto the paper's T_11 = {0,1,3,5,6,8} (EXACT), so the theorem is (*) for T_11 via the paper's symmetry (eq. (4.7), CITED). Scope: contacts are non-strict and every O is covered. Collisions, real circles and b = 0 are delegated to the branch documents and the global theorems.

## B. Dependencies

| dependency | status |
|---|---|
| Paper Lemma 4.2, Prop. 4.1 and (4.2) | CITED (lib-replayed facets) |
| Group covariance (4.7) | CITED (lib-tested) |
| final_referee_T10_A_branch_20261022 (A < 0) | audited in this batch: confirmed-modulo-citations |
| root_moc4_T10_B_branch_20261022 (B < 0), needing ab(c-d) > d^2(a+b+c+d) on (f4, f7, f3, f22) | CITED, not audited |
| SOC on (f2, f12, f26) | CITED (source unresolved). EXACT: it is the image of (f1, f11, f28) under ((0,1,3,2), (0,3,1,2), 0), the same orbit as the T9 SOC |
| four-coefficient inequality on (f4, f7, f3, f22) | CITED (source unresolved). EXACT: it is the image of (f4, f6, f2, f28) with roles (a, B, c, d) under (id, (0,1,3,2), 0) |
| Appendix B and Theorem C.1 (collisions and real v) | CITED |

## C. Steps

| # | claim | status |
|---|---|---|
| 1 | generators(T) = {0,1,2,3,4,5,7,11,12,21,22,26}; the non-minor ones are f22 and f26 | EXACT |
| 2 | A = g16 - g3 is 1 on f1, f4, f7, f22, -1 on f26, 0 elsewhere; B is 1 on f0, f2, f12, f26, -1 on f22, 0 elsewhere; so A + B >= 0 | EXACT |
| 3 | the seven minor replacements | EXACT |
| 4 | if A, B >= 0 the cost is in the minor cone: cancel min(c22, c26) with the first identity; then A >= 0 (or B >= 0) supplies enough f1, f4, f7 (or f0, f2, f12) | PROOF-READ (the cancellation adds only f3, f5, f11, f21, so the coefficients used next are untouched) |
| 5 | the A- and B-branch decompositions hold on the whole 8-dimensional space {Pc : (Pc)_T = 0} | EXACT |
| 6 | the three cases (A, B >= 0; A < 0; B < 0) are exhaustive | PROOF-READ |
| 7 | the A < 0 branch is safe | CITED (A-branch note, audited here) |
| 8 | the B < 0 branch is safe | CITED (B-branch note, not audited) |

## D. Falsification

None was run at the level of this assembly document: its content is finite integer identities, all replayed. The branch falsification is in the A-branch report.

## E. Verdict

**confirmed-modulo-citations.** The assembly is correct. Its validity depends on the B-branch note (not audited), the A-branch note (audited), the SOC and four-coefficient inequalities (source unresolved) and the paper's global results.
