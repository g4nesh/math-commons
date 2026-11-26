Evidence level: EXACT for the chart, curve, inverse, pivotal formulas, residue table, source-cone decomposition, nonvanishing identities and endpoint inventory (replay script below); PROOF-READ for the harmonic-minimum argument, the perturbation step and the two-circle chart lemma (argument supplied by auditor); CITED for Appendix B and the earlier contact cases.

# independent_moc4_T8_polar_reduction_20261022 (4 pages)

Supplement label T8 = paper T_9 = {0,1,2,3,11,12} (same set, no relabelling needed). This is a reduction lemma for T_9. It is not a primary case proof. Its text is repeated verbatim as Section 1 of `independent_moc4_T8_uniform_20261022`.

## Deviations and limits

- The document cites a "companion checker". I did not use it. Every identity below was re-derived from the definitions with my own sympy code: `replay/independent_moc4_T8_polar_reduction_20261022__check.py`, 57 EXACT checks, all pass.
- The two-circle chart is asserted without proof. I supplied a short proof (A.1 below).
- No numeric falsification of the main statement was run. It is a reduction, and every intermediate claim is an exact identity that was replayed.

## A. Statements

**Theorem 1.** Fix O in SO(4). Suppose some x in X has P c(x) >= 0, contacts containing T_9, and Phi_{c(x)}(O) < 0. Then some x' in X has P c(x') >= 0, contacts containing T_9 ∪ {18} or T_9 ∪ {20}, and Phi_{c(x')}(O) < 0. Degenerate spectra are handled by citing Appendix B. The document says explicitly that it does not prove either seven-contact family.

A.1 (auditor's lemma, used implicitly). Take 8 distinct points with no a_i + b_j = 0 and contacts {0,1,2,3} relative to pi* = 0. Normalize A2 = 0, A3 = 1, C4 = oo. Then A4 is real and C3/C2 is real. This is the two-circle chart A = (z,0,1,x), C = (w,v,tv,oo).

Proof: D_s = F_s - F_0 for s = 1,2,3, divided by (A1 - C1), equals p(1-al), p-q and p-q al (EXACT). Here al = A4, p = C2, q = C3. The contacts make D2/D1 and D3/D1 real, say lm and mu. Eliminating q/p leaves a real quadratic in al with root al = 1, so its other root -(1-mu)/lm is real. Since al ≠ 1, al is real. Then q/p = 1 - lm + lm al is real. When lm = 0, q = p, which is a collision. The contact phase gives kappa = i mu/(v(z-w)) with mu real.

## B. Dependencies

| ref | kind | used for | resolution |
|---|---|---|---|
| Appendix B (Thm B.1, Lemma B.2) | common-lemma | degenerate start points and degenerate endpoints (collisions) | CITED (audited by project, holds as written) |
| T_9 ∪ {9} | seven-contact | case u = 0 and endpoint 9 | = G-image of Table 5 T^(7)_1, Stage 1 (EXACT orbit check). T_9 ∪ {9} also contains paper T_8, but T_8 is Stage 4, so the citation must go to T^(7)_1 |
| T_9 ∪ {6} | paper-case | endpoint 6 | contains G-images of T_5 and T_6 (Stage 1) |
| T_9 ∪ {7} | paper-case | endpoint 7 | contains T_4 and T_7 (Stage 1) |
| T_9 ∪ {21} | paper-case | endpoint 21 | contains T_5 (Stage 1) |
| T_9 ∪ {23} | paper-case | endpoint 23 | contains T_4 (Stage 1) |
| T_9 ∪ {18}, T_9 ∪ {20} | seven-contact | conclusion (not used as true here) | = T^(7)_3 and T^(7)_4 (Stage 2). Proved in the uniform document, Thms 3 and 2 |
| double boundary T_9 ∪ {18,20} | other-doc | mentioned only. Logically redundant, since Thms 2 and 3 cover the intersection | — |

Labelling note: the document says the endpoints 6, 7, 21 and 23 contain "T2, T4, T5". Under the "paper index minus one" convention these would be paper T_3, T_5 and T_6. Paper T_3 lies in none of the four sets, and endpoints 7 and 23 need paper T_4. The coverage itself is correct (paper T_4..T_7, all Stage 1). Only the labels are inconsistent (presentation).

## C. Steps

| # | claim | kind | status | how |
|---|---|---|---|---|
| 1 | two-circle chart covers every 8-distinct T_9 point | identity/lemma | PROOF-READ (argument supplied by auditor) + EXACT algebra | A.1 |
| 2 | H0..H3 = 0, 1-x, P, 1-tx; r = P - H12 = z(v-1)(tv-w)/(v(z-w)) ≠ 0 | identity | EXACT | sympy on sigma-point products |
| 3 | inverse formulas give H12 = P-r, H18-H4 = d identically. Conversely, on any chart point d := H18-H4 recovers w and z. Transfer d = (v-1)(eta+b)/(Pv) | identity | EXACT | |
| 4 | H11 = u ⇔ eta^2 - a eta - bc = 0 (numerator = -1 × quadratic) | identity | EXACT | |
| 5 | contact values constant on the curve; pivotal H4..H18 = displayed formulas + real constants (constants printed) | identity | EXACT | reduction mod the quadratic |
| 6 | the curve lies in X, keeps T_9 (complex values held), and the gaps are -mu Im of meromorphic functions, so Phi is harmonic | analytic | PROOF-READ | |
| 7 | poles only over v = 0, 1, 1/t, x/t; regular at v = oo | analytic | EXACT (denominators) + PROOF-READ | |
| 8 | residue table, all 7 rows, and A(R) = 0, R21 or R18 | identity | EXACT | limits along each branch |
| 9 | v = x/t coalescent case: F_v = -b'c; ramified case a = c = 0 ⇔ r = x-t, u = -t(x-1), with a, b, c divisible by (tv-x) | identity | EXACT | |
| 10 | source-cone theorem: g in N_T and g18 >= g4 ⇒ nonnegative minor combination | inequality | EXACT | auditor's decomposition with f0 and f10 (λ-interpolation). Feasibility conditions are exactly g4, g10, g13, g14, g15, g17 >= 0 |
| 11 | lim inf Phi >= 0 at every puncture (real Laurent coefficients, sign of R) | analytic | PROOF-READ | uses steps 7, 8, 10 |
| 12 | min principle; the constant case is excluded by conjugation invariance (a real point over v = 0) | analytic | PROOF-READ | |
| 13 | for nonreal v: d, w, dw ≠ 0. z = w forces (v-1)(tv-x) = 0 | identity | EXACT | b+a-c and b+ta-t^2c identities, z-w factorization |
| 14 | remaining collisions at endpoints are covered by Appendix B | citation | CITED | |
| 15 | exclusions u(r+1)(x-t) ≠ 0: u = 0 is T_9 ∪ {9}. r = -1 and x = t are removed by a small real perturbation through the rational inverse | analytic | PROOF-READ | open conditions, and the chart covers a neighbourhood |
| 16 | endpoint inventory {6,7,9,18,20,21,23}; other j leave only minor generators | finite | EXACT | lib generators |
| 17 | endpoints 6, 7, 9, 21, 23 contain earlier cases | finite | EXACT containment + CITED cases | lib group |

## D. Falsification

Not needed beyond the exact replays: every displayed identity was replayed exactly, and the statement is a reduction. The SOC and seven-contact statements are tested in the uniform report.

## E. Verdict

**confirmed-modulo-citations.** Citations: Appendix B (B.1, B.2), T^(7)_1, T_4, T_5, T_6 and T_7 (all Stage 1, before T_9's Stage 3). Corollary B.3 is not invoked.

Presentation-level gaps:
- p.1: the two-circle chart is asserted without proof. Repaired by A.1.
- p.3: the labels "T2, T4, T5" do not match the paper's cases. The correct ones are T_4..T_7.
- Throughout: the "companion exact file" is the only cited evidence for the identities. Repaired by this replay.
- Endpoint 9: citing the paper's T_8 would be circular in the Table 6 order (T_8 is Stage 4). The text names a "u = 0 seven-contact theorem", which resolves to T^(7)_1 (Stage 1), so the order is fine.
