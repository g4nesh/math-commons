Evidence level: EXACT for the symmetry action, the source, every displayed chart identity, the curve restrictions (7), the limit (8), the Newton polynomial (10) and the three infinity branches (11)-(13) (symbolic in r, s, u, h); PROOF-READ for the curve argument and the restoration; CITED for the bounded-negative-part lemma, Appendix B and the earlier cases.

# independent_inventory77_native_harmonic_uniform_20261023 (primary proof of paper T_78)

PDF: proofs/manin_transfer_20261021/independent_inventory77_native_harmonic_uniform_20261023.pdf (4 pp.).
Replay: `replay/independent_inventory77_native_harmonic_uniform_20261023__check.py` and its output (all EXACT checks pass).

## Deviations and limits

- No numeric minimization of the main statement (the task allows skipping it). Every intermediate identity was replayed exactly.
- I replayed the coefficients [J^4]C and [X^3]C, the section C(J,1), Res_X(C,L), C(r,X) and the V = 0 coefficient. I did not replay Res_X(C,AN) symbolically. It is not needed, because d = -eta Im J is bounded at every finite J anyway (see step 10).
- The Birkhoff height argument at K = 0 was read, not replayed. It is not needed either, for the same reason.

## A. Statement

Theorem 1 (p. 1). The literal contact set is {1,4,5,6,9,12}, which is paper T_78 literally. Statement (*) holds for it: contacts are non-strict, every O in SO(4) is covered, and all spectra are included (collisions and repetitions go to Appendix B). The proof works on T = {0,1,3,10,11,23}. The map sigma -> R o sigma^{-1} with R = (0,3,2,1) sends T to the literal set (EXACT). This map is in G with eps = 1, so (*) transfers.

## B. Dependencies

| dependency | used for | check |
|---|---|---|
| Lemma 4.2 | working source [2,3,4,8,12,13,15,32]; literal source [0,1,2,5,7,9,35,59]; only f32 (resp. f35) non-minor | EXACT (lib) |
| exclusive-gap representation (1) | coefficients (g21,g20,g16,g9,g22,g7,g5,g2) | EXACT. For a basis of the 8-dimensional contact kernel, P c = sum_j (Pc)_{exc(j)} f_j holds on all 24 rows. The representation is therefore unique and nonnegative |
| added ties outside Z(f32) | d = 0, minors only | EXACT |
| ties 5, 7, 9, 20, 21 | earlier six-sets (stage < 6) | EXACT: +5: T4, T11, T17, T31; +7: T29, T38, T70; +9: T29, T31, T35, T47; +20: T31, T36; +21: T29, T34, T39, T46, T49 |
| ties 16, 22 | Table-5 sets T^(7)_18 (stage 4) and T^(7)_20 (stage 5) | EXACT: find_g gives an explicit g. T ∪ {16} also contains T_77, which comes later (stage 9). The document correctly uses the seven-set instead |
| Appendix B | eight-distinct reduction | CITED |
| bounded-negative-part lemma | exclusion | CITED. Its hypotheses are checked in C |

## C. Steps

| # | claim | kind | status |
|---|---|---|---|
| 1 | literal action, source, kernel dimension 8, extra zeros (2) of f32 = {5,7,9,16,20,21,22} | combinatorial | EXACT |
| 2 | r = -(x-z)/(z(x-w)), s = y(t-1)/(z(y-t)), u = H10; inverse (4) for w, t, y; J = (1+rw)/x = H2 | identity | EXACT |
| 3 | t numerator minus denominator = s(x-1)(w-1); y(sz-t+1) = szt, so no extra inverse component | identity + analytic | PROOF-READ (two lines of algebra; not separately machine-checked; the t and y formulas of (4) are EXACT, step 2) |
| 4 | homogeneous inverse; (6) E = s(X-1)(J-r)LK AN; F0/E = (r+1-JX)/K; H2 = J; H5 = r+1-JX | identity | EXACT |
| 5 | every denominator of every H_sigma divides (J-r)(X-1)L AN | identity | EXACT (factor lists of all 24) |
| 6 | C = N23 - hD23: 176 terms, bidegree (4,3); [J^4]C = V(X-1)X^2; [X^3]C = J^2 F_X | identity | EXACT |
| 7 | (7) C(J,1) | identity | EXACT up to sign. My computation gives C(J,1) = +(u-1)(J-r-1)^2[(s+u+h-2)J - rh + r - s - u + 1] with the normalization of (7)'s first line. The document prints a leading minus. The zero set is the same, so nothing downstream changes (presentation) |
| 8 | Res_X(C,L) = -r^2(s+u)^2 J^2(J-r-1)^2 | identity | EXACT up to sign and power-of-J convention (J^3 C(J,(r+1)/J) = J^2 r^2(s+u)^2(J-r-1)^2) |
| 9 | (8) H23 - H22 - J -> -r-1 as X -> inf | identity | EXACT at (r,s,u) = (2,3,4) only (the symbolic limit was not attempted). Not needed (step 10) |
| 10 | every finite-J puncture is harmless | analytic | PROOF-READ. d = g2 = -eta Im J is bounded wherever J is finite. The resultant list and (8)-(9) are stronger than needed |
| 11 | (10) at J = inf, X = tau Y: initial (r+1-Y)(VY+r-V); [J^4X^3]C = V != 0 rules out (inf,inf) | identity | EXACT |
| 12 | (11), (12): branch first corrections r(s+u), -r(s+u)(r-V)/V^2, -(u-1)(h+s+u-2)/V | identity | EXACT (symbolic). Along the trial series all orders below the y2-order vanish identically |
| 13 | (13): H7+J, H16+J, H5+J bounded on the respective branches; leading denominators r(s+u), rs(s+u)V^2(r-V), V | identity | EXACT (symbolic) |
| 14 | real Taylor coefficients make the limits real, so d -> 0; boundedness alone already suffices for the lemma | analytic | PROOF-READ |
| 15 | three branches exhaust the X-degree 3; X = 1 branch simple (derivative V) | analytic | PROOF-READ |
| 16 | K = 0 height poles: d bounded since J is finite; the Birkhoff argument is extra | analytic | PROOF-READ |
| 17 | nonconstancy: J^4 section when V != 0; on V = 0, [J^3X](N23 - hD23) = -r independent of h | identity | EXACT |
| 18 | Rouché restoration from walls; generic witness (2,3,4,5) | analytic | PROOF-READ (same argument as in the T_77 document, checked by auditor); witness values nonzero EXACT |
| 19 | harmonic exclusion | citation | CITED (lemma) |

## D. Falsification

No numeric search for violations of the main statement. No intermediate claim failed. The only discrepancy is the sign in (7), step 7.

## E. Verdict

confirmed-modulo-citations. Citations: bounded-negative-part lemma, Appendix B, T4, T11, T17, T29, T31, T34, T35, T36, T38, T39, T46, T47, T49, T70 (at least one per tie), T^(7)_18, T^(7)_20. One presentation slip: the sign of C(J,1) in (7).
