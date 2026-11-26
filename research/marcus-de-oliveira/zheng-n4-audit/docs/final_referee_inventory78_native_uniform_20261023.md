Evidence level: EXACT for the source and its seven frontiers, the first-five inverse (2), the formulas (3)-(7), the 282-term curve and its sections (8)-(12), both Newton models, the slope cubic, the residue cancellation, the infinity limits and the constant-polynomial exception (symbolic in r, s, t, u); PROOF-READ for the puncture analysis and the restoration; CITED for the bounded-negative-part lemma, Appendix B, the Hessian-kernel height exception and the earlier cases, including T^(7)_25.

# final_referee_inventory78_native_uniform_20261023 (primary proof of paper T_79)

PDF: proofs/manin_transfer_20261021/final_referee_inventory78_native_uniform_20261023.pdf (5 pp.).
Replay: `replay/final_referee_inventory78_native_uniform_20261023__check.py` and its output (all EXACT checks pass).

## Deviations and limits

- No numeric minimization of the main statement. Intermediate claims were replayed exactly instead.
- The height exception argument (H_kappa v0 = 0, kernel = singleton span, Birkhoff) was read but not replayed. It rests on (3.10)-(3.11) of the paper, which the library checks.
- The document labels earlier cases by inventory index (paper index minus one). I resolved every label by contact set.
- The statement labelled "corner" (p. 3) in the catalogue is the x = z = 0 Newton step, step 11 below.

## A. Statement

Theorem 1. Statement (*) holds for T = {1,3,4,8,11,12} = paper T_79 (literal), with reference pi* = pi_1. Contacts are non-strict, every O in SO(4) is covered, and all spectra are included. Repeated and cross-colliding spectra go to Appendix B.

## B. Dependencies

| dependency | check |
|---|---|
| source (1): generators f0, f1, f9, f10, f56, f53, f59, f2 with exclusive gaps 20, 17, 10, 9, 6, 14, 2, 15 | EXACT (lib generators; the kernel has dimension 8 and the exclusive-gap representation holds on all 24 rows for a kernel basis) |
| Z(f56) outside T = {2,9,10,14,15,17,20} | EXACT |
| added 2, 9, 10, 14, 17, 20 -> document labels T12, T69, T70, T77, T57, T39 = paper T13, T70, T71, T78, T58, T40 | EXACT. Each T ∪ {j} contains a G-image of the cited case. All of these have stage <= 6, and T_79 is stage 8 |
| added 15: T ∪ {15} = T^(7)_25 literally (stage 7) | EXACT. The seven-set is proved in final_referee_inventory78_boundary15_uniform_20261023, audited in this batch (verdict gap-found, presentation level only) |
| all other added ties make d = 0 (minors only) | EXACT |
| Appendix B (Thm B.1, Lemma B.2) | CITED |
| bounded-negative-part lemma | CITED (hypotheses compared) |
| (3.10)-(3.11): rank-4 H_kappa with kernel the 10-dimensional entry span; Birkhoff | CITED (paper, lib-checked) |

## C. Steps

| # | claim | kind | status |
|---|---|---|---|
| 1 | chart A = (w,0,1,z), C = (v,inf,kz,x); E = F3 - F1 = (w-v)z(1-x); H4 = (k-1)/(x-1), so k = 1+s(x-1) | identity | EXACT |
| 2 | H8 = r and H11 = t are equivalent to the two displayed linear equations in v and y = w - v | identity | EXACT (polynomial remainder) |
| 3 | Cramer solution (2) solves both; L = 0 is impossible on eight-distinct sources (Y != 0 because x, k, A0, B0, z-1 are pair differences) | identity + analytic | EXACT + PROOF-READ |
| 4 | homogeneous inverse; (3) F1/E = s(x-1); H3, H4, H8, H11 = 1, s, r, t; (4)-(7) for H2, H5, H6, H10 | identity | EXACT |
| 5 | g5 = d + b | identity | EXACT (f56(5) = f53(5) = 1, other source values at 5 are 0) |
| 6 | poles of H6 only at x = 0, k = 0, z = 1, infinity | identity | EXACT (denominator x k (z-1)) |
| 7 | C = N12 - uD12: 282 terms, bidegree (5,3); (10) [x^5z^3]C = s^3(r+s-u-1) (normalizing constant 1) | identity | EXACT |
| 8 | (8) C(x,1) = (x-1)^3 B(x)^2; B linear with leading coefficient s(r+s-t-1) | identity | EXACT |
| 9 | (9) C(1,z) = -(r+t)(z-1)^3; (11) C(0,z) = -(s-1)^2(s-t)^2 z^3 exactly, so x = 0 occurs only at z = 0 or infinity; (12) | identity | EXACT |
| 10 | generic witness (2,3,5,7): all displayed factors and the local discriminant are nonzero | finite | EXACT |
| 11 | x = z = 0: weighted initial form on 3 deg x + 2 deg z = 6 as displayed; branches x ~ z^{3/2}; H5 -> s, so 0 <= d <= g5 -> 0 | identity + analytic | EXACT + PROOF-READ |
| 12 | k = 0: section -s^2(r-1)^2(s-1)^2 z; graph z = (s+u-2)K/(r-1) + O(K^2); H6 + H10 has limit (rs - r + st - s - t - u + 2)/(s-1), real, independent of the second-order term; 0 <= d <= d + g10 | identity + analytic | EXACT + PROOF-READ |
| 13 | z = 1, x != 1: x real by (8); H5 = s there | identity | EXACT |
| 14 | x = z = 1: least degree 3; slope cubic (Yr+Yt+rs-t)(-Y^2+Y(t-u+1)+rs-t); H6 -> (Yr+Yt+rs-t)/Y, bounded for the generically nonzero finite slopes | identity | EXACT |
| 15 | z = inf: H6 -> r + (t+s(x-1))/x; x = inf: H6 -> r+s+(r+s-t-1)/(z-1); (x,z) = (inf,inf) excluded by (10); x = inf, z = 1 excluded by (12) | identity | EXACT |
| 16 | k = 0, z = inf: lowest degree 2; K^2 and tau^2 coefficients nonzero; discriminant s^4(r-1)^2(s-1)^3(...) generically nonzero; H6 along K = lambda tau has a finite limit for lambda != 0 | identity | EXACT |
| 17 | these exhaust the poles of d | case split | PROOF-READ (follows from steps 6, 9, 13, 15) |
| 18 | constant-coordinate exception via the explicit inverse and the height kernel | analytic | PROOF-READ + CITED (3.10-3.11) |
| 19 | C == 0 only at (r,s,t,u) = (0,1,0,0) or (1,-1,-1,-1); the first has [x^5z^2]C = 3; the second has C == 0 and H6 == 0 | identity | EXACT |
| 20 | harmonic exclusion with L_O >= d f56(O), d bounded at every puncture | citation | CITED (lemma) |
| 21 | Rouché restoration from walls | analytic | PROOF-READ (same argument as in the T_77 document, checked by auditor) |

## D. Falsification

No numeric search for violations of the main statement. No intermediate claim failed.

## E. Verdict

confirmed-modulo-citations. Citations: bounded-negative-part lemma, Appendix B, (3.10)-(3.11), paper cases T13, T40, T58, T70, T71, T78, and the seven-set T^(7)_25. The seven-set document's own verdict is gap-found, at presentation level only: its certificate polynomial G is not printed, and the auditor's replay of the claims about G is exact only on slices. This document's dependence on it is therefore conditional on that replay.
