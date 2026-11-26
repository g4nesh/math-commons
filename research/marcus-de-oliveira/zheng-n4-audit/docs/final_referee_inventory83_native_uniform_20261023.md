# Audit: final_referee_inventory83_native_uniform_20261023 (primary proof of paper case T_84)

Evidence level: EXACT for every displayed identity, the curve C, the fibre and infinity analysis, the source rays, the coefficient majorants and the boundary coverage (40 checks, `replay/final_referee_inventory83_native_uniform_20261023__check.py`, all pass); PROOF-READ for the pole classification and the compact-curve argument; CITED for the boundary cases, the degenerate-spectra theorems, Lemma 3.5 and Prop. 4.1/Lemma 4.2; NUMERIC for falsification.

Verdict: **confirmed-modulo-citations**, with one presentation gap (the list of boundary cases is not in the PDF; the auditor supplies it).

## Deviations and limits

- PDF: `proofs/manin_transfer_20261021/final_referee_inventory83_native_uniform_20261023.pdf`, 4 pages, read in full.
- Normalization used for the replay: A = (0, 1, z, w), C = (inf, v, x, t), [A_i, C_j] = A_i - C_j for finite C_j and [A_i, inf] = -1. With this choice E = F12 - F2 = (z-v)(w-t) and F12/E = x exactly as displayed; any other constant for [A, inf] only rescales all products by a common factor.
- The PDF says the boundary ties are "recorded in the certificate" and does not list them. The certificate was not used; the coverage below is the auditor's own computation.
- Cited theorems for repeated spectra / cross collisions are not named in the PDF ("retained full theorems"); the auditor reads them as paper Appendix B (Theorem B.1, Lemma B.2). Not audited here.

## A. Statements

**Theorem 1 (p. 1).** For every x in X with P c(x) >= 0 and contacts containing {3,4,8,10,12,13} (non-strict), Phi_{c(x)}(O) >= 0 for all real orthogonal O. Equivalently with the working set T = {2,5,9,11,12,13}.
- {3,4,8,10,12,13} is mapped onto T by g = (rho = id, rho' = 1243, eps = 0) (EXACT), and T is literally paper T_84.
- Strength: at least (*) for T_84. Eight-distinct spectra are treated by the curve argument; repeated eigenvalues and a_i + b_j = 0 are delegated to cited theorems (p. 1, last lines of the boundary paragraph). Reference permutation pi* = pi_12.

## B. Dependencies

| dependency | kind | used for | status |
|---|---|---|---|
| Prop. 4.1 / Lemma 4.2 (82-facet completeness, normal cone) | core-doc | source = 11 rays | CITED (lib replays the certificate; generator list EXACT) |
| boundary ties 0,4,7,8,15,16,19,23 | paper-case | new contact where f40 or f68 vanishes | auditor-supplied: covered by paper T_49 (stage 4), T_48, T_59 (stage 8) or Table 5 T^(7)_23 (stage 2) |
| boundary tie 3 | paper-case | idem | covered by T_22 or T_23 (stage 1) |
| boundary ties 6, 17 | paper-case | idem | covered by T_25 or T_28 (stage 1) |
| repeated spectra / cross collisions | common-lemma | reduction to eight distinct points | CITED (Appendix B) |
| bounded-negative-part subharmonic lemma (multi-profile) | common-lemma | compact exclusion | CITED (paper Lemma 3.5) |
| simultaneous two-real-response theorem (Thm 1.2) | core-doc | unitary corollary only | CITED, not needed for (*) |

All covering cases are at stages < 9 (T_84's stage), so the dependency is acyclic. The supplement catalogue lists this document as used in T_83; T_83 is in the same stage 9, but the auditor's T_83 report shows T_83 does not need T_84.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | source rays S = (0,1,7,12,39,40,41,45,67,68,72); only f40, f68 non-minor; rank 8 | identity | EXACT (+CITED completeness) | lib generators |
| 2 | (2) Theta = theta40 + theta68 <= g10, g14, g18 | inequality | EXACT | f40, f68 >= 1 at 10, 14, 18 |
| 3 | new ties outside {0,3,4,6,7,8,15,16,17,19,23}: both non-minor coefficients vanish, minor sum | identity | EXACT | profiles |
| 4 | new ties in that list are covered by earlier cases | citation | CITED; list not in PDF (presentation gap G1), inclusion EXACT | group search |
| 5 | (3) L_O >= -M Theta >= -M B, B in {g10, g14, g18}, any choice per puncture | inequality | PROOF-READ | |
| 6 | E = (z-v)(w-t), height x, H2 = -1, H12 = 0 | identity | EXACT | |
| 7 | elementary forms of F5/F13, F9/F13; (4) v, t solve them | identity | EXACT | |
| 8 | (5) nonzero factors; V != 0 from Vv = x(x+u)w; Td - Tn = x(x+s)(1-z) excludes Td = 0 | identity + argument | EXACT / PROOF-READ | |
| 9 | H13 = r is (x+r) C = 0, C 53 terms, bidegree (2,2) | identity | EXACT | |
| 10 | (6) H11 = K(x) on C = 0, so fixing x keeps the sixth contact | identity | EXACT (holds only on C = 0, as the PDF says) | pseudo-remainder |
| 11 | (7) [zw^2]C = -x(x+u)(x+r) != 0: a curve through the point | identity | EXACT | |
| 12 | (8) H10 = J on C = 0; poles of J only at z = x or w = inf | identity + analytic | EXACT / PROOF-READ | |
| 13 | (9), (10) H14, H18 identities, valid before the curve equation | identity | EXACT | |
| 14 | (11) fibre z = x: exactly w = x and w = w0, distinct by (12) | identity | EXACT | |
| 15 | at (x,x): H14 holomorphic, value V/(x+r) | identity | EXACT | |
| 16 | at (x,w0): v = t = x, H18 regular with value 0, cleared denominator x(x+s)(1-x)(x+r) | identity | EXACT | |
| 17 | (13) w = inf only over z = x(r-s)/(x+r) (!= x since x != -s); there H14 -> 0 | identity + analytic | EXACT / PROOF-READ | |
| 18 | z = inf, w finite: J -> -x | analytic | PROOF-READ | |
| 19 | corner z = w = inf: Ctilde(0,0) = 0, d_tau Ctilde = -x(x+u)(x+r) != 0, smooth graph tau(sigma), J bounded | identity + analytic | EXACT / PROOF-READ | |
| 20 | nonconstancy: x F0 = (x-z) F14 gives z, (8) gives w, (4) gives v,t | identity + argument | EXACT / PROOF-READ | |
| 21 | compact exclusion via Lemma 3.5 and the pole classification | analytic | PROOF-READ (lemma CITED) | see note |

Note on step 21. The PDF says only that a positive constant "would force all anchored coordinates to be constant". The missing line: if -L_O (extended by 0) is a positive constant, the negative component is the whole punctured curve, so every noncontact gap -eta Im H_i is > 0 there; a nonconstant meromorphic H_i on a compact Riemann surface takes real values on an infinite set, which is not contained in the finite puncture set, so a contradiction follows from step 20. Supplied by the auditor; presentation level only.

Every puncture is covered: J (g10) is bounded except at z = x (finite w) and at w = inf; those points are exactly the fibre points of step 14 and the point of step 17, where g14 or g18 is bounded instead. Spectral collisions and other inverse exceptions are punctures at which J, H14 or H18 is bounded by the same classification. No perturbation or genericity argument is used, and x need not be real.

## D. Falsification

NUMERIC: 6 points on the T_84 stratum (seeds 500-505, margin 1e-3), 10-start BFGS over SO(4): worst min Phi / max gap = 1.5e-12. No violation.

## E. Verdict and gaps

confirmed-modulo-citations.

- G1 (presentation; p. 1, paragraph "Every relevant additional tie"): the boundary ties and the cases covering them are not stated in the PDF. Believed true: the auditor computed that the eleven ties 0,3,4,6,7,8,15,16,17,19,23 are each covered by an earlier-stage paper case or Table 5 set (listed in B). Repair: add the table.
