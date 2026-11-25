# Audit: final_referee_inventory82_native_uniform_20261023 (primary proof of paper case T_83)

Evidence level: EXACT for every displayed identity, the curve C(alpha,k) and its normalization, the fibre alpha = u and both branch limits, the source rays, coefficient bounds and boundary coverage (53 checks, `replay/final_referee_inventory82_native_uniform_20261023__check.py`, all pass); PROOF-READ for the case split R = 0 / R != 0, the puncture analysis, the compact exclusion and the Rouche restoration of u = 0; CITED for the boundary cases, the degenerate-spectra theorems, Lemma 3.5 and Prop. 4.1/Lemma 4.2; NUMERIC for falsification.

Verdict: **confirmed-modulo-citations**.

## Deviations and limits

- PDF: `proofs/manin_transfer_20261021/final_referee_inventory82_native_uniform_20261023.pdf`, 5 pages, read in full; pages 3-4 checked against the page images.
- Normalization for the replay: A = (inf, 0, 1, z), C = (x, y, v, t), F_sigma = prod over the three finite A_i of (A_i - C_sigma(i)) (the factor at A_1 = inf is 1). With this choice E = F5 - F3 = (v-t)(z-y) and every displayed formula holds verbatim, including the sign of (20).
- The supplement's checker and JSON were not used.
- The degenerate-spectra inputs ("uniform within-spectrum pair law", "uniform shared-zero theorem") are read as paper Lemma B.2 and Theorem B.1; not audited here.

## A. Statements

**Theorem 1 (p. 1).** For every x in X with P c(x) >= 0 and contacts containing T = {3,5,6,9,14,18} (non-strict), Phi_{c(x)}(O) >= 0 for every real orthogonal O; with the two-real-response theorem (paper Thm 1.2) this gives the unitary version.
- T is literally paper T_83 (identity group element). Reference permutation pi* = pi_3.
- Spectra: eight distinct projective points by the curve argument; repeated points within a spectrum and cross collisions by the two cited theorems (p. 1). Strength: at least (*) for T_83.

Unlabelled claims: source (2) and minor events; bounds (3), (4); boundary list; Cramer system (8)-(11); coordinates (13)-(15); identities (16), (17); curve (18)-(20); coordinate functions (21)-(23); fibre (24)-(27); nonconstancy; u = 0 restoration.

## B. Dependencies

| dependency | kind | used for | status |
|---|---|---|---|
| Prop. 4.1 / Lemma 4.2 | core-doc | source = the 8 rays (2) | CITED (generator list EXACT) |
| new ties 0, 1, 7, 8, 19, 22 | paper-case | boundary where theta_28 need not vanish | covered by paper T_40 (stage 1) for every one of them (also T_36, T_38, T_45, T_76, T_79, T_90); PDF says only "previously proved contact-family list" |
| new tie 20 | seven-contact | H = {0,1,6,9,16,17,22} = paper Table 5 T^(7)_21 (stage 1) | CITED; g maps H onto T u {20} EXACT |
| within-spectrum repetitions, cross collisions | common-lemma | reduction to eight distinct points | CITED (Appendix B) |
| bounded-coefficient harmonic lemma | common-lemma | compact exclusion | CITED (paper Lemma 3.5) |

All covering statements are at stages < 9; T_83 does not use T_84 or T_85 (the catalogue's listing of the T_84 document as "used in T_83" is not reflected in this PDF).

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | source rays 1,2,5,9,28,29,39,59; only f28 non-minor; contact rank 6, source rank 8; minor events as listed | identity | EXACT (+CITED completeness) | lib |
| 2 | (3) theta <= g12, g13, g2 and (4) -L_O <= M min(g12,g13,g2) | inequality | EXACT / PROOF-READ | f28 = 1 at 12, 13, 2 |
| 3 | ties not forcing theta = 0 are exactly 0,1,7,8,19,20,22; others give a minor sum | identity | EXACT | |
| 4 | boundary ties covered (5) | citation | CITED; inclusions EXACT | group search |
| 5 | (7) E = (v-t)(z-y); F9 - F3 = v(1-t)(x-y) so r != 0 | identity | EXACT | |
| 6 | (8) four contact rows, linear in (vt, v, t, 1) | identity | EXACT | |
| 7 | (9) Cramer and (10) Delta = -xyz(z-1)(x-1)(x-y)/(v-t) at raw ratios (nonzero at eight distinct points) | identity | EXACT | |
| 8 | (11) det M = -xz(z-1)(x-y) Q, (12) Q 15 terms, (14) Q = 0 iff r + u beta - h alpha - alpha beta = 0 | identity | EXACT | |
| 9 | (15) inverse; denominators nonzero | identity | EXACT | |
| 10 | (16) H20 = u - alpha, H12 = h + beta, r - hu = (alpha-u)(beta+h) as raw product identities | identity | EXACT (all three hold identically in x,y,z,v,t) | |
| 11 | R = 0 layer: alpha != u gives g12 = 0 (minor sum); alpha = u gives g20 = 0 (T^(7)_21) | case split | PROOF-READ | |
| 12 | (17) beta, H12 = R/(alpha-u) | identity | EXACT | |
| 13 | (18) N0 Delta - N1 N2 = z(z-1)y(y-z)(x-1)(x-y) S64, 64 terms | identity | EXACT | |
| 14 | (19) after (15), (17): numerator = (unit) alpha k^2 C, C 179 terms, bidegree (5,2); (20) C(u,k) = -R^3 k(uk-1) with the stated sign | identity | EXACT | |
| 15 | Cramer coordinate functions reproduce the fixed values r,s,u,h,1 and (17), (21), (22), (23) | identity | EXACT (identically in (alpha,k)) | |
| 16 | H12 bounded off alpha = u, incl. alpha = inf | analytic | PROOF-READ | |
| 17 | u != 0: fibre alpha = u is {k=0, k=1/u}; leading coefficient -u R^3 != 0; (24) smooth points | identity | EXACT | |
| 18 | (25) first-order jets | identity | EXACT | implicit derivative |
| 19 | (26) H13 -> -1 and (27) H2 -> -1 along the two branches | analytic | EXACT | sympy limits with arbitrary second-order term |
| 20 | nonconstancy: H12 or H8 nonconstant (H12 nonreal at the point, alpha != 0) | argument | PROOF-READ | |
| 21 | compact exclusion for u != 0, R != 0 | analytic | PROOF-READ (lemma CITED) | see note 1 |
| 22 | u = 0 restoration by one-variable Rouche | analytic | PROOF-READ (argument made explicit by auditor) | see note 2 |

Note 1. The curve in X is x(p) = (i eta / E(p)) m(p), a rational image of the plane curve C(alpha,k) = 0 (the map (alpha,k) -> configuration is rational by (15), (17), (9)), so it is algebraic and inside the T-stratum with fixed real ratios. Gaps are -eta Im H_i. The bound (4) on the negative component D needs only one of g12, g13, g2 bounded at each puncture: g12 everywhere except over alpha = u, where by steps 17-19 the only points are the two smooth branch points and there g13 or g2 is bounded. Collisions, E = 0 and points at infinity are punctures covered by the same meromorphic-function analysis. New ties at ordinary points are steps 3-4.

Note 2. When u = 0 the fibre alpha = 0 meets k = inf, so the direct analysis fails; the PDF perturbs u. The details: R = r != 0, C(alpha,k) is a nonzero polynomial whose coefficients depend polynomially on (r,s,u,h); on a complex line through the point where the zero is isolated, Rouche gives zeros for every u' near 0; by (15), (17), (9) each zero is an eight-distinct configuration with ratios (r,s,u',h). The 18 strict gap signs, the real phase and Phi(O) < 0 are open conditions, so for small real u' != 0 (with R' = r - hu' != 0) one gets a counterexample of the already excluded kind. No closure or density claim is needed. PROOF-READ.

## D. Falsification

NUMERIC: 6 points on the T_83 stratum (seeds 500-505, margin 1e-3), 10-start BFGS over SO(4): worst min Phi / max gap = 1.4e-12. No violation.

## E. Verdict

confirmed-modulo-citations. No gap found. Minor presentation remark: the boundary paragraph (p. 1) does not name which earlier case covers each of the ties 0, 1, 7, 8, 19, 22; the auditor found T_40 covers all of them.
