# Audit: blind_moc4_T27_uniform_20261022 (primary proof of paper case T_28)

Evidence level: EXACT for the generator list, the three cone identities, the symmetry, the three orbit witnesses, C, Q, the radial factorization (G9, 2124 terms), the split G9 = k0 + k1 j + k2 j^2, k0 = C G6, the four-chamber certificate (8), the e = 0 claim, and all eight Bernstein charts (all 24 tables, with term counts that match the document exactly, including the 105 + 197 mixed-sign E-groups of chart row 4, each settled by an exact real-root count). PROOF-READ for the chamber logic and the boundary and continuity arguments. CITED for the pointwise capacities (3)-(6), the double SOC at ab = 1, the two four-coefficient inequalities, and paper Lemma 4.2. NUMERIC tests of every cited capacity on its boundary.

PDF: proofs/manin_transfer_20261021/blind_moc4_T27_uniform_20261022.pdf (5 pages). Auditor label C11-T27-T28.

## Deviations and limits

- Nothing from the supplement was used: no stored Hessian input, no Sturm JSON. The Sturm part was redone with sympy's exact `count_roots` on each grouped univariate polynomial in E. This is an exact real-root count, but it is not a replay of the document's own Sturm sequences.
- For three of the four chambers my denominator clearing in (8) differs from the document's, so the coefficient counts differ: (135,94,94,105), (135,98,164,161), (135,164,98,161), against the document's (135,182,182,342), (135,188,164,291), (135,164,188,291). The fourth chamber matches exactly: (130,180,180,268). Nonnegativity holds in all four, and a positive denominator does not affect the conclusion.
- The cited capacities are tested numerically only; their proofs are in other documents ("T24 core" = a core statement of paper case T_25, the four-coefficient orbit, the double SOC).

## A. Statements

**Theorem 1 (p. 1).** Let T = {2,3,5,6,9,12}, exactly the paper's T_28 (identity group element, EXACT). For every x in X with P c(x) >= 0 vanishing on T, and every O in SO(4), Phi_{c(x)}(O) >= 0. Contacts are non-strict. Every spectrum is covered, because the radial necessity holds for every x in X. The phrase "closed limiting spectral configurations covered by the general radial-necessity theorem" therefore does not restrict the statement. This is (*) for T_28 at full strength.

Intermediate claims (d = 1 after homogeneity), for the representation (1) g = q f7 + p f39 + j f5 + c f41 + e f2 + d f40 + b f0 + a f67:
- (R) every realizable g of form (1) satisfies G9 = 0;
- (K0) k0 = C G6 > 0 on the closed hard chamber {Le e < 1, Lp p < 1};
- (E0) G9 > 0 when e = 0;
- (B) in each chamber, G9 > 0 for 0 <= j < J, where J is the cited capacity endpoint.

## B. Dependencies

| dependency | kind | status |
|---|---|---|
| paper Lemma 4.2 / Prop. 4.1 | core paper | CITED |
| X contained in Sing C | core paper | EXACT in lib |
| G invariance (4.7) | core paper | EXACT in lib |
| "T24 core" capacities (3)-(6) on (f5,f2,f0,f67,f40), the image of (f59,f53,f55,f67,f54) under a group element (EXACT). Common zero set {3,4,7,8,10,12,13,17,22}, which contains paper T_25 = {3,4,7,8,10,12}, so supplement "T24" = paper T_25 | other-doc (core of T_25) | CITED; NUMERIC on each boundary: no violation |
| double SOC at ab = 1: j >= 1/e | other-doc | CITED; NUMERIC: no violation |
| four-coefficient inequality, base tuple (f0,f6,f10,f54). Maps in order to (f0,f67,f2,40) and to (f67,f41,f39,f40) (EXACT). The same orbit as input (2) of the T_27 proof | other-doc | CITED; NUMERIC: no violation |
| Theorem 1.2 (two real orthogonals) | core paper | mentioned only; not needed for the SO(4) statement |

Not used: Appendix B, Corollary B.3, single-collision compensation, curve or puncture arguments, seven-contact cases.

## C. Steps

| # | claim | kind | status | how |
|---|---|---|---|---|
| 1 | F_T = {0,1,2,5,7,9,39,40,41,67,68}; non-minor 40, 68 | finite | EXACT | lib |
| 2 | f40+f68 = f0+f7+f41+f67; f40+f1 = f0+f7+f39; f40+f9 = f2+f7+f41 | identity | EXACT | |
| 3 | a symmetry fixing T swaps 40 and 68 | finite | EXACT | g = ((2,0,1,3),(2,0,1,3),1), the document's description. It also swaps 0<->41, 1<->39, 2<->9, 7<->67 |
| 4 | cancellation gives (1) with d > 0, or an all-minor cost | case split | PROOF-READ | |
| 5 | radial necessity: C(Re x) = 0, det H(Re x) = 0 | analytic | PROOF-READ | as in the T_27 audit |
| 6 | C(y0 + z kappa) = C + zQ; Q = d^2+de+dp+ep-ej; C has 25 positive coefficients including d^3 | identity | EXACT | |
| 7 | C > 0, so Q = 0 has no realizable height | logic | PROOF-READ | |
| 8 | Q^4 det H(y0 - (C/Q) kappa) = -4 Q^2 G9^2, with G9 homogeneous of degree 9 and 2124 terms | identity | EXACT | determinant 17446 terms, degree 3 in z; exact square root. There is no C^2 factor here, unlike T_27, which is consistent with the document |
| 9 | G9 = k0 + k1 j + k2 j^2 (deg_j = 2); k0 = C G6 with G6 of 196 terms, bilinear in (e,p); g00 positive with constant 1 | identity | EXACT | |
| 10 | (8): the four functions are nonnegative for all a, b, c, q >= 0, on the four chambers | certificate | EXACT | own charts a = A, b = T/(A(1+T)) or (1+T)/A, c likewise with U. Le, Lp as rational functions with positive denominators |
| 11 | G6 >= g00 (1 - Le e)(1 - Lp p) > 0, so k0 > 0 on the closed hard chamber | inequality | PROOF-READ | the substitution e' = e/(1 - Le e) re-derived |
| 12 | k2 at e = 0 is 0; k1 at e = 0 has 375 nonnegative coefficients; so G9 > 0 for e = 0, j >= 0 | identity | EXACT | |
| 13 | Bernstein identity G9 = (1-t)^2 k0 + t(1-t)(2k0 + J k1) + t^2 (k0 + J k1 + J^2 k2), j = tJ | identity | PROOF-READ | expanded by hand |
| 14 | charts 1-6 (ab < 1; three sub-charts for a,b <= 1 / a >= 1 / b >= 1, each with ac < 1 or ac > 1) and chart 7 (ab > 1, ac < 1): all Bernstein data nonnegative | certificate | EXACT | own substitution and clearing. All 21 term counts equal the document's. Row 4 (a >= 1, ac > 1): 105 and 197 E-groups with a negative coefficient, all with positive constant and no root in [0, inf) (exact root count). All denominators have nonnegative coefficients |
| 15 | chart 8 (ab >= 1, ac >= 1): shared-resource split d = delta + r; piece 2 on the boundary of the 4-coefficient inequality #2; piece 1 gives (11) via (6) after scaling by 1/r; Bernstein tables 30504 / 63441 / 100982 terms, all nonnegative | inequality + certificate | EXACT (tables) + PROOF-READ (split) | I re-derived that the homogeneous form of #2 at (delta a, c, p, delta) is equivalent to (10), and that scaling piece 1 by 1/r turns (6) into (11). Also E_r <= E_ab, so the hard chamber condition Le e < 1 holds where e < E_r |
| 16 | exhaustiveness: hard region = {ab<1} x {ac<1 or ac>1} plus {ab>1, ac<1} plus {ab>=1, ac>=1} plus the boundaries ab = 1, ac = 1, zero coefficients | case split | PROOF-READ | all cases are covered. The boundaries are handled by continuity of the Bernstein data (polynomial in the original variables and in J, with J continuous up to the boundary) together with the strict k0 > 0, which is proved directly on the closed chamber. The argument is valid but terse |
| 17 | conclusion: a realizable root has j >= J, so the cited capacity applies; adding the minor terms gives Phi >= 0 | citation | CITED | |

## D. Falsification

- NUMERIC: each cited capacity at its boundary, 30 random parameter sets, 12 BFGS restarts. Worst min Phi / max gap, in the order (3), (4), (5), (6), double SOC, 4-coefficient #1, #2: 4e-18, 4e-16, 3e-17, 1e-16, 4e-16, 3e-17, 2e-16. Power check: (3) at 0.9 times the bound gives -1.7e-6, so the test detects violations just inside a bound.
- NUMERIC, main statement: 8 stratum points for T_28, worst -2.5e-15 (round-off).
- No violation was found.

## E. Verdict

**confirmed-modulo-citations.** Citations: the T_25-core capacities (3)-(6), the double SOC at ab = 1, the four-coefficient inequality (two orbit images), and paper Lemma 4.2. Every algebraic certificate the document states was replayed exactly, and the replay reproduced the stated term counts.

Presentation remarks (not gaps):
- The boundary paragraph (ab = 1, ac = 1, zero a, b, c, p, q) relies on continuity of the Bernstein data in the original variables. This is valid because the data are polynomials in (a,b,c,e,p,q) and J, and J extends continuously to each boundary used (1/e at ab = 1). The document states this only briefly.
- The four coefficient counts quoted for (8) depend on how denominators are cleared.
