Evidence level: EXACT for the cone facts, the endpoint polynomial N(h,Y) and the final factorization (symbolic); EXACT-pointwise for the spectral identities (1)-(4) and the later chart identities (exact rational arithmetic at 12 random points, with v and conj(v) independent, which is a randomized identity test and not a symbolic proof); PROOF-READ for the sign analysis; NUMERIC for the chart's surjectivity onto the contact stratum and for the cited four-coefficient lemma; CITED for that lemma and for the global inputs.

# root_moc4_T10_B_branch_20261022: audit (branch B < 0 of paper case T_11)

Document: `proofs/manin_transfer_20261021/root_moc4_T10_B_branch_20261022.pdf`, 4 pages, dated 22 Oct 2026. It is a branch document, not the primary proof of a paper case. It is one of the two branch propositions T_11^(1), T_11^(2) of the paper's Table 6, stage 2.

## Deviations and limits

- **A full symbolic proof of the spectral identities in (h,t,s,l,v,conj v) was attempted and timed out** in sympy. These identities are therefore replayed by exact rational evaluation at 12 random points, treating v and conj(v) as independent rationals. Every identity was checked this way and every check passed. This is strong evidence, but formally probabilistic (Schwartz-Zippel). The identities that matter for the conclusion, namely N(h,Y), the four endpoint values and the final factorization, are proved symbolically.
- I did not derive exactly that the six contacts force h, t, s, l to be real: that is, which contact differences have those quotients. Two pieces are exact: contact 2 gives kappa = i mu/((z-c0)(v-r)), since the normalised difference equals 1 identically, and contact 1 gives k real, since its normalised difference is -1/(k-1). For the rest I have only a NUMERIC check: 12 points on the stratum with contact set exactly T, where h, t, s, l were real to rounding. See gap G1.
- The four-coefficient inequality on (f4, f7, f3, f22) is cited from an unnamed supplement result. It was not audited, only tested numerically.
- The replay script is `replay/root_moc4_T10_B_branch_20261022__check.py`; its output file shows 0 failures.

## A. Statements

Contact set T = {0,1,2,4,7,9}. This is the image of paper T_11 = {0,1,3,5,6,8} under g = (id, (0,1,3,2), 0), i.e. sigma -> sigma o (3 4) (EXACT). The reference is pi_0. Generators: {0,1,2,3,4,5,7,11,12,21,22,26}, with non-minor f22 = h_{3412,+} and f26 = h_{4132,-} (EXACT). Set A = g16 - g3 and B = -g3 + g6 - g12 + g13 + g18. Then A + B >= 0 on N_T (EXACT: it holds on every generator).

**Theorem 1 (p. 1).** For every x in X with P c(x) >= 0, contacts containing T (non-strict, extra contacts allowed) and B < 0: Phi_{c(x)}(O) >= 0 for every O in SO(4), and hence the unitary statement by Thm 1.2. In the chart of eight distinct points with v nonreal and b = g8 > 0, the coefficient bound ab(c-d) >= d^2(a+b+c+d) (a = g15, b = g8, c = g18, d = -B) holds strictly.

Scope: the branch A < 0 is explicitly not covered; it is the separate companion branch. All spectra are covered: collisions go to App. B / Thm C.1, real v makes the A-spectrum real (Thm C.1), and b = 0 goes to paper T_7.

## B. Dependencies

| dependency | status |
|---|---|
| canonical profile g = g23 f1 + g18 f3 + g15 f4 + g21 f5 + g8 f7 + g6 f11 + g19 f21 + d f22 | EXACT on W_T |
| **four-coefficient inequality**: L_O(a f4 + b f7 + c f3 + d f22) >= 0 whenever a,b,c,d >= 0 and ab(c-d) >= d^2(a+b+c+d) | CITED (supplement; the document does not name its source); NUMERIC test passed (D) |
| b = 0: T∪{8} contains {0,1,2,4,7,8} = paper T_7 (stage 1, before this stage-2 branch) | EXACT resolution; T_7 CITED |
| v real: A_sp = (inf,0,1,z) with z = hv/(1+(h-1)v) real, so one whole spectrum lies on a generalized circle | paper Thm C.1, CITED |
| collisions within a spectrum or across | paper Thm C.1 / Thm B.1, CITED (App. B audited by the project) |
| projective covariance | PROOF-READ (common scalar absorbed in kappa) |
| paper Thm 1.2 (unitary to orthogonal) | CITED |

## C. Steps

| # | claim | status | how |
|---|---|---|---|
| 1 | canonical profile for d = -B > 0; A + B >= 0 | EXACT | basis of W_T; values on generators |
| 2 | chart (z, k, c0, r, w) in terms of h, t, s, l, v satisfies the six contacts; the ratio formulas for h, t, s, l are consistent | EXACT-pointwise | 12 exact points |
| 3 | contacts imply kappa = i mu/((z-c0)(v-r)) and h, t, s, l real | EXACT for kappa (D2 = (z-c0)(v-r)) and for k (D1/K = -1/(k-1)); NUMERIC for h, t, s, l | **G1 (presentation)** |
| 4 | the denominator of r has no zero at nonreal v; N0 != 0 because r != 0 | PROOF-READ + EXACT-pointwise (r numerator = -h v N0) | — |
| 5 | identities (1)-(4) | EXACT-pointwise (g19 = G-d and g21 = V-b-d are EXACT linear identities on W_T) | 12 points |
| 6 | a > 0, G > 0, g6 > 0 | PROOF-READ | needs M1, M2 != 0. The document's justification is terse; argument supplied: M1 = 0 at nonreal v would force s(k-1) = D and Dh = 0, which is impossible; J0 != 0 by J0 g6 = -hb with h, b != 0, so the constant term of M2 is nonzero |
| 7 | t > 1; h <= h* < 1; D > 0 | PROOF-READ + EXACT (q/a formula) | — |
| 8 | sign rules sgn J0 = -sgn h, sgn s = sgn h, rho l > 0, l N0 > 0 | PROOF-READ | — |
| 9 | h < 0 impossible (via delta = P - Y eta and (l-1)NP + (h-l)N0 = C) | PROOF-READ + EXACT-pointwise | — |
| 10 | h > 0 and l > 0 impossible (cases k >= 0 and k < 0) | PROOF-READ | k < 0 implies D > h, which I checked by hand |
| 11 | C < 0 and 0 < Y < 1/t | PROOF-READ | Y < 1/t is equivalent to lD < 0 |
| 12 | E = pZ + Y, W = -N0 Dn > 0, R = tYEW/(Z eta), the G/a formula, the k g21 formula, B0 > 0 | EXACT-pointwise + PROOF-READ | — |
| 13 | the g19 and g21 formulas, so LV - LG >= 0; the k = 0 singularity of LV is removable (denominator 2 eta (tY-1) B0) | EXACT-pointwise; singularity EXACT symbolic | — |
| 14 | LV - LG = -h N(h,Y)/(2 eta delta Z J B0), with N affine in h and in Y; the four endpoint identities | EXACT (symbolic) | — |
| 15 | N decreasing in h on 0 <= Y <= 1/t; K <= 0 implies N > 0, a contradiction; so K > 0 | PROOF-READ | interpolation argument correct; strictness from Y > 0 |
| 16 | ab(c-d) - d^2(a+b+c+d) = (a+d)K > 0 when c = (t-1)(a+d) | EXACT (symbolic) | — |
| 17 | conclusion from the four-coefficient lemma plus the nonnegative minor terms | CITED | — |

Analytic checks (i)-(vi): the argument is pointwise algebraic on the chart, with no curve or puncture argument. Extra ties are allowed throughout: only g8 > 0 is used strictly, and g8 = 0 is cited to T_7. Degenerate spectra go to the cited global theorems.

## D. Falsification

- The cited four-coefficient lemma (NUMERIC, seed 11): 120 samples, two thirds of them exactly on the boundary ab(c-d) = d^2(a+b+c+d). The worst normalised minimum of Phi over SO(4) (6 BFGS restarts) was +7e-17. As a control outside the region (a = b = 1, c = d = 0.9), the minimum of Phi is -0.0065, so the test can detect failure.
- Chart sampling (NUMERIC, seed 11): 200000 random parameter draws gave only 4 supported points with B < 0. At all of them, ab(c-d) - d^2(a+b+c+d) >= 1.0e-2 and the minimum of Phi was >= -1e-12.
- Surjectivity probe: 12 lib stratum points with contact set exactly T (margin 0.01), Moebius-normalised. The h, t, s, l computed from the ratio formulas were real to rounding.

## E. Verdict

**gap-found**, with one presentation-level gap.

- **G1** (p. 1, "the original product differences give the real ratios"): the document does not show which contact differences give the quotients t, s, l. So the claim that every configuration with contacts T lies in the real-parameter chart is asserted, not shown. Exact: kappa from contact 2, and k = 1 + t(h-1) real from contact 1. Numeric: realness of h, t, s, l on 12 stratum points. Believed true. Repair: display each parameter as a quotient of two contact differences, which is a short computation.

Otherwise everything was replayed (EXACT or EXACT-pointwise) or proof-read. One step rests on an auditor-supplied argument: that M1, M2 != 0. The conclusion is CITED to the unnamed four-coefficient lemma, which survived a boundary-heavy numeric test. If G1 is repaired, the verdict would be confirmed-modulo-citations: four-coefficient lemma, T_7, Thm B.1, Thm C.1, Thm 1.2. This document covers only B < 0. The A < 0 branch of paper T_11 must be audited separately.
