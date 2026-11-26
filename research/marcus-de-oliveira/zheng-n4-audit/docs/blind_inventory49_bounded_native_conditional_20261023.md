# Audit: blind_inventory49_bounded_native_conditional_20261023

Evidence level: EXACT for every algebraic identity and combinatorial claim (own sympy / exact-rational code, `replay/blind_inventory49_bounded_native_conditional_20261023__check.py`, 54 exact checks, all pass); PROOF-READ for the analytic curve argument and the wall-restoration argument (some steps supplemented by the auditor, marked below); CITED for the boundary cases, the seven-element cases F and G (explicit hypotheses), the bounded-negative-part lemma and Appendices B and C; NUMERIC for the falsification run.

PDF: `proofs/manin_transfer_20261021/blind_inventory49_bounded_native_conditional_20261023.pdf`, 4 pages. Role: branch document (conditional reduction) used by the primary proof of paper case T_50.

## Deviations and limits

- The supplement's companion checker (`..._exact_20261023.py/.json`) was not run or read; every identity was re-derived from the formulas printed in the PDF and from the definitions (sigma-points with b_4 = infinity, the generator vectors f_j of Appendix D as built in `lib/moc4lib.py`).
- The bounded-negative-part harmonic lemma (`final_referee_bounded_bad_harmonic_lemma_20261023`) was read to check that its hypotheses match the use made here; it was not audited (CITED).
- The numeric falsification is weak: 6 floating-point points on the contact stratum (lib `random_point_on_contact_stratum`, seeds 0-5, pi* = pi_1), 12 BFGS restarts each.
- The replay was run through the job limiter with the harness's background option (one job at a time), since the shared slots were busy.

## A. Statements

**Theorem 1 (p. 1).** Let T = {1,3,5,6,10,11} (= paper T_50 literally; the stated representative {1,3,4,6,7,8} maps to T under pi -> column o pi o row with row (0,2,3,1), column (1,0,2,3): EXACT). Assume (H_F) and (H_G): for every x in X with P c(x) >= 0 whose contact set contains F = {0,1,2,7,9,10,14} (resp. G = {0,1,2,7,9,10,23}), Phi_{c(x)}(O) >= 0 for all O in SO(4). Then for every x in X with P c(x) >= 0 and (P c(x))_pi = 0 for all pi in T (non-strict contacts), Phi_{c(x)}(O) >= 0 for all O in SO(4). F and G are Table 5 sets T_7^(7) and T_10^(7). The proof additionally uses the cases T_10, T_43, T_66 (six-element), A = T_9^(7), B = T_8^(7), and Appendices B, C for degenerate spectra; with those, all spectra (repeated eigenvalues, a_i + b_j = 0) are covered.

Sub-claims (all restated in the paper's notation):
- (1) Every c in N_T has a nonnegative expansion in the ten generators [f0,f1,f2,f3,f4,f8,f10,f27,f30,f59] with coefficients given explicitly in terms of the kernel coordinates x_j and a parameter zeta.
- (2)/(3) For every O in SO(4): L_O := Phi_c(O) >= d f30(O), 0 <= d <= g4, d <= g15, hence -L_O <= M min(g4,g15), M = max(0, -Phi_{f30}(O)).
- (4) Extra ties with f30 > 0 force d = 0; the remaining possible extra ties are {0,7,12,13,14,16,20} and each is covered by a cited case.
- (5)-(9) Explicit inverse parametrization of the eight-distinct stratum and the formulas for H_i := (F_i - F_1)/E.
- (10)-(12) For generic real (r,s,u,h) the fibre W = 0 has a single branch and H15 = -su/(r(r+s)) W + O(W^2) on it.
- (13) The Z-coefficient of the last-contact equation is a monic cubic in W, so the last germ is nonconstant at every regular point; small real perturbations of (r,s,u,h) are realized (Rouche).

## B. Dependencies

| dependency | kind | resolved to | check |
|---|---|---|---|
| F = {0,1,2,7,9,10,14} | seven-contact, hypothesis | Table 5 T_7^(7); tie 13 (identity row/column action (0,2,3,1),(0,1,2,3) maps F onto T+{13}) | EXACT action; inequality is a hypothesis |
| G = {0,1,2,7,9,10,23} | seven-contact, hypothesis | Table 5 T_10^(7); tie 20 | EXACT action |
| A = {0,1,2,7,9,10,17} | seven-contact | T_9^(7); tie 14 | EXACT action; CITED |
| B = {0,1,2,7,9,10,16} | seven-contact | T_8^(7); tie 12 | EXACT action; CITED |
| "+0 -> T9" | paper case | {0,1,3,5,7,8} = paper T_10 (stage 2) | EXACT action (row (0,1,2,3), col (0,1,3,2)) maps T_10 into T+{0}; CITED |
| "+7 -> T42" | paper case | {1,2,4,6,7,8} = paper T_43 (stage 1) | EXACT; CITED |
| "+16 -> T65" | paper case | {2,4,5,6,9,12} = paper T_66 (stage 7) | EXACT; CITED |
| bounded-negative-part harmonic lemma | common lemma | final_referee_bounded_bad_harmonic_lemma | hypotheses compared (match); CITED |
| one-circle theorem | common lemma | paper Thm C.1 (a spectrum with a repeated eigenvalue has its four points on a generalized circle) | CITED |
| shared-zero Birkhoff theorem | common lemma | paper Thm B.1 / Cor. B.3 (a_i + b_j = 0) | CITED |
| complete facet list | core | paper Prop. 4.1 / Lemma 4.2 | not even needed: completeness is proved directly here (below) |

All cited six-element cases are earlier than T_50 (stage 8) in Table 6; all seven-element sets are stage 1.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | T = paper T_50; representative action | identity | EXACT | lib `g_act_set` with the document's convention pi -> column o pi o row |
| 2 | contact rank 6; f0,f1,f2,f3,f8,f10,f27,f59 a basis of the contact kernel | identity | EXACT | rank of P rows at T; rank of the 8 profiles; all vanish on T |
| 3 | eleven gap formulas g0 = x59, ..., g22 = x0 + x1 | identity | EXACT | values of f_j at those indices |
| 4 | f4 = f1+f2+f3-f0, f30 = f0-f3+f8+f27 | identity | EXACT | vector identities in V_R |
| 5 | intervals [-S,D], [-x3,x0] meet (S = min(x1,x2), D = min(x8,x27)) | inequality | PROOF-READ | S,D >= 0 from g16,g13,g8,g4; four endpoint comparisons from g22,g19 (x0+S), g14,g12 (x3+D), g20 (x0+x3) |
| 6 | decomposition (1) reproduces c and is nonnegative | identity + inequality | EXACT (identity, both signs of zeta) + PROOF-READ (signs) | symbolic in x, zeta; the generator list equals lib `generators(T)` |
| 7 | L_O >= d f30(O), d <= g4, d <= g15, (3) | inequality | PROOF-READ + EXACT instances | minors have Phi >= 0 by (4.2); f30(4) = f30(15) = 1, all generators >= 0; 30 exact random instances of (3) |
| 8 | extra tie with f30 > 0 forces d = 0; other ties exactly {0,7,12,13,14,16,20} | case split | EXACT | zero set of f30; for ties with f30 > 0 every generator of N_{T+j} is a minor |
| 9 | ties 0,7,12,13,14,16,20 map to T_10, T_43, B, F, A, T_66, G | citation | EXACT actions, CITED results | literal row/column actions checked on all 24 entries |
| 10 | degenerate spectra (repeated, cross collision) | citation | CITED | Appendix C (Thm C.1), Appendix B |
| 11 | E, rE, sE, uE raw ratios | identity | EXACT | sigma-points with b = (t,0,1,inf) |
| 12 | inverse (5), denominators: p + q = 1 + W; both vanish together | identity + PROOF-READ | EXACT / PROOF-READ | symbolic; 20 random rational configurations recovered exactly |
| 13 | homogeneous inverse P, Da, Dt; W Dt + Da = (W+1)(r-W)X; E Da Dt = P (W+1)(r-W) X W Z | identity | EXACT | sympy |
| 14 | contact levels H1,H3,H5,H11,H10; (8) H6 - h = -C/((r-W)(W+1)); (9) H4, H15, H0; height F1/E | identity | EXACT | sympy, from raw products |
| 15 | only W = 0 needs a second bound (H4 bounded elsewhere incl. W = inf) | analytic | PROOF-READ | H4 = r(W+1)/W |
| 16 | C(0,Z) = r(h-r-s-1); all monomials have deg_W - deg_Z >= 0; weight-zero part (11); single branch V0 | identity + analytic | EXACT + PROOF-READ | C is linear in Z, so the curve is the graph Z = (W+1)Q_W/(W(r-W)(r+s-u-W)) and has exactly one branch over W = 0 |
| 17 | (12) H15 = -su/(r(r+s)) W + O(W^2) on that branch; denominator r^2(r+s)(r+s+1-h) nonzero under (10) | identity | EXACT | H15 restricted to the graph, expanded at W = 0 |
| 18 | generic conditions jointly proper, witness (2,3,4,5) | identity | EXACT | |
| 19 | nonconstant coordinate function on every component (H4 if W varies, H0 = (W-r)(Z-1) if W is constant) | analytic | PROOF-READ | W != r at a realizable point |
| 20 | bounded-negative-part lemma excludes a generic strictly supported negative component | analytic / citation | PROOF-READ (application) + CITED (lemma) | see note (a) |
| 21 | wall restoration by Rouche along a complex line | analytic | PROOF-READ, argument supplemented by auditor | see note (b) |

Note (a), puncture set (auditor's clarification). The document lists punctures only as "poles, zero spinors and collisions". The lemma needs hypothesis 1 only at regular eight-distinct points, where the boundary cases apply literally; every non-regular point of the component (poles of any H_i, collisions, points where E or the height F_1/E degenerate, W = 0, W = infinity) can be put into the finite puncture set, because d <= g4 = -eta Im H4 is bounded near every point with W != 0 (including W = infinity) and d <= g15 is bounded near W = 0 by step 17. So the curve argument does not need the collision theorems at boundary points; they are needed only for counterexamples whose own spectra are degenerate. The component is not contained in the degenerate locus (it contains a regular point), so that locus meets it in finitely many points.

Note (b), wall restoration. Given a strict counterexample at parameters p0 = (r,s,u,h) violating (10), with point (W0,Z0): C(.,.;p0) is not the zero polynomial (its Z-coefficient W(r-W)(r+s-u-W) is a nonzero cubic), so some complex line through (W0,Z0) meets {C = 0} in an isolated zero; by Hurwitz/Rouche, for p near p0 the line carries a zero of C(.,.;p) near (W0,Z0). The inverse (5) is regular there (all denominators nonzero at the realizable eight-distinct point), the contact levels equal the perturbed real r,s,u,h identically (step 14), eta is fixed, so the perturbed point is in X with contacts T, strictly positive other gaps and Phi < 0 by continuity. The generic set is the complement of finitely many proper real hypersurfaces in R^4, hence dense, so p can be chosen generic. This is the document's argument with the quantifiers made explicit; it is correct.

## D. Falsification

- EXACT finite instances: (3) holds on 30 random rational points of N_T times random rational SO(4) matrices.
- NUMERIC: 6 points on the T_50 contact stratum (seeds 0-5), min over SO(4) by 12-restart BFGS: worst value -1.0e-15 (zero up to rounding; the minimum 0 is attained at permutation matrices of contacts). No violation.

## E. Verdict

**confirmed-modulo-citations.** Every algebraic and combinatorial claim was replayed exactly and is correct; the analytic argument is complete once the puncture set is taken as in note (a). Citations: hypotheses F, G; cases T_10, T_43, T_66, A = T_9^(7), B = T_8^(7); the bounded-negative-part lemma; Appendices B and C.

Presentation remarks (not gaps): the puncture set is described loosely (note a); the companion checker is cited for identities that the PDF states (all re-derived here).
