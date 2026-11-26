Evidence level: EXACT for every algebraic claim (cone, barriers, ties and their labels, inverse, 355-term curve, all four pole branches, the constant-germ classification); PROOF-READ for the curve, compactness and perturbation argument (one step supplied by the auditor); CITED for earlier cases, Appendix B, the bounded-negative-part lemma, Lemma 3.4 / (3.11) (Birkhoff on E) and Theorem 1.2; NUMERIC for a brief main-statement search.

# Audit: final_referee_inventory81_native_uniform_20261023 (paper case T_82)

PDF: `proofs/manin_transfer_20261021/final_referee_inventory81_native_uniform_20261023.pdf`, 5 pages.
Replay: `replay/final_referee_inventory81_native_uniform_20261023__check.py`, output in `..._check_output.txt` (all EXACT checks pass).

## Deviations and limits

- None from the assigned scope. The "accompanying verifier" (p. 5) was not run or read.
- The certificate labels T33, T59, T30 and T78 on p. 1 are inventory indices, not paper indices. I resolved each as paper index = label + 1 (T_34, T_60, T_31, T_79) and checked that the case maps into T∪{j} and precedes T_82 in Table 6.
- Main-statement numerics are brief (3 points found).

## A. Statement

**Theorem 1 (p. 1).** Let T = {2,5,9,10,12,13} (equal to the paper's T_82; checked). For every x in X with P c(x) >= 0 and (P c(x))_pi = 0 on T (reference pi* = 12), and every real orthogonal O, Phi_{c(x)}(O) >= 0. The unitary form follows from Theorem 1.2. The contacts are non-strict. Collisions go to Appendix B (p. 1), and the parameter walls, including identically vanishing C, are treated on pp. 4-5. This is at least as strong as (*) for T_82.

Internal claims:
- (1): N_T is simplicial with generators f_0, f_1, f_12, f_39, f_45, f_67, f_70, f_72, and coefficients (q,p,j,c,e,b,d,a) = (g_18,g_16,g_8,g_4,g_3,g_14,g_11,g_7). The only non-minor generator is f_70 = h_{4231,+}, with coefficient d = g_11.
- (2): the barriers g_6 = j+d+a, g_15 = c+b+d and g_19 = q+j+d, each >= d.
- Ties: the zeros of f_70 outside T are 3,4,7,8,14,16,18, and these reduce to earlier cases. Every other tie forces d = 0.
- Normalization A = (0,1,z,w), C = (oo,v,x,t) with real H_13 = r, H_5 = s, H_9 = u, H_10 = h. Inverse (4) for (t,w,v). Curve C(x,z): 355 terms, bidegree (6,2). On the curve, H_11 = x(x+s)(x+u)/((x-1)(x+r)) - x, so d is bounded away from x in {1, -r}.
- Four branches: (1,1), (1,oo), (-r,-r), (-r,oo). On each, one barrier g_19, g_15 or g_6 has a real limit.
- Constant-germ classification: C is the zero polynomial only on parameter sets where H_11 is a real constant, so d = 0.

## B. Dependencies

| dependency | kind | resolved to | check |
|---|---|---|---|
| ties 0,1,6,11,15,17,19,20,21,22,23 force d = 0 | Lemma 4.2 | | EXACT |
| tie 3 | earlier case | T_34 (stage 1) [label T33]; also T_21, T_22, T_33 | EXACT |
| ties 4, 7, 16 | earlier case | T_60 (stage 9) [label T59]; also T_31, T_36, T_48, T_78 | EXACT |
| tie 8 | earlier case | T_31 (stage 1) [label T30] | EXACT |
| ties 14, 18 | earlier case | T_79 (stage 8) [label T78]; also T_30 (stage 9) | EXACT |
| T_82 is in stage 10, so all of the above precede it | Appendix G | | EXACT |
| collisions | Theorem B.1, Lemma B.2 | | CITED |
| constant-coordinate component: line x kappa + v_0 in X gives Hess C(kappa) v_0 = 0, so c lies in E, and Birkhoff | radial necessity (brief item 5), (3.11) (library-verified), Lemma 3.4 | | CITED / EXACT (library (3.11)) |
| compact pasting | Lemma 3.5 | | CITED |
| unitary form | Theorem 1.2 | | CITED |

## C. Steps

| # | claim (location) | kind | status | how checked |
|---|---|---|---|---|
| 1 | T = paper T_82 | finite | EXACT | |
| 2 | (1) generators, identity rows, coefficients equal gaps | finite | EXACT | |
| 3 | (2) barrier identities | finite | EXACT | the values of the generators at 6, 15, 19 |
| 4 | tie list, seven reductions, label resolution | finite | EXACT | orbit and stage computation |
| 5 | E = (z-v)(w-t), F_12 = xE, the three ratio formulas (p. 2) | identity | EXACT | sympy |
| 6 | inverse (4); Td - Tn = x(x+s)(1-z) | identity | EXACT | substituted r, s, u, h as functions of the original variables |
| 7 | nonvanishing of x, x-1, x+r, x+s, x+u, x+h, V at eight-distinct points | identity | EXACT | each factors into brackets of distinct points |
| 8 | numerator of H_13 - r is (x+r)C up to a constant; C has 355 terms, bidegree (6,2) | identity | EXACT | actual factorization: x (x-z)(x+u)(x+r) C. The extra factors are units at eight-distinct points. Presentation only |
| 9 | all six contacts hold modulo C | identity | EXACT | flint division |
| 10 | (5) H_11 formula; (9) H_19; (10) H_15; (11) H_6 | identity | EXACT | the first three modulo C, H_6 identically |
| 11 | (7), (8) sections; z-degree 2; four simple points | identity | EXACT | C_z != 0 at (1,1) and at (-r,-r); d_tau C~ != 0 at (1,oo) and at (-r,oo) |
| 12 | (-r,-r): H_19 = r(s+1)/(r+1); (-r,oo): H_19 -> s | identity / limit | EXACT | |
| 13 | (1,1): H_15 -> u+r+1 | limit | EXACT | computed for an arbitrary branch slope |
| 14 | (1,oo): tau/eps -> kappa; w -> (s+r+2)/(r+1); v -> w; H_6 -> v - 1 | limit | EXACT | |
| 15 | d bounded at every puncture: 0 <= d <= barrier on Omega | analytic | PROOF-READ | d = -eta Im H_11, which depends on x only. x = oo gives the limit s+u-r+1 (EXACT) |
| 16 | constant H_i on a component: radial necessity, so c is in E, and Birkhoff; constant x: (12) gives z, so a point | algebraic | PROOF-READ + CITED | (12) x F_0 = (x-z)F_14 EXACT |
| 17 | compact pasting at generic parameters (6) | citation | CITED (Lemma 3.5) | |
| 18 | C zero polynomial: (s+1)(u+1) = 0 and r(r-u)(r-s) = 0, then the case tree | identity | EXACT | the four specialized coefficients printed on p. 4 replayed. On every leaf, H_11 is a real constant (u, 0, s, 0, 0), so d = 0. In fact no leaf makes C identically zero, so the branch is vacuous (r = 0 is also impossible, since F_13 - F_12 = x(t-v)(z-w) != 0) |
| 19 | Rouche restoration of the walls of (6) | analytic | PROOF-READ (argument supplied by auditor) | as for T_80 and T_81. The counterexample lies on C = 0 because the other factors in step 8 are units |

## D. Falsification

NUMERIC, brief: seeds 8200-8205 gave 3 points with contact set exactly T. The worst scaled min_O Phi is +9.5e-13. No violation. Every intermediate claim (barriers, branch limits, classification) was replayed exactly.

## E. Verdict

**confirmed-modulo-citations.** No gap found. The citations are: the earlier cases T_21, T_22, T_30, T_31, T_33, T_34, T_36, T_48, T_60, T_78, T_79 (one per tie suffices); Appendix B; Lemma 3.4 with (3.11); Lemma 3.5; Theorem 1.2.

Presentation remarks:
- The numerator factorization omits the units x, x - z and x + u.
- The certificate labels on p. 1 use inventory numbering (paper index minus one).
- In the classification, "the possibilities r = 0" are vacuous, because r != 0 at every genuine point.
