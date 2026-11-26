Evidence level: EXACT for every algebraic claim (cone, ties, inverse, curve, fibre x = 0, branch expansions, four-factor identity and factor limits); PROOF-READ for the curve, compactness and perturbation argument (one step supplied by the auditor); CITED for earlier cases, Table 5 set T7_21, Corollary B.3, Appendix B, the bounded-negative-part lemma and Theorem 1.2; NUMERIC for a brief main-statement search.

# Audit: root_inventory80_native_uniform_20261023 (paper case T_81)

PDF: `proofs/manin_transfer_20261021/root_inventory80_native_uniform_20261023.pdf`, 4 pages.
Replay: `replay/root_inventory80_native_uniform_20261023__check.py`, output in `..._check_output.txt` (all EXACT checks pass).

## Deviations and limits

- None from the assigned scope.
- The "companion exact coefficient array" and the verifier mentioned on p. 2 were not used. C(x,t) was reconstructed independently, as the unique 74-term irreducible factor of the numerator of H_13 - r computed from my own 24 products.
- The equation numbers in the PDF repeat ((1), (2), (4), (5), (6) are each used twice). I cite by page and content.
- Main-statement numerics are brief (4 points found).

## A. Statement

**Theorem 1 (p. 1).** Let T = {1,2,10,11,12,13} (equal to the paper's T_81; checked). For every x in X with P c(x) >= 0 and (P c(x))_pi = 0 for pi in T (reference pi* = 12), and every real orthogonal O, Phi_{c(x)}(O) >= 0. Extra ties, zero coefficients and boundary configurations are included. The unitary form follows from Theorem 1.2. The contacts are non-strict, and all spectra are covered (collisions go to Appendix B, p. 1; parameter walls to p. 4). This is at least as strong as (*) for T_81.

Internal claims:
- (1)-(2), p. 1: N_T is simplicial with generators f_0, f_1, f_8, f_10, f_12, f_15, f_24, f_53. Rows 21,17,15,7,6,4,3,0 form the identity, so the coefficients are gaps. The only non-minor generator is f_24 = h_{4312,+}, with coefficient d = g_3 and -L_O <= M d.
- (3): the ties not annihilating d are 0,4,6,7,15,17,21.
- Normalization A = (0,1,z,w), C = (oo,v,x,t), with E = F_12 - F_2 = (z-v)(w-t), F_12/E = x, and the real parameters H_13 = r, H_10 = s, H_1 = u, H_11 = s + kr. There is an explicit inverse (x,t) -> (v,w,z). The curve is C(x,t), with 74 terms and bidegree (5,2). On the curve, H_3 = r - 1 - r/x.
- x = 0 fibre: C(0,t) = rs t(1-t), and two smooth branches. Branch (0,0) has primitive limit A' = (0,oo,1,1), C' = (oo,0,1,1), and there only [A'_3 C'_3] of E vanishes. Branch (0,1) has limit A = (0,1,1,s/(s+kr)), C = (oo,1,0,1), and there only [A_2 C_1] vanishes.

## B. Dependencies

| dependency | kind | resolved to | check |
|---|---|---|---|
| ties that kill d (3,5,8,9,14,16,18,19,20,22,23) | Lemma 4.2 / (4.2) | none needed | EXACT |
| tie 0 | earlier case | T_46 (stage 2), T_76 (6) | EXACT |
| tie 4 | earlier case | T_46, T_36, T_38, T_37 | EXACT |
| tie 6 | earlier case | T_34 (1), T_58 (6) | EXACT |
| tie 7 | earlier case | T_37, T_46, T_36, T_38 | EXACT |
| tie 15 | earlier case | T_46, T_76 | EXACT |
| tie 17 | earlier case | T_58, T_34 | EXACT |
| tie 21 | seven-contact | H = {0,1,6,9,16,17,22} = Table 5 set T7_21 (stage 1), with explicit g; no earlier six-set exists | EXACT |
| T_81 is in stage 7, and every reduction above is to an earlier stage | Appendix G | | EXACT |
| branch punctures (0,0), (0,1) | Corollary B.3(b) ("single-divisor compensation") | foundations/appendix-b.md | CITED; its use checked below |
| collisions at interior points | Theorem B.1, Lemma B.2 | | CITED |
| compact pasting | Lemma 3.5 | | CITED |
| unitary form | Theorem 1.2 | | CITED |

**Use of Corollary B.3 (checked as the brief requires).**
- (i) The gaps are g_sigma = Re(i eta (Pi_sigma - Pi_12)/E), with constant eta. E = Pi_12 - Pi_2 = -[A_0A_1][C_0C_2][A_2C_1][A_3C_3] for completely free homogeneous lifts (EXACT). E has multidegree one in each of the eight spinors, so the same holomorphic lifts are used in numerator and denominator. The common projective change q -> q/x on branch 1 multiplies every bracket by the same factor.
- (ii) Holomorphic primitive lifts exist, and R(0) != 0. On branch 1, t = x + O(x^2) and the orders and leads are VD 0/1, WD 0/r, WN 1/r, ZN 1/rs, ZD 0/rs, VN 2/k (EXACT, from the exact power series of the branch). The factor limits are [A'_0A'_1] = -1, [C'_0C'_2] = 1, [A'_2C'_1] = rs, [A'_3C'_3] = 0. On branch 2, t = 1 - (uk^2/s)x^2 + ...; VN and VD have order 1 with lead k; WN and WD have order 2 with leads uk and uk(s+kr)/s; ZN and ZD have order 2 with equal leads. The factor limits are -1, 1, 0 and -k^2 ru/s. Under (11), R(0) is a nonzero multiple of rs (branch 1) or of k^2 ru/s (branch 2).
- (iii) No puncture has two vanishing E-factors. These two points are the only points over x = 0: there is no point at t = oo, since [t^2]C(0,t) = -rs != 0. Both are smooth points, since C_t = ±rs != 0. Further coincidences among brackets outside E (A'_2 = A'_3 = C'_2 = C'_3, and A_1 = A_2 = C_1 = C_3) are allowed by B.3, because (B.9) is an ambient identity.

## C. Steps

| # | claim (location) | kind | status | how checked |
|---|---|---|---|---|
| 1 | T = paper T_81 | finite | EXACT | |
| 2 | (1) generator list; identity rows; d = g_3 | finite | EXACT | library values |
| 3 | (2) -L_O <= M d | inequality | PROOF-READ | the minor generators have Phi >= 0 |
| 4 | (3) tie list and reductions | finite | EXACT | orbit and stage computation |
| 5 | E, F_12 = xE, F_13 - F_12, F_11 - F_10 (p. 1, (1), (4)) | identity | EXACT | flint |
| 6 | inverse (2)-(5), p. 2: v = VN/VD, w = WN/WD, z = ZN/ZD | identity | EXACT | sympy. Parameters r, s, u, k substituted as functions of the original variables |
| 7 | raw identities (6), (7), p. 2 | identity | EXACT | identically zero in the free original variables |
| 8 | WN - xWD = x(x+u)(VN - xVD); VN - VD = (t-1)(1+kx); v_t = kx(1+kx)/VD^2 | identity | EXACT / PROOF-READ (derivative) | |
| 9 | nonvanishing of VD, WD, ZD, 1+kx, x+r, x+s, x+u, x+s+kr at eight-distinct points | identity | EXACT | each factors into brackets of distinct points |
| 10 | (9) numerator of H_13 - r is -r(x+r)C; C has 74 terms, bidegree (5,2) | identity | EXACT | actual factorization: -(t-1) x (VD - kt) r (x+r) C. The extra factors are nonzero at eight-distinct points: VD - kt = -t(t-1)(x-v)/(x(t-v)). Presentation only |
| 11 | (5) all six contacts hold modulo C | identity | EXACT | flint division |
| 12 | (10) H_3 = r - 1 - r/x on the curve; raw identity | identity | EXACT | |
| 13 | d bounded at every puncture with x != 0, including x = oo | analytic | PROOF-READ | d = -eta Im H_3 |
| 14 | (12) C(0,t) = rs t(1-t); no point at t = oo; both points simple | identity | EXACT | |
| 15 | (13), (16) local derivatives and series | identity | EXACT | exact series of t(x) to order 3 and 4 |
| 16 | (14)-(17): primitive lifts, limits, four-factor identity (15), single vanishing factor at each branch | identity | EXACT | see B |
| 17 | branch punctures bounded by Corollary B.3 | citation | CITED | hypotheses (i)-(iii) verified |
| 18 | nonconstancy: all H_i constant forces x, v, t constant, so a point | algebraic | PROOF-READ | uses (10), F_6/(F_6 - F_0) = v (EXACT) and v_t != 0 |
| 19 | compact pasting excludes a negative region | citation | CITED (Lemma 3.5) | |
| 20 | Rouche restoration of the walls of (11) when C_p is not identically zero | analytic | PROOF-READ (argument supplied by auditor) | as in the other cases: Hurwitz on a line, (11) is a proper real-algebraic condition, the inverse is continuous because steps 6 and 9 hold at every genuine point, and strict conditions persist. The counterexample lies on C = 0 because the other factors in step 10 are units |
| 21 | C_p identically zero: fix x, vary t on P^1; d constant | analytic | PROOF-READ | the contacts hold identically, since step 11 gives a numerator equal to C times something |

## D. Falsification

NUMERIC, brief: `random_point_on_contact_stratum(T, margin=1e-3)` with seeds 8100-8105 found 4 points, all with contact set exactly T. With 8 BFGS restarts, the worst scaled min_O Phi is +7e-14. No violation. Every intermediate algebraic claim was replayed exactly. For the behaviour of Corollary B.3 near single-factor germs, see the numeric tests in foundations/appendix-b.md.

## E. Verdict

**confirmed-modulo-citations.** No gap found. The citations are: the earlier cases T_34, T_36, T_37, T_38, T_46, T_58, T_76; Table 5 set T7_21; Corollary B.3 (with hypotheses (i)-(iii) verified here); Theorem B.1 and Lemma B.2; Lemma 3.5; Theorem 1.2.

Presentation remarks:
- (9) omits the unit factors (t-1), x and VD - kt.
- Equation numbers are reused.
- On p. 1, "the 67 preceding cases" are not named. Valid earlier-stage reductions exist (B).
