# Audit: independent_inventory40_native_uniform_20261023 (primary proof of paper case T_41)

Evidence level: EXACT for every finite and algebraic claim (contact set, generators, bound (2), lists (4) and (5), all twelve reductions of an extra contact, the parametrization (6)-(12), both fibres over v = 0 and their bracket limits, (17), (19), (20)); PROOF-READ for the curve argument (punctures, constant case, wall perturbation); CITED for Corollary B.3, Theorem B.1 / Lemma B.2, Lemma 3.5 (bounded negative part / maximum principle), Prop. 4.1 / Lemma 4.2, Theorem 1.2, and the earlier contact cases.

PDF: proofs/manin_transfer_20261021/independent_inventory40_native_uniform_20261023.pdf (5 pages). Auditor label C16-T40-T42.

## Deviations and limits

- The supplement's "inventory 40" is the paper's T_41 = {0, 1, 3, 6, 11, 12} (identical set).
- The document uses homogeneous "spinors" A_i, C_j in C^2 with F_pi = prod_i [A_i, C_pi(i)]. In the paper's chart this is v_i = A_i = (a_i, 1), v_{4+j} = C_j = (-b_j, 1), so [A_i, C_j] = a_i + b_j. I re-derived every identity from that definition with my own sympy code; nothing was taken from the supplement's scripts.
- The "69 preceding contact cases" are interpreted as the paper's Table 6 stages 1 to 8 (T_41 is in stage 9). I checked, for each extra contact, that some six-subset of T u {j} is either an all-minor orbit or the orbit of a case of stage <= 8. I did not count the list to 69 (stages 1-8 contain 71 cases).
- No numerical search on the main inequality (it is a special case of the conjecture).

## A. Statements

Theorem 1 (p. 1). For every x in X with P c(x) >= 0 and (P c(x))_pi = 0 for pi in T_41 = {0,1,3,6,11,12} (further contacts allowed), Phi_{c(x)}(O) >= 0 for every O in O(4) ("real orthogonal response"); by Theorem 1.2 (two real responses) the same holds for unitary responses. Spectra: all, including repeated eigenvalues and a_i + b_j = 0 ("collisions"), which are sent to Appendix B (Theorem B.1, Lemma B.2). This is at least as strong as (*) for T_41.

## B. Dependencies

| dependency | kind | used for | status |
|---|---|---|---|
| Prop. 4.1 / Lemma 4.2 | core | rays (1) of N_T | CITED; ray list EXACT |
| eq. (4.2) | core | minor rays have Phi >= 0 | CITED (library replays) |
| Theorem B.1, Lemma B.2 | common-lemma (Appendix B, audited by the project) | any counterexample with a repeated eigenvalue or a_i + b_j = 0; ordinary collision points on the curve | CITED |
| Corollary B.3 | common-lemma | the two punctures over v = 0 | CITED; the three hypotheses were checked here (see C, steps 14-15) |
| Lemma 3.5 (bounded negative part / maximum principle on a compact curve) | core | compact exclusion | CITED |
| Theorem 1.2 | core | unitary version | CITED |
| earlier contact cases for the 12 extra contacts of list (5) | paper-case | ordinary additional-tie boundary | CITED; group reductions EXACT: 2 -> T_5, 5 -> T_6, 7 -> T_19, 9 -> T_43, 10 -> T_44, 13 -> T_44, 14 -> T_37, 16 -> T_58, 18 -> T_57, 20 -> T_40, 21 -> T_5, 23 -> T_58 (stages 1, 1, 1, 1, 1, 1, 1, 6, 6, 1, 1, 6; all < 9) |

## C. Steps

| # | claim (location) | kind | status | how checked |
|---|---|---|---|---|
| 1 | T = paper T_41 | finite | EXACT | library |
| 2 | rays (1) = 0,1,2,7,8,9,10,15,18,29,30; non-minor exactly 18, 30; span dim 8; contact rows rank 6 (p. 1) | finite | EXACT | `generators`, `rank` |
| 3 | (2): theta = lam18 + lam30 <= g17 for any nonnegative representation | inequality | EXACT | f18 and f30 have value 1 at pi_17; all other rays have value >= 0 there |
| 4 | (3): -L_O <= M theta <= M g17 for fixed O | inequality | PROOF-READ | Phi_c = sum lam_k Phi_{f_k}, minor terms >= 0 |
| 5 | list (4): extra contact at 4,8,15,17,19,22 forces lam18 = lam30 = 0 | finite | EXACT | these are exactly the noncontact indices where f18 > 0 and f30 > 0 |
| 6 | list (5): each remaining extra contact reduces to an earlier case | citation | EXACT (group) + CITED | see table B |
| 7 | E = F6 - F0 = (z - x)(w - t); F0/E = v - 1 (p. 2) | identity | EXACT | sympy, from the spinor definition (E = (t-w)(x-z), same thing) |
| 8 | (6), (7): definitions of r, s, u, h; levels; r, s, u, h are nonzero at eight-distinct points (factor into brackets) | identity | EXACT | factorizations printed by the replay |
| 9 | (8) three edge equations; (9) the inverse x, t, w | identity | EXACT | substituted back, all three equations hold identically |
| 10 | inverse denominators nonzero at eight-distinct points (h + z = 0 forces v = z = -h; u + v = 0 forces t = 1) | algebra | PROOF-READ | reread from (9); both are collisions |
| 11 | (10): H1 = r is equivalent to the 28-term C(v,z) = 0, bidegree (2,2) | identity | EXACT | the numerator of r(v,z) - r after cancellation is exactly C(v,z); the remaining denominator is h (z - v) (z(uh + sv) - uh), i.e. collisions a2 = c1 and t = 1 |
| 12 | (11): H17 = 1 + u + h + uh/v | identity | EXACT | raw identity before the last cut |
| 13 | punctures with v != 0 (including v = infinity, collisions, E = 0, inverse poles): negative part bounded by (3) and (11) | analytic | PROOF-READ | u, h are constant on the curve, so g17 = -eta Im H17 is bounded wherever v stays away from 0 |
| 14 | (12), (14): C(v,0) and C(0,z); fibre over v = 0 is {z = 0, z = z*} with z* = h(r-1)/Delta, both simple; no point at z = infinity | identity | EXACT | C_z(0,0) = h^2 u (r-1), C_z(0,z*) = -h^2 u (r-1), both nonzero under (13); so v is a local parameter on each branch |
| 15 | first branch: slope alpha (15); limits x/v, t/v -> alpha, w/v -> alpha + 1/u; rescaled primitive lifts (16); E' factors -> (-1, 1, 0, 1/u) | identity | EXACT | sympy with an unknown second-order coefficient beta |
| 16 | second branch: limits (18); E factors -> (-1, 1, z*^2, 0) | identity | EXACT | sympy |
| 17 | (17): E = -[A0A1][C0C1][A2C2][A3C3] for free spinors | identity | EXACT | polynomial identity in 16 variables |
| 18 | Corollary B.3 applies at both branches | citation | CITED, hypotheses checked | (1) gaps g_i = -eta Im H_i = Re(i eta (Pi_i - Pi_0)/E) with the same lifts in numerator and E (all are multilinear, so rescaling lifts multiplies everything by one common factor); (2) E = d R with d = [A2'C2'] (branch 1) or [A3C3] (branch 2) and R(0) = (-1)(1)(1/u) resp. (-1)(1)(z*^2), nonzero since u, z* != 0; the lifts are holomorphic in v (z/v, x/v, t/v, w/v holomorphic); (3) exactly one factor of E vanishes at each branch, so no two-factor puncture arises over v = 0; at all other punctures the argument of step 13 is used instead |
| 19 | compact exclusion: pasted -L_O subharmonic, bounded above at all punctures, hence constant (Lemma 3.5); positive constant forces all H_i constant | analytic | PROOF-READ + CITED (Lemma 3.5) | boundary points of the negative component are either L = 0, a new contact (steps 5-6) or a collision (Appendix B) |
| 20 | constant case: (11) gives v constant; (19) x(F12 - F2) = F12 gives x constant; (20) z = xh/(h + v - x) gives z constant; contradiction | identity + logic | EXACT + PROOF-READ | (19) is a polynomial identity of products (F12 - F2 = (t - w)(z - v) up to sign, nonzero at eight-distinct points) |
| 21 | walls r = 1, Delta = 0 restored by Rouche: a counterexample at a wall parameter perturbs to a counterexample at generic real (r,s,u,h) | analytic | PROOF-READ | the perturbed point is an honest eight-point configuration with H1, H3, H11, H12 real (so the six contacts hold), and strict support, negativity and eight-distinctness are open conditions; C is a nonzero polynomial by (12), so a generic complex line through the point meets the curve in an isolated zero |

## D. Falsification

All intermediate identities were replayed exactly, so no numerical test was needed for them. The inequality steps (2)-(3) are exact finite facts. No numerical search on the main statement.

## E. Verdict

confirmed-modulo-citations. No gap found. Citations: Prop. 4.1/Lemma 4.2, (4.2), Theorem B.1, Lemma B.2, Corollary B.3 (hypotheses checked), Lemma 3.5, Theorem 1.2, and the earlier cases T_5, T_6, T_19, T_37, T_40, T_43, T_44 (stage 1) and T_57, T_58 (stage 6) for the extra contacts.

Minor presentation points (not gaps): (a) equation (10) is labelled "(1)" in the PDF for the D, B line; (b) the statement "with no discarded factor" is true only after cancelling a factor v that appears in both the numerator and the denominator of r(v,z) - r; (c) the inverse-chart exceptional points (z = v and z(uh + sv) = uh) are collisions, and the document correctly adds them to the puncture set, where step 13 applies because v != 0 there (at v = 0 the first is the branch point (0,0); the second does not meet the curve, since C(0,1) = -uhs != 0).

Replay: `replay/independent_inventory40_native_uniform_20261023__check.py`; output `replay/independent_inventory40_native_uniform_20261023__check_output.txt` (all EXACT checks pass).
