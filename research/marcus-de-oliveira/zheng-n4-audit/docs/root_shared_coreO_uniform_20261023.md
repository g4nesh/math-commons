# Audit: root_shared_coreO_uniform_20261023 (seven-contact set O, Table 5 no. 15)

Evidence level: decomposition, circuit, group actions, cubic, Hessian identity (7), chart positivity and the zero-leaf identity (10) are EXACT (own replay). Chart coverage and the conic addition (5) are PROOF-READ. The SOC reference inequality is CITED, and NUMERIC-tested at its sharp form. The unitary sentence is CITED to Theorem 1.2.

PDF: `proofs/manin_transfer_20261021/root_shared_coreO_uniform_20261023.pdf`, 3 pages. Replay: `replay/root_shared_coreO_uniform_20261023__check.py`, with output in `__check_output.txt`.

## Deviations and limits

- The document's certificate (`root_shared_coreO_uniform_exact_20261023.py`) was neither executed nor loaded. Every polynomial was recomputed from the paper's cubic (A.2), P and Appendix D through `lib/moc4lib.py`.
- The document does not say which exponents it uses to clear denominators in the charts ("full individual multidegrees"). I replayed both the minimal exponents and the full ones (D^6 for chart (8), L^6 (1+Y)^2 (1+Z)^4 for chart (9)). The full exponents reproduce the stated coefficient counts exactly, 2259 and 3096, and both choices give all-positive polynomials. Any positive multiplier preserves the sign, so the conclusion does not depend on this choice.
- I did not audit the SOC reference inequality for (f1, f9; f22). It is attributed to F03 (`final_referee_T2_scalar_compensation_20261022`). That document itself says that identifying its quartics with the atoms f1, f9, f22 needs a separate exact coordinate bridge, which I did not locate or check. I tested the atom inequality numerically only.

## A. Statements

Contact set O = {0,1,2,7,14,16,23}, which is paper Table 5 set no. 15 (EXACT, `AUX7[15]`). Here "O" names a set, not the orthogonal matrix.

**Theorem 1 (p. 1).** Let x lie in X and pi* be such that c = c(x) satisfies Pc >= 0 and (Pc)_sigma = 0 for every sigma in O. Further contacts, collisions and zero coefficients are allowed. Then Phi_c(M) >= 0 for every M in SO(4), and on U(4) by Theorem 1.2. The argument is purely algebraic on X: realizability is used only through "Re x is a singular point of the Hessian of C". It therefore covers every spectrum, degenerate ones included. `covers_all_spectra` = true.

Internal intermediate claims, each stronger than the theorem:
- (4) For {i,j} ⊂ {2,12,15} and all x, y >= 0 with xy >= delta^2, Phi_{x f_i + y f_j + delta f26}(M) >= 0 for every M in SO(4), with no realizability needed.
- (5) a f2 + c f12 + e f15 + d f26 is safe whenever 0 <= d <= H3(a,c,e).
- If d = 1, H3(a,c,e) < 1 and a, b, q, c, e, v >= 0, then G(a,b,q,c,e,v,1) > 0.

## B. Dependencies

| dependency | kind | used for | match |
|---|---|---|---|
| SOC reference inequality: x f1 + y f9 + delta f22 is safe whenever xy >= delta^2 | common-lemma (F03, final_referee_T2_scalar_compensation_20261022) | atoms (4) | EXACT: for each pair {i,j} ⊂ {2,12,15} a group element maps {f1,f9} to {f_i,f_j} and f22 to f26. Elements found: ((0,1,3,2),(0,2,1,3),0), ((0,3,2,1),(3,2,0,1),0), ((0,2,1,3),(1,2,3,0),0) in lib convention. The inequality itself is CITED, with the bridge caveat above. NUMERIC: min over SO(4) of Phi_{f26} + 2 sqrt(Phi_{f_i} Phi_{f_j}) is 0 to within 1e-12 for all three pairs, and also for the reference (f1,f9;f22). The inequality is sharp and was not violated. |
| X contained in Sing C | paper A.1 | the Hessian of a realizable lift is singular | EXACT in lib |
| G action preserves realizability and Phi (4.7) | paper Sect. 4.2 | the transposition in (3) and the six permutations of (a,c,e) | EXACT in lib |
| (4.2) Phi of a minor vector = det^2 | paper | d = 0; unused minor leaves | EXACT in lib |
| Theorem 1.2 | paper | unitary sentence only | CITED |

Prop. C.6 is not needed: this document uses only that the Hessian is singular, never its inertia.

## C. Steps

| # | claim (location) | kind | status | how checked |
|---|---|---|---|---|
| 1 | rank P = 14; rows at O have rank 7; the kernel is spanned by [f2,f5,f11,f12,f15,f23,f26], giving the unique decomposition (1) (p. 1) | identity | EXACT | lib |
| 2 | rows [15,10,17,8,4,9] = first six rows of I_7, so a,b,q,c,e,v >= 0; rows [13,20,12,22,18,3] read a+d, ..., v+d >= 0 (p. 1) | identity | EXACT | integer table |
| 3 | vanishing-facet list = [2,5,11,12,15,20,23,26]; circuit (2) f20 + f26 = f2+f5+f11+f12+f15+f23 (p. 1) | identity | EXACT | |
| 4 | d = 0 is safe (six minor vectors) | identity | EXACT | |
| 5 | d < 0: (3) has nonnegative coefficients; transposition (rho = rho' = id, eps = 1) stabilizes O, swaps f20 and f26, and acts by 2<->5, 11<->12, 15<->23. So the cost becomes a realizable form (1) with d' = -d > 0 (p. 1) | identity + group | EXACT | own check plus lib (4.7) |
| 6 | Six elements of Stab(O) fix f26, induce all permutations of (f2,f12,f15) and permute the six minor vectors (p. 2) | group | EXACT | \|Stab(O)\| = 12; all 6 permutations occur |
| 7 | Atoms (4) | citation | CITED (+ EXACT bridge, NUMERIC sharpness) | see B |
| 8 | Conic addition (5): the balanced allocation (atoms ac, ae, ce with x = y = delta = (a+c-e)/2, and so on) sums to minor coefficients exactly (a,c,e) and non-minor coefficient s/2. In the dominant case the atoms (ac/(c+e), c) and (ae/(c+e), e) carry delta_1 + delta_2 = sqrt(a(c+e)). Lowering delta and adding unused minor terms keeps Phi >= 0 (p. 2) | inequality | PROOF-READ | Algebra rechecked: x y = delta^2 holds in each atom. x f_i + y f_j + delta' f26 with delta' < delta is a convex combination of the atom and x f_i + y f_j. NUMERIC: 25 random (a,c,e), 5 of them dominant, gave worst normalized min Phi = +7e-16. |
| 9 | F(g + z kappa) = C + zQ with C = ace + (b+d)(q+d)(v+d) and Q as in (6); C > 0 when d > 0 (p. 2) | identity | EXACT | polynomial equality with the displayed C and Q |
| 10 | A realizable lift has singular H, so Q ≠ 0, z = -C/Q and det H = 0 (p. 2) | analytic | PROOF-READ | Same polarization/Euler argument as in coreN. Q ≠ 0 because F = C + zQ = 0 and C > 0. |
| 11 | (7): Q^4 det H(g - C kappa/Q) = -4 d^2 Q^2 G^2, where G has 415 terms, degree 8, individual degrees (2,2,2,2,2,2,8) and [d^8]G = 1 (p. 2) | identity | EXACT | 14 x 14 symbolic Bareiss; exact division by -4 d^2 Q^2; exact square root |
| 12 | Balanced region (triangle inequalities and s < 2), walls included, equals chart (8) with X, Y, Z >= 0; the inverse is as stated (p. 2) | case split | PROOF-READ + EXACT | sympy inverse; a+c+e = 2(X+Y+Z)/D < 2, and c+e-a = 2X/D (cyclically) |
| 13 | D^6 G on chart (8) at d = 1: 2259 coefficients, all positive, constant term 1, so G > 0 (p. 3) | inequality | EXACT | flint. The minimal exponent also gives all-positive coefficients. |
| 14 | Dominant region (a = max, 2a >= s, a(c+e) < 1, S = c+e > 0, e > 0 after a swap) equals chart (9) with X > 0, Y >= 0, Z >= 0. S^2 <= aS < 1 gives X > 0, and aS - S^2 = S(a-S) >= 0 gives Y >= 0 (p. 3) | case split | PROOF-READ + EXACT | sympy inverse |
| 15 | L^6 (1+Y)^2 (1+Z)^4 G on chart (9): 3096 coefficients, all positive, including 8X; constant term 0, so G > 0 because X > 0 (p. 3) | inequality | EXACT | flint |
| 16 | S = 0: G(a,b,q,0,0,v,d) = (b+d)(a^2 b + abd + ad^2 + bd^2 + d^3)(v+d)^2(q+d)^2 > 0 for d > 0 (10) (p. 3) | identity | EXACT | flint |
| 17 | Exhaustion: d = 0 is safe; d < 0 reduces to d > 0 by step 5; d > 0 with H3 >= d is safe by step 8; d > 0 with H3 < d is either balanced (step 13) or dominant with S > 0 (step 15) or S = 0 (step 16), and in each case G ≠ 0, so the cost is not realizable (p. 3) | case split | PROOF-READ | Each dominant subcase is brought to "a = max, e > 0" by the actions of step 6, which keep d and keep the form (1). Positivity in steps 13/15 holds for all b, q, v >= 0, so the permutation of (b,q,v) does not matter. |
| 18 | Unitary extension | citation | CITED | Theorem 1.2 |

## D. Falsification (NUMERIC)

- SOC atoms and the H3 bound: see steps 7 and 8. Both are sharp (minimum 0) and neither was violated.
- Realizable points with contacts containing O (lib least squares, seeds 0..19, 20 points found). After the sign normalization of step 5, every point satisfies H3(a,c,e) >= d up to rounding. Two points have a genuinely nonzero non-minor coefficient (d/mean gap = 0.037 and -0.344) with H3 - d > 0, so the safe region (5) is actually used. A first version of this test wrongly read the main coefficients after transposition as (a+d, c+d, e+d). The correct images are (b+d, q+d, v+d), because the transposition sends f5 to f2, f11 to f12 and f23 to f15. The corrected test is the one in the replay.
- No violation of the main inequality was searched beyond this, since the main statement is itself a special case of the conjecture.

## E. Verdict

**confirmed-modulo-citations.** Citations: the sharp SOC reference inequality for (f1, f9; f22). It is attributed to F03, and its link to that document's Theorem 2 needs a coordinate bridge that I did not audit. The group bridge from (f1,f9;f22) to the three atoms used here is EXACT. Theorem 1.2 is cited for the unitary sentence. No gaps found in this document.
