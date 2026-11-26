# Audit: final_referee_shared_coreK_uniform_20261023 ("Three Lorentz-cone inequalities on a seven-contact source cone")

Evidence level: every identity, the derivation of the degree-nine polynomial G, both chart positivity tables (coefficient counts included), the zero-leaf identity, the symmetry claims and the resource allocations were replayed EXACT by our own code. The covering and exhaustion logic is PROOF-READ. One input is CITED, the "sharp SOC input inequality" for (f1, f9; f22); we have NUMERIC support for it.

PDF: `proofs/manin_transfer_20261021/final_referee_shared_coreK_uniform_20261023.pdf`, 3 pages. Replay: `replay/final_referee_shared_coreK_uniform_20261023__check.py`, output in `__check_output.txt`. On the last run 20 EXACT checks passed and none failed.

## Deviations and limits

- G is not displayed in the document; it says G is "stored in the certificate". We did not read that certificate. We re-derived G ourselves as the square root of Q^4 det H / (-4Q^2), and the result matches every property the document states: 746 terms, multidegree (3,2,2,3,2,3,9), d^9 coefficient +1, and coefficient counts of 5373 and 7297 in the two charts with an 8X monomial.
- The SOC atom (f1 = n12,34, f9 = n14,23; f22 = h3412,+), that x Phi_f1 + y Phi_f9 + delta Phi_f22 >= 0 whenever x, y, delta >= 0 and xy >= delta^2, is CITED and was not proved here. NUMERIC: the minimum over SO(4) of Phi_f22 + 2 sqrt(Phi_f1 Phi_f9) was -1.3e-28 over 40 Nelder-Mead restarts (seed 7). It is sharp, with minimum 0, and no violation was found.
- The unitary remark uses Theorem 1.2 (CITED). It is not needed for (*).

## A. Statement

Let K = {0, 1, 2, 9, 10, 14, 17}, the paper's Table 5 seven-set #16 (EXACT).

Theorem 1: let x be in X and pi* in K, and suppose P c(x) >= 0 with (P c(x))_pi = 0 on K; extra contacts are allowed. Then Phi_{c(x)}(O) >= 0 for every O in SO(4) (equivalently O(4)).

All spectra are covered: the necessity step uses only X contained in Sing C, the conic step is pointwise in c, and no genericity is used. The statement is stronger than (*) for K. As the document itself says, using it for a six-contact case requires a separate reduction to K; that reduction is not in this document.

## B. Dependencies

| dependency | kind | status |
|---|---|---|
| generators(K) = [0,2,12,13,15,21,26], simplicial, with f26 = h4132,- the only non-minor | paper App. D / lib | EXACT |
| X is contained in Sing C, and X is closed under conjugation | paper / lib | EXACT in lib; PROOF-READ |
| the group action preserves X and Phi (4.7) | paper / lib | EXACT in lib |
| SOC atom (f1, f9; f22) | supplement lemma ("sharp SOC input inequality") | CITED; NUMERIC support |
| Phi of a minor vector is >= 0 (4.2) | paper | EXACT/PROOF-READ in lib |
| Theorem 1.2 | paper | CITED (unitary remark only) |

## C. Steps

| # | claim | kind | status | how |
|---|---|---|---|---|
| 1 | rank P[K] = 7; the 7 profiles evaluated at rows [23,15,8,7,4,16,11] give I7, so the coefficients are gap readings and are >= 0; this is (1) | identity | EXACT | |
| 2 | (2): for each pair {i,j} in {2,12,15} there is a group element sending {f1, f9} to {fi, fj} and f22 to f26 | identity | EXACT | 8 elements per pair |
| 3 | the stabilizer of K has 6 elements and acts as S3 on (a,b,c) with the same permutation on (q,v,w), fixing d | identity | EXACT | also G is invariant under the induced (a b)(q v) and (b c)(v w) |
| 4 | (3): the balanced allocation gives non-minor total s/2 with resource sums exactly a, b, c; the dominant split gives sqrt(a(b+c)); reducing delta keeps xy >= delta^2; so Phi >= 0 whenever d <= H3 | identity + logic | EXACT (sympy) + PROOF-READ | also NUMERIC below |
| 5 | (5): C(g + z kappa) = C + zQ with C = (q+d)(v+d)(w+d) + abc and Q = ab+ac+bc+d^2 | identity | EXACT | |
| 6 | realizable implies F = 0 and det H = 0 (polarization, as in core J) | analytic | PROOF-READ | |
| 7 | (6): Q^4 det H(g - C kappa/Q) = -4 Q^2 G^2 with G a polynomial having 746 terms, degree 9, multidegree (3,2,2,3,2,3,9), d^9 coefficient 1 | identity | EXACT | 14x14 Bareiss with 819 exact divisions; the z-degree is 3; exact sqrt |
| 8 | (8): chart coverage of the closed balanced strict-failure region, with inverse (X,Y,Z) = (x,y,z0)/(1-x-y-z0) | identity + logic | EXACT + PROOF-READ | |
| 9 | (8): D^6 G(chart) has 5373 coefficients, all > 0, constant term 1, so G > 0 | identity | EXACT | |
| 10 | (9): dominant chart coverage of {0 < S < 1, S <= a < 1/S, c > 0}; inverse checked | identity + logic | EXACT + PROOF-READ | |
| 11 | (9): L^6 (1+Y)^2 (1+Z)^4 G(chart) has 7297 coefficients, all > 0, including the monomial 8X; with X > 0 this gives G > 0 | identity | EXACT | |
| 12 | (10): G(q,a,0,v,0,w,d) = {a(v+w+d)+(v+d)(w+d)}(v+d)^2(w+d)^2(q+d)^3 > 0 when d > 0 | identity | EXACT | |
| 13 | exhaustion: H3 >= 1 (SOC), else balanced (2m <= s), or dominant with S > 0, or S = 0; the symmetry from step 3 puts the dominant coordinate at a, and swapping b with c gives c > 0 | case split | PROOF-READ | |

## D. Falsification

- NUMERIC, SOC atom: 40 restarts gave a minimum of -1.3e-28, i.e. 0 (sharp).
- NUMERIC, conic claim (3) on 12 random (a,b,c), a third of them dominant, with d = H3(a,b,c) and 10 BFGS restarts each: the minimum of Phi/(a+b+c) was 1.8e-13, so no violation.
- No search was made on the realizable stratum. G > 0 is exact on the whole failure region.

## E. Verdict

**confirmed-modulo-citations.** Every computation was replayed exactly and the case split is exhaustive. The only CITED input is the SOC atom inequality for (f1, f9; f22); it is numerically sharp and showed no violation. No gaps were found.
