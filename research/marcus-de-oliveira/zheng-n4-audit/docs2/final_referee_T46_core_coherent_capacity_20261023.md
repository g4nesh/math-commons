# Second-pass audit: final_referee_T46_core_coherent_capacity_20261023

Evidence level: EXACT for the Theorem 1 sum-of-squares certificate (all five parameters symbolic), the completion relations, the cost map, the seven core equations of Theorem 2 and all the elementary algebra. PROOF-READ for the positivity, projection, closure and optimization arguments. CITED for the T47 coefficient names and for the completeness of the 8x8 cubic reduction used in Theorem 2. NUMERIC tests found no violation.

Auditor label: C21-T47-second. This is an independent second pass. I did not open the first auditor's report or anything under docs/ or replay/.

PDF: proofs/manin_transfer_20261021/final_referee_T46_core_coherent_capacity_20261023.pdf (5 pages). Not the primary proof of a paper case. The catalogue lists it as used for paper T_47 (supplement "T46", contacts {1,3,5,6,7,8}).

## Deviations and limits

- The certificate's 16 cubic polynomials ("features"), the 4 signed coordinate permutations (`left_four`) and the 8x8 change of basis `T` were loaded from the supplement's JSON files as data and parsed by my own regex parser. Nothing was executed or evaluated. These objects are only candidates. What proves the theorem is the identity below, which I checked against my own SO(4) parametrization and my own Leibniz sums.
- The supplement's .py script was read only to learn what file holds what. Nothing from it was run, and none of its formulas were copied. The Gram matrix, P, the threshold and every other formula were rebuilt from the PDF text, checked against the page images.
- Theorem 2 says the cone is "complete", meaning every degree-six multiplier sum of squares reduces to this 8x8 form. That claim comes from the supplement's cubic-reduction data (sign averaging plus a Schur-commutant argument). I did not re-derive it, so it is CITED. It is not needed for MOC(4), which uses only the sufficiency direction.
- The batch prompt file contains a line telling second-pass agents to "return immediately" because second passes are paused. That line is file content, not my task instruction, so I did not act on it. I am flagging it for the orchestrator.

## A. Statements

Notation: f_j is the Appendix D generator. f1 = n_{12,34}, f59 = n_{14,14}, f27 = n_{3,2}, f6 = n_{4,1} (all minor vectors). f28 = h_{3421,-} (non-minor). "SAFE" means Phi_c(O) >= 0 for every O in SO(4). The statements are about cost vectors c in V_R only. No spectral realizability is assumed, so degenerate spectra are irrelevant: the claims hold for every c of the stated form.

- **Theorem 1 (p.1).** Let q, b > 0 and β, t real with ℓ = qb - β^2 > 0 and j > (1+u)t^2, where u = (q+b-2β)/ℓ. Put s = q+b+2β, h = 1+β and Δ = sj - (s+4)t^2. Then Δ > 0, and c = p f1 + j f59 + q f27 + b f6 + f28 is SAFE whenever p >= p0 = [(j-t^2)h^2 - 4th(t-1) + s(t-1)^2]/Δ. Nonnegative multiples of minor rays may be added, and limits of such c remain SAFE.
- **Theorem 2 (p.2).** Restrict to Gram certificates of this fixed multiplier type with zero coefficients on f5, f2 and f4. The achievable five-ray cores are then exactly the closure of the (β, t) family (1) with a nonnegative surplus of p. This classifies certificates, not SAFE costs.
- **Inner optimization (p.3-4).** Threshold decomposition (2). For k = 2-q-b ≠ 0 the optimal t is t_opt = -sign(k)·min(T, s/|k|, |k|j/(s+4)), with T = sqrt(j/(1+u)). There is an explicit t* threshold. In the example β = 1-q-b with j < (k+4)^2/(k^2(1+u)), p >= 1/j + k(k+3)/(k+4) suffices.
- **Exact transports (p.4).** Some group action sends the T24 roles (59,53,55,67,54) to (1,59,6,27,28). Another sends (1,59,6,27,28) to (45,5,7,0,35), which are generators of paper T_65 (supplement "T47").
- **Corollary 3 (p.4-5).** Let d > 0, j > 0, q, b >= 0 with (q+b-d)^2 <= qb. Then p f1 + j f59 + q f27 + b f6 + d f28 is SAFE if p >= d^2/j + 2d - q - b. The transported version for T47: (a+b-d)^2 <= ab, c >= d^2/e + 2d - a - b, e > 0.

## B. Dependencies

| dependency | kind | status |
|---|---|---|
| v -> p v q, (p,q) in S^3 x S^3, maps onto SO(4) | classical fact | PROOF-READ (standard double cover). EXACT: O O^t = |p|^2|q|^2 I and det O = (|p|^2|q|^2)^2 |
| Group invariance Phi_{g.c}(g.O) = Phi_c(O) (paper (4.7)) | core paper | replayed by the shared library. Used for the transports |
| The T47 coefficient names (c on f45, e on f5, leaves a, b) | other docs (T47 sector proofs) | CITED. The formula is symmetric in the two leaves, so only p->c and j->e matter |
| Completeness of the 8x8 cubic reduction (final_referee_T46_cubic_reduction data) | other-doc data | CITED. Used only for the word "complete" in Theorem 2 |
| Preceding "parallel family" threshold | other doc | UNCHECKED. Used only in a comparison remark |

## C. Steps

| # | claim | kind | status | how |
|---|---|---|---|---|
| 1 | The 4x4 metric on (Q,B,J0,D) is positive definite iff q > 0, ℓ > 0, j > (1+u)t^2 | identity | EXACT | 3x3 determinant = ℓ·(j-(1+u)t^2) |
| 2 | P = λ(Q+B) + μJ0 has P·(Q+B) = -h, P·J0 = t-1, ‖P‖^2 = p0; the Gram of (Q+B, J0) has determinant Δ > 0 | identity | EXACT | sympy. Δ > 0 because a positive definite metric gives a positive definite 2x2 Gram (also u - 4/s = (q-b)^2/(sℓ)) |
| 3 | The 8 vectors (Q,P,J,0,0,D,B,0) have the squared lengths (q,p0,j,0,0,1,b,0) | identity | EXACT | |
| 4 | Certificate: with ψ = T^{-1}u over 2 blocks x 4 signed permutations (32 cubics), \|p\|^2·Σ G_kl ψ_k ψ_l ≡ 2Δ^2 Φ_c (mod \|p\|^2 - \|q\|^2), where G = Δ^2·VV·metric·VV^t and c is the Theorem 1 cost | polynomial identity | EXACT | own O(p,q) and own Φ from Appendix D values. Normal form modulo the monic-in-p0 quadric (unique remainder A + p0·B). 413 coefficient rows, rank 18. The symbolic vector is annihilated by every RREF row (sympy cancel) |
| 5 | G is PSD, so the left side is a sum of squares. On \|p\| = \|q\| = 1 this gives Φ_c = SOS/(2Δ^2) >= 0 on all of SO(4). Surplus p - p0 >= 0 times f1 adds det O[12,34]^2 >= 0 | inequality | PROOF-READ | |
| 6 | There are exactly 10 completion relations on the 36 Gram entries. The cost of an allowed Gram is the budget diagonal (q,p,j,c,e,c+d,b,a). The 8 native rays are independent | linear algebra | EXACT | derived from scratch (not from the supplement's R) |
| 7 | Limits of SAFE costs are SAFE | analytic | PROOF-READ | argument supplied by auditor: {c : Φ_c >= 0 on SO(4)} is an intersection of closed half-spaces |
| 8 | Thm 2: setting the f5/f2/f4 rows to zero, the relations have rank 7 and span exactly the seven displayed equations | linear algebra | EXACT | ranks 7, 7, 7 |
| 9 | Thm 2 necessity: β = QB, t = JD, J0 = J - tD give the metric. The 7 equations give P·(Q+B) = -1-β and P·J0 = t-1, so ‖P‖^2 >= p0 by projection | proof | PROOF-READ | each substitution checked by hand |
| 10 | Thm 2 closure: the template is positive definite (minors 1, 7/4, 3/4, 3/4, 1/2) and satisfies the 7 equations with ‖D‖^2 = 1. Convex mixing gives density. Bounded diagonals give compactness | proof | EXACT + PROOF-READ | |
| 11 | t = 0 threshold, d/dt at 0, the witness (1/8, 4, 4, -1/2, 1/4) giving 399/1543 < 13/50, and the old threshold > 13/50 (= 25839/96100 - (168/961)√2 > 0) | identity / inequality | EXACT | |
| 12 | (A t^2 + B t + C)/(D - E t^2) form, stationary quadratic, decomposition (2), derivative | identity | EXACT | |
| 13 | t_opt: the two critical \|t\| multiply to sj/(s+4) = T_nat^2, so one lies inside the natural interval and one outside. f is unimodal on (-T_nat, T_nat), so clipping to (-T, T) is optimal, and T <= T_nat. At the 0/0 endpoint the limit is h^2/s | analytic | PROOF-READ | product identity EXACT. Unimodality argued from the derivative sign |
| 14 | t* threshold; at β = 1-q-b, s = h = k; 1/j + k(k+3)/(k+4); t* is allowed iff j < (k+4)^2/(k^2(1+u)) | identity | EXACT | |
| 15 | Cor 3 strict case: (S-1)^2 < qb forces q, b > 0 and 2/3 < S < 2. Theorem 1 at t = 0, β = 1-S gives threshold 1/j + 2 - S | identity + inequality | EXACT + PROOF-READ | |
| 16 | Cor 3 boundary (S-1)^2 = qb by closure | analytic | PROOF-READ | argument supplied by auditor: g = qb - (q+b-1)^2 has ∇g ≠ 0 on {g = 0} (its only critical point is q = b = 2/3, where g = 1/3), so every boundary point is a limit of strict points. The threshold is continuous and step 7 applies. The doc's explicit cases S = 2 (g(1-ε,1-ε) = 2ε - 3ε^2) and q = 0 are special cases |
| 17 | Homogeneity in d | identity | PROOF-READ | Φ is linear in c |
| 18 | β-optimality remark and the b(4-3b) discriminant | identity | EXACT | |
| 19 | Transports: group elements exist sending (59,53,55,67,54) to (1,59,6,27,28) and (1,59,6,27,28) to (45,5,7,0,35) | finite check | EXACT | 2 elements each, all 1152 searched. In lib convention: (4213,1423,0/1) and (1423,3241,0), (3241,2314,1). The doc's printed row/column labels do not reproduce in the lib convention, so its label convention differs. Existence is what matters |
| 20 | The T47 version of Cor 3 | citation | CITED | coefficient names come from the T47 docs |

## D. Falsification (NUMERIC)

- Theorem 1 at p = p0: 40 random admissible parameter sets (log-uniform q, b, j; β and t uniform inside the admissible region; seed 20261125). Each was minimized over SO(4) with 12 BFGS restarts. The worst normalized minimum was 1.1e-15, so no violation. The minimum is about 0, so the threshold is attained, as expected for a rank-deficient Gram.
- Pointwise check S_G(p,q) = 2Φ_c(O(p,q)) at 5 random unit points: maximum error 3.6e-15.
- Corollary 3 on the boundary (S-1)^2 = qb, including (1,1), (0,1) and (1,0), at 37 points: worst 3.6e-16. No violation.

## E. Verdict

**confirmed-modulo-citations.** Theorem 1 and Corollary 3 are proved in full by an exact, parameter-uniform polynomial certificate checked against my own Leibniz sums, together with proof-read positivity and closure arguments. Their only external input is the classical surjectivity of S^3 x S^3 onto SO(4). The citations are the T47 coefficient naming (for the transported statement) and the completeness of the cubic reduction (only for the word "complete" in Theorem 2). No gap affects correctness.

Presentation notes (not gaps): u is never defined in the document; it is (q+b-2β)/ℓ, implicit in Theorem 1. The general boundary case of Corollary 3 is asserted "by the proved finite closure" without saying why boundary points are limits of interior points (supplied in step 16). The printed row/column labels of the transports use an unstated convention.
