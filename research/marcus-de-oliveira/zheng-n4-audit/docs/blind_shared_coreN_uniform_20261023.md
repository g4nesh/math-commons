# Audit: blind_shared_coreN_uniform_20261023 (seven-contact set N, Table 5 no. 24)

Evidence level: every algebraic identity and coefficient claim EXACT (own flint/sympy replay from the paper's definitions); the analytic root-graph, rank and inertia argument PROOF-READ; the final contradiction CITED to the paper's Proposition C.6; the unitary sentence CITED to Theorem 1.2. Numerical tests are NUMERIC only.

PDF: `proofs/manin_transfer_20261021/blind_shared_coreN_uniform_20261023.pdf`, 4 pages. Replay: `replay/blind_shared_coreN_uniform_20261023__check.py` (output in `__check_output.txt`, about 30 s).

## Deviations and limits

- The document's certificate files (`..._exact_20261023.py/json`) were not executed or loaded. Every polynomial (C, Q, the 14 x 14 determinant, G, R, the three 13 x 13 borders) was recomputed from scratch from the paper's cubic (A.2), the matrix P and the vectors f_j of Appendix D, through the shared library `lib/moc4lib.py`.
- My Bareiss uses full pivoting and does 650 exact divisions for the 14 x 14 determinant and 506 for each 13 x 13 border. The document's counts (819, 1950) belong to its own elimination order. The counts are not mathematical claims, so I did not try to reproduce them.
- Proposition C.6 (Appendix C) is used as a black box. This batch did not audit it.

## A. Statements

Contact set N = {0,3,7,11,12,16,20} = paper Table 5 set no. 24 (EXACT, `AUX7[24]`).

**Theorem 1 (p. 1).** Let x be in X and pi* be such that c = c(x) satisfies Pc >= 0 and (Pc)_sigma = 0 for every sigma in N. Further contacts, collisions, a_i + b_j = 0 and zero coefficients are allowed. Then c is a nonnegative combination of the six minor vectors f2 = n_{13,34}, f8 = n_{12,23}, f12 = n_{14,24}, f15 = n_{12,14}, f29 = n_{13,12} and f39 = n_{14,13}. Equivalently, the coefficient d of the only non-minor generator f30 = h_{4321,-} is 0. Hence Phi_c(O) >= 0 for every O in SO(4), and by Theorem 1.2 the same holds on U(4).

Coverage of spectra: the argument uses only that x lies in X, that X lies in Sing C, and Prop. C.6, which is stated for every x in X. All spectra are therefore covered, degenerate ones included. `covers_all_spectra` = true.

## B. Dependencies

| dependency | kind | used for | match |
|---|---|---|---|
| Paper Prop. C.6 (Hessian inertia of real parts of matching points): for every x in X, every 12 x 12 principal minor of H(Re x) is >= 0 | core-doc (paper App. C) | final contradiction | Checked. The document's "full real lift" w + z kappa equals Re x. The reason is that c(x) = Re(x - F_{pi*} kappa), so Re x = c + (Re F_{pi*}) kappa. The document's "original cubic" is the paper's C in the paper's coordinates: my replay with (A.2) reproduces every number the document states (38/19 terms, 498-term G, det A_I = -49152 at the sample point). A real rescaling of the cubic would not change the sign of an even-order minor. CITED, not audited here. |
| X contained in Sing C (paper A.1) | paper | singular Hessian at realizable lifts | EXACT in lib (grad C vanishes identically on matching vectors) |
| (4.2): Phi_{n_{I,J}}(O) = det O[I,J]^2 | paper | d = 0 gives Phi >= 0 | EXACT in lib |
| Theorem 1.2 / Corollary 2.2 | paper | the "unitary response" sentence only | CITED |

The document cites no other contact case, branch or Appendix B result.

## C. Steps

| # | claim (location) | kind | status | how checked |
|---|---|---|---|---|
| 1 | rank P = 14; rows at N have rank 7; generators(N) = {f2,f8,f12,f15,f29,f39,f30}, independent, so they span the kernel (p. 1) | identity | EXACT | lib rank and `generators` |
| 2 | rows [13,14,6,5,10,1,23] of these generators = diag(1,1,1,1,1,1,2), so the decomposition (1) is unique, all coefficients are >= 0, and d = g_23/2 (p. 1) | identity | EXACT | integer table |
| 3 | f2, ..., f39 are minor vectors, so d = 0 gives Phi >= 0 | identity | EXACT | lib (4.2) |
| 4 | F(w + z kappa) = C + zQ (2); C has 38 positive terms with [d^3]C = 1; Q has 19 positive terms with [d^2]Q = 1; so C, Q > 0 for x >= 0, d > 0 (pp. 1-2) | identity | EXACT | flint expansion of (A.2) along w. There is no z^2 or z^3 term because kappa = e_13 and m_13 enters C only linearly. |
| 5 | A realizable lift has F(Re x) = 0 and singular H(Re x). Hence z = -C/Q and det A = 0 (p. 2) | analytic | PROOF-READ | Re x = u + ubar with u = x/2 in X ⊂ Sing C (C has real coefficients). With T the symmetric trilinear form of C: C(u + ubar) = 0, and H(u + ubar)(u - ubar) = 0. If u = ubar, Euler gives H(2u)(2u) = 0 with 2u = Re x ≠ 0 (d > 0). |
| 6 | Identity (4): Q^4 det A = -4 Q^2 G^2 with G homogeneous of degree 9, 498 terms, [d^9]G = 1, individual degrees (2,2,2,2,2,2,9) (p. 2) | identity | EXACT | det H(w + z kappa) in 8 variables (z-degree 3). Then sum_k D_k (-C)^k Q^{4-k} is divided exactly by -4Q^2, and the quotient is an exact square (flint `sqrt`). |
| 7 | Coefficient signs: d^k with k >= 4 positive, k <= 2 negative; R = dG_d - 3G has 427 positive terms including 6d^9 (6), so d(G/d^3)/dd = R/d^4 > 0 for x >= 0, d > 0 (7) (p. 2) | inequality | EXACT | flint. The d^3 coefficient has 70 positive terms and 1 negative term, and it cancels in R. |
| 8 | A positive root is unique and simple, with G_d = R/d > 0 there | analytic | PROOF-READ | from 7 |
| 9 | G(x,0) = -16 abcefg C(x,0), where C(x,0) has 16 positive cubic terms (8). G -> +inf as d -> inf. So each x in (0,inf)^6 has exactly one root rho(x), rho is real analytic, and its graph is connected (pp. 2-3) | identity + analytic | EXACT + PROOF-READ | flint, IFT |
| 10 | Every root (x0,d0) with x0 >= 0, d0 > 0 is a limit of graph points (x0 + eps 1, rho(x0 + eps 1)) (p. 3) | analytic | PROOF-READ | IFT at (x0,d0) uses G_d > 0. The root is unique by 8. |
| 11 | For each of the three borders ({1}∪I,{1}∪I), ({0}∪I,{0}∪I), ({1}∪I,{0}∪I) with I = {2..13}: the z-degree is <= 2, and Q^2 det(...) at z = -C/Q is divisible by G (9) (p. 3) | identity | EXACT | 13 x 13 Bareiss on the symbolic Hessian; exact division by G |
| 12 | Where det A_I ≠ 0, the 2 x 2 Schur complement vanishes, so rank A = 12 (p. 3) | linear algebra | PROOF-READ | The Schur entries are the three border minors divided by det A_I. |
| 13 | Sample (10): at a = ... = d = 1, C = 38, Q = 19, G = 0, z = -2, det A_I = -49152 (p. 3) | identity | EXACT | exact Fractions; also rank A = 12 there |
| 14 | det A_I(x, rho(x)) is real analytic and not identically 0, so it is nonzero on a dense set. Continuity then gives rank A <= 12 on the graph and at every boundary root (p. 3) | analytic | PROOF-READ | |
| 15 | (11): d^2/dd^2 det A at a root = -8 G_d^2/Q^2 < 0. Rank <= 11 would force vanishing order >= 3, so rank >= 12 (p. 3) | identity + linear algebra | EXACT (calculus, sympy) + PROOF-READ | If rank M(0) = r, then det M(t) = O(t^{14-r}). A is analytic in d because Q >= d^2 > 0. |
| 16 | Rank 12 is constant on the connected graph, so the inertia is constant. A nonzero principal 12-minor has the sign of the product of the nonzero eigenvalues (Cauchy-Binet). By 13 that product is < 0 everywhere, so some principal 12-minor is < 0, also at boundary roots (p. 4) | analytic | PROOF-READ | The limit argument at boundary roots is valid: rank stays 12 in the limit, so no nonzero eigenvalue can reach 0. |
| 17 | This contradicts Prop. C.6, so d = 0 (p. 4) | citation | CITED | Prop. C.6 |
| 18 | Unitary extension | citation | CITED | Theorem 1.2 |

The logic is complete: realizable with d > 0 implies G = 0 (5, 6); then x >= 0 and d > 0 is a root, on the graph or at a boundary root (9, 10); then a negative principal 12-minor (11-16); this contradicts C.6. I found no unjustified step.

## D. Falsification (NUMERIC)

- 40 random x >= 0 (exponential, seed 20261125; every fourth with one zero coordinate). Univariate roots of G(x,·): at most one positive root in every case. At every root, the floating-point inertia of A was (7,5,2) or (5,7,2), never (6,6,2). This is consistent with steps 8 and 16.
- Realizable points with contacts containing N, found by lib least squares (seeds 0..19): all 20 seeds returned a point, and the predicted consequence d = g_23/2 = 0 held at every one (max |g_23| / mean off-N gap = 2.3e-10, at the solver tolerance). So the stratum appears to force the extra contact 23, as the theorem says.
- No search for violations of the main inequality was done. The theorem reduces it to minor vectors, whose Leibniz sums are exact squares.

## E. Verdict

**confirmed-modulo-citations.** Citations: paper Prop. C.6 (the only mathematical input beyond the replayed identities) and Theorem 1.2 (only for the unitary sentence). No gaps found. Minor presentation remarks: the division counts (819/1950) are tied to the document's own elimination order, and the d^3 coefficient of G does have mixed sign, 70 positive terms and 1 negative.
