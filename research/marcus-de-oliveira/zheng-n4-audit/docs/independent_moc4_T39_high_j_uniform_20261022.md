# Audit of independent_moc4_T39_high_j_uniform_20261022 (high-j chamber of paper T_40)

Evidence level: EXACT for C, Q, the N bound and the chart inverse. PROOF-READ for the sign logic of Section 1, the trace-norm argument, and the rank-12 and Schur argument of Sections 2-3 (each re-derived by the auditor). NUMERIC for the four strict inequalities at 15 points and for the Hessian structure at the positive root at 7 points. UNCHECKED for the eight coefficient tables, the pseudo-remainder identity, the L16/M16 factorization and the bordered-minor divisibility as polynomial identities.

Auditor label C15-T38-T40. Document: proofs/manin_transfer_20261021/independent_moc4_T39_high_j_uniform_20261022.pdf, 3 pages. Replay: `replay/independent_moc4_T39_high_j_uniform_20261022__check.py` (0 failed checks).

## Deviations and limits

- I did not rebuild G (1313 terms), P6, F±, P∗ (5188 terms), L16 or M16. None of the 8 strict tables was replayed.
- At each rational test point I computed G as a function of c exactly. I interpolated -(Q^2/4) det Hess along c, certified the degree with 4 extra points, and took the square-free part. G is defined only up to sign by G^2, so I fixed the sign by the claimed gamma < 0 and then tested whether ABβ - N_j α < 0.
- At the positive root c+ (to 40 digits, computed with 50-digit mpmath) I computed the following: the singular values of Hess C(w(c+)); P∗ = the principal minor that deletes indices 0 and 1; the two bordered 13-minors; and the determinant of the binary form K^T Hess(w_c) K on the numerical kernel.
- First observation: in 20 random chamber points α was <= 0 every time, so G has no root at all and the Section 2-3 machinery is not needed there. A wider search (23668 random points, many near j = 1 with large q, b, a, e) found 558 points with α > 0; 7 of them are used in the replay.

## A. Statements

Theorem 1 (p. 1). With the canonical source g = q f5 + p f1 + j f12 + c f9 + e f2 + d f35 + b f33 + a f29 (all coefficients >= 0), every spectrally realizable cost with j >= d satisfies Phi >= 0 on O(4). More precisely, with d = 1, no point of the chamber (1) {j >= 1, pj < 1, (j-1)(a+e) < 1} is Re x - (height) kappa for x in X. j = 1 and all zero resources are included. The scope is the stated chamber of T_40, used as a branch by independent_moc4_T39_uniform_20261022.

## B. Dependencies

| kind | ref | used for | status |
|---|---|---|---|
| core-doc | SOC and old-H atoms (capacities pj >= d^2, j >= d + d^2/(a+e)) | reduction to (1) | CITED |
| common-lemma | X ⊂ Sing C (paper App. A) | radial necessity | EXACT in lib |

## C. Steps

| # | claim | kind | status |
|---|---|---|---|
| 1 | C = -(ABc + Nj), -Q formula | identity | EXACT |
| 2 | N - (1+S+ea) has nonnegative coefficients; Nj > 1 + 3S/4 (S > 0), Nj >= 1; C < 0 and Q < 0 on (1) | identity + inequality | EXACT + PROOF-READ |
| 3 | det Hess(w(c)) = -4 G^2 / Q^2; (3) dh-identity | identity | UNCHECKED as polynomials; square property EXACT at 15 points (via interpolation) |
| 4 | chart (5) covers (1) for j > 1 | analytic | EXACT (sympy) + PROOF-READ |
| 5 | the four strict inequalities (4), from 3648 / 14147 / 14123 / 14123 negative coefficients, plus the j = 1 wall tables | certificate | UNCHECKED; NUMERIC: gamma < 0 and ABβ - Njα < 0 at 15 points |
| 6 | the sign logic: α <= 0 means no root on c >= 0; α > 0 means a unique simple positive root; F± < 0 | logic | PROOF-READ |
| 7 | pseudo-remainder and trace-norm identity: U+ + U- = T/α <= 0, U+ U- = Z^2/(16α^2), so P∗(c+) <= 0 | identity + logic | algebra from the displayed relations PROOF-READ (auditor derived U+U- = Z^2/(16α^2) from r1 = LM and 2T = F+M^2 + F-L^2); the 13576-term identities UNCHECKED |
| 8 | the bordered 13-minors are divisible by G; on Z != 0 rank is 12 and all principal 12-minors are <= 0 (by congruence) | identity + logic | divisibility UNCHECKED; logic PROOF-READ |
| 9 | extension to the closed chamber by density of Z != 0 and continuity of the simple root | analytic | PROOF-READ (closed conditions; α > 0 is open near a root point) |
| 10 | d^2/dc^2 det = -8 G_c^2/Q^2 < 0 excludes corank >= 3, so the rank is exactly 12 and some principal 12-minor is < 0 | analytic | PROOF-READ |
| 11 | real-lift contradiction: H(w)δ = 0, grad F(δ) = grad F(w), Q != 0 gives δ != 0; Q_v(δ) = 0; the Schur expansion gives d^2 det = 8 det A' det(K^T Q_v K); definiteness gives a contradiction | analytic | PROOF-READ (factor 8 = 2·4 checked; curvature irrelevant at corank 2) |

## D. Falsification

- NUMERIC (seed 40, 8 random points with j = 1 included, plus 7 α > 0 points from an auditor search): G is quadratic in c with gamma < 0 and ABβ - Njα < 0 at all 15 points.
- At the 7 α > 0 points there is exactly one positive root c+. Hess C(w(c+)) has two singular values below 1e-40 (the next is about 1). P∗ < 0. Both bordered 13-minors vanish to relative precision 1e-25. det(K^T Hess(w_c) K) > 0 and Q < 0. This matches every claim of Sections 2-3, and no realizable point was found.

## E. Verdict

could-not-check. No gap found. The analytic argument of Sections 1-3 is sound as written, and I re-derived each step. The numerical tests at α > 0 points match every intermediate claim. The polynomial certificates (strict tables, pseudo-remainder, L16/M16, divisibility) were not replayed.
