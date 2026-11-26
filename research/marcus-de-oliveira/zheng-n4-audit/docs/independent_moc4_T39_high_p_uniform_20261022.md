# Audit: independent_moc4_T39_high_p_uniform_20261022 (T_40 source, region p >= d, pj < d^2)

Evidence level: PROOF-READ for the whole logical chain (Sections 1-4); EXACT for the restricted cubic C + hQ, for N - ea >= 0 and the substitution; EXACT pointwise (at sampled rational base points) for identity (1) and for the first two certificates of (2); NUMERIC (60-digit) at the positive root for P*(c+) < 0, the bordered minors, rank 12, the signs of all 91 maximal principal minors and the definiteness used in Section 4; UNCHECKED: the global coefficient certificates (2) (832, 3084, 1821 terms), the octic factorization beta^2 - 4 alpha gamma = F+ F-, and the trace-norm identity with L, M (13576 terms each), because the 1313-term G and the other large polynomials are not displayed and could not be rebuilt within the budget.

PDF: proofs/manin_transfer_20261021/independent_moc4_T39_high_p_uniform_20261022.pdf (3 pages). Auditor label C16-T40-T42.

## Deviations and limits

- "T39" is the supplement's historical label. The eight rays f5, f1, f12, f9, f2, f35, f33, f29 vanish exactly on {0,1,5,9,14,18} = paper T_40 (EXACT), and all eight are generators of N_{T_40}; f35 is the only non-minor one. The two remaining generators of N_{T_40}, f11 and f28, are absent; reducing a general T_40 cost to this canonical source is the job of the T_40 assembly document, not of this one.
- I did not rebuild G as a polynomial in the seven base parameters. Instead, for each sampled rational base point (q, p, j, e, b, a; d = 1) I computed the 14 x 14 Hessian of the paper's cubic (A.2) at w(c) = g + h(c) kappa, h = -C/Q, at twelve rational values of c, exactly, and interpolated. This tests identity (1) and the sign claims pointwise but is not a proof of them on the region.
- The sign of G is fixed only up to +-1 by (1); I normalized gamma < 0 and then tested the second certificate. The third certificate needs F+ and F- separately and was not testable; its consequence (P*(c+) <= 0) was tested instead.

## A. Statements

Theorem 1 (p. 1). Let d > 0 and g = q f5 + p f1 + j f12 + c f9 + e f2 + d f35 + b f33 + a f29 with all coefficients >= 0, p >= d and pj < d^2. Then there is no x in C^14 with grad C(x) = 0 and Re x = g + h kappa for some real h. Since X is contained in {grad C = 0} (paper App. A.1, replayed by the shared library), no x in X has gap vector c(x) = g (any reference permutation in T_40, h = Re F_{pi*}). Consequently a spectrally realizable cost of this form with p >= d has pj >= d^2. Zero coefficients and all spectra are included (the statement is about vectors, not about a spectral chart). This is an intermediate exclusion lemma, stronger than the conjecture; it is not itself a case of (*).

## B. Dependencies

| dependency | kind | status |
|---|---|---|
| X subset of {grad C = 0} (paper A.1) | core | CITED; replayed EXACT by the shared library |
| facet vectors f_j (App. D) | core | EXACT via library |
| "already proved SOC capacity pj >= d^2" (the consequence used by the T_40 proof) | other-doc | not used inside this proof; CITED by the users of this document |
| companion script independent_moc4_T39_high_p_uniform_exact_20261022.py | supplement script | not executed (hard rule); its claims are treated as unverified |

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | C(g + h kappa) = C + hQ with C = -(ABc + Nj), Q as displayed | identity | EXACT | sympy from (A.2) |
| 2 | N - ea has nonnegative coefficients, constant 1; so Nj > 0, C < 0 for j < 1 (j < 1 follows from p >= 1, pj < 1) | identity | EXACT | |
| 3 | singular x forces C(w) = 0, hence Q != 0 and h = -C/Q | algebra | PROOF-READ | real part of C(w + i delta) with H(w) delta = 0 (auditor: C(w) = (1/2) delta^T H(w) delta = 0) |
| 4 | (1): det H(w(c)) = -4 G(c)^2 / Q(c)^2 with G quadratic in c | identity | EXACT at sampled points only; global identity UNCHECKED | see D |
| 5 | substitution p = 1 + P, j = J/((1+P)(1+J)) covers p >= 1, 0 <= pj < 1 | algebra | EXACT | pj = J/(1+J) |
| 6 | (2): gamma < 0, AB beta - Nj alpha < 0, AB(F+ + F-) + 2 N alpha < 0 on the region (coefficient-sign certificates) | inequality | UNCHECKED globally; first two EXACT at sampled points | see D |
| 7 | consequences: alpha > 0 at any nonnegative root, unique simple positive root c+, F+ F- > 0, F+ + F- < 0, so (3) F+, F- < 0 | logic | PROOF-READ | |
| 8 | pseudo-division alpha^3 P* = r1 c + r0 mod G; r1 = LM, 2T = F+ M^2 + F- L^2 | identity | UNCHECKED | polynomials not displayed |
| 9 | U+ + U- = T/alpha, U+ U- = Z^2/(16 alpha^2) | algebra | PROOF-READ (auditor re-derived both from Vieta) | |
| 10 | P*(c+) <= 0, strictly off Z = 0 | consequence | PROOF-READ given 8; NUMERIC at sampled roots | |
| 11 | bordered minors on I u {0}, I u {1} vanish on G = 0 | identity | UNCHECKED globally; NUMERIC at sampled roots | |
| 12 | rank 12 and all maximal principal minors <= 0 where Z != 0 | linear algebra | PROOF-READ (auditor: H = W^T A'^{-1} W with W = [B, A'], so every 12-minor is det(W_J)^2/det A' <= 0) | |
| 13 | extension to Z = 0 and to walls/zero coefficients by limits inside the open domain alpha > 0, Q(c+) != 0 | analytic | PROOF-READ | closed conditions; Z not identically zero is confirmed by the strictly negative P*(c+) values found |
| 14 | (4): second derivative -8 G_c^2/Q^2 < 0, so corank exactly 2 | analytic | PROOF-READ | from (1) at a simple root |
| 15 | Section 4: H(w) delta = 0, grad F(delta) = grad F(w), delta != 0; Q_v(delta) = 0; Schur expansion d^2/dc^2 det = 8 det A' det(K^T Q_v K); hence Q_v definite on ker H, contradiction | analytic | PROOF-READ (auditor checked that the t^2 term of det(H0 + tH1 + O(t^2)) is det A' det(K^T H1 K) for corank 2 and graph basis K) + NUMERIC | |

## D. Falsification

Pointwise exact tests (script `replay/independent_moc4_T39_high_p_uniform_20261022__check.py`, seed 7, base points with q, e, a, b, P, J random rationals including zeros):

- 60 of 60 base points: -det H(w(c)) Q(c)^2/4 is a polynomial of degree 4 in c (interpolated on 9 values of c, confirmed on 3 more) and is a perfect square G(c)^2 with G quadratic. EXACT at each point.
- 60 of 60: with the sign of G chosen so that gamma < 0 (gamma was never 0), AB beta - Nj alpha < 0. EXACT at each point. This is consistent with the first two certificates of (2) but does not prove them on the region.
- 60 of 60: P*(c) has degree <= 4 in c and does not depend on h (EXACT: H(kappa) only has entries (0,6), (1,5) and their transposes).
- 8 base points had alpha > 0 (a positive root c+). At each, with 60-digit arithmetic: P*(c+) < 0 (largest value -3.6e7), both bordered minors vanish (|.| < 2e-42), H has exactly two singular values < 1e-40, all 91 principal 12 x 12 minors are strictly negative (largest -223), and the binary form K^T Q_v K on the kernel has positive determinant (smallest 4.5e-16; working precision about 1e-40). NUMERIC.

No violation of any intermediate claim was found. Limitation: only 8 positive-root samples; the sampling was not aimed at the boundary p = 1 or pj -> 1.

## E. Verdict

could-not-check (no gap found). The logical chain is sound as written; every inference was checked line by line, and every large identity that I could test pointwise held at every sampled point. What keeps the verdict below "confirmed" is that the global certificates of (2), the octic factorization, the trace-norm identity and the bordered-minor divisibility rest on undisplayed polynomials of 300 to 13576 terms that I did not rebuild; they are supported only by exact pointwise tests and 60-digit evaluations. A future replay should rebuild G (for example by multivariate interpolation of -det H Q^2/4) and check the three coefficient tables.
