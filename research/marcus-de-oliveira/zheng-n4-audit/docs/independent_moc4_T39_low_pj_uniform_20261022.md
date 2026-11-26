# Audit: independent_moc4_T39_low_pj_uniform_20261022 (T_40 source, region 0 <= p, j < d)

Evidence level: PROOF-READ for the logical chain (Sections 1-4); EXACT for the restricted cubic, C < 0 and -Q > 0 on the region, the substitution; EXACT pointwise (sampled rational base points) for identity (1) and the two certificates (3); NUMERIC (60-digit) at the positive root for (8), the bordered minors, rank 12, all maximal principal minors and the Section 4 definiteness; UNCHECKED: the global coefficient certificates (3) and (5) (696, 1985, 1457, 1457 terms), the octic factorization (4), the pseudo-division and trace-norm identities (6)-(7), the bordered-minor divisibility, because G and the other large polynomials are not displayed and were not rebuilt.

PDF: proofs/manin_transfer_20261021/independent_moc4_T39_low_pj_uniform_20261022.pdf (3 pages). Auditor label C16-T40-T42.

## Deviations and limits

- Same source and the same method as the companion document independent_moc4_T39_high_p_uniform_20261022 (see that report for the identification "T39" = paper T_40 = {0,1,5,9,14,18}; EXACT). The two documents share one replay design; this document's replay is `replay/independent_moc4_T39_low_pj_uniform_20261022__check.py` (region taken from the file name, seed 8).
- G is determined by (1) only up to sign; I normalized gamma < 0. Certificate (5) needs F+ and F- separately and was not testable; its consequence P*(c+) <= 0 was tested.

## A. Statements

Theorem 1 (p. 1). Let d > 0 and g = q f5 + p f1 + j f12 + c f9 + e f2 + d f35 + b f33 + a f29 with all coefficients >= 0, 0 <= p < d, 0 <= j < d. Then no x with grad C(x) = 0 has Re x = g + h kappa for a real h. Hence (via X subset of {grad C = 0}) every spectrally realizable cost of this form with d > 0 has p >= d or j >= d. All zero coefficients included; no spectral chart is used. This is an exclusion lemma for a branch of T_40, not a case of (*) itself. The document says the high-p and high-j regions are separate tasks; together with the high-p document, the remaining T_40 region is {p < d, j >= d} plus {p >= d, pj >= d^2} (the latter covered by a cited "SOC capacity").

## B. Dependencies

| dependency | kind | status |
|---|---|---|
| X subset of {grad C = 0} (paper A.1) | core | CITED; replayed EXACT by the shared library |
| facet vectors f_j (App. D) | core | EXACT via library |
| companion script independent_moc4_T39_low_pj_uniform_exact_20261022.py | supplement script | not executed; claims unverified |

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | C = -(ABc + Nj) < 0 and -Q > 0 on the region, so Q never vanishes | identity + sign | EXACT | sympy; -Q - (1-p)c - (1-j)(e+a) - pj - 1 has positive coefficients |
| 2 | singular x forces C(w) = 0, so h = -C/Q | algebra | PROOF-READ | as in the high-p report |
| 3 | (1) det H = -4G^2/Q^2; (2) the h-derivative formula | identity | (1) EXACT at sampled points, global UNCHECKED; (2) not used and UNCHECKED |
| 4 | (3): gamma < 0, AB beta - Nj alpha < 0 under p = P/(1+P), j = J/(1+J) | inequality | UNCHECKED globally; EXACT at sampled points | |
| 5 | alpha > 0 at any nonnegative root; unique simple positive root; G_c(c+) > 0 | logic | PROOF-READ | |
| 6 | (4) disc = F+ F-; (5) AB F+- + N alpha < 0, so F+- < 0 where alpha > 0 | identity + inequality | UNCHECKED | |
| 7 | (6)-(7) pseudo-division and trace-norm identity; U+ + U- = T/alpha <= 0, U+ U- = Z^2/(16 alpha^2) | identity + algebra | identities UNCHECKED; the Vieta algebra re-derived by the auditor (PROOF-READ) | |
| 8 | (8) P*(c+) <= 0, < 0 off Z = 0 | consequence | PROOF-READ given 6-7; NUMERIC at sampled roots | |
| 9 | bordered minors divisible by G; rank 12 on Z != 0; all 12-minors <= 0 | identity + linear algebra | divisibility UNCHECKED globally (NUMERIC at roots); the linear algebra PROOF-READ (H = W^T A'^{-1} W, each 12-minor = det(W_J)^2/det A') | |
| 10 | extension by limits to Z = 0 and to zero coefficients (open domain alpha > 0 inside p, j < 1) | analytic | PROOF-READ | |
| 11 | (9) second derivative < 0, corank exactly 2, some 12-minor < 0 | analytic | PROOF-READ | |
| 12 | Section 4 contradiction (Q_v definite on ker H, Q_v(delta) = 0, delta != 0) | analytic | PROOF-READ + NUMERIC | |

## D. Falsification

Replay output `replay/independent_moc4_T39_low_pj_uniform_20261022__check_output.txt` (seed 8, 60 base points with random rationals including zeros):

- 60/60: -det H Q^2/4 is a degree-4 perfect square G^2 in c (EXACT per point); gamma < 0 and AB beta - Nj alpha < 0 for that sign (EXACT per point); P* has degree <= 4 in c.
- 14 base points had a positive root. At each (60 digits): P*(c+) < 0 (largest -1.7e6), both bordered minors vanish (< 2e-45), rank exactly 12, all 91 maximal principal minors < 0 (largest -21.2), det(K^T Q_v K) > 0 (smallest 3.3e-10). NUMERIC.

No violation found.

## E. Verdict

could-not-check (no gap found). The written argument is logically complete, and every large identity and sign claim held at every sampled point, but the global certificates rest on undisplayed polynomials (G with 1313 terms, two 306-term octics, L and M with 13576 terms) that I did not rebuild. A follow-up replay needs only G: from it F+-, r0, r1, L, M can be recomputed and the four coefficient tables checked.
