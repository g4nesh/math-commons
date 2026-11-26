Evidence level: EXACT for every polynomial identity, factorization and coefficient-sign certificate of the document (rebuilt by own python-flint code from the paper's cubic (A.2) and the Appendix D facet vectors, on the whole chamber including its walls); PROOF-READ for the analytic steps that use them (first auditor's reading re-examined against the replayed polynomials); CITED: X contained in {grad C = 0} (paper App. A.1).

# Certificate replay: independent_moc4_T39_high_p_uniform_20261022 (paper T_40, chamber p >= d, pj < d^2)

Auditor label cert-T39. Continues `docs/independent_moc4_T39_high_p_uniform_20261022.md` (first verdict could-not-check because the global certificates were not rebuilt). Replay: `replay2/independent_moc4_T39_high_p_uniform_20261022__certificates.py`, output `..._certificates_output.txt` (0 failures, 113 s). This script also holds the shared core used by the low_pj and high_j replays.

## Deviations and limits

- None in scope. Every polynomial was rebuilt from definitions: w = q f5 + p f1 + j f12 + c f9 + e f2 + f35 + b f33 + a f29 + h e13 in the paper's coordinates m_0..m_13 (d = 1, kappa = e_13), the 14 x 14 Hessian of the cubic (A.2) at w, its determinant (fraction-free elimination, 26218 terms), and from there G, P6, F+-, P*, r0, r1, L, M, Z and the bordered minors. The supplement's programs were read (not run) only to learn which polynomial each certificate is; its JSON was loaded as data for comparison.
- Normalization of G: G is fixed by (1) only up to sign; I take constant term -1 (it is +-1). The certificates hold for this one fixed sign.
- Section 4 (the Schur-expansion contradiction) and the limit argument of Section 3 contain no certificate; I re-examined them (below) but did not add a machine check.

## What each certificate asserts and how it was rebuilt

| # | claim in the PDF | rebuilt as | result |
|---|---|---|---|
| 1 | C(g + h kappa) = C + hQ with the displayed C, Q | cubic (A.2) composed with w | EXACT, equal |
| 2 | N - ea has nonnegative coefficients, constant 1 | direct | EXACT |
| 3 | (1) det H_F(w(c)) = -4 G^2/Q^2, G = alpha c^2 + beta c + gamma, 1313 terms | det H has height degree 3 (H(kappa) has entries only at (0,6),(1,5)); Q^3 det H(-C/Q) = sum_i D_i (-C)^i Q^(3-i) is divisible by -4Q and the quotient is a perfect square (flint sqrt) | EXACT; G has 1313 terms, degree 2 in c, total degree 9 |
| 4 | D_h det H at h = -C/Q equals 4 G P6/Q, P6 464 terms | Q^2 D'(-C/Q) = 4 G P6 Q, exact division | EXACT, 464 terms (not used in the logic) |
| 5 | beta^2 - 4 alpha gamma = F+ F-, two 306-term octics | flint factorization; unit 1, two simple factors, product rechecked | EXACT |
| 6 | (2a) gamma < 0 | substitute p = 1+P, j = J/((1+P)(1+J)), multiply by ((1+P)(1+J))^deg_j | EXACT: 832 terms, all < 0, constant -2 |
| 7 | (2b) AB beta - N_j alpha < 0 | same | EXACT: 3084 terms, all < 0, constant -2 |
| 8 | (2c) AB(F+ + F-) + 2N alpha < 0 | same | EXACT: 1821 terms, all < 0, constant -4 |
| 9 | P* (delete indices 0,1) height free, degree 4 in c, 5188 terms | 12 x 12 determinant | EXACT |
| 10 | alpha^3 P* = r1 c + r0 mod G | pseudo-division with the quotient tracked; alpha^3 P* = Quot G + r1 c + r0 checked as a full polynomial identity | EXACT; r0 578728, r1 527140 terms |
| 11 | r1 = L M, each 13576 terms | flint factorization, unit 1, product rechecked | EXACT |
| 12 | 2T = F+ M^2 + F- L^2, T = 2 alpha r0 - beta r1 | full identity (labelling of F+- fixed by it) | EXACT |
| 13 | Z = F+ M^2 - F- L^2 nonzero; U+ + U- = T/alpha, U+ U- = Z^2/(16 alpha^2) | Z has 1619914 terms; 4(T^2 - disc r1^2) = Z^2 and T^2 - disc r1^2 = 4 alpha(alpha r0^2 - beta r0 r1 + gamma r1^2) as full identities (the second is the Vieta product numerator) | EXACT |
| 14 | bordered minors on I u {0}, I u {1} (I = {2..13}) have height degree 2; Q^2 minor(-C/Q) divisible by G, quotients 2021, 1091 terms | 13 x 13 determinants (raw 20639, 16186 terms), projection, division, product rechecked | EXACT |

Data comparison (supplement JSON as data): my three high-p certificate tables are coefficient-for-coefficient equal to the stored ones; my octics and my L, M equal the stored pairs up to order and sign.

Region cover. The chart is a bijection from P, J >= 0 onto {p >= 1, 0 <= pj < 1, j >= 0}: inverse P = p - 1, J = pj/(1 - pj); then 1 + P = p and J/((1+P)(1+J)) = pj/p = j. The multiplier is positive. Each numerator has all coefficients strictly negative and a nonzero (negative) constant, so it is < 0 at every point of the closed orthant; hence (2a)-(2c) hold strictly on the entire chamber including p = 1, j = 0 and every zero resource. EXACT.

## Re-examination of the dependent PROOF-READ steps

- Root structure (p. 2): N_j > 0 (item 2 and j < 1), so alpha <= 0 with (2b) gives beta < 0 and G < 0 on c >= 0. Hence any root c >= 0 has alpha > 0; with gamma < 0 there is one positive and one negative simple root and G_c(c+) > 0. disc = beta^2 - 4 alpha gamma > 0 holds already from alpha > 0 > gamma, so F+ F- > 0; (2c) with AB > 0, N alpha > 0 gives F+ + F- < 0, hence (3) F+, F- < 0. Correct.
- Section 2: at c+, alpha^3 P*(c+) = U+ (item 10); U+- real with U+ + U- = T/alpha = (F+ M^2 + F- L^2)/(2 alpha) <= 0 and U+ U- = Z^2/(16 alpha^2) >= 0 (items 12-13), so U+ <= 0, < 0 if Z != 0. Correct.
- Rank 12: with P* < 0 the block on I is invertible; items 14 and 3 make the 2 x 2 Schur complement zero. Correct. Z depends only on the base parameters (q, p, j, e, b, a) and is a nonzero polynomial, so Z != 0 is dense in the open set {alpha > 0, Q(c+) != 0}; the extension to Z = 0, walls and zero resources is by closed conditions along the continuous simple root. Correct.
- (4) follows from (1) at a simple root with Q != 0. Section 4: Re/Im of grad F(w + i delta) = 0 give grad F(w) = grad F(delta) and H(w) delta = 0; F is linear in m13 with coefficient Q, so Q != 0 forces delta != 0; then F(w) = 0 gives h = -C/Q and det H(w) = 0 gives G(c) = 0. The remaining Schur argument is as read by the first auditor. No step depends on anything not now replayed.

## Verdict

confirmed-modulo-citations. Every certificate, factorization and identity the proof uses is EXACT on the whole chamber including its boundary; the analytic steps are PROOF-READ. Citation: X subset of {grad C = 0} (paper App. A.1; also replayed by the shared library). The consequence "a realizable source with p >= d has pj >= d^2" is what the assembly uses.
