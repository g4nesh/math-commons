Evidence level: EXACT for every polynomial identity, factorization and coefficient-sign certificate of the document (own python-flint rebuild from the paper's cubic (A.2) and facet vectors, whole chamber 0 <= p, j < d including walls); PROOF-READ for the analytic steps that use them; CITED: X contained in {grad C = 0} (paper App. A.1).

# Certificate replay: independent_moc4_T39_low_pj_uniform_20261022 (paper T_40, chamber 0 <= p < d, 0 <= j < d)

Auditor label cert-T39. Continues `docs/independent_moc4_T39_low_pj_uniform_20261022.md` (first verdict could-not-check). Replay: `replay2/independent_moc4_T39_low_pj_uniform_20261022__certificates.py`, output `..._certificates_output.txt` (0 failures, 97 s). It imports the shared core from the high_p replay in the same directory, which it re-runs in full.

## Deviations and limits

- None in scope. Same construction as the high_p certificate report (same source, same G, F+-, P*, r0, r1, L, M, bordered minors; G normalized to constant term -1). Supplement programs read, not run; JSON used only for comparison.
- Section 4 and the limit argument contain no certificate; re-examined, not machine-checked.

## What each certificate asserts and how it was rebuilt

| # | claim in the PDF | rebuilt as | result |
|---|---|---|---|
| 1 | C = -(ABc + N_j) < 0; -Q equals the displayed positive expression | restriction of (A.2); N - ea >= 0 coefficientwise with constant 1; -Q - (1-p)c - (1-j)(e+a) - pj - 1 has nonnegative coefficients; displayed -Q equal | EXACT; so Q < 0 on the chamber |
| 2 | (1) det H = -4G^2/Q^2, G 1313 terms; (2) h-derivative = 4 G P6/Q, P6 464 terms | shared core | EXACT |
| 3 | (3) gamma < 0 | substitute p = P/(1+P), j = J/(1+J), multiply by (1+P)^deg_p (1+J)^deg_j | EXACT: 696 terms, all < 0, constant -1 |
| 4 | (3) AB beta - N_j alpha < 0 | same | EXACT: 1985 terms, all < 0, constant -1 |
| 5 | (4) beta^2 - 4 alpha gamma = F+ F- (306 terms each) | flint factorization | EXACT |
| 6 | (5) AB F+ + N alpha < 0 and AB F- + N alpha < 0 | same chart | EXACT: 1457 and 1457 terms, all < 0, constants -1 |
| 7 | P* 5188 terms, degree 4 in c, height free; (6) pseudo-division; (7) r1 = LM (13576 terms each), 2T = F+ M^2 + F- L^2; Z nonzero; U+ U- = Z^2/(16 alpha^2) | shared core, all as full polynomial identities | EXACT |
| 8 | bordered minors on I u {0}, I u {1}: raw 20639, 16186 terms, height degree 2, projected numerators divisible by G with quotients 2021, 1091 | shared core | EXACT, all counts as claimed |

Data comparison: my gamma and AB beta - N_j alpha tables equal the stored ones coefficient for coefficient; my two F+- certificate tables, the octics and L, M equal the stored pairs up to order and sign.

Region cover. p = P/(1+P), j = J/(1+J) maps P, J >= 0 bijectively onto 0 <= p, j < 1 (inverse P = p/(1-p), J = j/(1-j)); the multiplier is positive. All coefficients are strictly negative with a nonzero constant, so the four inequalities are strict on the entire closed chamber, every zero resource included. EXACT.

## Re-examination of the dependent PROOF-READ steps

- With N_j > 0, (3) gives: alpha <= 0 forces beta < 0 and no root c >= 0; otherwise alpha > 0, one simple positive root c+, G_c(c+) > 0. (5) with alpha > 0 gives F+ < 0 and F- < 0 individually. Correct.
- (8): alpha^3 P*(c+) = U+, U+ + U- = T/alpha <= 0, U+ U- = Z^2/(16 alpha^2) >= 0, so P*(c+) <= 0 and < 0 off Z = 0 (Z is a nonzero polynomial in the base parameters only). Correct.
- Rank 12 on Z != 0 from P* < 0, the two vanishing bordered minors and (1); extension to Z = 0 and zero coefficients by closed conditions along the continuous simple root (Q never vanishes here, item 1); (9) from (1); Section 4 contradiction as read by the first auditor (delta != 0 because the m13-coefficient of F is Q != 0). No step depends on anything not now replayed.

## Verdict

confirmed-modulo-citations. All four strict certificates, the discriminant factorization, pseudo-division, trace-norm identity and bordered-minor divisibility are EXACT on the whole chamber; analytic steps PROOF-READ. Citation: X subset of {grad C = 0} (paper App. A.1).
