Evidence level: EXACT for every polynomial identity, factorization and all eight coefficient-sign certificates (own python-flint rebuild from the paper's cubic (A.2) and facet vectors; chart (5) for j > d and the wall j = d, boundaries included); EXACT for the T_40 cone data used by the assembly; PROOF-READ for the analytic steps and the four-region cover; CITED: X contained in {grad C = 0} (paper App. A.1), and the SOC and old-H capacities (only for the "support inequality" wording of Theorem 1, not for the exclusion).

# Certificate replay: independent_moc4_T39_high_j_uniform_20261022 (paper T_40, chamber j >= d, pj < d^2, (j-d)(a+e) < d^2)

Auditor label cert-T39. Continues `docs/independent_moc4_T39_high_j_uniform_20261022.md` (first verdict could-not-check). Replay: `replay2/independent_moc4_T39_high_j_uniform_20261022__certificates.py`, output `..._certificates_output.txt` (0 failures, 98 s). It imports the shared core from the high_p replay in the same directory and re-runs it in full; it also contains the chamber-cover checks for the assembly document.

## Deviations and limits

- None in scope. Shared core as in the high_p certificate report (G normalized to constant term -1). Supplement programs read, not run; JSON used only as data.
- The four-region cover is proved below on paper; the script additionally tests it on a rational grid (a sanity check, not the proof).
- The SOC and old-H capacities (assembly eqs. (3), (4)) are CITED; they are not part of this document's exclusion theorem.

## What each certificate asserts and how it was rebuilt

| # | claim in the PDF | rebuilt as | result |
|---|---|---|---|
| 1 | C = -(ABc + N_j); -Q = (q+b+1+S-p)c + 1 - JS + (q+b)(1+S) + pj, S = a+e, J = j-1 | restriction of (A.2) | EXACT |
| 2 | N - (1 + S + ea) has nonnegative coefficients | direct | EXACT; with ea <= S^2/4 and JS < 1 this gives N_j >= 1 + S - J ea > 1 + 3S/4 (S > 0), N_j >= 1 (S = 0): PROOF-READ. Under (1) p < 1 and JS < 1, so -Q > 0 |
| 3 | (2) det H = -4G^2/Q^2 (G 1313 terms), (3) h-derivative 4 G P6/Q (464 terms), beta^2 - 4 alpha gamma = F+ F- (306 each) | shared core | EXACT |
| 4 | (4) for j > 1, chart (5): gamma, AB beta - N_j alpha, AB F+ + N alpha, AB F- + N alpha; multiplier [j(1+P)]^2 [J(1+V)]^D (L+M)^D, D = max total degree in (a, e) | each a = V L/(J(1+V)(L+M)), e = V M/(J(1+V)(L+M)) monomial cleared jointly, j = 1+J, p = P/(j(1+P)) | EXACT: D = 4, 6, 6, 6; deg_p 2; 3648, 14147, 14123, 14123 terms, all coefficients < 0 |
| 5 | the q = b = V = 0 part of each is -j(j+1)[j(1+P)+P] J^D (L+M)^D | extracted and compared | EXACT, equal in all four |
| 6 | (4) on the wall j = 1, p = P/(1+P), multiplier (1+P)^2 | substitution | EXACT: 228, 491, 479, 479 terms, all < 0, constant -2 |
| 7 | P* (5188 terms), pseudo-division, r1 = L16 M16 (13576 each), 2T = F+ M16^2 + F- L16^2, Z nonzero, norm identity, bordered minors with quotients 2021, 1091 | shared core, full polynomial identities | EXACT |

Data comparison: my gamma and AB beta - N_j alpha tables (chart and wall) equal the stored ones coefficient for coefficient; the F+- certificate tables (chart and wall), the octics and L16, M16 equal the stored pairs up to order and sign.

Region cover. For j > 1: the chart sends J > 0, P, V, L, M, q, b >= 0, L + M > 0 into (1), and onto it: P = pj/(1-pj), V = JS/(1-JS), L : M = a : e (any L + M > 0 if S = 0) reproduce (p, S, a, e) (sympy identity in the script). Positivity: at such a point every monomial is >= 0 and the base part is a sum of monomials with negative coefficients that is strictly negative (J > 0, L + M > 0, j(j+1)(j(1+P) + P) > 0), so the numerator is < 0; the multiplier is > 0 (J > 0, L + M > 0). This covers P = 0, V = 0, L = 0, M = 0, q = 0, b = 0. The wall j = 1 (where (1) reduces to p < 1, a, e arbitrary) is covered by item 6: negative constant, all coefficients negative. So the four inequalities hold strictly on all of (1). EXACT.

## Re-examination of the dependent PROOF-READ steps

- N_j > 0 and (4) give alpha > 0 at every root c >= 0, a unique simple positive root with G_c(c+) > 0, and F+, F- < 0 individually. Correct.
- P*(c+) <= 0 (< 0 off Z = 0), rank 12 via the bordered minors, extension through Z = 0, j = 1 and zero resources (approach inside j > 1 keeping alpha > 0, the strict inequalities and Q < 0), (6) from (2), and the Section 3 contradiction: as in the high_p report; every polynomial input is now EXACT.
- The first auditor's observation that alpha <= 0 at most chamber points (no root at all) is consistent with the certificates; the machinery is only needed where alpha > 0.

## T_40 cone and the four-region cover (assembly `independent_moc4_T39_uniform_20261022`)

EXACT (script): the generators of N_{T_40} are f1, f2, f5, f9, f11, f12, f28, f29, f33, f35 (library, Lemma 4.2), the non-minor ones exactly f28, f35; the three rewriting identities f28 + f12 = f1 + f2 + f11 + f29, f35 + f11 = f1 + f5 + f12 + f33, f28 + f35 = 2f1 + f2 + f5 + f29 + f33 hold on all 24 values; sigma -> sigma^{-1} fixes T_40 and maps the g28 source termwise onto the g35 source. PROOF-READ: the cancellation (subtract min of the two non-minor coefficients by the third identity, then absorb min(d, coefficient of f11) by the second, or min(d, coefficient of f12) by the first) keeps all coefficients nonnegative and ends in a minor sum or in g35 / g28 with d > 0; scaling d = 1 is harmless because {grad C = 0} and the cone are homogeneous.
Cover, PROOF-READ, for nonnegative coefficients with d = 1: if pj >= 1, SOC (CITED). If pj < 1 and j >= 1: either (j-1)(a+e) < 1, excluded by high_j, or (j-1)(a+e) >= 1, which forces a + e > 0 and j >= 1 + 1/(a+e), the old-H capacity (CITED). If pj < 1, j < 1, p < 1: excluded by low_pj. If pj < 1, j < 1, p >= 1: excluded by high_p (p >= 1, pj < 1). The four cases are exhaustive (split on pj >= 1, then on j >= 1, then on p >= 1), including every wall and zero resource. d = 0 is a minor sum. So the three chambers plus the SOC and old-H capacities cover every coefficient configuration of the T_40 cone; this rests on the two CITED capacities, the transpose symmetry (paper Sec. 4.2 group action) and the two-orthogonal reduction, which are outside this replay.

## Verdict

confirmed-modulo-citations. All eight strict certificates, their base parts, the factorizations and identities are EXACT on the whole chamber including j = 1 and every zero resource; analytic steps PROOF-READ. Citations: X subset of {grad C = 0} (paper App. A.1); the SOC and old-H capacities are needed only to turn the exclusion into the support inequality and for the assembly cover.
