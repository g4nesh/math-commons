Evidence level: EXACT for the atom orbit and role checks, the capacity algebra (2), Q, the two sector-I ab <= 1 charts, and the ab = 1 wall (a < 1) negative-fibre structure; PROOF-READ for the ab = 1 limit argument and the comparison-chart logic; CITED for the T24 coherent core (3) (foundations part 2a), the negative-Phi theorem and the unitary transfer; UNCHECKED for the remaining 16 positive charts, the last-II comparison tables, 18 marker tables, 20 projective corners and the ab = 1, a = 1+A tables (time).

# Audit of `independent_moc4_T47_I_II_uniform_20261022`

Auditor: C22-T47. PDF: `proofs/manin_transfer_20261021/independent_moc4_T47_I_II_uniform_20261022.pdf` (4 pages). Replay: `../replay/independent_moc4_T47_I_II_uniform_20261022__check.py` (output `../replay/independent_moc4_T47_I_II_uniform_20261022__check_output.txt`).

**Verdict: could-not-check.** Support inequality on the T_65 face for sector I (e <= j) and sector II (ab >= 1, e >= j). The capacity inputs, Q, 2 of 18 positive charts (19536/43392, matching) and the ab = 1 wall structure (9233 terms, 4 negatives, the four listed fibre polynomials) were replayed EXACT. The remaining 16 charts, the last comparison chart and the marker and corner tables are UNCHECKED, so no overall judgement; no discrepancy found.

## Deviations and limits

- Nothing from the supplement was executed or imported. The radial octic G, the principal minor P and every certificate polynomial were rebuilt from the paper's definitions (f_j of Appendix D via lib moc4lib, the cubic C of (A.2), the height chart w(g,h) = c(g) + h kappa) by the replay script.
- The source vector is c(g) = sum of the eight coefficients times the paper's f_j (0-based Appendix D indices, as the document prints them); d is set to 1 only where the document normalises d = 1.
- Chart polynomials are cleared by multiplying with prod_i den_i^(deg_i G) (each coordinate's denominator to the coordinate degree of G). Where my term counts equal the document's, the cleared polynomial is very likely the same; where they differ, my polynomial is the document's times extra positive denominators, which is enough for the sign argument.
- The sign of G is fixed by the square root; it agrees with the documents' convention, since the documents' positivity tables (low-j, large-pc, stationary) come out positive with my G, with identical term counts. The signs of H6 and H8 are then fixed by the displayed identities (4LG + H6 Phi = K5^2 and the S+- forms), so the hypotheses H6 < 0, H8 < 0 have a well-defined meaning.
- Out of time: I replayed 2 of the 18 positive-coefficient charts and one of the two ab = 1 wall charts. Every other chart table is UNCHECKED. The verdict could-not-check reflects incompleteness, not a found problem.

## A. Statements

- **Theorem 1** (p. 2; contact set [[0, 4, 5, 6, 9, 12]]). T_65 face, d = 1, a, p, e, c > 0, q, b, j >= 0. Every realizable cost with I: e <= j, or II: ab >= 1 and e >= j, has Phi >= 0 on SO(4). It includes q = 0, b = 0 where compatible, and j = 0 in II; it excludes the faces a, p, e, c = 0. The region ab < 1, e > j is not covered.

The contact set {0,4,5,6,9,12} is literally paper T_65 (lib CONTACT_SETS[65]; this is historical T47). Its cone N_T is simplicial and its generators are exactly the eight f_j used (lib `generators`). The statements concern realizable costs, that is, x in X with any spectra. No genericity of the spectra is used: radial necessity and Prop. C.6 hold for every x in X.

## B. Dependencies

- [common-lemma] radial necessity (paper p.27): for x = u + iv in X, C(u) = 0 and det Hess C(u) = 0: G = 0 at every spectrally realizable point, via det H(w(g,-C/Q)) = -4 d^2 G^2/Q^2 (replayed EXACT) and C < 0 (EXACT), so Q != 0 and h = -C/Q. Match with the use: checked.
- [core-doc] Theorem 1.2 / Corollary 2.2 (audited in foundations): transfer of the orthogonal support inequality to unitary matrices. Not audited here (CITED).
- [common-lemma] SOC (45,1;35), old-H (9,1;35), chain (1;7,9;35), four (0,7;5;35): foundations/so4-atoms.md; orbit and roles re-checked EXACT here: capacities (1), (2). Match with the use: checked.
- [common-lemma] T24 coherent core transported by (59,53,55,67,54,2) -> (45,5,7,0,35,1) (foundations/so4-atoms-part2a.md; transport realised by a group element, EXACT here): capacity (3). Match with the use: checked.
- [other-doc] root_moc4_T47_negative_Phi_20261022: Phi < 0 in the comparison charts. Not audited here (CITED).

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | C < 0 and Q formula; G = 0 necessary | identity | EXACT | replayed |
| 2 | capacities (1): SOC, old-H, chain tuples are in the proved orbits with the printed roles (old-H reversed roles not in orbit) | citation | EXACT | group orbit computation with lib |
| 3 | capacity (2) is the four-coefficient atom with (a,b,c,d) = (b,a,e,1) | identity | EXACT | orbit check and sympy |
| 4 | capacity (3) is the transported T24 coherent core | citation | CITED | part 2a; transport map EXACT |
| 5 | ab = 1 closure c >= 1/e as a limit of valid pointwise inequalities | analytic | PROOF-READ | b -> 1/a from above, Emax -> infinity |
| 6 | sector I, ab <= 1, a = A/(1+A) and a = 1+A: cleared G positive (19536, 43392) | inequality | EXACT | flint; counts match |
| 7 | sector I, ab > 1 (8 charts) and sector II small/large p (8 charts) | inequality | UNCHECKED | time |
| 8 | last II chart: HG = alpha Phi + U0 + c U1 interpolation (the identity is the one replayed for the low-j document); tables 6673/50347/107703 | inequality | UNCHECKED | logic PROOF-READ, tables not replayed |
| 9 | ab = 1 wall, a = A/(1+A): 9233 terms, 4 negative; their U-fibres are exactly 4U^2-U+13, 8U^2-2U+6, 5U^2-U+1, U^3+7U^2-U+2, each positive on U >= 0 | inequality | EXACT | flint fibres; sympy real roots |
| 10 | ab = 1 wall, a = 1+A, with its three tables; 18 marker tables; 20 projective corners | inequality | UNCHECKED | time |

## D. Falsification

Not run (time). The realizable samples from the doc 2 test all satisfy G = 0, consistent with the charts.

## E. Verdict and gaps

Verdict: **could-not-check**.

Gaps: none found.

Exact replays: atom orbit and role checks (SOC, old-H, chain, four, T24 transport); four-atom capacity algebra = (2); Q formula; sector I ab <= 1 charts (19536, 43392 positive); ab = 1 wall (a < 1) fibre structure.

Unchecked: 16 of the 18 positive charts; last-II comparison tables; 18 strict-marker tables and the j = 0 markers; 20 projective corner endpoints; ab = 1, a = 1+A tables.
