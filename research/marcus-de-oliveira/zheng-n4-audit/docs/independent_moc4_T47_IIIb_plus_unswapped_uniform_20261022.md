Evidence level: EXACT for the cross-range chart positivity (own clearing), the stationary chart and its two-guard certificate (weights transcribed from the companion verifier as data, identity and signs recomputed), and the coverage identities; PROOF-READ for the assembly and exhaustiveness; CITED for the merged-C high-j theorem, the negative-Phi theorem, and the SOC, old-H and merged-C capacities (foundations).

# Audit of `independent_moc4_T47_IIIb_plus_unswapped_uniform_20261022`

Auditor: C22-T47. PDF: `proofs/manin_transfer_20261021/independent_moc4_T47_IIIb_plus_unswapped_uniform_20261022.pdf` (3 pages). Replay: `../replay/independent_moc4_T47_IIIb_plus_unswapped_uniform_20261022__check.py` (output `../replay/independent_moc4_T47_IIIb_plus_unswapped_uniform_20261022__check_output.txt`).

**Verdict: confirmed-modulo-citations.** Assembles the support inequality on the whole unswapped IIIb-plus region of T_65 (b > 1, a < 1/b, p, c > 0, 0 < j <= e, j > eb/(e+b)). The new cross-range and stationary-range certificates were replayed EXACT, with stationary counts matching. Correct modulo the cited high-j theorem, the negative-Phi theorem and foundation atoms.

## Deviations and limits

- Nothing from the supplement was executed or imported. The radial octic G, the principal minor P and every certificate polynomial were rebuilt from the paper's definitions (f_j of Appendix D via lib moc4lib, the cubic C of (A.2), the height chart w(g,h) = c(g) + h kappa) by the replay script.
- The source vector is c(g) = sum of the eight coefficients times the paper's f_j (0-based Appendix D indices, as the document prints them); d is set to 1 only where the document normalises d = 1.
- Chart polynomials are cleared by multiplying with prod_i den_i^(deg_i G) (each coordinate's denominator to the coordinate degree of G). Where my term counts equal the document's, the cleared polynomial is very likely the same; where they differ, my polynomial is the document's times extra positive denominators, which is enough for the sign argument.
- The sign of G is fixed by the square root; it agrees with the documents' convention, since the documents' positivity tables (low-j, large-pc, stationary) come out positive with my G, with identical term counts. The signs of H6 and H8 are then fixed by the displayed identities (4LG + H6 Phi = K5^2 and the S+- forms), so the hypotheses H6 < 0, H8 < 0 have a well-defined meaning.
- Cross range: my cleared polynomial has 2,009,232 terms, against the document's 1,232,529. I clear the nested rational map with a different common denominator, so mine is the document's polynomial times extra positive factors. All its coefficients are positive, with 69552 positive monomials at U = Q = 0. That proves G > 0 on the chart independently of their table.
- Stationary certificate: the 11 positive rational weights and monomials (common denominator 42) are not printed in the PDF. I transcribed them as data from the companion script by reading it, without executing it. With my own G~ (134717 terms, 24 negative, as the document states) and my own reduced guard numerators, R0 = 42 G~ - M_pc F_pc - M_Phi F_Phi has 134706 coefficients, all positive, and 16827 monomials free of U, Q. Both counts match the document.
- The projective endpoint table at je = A0^2 (142167 terms) was not replayed. It is not needed: the stationary chart, at U = 0, contains je = A0^2, and there c = theta(S)/(j + A0) agrees with C* = 1/(j + A0).

## A. Statements

- **Theorem 1 (whole)** (p. 1; contact set [[0, 4, 5, 6, 9, 12]]). T_65 face, d = 1: b > 1, 0 < a < 1/b, p, c > 0, q >= 0, 0 < j <= e, j > eb/(e+b). Every x in X with c(x) of this form has Phi_{c(x)} >= 0 on SO(4). The boundary j = eb/(e+b) and the faces a, p, c, e = 0, b = 1, ab = 1 are excluded and stated to be separate inputs.

The contact set {0,4,5,6,9,12} is literally paper T_65 (lib CONTACT_SETS[65]; this is historical T47). Its cone N_T is simplicial and its generators are exactly the eight f_j used (lib `generators`). The statements concern realizable costs, that is, x in X with any spectra. No genericity of the spectra is used: radial necessity and Prop. C.6 hold for every x in X.

## B. Dependencies

- [common-lemma] radial necessity (paper p.27): for x = u + iv in X, C(u) = 0 and det Hess C(u) = 0: G = 0 at every spectrally realizable point, via det H(w(g,-C/Q)) = -4 d^2 G^2/Q^2 (replayed EXACT) and C < 0 (EXACT), so Q != 0 and h = -C/Q. Match with the use: checked.
- [core-doc] Theorem 1.2 / Corollary 2.2 (audited in foundations): transfer of the orthogonal support inequality to unitary matrices. Not audited here (CITED).
- [common-lemma] SOC (45,1;35): cj >= 1, and old-H (aux 9, main 1; 35): j >= 1 + 1/p (foundations/so4-atoms.md; orbit and roles re-checked in the I_II replay): pointwise covers. Match with the use: checked.
- [common-lemma] merged coherent-C capacity C*(A0,j,e) (foundations/so4-atoms-part2b.md G4: PROVED IN THE DOCUMENT AND CHECKED): c >= C* cover. Match with the use: checked.
- [other-doc] independent_moc4_T47_IIIb_plus_unswapped_low_j_20261022 (audited in this batch): 0 < j <= 1. Match with the use: checked.
- [other-doc] independent_moc4_T47_IIIb_plus_large_pc_20261022 (audited in this batch): stationary range with pc >= a^2. Match with the use: checked.
- [other-doc] independent_moc4_T47_merged_C_high_j_uniform_20261022: j >= b. Not audited here (CITED).
- [other-doc] root_moc4_T47_negative_Phi_20261022: Phi < 0. Not audited here (CITED).

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | covers: old-H j >= 1 + 1/p, SOC cj >= 1, merged-C c >= C* | citation | CITED | foundations atoms (proved) |
| 2 | 0 < j <= 1: after SOC (c >= 1/j), the low-j theorem applies; its region matches (j <= 1 < b makes j > eb/(e+b) equivalent to e < bj/(b-j)) | citation | CITED | doc 4 of this batch (EXACT); match PROOF-READ |
| 3 | j >= b: merged-C high-j theorem | citation | CITED | not in this batch |
| 4 | 1 < j < b hard part: j = 1+T, b = j+A, a = theta(V)/b, p = theta(P)/T covers p < 1/(j-1), a < 1/b | identity | PROOF-READ | direct |
| 5 | cross range je <= A0^2: R = 1 + (A0/j - 1) theta(U), e = jR^2, c = C* theta(S) with C* = (A0+e)/(A0 j (1+R)^2) rational | identity | EXACT | sympy: sqrt(je) = jR; R -> A0/j |
| 6 | cross range: cleared G has only positive coefficients with strict markers at U = Q = 0, so G > 0 | inequality | EXACT | flint, 2,009,232 terms |
| 7 | stationary range je >= A0^2: e = A0^2/j + U, c = theta(S)/(j + A0); cleared G has 134717 terms, 24 negative | identity | EXACT | flint; counts match |
| 8 | stationary, pc >= a^2: the large-pc theorem (its hypotheses j > 1, j < 1 + 1/p, plus condition, pc >= a^2 hold) | citation | CITED | doc 3 of this batch |
| 9 | stationary, pc <= a^2, Phi >= 0: 42 G~ = R0 + M_pc F_pc + M_Phi F_Phi with R0 > 0 and positive multipliers and guard denominators | inequality | EXACT | recomputed; weights transcribed as data |
| 10 | stationary, Phi < 0: negative-Phi theorem | citation | CITED | other document |
| 11 | exhaustiveness of the assembly, including the interfaces j = 1, j = b, je = A0^2, e = j, q = 0, pc = a^2 | case split | PROOF-READ | each interface lies in a closed chart or a theorem with inclusive bound |

## D. Falsification

The stationary branch was also checked numerically (3000 samples): G stays well away from 0. No violation.
- NUMERIC: stationary range, pc <= a^2, Phi >= 0: 3000 log-uniform chart samples (seed 11), min G/scale = 0.31 > 0

## E. Verdict and gaps

Verdict: **confirmed-modulo-citations**.

Gaps: none found.

Exact replays: G (965 terms); cross-range chart: all 2,009,232 cleared coefficients positive, markers at U = Q = 0; stationary chart: 134717 terms, 24 negative; stationary two-guard certificate R0 (134706 positive, 16827 markers); chart identities for C* and R.

Unchecked: projective endpoint table at je = A0^2 (not needed, see deviations); merged-C high-j theorem (other document).
