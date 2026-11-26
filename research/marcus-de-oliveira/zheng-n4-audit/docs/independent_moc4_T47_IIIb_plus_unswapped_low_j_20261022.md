Evidence level: EXACT for the quadratic-in-c decompositions, the two comparison identities, the chart and all four positivity tables with strict markers; PROOF-READ for the interpolation argument and coverage; CITED for the negative-Phi exclusion (corollary sentence only).

# Audit of `independent_moc4_T47_IIIb_plus_unswapped_low_j_20261022`

Auditor: C22-T47. PDF: `proofs/manin_transfer_20261021/independent_moc4_T47_IIIb_plus_unswapped_low_j_20261022.pdf` (2 pages). Replay: `../replay/independent_moc4_T47_IIIb_plus_unswapped_low_j_20261022__check.py` (output `../replay/independent_moc4_T47_IIIb_plus_unswapped_low_j_20261022__check_output.txt`).

**Verdict: confirmed-modulo-citations.** For the T_65 face on b > 1, a < 1/b, 0 < j <= 1, j <= e < bj/(b-j), 0 <= c < 1/j: Phi >= 0 implies G > 0. Fully replayed EXACT, with all four positive tables (1320/6332/28758/28038) matching; the no-lift conclusion cites the negative-Phi theorem.

## Deviations and limits

- Nothing from the supplement was executed or imported. The radial octic G, the principal minor P and every certificate polynomial were rebuilt from the paper's definitions (f_j of Appendix D via lib moc4lib, the cubic C of (A.2), the height chart w(g,h) = c(g) + h kappa) by the replay script.
- The source vector is c(g) = sum of the eight coefficients times the paper's f_j (0-based Appendix D indices, as the document prints them); d is set to 1 only where the document normalises d = 1.
- Chart polynomials are cleared by multiplying with prod_i den_i^(deg_i G) (each coordinate's denominator to the coordinate degree of G). Where my term counts equal the document's, the cleared polynomial is very likely the same; where they differ, my polynomial is the document's times extra positive denominators, which is enough for the sign argument.
- The sign of G is fixed by the square root; it agrees with the documents' convention, since the documents' positivity tables (low-j, large-pc, stationary) come out positive with my G, with identical term counts. The signs of H6 and H8 are then fixed by the displayed identities (4LG + H6 Phi = K5^2 and the S+- forms), so the hypotheses H6 < 0, H8 < 0 have a well-defined meaning.

## A. Statements

- **Theorem 1** (p. 1; contact set [[0, 4, 5, 6, 9, 12]]). d = 1; b > 1, 0 < a < 1/b, p > 0, q >= 0, 0 < j <= 1, j <= e < bj/(b-j), 0 <= c < 1/j. If Phi = [(p-c)j + (b-c)(p+1+a) - pa]^2 - 4ac(p-b)(p+1+a) >= 0, then G > 0. With the negative-Phi theorem and G = 0 necessity, no realizable cost of the T_65 face lies in this region.

The contact set {0,4,5,6,9,12} is literally paper T_65 (lib CONTACT_SETS[65]; this is historical T47). Its cone N_T is simplicial and its generators are exactly the eight f_j used (lib `generators`). The statements concern realizable costs, that is, x in X with any spectra. No genericity of the spectra is used: radial necessity and Prop. C.6 hold for every x in X.

## B. Dependencies

- [common-lemma] radial necessity (paper p.27): for x = u + iv in X, C(u) = 0 and det Hess C(u) = 0: G = 0 at every spectrally realizable point, via det H(w(g,-C/Q)) = -4 d^2 G^2/Q^2 (replayed EXACT) and C < 0 (EXACT), so Q != 0 and h = -C/Q. Match with the use: checked.
- [other-doc] root_moc4_T47_negative_Phi_20261022 (no realizable point with Phi < 0, coefficients >= 0, e > 0): the no-lift corollary for Phi < 0. Not audited here (CITED).

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | G = alpha c^2 + beta c + gamma; Phi = H c^2 + J c + K with H = (p+j+1+a)^2, K = [pj + b(p+1+a) - pa]^2 | identity | EXACT | flint |
| 2 | HG = alpha Phi + U0 + c U1 and KG = gamma Phi + c(V0 - c U0) | identity | EXACT | flint |
| 3 | alpha >= 0, gamma > 0, V0 >= 0, EJ = jU0 + U1 > 0 imply G > 0 on 0 <= c < 1/j when Phi >= 0 (two cases on the sign of U0) | inequality | PROOF-READ | interpolation of the affine remainder; K > 0 when U0 < 0 |
| 4 | chart b = 1+B, a = V/(b(1+V)), p = P, q = Q, j = T/(T+W), e = T[Z(1+U)+TU]/((T+W)Z(1+U)) covers the region exactly (j = 1 at W = 0) | identity | EXACT | sympy: e = j + j^2 U/((b-j)(1+U)), limit bj/(b-j) |
| 5 | cleared alpha, gamma, V0, EJ have 1320, 6332, 28758, 28038 coefficients, all positive | inequality | EXACT | flint; counts identical to the document |
| 6 | gamma and EJ keep positive monomials free of W, U, Q (84 and 270) | inequality | EXACT | flint |
| 7 | no-lift: G = 0 necessity plus negative-Phi exclusion | citation | CITED | root_moc4_T47_negative_Phi |

## D. Falsification

Not needed: the theorem's implication is replayed EXACT. Realizable samples (doc 2 test) satisfy G = 0, consistent with this.

## E. Verdict and gaps

Verdict: **confirmed-modulo-citations**.

Gaps: none found.

Exact replays: alpha, beta, gamma and H, J, K; both comparison identities; chart map identities; four positive tables 1320/6332/28758/28038 and strict markers.
