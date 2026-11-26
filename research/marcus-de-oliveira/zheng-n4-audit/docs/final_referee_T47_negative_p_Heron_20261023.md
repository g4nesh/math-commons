Evidence level: EXACT for every displayed polynomial identity (G, H6, H8, P, W18, R8, T8, both congruences, K5, the Phi identities, both gcds); PROOF-READ for the density, rank and sign argument; CITED for the HMSV set equality (C.2) behind the rank bound at non-realizable points, Corollary C.5, Proposition C.6 and the negative-Phi exclusion; NUMERIC support for the rank bound and a falsification sample.

# Audit of `final_referee_T47_negative_p_Heron_20261023`

Auditor: C22-T47. PDF: `proofs/manin_transfer_20261021/final_referee_T47_negative_p_Heron_20261023.pdf` (3 pages). Replay: `../replay/final_referee_T47_negative_p_Heron_20261023__check.py` (output `../replay/final_referee_T47_negative_p_Heron_20261023__check_output.txt`).

**Verdict: confirmed-modulo-citations.** Excludes realizable points of the T_65 face with H6 < 0 and H8 < 0 (Thm 1), and with H6 < 0 when p >= b, a, c, d, e > 0 (Cor 2). Every polynomial identity was replayed EXACT, including both full congruences; the argument is sound modulo the HMSV set equality, Prop. C.6 and the negative-Phi theorem (CITED).

## Deviations and limits

- Nothing from the supplement was executed or imported. The radial octic G, the principal minor P and every certificate polynomial were rebuilt from the paper's definitions (f_j of Appendix D via lib moc4lib, the cubic C of (A.2), the height chart w(g,h) = c(g) + h kappa) by the replay script.
- The source vector is c(g) = sum of the eight coefficients times the paper's f_j (0-based Appendix D indices, as the document prints them); d is set to 1 only where the document normalises d = 1.
- Chart polynomials are cleared by multiplying with prod_i den_i^(deg_i G) (each coordinate's denominator to the coordinate degree of G). Where my term counts equal the document's, the cleared polynomial is very likely the same; where they differ, my polynomial is the document's times extra positive denominators, which is enough for the sign argument.
- The sign of G is fixed by the square root; it agrees with the documents' convention, since the documents' positivity tables (low-j, large-pc, stationary) come out positive with my G, with identical term counts. The signs of H6 and H8 are then fixed by the displayed identities (4LG + H6 Phi = K5^2 and the S+- forms), so the hypotheses H6 < 0, H8 < 0 have a well-defined meaning.
- The document labels this T47, which is historical T47 = paper T_65 (contacts {0,4,5,6,9,12}, literally T_65 in Appendix E).
- I could not compare the sign of H6, H8 with the companion JSON, which I did not read. I fixed the signs by the document's own identities: 4LG + H6 Phi = K5^2 holds with H6 as factored, and S+ = H8 R8^2, S- = d^2 H6 T8^2 then hold for one sign of W18.

## A. Statements

- **Theorem 1** (p. 1; contact set [[0, 4, 5, 6, 9, 12]]). T = {0,4,5,6,9,12} = T_65. c = q f2 + p f9 + j f1 + c f45 + e f5 + d f35 + b f0 + a f7 (N_T is simplicial with exactly these 8 generators), coefficients >= 0, d > 0. If H6 < 0 and H8 < 0, where B^2 - 4AC0 = H6 H8 for G = Ap^2 + Bp + C0, then no x in X has c(x) = c with the restored real height. Zero coefficients and Phi = 0 are allowed.
- **Corollary 2** (p. 2; contact set [[0, 4, 5, 6, 9, 12]]). Same face, coefficients >= 0, a, c, d, e > 0, p >= b. If H6 < 0, there is no realizable lift. This includes H8 = 0, A = 0 and Phi = 0.

The contact set {0,4,5,6,9,12} is literally paper T_65 (lib CONTACT_SETS[65]; this is historical T47). Its cone N_T is simplicial and its generators are exactly the eight f_j used (lib `generators`). The statements concern realizable costs, that is, x in X with any spectra. No genericity of the spectra is used: radial necessity and Prop. C.6 hold for every x in X.

## B. Dependencies

- [common-lemma] radial necessity (paper p.27): for x = u + iv in X, C(u) = 0 and det Hess C(u) = 0: G = 0 at every spectrally realizable point, via det H(w(g,-C/Q)) = -4 d^2 G^2/Q^2 (replayed EXACT) and C < 0 (EXACT), so Q != 0 and h = -C/Q. Match with the use: checked.
- [external] HMSV via paper (C.2): {C = det H = 0} = Sigma2-hat (class (b); this document is one of the four listed in foundations/appendix-c-open-items.md Sec. 4.2): rank <= 12 at points of the perturbed radial graph that are not known to be realizable. Match with the use: checked.
- [core-doc] Corollary C.5 and Proposition C.6 (foundations/appendix-c-part1.md: proved, C.6 modulo standard GIT): rank <= 12 on Sigma2-hat; nonnegative 12x12 principal minors at H(Re x), one positive at rank 12. Match with the use: checked.
- [other-doc] root_moc4_T47_negative_Phi_20261022 (no realizable point with Phi < 0 for coefficients >= 0, e > 0): Corollary 2: excludes Phi < 0. Not audited here (CITED).

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | C(w(g,h)) = C + hQ with C = -d{ab+ad+aq+bd+bp+d^2+de+dp+dq+ep+eq+pq} < 0 | identity | EXACT | own expansion |
| 2 | det H(w(g,-C/Q)) = -4d^2 G^2/Q^2; G has 965 terms, degree 8, is irreducible | identity | EXACT | 14x14 Bareiss, exact square root |
| 3 | G is quadratic in p; B^2 - 4AC0 = H6 H8 (123 and 260 terms, primitive) | identity | EXACT | flint factor |
| 4 | principal minor P deleting coordinates 1, 6 is height-free; 3376 terms; degree 2 in p | identity | EXACT | 12x12 Bareiss |
| 5 | u0^2 - u1(B gamma - beta C0) = W18^2 (norm 620533 terms, W18 15999 terms, degree 18) | identity | EXACT | exact square root |
| 6 | S+ = H8 R8^2, S- = d^2 H6 T8^2 (465, 595 terms); W18, R8, T8 free of p | identity | EXACT | exact division and square roots |
| 7 | P S+- = (u1 p + u0 +- W18)^2 + lambda+- G with p-free lambda (150891 and 159898 terms) | identity | EXACT | full polynomial congruences |
| 8 | on G = 0 with W18 R8 T8 != 0 and H6, H8 < 0: P < 0 | inequality | PROOF-READ | S+- < 0 and two squares; P = 0 would force W18 = 0 |
| 9 | H6 H8 > 0 gives G_p = +-sqrt(disc) != 0 at real roots (including A = 0); local graph over the other 7 coordinates | analytic | PROOF-READ | implicit function theorem |
| 10 | density of W18 R8 T8 != 0 on the graph | analytic | PROOF-READ | argument supplied by auditor: all three are nonzero polynomials free of p, and the graph projects onto an open set of the 7 free coordinates |
| 11 | d^2/dp^2 det = -8 d^2 G_p^2/Q^2 < 0 gives corank <= 2; with rank <= 12 the rank is exactly 12 | analytic | CITED | PROOF-READ; rank <= 12 at non-realizable points needs the HMSV set equality |
| 12 | a rank-12 real symmetric matrix has all nonzero 12x12 principal minors of one sign (X B X^t and Cauchy-Binet) | identity | PROOF-READ | standard; checked |
| 13 | continuity: all 12x12 principal minors <= 0 at the putative point; Prop. C.6 gives one > 0 at rank 12; contradiction | citation | CITED | Prop. C.6 |
| 14 | Cor 2: 4LG + H6 Phi = K5^2 (K5 193 terms); L > 0 | identity | EXACT | exact square root |
| 15 | Cor 2: on G = 0 and H6 < 0, Phi <= 0; Phi < 0 is excluded by the cited theorem (e > 0), so Phi = K5 = 0 and 4L grad G = -H6 grad Phi | citation | CITED | PROOF-READ plus root_moc4_T47_negative_Phi |
| 16 | Cor 2: grad Phi != 0 (cases p > b with p != c, p = c, p = b): Phi_j = 2X(p-c), Phi/c=p = R^2, Phi_c/c=p, Phi/p=b = Xb^2, Phi_p/p=b | identity | EXACT | sympy identities; the case logic PROOF-READ |
| 17 | Cor 2: gcd(G, G_p) = 1, so G_p != 0 is dense; H6 H8 = G_p^2 there gives H8 < 0; Theorem 1 mechanism; G_x != 0 gives rank 12; contradiction | analytic | PROOF-READ | as steps 8-13 |

The phrase 'the exact nonzero quadratic resultant W18^2 prevents a common principal/radial component' is not needed. Density follows because W18, R8 and T8 are free of p (EXACT), while the local radial graph is a graph over the remaining seven coordinates. I record this as an argument supplied by the auditor, not as a gap.

## D. Falsification

Realizable samples satisfy G = 0 to 2e-13. This is an independent check of the radial factor and of radial necessity. None of the 15 samples falls in the excluded regions. The rank bound used at non-realizable points held numerically at 8 random points of G = 0. No violation found.
- NUMERIC: 15 realizable points on the T_65 stratum (lib sampler, seeds 5000+, margin 1e-2): |G|/scale <= 2e-13; none has H6 < 0 and H8 < 0; none lies in the Cor. 2 region with H6 < 0
- NUMERIC: 8 random real points of G = 0 (seed 7): sigma_13/sigma_1 <= 1.6e-16, i.e. numerical corank >= 2, consistent with the cited rank bound

## E. Verdict and gaps

Verdict: **confirmed-modulo-citations**.

Gaps: none found.

Exact replays: C constant; det identity with prefactor d^2, G (965 terms); disc_p G = H6 H8 (123/260); P (3376 terms, height-free, quadratic in p); W18^2 identity (620533-term norm, 15999-term W18); S+ = H8 R8^2, S- = d^2 H6 T8^2 (465/595); both full congruences mod G; 4LG + H6 Phi = K5^2 (193); Phi identities of Cor 2; gcd(G,G_p) = 1, G irreducible.

Unchecked: the HMSV set equality itself (external); the negative-Phi exclusion theorem (another document).
