Evidence level: EXACT for all displayed identities and the rebuilt radial factor; PROOF-READ for the smoothness and contradiction argument; CITED for the predecessor theorem and its sign mechanism, the HMSV set equality (C.2) behind the rank bound, Corollary C.5 and Proposition C.6.

# Audit of `final_referee_T46_two_Heron_necessity_20261023`

Auditor: C22-T47. PDF: `proofs/manin_transfer_20261021/final_referee_T46_two_Heron_necessity_20261023.pdf` (2 pages). Replay: `../replay/final_referee_T46_two_Heron_necessity_20261023__check.py` (output `../replay/final_referee_T46_two_Heron_necessity_20261023__check_output.txt`).

**Verdict: confirmed-modulo-citations.** Proves H6 >= 0 and H8 >= 0 at realizable points on one 8-generator face of the T_47 cone. Every identity replayed EXACT; the logic is sound; it rests on the preceding negative-j-Heron theorem and its principal-minor mechanism, the HMSV set equality and Prop. C.6, all CITED.

## Deviations and limits

- Nothing from the supplement was executed or imported. The radial octic G, the principal minor P and every certificate polynomial were rebuilt from the paper's definitions (f_j of Appendix D via lib moc4lib, the cubic C of (A.2), the height chart w(g,h) = c(g) + h kappa) by the replay script.
- The source vector is c(g) = sum of the eight coefficients times the paper's f_j (0-based Appendix D indices, as the document prints them); d is set to 1 only where the document normalises d = 1.
- Chart polynomials are cleared by multiplying with prod_i den_i^(deg_i G) (each coordinate's denominator to the coordinate degree of G). Where my term counts equal the document's, the cleared polynomial is very likely the same; where they differ, my polynomial is the document's times extra positive denominators, which is enough for the sign argument.
- The sign of G is fixed by the square root; it agrees with the documents' convention, since the documents' positivity tables (low-j, large-pc, stationary) come out positive with my G, with identical term counts. The signs of H6 and H8 are then fixed by the displayed identities (4LG + H6 Phi = K5^2 and the S+- forms), so the hypotheses H6 < 0, H8 < 0 have a well-defined meaning.
- In this frame the Hessian determinant identity is det H(w(g,-C/Q)) = -4 G^2/Q^2 with prefactor 1. G has degree 9 and 1556 terms. The document prints neither number; the identity agrees with its second-derivative formula -8 G_x^2/Q^2.
- The 'signed radial matrix mechanism' (principal norm congruences in the T46 frame) is in the predecessor final_referee_T46_negative_j_Heron_uniform_20261023. I read that document's text only far enough to confirm that its mechanism's hypotheses (G = 0, H6 < 0, H8 < 0, q + d > 0, C < 0, Q != 0) match the use made here. I did not replay its congruences.
- No numerical falsification was run in this frame (time).

## A. Statements

- **Theorem 1** (p. 1; contact set [[1, 3, 5, 6, 7, 8]]). Let T = {1,3,5,6,7,8} = paper T_47 (literal). Let c = q f27 + p f1 + j f59 + c f5 + e f2 + d f28 + b f6 + a f4 with all coefficients >= 0 and b, j, d > 0. This is one simplicial 8-generator face of N_T, which has 10 generators; it omits f0 = n1,4 and f3 = n14,34. If c = c(x) for some x in X (with the restored real height), then H6 >= 0 and H8 >= 0, where disc_j G = (q+d)^2 H6 H8 and H6 = T0^2 + 4bDM(q-c)(q-e).

The contact set {1,3,5,6,7,8} is literally paper T_47 (historical T46). Its cone has 10 generators; the document works on the simplicial face spanned by 8 of them.

## B. Dependencies

- [common-lemma] radial necessity (paper p.27): for x = u + iv in X, C(u) = 0 and det Hess C(u) = 0: G = 0 at every spectrally realizable point, via det H(w(g,-C/Q)) = -4 d^2 G^2/Q^2 (replayed EXACT) and C < 0 (EXACT), so Q != 0 and h = -C/Q. Match with the use: checked.
- [other-doc] final_referee_T46_negative_j_Heron_uniform_20261023 (Theorem 1: no realizable point with H6 < 0; its signed principal-minor mechanism): excludes H6 < 0; makes all 12x12 principal minors <= 0 at radial points with H6, H8 < 0. Not audited here (CITED).
- [external] HMSV via paper (C.2): {C = det H = 0} = Sigma2-hat (class (b) in foundations/appendix-c-open-items.md, which lists this document): rank <= 12 at non-realizable points of the radial graph. Match with the use: checked.
- [core-doc] Corollary C.5 (rank H <= 12 on Sigma2-hat) and Proposition C.6 (12x12 principal minors of H(Re x) >= 0; one positive if rank 12), foundations/appendix-c-part1.md: the final contradiction. Match with the use: checked.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | the 8 sources vanish exactly on {1,3,5,6,7,8} = T_47; they are independent; generators(T_47) has 10 elements | identity | EXACT | lib F_ZERO, generators, rank |
| 2 | height equation affine in h; constant = -(q+d)L - jM < 0 | identity | EXACT | own expansion of C(c+h kappa) |
| 3 | det H(w(g,-C/Q)) = -4G^2/Q^2, G irreducible (1556 terms, degree 9) | identity | EXACT | 14x14 Bareiss determinant, exact square root |
| 4 | disc_j G = (q+d)^2 H6 H8 (150- and 279-term factors) | identity | EXACT | flint factor |
| 5 | (1) H6 = T0^2 + 4bDM(q-c)(q-e); (H6)_p = 2LT0; (H6)_q = 2(bD-M)T0 + 4bDM(2q-c-e) | identity | EXACT | polynomial identities |
| 6 | case analysis: grad H6 != 0 on {H6 = 0} except on p=q=c=e=t | case split | PROOF-READ | b, D, M, L > 0 from d > 0, b > 0 |
| 7 | on p=q=c=e=t: H8 = W0^2 + 4dDb^2 L t(t+d+b) with W0(0) = -bdL, so H8 > 0 there | identity | EXACT | square root of the restricted difference |
| 8 | at a realizable point with H8 < 0: G = 0, so G_j^2 = (q+d)^2 H6 H8 gives H6 <= 0; H6 < 0 is excluded by the predecessor; so H6 = 0 and G_j = 0 | inequality | CITED | PROOF-READ, plus the predecessor theorem |
| 9 | grad disc = (q+d)^2 H8 grad H6 != 0, hence A != 0 and grad G = -grad(disc)/(4A) != 0 | analytic | PROOF-READ | argument checked line by line |
| 10 | gcd(G, G_j) = 1, so G_j != 0 is dense on the local smooth real piece of {G = 0} | analytic | PROOF-READ | gcd and irreducibility EXACT; density argument completed by the auditor (a real 7-dim analytic piece is Zariski dense in the irreducible G = 0) |
| 11 | at dense nearby points H6, H8 < 0; the predecessor mechanism gives all 12x12 principal minors <= 0 there, hence at the putative point by continuity | citation | CITED | predecessor mechanism, not replayed |
| 12 | some G_x != 0, so d^2/dx^2 det = -8 G_x^2/Q^2 < 0 gives corank <= 2; rank <= 12 gives rank exactly 12 | analytic | CITED | PROOF-READ; rank <= 12 needs the HMSV set equality |
| 13 | Prop. C.6 gives a positive 12x12 principal minor, a contradiction | citation | CITED | foundations Appendix C |

Scope note: the theorem is about one simplicial face (8 of the 10 generators of N_{T_47}). Whether the T_47 assembly needs other faces is a question for the T_47 case audit. Here, as in the predecessor, Theorem 1 is a necessity statement (H6, H8 >= 0 at realizable points), not a support inequality, as the document says.

## D. Falsification

Not run in this frame (time). The analogous tests in the T_65 frame, for the companion negative-p Heron document, found no violation: realizable points satisfy G = 0 to 1e-13, and the numerical corank of H at random points of G = 0 is 2.

## E. Verdict and gaps

Verdict: **confirmed-modulo-citations**.

Gaps: none found.

Exact replays: contact set and generator check; C(c + h kappa) = C + hQ with C = -(q+d)L - jM; det identity and G (1556 terms, degree 9), G irreducible; disc_j G = (q+d)^2 H6 H8; identity (1) and both gradient formulas; exceptional-locus H8 formula and W0(0) = -bdL; gcd(G,G_j) = 1.

Unchecked: the predecessor's principal norm congruences in the T46 frame (cited, not replayed); no numerical falsification in this frame.
