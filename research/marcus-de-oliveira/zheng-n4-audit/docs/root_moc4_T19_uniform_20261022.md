# Audit: root_moc4_T19_uniform_20261022 (paper case T_20)

Evidence level: EXACT for every displayed identity, the node parametrization, the generator list, the minor-cone dual, the sign-cell analysis and the root-location values. PROOF-READ for the exhaustiveness of the chart, the H = 0 argument and the interlacing/residue sign argument. CITED for Appendix B, the two-orthogonal reduction and (4.2). The real-v layer rests on an argument the auditor supplies (PROOF-READ plus an EXACT generator check). The falsification run is NUMERIC.

Document: `proofs/manin_transfer_20261021/root_moc4_T19_uniform_20261022.pdf` (3 pages). Supplement label T19 = paper T_20 = {0, 1, 5, 6, 7, 8}, the same set literally. The supplement catalogue also lists this document as used for T_24; no T_24 claim is made in this document, and none was audited here.

## Deviations and limits

- These follow the conventions of the T_18 audit (root_moc4_T17_uniform_20261022.md): own code, the shared library, gaps Re{kappa(F_j - F_0)}, and imaginary parts as identities in (q, qb).
- The citation "three-A/three-C common-circle theorem" for the case of real v is not a result of the paper. The paper's Theorem C.1 needs a whole spectrum on one circle. The auditor supplies the argument (S9).
- The companion checker that the document names was not run or read.

## A. Statements

**Theorem 1 (p. 1).** Take x in X with P c(x) >= 0 and (P c(x))_pi = 0 for pi in T_20 = [0, 1, 5, 6, 7, 8]. Extra contacts are allowed, and so are all finite spectra, with or without repeats and with or without a_i + b_j = 0. Then Phi_{c(x)}(O) >= 0 for all O in SO(4) (hence in O(4) and U(4)). Stronger claim: on the eight-distinct layer with v nonreal, c(x) lies in the cone of the nine minor generators. This is at least (*) for T_20, in full.

Intermediate claims:
- (1) The minor generators {0,...,6,27,33} span the same 8-dimensional space as all eleven generators. Their cone has exactly 10 facets: g_j >= 0 for j in {3,4,9,10,14,16,19,22}, A = g3 - g2 + g18 and B = g13 - g3.
- (2)-(3) The l-equation is the curve (2). Its node at (v,w) = (1,1) gives a rational parametrization by the slope q. The coefficient qT(q) is nonzero on original configurations.
- (4) Four sign kernels. (5) The Pick identity. (6) The formula for Im w / Im q. (7) The supported nonreal layer lies in 0<h<1, k<0, l>1, x<0 with mu Im q > 0.
- (8) A_R and B_R. Interlacing alpha < z0 < p < 1 < beta. Both residues are positive.

## B. Dependencies

| dependency | kind | used for | status |
|---|---|---|---|
| Prop. 4.1 / Lemma 4.2, generators | core | F_T | EXACT (library) |
| (4.2) | core | minor cone gives Phi >= 0 | EXACT (library) |
| Appendix B (Thm B.1, Lemma B.2) | common lemma | collision layer | CITED |
| "three-A/three-C common-circle theorem" | supplement, unstated | real v | replaced by the auditor's argument (S9). The face {sigma(4)=4} together with T has generators [0,1,2,5,6], all minor (EXACT) |
| two-orthogonal reduction | core | U(4) | CITED |

No other contact case is used.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| S1 | F_T = {0,1,2,3,4,5,6,27,28,33,34}; nonminor {28,34} | identity | EXACT | |
| S2 | span(minor) = span(all) = im P intersected with {g_T = 0}, dimension 8 | identity | EXACT | |
| S3 | the minor-cone dual is exactly (1), with 10 rays | finite | EXACT | own facet enumeration. INFO: A = -1, B = 1 on f28; A = 1, B = -1 on f34 |
| S4 | collisions and kappa = 0 | citation | CITED | Appendix B |
| S5 | F0 F7 = F1 F6 implies H = 0 and all contact ratios are real | algebra | PROOF-READ | as for T_18 |
| S6 | F6/F0 = x, F8/F6 = 1/(1-k), F5/F0 = (1-w)A3/(A3-w), A3 = hw/(h-1+w) | identity | EXACT | |
| S7 | F1/F0 = l is equivalent to curve (2) (up to the factor -1) | identity | EXACT | |
| S8 | with w = 1 + q(v-1), the curve equals (v-1)^2 times a polynomial linear in v; its root is v(q) of (3); the v-coefficient is -qT(q); the numerator at q = k/c is kl(h-1)(k-1)^2(l-1)/c^2, nonzero | identity | EXACT | |
| S8b | the chart is exhaustive: v != 1 by distinctness, so q is defined; qT(q) != 0 on original configurations; q real implies v real | analytic | PROOF-READ | |
| S9 | real v: A0,A1,A2 and C0,C1,C2 are real. The face {sigma(4)=4} = {0,2,6,8,12,14} meets T in {0,6,8} with F0 != F6, so Re(kappa(A3 - C3)) = 0 and the whole face consists of contacts | analytic + finite | PROOF-READ (argument supplied by auditor) + EXACT | generator list [0,1,2,5,6], all minor |
| S10 | contact ratios; R14, R16 = h R14, R19 | identity | EXACT | |
| S11 | g14 != 0 and h > 0 | algebra | PROOF-READ | Im R14 = -(x-1) Im v / ((k-1)|v-1|^2) != 0, and g16 = h g14 |
| S12 | the four sign kernels (4) | identity | EXACT | |
| S13 | the preliminary table is exactly the feasible set | case split | EXACT | enumeration over cells, using h > 0. The six feasible cells are the document's four rows; row 1 covers two x-cells |
| S14 | Pick identity (5) and its positivity | identity + inequality | EXACT + PROOF-READ | Positivity: h > 0 and -l(h-1)(l-1) > 0 come from g4/g3 > 0. No division by c is needed |
| S15 | g14/g3 > 0 gives (x-1)/k > 0, leaving 3 cells | case split | EXACT | |
| S16 | (6) and the g19 formula; in the cells {k>1, x>1} and {k<0, 0<x<1} we get c > 1, bracket > 1 and g19 < 0 | identity + sign | EXACT | the signs of c - 1 = l(k-1) and of -kl(l-1) are fixed by the cell |
| S17 | region (7), with c < 0 | case split | EXACT + PROOF-READ | c = 1 + l(k-1) < k < 0 |
| S18 | identities (8) for A_R and B_R, and v - 1 = J/(qT) | identity | EXACT | |
| S19 | J(0), J(p), J(1), the slope of S, S(0), S(p), J(z0) | identity | EXACT | |
| S20 | in (7): 0 < p < 1, J has roots 0 < alpha < p < 1 < beta, alpha < z0 < p, J'(alpha) > 0 > J'(beta) | inequality | PROOF-READ | p < 1 iff k > c iff l > 1. The leading coefficient h - 1 is negative |
| S21 | the residues of A_R and B_R at alpha and beta are positive. Each function is a real constant plus two simple-pole terms, so A, B > 0 when mu Im q > 0 | inequality | PROOF-READ | Each residue is (numerator)/J' with the signs listed. S(alpha) != 0 and T(alpha) != 0, so there is no cancellation |
| S22 | A, B > 0 puts g in the minor cone, so Phi >= 0 | conclusion | EXACT (S2, S3) + (4.2) | |

## D. Falsification

4000 random chart points were drawn (seed 20261125); half were biased toward region (7). Each was tried with both signs of mu. 2013 were supported (all gaps >= -1e-9 max|g|), and all 2013 lay in region (7). The smallest value of min(A, B)/max|g| was 1.6e-5 > 0. No violation was found.

## E. Verdict

**confirmed-modulo-citations.** Citations: Appendix B (Theorem B.1, Lemma B.2), the two-orthogonal reduction, and (4.2).

Presentation gap (1): p. 1, "so are separately safe" for real v. The cited three-A/three-C common-circle theorem is not stated, and it is not the paper's Theorem C.1. The claim is true, and the auditor supplies the argument (S9).

Presentation gap (2): p. 2. The doc says "interpret divided differences below as Im v/Im q": the ratios in (4) are written with |q|^2 and so on, while g3 carries y = Im q. The identities hold exactly as replayed, and only the reading is ambiguous.

Corollary B.3 is not invoked.
