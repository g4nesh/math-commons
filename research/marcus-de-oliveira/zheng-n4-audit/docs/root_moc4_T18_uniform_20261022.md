# Audit: root_moc4_T18_uniform_20261022 (paper case T_19)

Evidence level: EXACT for every identity, the chart solves, the generator list, the minor-cone dual and the sign-cell analysis. PROOF-READ for exhaustiveness of the chart, the H = 0 argument and the region inequalities. CITED for Appendix B, the two-orthogonal reduction and (4.2). The real-v layer rests on an argument supplied by the auditor (PROOF-READ plus an EXACT generator check). NUMERIC for the falsification run.

Document: `proofs/manin_transfer_20261021/root_moc4_T18_uniform_20261022.pdf` (3 pages). Supplement label T18 = paper T_19 = {0, 1, 4, 6, 7, 8}, literally the same set.

## Deviations and limits

- Conventions follow the T_18 audit (root_moc4_T17_uniform_20261022.md).
- Real v is handled by the "three-A/three-C common-circle theorem", which is not a paper theorem (Theorem C.1 needs a whole spectrum on one circle). The auditor supplies the argument (S9).
- The check output was assembled from two runs of the same script. The exact stages come from one full run, and the numeric stage comes from a rerun after the sampler was biased toward the regions. The code of the exact stages did not change between the runs.

## A. Statements

**Theorem 1 (p. 1).** Let x be in X with P c(x) >= 0 and (P c(x))_pi = 0 for pi in T_19 = [0, 1, 4, 6, 7, 8]. Extra contacts are allowed, and so are all finite spectra, including repeats and a_i + b_j = 0. Then Phi_{c(x)}(O) >= 0 for all O in SO(4). On the eight-distinct layer with v nonreal, c(x) moreover lies in the cone of the ten minor generators. This is at least (*) for T_19, in full.

Intermediate claims:
- (1) The minor generators {0,...,6,21,23,33} span the same 8-dimensional space as all twelve generators, and their cone has exactly 14 facets. Twelve are coordinates g_j >= 0 for j in {2,5,9,10,14,15,16,17,18,19,22,23}. The other two are A = g13 - g2 and B = g2 - g3 + g9 - g12 + g13 + g18.
- Chart (2): A = (x, 1, a2, A3), C = (inf, 0, v, w), with a2 = 1 + k(v-1), w = (1-h)a2/(a2-h), and A3 as displayed. The contact ratios are (1, l, hl, x, lx, -x/(k-1)).
- (3) Six sign kernels. (4) Region I is k>1, 0<h,l<1, x>1; region II is h,k<0, l>1, 0<x<1. (5) Two identities, each a positive sum.

## B. Dependencies

These are the same as for T_18. Generators and dual: EXACT (library). (4.2): EXACT (library). Appendix B (Theorem B.1, Lemma B.2) covers the collision layer: CITED. The common-circle step is replaced by the auditor's argument; the face {sigma(4)=4} together with T has generators [0,1,2,5,6,21,23], all minor (EXACT). The two-orthogonal reduction is CITED. No contact case is cited.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| S1 | F_T = {0,...,6,21,22,23,33,34}; nonminor {22,34} | identity | EXACT | |
| S2 | span(minor) = span(all) = im P intersected with {g_T = 0}, dimension 8 | identity | EXACT | |
| S3 | the dual is exactly (1), with 14 rays | finite | EXACT | own facet enumeration. INFO: A = 1, B = -1 on f22; A = -1, B = 1 on f34 |
| S4 | collisions and kappa = 0 | citation | CITED | Appendix B |
| S5 | F0 F7 = F1 F6 forces H = 0 | algebra | PROOF-READ | |
| S6 | F6/F0 = x; F8/F6 = 1/(1-k); F4/F1 = (1-C3)a2/(a2-C3); unique solves for w and A3 | identity | EXACT | |
| S7 | the chart is exhaustive and the exclusions are collisions | analytic | PROOF-READ | |
| S8 | contact ratios (1, l, hl, x, lx, -x/(k-1)) | identity | EXACT | |
| S9 | real v: the face {sigma(4)=4} meets T in {0,6,8}, so it joins the contact set and gives a minor-only cone | analytic + finite | PROOF-READ (argument supplied by auditor) + EXACT | |
| S10 | the six sign kernels (3); D and a2 have only real zeros | identity | EXACT | |
| S11 | the cell analysis gives exactly regions I and II (4) | case split | EXACT | 81 cells |
| S12 | e > 1, P > e, L > 0 in both regions | inequality | PROOF-READ | h(k-1) > 0 and k(x-1) > 0. For L: l > 0 and k(1-l) > 0 |
| S13 | the identities (5) | identity | EXACT | |
| S14 | every summand of (5) is positive | inequality | PROOF-READ | The signs of h(h-1)(l-1), (x-1)/k, -hk(k-1)(l-1) and -(h-1)(k-1)(x-1) were checked in each region |
| S15 | A, B > 0 puts g in the minor cone, so Phi >= 0 | conclusion | EXACT + (4.2) | |

## D. Falsification

4000 chart samples were drawn (seed 20261125); half were biased into regions I and II, and both signs of mu were tried. None of them gave a supported gap vector. The largest min-gap/max|g| among the 1982 region points was -5.8e-6. A short Nelder-Mead search that maximized the smallest gap reached 0 only at a collision boundary (x -> 1, h -> 0, l -> 1). This is NUMERIC evidence that the supported nonreal eight-distinct layer of T_19 may be empty. The document's argument does not need this. The intermediate claim (5) was therefore tested directly on the 1982 region points, with mu chosen so that g2 > 0: min(A, B)/max|g| = 6.2e-5 > 0 and no violation occurred. The main-statement falsification is vacuous here.

## E. Verdict

**confirmed-modulo-citations.** Citations: Appendix B (Theorem B.1, Lemma B.2), the two-orthogonal reduction, and (4.2).

Presentation gap (1): p. 1, real v. The common-circle theorem is cited without a statement. The auditor supplies the argument (S9).

Presentation gap (2): p. 1 and p. 2. Two garbled sentences, "Thus A, B >= 0 suffices for a nonnegative sum of minor responsesty of a supported source" and "proves a nonnegative sum of minor responsesty of this entire layer". The intended meaning is clear.

Corollary B.3 is not invoked.
