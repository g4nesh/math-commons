Evidence level: EXACT for every identity, certificate and coefficient-sign table in the document (own code, `replay/final_referee_inventory87_uniform_20261023__check.py`). PROOF-READ for the SOC resource allocation, the chart exhaustion and the radial-necessity argument. CITED for the SOC atom (proved in `foundations/so4-atoms.md`) and Theorem 1.2. NUMERIC for the falsification runs.

# Audit: final_referee_inventory87_uniform_20261023 (primary proof of paper case T_88)

PDF: `proofs/manin_transfer_20261021/final_referee_inventory87_uniform_20261023.pdf`, 4 pages. Auditor label C44-T87-T89.

## Deviations and limits

- None from DOC_AUDIT. I did not open the certificate `final_referee_inventory87_uniform_exact_20261023.py` or its tables. G was recomputed from scratch: the Hessian of the paper's C, an exact Bareiss determinant, division by -4Q^2 and an exact square root. The ten coefficient tables were rebuilt from my own G.
- The Hessian identity (8) was computed at d = 1, and G was then re-homogenized. This is lossless, because every polynomial involved is homogeneous: the Hessian entries are linear forms, N = Q^4 det has degree 22, and Q^2 G^2 has degree 22. I used Bareiss with full pivoting on the sparsest entry. It makes the same number of exact divisions (819) as the document's.
- The document's cubic is the paper's C, with sign +1 in (7). Its profile numbering is the paper's Appendix D numbering (B1).
- The chart inverses were checked exactly on 500 rational samples. The algebra of the inverses was also proof-read.

## A. Statements

**Theorem 1 (p.1).** For every x ∈ X with P c(x) ≥ 0 and (P c(x))_σ = 0 on T = {1, 5, 6, 10, 13, 14}, Φ_{c(x)}(O) ≥ 0 for all O ∈ O(4). Contacts are non-strict, and all spectra are included. The proof never uses a spectral lift, only x ∈ X and radial necessity, so it covers repeated spectra and cross-collisions automatically. T is literally the paper's T_88. This is (*) for T_88. The final sentence, which carries the result to unitary U, is Theorem 1.2 (foundations); it is not needed for (*).

Intermediate claims stronger than the conjecture:

- (I1) SOC capacity: the source is safe whenever d ≤ H3(a, b, c).
- (I2) Every realizable source with d > 0 has height H = -C/Q and det Hess C(u) = 0.
- (I3) G > 0 on {H3 < d} ∩ (nonnegative orthant).

Together, (I2) and (I3) say that realizable sources never have H3 < d.

## B. Dependencies

| dependency | kind | resolution | status |
|---|---|---|---|
| SOC atom: x f_i + y f_j + δ f30 ≥ 0 on SO(4) if xy ≥ δ², {i, j} ⊂ {1, 10, 59} | common lemma (SO(4) atom) | family SOC, representative (1, 9; 22), orbit 144 (pair unordered). The tuples (1,10;30), (1,59;30), (10,59;30) all lie in the orbit, with f30 as the compensated term (EXACT, B9 and B10) | established (`foundations/so4-atoms.md`) |
| radial necessity | common lemma | re-proved in the document by polarization; I proof-read that argument. It matches the brief's item (5) | established + PROOF-READ |
| Theorem 1.2 (two-real-response) | core | only used for the unitary restatement | CITED (foundations) |

No other contact case is used, which is consistent with Stage 1 in Appendix G.

## C. Steps

| # | claim (location) | kind | status | how checked |
|---|---|---|---|---|
| 1 | P[T,:] has rank 6; the eight profiles (1) are a kernel basis | identity | EXACT | B3 |
| 2 | (2): q, a, v, b, w, t, c = g20, g16, g11, g7, g12, g3, g0 and d = g19 - g20 | identity | EXACT | B4, B5 |
| 3 | generators = {0,1,6,10,27,30,34,36,59} (the nine "crosscheck" profiles) | finite | EXACT | B1, B2 |
| 4 | (3) f30 + f34 = f0 + f6 + f27 + f36; g19, g15, g4, g8 = q+d, v+d, w+d, t+d; rewrite for d < 0 | identity | EXACT | B6, B7 |
| 5 | the transpose fixes T, swaps f30 and f34, and permutes the seven minors | finite | EXACT | B8 |
| 6 | (4) SOC atoms for the pairs in {1, 10, 59} with f30 | citation | EXACT (orbit) + established | B9, B10 |
| 7 | (5) H3 allocation: balanced case; dominant case; scaling for smaller d; S = 0 | inequality (allocation) | PROOF-READ + EXACT samples | C1. Atom resources sum to a, b, c exactly. The capacities are (a+b-c)/2 and so on in the balanced case, and b·sqrt(a/S) and c·sqrt(a/S) in the dominant case, with sum sqrt(aS) |
| 8 | (7): C(L + Hκ) = C + HQ with the displayed C and Q; C ≥ d³ > 0, Q ≥ d² > 0, H = -C/Q | identity | EXACT | C2 (sign +1) |
| 9 | radial necessity: u = Re x; C(u) = 0 and det Hess C(u) = 0 (including the case x = y) | analytic | PROOF-READ | polarization identities; Euler's relation when Im x = 0 |
| 10 | (8): Q^4 det Hess C(L - Cκ/Q) = -4 Q² G², with G of degree 9, 2285 terms, d^9 coefficient +1 | identity | EXACT | C3 to C10: 819 divisions; det degree 3 in H; exact division by -4Q²; exact square root; 2285 terms; multidegree (3,2,3,2,3,3,2,9); G has 795 negative coefficients in the original variables |
| 11 | balanced chart (10): cleared G has 20241 coefficients, all positive, constant 1 | coefficient signs | EXACT | D1 (D^6, 20241 terms, matching the count in the document) |
| 12 | six dominant charts (12)-(13): 27372 positive coefficients each, with the strict monomial 8X | coefficient signs | EXACT | D2 (27372 terms each; monomial 8X present) |
| 13 | three solo charts: 450 positive coefficients each, constant 1 | coefficient signs | EXACT | D3 |
| 14 | the charts and inverses cover {H3 < 1} ∩ orthant exactly, including all zero faces | case split | PROOF-READ + EXACT samples | D4. Balanced: σ = x+y+z < 1, X = x/(1-σ). Dominant: X = 1/S - 1 > 0, Y ≥ 0 iff M ≥ S, Z = L1/L2 or L2/L1. Solo: S = 0 |
| 15 | conclusion: G > 0 on H3 < d contradicts G = 0; d ≤ H3 is safe by (5); d = 0 is a sum of minors; d < 0 reduces to d > 0 by step 5 | logic | PROOF-READ | |

## D. Falsification

- NUMERIC E1. G evaluated at 3000 random nonnegative sources with H3 < d = 1 (seed 11), with about 20% of the coordinates set to 0. Minimum 0.305 > 0.
- NUMERIC D8 (uninformative). The least-squares stratum search (seeds 0-5) found 6 points of the T_88 stratum, but all had d ≈ 0. So the test of "d ≤ H3 at realizable sources" is vacuous here. A targeted search for realizable points with d > 0 was not done.
- NUMERIC D9. At the same 6 points, the minimum of Φ over SO(4) (8-start BFGS) is 4.5e-15 ≥ 0. This is a brief test only.

## E. Verdict

**confirmed-modulo-citations.** The citations are the SOC atom (established, with the tuples verified in the proved orbit and the roles as printed) and Theorem 1.2 (only for the unitary rephrasing). Every identity was reproduced exactly from the paper's definitions, with the document's term counts matched exactly. These are (1) to (3), the d < 0 rewrite, the transpose, (7), (8) with G (2285 terms), and all ten positivity tables (20241, 6 × 27372 and 3 × 450 coefficients). The argument has no gap.

Presentation remark: the display (1) is labelled twice on p.2. Both "(1)" and "(7)" are used for the cubic identity.
