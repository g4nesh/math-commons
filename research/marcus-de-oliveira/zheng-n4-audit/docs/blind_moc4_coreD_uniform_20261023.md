# Audit: blind_moc4_coreD_uniform_20261023 (seven-contact set {0,1,2,7,9,16,23}, whole cone)

Evidence level: identities, chart inverses, the dominant and zero-main certificates, the vertex enumeration and the triangulation facet pairing EXACT (own code); case assembly PROOF-READ; the low-pair exclusion relied on through its own audit (blind_moc4_coreD_low_pair_20261023: confirmed-modulo-citations); atoms and the Hessian theorems CITED.

PDF: `proofs/manin_transfer_20261021/blind_moc4_coreD_uniform_20261023.pdf` (4 pages). Replay: `replay/blind_moc4_coreD_uniform_20261023__check.py`, output `replay/blind_moc4_coreD_uniform_20261023__check_output.txt`.

## Deviations and limits

- G (385 terms) is the one validated in the low-pair replay by the exact Hessian identity (3). The uniform script recomputes that identity.
- The balanced polytope is the region a + c > 1, b <= 1. I enumerated its vertices myself (27, which equals the certificate's list). For the 81 simplices I took only the vertex index lists from the JSON, and I checked them myself: each is nondegenerate, and the facet pairing holds (each boundary facet lies in one simplex; each interior facet lies in exactly two simplices, on opposite sides). The coverage then follows from a degree argument (PROOF-READ). I recomputed all 81 positivity polynomials myself (about 52k terms each), in a second run with `--polytope-tables`. That run is appended to the output file.
- Capacity atoms (sharp SOC on {2,12,15}; old-H (3,15), (5,2), (11,12)) are CITED and tested numerically over SO(4).

## A. Statements

**Theorem 1 (p. 1).** Let T = {0,1,2,7,9,16,23} (Table 5 no. 12). For every x in X with P c(x) >= 0 vanishing on T (non-strict), Phi_{c(x)}(O) >= 0 for all O in SO(4). This covers all spectra, through the closed conditions on Re X.

## B. Dependencies

| dependency | kind | status |
|---|---|---|
| canonical cone, I7 rows [15,14,10,17,8,4,3]; source S3 | own | EXACT |
| (3): det H = -4 G^2/Q^2 | identity | EXACT (own Hessian) |
| det H(Re x) = 0 and principal 12-minors >= 0 on X | common lemma | CITED (randomized exact support) |
| blind_moc4_coreD_low_pair_20261023 (region (4)) | other-doc | audited in this batch: confirmed-modulo-citations; its statement is exactly the region (4) |
| sharp SOC atoms (f_i, f_j; f26), {i,j} in {2,12,15} | other doc | CITED; NUMERIC no violation |
| old-H atoms t(c-d) >= d^2, u(a-d) >= d^2, v(b-d) >= d^2 | other doc | CITED; NUMERIC no violation |

## C. Steps

| # | claim | kind | status |
|---|---|---|---|
| 1 | (1)-(2): cone, C, Q; ordering a <= c <= b by S3 | identity | EXACT |
| 2 | (5), (6): H3 bound; failure gives a(b+c), b(a+c), c(a+b) < 1 | inequality | PROOF-READ |
| 3 | (7)-(9): cap(X,m) solves lambda^2 = X(m - lambda); the three old-H pairs use disjoint coefficients, so lambda_a + lambda_c + lambda_b >= 1 gives safety | identity + logic | EXACT + PROOF-READ |
| 4 | assembly: d > 0, ordered positive mains fall in exactly one of (4) [b <= 1, a+c <= 1], dominant [b > 1], balanced [b <= 1, a+c > 1] | case split | PROOF-READ (b > 1 forces b dominant, because the balanced case would give H3 >= b > 1; a + c > 1 >= b forces the balanced case with s < 2) |
| 5 | (10)-(12): the dominant chart is onto {b > 1, a, c > 0, b(a+c) < 1, lambda-sum < 1}, with D = 1/theta | identity | EXACT (sympy, full inverse) |
| 6 | dominant table: X^2 Y^2 VD^2 G(list) has nonnegative coefficients, with a strict monomial | certificate | EXACT (see output for the term count) |
| 7 | (13)-(15): the balanced polytope, 27 vertices, 81-simplex triangulation | finite geometry | EXACT (own enumeration; facet pairing) + PROOF-READ (degree argument) |
| 8 | 81 simplex polynomials nonnegative and strict under (14) | certificate | see E |
| 9 | two zero mains: G at a = c = 0 has positive coefficients and contains d^9 | certificate | EXACT |
| 10 | one zero main, b <= 1, c < 1: the low-pair a+c<1 tables have nonzero constant terms, so they stay strict at A = 0; the low-pair root and rank argument goes through unchanged (its density step is about the algebraic root graph, which is analytic across a = 0) | extension | EXACT (constant terms) + PROOF-READ |
| 11 | a = 0, c = b = 1: H3(0,1,1) = 1 = d, so the cost is safe | inequality | PROOF-READ |
| 12 | one zero main, b > 1: chart (16) with u = U free; Y^2 VD^2 G has nonnegative coefficients and a strict monomial | certificate | EXACT (see output) |

## D. Falsification

Capacity atoms over SO(4) (40 Nelder-Mead starts each): no violation. Realizable stratum points (seeds 2026-2028): G is about 0, as expected, and min Phi >= about 0.

Results of the full run. The output file holds two runs, both exit 0 with 0 EXACT failures. The second run is the polytope tables, about 330 s.

- (3) holds exactly with my own Hessian.
- Dominant table: 180265 terms, exactly the document's count, all coefficients nonnegative. The coefficient of X^2Y^2 is 2, as the document states.
- Zero-main chart (16): 21497 terms, again equal to the document's count.
- G at a = c = 0: 63 terms, all positive, with a d^9 term.
- Vertex enumeration gives 27 vertices, the same set as the certificate's.
- The triangulation facet pairing holds: 259 boundary facets and 154 interior facets (81 x 7 = 259 + 2 x 154).
- All 81 simplex polynomials have nonnegative coefficients (52-54k terms each). The support test passes on all 127 subsets of every simplex.
- The constant terms of the low-pair a+c<1 tables are nonzero (+-1, +-2).
- Capacity atoms (NUMERIC): minima are about -1e-29 to 0, with no violation.
- Three realizable points with contact set exactly T: G is about 0 there, as it should be.

## E. Verdict

**confirmed-modulo-citations.** I found no gap. The citations are:

- the low-pair exclusion lemma, audited in this batch as confirmed-modulo-citations;
- the sharp-SOC and old-H atoms;
- the global Hessian determinant and inertia theorems for Re X;
- (4.2).

The extension of the low-pair argument to a = 0 is written tersely in the document ("remain valid verbatim"). It is correct: the a+c<1 tables have nonzero constant terms, so they stay strict at A = 0, and the density argument concerns the algebraic root graph, which is analytic across a = 0. I record this as PROOF-READ with the argument supplied by the auditor. It is not a gap.
