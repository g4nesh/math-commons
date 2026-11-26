Evidence level: EXACT for G9, the scalar capacity algebra, chart coverage and all twelve sign tables (eleven all-negative, and one by an auditor-built AM-GM certificate); PROOF-READ for the case logic; CITED for the zero-coefficient transfers, the a = 0 face, the core and four-coefficient capacities with finite Gram closure, and radial necessity.

# final_referee_T46_positive_a_large_product_20261023 ("The large-product sector of the T46 contact cone", 3 pp.)

Replay: `replay/final_referee_T46_positive_a_large_product_20261023__check.py`. It was run in three chunks (`scalar B`, `A`, `C`), and the output is concatenated in `..._check_output.txt`. All EXACT checks pass.

## Deviations and limits

- G9 is rebuilt from the paper's definitions (see the negative_j_Heron report). The document's square list and verifier were not used.
- My clearing of denominators is not the document's minimal clearing: I multiply by D_x^{deg_x G} for each variable and then divide out the per-variable denominators. My polynomials are therefore the document's numerators times extra factors with positive coefficients, and my term counts are larger: for example 265584 against the document's 196914. The sign conclusions are the same, since every factor involved has positive coefficients on the positive chart.
- In the one mixed-sign table my polynomial has 86 positive coefficients; the document's has 80. I built my own certificate: 86 AM-GM pairs X^{N1} - 2X^M + X^{N2} with M = (N1+N2)/2, a nonnegative remainder, and strictly positive remainder monomials free of the optional variables (and with J exponent 0). It was verified in exact rational arithmetic. This is an auditor-supplied certificate in place of the document's 82 binomial squares.
- c and a (the variables C, H) were treated as possibly zero in the strictness checks.

## A. Statements

Contact set {1,3,5,6,7,8} = paper T_47. The cost is q f27 + p f1 + j f59 + c f5 + e f2 + d f28 + b f6 + a f4 with all coefficients >= 0.

**Theorem 1.** Every such cost with d > 0 (normalized d = 1), qb > d^2 and q != b that equals c(x) for some x in X satisfies Phi_{c(x)}(O) >= 0 for all O in SO(4). This holds for arbitrary a >= 0.

The proof shows that the coefficient set outside two sufficient capacities contains no realizable cost. Zeros among q, p, j, c, e, b, and the face a = 0, are delegated to earlier results. Degenerate spectra are covered because the only spectral input is radial necessity. The line q = b and the wall qb = 1 are excluded, and are handled elsewhere: q = b in the small_or_balanced note, qb = 1 in the T46 mixed note, which is not in this batch.

## B. Dependencies

| dependency | status |
|---|---|
| zero-coefficient transfers for q, p, j, c, e, b = 0 and the closed a = 0 face ("previously established") | CITED (not resolved to specific documents here) |
| beta = -1 coherent core: p >= (t-1)^2/(j - lambda t^2), j > K t^2, with finite closure | the core document of this batch (Thm 1 at beta = -1) and foundations/so4-atoms-part2b G3 ("the beta = -1 forms (lambda, K, K - lambda, t = j/lambda, balanced case): EXACT") |
| four-coefficient capacity ab(e-1) >= a+b+e+1, i.e. four(4,6;2;28) | foundations/so4-atoms.md (proved family); orbit membership (0,6;10;54) -> (4,6;2;28) re-checked EXACT with moc4lib |
| radial necessity G9 = 0 | audit item (5) + the identity det Hess = -4G^2/Q^2 (EXACT here) |

## C. Steps

| # | claim | status |
|---|---|---|
| 1 | C0 < 0, so the chart Q = 0 contains no realizable point | EXACT (all coefficients negative) |
| 2 | K - lambda = (q-b)^2/((qb-1)(q+b-2)) > 0 | EXACT |
| 3 | the core at beta = -1 gives threshold (t-1)^2/(j - lambda t^2) on j > K t^2 | EXACT |
| 4 | the optimum t = j/lambda gives 1/j - 1/lambda when j < lambda^2/K; the endpoint t = sqrt(j/K) gives (1 - sqrt(j/K))^2/(j(1 - lambda/K)); j >= K is safe | EXACT (stationarity and values); the endpoint by finite closure (CITED) |
| 5 | the charts (big, small) = (r+A, r^2/(r+A)) with r = 1+V cover qb > 1, q != b, in both orientations | PROOF-READ |
| 6 | the low and middle j-charts cover (0, lambda^2/K) and [lambda^2/K, K); the cap formulas are right | EXACT |
| 7 | the e/a charts cover the complement of the four-coefficient capacity (e0 = 1 + 1/b; a4 > 0 iff e > e0) | EXACT/PROOF-READ |
| 8 | 11 tables have every coefficient negative, positive denominators, a strictly negative monomial free of c and a, and one with J exponent 0 in the middle range | EXACT (term counts differ, see the limits) |
| 9 | the last table (b > q, middle, e > e0) has mixed signs and a certificate is needed | EXACT, by the auditor's AM-GM certificate |
| 10 | G9 < 0 on every complementary chart, contradicting radial necessity | PROOF-READ |

## D. Falsification

NUMERIC. In 77 stratum samples (see the negative_j_Heron report), the one point with d > 0 lies in this sector, and it satisfies both the core and the four-coefficient capacity, as the theorem requires. No Leibniz minimization of the main statement was run; it is a special case of the conjecture.

## E. Verdict

**confirmed-modulo-citations.** The citations are the zero-coefficient transfers and the a = 0 face, the core and four-coefficient capacities with finite Gram closure (foundations), and radial necessity. Every chart table replays. The mixed-sign table needed the auditor's own certificate, because the document's square list was not used.
