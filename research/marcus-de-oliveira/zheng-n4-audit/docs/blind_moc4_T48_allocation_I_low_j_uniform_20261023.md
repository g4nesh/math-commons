Evidence level: EXACT for G, the both-chain capacity algebra, every chart identity and all three sign certificates (245193, 166518 and 457338 terms, every coefficient of one sign); PROOF-READ for chart exhaustiveness and the B0 = 0 hand-off; CITED for the zero-face cases and for allocation II (audited in this batch); NUMERIC for falsification only.

# blind_moc4_T48_allocation_I_low_j_uniform_20261023 (2 pp.): large-lambda, low-j sector of the both-chain cell

Auditor label C26-T49. Script: `replay/blind_moc4_T48_allocation_I_low_j_uniform_20261023__check.py` (0 EXACT failures).

## Deviations and limits

- G is rebuilt from the definitions; my G = -(supplement's G). The scope is the eight-role subcone, as in the allocation-II report.
- My cleared polynomials are G(chart) times prod den^deg, with the factors (2+X)^2, (2+X)^4 and (XB+X+B)^2 removed for charts 1, 2 and 3. Their term counts equal the document's exactly (245193, 166518, 457338), and so do the markers (constant 10616832, 131072 X^3, 1296 X^13, with the opposite sign).

## A. Statement

Theorem 1. Normalize d = 1 and take a, e, j, q, b > 0, c >= 0, Q = q + b, Ea = ae/(a+e), JQ = jQ/(j+Q), lambda = 1/(Ea + JQ), u = lambda Ea, v = lambda JQ. Suppose the source is realizable (any x in X) and satisfies 2u <= a, 2v <= Q, lambda >= 1/2 and j <= 1. Then Phi >= 0 on SO(4). Finite walls and c = 0 are included; zeros of q, p, j, e, b, a are delegated.

## B. Dependencies

| dependency | check |
|---|---|
| chains (f5; f0, f2; f40), (f5; f7, f29; f40), (f5; f29, f41; f40) | EXACT orbit membership |
| radial necessity, det Hess = -4G^2/Q0^2, the 42-term identity | established + EXACT |
| allocation II (this batch) for B0 = 0, L > 1 | EXACT: u(1/a+1/e) - v/j = v/Q holds identically, so at Q = 2v it equals 1/2. Allocation II's u is the unique solution of that linear equation, so the two u agree, and the point satisfies II's hypotheses 0 < u < 1, 2u <= a, Q = 2v. II is confirmed-modulo-citations in this batch |
| zero faces -> T_39, T_1, T_3, T_24, T_43, T_44 | embeddings EXACT, proofs CITED |

## C. Steps

| # | claim | status | how |
|---|---|---|---|
| 1 | three chains with currents u, vq/Q, vb/Q, j shared in proportion: requirement = lambda, currents sum to 1 | EXACT | sympy |
| 2 | base chart: a = 2u(1+A), Q = 2v(1+B0), e = 2u(1+A)/(A+T+AT), j = 2v(1+B0)/(B0+T+B0T), p = (1+T)/2 * P/(1+P) | EXACT | lambda = (1+T)/2, u = lambda Ea, v = lambda JQ, rechecked in each chart with reduced flint fractions |
| 3 | j <= 1 iff T >= Tj = 2/(1+L) - B0/(1+B0) | EXACT | (1-j)(B0+T+B0T) = (1+B0)(T - Tj), with a positive denominator |
| 4 | charts 1-3 exhaust B0 > 0 (L <= 1; L > 1 with B0(L-1) >= 2; L > 1 with 0 < B0(L-1) < 2), and in each the constraint T >= max(0, Tj) is exactly the chart's range | EXACT + PROOF-READ | EXACT: the formula for Tj in chart 1, Tj = (2 - B0 X)/((2+X)(1+B0)) <= 0 in chart 2, B0(L-1) = 2B/(1+B) in chart 3; the inverses were checked by hand |
| 5 | B0 = 0 with L <= 1 is chart 1 at B = 0; B0 = 0 with L > 1 goes to allocation II | PROOF-READ + EXACT | step B |
| 6 | three sign certificates: every coefficient of one sign, with a strictly positive marker for X > 0 (charts 2 and 3) and a constant marker (chart 1) | EXACT | flint |
| 7 | finite walls (A = 0, Z = 0, B = 0, c = 0, lambda = 1/2) | PROOF-READ | all chart variables >= 0 and all coefficients of one sign; the markers involve only X or the constant; e finite requires A + T + AT > 0, which is what excludes A = Z = 0 in chart 2 |

## D. Falsification

NUMERIC: 12 random cell points with p = lambda (seed 13), min Phi over SO(4) (8 BFGS restarts): worst -1.1e-16 relative, so the capacity is sharp and no violation was found.

## E. Verdict

**confirmed-modulo-citations.** Citations: the zero-face cases, and allocation II (audited in this batch). No gap.
