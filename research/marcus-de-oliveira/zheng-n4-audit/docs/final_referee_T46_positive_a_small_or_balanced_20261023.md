Evidence level: EXACT for G9, all six sign tables (term counts match the document exactly), the nine wall evaluations and the chart identities; PROOF-READ for the case logic; CITED for the zero-coefficient transfers, the a = 0 face, the four-coefficient and beta = -1 core capacities with finite closure, and radial necessity.

# final_referee_T46_positive_a_small_or_balanced_20261023 ("Two coefficient regions of the T46 support inequality", 2 pp.)

Replay: `replay/final_referee_T46_positive_a_small_or_balanced_20261023__check.py` (output `..._check_output.txt`: 20 EXACT checks pass).

## Deviations and limits

- G9 is rebuilt from the paper's definitions (see the negative_j_Heron report). The document's verifier was not used.
- The six tables reproduce the document's term counts exactly: 4916, 1815, 10809 and 4095, 1377, 8398.
- For the walls q = 1 and/or b = 1 I ran my own nine evaluations. Each is all-negative with a strict monomial. I did not compare them with the document's "leading coefficient" remark.

## A. Statements

Contact set {1,3,5,6,7,8} = paper T_47. The cost has nonnegative canonical coefficients (q,p,j,c,e,d,b,a).

**Theorem 1.** If the cost equals c(x) for some x in X and (q <= d and b <= d) or q = b, then Phi_{c(x)}(O) >= 0 for all O in SO(4). This holds for every a >= 0.

The case d = 0 is a sum of minors. Zeros of q, p, j, c, e, b, and the face a = 0, are delegated to earlier results. Degenerate spectra are covered through radial necessity.

## B. Dependencies

| dependency | status |
|---|---|
| zero-coefficient transfers and the a = 0 face | CITED (not resolved here) |
| four-coefficient capacity ab(e-1) >= a+b+e+1 = four(4,6;2;28) | foundations/so4-atoms.md; orbit membership re-checked EXACT |
| beta = -1 core capacity p >= max(1/j - (x-1)/(x+1), 0) at q = b = x > 1 | core document of this batch / foundations part2b G3 (balanced case, EXACT); at q = b, K = lambda = (x+1)/(x-1) (EXACT) |
| radial necessity G9 = 0 | audit item (5) + det Hess = -4G^2/Q^2 (EXACT) |

## C. Steps

| # | claim | status |
|---|---|---|
| 1 | C0 < 0 (d, b > 0), so Q != 0 at realizable points | EXACT |
| 2 | the complement of the four-coefficient capacity is {e <= e0} union {e > e0, a < a4}: for e <= e0, b(e-1) <= 1 gives ab(e-1) <= a < a+b+e+1 | PROOF-READ; a4 > 0 exactly when e > e0 (EXACT) |
| 3 | sector q, b < 1: three tables, all coefficients negative, strict monomial free of A and C; p, j, c are unrestricted | EXACT (4916 / 1815 / 10809 terms) |
| 4 | walls q = 1 and/or b = 1: nine direct evaluations, all negative and strict | EXACT |
| 5 | balanced q = b = 1 + V: the core capacity, its complement chart j = lambda J/(1+J), p = VP/((V+2)J(1+P)), and e0 = (V+2)/(V+1) | EXACT |
| 6 | the three balanced tables are all negative and strict | EXACT (4095 / 1377 / 8398 terms) |
| 7 | the equality walls are covered by the capacities with finite Gram closure | CITED |

## D. Falsification

NUMERIC. No stratum sample with d > 0 fell in this sector (see the negative_j_Heron report). In this sector the claim is that every realizable cost satisfies the four-coefficient capacity, which is a strong claim. The exact tables prove it, and I found no counter-evidence. No Leibniz minimization was run.

## E. Verdict

**confirmed-modulo-citations.** The citations are the zero-coefficient transfers and the a = 0 face, the four-coefficient and core capacities with closure, and radial necessity. All tables replay exactly, with matching term counts.
