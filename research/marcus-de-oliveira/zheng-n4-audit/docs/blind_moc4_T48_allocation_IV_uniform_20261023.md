Evidence level: EXACT for G, the allocation algebra, the 111608-term certificate with its 45 square groups, the closed-wall homogenization in (A, WA) and (B, WB), and the zero-face embeddings; PROOF-READ for the logic; CITED for the six zero-face contact cases; NUMERIC for falsification only.

# blind_moc4_T48_allocation_IV_uniform_20261023 (3 pp.): double old-H allocation region of T_49

Auditor label C26-T49. Script: `replay/blind_moc4_T48_allocation_IV_uniform_20261023__check.py` (0 EXACT failures).

## Deviations and limits

- G is rebuilt from the paper's definitions. See the allocation-II report for the construction. It is hash-checked here (sha256 c89c8cbc...0a49), and my G = -(supplement's G).
- The scope is the eight-role subcone of the ten-generator cone N_{T_49}: f1 and f9 are absorbed elsewhere. The same remark applies as in the allocation-II report.
- I did not check that the T_49 cells cover the eight-coefficient face. That is the assembly's job.

## A. Statement

Theorem 1. Take an eight-role source with nonnegative coefficients, normalized to d = 1, with S = e + j > 0, a <= 2e/S and Q = q + b <= 2j/S, coming from any x in X (collisions included). Then Phi_{c(x)}(O) >= 0 for all O in SO(4). Both equality walls, their intersection and c = 0 are included. The faces q, p, j, e, b, a = 0 are delegated, and d = 0 is a sum of minors.

## B. Dependencies

| dependency | check |
|---|---|
| old-H (aux f2, main f5; f40) and (aux f29, main f5; f40) | EXACT: both are in the proved 288-member orbit with these roles (the reversed roles are not). The coefficient form e(p - h) >= h^2 matches t(c - d) >= d^2 |
| chains (f5; f0, f2; f40), (f5; f7, f29; f40), (f5; f29, f41; f40) | EXACT: all are in the proved 576-member orbit |
| radial necessity + det Hess = -4G^2/Q0^2 + the 42-term identity | established + EXACT |
| zero faces -> paper T_39, T_1, T_3, T_24, T_43, T_44 | embeddings EXACT; proofs CITED |

## C. Steps

| # | claim | status | how |
|---|---|---|---|
| 1 | currents he = e/S - a/2, hj = j/S - Q/2, alpha = a/2, beta = Q/2 are >= 0 on the region and sum to 1; the e and j allocations S he + S alpha = e and S hj + S q/2 + S b/2 = j are exact | EXACT | sympy |
| 2 | requirement sum = p* = 1 + 1/S - (a+Q)/4 > 0 | EXACT | sympy (each old-H needs h + h^2/(S h); each chain needs x^2/leaf + x^2/(S x)) |
| 3 | realizable implies G = 0 | EXACT + established | as in allocation II |
| 4 | the chart e = SL/(1+L), j = S/(1+L), a = 2LA/((1+L)(1+A)), q + b = 2B/((1+L)(1+B)), q/b = C1, p = P/(1+P) p* is onto the strict region | PROOF-READ | each fraction inverts uniquely; L = e/j, S = e + j |
| 5 | D G, with D = 4S^2(1+L)^9(1+A)^4(1+B)^6(1+C1)^6(1+P)^2, has 111608 terms; its 135 coefficients of the wrong sign lie in 45 C1-groups, each (C1-1)^2 R(C1) with R >= 0; marker 108 A^4 B^6 S P with a nonnegative group | EXACT | flint |
| 6 | walls: homogenize to degree 4 in (A, WA) and degree 6 in (B, WB); the pattern holds for WA, WB >= 0; the marker survives WA = 0, WB = 0, both, and c = 0 | EXACT | direct substitution with A + WA and B + WB, then the four specializations |
| 7 | p = p*: capacity; zeros: delegated | PROOF-READ / CITED | |

## D. Falsification

NUMERIC: 12 random region points with p = p* (seed 12), min Phi over SO(4) by BFGS (8 restarts): worst -4.1e-16 relative, so the capacity is sharp and no violation was found.

## E. Verdict

**confirmed-modulo-citations.** Citations: the six zero-face cases (embeddings verified). The atoms and radial necessity are established by this audit. No gap.
