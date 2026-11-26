Evidence level: EXACT for G, the capacity algebra, every chart identity and the full 374775-term certificate (with the document's denominator); PROOF-READ for exhaustiveness; CITED for zero faces and the large-lambda note; NUMERIC for falsification only.

# blind_moc4_T48_allocation_I_small_lambda_low_j_uniform_20261023 (2 pp.): small-lambda, low-j both-chain region

Auditor label C26-T49. Script: `replay/blind_moc4_T48_allocation_I_small_lambda_low_j_uniform_20261023__check.py` (0 EXACT failures).

## Deviations and limits

- G is rebuilt from the definitions; my G = -(supplement's G). The scope is the eight-role subcone (see the allocation-II report).

## A. Statement

Theorem 1. Normalize d = 1 and take a, e, j, q, b > 0, c >= 0, with lambda, u, v as in allocation I. If 0 < lambda < 1/2 and j <= 1, every realizable source satisfies the support inequality. The combination sentence ("together with the large-lambda note, the whole both-chain cell with j <= 1") is a remark. It is correct: the large-lambda note covers lambda >= 1/2 with 2u <= a and 2v <= Q, and for lambda < 1/2 those two conditions hold automatically (step 3).

## B. Dependencies

Chains (f5; f0, f2; f40), (f5; f7, f29; f40) and (f5; f29, f41; f40): EXACT orbit membership. Radial necessity, the Hessian identity and the 42-term identity: established + EXACT. Zero faces (T_39, T_1, T_3, T_24, T_43, T_44): embeddings EXACT, proofs CITED. Large-lambda note (allocation I low-j, this batch) only for the remark and for lambda = 1/2.

## C. Steps

| # | claim | status |
|---|---|---|
| 1 | three chains need exactly lambda; currents sum to 1 | EXACT |
| 2 | T = 2 lambda/(1 - 2 lambda), A = a/e, B = Q/j, L = u/v give a = 2u(1+T)(1+A)/T and Q = 2v(1+T)(1+B)/T | EXACT (lambda = T/(2(1+T)), u = L/(1+L), v = 1/(1+L), a/e = A, Q/j = B) |
| 3 | a - 2u and Q - 2v are positive/positive, so the cell is automatic | EXACT |
| 4 | j <= 1 iff L >= 1 + 2/T + 2/B + 2/(TB); chart L = that + X, LN = BTL, K = LN + BT, K - 2(1+T)(1+B) = BTX, 1 - j = BTX/K | EXACT |
| 5 | onto the strict-failure region (inversion C1 = q/b, P = p/(lambda - p)) | PROOF-READ |
| 6 | D G with D = 4K^9 T^4 A^2 (1+T)^2 (1+C1)^6 (1+P)^2 is a polynomial with 374775 coefficients, all of one sign; the pure monomial T has coefficient 65536 (the document says -65536) | EXACT |
| 7 | X = 0 (j = 1) and c = 0 keep the marker; lambda = 1/2 goes to the large-lambda note; lambda = 0 is impossible | PROOF-READ |

## D. Falsification

NUMERIC: 12 random points with lambda < 1/2, j <= 1 and p = lambda (seed 14), min Phi (8 restarts): worst -1.1e-16 relative, so the capacity is sharp and no violation was found.

## E. Verdict

**confirmed-modulo-citations** (the zero-face cases). No gap.
