Evidence level: EXACT for the Gram algebra of Theorem 1, the optimization formulas, Corollary 3, the numeric example and the group transports (own code); CITED for the face's coefficient-relation system (foundations/so4-atoms-part2b.md G3, which proves Theorem 1 and Corollary 3); PROOF-READ for the closure and projection arguments; UNCHECKED for the completeness ("exactly") half of Theorem 2.

# final_referee_T46_core_coherent_capacity_20261023 ("A rank-four Gram bound for a five-generator source cone", 5 pp.)

Replay: `replay/final_referee_T46_core_coherent_capacity_20261023__check.py` (output `..._check_output.txt`: 22 EXACT checks pass, plus 3 INFO lines).

## Deviations and limits

- The face's degree-six coefficient-relation system (the "ten completion relations") was not rebuilt here. Foundations part 2b (G3) checked Gram G0 against the full linear system: the relations hold, |P|^2 equals threshold (1), and G0 is positive definite iff the Schur condition holds. I cite that and replay the remaining algebra.
- The literal (row, column, transpose) tuples printed for the transports do not reproduce the role maps in moc4lib's action convention. Some element of G does realise each ordered role map (2 elements each). This is a convention point, not a gap; foundations reports the same for other faces.
- Theorem 2's claim that the seven printed equations are the complete restricted system ("both of rank seven") is not checked. It is a classification statement. The sufficiency results (Theorem 1, Corollary 3) do not need it, and no document in this batch uses the necessity direction.
- Garbled extract: the example reads "q = 1/8, b = 4, j = 4, beta = -1/2, t = 1/4". This reading reproduces 399/1543 exactly.

## A. Statements

The face is historical T46, the contact set {1,3,5,6,7,8} = paper T_47. The core is p f1 + j f59 + q f27 + b f6 + f28; only f28 = h_{3421,-} is non-minor.

- **Thm 1.** Let q, b > 0, l = qb - beta^2 > 0, j > (1 + (q+b-2beta)/l) t^2, s = q+b+2beta, h = 1+beta, Delta = sj - (s+4)t^2. If p >= [(j-t^2)h^2 - 4th(t-1) + s(t-1)^2]/Delta, then Phi_c(O) >= 0 for every O in SO(4), for the cost c = p f1 + j f59 + q f27 + b f6 + f28 plus any nonnegative minor rays. This is a universal SO(4) inequality; it needs no spectral hypothesis. Finite limits are included.
- **Thm 2.** The PSD completion cone of this fixed-multiplier core equals the closure of the family in Thm 1.
- **Optimization remarks.** Form (2), the optimal t, and the formulas at t* = -kj/(s+4) and at beta = 1-q-b.
- **Cor 3.** If d > 0, j > 0, q, b >= 0 and (q+b-d)^2 <= qb, then the core is safe when p >= d^2/j + 2d - q - b. By transport, the T47 core is safe when (a+b-d)^2 <= ab, c >= d^2/e + 2d - a - b and e > 0.
- **Transports.** T24 roles (59,53,55,67,54) map to (1,59,6,27,28), which map to the T47 roles (45,5,7,0,35).

## B. Dependencies

| dependency | status |
|---|---|
| the face Gram relation system and the PSD-implies-SAFE principle (sum of cubic squares modulo the quaternion norm) | CITED: foundations/so4-atoms-part2b.md G3 ("Core coherent Thm 1 ... PROVED HERE") |
| finite Gram closure | foundations G3 (PROOF-READ there); re-read here, correct: fixed-diagonal slices are compact |
| T24 core (part 2a) as the source of the transport | foundations part 2a; only the transport was re-checked |

## C. Steps

| # | claim | status |
|---|---|---|
| 1 | Gram of (Q+B, J0) is [[s,-2t],[-2t,j-t^2]], with det Delta; P = lambda(Q+B) + mu J0 has P.(Q+B) = -h, P.J0 = t-1 and |P|^2 = (1) | EXACT |
| 2 | the seven printed completion relations hold for (Q,P,J,0,0,D,B,0) | EXACT (that they are the full system: CITED G3) |
| 3 | leaf Schur complement condition equals the hypothesis j > (1+u) t^2 | EXACT |
| 4 | t = 0 threshold and its derivative 2(2-q-b)/(sj) | EXACT |
| 5 | example: 399/1543 < 13/50 < the parallel value | EXACT |
| 6 | Thm 2: the template is positive definite, has |D|^2 = 1 and satisfies the seven equations; the convex-mix density and projection argument | EXACT / PROOF-READ |
| 7 | Thm 2: the seven equations are the complete restricted system (row-span equality, rank 7) | UNCHECKED |
| 8 | the stationary quadratic, form (2), its derivative, u - 4/s = (q-b)^2/(sl), the clipped minimizer | EXACT / PROOF-READ (sign analysis of the derivative re-done) |
| 9 | values at t* and at beta = 1-q-b | EXACT |
| 10 | Cor 3: 2/3 < S < 2; beta = 1-S, t = 0 gives 1/j + 2 - S; closure at the boundary and at S = 2; minimality over beta; discriminant b(4-3b) | EXACT / PROOF-READ |
| 11 | transports realised by group elements; T47 images are T47 generators | EXACT |

## D. Falsification

Not repeated. Foundations G3 attacked the core on its threshold boundary (minimum 1.1e-18; with the p-term removed, a control value of -0.054).

## E. Verdict

**confirmed-modulo-citations** for Theorem 1, the optimization formulas, Corollary 3 and the transports. The citation is the face relation system of foundations G3. Theorem 2's completeness claim is UNCHECKED; it is not load-bearing for the sufficiency statements used downstream.
