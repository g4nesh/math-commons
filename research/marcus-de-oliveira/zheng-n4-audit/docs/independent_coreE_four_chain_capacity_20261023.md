# Audit: independent_coreE_four_chain_capacity_20261023 (four-coefficient plus chain capacity on E)

Evidence level: EXACT for every displayed identity, the capacity matrix, the endpoint reductions and the worked example; PROOF-READ for the positivity and allocation logic; CITED for the four-coefficient and chain reference inequalities; NUMERIC for boundary falsification.

PDF: proofs/manin_transfer_20261021/independent_coreE_four_chain_capacity_20261023.pdf (2 pages). Replay: `replay/independent_coreE_four_chain_capacity_20261023__check.py` (output in `__check_output.txt`, 20 PASS, 0 FAIL).

## Deviations and limits

- The two reference inequalities ("four-coefficient" and "chain") are CITED. Their ordered group actions onto the literal families used here were replayed in the core E replay (`independent_coreE_uniform_20261023__check.py`). Their source proofs were not audited.
- This note is a pointwise (non-realizability) sufficient condition. It proves no case of (*) by itself. Core E reproves its special section delta = dv/(c+v) as its Lemma 1.

## A. Statements

Source on E = {0,1,2,9,10,16,17}: c = q f0 + p f2 + t f3 + b f12 + c f15 + v f13 + d f26. Here generators(E) = {f0,f2,f3,f12,f13,f15,f26} (EXACT).

- **Theorem 1 (p. 1).** Let q, v, c, d > 0 and p, b, t >= 0. Suppose some delta in [0, d] satisfies P(delta) = cvt(q+delta) - (t+q+delta){cv delta + c delta^2 + v(d-delta)^2} >= 0. Then Phi_c(O) >= 0 for every O in SO(4), for every such vector c in V_R, with no realizability assumed. The same holds with q and v interchanged (using the chain with v f13).
- **Capacity matrix (p. 2).** For the same parameters, P(delta) >= 0 is equivalent to M(delta) ⪰ 0, where M(delta) is the displayed 4x4 matrix: the Schur complement in the last two diagonals is [[W, K], [K, t-K]] and det M = P(delta).
- **Simple section.** delta = dv/(c+v) gives the sufficient condition q(c+v) > d^2 and t >= d(v+d){q(c+v)+dv} / ((c+v){q(c+v)-d^2}). This is core E's Lemma 1.
- **Optimality within the section.** For fixed delta, the allocation minimizes A + B + A^2/X + B^2/(W-X) over 0 < X < W.
- **Example.** d = 1, q = 100, v = c = 1/50: every such source with t >= 680/19 is safe.

## B. Dependencies

| dependency | status | match |
|---|---|---|
| four-coefficient reference (f0,f6,f10;f54) mapped to (f0,f13,f3;f26): q f0 + v f13 + t f3 + d f26 safe if qv(t-d) >= d^2(q+v+t+d) | CITED | ordered action EXACT (core E replay) |
| chain reference (f2;f0,f12;f70) mapped to (f15;f0,f3;f26) and (f15;f13,f3;f26): c >= d^2/q + d^2/t | CITED | ordered actions EXACT (core E replay) |
| minor vectors have Phi = det^2 >= 0 (paper (4.2)) | EXACT in moc4lib | |

## C. Steps

| # | claim | status | how checked |
|---|---|---|---|
| 1 | (1) is exactly t W >= (q+delta) K, with P = cv(tW - (q+delta)K) | EXACT | sympy |
| 2 | for 0 < delta < d: K > 0 and q + delta > 0, so P >= 0 forces t > 0 and W > 0; then X = WA/K in (0, W) | PROOF-READ | elementary |
| 3 | q_F v (t_F - delta) = delta^2(q_F + v + t_F + delta), lambda^2/q_C + lambda^2/t_C = c, and t_F + t_C = (q+delta)K/W <= t | EXACT | sympy |
| 4 | q_F v - delta^2 = Xv > 0 and q_C - B = W - X = WB/K > 0, so the allocations are finite and positive | EXACT and PROOF-READ | |
| 5 | Sum: a four atom (non-minor weight delta), a chain atom (weight lambda), and minor remainders (t - t_F - t_C, p, b) | PROOF-READ, using CITED (3) and (4) | linearity of Phi |
| 6 | delta = 0: P(0) = v{cqt - d^2(q+t)}, which is the chain condition. delta = d: P(d) = c{qv(t-d) - d^2(q+v+t+d)}, which is the four-coefficient condition | EXACT | sympy |
| 7 | Schur complement and det M = P(delta); the equivalence with PSD follows because v, c > 0 and K > 0 | EXACT and PROOF-READ | sympy |
| 8 | the section delta = dv/(c+v) gives W = q - d^2/(c+v) and the threshold R_d | EXACT | sympy |
| 9 | A^2/X + B^2/(W-X) >= (A+B)^2/W, with equality at X = WA/(A+B) | EXACT | explicit square identity |
| 10 | Example: the four inequality is strict (78/5 > 1513/100); (1/2)^2/60 + (1/2)^2/(300/19) = 1/50; the allocations exhaust q = 100; the total t is 680/19. 680/19 is below both separate thresholds (four alone 5151/50, chain alone 100) | EXACT | Fractions |

## D. Falsification (NUMERIC)

- 60 random (q, v, c, delta) with d = 1, on the boundary t = (q+delta)K/W (P(delta) = 0), p = b = 0. Seed 11, 8 BFGS restarts. Worst min Phi / max|Pc| = 1.7e-18.
- The worked example (t = 680/19) with 30 restarts gives 2.2e-17.
- No violation was found.

## E. Verdict

**confirmed-modulo-citations.** Every algebraic step is EXACT and the logic is complete. The result rests only on the CITED four-coefficient and chain reference inequalities, whose literal images were verified exactly. One remark, not a gap: the claimed equivalence with the capacity matrix being PSD is correct but is not used in the proof.
