Evidence level: EXACT for the Gram algebra of Theorem 1, the density template, Corollaries 2-4 and the T47 transport (own code); CITED for the face's coefficient-relation system (foundations/so4-atoms-part2b.md G3, "Six-ray Thm 1 ... PROVED HERE"; "Singular leaf ... PROVED HERE"); PROOF-READ for the closure and projection arguments; UNCHECKED for the "exactly the finite closure" classification.

# final_referee_T46_six_ray_coherent_capacity_20261023 ("A three-scalar criterion for a six-generator Gram cone", 4 pp.)

Replay: `replay/final_referee_T46_six_ray_coherent_capacity_20261023__check.py` (output `..._check_output.txt`: 21 EXACT checks pass, plus 1 INFO line).

## Deviations and limits

- As in the core note, the face relation system (the "complete original completion equations ... rank nine" and the "156 coefficient rows") is CITED from foundations G3. I did not rebuild it.
- The literal action tuple does not reproduce the role map in moc4lib's convention. Exactly one g in G realises the full ordered 8-tuple map (27,1,59,5,2,28,6,4) -> (0,45,5,9,29,35,7,36). Phi_{g.c}(g.O) = Phi_c(O) was checked exactly at a rational O. This g does not map the T46 contact set to the historical T47 contact set, and it need not: the transported statement is a universal SO(4) inequality for a fixed cost.
- The classification ("the complete ... PSD8 completion cone on these six rays is exactly the finite closure") requires the nine equations to be the complete system. That is UNCHECKED. It is not used for sufficiency.

## A. Statements

- **Thm 1.** Let q, b > 0, l = qb - beta^2 > 0 and j - t^2 > (bt^2 - 2beta tv + qv^2)/l. Then q f27 + p f1 + j f59 + c f5 + f28 + b f6 satisfies Phi >= 0 on SO(4) whenever p >= p* and c >= c*. Here p* = [(j-t^2)h^2 - 2h(t+v)(1+v) + s(1+v)^2]/Delta with Delta = s(j-t^2) - (t+v)^2, and c* = (t-v)^2/(j - t^2 - v^2/b). Finite limits are included.
- **Classification.** The completion cone equals the closure of this family.
- **T47 transport.** Roles (q,b,j,p,c)_46 = (b,a,e,c,p)_47 (historical T47 = paper T_65).
- **Cor 2 (singular leaf, T47 roles).** For b > 1, 0 < a < 1/b and e > kappa t^2, the cost is safe when c >= k + (1-eta t)^2/(e - kappa t^2) and p >= rho^2 t^2/(e - kappa t^2). In particular it is safe when e <= p kappa^2/(eta^2(p kappa + rho^2)) and c >= 1/e + k - eta^2/kappa.
- **Cor 3.** For every p, e > 0, safe when c >= k + D0^2/(eN).
- **Cor 4.** For every p, e > 0, safe when c >= k + (D0-A)^2/(e(D0^2 - kappa A)) for A <= D0, and when c >= k for A >= D0, where A = eta^2 e.

## B. Dependencies

| dependency | status |
|---|---|
| T46 face relation system; the PSD-implies-SAFE principle | CITED: foundations G3 |
| singular-leaf Gram (rank-one leaf block, PSD only) | CITED: foundations G3 "Singular leaf (q != b): explicit PSD vectors, p* and c* formulas and the rational section p >= k + D0^2/(jN). PROVED HERE" |
| finite closure | foundations G3; re-read, correct |

## C. Steps

| # | claim | status |
|---|---|---|
| 1 | P.S = -h, P.J0 = -1-v, |P|^2 = p* | EXACT |
| 2 | C.B = C.Z = 0, C.J = t - v, |C|^2 = c*, |D|^2 = c*+1 | EXACT |
| 3 | the nine displayed relations hold for the construction | EXACT (that they are the full system: CITED) |
| 4 | the Schur condition equals the hypothesis | EXACT |
| 5 | surplus vectors orthogonal to everything add excess p and c | PROOF-READ |
| 6 | converse projection argument | PROOF-READ, given the nine equations are complete; completeness UNCHECKED |
| 7 | density template: positive definite, |D|^2 - |C|^2 = 1, satisfies the nine relations | EXACT |
| 8 | transport realised by an element of G; the images of the six roles are T47 generators, with f35 the non-minor one | EXACT |
| 9 | Cor 2: the substitution beta = -r, q = b47, b = a47, v = -rt/b47 gives l = 0 and both thresholds | EXACT |
| 10 | Cor 2: at t = eta e/kappa, the values and the equivalence p* <= p iff e <= p kappa^2/(eta^2(p kappa + rho^2)) | EXACT |
| 11 | Cor 3: the identities for N and N - rho^2 eta^2 e/p; e - kappa t^2 = eN/M0^2; the consumed c and p | EXACT |
| 12 | Cor 4: both branches, the consumed p, the comparison identity, and the brace at both endpoints with concavity | EXACT / PROOF-READ |

## D. Falsification

Not repeated. Foundations G3 found a minimum of 2.1e-18 on the threshold boundary.

## E. Verdict

**confirmed-modulo-citations** for Theorem 1 and Corollaries 2-4 and the transport. The citations are the face relation system and the singular-leaf Gram of foundations G3. The "exactly the finite closure" classification is UNCHECKED and not load-bearing.
