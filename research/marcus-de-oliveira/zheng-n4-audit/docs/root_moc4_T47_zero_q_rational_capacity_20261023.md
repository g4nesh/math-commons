Evidence level: EXACT for the capacity-curve identities (c* = Crat, p* < p, positivity of the denominators); NUMERIC-level (exact rational evaluation at sampled points) for the main inequality G(0, ...) > 0; CITED for the six-ray Gram construction and for the 1 675 982-term positivity certificate, which I did not rebuild.

# Audit: root_moc4_T47_zero_q_rational_capacity_20261023 (branch lemma for paper T_65)

PDF: proofs/manin_transfer_20261021/root_moc4_T47_zero_q_rational_capacity_20261023.pdf, 2 pages. Replay: `replay/root_moc4_T47_zero_q_rational_capacity_20261023__check.py`.

## Deviations and limits

- The theorem is a polynomial positivity statement whose only proof is a cleared polynomial with 1 675 982 positive coefficients produced by the supplement's verifier. Rebuilding it needs the Crat substitution into the 431-term octic with high denominator powers; I judged this beyond the batch budget and did NOT replay it. In its place I evaluated G exactly (rational arithmetic) at 2 500 sampled points of the closed region, about half of them on one or more walls. This is evidence, not proof.
- G is the supplement's 965-term polynomial, read from `root_moc4_T47_full_radial_20261022.json` as data. Its link to the paper's cubic (det HF(w(g,-C/Q)) Q^2 = -4 d^2 G^2) is replayed exactly at 25 random rational points in the companion report root_moc4_T47_negative_Phi_20261022.
- "Historical T47" = paper T_65 = {0,4,5,6,9,12}; the roles (q,p,j,c,e,d,b,a) = (f2,f9,f1,f45,f5,f35,f0,f7) are exactly the generators of T_65 (EXACT, companion replay).

## A. Statements

**Theorem 1 (p.1).** With d = 1: if b > 1, 0 < a < 1/b, p >= b, 0 < j <= 1, e >= bj/(b-j), 0 <= c < Crat(p,e,b,a), then G(0,p,j,c,e,1,b,a) > 0. A pure algebraic inequality on a real region; it asserts nothing about realizability, and the document says so (a separate monotonicity argument is needed to use it at q > 0).

**Unnumbered claim (p.1).** For every p, e > 0 (with the other parameters as above) the cost with c >= Crat satisfies the support inequality, by the six-ray cubic-Gram construction: at t = eta e/(D0 + eta^2 e) the singular-leaf family gives c*(t) = Crat and p*(t) < p, and the construction proves safety for c >= c*, p >= p*.

Use downstream: on the q = 0 face of the T_65 source, below the capacity Crat a realizable point would need G = 0 (radial necessity); Theorem 1 excludes it; at or above Crat the six-ray atom gives safety.

## B. Dependencies

| ref | kind | used for | status |
|---|---|---|---|
| final_referee_T46_six_ray_coherent_capacity_20261023 | common lemma (face Gram family) | safety for c >= c*(t), p >= p*(t) | CITED; the so4-atoms foundation lists the six-ray family as tested numerically only, not proved |
| root_moc4_T47_full_radial (G) | data | the octic | link to the cubic replayed (companion report) |
| root_moc4_T47_zero_q_rational_capacity_exact (1 675 982-term certificate) | supplement verifier | Theorem 1 | CITED (not executed, not rebuilt) |

## C. Steps

| # | claim | kind | status | how |
|---|---|---|---|---|
| 1 | c*(t) = Crat at the chosen t | identity | EXACT | sympy |
| 2 | e - kappa t^2 = e N/(D0+eta^2 e)^2 > 0 | identity | EXACT | sympy; N > 0 by step 4 |
| 3 | N - (rho^2/p) eta^2 e = D0^2 + D0 eta^2 e + eta^4 e^2 > 0, hence p* < p | identity | EXACT | sympy (uses D0 - kappa = rho^2/p) |
| 4 | N > D0^2 > 0 | identity | EXACT | N - D0^2 = eta^2 e(kappa + 2 rho^2/p + eta^2 e) |
| 5 | safety at c >= Crat (larger c, p only add minors f45, f9) | citation | CITED | six-ray family; monotonicity in minor coefficients PROOF-READ |
| 6 | chart b = 1+A, r = V/(1+V), ..., covers the closed region including j = 1 via T/(T+Wj) | parametrization | PROOF-READ | each map is onto its interval; c = Crat S/(1+S) covers [0, Crat) |
| 7 | cleared polynomial has 1 675 982 positive coefficients; monomial 18 T^10 survives all walls | certificate | CITED (UNCHECKED here) | not rebuilt |
| 8 | G(q=0, d=1) has 431 terms | identity | EXACT | count matches |

## D. Falsification

Exact rational evaluation of G(0,p,j,c,e,1,b,a) at 2 500 points of the region (seed 11), 1 316 of them on at least one of the walls p = b, j = 1, e = bj/(b-j), c = 0, with b down to 1 + 10^-5 and r down to 10^-6 and up to 1 - 10^-6: every value positive (minimum 67.46). Just outside the region (c = 1.01, 2, 10 times Crat, seed 12) 14 of 1 200 points give G < 0, so the inequality is not vacuous and the capacity matters.

## E. Verdict

**confirmed-modulo-citations.** The identities that the document states in its text are EXACT. The theorem itself rests on the supplement's 1.7-million-coefficient certificate, which I did not rebuild; my exact pointwise sampling found no counterexample. The "c >= Crat is safe" claim rests on the six-ray Gram construction, which this audit has only tested numerically.
