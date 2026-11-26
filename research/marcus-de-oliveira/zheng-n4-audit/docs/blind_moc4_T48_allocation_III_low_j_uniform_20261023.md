Evidence level: EXACT for the allocation algebra, stationarity, the chart formulas and their coverage, the forced-height exclusion, the zero-face reductions and the atom-orbit membership; EXACT pointwise for det H = -4G^2/Q^2 (own Hessian, supplement G as data); NUMERIC-level (exact rational evaluation at sampled chart points) for the strict negativity of G below p*; PROOF-READ for Corollary 2; CITED for the three large negativity certificates and for allocations I, II, IV.

# Audit: blind_moc4_T48_allocation_III_low_j_uniform_20261023 (branch of paper T_49)

PDF: proofs/manin_transfer_20261021/blind_moc4_T48_allocation_III_low_j_uniform_20261023.pdf, 3 pages. Replay: `replay/blind_moc4_T48_allocation_III_low_j_uniform_20261023__check.py`.

## Deviations and limits

- The exclusion below p* rests on three cleared integer polynomials with 163 753, 189 096 and 309 412 coefficients, all claimed negative. I did not rebuild them (budget). Instead I evaluated the supplement's 1772-term G exactly at random points of all three charts, including c = 0, A = 0 and j = 1 (W_V = 0); see D. This is evidence, not proof.
- G is read as JSON data from `blind_moc4_T48_native_full_radial_20261023.json`; its identity with the paper's cubic, det H(g - C kappa/Q) Q^2 = -4 G^2, was checked exactly at random rational points with my own Hessian (pointwise exact; as a polynomial identity only with high probability).
- The supplement's labels "T38, T0, T3, T23, T42, T43" for the zero faces are historical labels; they are not all the paper index minus one. I did not try to decode them; instead I checked directly that every zero face lands in a six-set of an earlier paper stage.

## A. Statements

Source roles (q,p,j,c,e,d,b,a) = (f7,f5,f29,f39,f2,f40,f41,f0) on T = {0,3,5,6,9,12} = paper T_49 (EXACT). This is the canonical form (1) of root_true_T48_whole_assembly after the substitutions f40+f1 = f7+f39+f0, f40+f9 = f7+f2+f41.

**Theorem 1 (cell III).** For realizable supported costs of this form with d = 1, a,e,j,q,b > 0, c >= 0, Q = q+b, u = (1+2/e)/(2(1/Q+1/j+1/e)), v = 1-u: if 0 < u < 1, Q >= 2u, a <= 2v, j <= 1, then Phi_{c(x)}(O) >= 0 on SO(4). Including all equality walls and c = 0.

**Corollary 2.** Every realizable canonical source with 0 <= j <= d is safe, given allocation I, II and IV theorems.

## B. Dependencies

| ref | kind | used for | status |
|---|---|---|---|
| chain atom, tuples (5;7,29;40), (5;41,29;40), (5;0,2;40) | so4-atoms (proved) | pointwise bound | EXACT: all three in the proved 576-orbit with centre and leaves as used |
| old-H atom, (aux f2 = e, main f5 = p; f40) | so4-atoms (proved) | pointwise bound | EXACT: in the 288-orbit; the reversed roles are not, and the document uses the correct order |
| radial necessity (established item 5) + det identity | core | realizable implies G = 0 | established; identity pointwise EXACT |
| allocation-I "both-chain" j <= 1 cell (two lambda notes) | branch | wall a = 2v; cell I | CITED |
| allocations II, IV | branch | Corollary 2 | CITED |
| zero faces q,p,j,e,b,a | earlier cases | boundary | EXACT: extra contacts 16, 11, 2, 13, 7, 23; each T_49+{j} contains a six-set of paper stage <= 3 |
| c = 0 face | seven-contact frontier | included in the theorem | (handled inside the theorem; T_49+{17} = Table 5 set 23) |
| three negativity certificates | supplement verifier | p < p* exclusion | CITED (not rebuilt) |

## C. Steps

| # | claim | kind | status | how |
|---|---|---|---|---|
| 1 | stationarity u(1/Q+1/j) - v/e = 1/2 is the definition of u | identity | EXACT | sympy |
| 2 | allocation: (q,j),(b,j) chains with currents uq/Q, ub/Q need u^2(1/Q+1/j); old-H current h = v - a/2 plus (a,e) chain current a/2, e shared in ratio h : a/2, need v - a/4 + v^2/e | inequality bookkeeping | EXACT | sympy on the stated shares; each atom with roles in the proved orbit |
| 3 | p* = 1 - u/2 + v/e - a/4 > 0 | identity | EXACT | sympy |
| 4 | a = 2v wall is the allocation-I cell | citation | CITED | stationarity identity reduces correctly (PROOF-READ) |
| 5 | j = 2ue(1+A0)/(A0 e + 2v(1+A0)) with Q = 2u(1+A0) | identity | EXACT | sympy |
| 6 | three charts cover j <= 1: K <= 0 (chart 1, j < 1 automatically), L >= 1 (chart 2), L < 1 & K > 0 (chart 3) | case split | EXACT + PROOF-READ | chart 1: den - num of j is a positive polynomial; chart maps K(X/(1+X), 2X+A) = -A/(1+X) and K(.., 2XA/(1+A)) = 2X/((1+X)(1+A)); j(e) is increasing with supremum 2u(1+A0)/A0 >= 1 exactly when K > 0 or L >= 1, so the j-parametrized charts reach every e with j <= 1 |
| 7 | chart e-formulas | identity | EXACT | sympy solve |
| 8 | realizable implies G = 0 (forced-height exception excluded by (e+j)C - ejQ > 0) | identity | EXACT | own cubic: all coefficients of (e+j)C - ejQ positive (term count in output) |
| 9 | below p*, G < 0 on all three charts including walls and j = 1 | certificate | CITED (UNCHECKED here) | pointwise sampling in D |
| 10 | p = p* handled by the pointwise bound | inequality | EXACT (steps 2-3) | |
| 11 | zero faces | reduction | EXACT | see B |
| 12 | Corollary 2: T(s) strictly convex, T'(0) < 0 < T'(1), four cells by s <= a/2, 1-s <= Q/2; cell III minimizer is the stationarity u | analytic | PROOF-READ | I re-derived T(s) from the atoms: the e-side cost with current s is s^2/e + s - phi_a(s), the Q-side cost with current 1-s is (1-s)^2/j + (1-s) - phi_Q(1-s) |

## D. Falsification

Exact rational evaluation of G at random points of the three charts (seed 23), with c = 0 in 20% of draws, A = 0 in 15%, and j = 1 in 15% of chart-2/3 draws; every sampled point was first checked to satisfy the Theorem-1 hypotheses and the stationarity identity exactly. All 1 800 sampled points (600 per chart) gave G < 0; the largest value was -2.39 (chart 1). No positive or zero value of G was found.

## E. Verdict

**confirmed-modulo-citations.** All algebra stated in the document is EXACT, the atoms are used in their proved orbits with the printed roles, and the charts cover the region. The core exclusion rests on the uncited-in-PDF negativity certificates (CITED; no counterexample in sampling), and Corollary 2 on allocations I, II, IV (CITED).
