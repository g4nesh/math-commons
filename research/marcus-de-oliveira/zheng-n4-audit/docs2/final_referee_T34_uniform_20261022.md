Evidence level: EXACT for every identity and every positivity table, PROOF-READ for the logic that joins them (two short arguments supplied by the auditor), CITED for the four pointwise coefficient inequalities (3)-(6). Verdict: confirmed-modulo-citations.

# Second-pass audit: final_referee_T34_uniform_20261022 (primary proof of paper case T_35)

Document: `proofs/manin_transfer_20261021/final_referee_T34_uniform_20261022.pdf`, 4 pages, dated 22 October 2026. Auditor label: C14-T35-T37-second. This audit was done independently. The first auditor's report was not opened.

Replay: `replay2/final_referee_T34_uniform_20261022__check.py`, with its output in `replay2/final_referee_T34_uniform_20261022__check_output.txt`. The run passes 61 EXACT checks and 1 DATA comparison and reports the NUMERIC results; it takes about 200 s, almost all of it numeric. Everything is built from the paper's definitions through `lib/moc4lib.py`.

## Deviations and limits

- The four pointwise inequalities (3)-(6) are only CITED. They are proved in other supplement documents, which I did not audit. I checked two things about them. EXACT: a group element carries each cited base family onto the roles the document uses. NUMERIC: on their boundaries, a multistart minimization over SO(4) found no violation.
- I derived G myself: from the paper's cubic (A.2), through the full 14 x 14 Hessian determinant, to G (1475 terms). Only afterwards did I load the supplement's stored G, as JSON data, to compare. The two are identical, sign included. I did not run or import any supplement script.
- I rebuilt all nineteen positivity tables from my own G with my own substitutions. Every term count the document states is reproduced exactly.
- The stratum sampler is best effort. The final run found 6 points on the T_35 stratum and 1 of them has d > 0. An earlier scratch run with the same code found 12 points, 5 with d > 0.

## A. Statements

**Theorem 1 (p. 1), in the paper's notation.**
- Let T = {2, 3, 4, 6, 9, 12} (0-based lexicographic indices). This is literally the paper's T_35 of Appendix E; EXACT, no group element is needed.
- Take any x in X: any spectra a, b in C^4 and any alpha in C, including repeated eigenvalues and a_i + b_j = 0.
- Suppose P c(x) >= 0 and (P c(x))_pi = 0 for every pi in T. The contacts are non-strict, so further contacts are allowed.
- Then Phi_{c(x)}(O) >= 0 for every O in SO(4).

This is exactly (*) for T_35, that is, Proposition 4.3 for nu = 35.

The supplement writes gaps as Re{zeta(F_sigma - F_tau)} with tau in T. Absorbing zeta into x, and using Re F_tau = Re F_pi* on T, the gap vector is the paper's c(x).

**Covers all spectra.** Yes. The proof uses only two facts about x: x lies in X, and X is contained in {grad C = 0}. That inclusion holds on the whole Zariski closure, so degenerate spectra are included. Nothing in the proof is generic, perturbative or continuous in the spectra. The continuity steps are all in the finite coefficients of the cone.

**Internal claims strong enough to be false, all restated as replayed.**
- (2): Q^4 det Hess C(g - (C*/Q) kappa) = -4 Q^2 G^2.
- (7): k0 >= G00 (1 - L_e e)^2 (1 - L_j j)^2.
- (9): the analogous bound for G at c = 0, for every p >= 0.
- The Bernstein tables for b <= 1 and for b >= 1.

## B. Dependencies

| ref | kind | used for | status |
|---|---|---|---|
| Lemma 4.2 / Appendix D (generators of N_T) | core | the source {0,1,2,5,7,9,53,67,68} | EXACT, lib `generators` |
| eq. (4.2): Phi of a minor vector = det^2 >= 0 | core | d = 0 case; unused minor coefficients | PROOF-READ (lib README proof) |
| X contained in {grad C = 0} (Appendix A.1) | core | realizable => C(Re x) = 0, det Hess C(Re x) = 0 | EXACT in lib (16-variable identity) |
| (4.7) group invariance Phi_{g.c}(g.O) = Phi_c(O) | core | transporting the cited inequalities | EXACT in lib |
| four-coefficient Theorem 1 of final_referee_T5_four_budget_capacity_20261022: a f0 + b f6 + c f10 + d f54 safe if ab(c-d) >= d^2(a+b+c+d) | other-doc (T_5 = {1,2,3,4,6,8}) | (5) and (6) | CITED. EXACT: a group element maps (f0,f6;f10;f54) to ({7,0};f1;f68) and to ({7,67};f53;f68), and with d = 1 the condition is literally (5)/(6) |
| old-H family: t f3 + c f5 + d f54 safe if t(c-d) >= d^2 (the "[3,5,54]" family, as quoted in final_referee_T41_bad54 / blind_moc4_coreB) | other-doc | (3) | CITED. EXACT: a group element maps (f3,f5;f54) to (f9,f5;f68) in that order, so it gives c(p-1) >= 1. I did not locate the proof document itself. |
| three-coefficient chain: reference (f2; f0, f12; f70), center coefficient >= d^2/q + d^2/t | other-doc | (4) and its rescaling (14) | CITED. EXACT: a group element maps (f2;f0,f12;f70) to (f5;f7,f9;f68). I did not locate the proof document. |

No other contact case, branch document, Appendix B or C, or curve argument is used.

## C. Steps

| # | claim (page) | kind | status | how checked |
|---|---|---|---|---|
| 1 | The document's T equals paper T_35 (p. 1) | finite | EXACT | lib CONTACT_SETS[35] |
| 2 | The source is 0,1,2,5,7,9,53,67,68; rank 8; only f68 is non-minor (p. 1) | finite | EXACT | lib `generators`, `rank`, `F_IS_MINOR` |
| 3 | f68 + f2 = f0 + f9 + f67 is the only relation. So either d = 0 or there is the unique f2-free representation (1) with d > 0 (p. 1) | identity + case split | EXACT | 1-dimensional nullspace. Moving min(s2, d) of (f68 + f2) to (f0 + f9 + f67) leaves s2 = 0 or d = 0. |
| 4 | d = 0 gives a sum of minors, so Phi >= 0 | citation | PROOF-READ | (4.2) |
| 5 | For realizable cost g: Re x = g + h kappa satisfies C = 0 and det Hess C = 0 (p. 1, "therefore") | analytic | PROOF-READ, argument supplied by auditor | The document gives no argument; see note (a). |
| 6 | C(g + h kappa) = C* + hQ; C* has positive coefficients and contains d^3; Q as stated. So Q != 0 and h = -C*/Q (p. 1) | identity | EXACT | flint, from the paper's (A.2) |
| 7 | (2): Q^4 det Hess = -4 Q^2 G^2, with G of degree 9 and 1475 terms. So realizable implies G = 0 (p. 1) | identity | EXACT | Exact 14x14 determinant by subset DP, cross-checked at 3 integer points; the quotient is a perfect square; equals the stored G |
| 8 | Pointwise bounds (3)-(6) (pp. 1-2) | citation | CITED | group transport EXACT; NUMERIC (section D) |
| 9 | (5), (6) fail exactly when L_e e < 1 and L_j j < 1 (p. 2) | identity | EXACT | sympy |
| 10 | (7) for b > 0, c >= 0, in all four product cases (8); b = 0 factorization; walls ab = 1, bq = 1 (p. 2) | inequality certificate | EXACT tables + PROOF-READ | 2227/3268/3268/4945 terms, all >= 0. The (E,J)-free part equals the cleared G00 (G00 factorization EXACT). This gives k0 (1+E)^2(1+J)^2 >= G00 in the "above" cases and k0 >= G00 in the "below" cases. The walls follow by closedness of >=, with right-hand side > 0. |
| 11 | (9): G at c = 0 is > 0 for all p >= 0 in the hard region, so c = 0 is not realizable (p. 2) | inequality certificate | EXACT tables + PROOF-READ | 1419/2110/2110/3211 terms; the (E,J,p)-free part equals the cleared (1+b)^3(a+1)^3(q+1)^3; b = 0 factorization EXACT |
| 12 | Bernstein criterion (10): k0 > 0 and middle, endpoint >= 0 give G > 0 on [0, P*) (p. 2) | identity + inequality | EXACT + PROOF-READ | sympy identity; on [0,1) the weight (1-t)^2 is > 0 |
| 13 | b <= 1 with P* = 1 + 1/c: eight tables >= 0 (p. 3) | inequality certificate | EXACT | (3612,3612), (7230,7230), (7230,7230), (12258,12258). b = 0 and the product-one walls follow by continuity; see note (b). |
| 14 | Allocation (11)-(14): E_a increasing, its range, unique delta_i in [0,1), (13) is (5) rescaled, the delta1 + delta2 >= 1 case, the residual chain threshold (14) (p. 3) | identity + logic | EXACT identities + PROOF-READ | sympy. Derivative numerator a^2 b + 2ab delta - delta^2 > a^2 b > 0 for 0 <= delta < ab. The allocations are disjoint and sum to the f7 and f68 coefficients. |
| 15 | Chart (15)-(16) covers every point with b >= 1, a, q > 0, c > 0, delta1 + delta2 < 1 (pp. 3-4) | identity + logic | EXACT + PROOF-READ | sympy identities. Inverse map: U = delta1/r, V = delta2/r, X = ab - delta1 > 0, Y = bq - delta2 > 0, B = b - 1. |
| 16 | b >= 1 tables: base, middle, endpoint all >= 0 (p. 4) | inequality certificate | EXACT | 117110 / 209489 / 289848 terms |
| 17 | Finite boundaries (p. 4) | logic | PROOF-READ | Cases: a -> 0 or q -> 0 (0 <= delta1 < ab forces delta1 -> 0, so P* is continuous and the beta's >= 0 by limit); e = 0 or j = 0 (U or V = 0, inside the chart); r = 0 (the double-four case); b = 1 (both charts); equality in (5)/(6); p = 0 ((7) strict); c -> 0 (excluded by (9)). |
| 18 | Degenerate spectra | logic | PROOF-READ | Only x in X and grad C = 0 on X are used |

**Note (a), step 5, supplied argument.**
- Write grad C(w) = (1/2) H(w) w, where H(w) is the Hessian; it is linear in w and symmetric as a trilinear form. For x = u + iv in X, grad C(x) = 0 splits into H(u)u = H(v)v and H(u)v = 0.
- If v != 0, then H(u) is singular.
- Also 6 C(u) = u^T H(u) u = u^T H(v) v = v^T H(u) v = 0.
- If v = 0, then u lies in Sing C, so C(u) = 0 and H(u)u = 0. Here u != 0, because the eight generators are independent and d > 0.
- In both cases C(Re x) = 0 and det H(Re x) = 0, with Re x = c(x) + Re F_tau(x) kappa.

**Note (b), step 13, supplied argument.** The chart b = B/(1+B) covers 0 < b < 1, ab != 1 and bq != 1. The Bernstein coefficients are continuous in the original coefficients for fixed c > 0. Each excluded boundary point (b = 0, ab = 1, bq = 1) is a limit of hard chart points with the same remaining coefficients, because L_e = 0 whenever ab <= 1. So the conditions beta >= 0 extend to those points. b = 1 is also covered by the b >= 1 chart.

Neither supplied argument changes the document's method. Both are recoverable in a few lines; I record them here because the document asserts these steps without proof.

## D. Falsification (NUMERIC)

- **Pointwise (3)-(6) on their boundaries.** 12 random coefficient tuples per inequality, 12 BFGS restarts each over SO(4), seed 7. Worst minimum of Phi: 2.3e-14 for (3), 3.9e-14 for (4), 1.8e-15 for (5), 5.3e-16 for (6). No violation was found. At 0.97 times the threshold, (3) and (5) go slightly negative (-4.2e-8 and -7.8e-11), which is consistent with those thresholds being sharp or nearly so.
- **Realizable points on the T_35 stratum.** lib `random_point_on_contact_stratum` with margin 0.01, so the contact set is exactly T.
  - Final run: 6 points, 1 with d > 0. Earlier scratch run: 12 points, 5 with d > 0.
  - At every d > 0 point, |G| / sum|terms| is at most 5.7e-14. This confirms the necessity G = 0 numerically.
  - At the final-run d > 0 point, one of the document's safe certificates applies (the four-coefficient one).
  - The minimum of Phi over SO(4) at every point is >= 2.6e-13 (10 restarts). No counterexample was found.

## E. Verdict

**confirmed-modulo-citations.**
- Every identity and polynomial-positivity certificate in the document is replayed exactly from the paper's definitions, with the stated term counts reproduced.
- The case logic is exhaustive. Two short steps (note (a), note (b)) are justified by arguments the auditor supplied.
- The proof rests on three cited pointwise families, none audited here:
  - the four-coefficient Theorem 1 of final_referee_T5_four_budget_capacity_20261022;
  - the old-H family (f3,f5;f54);
  - the chain family (f2;f0,f12;f70).
- Their group transports to (3)-(6) are EXACT, and they pass the numeric boundary tests. No GAP was found.
