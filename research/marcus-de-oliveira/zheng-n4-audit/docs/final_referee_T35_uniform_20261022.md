Evidence level: EXACT for all identities, both cubic restrictions and their guards, and all 10 positivity tables (own code; the Hessian factorizations (2) by an exact randomized identity test); PROOF-READ for the reduction, chart coverage and boundary continuation; CITED for the pointwise capacities (3), (4) and for X in Sing C; NUMERIC for falsification.

# Audit of `final_referee_T35_uniform_20261022` (primary proof of paper case T_36)

PDF: `proofs/manin_transfer_20261021/final_referee_T35_uniform_20261022.pdf`, 4 pages. Replay: `replay/final_referee_T35_uniform_20261022__check.py` (46 EXACT checks pass, exit 0).

## Deviations and limits

- G68 (1864 terms) and G69 (1367 terms) were read as data from `final_referee_T35_bad68/bad69_full_radial_20261022.json`. Their stored sign has d^9 coefficient -1, and I negated them as the document says.
- Each was validated against my own 14x14 Hessian of (A.2) at 10 random integer points, checking Q^2 det Hess = -4 G^2 exactly (randomized identity test, as in the T34 report).
- My C* and Q equal the JSON restrictions.
- The capacities (3) and (4) and their group witnesses are not audited: CITED, with numeric tests only.
- Budget: about 20 minutes.

## A. Statement

Theorem 1 (p. 1). T = {2,3,4,6,11,12}. This is literally paper T_36. The statement is (*) for T_36: non-strict contacts, all O in SO(4), all spectra. The method uses only the necessary conditions C = det Hess C = 0 on X, so degenerate spectra are included.

## B. Dependencies

- Lemma 4.2. The source 0,1,2,7,8,9,10,53,67,68,69 equals lib `generators(T_36)` (EXACT).
- X subset Sing C, plus the secant/Hessian necessity argument. PROOF-READ, with the argument supplied in the T34 report.
- Branch (68), CITED:
  - four-coefficient qj(p-1) >= q+j+p+1, roles (0,7;1;68)
  - four-coefficient jb(e-1) >= j+b+e+1, roles (7,67;53;68)
- Branch (69), CITED:
  - SOC pj >= 1 (2,8;69)
  - old-H c(p-1) >= 1 (9,2;69)
  - old-H atoms (10,8;69) and (53,8;69). Their H-sum combination gives the central requirement 1 + 1/(a+e). I re-derived this: eta = (a, e)/(a+e) gives eta1 + eta2 + eta1^2/a + eta2^2/e = 1 + 1/(a+e). PROOF-READ.
  - chain j >= 1/e + 1/b (center 8, leaves 53, 67)
- Not invoked: Appendix B and Corollary B.3.

## C. Steps

| # | claim | kind | status | how |
|---|---|---|---|---|
| 1 | Source, rank 8, nonminor 68 and 69 | identity | EXACT | |
| 2 | The five cancellation identities | identity | EXACT | 24-value checks |
| 3 | Cancellation order. If 69 is exhausted first, consume 2 and 8 against 68, leaving the basis {0,1,7,9,53,68,67,10}. Otherwise consume 0 and 7 against 69, leaving {1,2,8,9,53,69,67,10}. Both bases are independent | case split | EXACT + PROOF-READ | ranks; the coefficient bookkeeping only increases minor coefficients |
| 4 | (68): C* and Q have positive coefficients (d^3, d^2) | identity | EXACT | |
| 5 | (69): the only negative monomial of C* is -aej, and C* >= 1 + a + e + ae(1-j). Under H-sum failure, C* > 0 | identity + inequality | EXACT + PROOF-READ | the bound 1 + S - ae/S >= 1 + 3S/4 re-derived |
| 6 | (2) for both branches; d^9 normalization | identity | EXACT (randomized) | 10 points each |
| 7 | (5) baseline G68 at p=e=a=c=0 equals (j+1)^3(b+1)^3(q+1)^3; degrees (3,2,3,2,2,3,2) | identity | EXACT | |
| 8 | Charts (6) cover the strict hard set of branch (68) away from the walls | coverage | PROOF-READ | (3a) <=> p >= (q+1)(j+1)/(qj-1); the upper chart realizes this as p = threshold * P/(1+P) |
| 9 | Chamber tables with Z counts 6190, 8949, 8705, 12375 | tables | EXACT | identical counts; every multiplier M has nonnegative coefficients |
| 9a | Chamber 00: Z - M H^2 has exactly three -1 monomials at the N of (7), and 2(Z - MH^2) - sum x^(nu - e_J)(J-1)^2 is nonnegative (5893) | tables | EXACT | |
| 9b | Chambers 01, 10, 11: Z - MH^2 nonnegative (8802, 8439, 12197) | tables | EXACT | |
| 10 | Every subtracted term contains P*A*C, so (5) follows with the factors (1+P)^-2 = (1 - Lp p)^2 and (1+E)^-2 = (1 - Le e)^2 | argument | PROOF-READ | |
| 11 | (9) B0 = G69 at a=c=b=q=0; B0 > 0 when pj < 1 | identity | EXACT + PROOF-READ | |
| 12 | Charts (11) cover the open hard chambers of branch (69) | coverage | PROOF-READ + EXACT | The 1 - L x factors equal 1/(1+A), 1/(1+C), 1/(1+B) (EXACT). Exhaustiveness: j > 1 forces p < 1, and p > 1 forces j < 1 |
| 13 | Six tables: (3226 -> 3199 after the square (12)), 41018, 14406, 26220, 6760, 22059 | tables | EXACT | identical counts |
| 14 | Continuation to walls and zero coefficients: Lp, Le, La, Lb, Lc are continuous on the hard set, and the bounds (5), (10) extend by continuity of polynomials | continuity | PROOF-READ | |

## D. Falsification (NUMERIC)

- All six capacity statements at equality: 6 samples each. Worst normalized min Phi lies between 3e-16 and 1.7e-13, i.e. zero to rounding.
- G68/bound - 1: min 0.070 over 20000 hard points.
- G69/bound - 1: min 0.27 over 20000 hard points.
- Main statement: 3 stratum points, min Phi >= -2e-17.

## E. Verdict

**confirmed-modulo-citations.** Everything checkable was replayed exactly with the document's counts. The citations are the capacity inequalities (3)-(4), Lemma 4.2, X subset Sing C and the reduction to SO(4). No gaps found.

Presentation note: the source sign convention for G (stored with d^9 coefficient -1) is stated in the text, and the replay follows it.
