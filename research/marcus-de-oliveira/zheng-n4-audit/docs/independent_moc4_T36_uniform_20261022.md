Evidence level: EXACT for the three-way reduction, all identities and guards, both Hessian factorizations (exact randomized identity test), the existence of every stated group transfer, and all 85 positivity tables plus 4 R90 tables; PROOF-READ for case coverage, the KKT chart coverage and the zero-coefficient limits; CITED for the transferred T24 core bounds, the old-H, chain, SOC and four-coefficient inequalities, and X in Sing C; NUMERIC for falsification.

# Audit of `independent_moc4_T36_uniform_20261022` (primary proof of paper case T_37)

PDF: `proofs/manin_transfer_20261021/independent_moc4_T36_uniform_20261022.pdf`, 6 pages. Replay: `replay/independent_moc4_T36_uniform_20261022__check.py` (59 EXACT checks pass, exit 0).

## Deviations and limits

- G35 (983 terms) and G40 (1385 terms) were read as data from the supplement's `bad35/bad40_full_radial` JSON files and negated, as stated (d^9 coefficient +1). Each was validated against my own Hessian of (A.2) at 10 random integer points.
- Row II of the seven-interval table: my middle and endpoint counts are 6832 and 9514, against the document's 9616 and 16507.
  - The cause: I write the threshold J = 1/e - 1/Emax in lowest terms before clearing denominators. The document clears a larger, unreduced denominator.
  - All my coefficients are nonnegative. Since N and D have positive coefficients, this is an equally valid certificate. Not a gap.
  - All other counts agree with the document.
- The cited T24 core bounds (6)-(7) and the other capacities were tested numerically only.
- Budget: about 30 minutes.

## A. Statement

Theorem 1 (p. 1). T = {0,5,6,7,9,12}. This is literally paper T_37. The statement is (*) for T_37 at full strength: non-strict contacts, all O in O(4) (so also SO(4)), all spectra including coincidences. The method uses only necessary conditions valid on all of X.

## B. Dependencies

- Lemma 4.2. The generators of T_37 are {0,1,2,5,7,29,35,39,40,45}, rank 8 (EXACT). The document proves its decompositions directly as linear identities on the 8-dimensional contact nullspace, and I replayed those exactly, so Lemma 4.2 is not even needed.
- X subset Sing C, and the necessity argument of Section 2. The polarization identity H(x+y)(x-y) = 0 was checked by hand: PROOF-READ.
- Branch delta35. The T24 core bounds (6) cj >= 1, c >= 1 + 1/e and K >= 0, and the coherent bounds (7), all transferred by the core action. CITED.
  - The group element exists and is unique (EXACT).
  - The document's printed row and column labels match lib's convention only after reordering/inverting. The other actions match directly. Presentation only.
- Branch delta40. The two old-H and two chain families (transposed core action: element exists, EXACT) and SOC as >= 1 (element exists, EXACT). CITED.
- The exchange stabilizer (a,b,q,p) <-> (b,a,p,q). The group element exists, and the symmetry of G35 is exact.
- Theorem 1.2 (two-orthogonal law), for the return to U(4). Paper.
- Not invoked: Appendix B and Corollary B.3.

## C. Steps

| # | claim | kind | status | how |
|---|---|---|---|---|
| 1 | delta35 + delta40 = -g17; (1) is a minor sum; (2) and (3) are exact operator identities | identity | EXACT | on a basis of the contact nullspace |
| 2 | The three generating identities | identity | EXACT | |
| 3 | The branches are mutually exclusive and exhaustive | case split | PROOF-READ | follows from step 1 and g17 >= 0 |
| 4 | Section 2: C restrictions match the JSON; in (2) every C coefficient is negative (incl. -d^3); (5) C = bc(U-a) + U T0 and the M C identity, whose last parenthesis has nonnegative coefficients; C at M = 0 | identity | EXACT | |
| 5 | Guard for (3): energy failure plus J <= 1 + 1/(b+c) gives a < d + d^2/M, so C > 0 | inequality | PROOF-READ | |
| 6 | (4) for both branches | identity | EXACT (randomized) | |
| 7 | k0 = Cbar * R90; G35 at e = j = 0 factorization; exchange symmetry | identity | EXACT | |
| 8 | (8): four R90 tables (280, 404, 280, 404) nonnegative | tables | EXACT | identical counts |
| 9 | The seven-interval rows cover the hard set of branch delta35 | coverage | PROOF-READ | j >= e (I+/-); ab > 1, j <= e (II); ab < 1, j <= e with b >= a by exchange, split at b = 1 and at e/(1+e) or eb/(e+b) (III). Each J is a cited valid threshold |
| 10 | 21 Bernstein tables nonnegative; each J = N/D with positive N, D | tables | EXACT | counts identical except row II (see Deviations) |
| 11 | Closed layers: e = j = 0 identity; zero leaves; ab = 1; b = 1 | continuity | PROOF-READ | |
| 12 | (10) allocation | identity | PROOF-READ | The required central coefficients add to the energy. Any feasible (beta, delta) gives a valid sufficient threshold, so optimality is not even needed |
| 13 | Five KKT charts: J equals the energy of an explicit feasible allocation (beta = (theta, 1-theta), 0 <= delta_i <= beta_i) | identity | EXACT | |
| 14 | Coverage of every positive (p,b,c,q) by the charts | coverage | PROOF-READ, argument completed by auditor | The reduced objective is strictly convex in beta1. I re-derived that the HC and CC stationarity equations are exactly the chart relations |
| 15 | 60 tables (L0, B1, B2 by Y-power) nonnegative | tables | EXACT | identical counts |
| 16 | L0(0,r) has positive coefficients and constant term 1, giving (12) | identity | EXACT | |
| 17 | Identity (13), giving (14) | identity | EXACT | |
| 18 | Zero core coefficients: J_eps increases to J; limit of (14) | analytic | PROOF-READ | |

## D. Falsification (NUMERIC)

- Capacities at equality: 8 statements x 6 samples, min Phi over SO(4) with 6 restarts. Worst normalized value is -2.4e-16, i.e. zero to rounding. This covers the cited coherent bounds (7) and the energy threshold J (numerically minimized).
- (14): G40 - (1-as)^2(1-a/J)^2 has min 1.71 over 20000 random points, with J the numerically minimized energy (10). This tests the tables and the chart coverage together.
- Main statement: 3 stratum points, min Phi >= -2e-17.

## E. Verdict

**confirmed-modulo-citations.** The citations are the T24 core bounds (6)-(7) (including the coherent bounds), old-H, chain, SOC and four-coefficient inequalities, X subset Sing C, and Theorem 1.2. No gaps found.

Presentation notes:
- The row/column labels of the core action do not follow one convention throughout. A valid group element exists in every case.
- The row II term counts depend on how the denominator is cleared.
