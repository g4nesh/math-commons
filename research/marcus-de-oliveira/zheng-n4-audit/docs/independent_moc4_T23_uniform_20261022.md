# Audit: independent_moc4_T23_uniform_20261022 (paper case T_24)

Evidence level: EXACT for all finite and rational-function identities (source, decomposition, boundary list, five group witnesses, chart, inverse chart, pole structure, residue ratio); PROOF-READ for the harmonic minimum-principle argument, with one missing line supplied by the auditor; CITED for five boundary theorems, the "old SOC" inequality, the one-circle theorem and App. B; NUMERIC falsification.

PDF: proofs/manin_transfer_20261021/independent_moc4_T23_uniform_20261022.pdf (4 pages). Replay: `replay/independent_moc4_T23_uniform_20261022__check.py`.

## Deviations and limits

- Chart convention: in our replay the point C_3 = infinity contributes the factor 1 to every F_sigma. With the factor -1 every H_j changes sign, a global sign absorbed in kappa (mu). Both conventions give the same theorem. The document's formulas hold exactly with the factor 1.
- The supplement's assembly and harmonic scripts were not run. The two seven-contact inputs are documents 4 and 5 of this batch, audited separately; here they are CITED.
- The chart-sampling falsification (40000 draws) found 383 supported points but none in the SOC-failure region, so it gives little evidence there.

## A. Statements

Theorem 1 (p. 1). For every x in X with P c(x) >= 0 vanishing on T_24 = {2,3,4,6,7,8} (extra contacts allowed), Phi_{c(x)}(O) >= 0 for all O in SO(4). All spectra are included: real x (one-circle theorem), collisions (App. B), extra ties (boundary list) and t = 1 (perturbation). The contact set is literally paper T_24.

Inline claims: (1) L_O >= -M_O delta on 0 < delta <= min(g5, g9); "old SOC": delta^2 <= g10 g13 implies safety; the chart (s,t,r,u real, X complex) is exhaustive for eight-distinct points with X nonreal.

## B. Dependencies

| tie / input | cited as | resolved by own code | status |
|---|---|---|---|
| 1 | supplement T14 | {1,2,3,4,6,7} = paper T_15 (identity element) | CITED (paper T_15) |
| 10 | T19 | {2,3,4,6,8,10} = g.T_20, g = ((0,1,2,3),(0,2,1,3),0) | CITED (paper T_20) |
| 12 | T12 | {2,3,4,6,8,12} = g.T_13 | CITED (paper T_13) |
| 22 | T11 | {2,3,4,6,8,22} = g.T_12 (with transposition) | CITED (paper T_12) |
| 16 | T24 | {2,3,4,7,8,16} = g.T_25 | CITED (paper T_25) |
| 13 | independent_moc4_T23_g13 | T_24+13 = Table-5 seven-set 6 | other-doc (this batch: confirmed-modulo-citations) |
| 17 | independent_moc4_T23_g17 | T_24+17 = Table-5 seven-set 5 | other-doc (this batch: confirmed-modulo-citations) |
| other ties | f54(j) > 0 forces delta <= 0 | exact: the ties with f54(j) = 0 are exactly {1,10,12,13,16,17,22} | EXACT |
| old SOC: Phi_{f54} + 2 sqrt(Phi_{f5} Phi_{f2}) >= 0 | supplement lemma | not proved here | CITED; NUMERIC: min over SO(4) is -3.9e-30 (60 starts) |
| one-circle theorem | paper Thm C.1 | real x | CITED |
| repeated / shared values | paper App. B | collisions | CITED |
| projective covariance | elementary | normalization A = (z,0,1,x), C = (c,v,w,inf) | PROOF-READ (every F_sigma scales by the same factor) |

## C. Steps

| # | claim | status | how |
|---|---|---|---|
| 1 | generators {0,1,2,3,4,5,6,53,54,67}, rank 8, only f54 non-minor | EXACT | lib |
| 2 | canonical decomposition g = g22 f0 + g13 f2 + g12 f3 + g10 f5 + g17 f6 + g16 f53 + delta f54 + g1 f67, delta = g5-g16 = g9-g17 | EXACT | checked on every generator (the generators span the 8-dimensional contact subspace) |
| 3 | f54+f1 = f0+f6+f53, f54+f4 = f2+f3+f6+f53; delta = w54 - w1 - w4 | EXACT | so delta <= 0 gives the minor cone |
| 4 | extra ties with f54(j) = 0 are exactly 1,10,12,13,16,17,22; any other tie forces delta <= 0 | EXACT | |
| 5 | five group witnesses | EXACT | find_g |
| 6 | (1): L_O >= -M_O delta, with M_O = max(0, -Phi_{f54}(O)) | PROOF-READ | from step 2 |
| 7 | old SOC; failing region has delta^2 > g10 g13 | CITED | numeric support |
| 8 | F3-F2 = (z-c)wX; (H2,H3,H4,H6,H7,H8) = (0,1,s,r-u,r+t-u,r); H9 = r+t-uX | EXACT | sympy, in (z,x,c,v,w) |
| 9 | inverse chart returns s,t,r,u; w = ch, z-v = th(z-c), r = -d(z-1)/X | EXACT | |
| 10 | z - c = XJ/((t-uX) ell d) | EXACT | |
| 11 | every finite nonreal pole of every H_j is a zero of J | EXACT + PROOF-READ | label denominators d, ell, ell h and the factor 1/((z-c)wX) involve only real-linear polynomials in X besides J |
| 12 | H1 = s + A1/L1, H10 = (r(q+t)-su)/q - (tr/q)A1/L1 (so g10 = -(tr/q) g1), H5 -> 1 as x -> inf; H1, H5, H9, H10 have no J pole; H13, H16 do | EXACT | |
| 13 | res H13 / res H16 = z/(z-1) at zeros of J | EXACT | remainder modulo J |
| 14 | F16-F2 at c = v = z equals (z-1)(x-z)(z-w) up to sign | EXACT | |
| 15 | a failing sequence cannot approach a nonreal pole | PROOF-READ, line supplied | see gap G1 |
| 16 | liminf L_O >= 0 at real X (g9 = mu u Im X -> 0) and at X = inf (g5 -> 0) | PROOF-READ | uses (1) and delta <= min(g5, g9) |
| 17 | harmonic minimum principle on a component of {noncontact gaps > 0, L_O < 0} | PROOF-READ | L_O = -mu Im(sum_j w_j(O) H_j(X)) is harmonic off the poles |
| 18 | t = 1 by perturbing t for strict six-contact points | PROOF-READ | open conditions; the chart keeps contacts exactly |
| 19 | qA1 = 0 gives g1 = 0 identically (tie 1) | PROOF-READ | H1 is then real |

## D. Falsification

- SOC lemma: min over SO(4) of Phi_{f54} + 2 sqrt(Phi_{f5} Phi_{f2}) = -3.9e-30 (60 Nelder-Mead starts, seed 7). No violation; it appears tight (minimum 0).
- Chart sampling with our own products: 40000 draws (seed 11), 383 supported points, none in the SOC-failure region.
- Library stratum points for T_24: 10/10 found, worst min Phi = -3.3e-12.

## E. Verdict

confirmed-modulo-citations. Citations: paper cases T_12, T_13, T_15, T_20, T_25; the seven-contact documents independent_moc4_T23_g13 and _g17 (this batch); the old SOC lemma; paper Thm C.1 and App. B.

Gaps:
- G1 (presentation; p. 3, "g10 has a strictly positive finite limit because g1 > 0"). Positivity of g1 inside the region only gives limit >= 0. Supplied argument: H1 is a real Mobius function of X, nonconstant when qA1 != 0, so H1(X0) is nonreal at a nonreal pole X0. Hence g1(X0) != 0, so g1(X0) > 0, and g10 = -(tr/q) g1 has a positive limit. With that, the bound g13 < delta^2 / g10 is uniform near X0 and the residue argument goes through.
- G2 (presentation; p. 3). The norm-kernel display "(X,1) M (X,1)^T" in the companion documents must be read as the Hermitian form (X,1) M (Xbar,1)^T; here it does not matter.
