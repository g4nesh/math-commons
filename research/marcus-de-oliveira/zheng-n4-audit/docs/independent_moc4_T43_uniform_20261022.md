# Audit of independent_moc4_T43_uniform_20261022

Evidence level: EXACT for (1), the d<=0 reduction, the core action, (4) and the Q=0 guard; NUMERIC (exact-at-points) for (5), k0, (10) and the Bernstein numerators (8); CITED for the coefficient bounds (2),(3) and the 33-table certificate.

Auditor batch: C18-T43-T45. Source: `proofs/manin_transfer_20261021/independent_moc4_T43_uniform_20261022.pdf` (5 pages). Replay: `replay/independent_moc4_T43_uniform_20261022__check.py` (output in `replay/independent_moc4_T43_uniform_20261022__check_output.txt`).

## Deviations and limits

- Batch budget about one hour for five documents; depth was put on the three primary proofs.
- Nothing from the supplement was executed or imported; no supplement JSON was loaded. All polynomials were rebuilt from the paper's definitions with the shared library (cubic (A.2), P, f_j of Appendix D).
- The large stored certificates (1000+-term radial factors, Bernstein tables) were NOT re-expanded symbolically. Instead the radial factor was reconstructed EXACTLY at sample points: G^2 = -Q^2 det Hess C(g - (C/Q) kappa)/4 by exact 14x14 rational determinants, interpolated along lines and square-rooted exactly; sign fixed by the d^9 coefficient = +1. Polynomial identities checked this way at finitely many random rational points are evidence (exact arithmetic, probabilistic), not proofs; they are labelled NUMERIC (exact-at-points).
- Library value check: the library cubic gives C, Q and the base normalisation G(d=1, rest 0) = 1 exactly as the documents state (constant factor 1); no disagreement with the library was found.

## A. Statements

- **main (Theorem 1)** (p. 1) - paper case(s) T_44: For every x in X, pi* with P c(x) >= 0 and contacts containing {1,4,5,6,7,8} (extra ties, zero coefficients, coincident spectral values included), Phi_{c(x)}(O) >= 0 on O(4). The argument uses only C(Re x)=0 and det Hess C(Re x)=0, valid for all x in X.

The contact set equals the paper's Appendix E set literally (no group element needed), and the generator list equals lib `generators(T)`. The statement is non-strict, covers all O in SO(4) (O(4) equivalently) and all spectra; it is at least as strong as (*).

## B. Dependencies

- [core-doc] 'original T24 core' (supplement label; central c, extra e, leaves a,b, SOC endpoint j) transported by g=((2,3,0,1),(1,3,0,2),0) in lib convention - used for: coefficient bounds (2),(3). Statement match checked: yes.
- [common-lemma] SOC bound cj>=1, old-H bound c>=1+1/e, four-coefficient bound K>=0, coherent bound (ab>1), parallel bound (ab<1,b>=1), closed coherent bound at ab=1 - used for: safe region. Statement match checked: no (CITED).
- [other-doc] supplement verifier: 33 Bernstein tables (8) with nonnegative integer coefficients - used for: (9). Statement match checked: no (CITED).
- [common-lemma] radial necessity C=0, det Hess C=0 on Re x - used for: G9=0 for hard realizable costs. Statement match checked: no (CITED).
- [core-doc] Theorem 1.2 (two-orthogonal law) - used for: extension to U(4) (not needed for Prop 4.3). Statement match checked: no (CITED).

## C. Steps

| step | claim | kind | status | how checked |
|---|---|---|---|---|
| S1 | generators {0,1,2,3,4,5,6,33,34,59}, nonminor 34 | identity | EXACT | lib |
| S2 | (1) with coefficients (g14,g13,g16,g0,g10,g3-g13,g15,g19) on the whole 8-dim contact nullspace | identity | EXACT | all 8 nullspace basis vectors |
| S3 | d<=0: gap identities g3=p+d, g9=b+d, g18=q+d+a, g22=j+d+a; two f-identities; replacement of -m f34; all coefficients nonnegative | identity/inequality | EXACT + PROOF-READ | sympy; x=min(a,m), y=m-x bookkeeping checked by auditor |
| S4 | profile action (59,53,55,67,54,2)->(59,5,6,0,34,1) is in G | identity | EXACT | element found |
| S5 | bounds (2),(3) sufficient | citation | CITED | other documents |
| S6 | (4): C = p{e(d-c)+d^2}+R, R >= 0 coefficientwise with d^3 coefficient 1; hence C>0 in the hard region and Q != 0 | identity/inequality | EXACT | sympy from library cubic |
| S7 | (5) Q^4 det H = -4Q^2 G9^2; G9 quadratic in c | identity | NUMERIC | exact-at-points: 6 points, perfect squares, degree <= 2 in c, G9(base)=1 |
| S8 | k0 = (b+1)(p+1)^2 P2 H4, H4 = (a+q+1)(-K)+D, D the 16 listed monomials; (7) | identity | NUMERIC | exact equality at 6 random rational points; (7) then follows since D >= 0 |
| S9 | seven chart rows cover the hard chamber (ab<1/ab>1, j vs e, leaf order, j-cuts) | case split | PROOF-READ | auditor checked each row's range and the cut values 1/j = 1+1/e, 1/e+1/b |
| S10 | 33 tables nonnegative | certificate | CITED | supplement verifier; sampled exactly at 28 chart points (S11) |
| S11 | Bernstein numerators (8) >= 0 | inequality | NUMERIC | exact values at 28 chart points, 2 per row and leaf orientation: all >= 0 |
| S12 | (9) from (8) and k0>0 | identity | PROOF-READ | standard quadratic Bernstein form |
| S13 | (10) at e=j=0 | identity | NUMERIC | exact equality at 3 points |
| S14 | Section 5 boundary layers (zero leaves, ab=1, b=1, j=e) | analytic | PROOF-READ | limits of polynomial inequalities with finite endpoints; k0 strict by (7) |

## D. Falsification

- NUMERIC: G9 degree <= 2 in c and perfect square at 6 exact points
- NUMERIC: k0 formula exact at 6 random rational points
- NUMERIC: (10) exact at 3 points
- NUMERIC: Bernstein numerators (8) >= 0 at 28 exact chart points (seed 44)
- NUMERIC: 11 realizable T_44 points (seeds 0-11): 2 with d>0, none fails every bound; min Phi ~ 8e-11

Exact replays:
- EXACT: generator inventory
- EXACT: operator identity (1) on the contact nullspace
- EXACT: gap identities and the d<=0 absorption
- EXACT: core profile action in G
- EXACT: (4) with R >= 0 and d^3 coefficient 1
- EXACT: G9(base)=1

## E. Verdict

**confirmed-modulo-citations**. Primary proof of paper T_44={1,4,5,6,7,8} (literal match): operator identity (1), the d<=0 minor-sum reduction, (4) and the C>0 guard replay exactly; k0 = (b+1)(p+1)^2 P2 H4, (10) and the Bernstein numerators check exactly at sampled points; rests on cited SOC/old-H/four-coefficient/coherent/parallel bounds and the unreplayed 33 tables.

No gap found. k0's factorization and (10) are checked only at sample points (exact arithmetic).


Unchecked:
- symbolic replay of G9 (1232 terms) and of the 33 tables
- the cited coefficient bounds themselves
