# Audit of blind_moc4_T44_conic_branch_exclusion_20261022

> **Second pass:** the polynomial identities, whose polynomials the document does not print, were rebuilt from the paper's definitions and verified exactly; the statement holds. See [the second-pass note](../docs2/blind_moc4_T44_conic_branch_exclusion_20261022__gapcheck.md).

Evidence level: EXACT for the cubic guard and the explicit algebra of Lemma 2; NUMERIC (exact-at-points) for the radial substitution identity; GAP for the square-class identities (2), whose polynomials are not printed; UNCHECKED for (3); CITED for the inertia theorem and secant rank bound.

Auditor batch: C18-T43-T45. Source: `proofs/manin_transfer_20261021/blind_moc4_T44_conic_branch_exclusion_20261022.pdf` (3 pages). Replay: `replay/blind_moc4_T44_conic_branch_exclusion_20261022__check.py` (output in `replay/blind_moc4_T44_conic_branch_exclusion_20261022__check_output.txt`).

## Deviations and limits

- Batch budget about one hour for five documents; depth was put on the three primary proofs.
- Nothing from the supplement was executed or imported; no supplement JSON was loaded. All polynomials were rebuilt from the paper's definitions with the shared library (cubic (A.2), P, f_j of Appendix D).
- The large stored certificates (1000+-term radial factors, Bernstein tables) were NOT re-expanded symbolically. Instead the radial factor was reconstructed EXACTLY at sample points: G^2 = -Q^2 det Hess C(g - (C/Q) kappa)/4 by exact 14x14 rational determinants, interpolated along lines and square-rooted exactly; sign fixed by the d^9 coefficient = +1. Polynomial identities checked this way at finitely many random rational points are evidence (exact arithmetic, probabilistic), not proofs; they are labelled NUMERIC (exact-at-points).
- Library value check: the library cubic gives C, Q and the base normalisation G(d=1, rest 0) = 1 exactly as the documents state (constant factor 1); no disagreement with the library was found.

## A. Statements

- **thm (Theorem 1)** (p. 1): For the T_45 cost q f2+p f1+j f53+c f12+e f8+d f18+b f7+a f15 with q,p,j,c,e,b,d>0, a>=0, pj<d^2, j<d+d^2/(a+e), and Phi<0: there is no x in X whose restored real point is c + h kappa. Hence realizable costs outside the SOC/Hsum capacities satisfy Phi>=0. A conditional exclusion, not a case of Prop 4.3.

## B. Dependencies

- [common-lemma] global Hessian inertia theorem: order-12 principal minors of Hess C(Re x) are >= 0 - used for: final contradiction. Statement match checked: no (CITED).
- [external] HMSV set equality / secant rank bound rank Hess <= 12 - used for: rank exactly 12. Statement match checked: no (CITED).
- [other-doc] square-class verifier (polynomials Psi, S, T, Omega, U) and earlier T44 height-pencil verifier - used for: (2), provenance of P. Statement match checked: no (CITED).

## C. Steps

| step | claim | kind | status | how checked |
|---|---|---|---|---|
| S1 | C>0 on the chamber: (a+e)C - ae{d(a+e)+d^2-j(a+e)} >= 0 coefficientwise with (a+e)d^3 | identity | EXACT | sympy from library cubic |
| S2 | (1): det H = -4G^2/Q^2, G quadratic in q | identity | NUMERIC | implicit in S7 sampling |
| S3 | (2): disc = e^2 Phi Psi; P(e^2 Phi S^2) = (U+Omega)^2, P(Psi T^2) = (U-Omega)^2 mod G | certificate | GAP | Psi, S, T, Omega, U not given in the document |
| S4 | conic form and (3) det H = 2jXY Phi^2 | identity | UNCHECKED | needs G as bivariate polynomial; time |
| S5 | Lemma 2: Y=0 impossible (DR identity, D^2 Phi = (D^2-X')^2) | identity | EXACT | sympy |
| S6 | X=0 case: Psi substitution | identity | GAP | Psi not printed (part of S3) |
| S7 | W^2E1^2 G|_{a=pc/W,e=-E0/E1} = d(d+b)Z T3 (qL1+L0)^2; L1+WE1 identity; signs Z,T3,L1,L0<0 | identity | EXACT + NUMERIC | L1+WE1 symbolic; substitution identity exact at 3 points (ratio 1); signs sampled |
| S8 | rank >= 12 from d^2/dt^2 det < 0; Schur congruence; all maximal principal minors <= 0, one < 0 | analytic | PROOF-READ | standard; argument re-derived (compound/Cauchy-Binet sign) |
| S9 | density of S Omega != 0 on the real hypersurface | analytic | PROOF-READ | valid given gcd(A,B,C0)=1 and S, Omega q-free |

## D. Falsification

- NUMERIC: radial substitution identity exact at 3 rational points (G normalisation matches, ratio 1)
- NUMERIC: signs Z,T3,L1+WE1 < 0 under pj<d^2 on 2000 samples
- NUMERIC: 11 realizable T_45 points: none in the chamber (test is weak)

Exact replays:
- EXACT: cubic guard (a+e)C - ae{...} >= 0
- EXACT: DR identity and Y=0 square identity
- EXACT: a=0: Y=-pcbd
- EXACT: L1+WE1 identity

## E. Verdict

**gap-found**. Branch lemma for paper T_45 (no realizable lift in an Hsum-hard chamber with Phi<0; not used by the primary T_45 proof): the C>0 guard, the Y=0 and X=0 algebra and the radial substitution identity check, the rank/inertia logic is sound, but (2) involves unprinted polynomials and (3) was not replayed.

Not used by the primary T_45 proof.

Gaps:
- p.1, display (2) and p.2 (X=0 octic identity) [substantive]: three square-class identities modulo G and the identity W^2 Psi|_{a=pc/W} = c^2(E1 e+E0)^2. Why: the polynomials Psi (137 terms), S, T, Omega, U are not printed, so the identities cannot be checked from the document; they are asserted via a companion verifier. Believed true: True. Repair: publish the polynomials or derive the sign of the 12x12 minor directly.

Unchecked:
- (3) conic determinant
- G itself (1188 terms) symbolically
