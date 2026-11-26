# Audit of blind_moc4_T44_coherent_chain_capacity_20261022

> **Second pass:** the Gram representation, for which the document cites only a companion program, is proved by the audit with exact code; the statement holds. See [the second-pass note](../docs2/blind_moc4_T44_coherent_chain_capacity_20261022__gapcheck.md).

Evidence level: EXACT for the reductions Prop 1 -> Prop 2 -> scalar test (1); NUMERIC for Props 1 and 2; GAP for the Gram identity behind Prop 1.

Auditor batch: C18-T43-T45. Source: `proofs/manin_transfer_20261021/blind_moc4_T44_coherent_chain_capacity_20261022.pdf` (2 pages). Replay: `replay/blind_moc4_T44_coherent_chain_capacity_20261022__check.py` (output in `replay/blind_moc4_T44_coherent_chain_capacity_20261022__check_output.txt`).

## Deviations and limits

- Batch budget about one hour for five documents; depth was put on the three primary proofs.
- Nothing from the supplement was executed or imported; no supplement JSON was loaded. All polynomials were rebuilt from the paper's definitions with the shared library (cubic (A.2), P, f_j of Appendix D).
- The large stored certificates (1000+-term radial factors, Bernstein tables) were NOT re-expanded symbolically. Instead the radial factor was reconstructed EXACTLY at sample points: G^2 = -Q^2 det Hess C(g - (C/Q) kappa)/4 by exact 14x14 rational determinants, interpolated along lines and square-rooted exactly; sign fixed by the d^9 coefficient = +1. Polynomial identities checked this way at finitely many random rational points are evidence (exact arithmetic, probabilistic), not proofs; they are labelled NUMERIC (exact-at-points).
- Library value check: the library cubic gives C, Q and the base normalisation G(d=1, rest 0) = 1 exactly as the documents state (constant factor 1); no disagreement with the library was found.

## A. Statements

- **Proposition 1** (p. 1): With d=1, real r,k,v, w>0, r!=k, x=(v(2r-k)-1)/(r-k), y=k(1-rv)/(w(r-k)): the cost q f2+p f1+j f53+c f12+e f8+f18+b f7+a f15 is >= 0 on SO(4) if p>=x^2+y^2, j>=r^2(1+v^2), c>=k^2(1+v^2), a,e>=v^2, b>=w^2, q>=0. A pointwise inequality on coefficients, no realizability.
- **Proposition 2** (p. 1): For b,j>0, c>=0: p f1+j f53+c f12+b f7+f18 >= 0 on SO(4) if p >= (b+c)/(b(sqrt j+sqrt c)^2) when jc<=b^2, p >= 1/(j+b) when jc>=b^2; boundary cases j=0 (chain 1/c+1/b), c=0 (SOC 1/j), b=0 (SOC 1/j).
- **scalar test (1)** (p. 2): With m=min(a,e), J=sqrt j, C=sqrt c, b>0, J+C>0: if some t in [0, sqrt(m/(1+m))] gives p >= [((2J+C)t-1)^2/(1-t^2) + c(1-Jt)^2/b]/(J+C)^2, the cost is >= 0 on SO(4).

## B. Dependencies

- [other-doc] 'already independently verified diagonal-coefficient Gram coordinates (J,C,E,A,B,Q,P,D)' and their ten coefficient relations / 160-row quaternion identity - used for: Prop 1. Statement match checked: no (CITED).

## C. Steps

| step | claim | kind | status | how checked |
|---|---|---|---|---|
| S1 | Gram of the eight listed R^3 vectors yields a nonnegative representation of Phi with the stated diagonal | certificate | GAP | the map from Gram entries to Phi_O and the ten relations are not given in the document |
| S2 | Prop 2 from Prop 1 (v=0, r=sqrt j, k=-z, w=sqrt b), derivative, optimal z, branch values | identity | EXACT | sympy |
| S3 | b=0 direct Gram | certificate | GAP | same missing framework as S1 |
| S4 | scalar test (1) from Prop 1 with the stated substitution | identity | EXACT | sympy, incl. j=J^2, c=C^2, v^2<=m |
| S5 | Props 1 and 2 true | inequality | NUMERIC | 24 + 20 random threshold costs, min Phi/sum(coeffs) >= -1e-16; example 1/2 f1+f53+f12+f7+f18 min ~ 2e-13 |

## D. Falsification

- NUMERIC: Prop 1 at threshold: 24 random (r,k,v,w), 8 restarts each, worst min Phi/sum = 7e-17
- NUMERIC: Prop 2 at threshold: 20 random (b,j,c), worst 8e-16
- NUMERIC: example cost min ~ 2e-13 (20 restarts)

Exact replays:
- EXACT: Prop 2 from Prop 1
- EXACT: derivative formula
- EXACT: both branch values
- EXACT: scalar test (1) from Prop 1
- EXACT: generator set of the cost

## E. Verdict

**gap-found**. Branch lemma for paper T_45 (sufficient coefficient bounds; not used by the primary T_45 proof): Prop 2 and test (1) follow exactly from Prop 1, and Props 1-2 pass numerical SO(4) minimisation, but Prop 1's Gram identity is not stated in the document (coordinates and ten relations only cited), so it cannot be checked here.

Not used by the primary T_45 proof (the harmonic document cites no capacity), so this gap does not affect the T_45 case.

Gaps:
- p.1, proof of Proposition 1 (and the b=0 Gram in Prop 2) [substantive]: the Gram matrix of the eight vectors gives a sum-of-squares representation of the Leibniz sum with the stated diagonal coefficients. Why: the 'diagonal-coefficient Gram coordinates' and the 'ten original coefficient relations' are cited from elsewhere and not stated; the identity is only asserted to be verified by a companion program. Believed true: True. Repair: state the bilinear map from the 3x3 Gram data to Phi_O (or cite the document that defines it) so the 160-row identity can be replayed.

Unchecked:
- the Gram construction itself
