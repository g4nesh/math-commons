# Audit of independent_moc4_T47_PCBE_energy_20261022

Evidence level: EXACT (scalar identities), CITED (Gram identity: audit foundation so4-atoms-part2b G4), NUMERIC (threshold attack)

PDF: proofs/manin_transfer_20261021/independent_moc4_T47_PCBE_energy_20261022.pdf (2 pages). Auditor label C23-T47. Replay: `replay/independent_moc4_T47_PCBE_energy_20261022__check.py`, output `replay/independent_moc4_T47_PCBE_energy_20261022__check_output.txt`.

## Deviations and limits

No deviation. The Gram matrix identity against the 160 quaternion coefficients was not re-derived here: it is proved in the audit foundation so4-atoms-part2b.md (G4, "PCBE Thm 1 ... PROVED IN THE DOCUMENT AND CHECKED"), which identified the documents' spectral Gram coordinates. This report checks that the foundation's statement is the document's statement (same face, same coefficient roles, same region) and re-derives the scalar steps.

- Supplement scripts and JSON were not executed or read; every polynomial was rebuilt from the paper (cubic C, Appendix D generators) through the audit library.
- Not the primary proof of a paper case; the catalogue lists it as used for paper T_47, T_48, T_65 (T_65 is the historical T47 face; uses in T_47/T_48 were not traced here).

## Common setting (applies to every T47 document of batch C23-T47)

Historical "T47" is the paper's case T_65: contact set {0,4,5,6,9,12} (0-based lexicographic indices). The supplement's facet numbers coincide with the paper's Appendix D numbering here: lib `generators(T_65)` = {0,1,2,5,7,9,35,45} = the document's facets (q,p,j,c,e,d,b,a) = (f2,f9,f1,f45,f5,f35,f0,f7). Only f35 = h_{4312,-} is non-minor; the eight vectors are linearly independent (rank 8), so a source's coefficients are unique. (EXACT, lib.)

The "original signed radial octic G" was re-derived by the auditor from the paper's cubic C (A.2), not taken from the supplement. With c the source cost and kappa = e_13:
- C(c + h kappa) = C(c) + h Q(c), Q = m0 m6 - m1 m5 (EXACT), and C(c) = -d (qp+qe+qd+qa+pe+pd+pb+ed+d^2+db+da+ba) (EXACT). This is the documents' height constant "C < 0" (they normalise d = 1). For d > 0 and the other coefficients >= 0, C(c) < 0, so a realizable source (which has C(u) = 0 with u = c + h kappa, radial necessity) has Q(c) != 0 and forced height h = -C(c)/Q(c).
- det Hess C(c + h kappa) has degree 3 in h, and at the forced height it equals -4 d^2 G^2 / Q^2 with G a 965-term integer polynomial, homogeneous of degree 8, multidegree (2,2,2,2,2,8,3,3) in (q,p,j,c,e,d,b,a) (EXACT, own 14x14 symbolic determinant and factorisation). This matches every description in the documents (degree eight, multidegree (2,2,2,2,2,8,3,3,0), "forced-height determinant -4G^2/Q^2"). The sign of G is fixed by the documents' positivity claims; every positivity certificate below came out with this one orientation, which is a strong consistency check that the auditor's G is the documents' G.
- Radial necessity (audit item (5): C(u) = 0 and det Hess C(u) = 0 for u = Re x, x in X) then gives G = 0 at every realizable source with d > 0 (d = 0 means all generators are minor, where Phi >= 0 is trivial).

Charts are replayed by substituting the displayed rational chart into G (and the other displayed polynomials) and clearing with the product of each variable's denominator to the power of its degree; signs are then read off coefficient by coefficient. A positive-coefficient numerator over a positive-coefficient denominator, with a positive monomial in the strictly positive chart variables only, gives strict positivity on the closed chart, walls included.


## A. Statements

- **Theorem 1** (p. 1). For the cost qf2+pf9+jf1+cf45+ef5+df35+bf0+af7 with d=1, p,b,j,e>0, all coefficients >=0, and 0<L<=min(p,b): if c >= min_x max{K(1-x)^2-1, Bx^2}, K=1/j+1/(L+1), B=max{1/e+1/(L+1),1/L} (equivalently sqrt((c+1)/K)+sqrt(c/B) >= 1), then Phi_c(O) >= 0 for every O in SO(4). Pointwise coefficient statement: no spectral hypothesis, so it covers all spectra. General d>0 by scaling.

Contact set: {0,4,5,6,9,12} = paper T_65 (identity element; the document works directly on this face).

## B. Dependencies

- [core-doc] audit foundation so4-atoms-part2b.md G4 (PCBE Gram identity, Schur complement, completion, closure): Gram matrix M satisfies the ten spectral relations / 160 coefficients; PSD completion argument. Statement match checked.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | Gram matrix M with diagonal (0,L,j,c,e,1,L,0) and displayed off-diagonals represents the cost | identity | CITED | foundation so4-atoms-part2b G4 (EXACT there, in the identified spectral coordinates) |
| 2 | Schur complement over D,B gives the displayed 3x3 matrix on J,E,C | identity | CITED | foundation G4 (EXACT there) |
| 3 | free JE entry w makes Schur complement of C0 diagonal (path-graph completion) | analytic | PROOF-READ | foundation G4 and re-read: AJ, AE as displayed; completion valid iff both 2x2 blocks PSD |
| 4 | ratio identity (3) | identity | EXACT | sympy in replay |
| 5 | three expressions simplify to K(1-x)^2-1 and Bx^2 in (1) | identity | EXACT | sympy in replay |
| 6 | x* is the crossing point when K>1 | identity | EXACT | sympy in replay |
| 7 | equality cases by closure of a fixed-diagonal slice; clamping x to [0,1]; L=min(p,b) optimal | analytic | PROOF-READ | bounded fixed-diagonal PSD slices are compact; K,B decreasing in L |

## D. Falsification

- 300 random (p,b,j,e) ~ lognormal, c at the threshold of (1), q=a=0, lib min_phi_over_SO4 with 4 restarts each, seed 7: worst min Phi = -4.2e-16 (no violation)

## E. Verdict

**confirmed-modulo-citations**. PCBE Gram capacity: sufficient coefficient bound for Phi >= 0 on the T_65 face; Gram identity is the one proved in the audit foundation, scalar steps re-derived exactly; confirmed modulo that citation.


The document is explicitly a sufficient coefficient bound and claims nothing about realizability; its use as the "PCBE network" in docs 2 and 5 is with L = b (p >= b), which is inside the stated range 0 < L <= min(p,b).

Unchecked: 160-coefficient Gram identity not re-derived in this report (cited from the foundation).
