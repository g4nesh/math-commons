# Audit of independent_moc4_T47_high_e_below_core_20261023

Evidence level: EXACT (factorisation of disc_p G, orientation, full nine-variable remainder certificate), PROOF-READ (logic)

PDF: proofs/manin_transfer_20261021/independent_moc4_T47_high_e_below_core_20261023.pdf (2 pages). Auditor label C23-T47. Replay: `replay/independent_moc4_T47_high_e_below_core_20261023__check.py`, output `replay/independent_moc4_T47_high_e_below_core_20261023__check_output.txt`.

## Deviations and limits

No deviation. The remainder was checked directly in the nine-variable chart (Wj, Ws free) instead of via the document's bihomogenisation argument; this is the stronger check and covers all walls at once.

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

- **Theorem 1** (p. 1). With d=1, b>1, 0<a<1/b, p>=b, q>=0, e>=b(b+1), 0<j<=1, r=sqrt(ab), k=b(1-r)^2/(b-r)^2: if 0<c<=k, Phi>=0 and F6>=0 then G>0 (all non-strict walls included). Purely algebraic; Phi and F6 >= 0 must be supplied by other results in any application.

Contact set: {0,4,5,6,9,12} = paper T_65 (identity element; the document works directly on this face).

## B. Dependencies

- [core-doc] own derivation of G and of F6 (common setting): definitions. Statement match checked.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | disc_p G has a 123-term degree-6 irreducible factor F6; orientation [c^2]F6=(e-j)^2(q+d+b)^2 | identity | EXACT | replay: disc_p G = F6 * (260-term octic) |
| 2 | chart covers the stated domain exactly | case split | PROOF-READ | r=V/(1+V) in (0,1); j=T/(T+Wj) in (0,1]; c=kS/(S+Ws) in (0,k] |
| 3 | 160 D_G [G - (p-b)F6 - q(e-e0)(q+e-e0+r(177+573r)/(160b))Phi] has nonnegative coefficients (26004 terms at Wj=Ws=1) with 61440 T^2 S^2 | identity + coefficient-sign | EXACT | replay: remainder N/D with N, D nonnegative in all nine variables; count 26004 and coefficient 61440 reproduced |
| 4 | multipliers nonnegative on the domain | inequality | EXACT | P=p-b>=0, U=e-e0>=0, r(177+573r)/(160b) = V(177+750V)/(160b(1+V)^2) |
| 5 | strict at all walls including P=U=Q=Wj=Ws=0 | inequality | EXACT | positive monomial in A,T,V,S only, in numerator and denominator |

## D. Falsification

- The statements are exact coefficient-sign certificates; replayed exactly, so no numerical search was needed beyond the exact checks. No violation found.

## E. Verdict

**confirmed**. Algebraic exclusion G > 0 on the high-e, below-core region given Phi >= 0 and F6 >= 0; F6, the identity and the remainder certificate replayed exactly; confirmed.


Interface note: the theorem is conditional on Phi >= 0 and F6 >= 0; the document says so. Its strength in an application depends on the separate Heron (Phi) and F6 necessity results, which are outside this batch.
