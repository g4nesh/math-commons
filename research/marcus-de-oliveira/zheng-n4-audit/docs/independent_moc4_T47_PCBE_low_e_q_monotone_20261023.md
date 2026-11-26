# Audit of independent_moc4_T47_PCBE_low_e_q_monotone_20261023

Evidence level: EXACT (all certificates re-derived from the auditor's own G), PROOF-READ (logic)

PDF: proofs/manin_transfer_20261021/independent_moc4_T47_PCBE_low_e_q_monotone_20261023.pdf (2 pages). Auditor label C23-T47. Replay: `replay/independent_moc4_T47_PCBE_low_e_q_monotone_20261023__check.py`, output `replay/independent_moc4_T47_PCBE_low_e_q_monotone_20261023__check_output.txt`.

## Deviations and limits

No deviation. G was rebuilt from the paper's cubic C (see common setting); the document's polynomial was never read from the supplement. The term counts 40647/41139/38301 and 38304 were reproduced exactly with the canonical multidegree clearing.

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

- **Lemma 1** (p. 1). On the rational chart b=1+A, p=b+P, r=V/(1+V), a=r^2/b, s=S/(1+S), e=b(b+1)U/(1+U), c=s^2(1/e+1/(b+1)), j=J* T/(1+T), J*=(b+1)(1-s)^2/(b+2s+(b+1)s^2/e), with A,T,U,V,S>0, P>=0, d=1: if Phi>=0 then G is strictly increasing in q>=0. Purely algebraic (no realizability claim, no claim about G(0)).

Contact set: {0,4,5,6,9,12} = paper T_65 (identity element; the document works directly on this face).

## B. Dependencies

- [core-doc] own derivation of G (common setting); radial necessity is not used by the lemma itself: definition of G. Statement match checked.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | canonical clearing denominator is positive and independent of q | identity | EXACT | replay: nonnegative coefficients, free of Q, positive monomial free of P |
| 2 | cleared G = g0 + Q g1 + Q^2 g2 with 40647, 41139, 38301 terms | identity | EXACT | replay, counts match |
| 3 | every coefficient of g1 is a positive integer; survives P=0 | coefficient-sign | EXACT | replay |
| 4 | g2 = R + 2AU^2V^2 N_Phi, R with 38304 positive integer coefficients, D_Phi positive | identity + coefficient-sign | EXACT | replay (N_Phi = reduced numerator of Phi) |
| 5 | Phi>=0 => dG/dQ = g1 + 2Q g2 > 0 for Q>=0 | inequality | PROOF-READ | immediate |
| 6 | chart origin: failure of (2) at L=b, low e <=> s in (0,1), j<J* | identity | EXACT | sympy: (c+1)/(1-s)^2 - 1/(b+1) = 1/J*; B=1/e+1/(b+1) when e<b(b+1) |

## D. Falsification

- The statements are exact coefficient-sign certificates; replayed exactly, so no numerical search was needed beyond the exact checks. No violation found.

## E. Verdict

**confirmed**. Strict q-monotonicity of the radial octic G on the low-e PCBE-failure chart when Phi >= 0; every certificate replayed exactly from the auditor's own G; confirmed.


Remark for users of this lemma: it is stated on the open chart (A,T,U,V,S > 0); the endpoints e = b(b+1), s = 1 (c = B) and j = J* are not covered, as the document itself says for e = b(b+1). Whoever cites it must handle those walls elsewhere.
