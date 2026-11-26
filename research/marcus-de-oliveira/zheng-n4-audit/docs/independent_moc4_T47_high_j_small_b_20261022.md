# Audit of independent_moc4_T47_high_j_small_b_20261022

Evidence level: EXACT (conic identity and all eleven sign tables), PROOF-READ (Lemma 1 and the logic), CITED (radial necessity, audit item (5))

PDF: proofs/manin_transfer_20261021/independent_moc4_T47_high_j_small_b_20261022.pdf (3 pages). Auditor label C23-T47. Replay: `replay/independent_moc4_T47_high_j_small_b_20261022__check.py`, output `replay/independent_moc4_T47_high_j_small_b_20261022__check_output.txt`.

## Deviations and limits

No deviation. The square parts were checked with multiplier 1 (cleared polynomial = square + nonnegative remainder); remainder term counts 2195, 1360, 3234, 3024 and sector counts 20490, 25008, 21798, 4032, 6057, 6057 were reproduced exactly.

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

- **Theorem 1** (p. 1). In the chamber j>1, b>j, 0<ab<1, p>0, p(j-1)<1, 0<c<1/j-a (d=1), no spectrally realizable source has 1<b<=b* or 0<c<=k (b*=(1+r+sqrt(1-r^2))/2, k=b(1-r)^2/(b-r)^2, r=sqrt(ab)), for any q,e>=0. Mechanism: G>0 there, contradicting G=0.
- **Lemma 1** (p. 2). In the chamber, with F(u)=A0u^2+B0u+K, any root of G with q,e>=0 has F(u)<=0 and A0<0; if moreover c>k and e<1/(c-k) (old core failure), no root exists when A0+B0(b+1)(c-k)+K(b+1)^2(c-k)^2>=0.

Contact set: {0,4,5,6,9,12} = paper T_65 (identity element; the document works directly on this face).

## B. Dependencies

- [common-lemma] radial necessity (audit item (5)) + auditor's forced-height derivation det Hess = -4d^2G^2/Q^2: G = 0 at realizable sources. Statement match checked.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | height constant C<0, so Q != 0 and forced-height determinant is -4G^2/Q^2 | identity | EXACT | common setting |
| 2 | conic identity G = e^2(A0+A1q+A2q^2)+e(q+b+1)(B0+B1q)+(q+b+1)^2K | identity | EXACT | replay (exact division, reassembly) |
| 3 | A1, A2, B1 > 0 on the broad chart (square part + nonnegative remainder) | coefficient-sign | EXACT | replay; square parts P^2S^2(V^2-1)^2, PS^2(V^2-1)^2, P^2S^2(1+V)(V^2-1)^2 with multiplier 1 |
| 4 | B0-A0 > 0 and K > 0 on the broad chart | coefficient-sign | EXACT | replay (3024 terms each) |
| 5 | broad chart covers chamber (1) | case split | PROOF-READ | 1/j-a = (j+A+AV)/(...) > 0; all ranges bijective |
| 6 | small-b sector: A0, B0, K > 0 including b=b* | coefficient-sign | EXACT | replay; b*-1 = V/((1+V)^2+V^2) > 0 |
| 7 | below-core sector: A0, B0, K > 0 including c=k | coefficient-sign | EXACT | replay; chart covers 0<c<=k, a superset of the chamber part |
| 8 | all six coefficients positive and u,w>=0 => G>0 | inequality | PROOF-READ |  |
| 9 | Lemma 1 (concavity argument) | inequality | PROOF-READ | F concave, F(0)=K>0, F(U)>=0 => F>0 on [0,U); u<=e/(b+1)<U |

## D. Falsification

- The statements are exact coefficient-sign certificates; replayed exactly, so no numerical search was needed beyond the exact checks. No violation found.

## E. Verdict

**confirmed-modulo-citations**. Conic decomposition of G and exclusion of the small-b and below-core high-j sectors; every identity and sign table replayed exactly from the auditor's own G; confirmed modulo radial necessity (audit item 5).

