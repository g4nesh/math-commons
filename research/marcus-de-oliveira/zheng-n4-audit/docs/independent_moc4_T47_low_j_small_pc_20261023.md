# Audit of independent_moc4_T47_low_j_small_pc_20261023

Evidence level: EXACT (chart and full positivity certificate, all walls), PROOF-READ (coverage), CITED (radial necessity, audit item (5))

PDF: proofs/manin_transfer_20261021/independent_moc4_T47_low_j_small_pc_20261023.pdf (2 pages). Auditor label C23-T47. Replay: `replay/independent_moc4_T47_low_j_small_pc_20261023__check.py`, output `replay/independent_moc4_T47_low_j_small_pc_20261023__check_output.txt`.

## Deviations and limits

No deviation. Positivity was checked directly in the nine-variable chart (Wj, Ws free).

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

- **Theorem 1** (p. 1). No spectrally realizable T_65 source (d=1) has b>1, 0<a<1/b, p>=b, q>=0, e>0, 0<j<=1, j<=eb/(e+b), 0<c<=a^2/p. Hence the support inequality holds vacuously in that sector. Zero layers c=0, j=0 excluded (stated).

Contact set: {0,4,5,6,9,12} = paper T_65 (identity element; the document works directly on this face).

## B. Dependencies

- [common-lemma] radial necessity (audit item 5) + forced-height derivation: G=0. Statement match checked.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | height constant C<0 (the -d^3 term) | identity | EXACT | common setting |
| 2 | chart covers the sector exactly, walls j=1, pc=a^2, p=b, j=eb/(e+b) at Wj=Ws=P=U=0 | case split | PROOF-READ | e>=bj/(b-j) <=> j<=eb/(e+b); a=V^2/(b(1+V)^2) in (0,1/b) |
| 3 | e = En/Ed with Ed=bWj+AT, En=bT+U Ed | identity | EXACT | replay |
| 4 | cleared G has 137198 positive integer coefficients; T^4S^2 coefficient 18 | coefficient-sign | EXACT | replay, counts match |
| 5 | strict at all five simultaneous walls | inequality | EXACT | positive monomial free of P,U,Q,Wj,Ws; denominators positive likewise |
| 6 | G>0 contradicts G=0 | citation | CITED | radial necessity |

## D. Falsification

- The statements are exact coefficient-sign certificates; replayed exactly, so no numerical search was needed beyond the exact checks. No violation found.

## E. Verdict

**confirmed-modulo-citations**. Exclusion of the low-j minus, p>=b, small-pc sector by G>0; the 137198-term certificate replayed exactly from the auditor's own G; confirmed modulo radial necessity.

