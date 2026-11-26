# Audit of independent_moc4_T47_merged_C_high_j_uniform_20261022

Evidence level: EXACT (all three positive-coefficient tables of ranges 1 and 2, including the je=A0^2 endpoint), CITED (Lemma 1 via foundation part2b G4; large-pc theorem; negative-Heron), NUMERIC (range 3), GAP (range 3 multipliers absent)

PDF: proofs/manin_transfer_20261021/independent_moc4_T47_merged_C_high_j_uniform_20261022.pdf (3 pages). Auditor label C23-T47. Replay: `replay/independent_moc4_T47_merged_C_high_j_uniform_20261022__check.py`, output `replay/independent_moc4_T47_merged_C_high_j_uniform_20261022__check_output.txt`.

## Deviations and limits

The 1810026-term certificate of range b<=j<A0, je<=A0^2 was replayed exactly (mode r2, long job): 1810026 positive coefficients, as stated. The endpoint je=A0^2 was replayed in the auditor's own reduced clearing (80874 positive terms, not the document's 201114 homogenised terms). Range 3 multipliers are not in the document.

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

- **Lemma 1** (p. 1). For d>0, A0=a+b: c >= d^2(A0+e)/(A0(sqrt j+sqrt e)^2) if je<=A0^2, c >= d^2/(j+A0) if je>=A0^2, and (e>=j) c >= d^2/(4j)+d^2/(4A0), are each sufficient for Phi_c >= 0 on SO(4) (pointwise, all spectra).
- **Theorem 2** (p. 2). With d=1, b>1, 0<a<1/b, p,c>0, q>=0, e>=j>=b: every spectrally realizable source satisfies the support inequality (including j=b, q=0, e=j, je=(a+b)^2).

Contact set: {0,4,5,6,9,12} = paper T_65 (identity element; the document works directly on this face).

## B. Dependencies

- [core-doc] audit foundation so4-atoms-part2b G4 (merged-C Lemma 1: six Grams and the allocation): Lemma 1. Statement match checked.
- [core-doc] old-H atom: j>=1+1/p. CITED, not audited here.
- [other-doc] independent_moc4_T47_IIIb_plus_large_pc_20261022 (high-j large-pc theorem for j>=b; not audited here): range 3, pc>=a^2. CITED, not audited here.
- [branch] negative-Heron exclusion (Phi<0): range 3. CITED, not audited here.
- [common-lemma] radial necessity (audit item 5): G=0. Statement match checked.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| L1 | six Grams and proportional allocation | identity | CITED | foundation G4 (PROVED IN THE DOCUMENT AND CHECKED) |
| T.1 | case split j>=A0 / b<=j<A0 with je<=A0^2 / je>=A0^2 is exhaustive on e>=j>=b | case split | PROOF-READ | for j<A0, A0^2/j>j so range 3 has e>=j |
| T.2 | range j>=A0: cleared G has 109050 positive coefficients, strict at T=U=Q=0 | coefficient-sign | EXACT | replay, count matches |
| T.3 | range b<=j<A0, je<=A0^2: 1810026 positive coefficients | coefficient-sign | EXACT | replay mode r2: 1810026 positive coefficients, strict monomial free of T,U,Q |
| T.4 | endpoint R=A0/j positive | coefficient-sign | EXACT | own reduced clearing: 80874 positive terms |
| T.5 | range 3: pc>=a^2 -> large-pc theorem; Phi<0 -> negative-Heron | citation | CITED |  |
| T.6 | range 3, pc<=a^2, Phi>=0: 9990288000 G~ = R0 + Mpc Fpc + MPhi FPhi | identity | GAP | multipliers "embedded in the standalone verifier", not in the document; NUMERIC only |

## D. Falsification

- range 2 interior (b<=j<A0, 1<=R<A0/j, c below crossing bound, p<1/(j-1)): 20000 samples, seed 11, min G = 89.8 > 0
- range 3 (e>=A0^2/j, c<1/(j+A0), p<1/(j-1)): 20000 samples, seed 11, min G = 86.7 overall and 109.7 on the 11790-sample pc<=a^2, Phi>=0 subset (no violation)

## E. Verdict

**gap-found**. Merged coherent-C capacity (cited, proved in the foundation) and the j >= b high-j sector; the two positive-coefficient range certificates (109050 and 1810026 terms) replayed exactly, the third range rests on an identity whose multipliers are not in the document: gap-found (substantive, believed true).

Gaps:
- (substantive) p.2-3, range b<=j<A0, je>=A0^2, identity 9990288000 G~ = R0 + Mpc Fpc + MPhi FPhi. Claim: G>0 when pc<=a^2 and Phi>=0 in range 3. Why: The twelve- and fourteen-monomial multipliers are not printed, so the certificate cannot be checked from the document. Believed true: yes (numeric evidence). Repair: print the multipliers and replay the identity.

Unchecked: large-pc theorem for j>=b (other document); negative-Heron exclusion; range-3 identity (multipliers not in the document).
