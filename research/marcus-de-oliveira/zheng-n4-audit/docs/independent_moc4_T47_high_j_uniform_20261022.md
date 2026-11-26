# Audit of independent_moc4_T47_high_j_uniform_20261022

Evidence level: EXACT (Theorem 1 side identities and both endpoint certificates, 144810 and 186813 terms), PROOF-READ (Theorem 1 logic), CITED (capacities, negative-Phi exclusion, cross/stationary certificates), GAP (Theorem 2 certificates not given), NUMERIC (Theorem 2 claims)

PDF: proofs/manin_transfer_20261021/independent_moc4_T47_high_j_uniform_20261022.pdf (3 pages). Auditor label C23-T47. Replay: `replay/independent_moc4_T47_high_j_uniform_20261022__check.py`, output `replay/independent_moc4_T47_high_j_uniform_20261022__check_output.txt`.

## Deviations and limits

Theorem 2's middle-range chart is not displayed and the multipliers of its 42G identity are not given in the document, so neither could be replayed from the PDF; those steps were tested numerically only. The capacities old-H, coherent-J stationary and "old stated core" c >= 1/e + k are CITED (the audit foundation part2b lists the old-core specialisation as UNCHECKED).

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

- **Theorem 1** (p. 1). Every spectrally realizable T_65 source (d=1) with q>=0, a,p,c,e>0, 0<ab<1, 1<j<b, pc>=a^2 satisfies Phi_c(O)>=0 on SO(4).
- **Theorem 2** (p. 2). Same with pc>=a^2 replaced by e>=j (1<j<b): every realizable source satisfies the support inequality.

Contact set: {0,4,5,6,9,12} = paper T_65 (identity element; the document works directly on this face).

## B. Dependencies

- [core-doc] old-H atom j >= 1+1/p (audit foundation so4-atoms): capacity. CITED, not audited here.
- [core-doc] coherent-J stationary family j >= 1/(a+c) when pc >= a^2: capacity. CITED, not audited here.
- [core-doc] old stated core capacity c >= 1/e + k (foundation part2b: UNCHECKED): capacity. CITED, not audited here.
- [branch] independent_moc4_T47_high_j_small_b_20261022 (this batch, confirmed-modulo-citations): conic identity, eleven signs, exclusion of c<=k and b<=b*. Statement match checked.
- [core-doc] merged coherent-C capacity (independent_moc4_T47_merged_C_high_j_uniform Lemma 1; foundation part2b G4): Theorem 2 capacity. Statement match checked.
- [branch] IIIb_low_j_cross_exact and IIIb_low_j_stationary_exact verifiers (not audited): Theorem 2 failure region. CITED, not audited here.
- [branch] negative-Phi (negative-Heron) exclusion theorem (not in this batch): Theorem 2, Phi<0. CITED, not audited here.
- [common-lemma] radial necessity (audit item 5): G=0. Statement match checked.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| T1.1 | capacities cover j>=1+1/p, j>=1/(a+c) (pc>=a^2), c>=1/e+k | citation | CITED | see dependencies |
| T1.2 | hard part = chamber (1) of high_j_small_b; that theorem excludes c<=k and b<=b* | citation | EXACT | same chamber; replayed in this batch |
| T1.3 | 1/b-a-k = r(1-r){2b^2-2b(r+1)+r^2+r}/(b(b-r)^2), b* the larger root; hence 1/j-a-k>0 on the chart | identity | EXACT | sympy, and chart numerator/denominator nonnegative |
| T1.4 | root => A0<0, F concave; old core failure => u<U | inequality | PROOF-READ | as Lemma 1 of high_j_small_b |
| T1.5 | E>0 on p<=b chart (p=bP/(P+W)), 144810 nonnegative coefficients | coefficient-sign | EXACT | replay mode le: 144810 terms reproduced, strict monomial free of W |
| T1.6 | E>0 on b<=p<H chart (p=(b+PH)/(1+P)), 186813 nonnegative coefficients | coefficient-sign | EXACT | replay mode ge: 186813 terms reproduced, strict monomial free of P |
| T2.1 | merged coherent-C capacity with L=a+b | citation | CITED | merged_C Lemma 1 (foundation G4) |
| T2.2 | failure region handled by cross and stationary certificates | citation | CITED | other documents, not audited |
| T2.3 | j<=e<=L^2/j: cleared G has 1232529 positive coefficients; endpoint 142167 | coefficient-sign | GAP | chart not displayed in the document; NUMERIC only |
| T2.4 | e>=L^2/j, pc<=a^2, Phi>=0: 42 G~ = R + M_pc num(a^2-pc) + M_Phi num(Phi) | identity | GAP | multipliers not given; NUMERIC only |
| T2.5 | pc>=a^2 -> Theorem 1; Phi<0 -> negative-Phi exclusion | citation | CITED |  |

## D. Falsification

- Thm 2 middle range (j<=e<=L^2/j, c below the merged-C crossing bound, p<1/(j-1)): 20000 samples, seed 11, min G = 107.9 > 0
- Thm 2 high range (e>=L^2/j, c<1/(L+j)): 20000 samples, seed 11, min G = 104.3 > 0 (144.1 on the pc<=a^2, Phi>=0 subset)
- Thm 1 endpoint E (b>b*, 1<j<b, k<c<1/j-a, 0<p<1/(j-1)): 20000 samples, seed 12, min E = 4.77 > 0 (consistent with the exact certificates)

## E. Verdict

**gap-found**. Theorem 1 (high-j large-pc) reduces exactly to the endpoint inequality E>0, replayed exactly; Theorem 2 rests on two certificates whose chart or multipliers are not given in the document (numeric evidence only) plus uncited-here companion results: gap-found (substantive, believed true).

Gaps:
- (presentation) p.3, Thm 2 proof, "For j <= e <= L^2/j, the original cleared G has 1232529 positive coefficients". Claim: G>0 on the merged-C failure region j<=e<=L^2/j, 1<j<b. Why: The chart (how e and c are parametrised, which bounds on p) is not stated, so the certificate cannot be replayed from the document; it is a pointer to a companion verifier. Believed true: yes (numeric evidence). Repair: state the chart (presumably j=1+(b-1)theta(T), e=jR^2, c=theta(S)(L+e)/(L j (1+R)^2), p=theta(P)/(j-1) as in merged_C Thm 2) and replay.
- (substantive) p.3, Thm 2 proof, identity 42 G~ = R + M_pc num(a^2-pc) + M_Phi num(Phi). Claim: G>0 when e>=L^2/j, c<1/(L+j), pc<=a^2, Phi>=0. Why: The multipliers M_pc and M_Phi are not printed; the identity cannot be checked from the document. Believed true: yes (numeric evidence). Repair: print the multipliers (5 and 6 monomials) so the identity can be replayed.

Unchecked: cited capacities (old-H, coherent-J stationary, old core); cross and stationary certificates; negative-Phi exclusion.
