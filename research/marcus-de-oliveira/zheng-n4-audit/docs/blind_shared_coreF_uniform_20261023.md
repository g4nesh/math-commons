# Audit: blind_shared_coreF_uniform_20261023 ("core F", seven-contact set F)

Evidence level: EXACT for the kernel basis, the relation (2), the d < 0 decomposition, the SOC group actions, C and Q, the determinant identity (7), G and all four positivity tables (own code from the paper's definitions); PROOF-READ for the realizability necessity and the chart coverage; CITED for the sharp SOC and old-H reference inequalities; NUMERIC for falsification.

PDF: proofs/manin_transfer_20261021/blind_shared_coreF_uniform_20261023.pdf (4 pages). Replay: `replay/blind_shared_coreF_uniform_20261023__check.py` (output in `__check_output.txt`, 39 PASS, 0 FAIL).

## Deviations and limits

- (7) was replayed at d = 1 and re-homogenized, as in core E: both sides are homogeneous of degree 22. G was derived independently as an exact square root and equals the document's 441 coefficients, which were read from its JSON as data for comparison only.
- The old-H reference inequality is CITED and its source was not audited. The document identifies it only as "one of the three explicit pairs in the proved scalar-compensation/principal-minor result", with the pair index in its checker output. The data's "mother minor indices [7, 6]" are not paper f-indices: no element of G maps paper f7, f6 to f11, f12. They index a local nine-minor list (a sibling document, core B, names that list). Under the action the data names (row (0,1,3,2), column (3,1,0,2), no transpose), the mother triple in paper numbering is (f12, f29; f80 = h_{1324,+}). This was checked exactly. I could not locate the statement "t f12 + c f29 + d f80 safe if t(c - d) >= d^2" (or its variant with roles swapped) in a document read for this audit.
- The unitary remark is not needed for (*) and was not checked.

## A. Statements

F = {0,1,2,7,9,10,14}. This is literally seven-set no. 7 of the paper's Table 5. Its generators are generators(F) = {f0,f1,f2,f11,f12,f15,f21,f26}: eight rays (a nonsimplicial cone), all minor vectors except f26 = h_{4132,-}.

- **Theorem 1 (p. 1).** For every x in X (all spectra, including degenerate ones) and every reference pi*: if P c(x) >= 0 and (P c(x))_s = 0 on F (non-strict, so extra contacts are allowed), then Phi_{c(x)}(O) >= 0 for every O in SO(4). This is (*) for F. The proof is pointwise and algebraic, so all spectra are covered.
- **Eq. (1)-(2).** c(x) = q f0 + p f2 + v f11 + b f12 + c f15 + w f21 + d f26 uniquely. q, p, v, b, c, w are the gaps at [23,15,17,8,4,16], so they are >= 0, while d is any real number. Also f26 + f1 = f0 + f11 + f21.
- **d <= 0.** The gaps at 20, 11 and 3 are q+d, v+d and w+d. So q, v, w >= -d, and c(x) is a nonnegative combination of minor vectors.
- **(3)-(4), pointwise sufficient tests.**
  - (3) p(b+c) >= 1, b(p+c) >= 1 and c(p+b) >= 1 (with d = 1), from three SOC atoms with leaves in {f2,f12,f15} and non-minor f26.
  - (4) v(b-1) >= 1 (old-H on (f11, f12; f26)).
- **(5)-(7).**
  - C(c + z kappa) = Cst + Q z with Cst = (q+d)(v+d)(w+d) + pbc and Q = pb + pc + bc + v(d-b) + d^2.
  - Q^4 det H(w(-C/Q)) = -4 (v+d)^2 Q^2 G^2, where G has degree 8 and 441 terms.
- **(8).** G is invariant under p <-> c and under q <-> w.
- **(9)-(10).** Four charts S, B>, B=, P>= cover every unsafe nonnegative tuple with d = 1. On each chart the cleared G has only positive coefficients.

## B. Dependencies

| dependency | kind | status | match |
|---|---|---|---|
| X ⊂ Sing C | paper | EXACT in moc4lib | as needed |
| sharp SOC (f1, f9; f22) | external doc (F03 family) | CITED | ordered actions to (f2,f15;f26), (f2,f12;f26), (f12,f15;f26) EXACT |
| old-H pair, mapped to (f11, f12; f26) with v(b-d) >= d^2 | external doc | CITED; source statement not pinned down by the document | action and mother triple (f12, f29; f80) EXACT. Numeric test at equality passes |

No paper case, Appendix B or Appendix C result, curve argument or puncture analysis is used.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | generators(F) are 8 rays; rank of P at F is 7; the 7 listed vectors are independent; rows [23,15,17,8,4,16] give I_6 with a zero f26 column | identity | EXACT | replay |
| 2 | f26 + f1 = f0 + f11 + f21 | identity | EXACT | all 24 values |
| 3 | gaps at 20, 11, 3 are q+d, v+d, w+d | identity | EXACT | replay |
| 4 | d = -delta < 0 gives an explicit nonnegative minor decomposition; d = 0 is immediate | identity | EXACT | symbolic identity on all 24 values |
| 5 | (3) from three SOC atoms | citation plus splitting | CITED (SOC); splitting EXACT; actions EXACT | |
| 6 | (4) old-H | citation | CITED | see Deviations. Numeric test: 60 draws at equality give worst -8.9e-17; with d doubled, -7.7e-4 |
| 7 | realizable w satisfies C(w) = 0 and det H(w) = 0 (same polarization argument as core E) | analytic-algebraic | PROOF-READ | |
| 8 | (5)-(6): C and Q; Cst >= d^3 > 0 when d > 0, so Q != 0 | identity | EXACT | |
| 9 | (7) with the 441-term G; v + d > 0, so G = 0 | identity | EXACT | own Bareiss (819 exact divisions), division and square root |
| 10 | (8) symmetries | identity | EXACT | |
| 11 | Coverage. If p, b, c < 1: S. If b > 1: B> (failure of b(p+c) >= 1 and of (4); onto, inverse EXACT). If b = 1: B= (p + c < 1). If b < 1 and max(p,c) >= 1: swap p <-> c to get p >= 1, then P>= (inverse EXACT). q and w are free | case split | PROOF-READ, inverse maps EXACT | the swap preserves G and the set of failed tests |
| 12 | The four cleared tables (852, 3075, 478, 2855 terms; 7260 in total) have only positive coefficients, a strictly positive marker (1, B^6, 2, 2), and denominators with nonnegative coefficients and a positive monomial in the strict parameters | certificate | EXACT | own flint computation. Term counts equal the document's |
| 13 | Conclusion | logic | PROOF-READ | |

## D. Falsification (NUMERIC)

Seed 2026, 8 BFGS restarts per point.

- SOC at equality (40 draws per family): worst 1.7e-16.
- Old-H at equality (60 draws): worst -8.9e-17.
- With d doubled, the same tests detect violations (about -1e-3).
- Realizable points with contact set exactly F (stratum search with margin 0.01, 2 seeds): both have d < 0 (-0.042 and -0.0069), matching the minor-sum branch. Min Phi / max|Pc| = 4.0e-13 and 4.4e-12. G / sum|terms| is about 1e-14 there as well.
- No violation was found.

## E. Verdict

**confirmed-modulo-citations.** All algebraic content is replayed exactly, including the 441-term G and the four tables, and the coverage argument is complete. It rests on two CITED pointwise inequalities: sharp SOC (f1,f9;f22) and old-H. For old-H the document does not name the source statement precisely. Its named group action is exact: the mother triple is paper (f12, f29; f80). The data's "[7, 6]" is a local index list, not paper f-numbering. This is a citation-precision issue, not a gap in the argument: the inequality passes numeric tests at its boundary. A referee should pin down the old-H source statement.
