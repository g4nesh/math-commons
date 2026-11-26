Evidence level: EXACT for every algebraic step, including my own derivation and replay of G and of the whole 149933-term coefficient certificate. PROOF-READ for the logic, with two short arguments I had to supply. CITED for the four capacity atoms and the six zero-face contact cases. NUMERIC for the falsification tests.

# blind_moc4_T48_allocation_II_uniform_20261023 (second-pass audit, label C26-T49-second)

PDF: `proofs/manin_transfer_20261021/blind_moc4_T48_allocation_II_uniform_20261023.pdf`, 2 pages. This is a cell lemma for paper case T_49. It is not a primary case proof.

## Deviations and limits

- I did not open the first auditor's report or anything under `docs/` or `replay/`.
- The supplement's scripts were not run or imported. I loaded two of its JSON files as data (`blind_moc4_T48_retained_atoms_exact_20261023.json` and `blind_moc4_T48_allocation_II_exact_20261023.json`), but only to learn which generator roles the four cited "atoms" use and what form their capacities take. No exact check depends on them. G, D·G and the groups were all re-derived from the paper's definitions (lib `moc4lib`).
- I tested the capacity atoms (two-leaf chains and the old-H inequality) numerically only. They are proved in other documents, which I did not audit.
- I tried to sample spectrally realizable sources with d > 0 by least squares, but they are rare: 3 points from about 590 starts, and none of them inside the cell. So the radial exclusion was not tested numerically on realizable points inside the cell. It rests on the exact certificate, which I replayed.
- The replay script `replay2/..._check.py` holds only the exact checks. It passed 35 of 35 in about 500 s on a loaded machine; its output is saved next to it. The numeric tests in D were run from separate scratch scripts, with the seeds stated in D.
- I did not check whether the cited cases T_1, T_3, T_24, T_39, T_43 and T_44 come before T_49 in the Appendix G order. That is an assembly-level question.

## A. Statements

Contact set (0,3,5,6,9,12), 0-based lexicographic indices. This is exactly the paper's T_49 (lib `CONTACT_SETS[49]`, identity group element). The roles are q,p,j,c,e,d,b,a = f7,f5,f29,f39,f2,f40,f41,f0, using the Appendix D numbering. All eight vanish on T. Only f40 = h_{4213,+} is non-minor. The source is c = q f7 + p f5 + j f29 + c f39 + e f2 + d f40 + b f41 + a f0.

**Theorem 1.** Normalize d = 1 and take a, e, j > 0 and q, p, c, b ≥ 0. Put Q = q + b, u = (1+2/j) / (2(1/a+1/e+1/j)) and v = 1 − u. Assume 0 < u < 1, 2u ≤ a and Q ≤ 2v. Then every x in X with c(x) equal to this source satisfies Phi_{c(x)}(O) ≥ 0 for all O in SO(4). This covers every spectrum (collisions, repeated eigenvalues), the walls a = 2u and Q = 2v, and c = 0. Contacts are non-strict.

Two parts of the theorem's region are handed to other results:
- the faces where one of q, p, j, e, b, a is zero go to other contact cases;
- the case d = 0 is a nonnegative combination of minors.

Scope remark (auditor). N_{T_49} has 10 generators: f0, f1, f2, f5, f7, f9, f29, f39, f40, f41. The document uses only 8. The exact relations are f1 = q + c + a − d and f9 = q + e + b − d, in role vectors. From these I proved:
- On N_T ∩ {d > 0}, all 8 role coordinates are nonnegative. So the document's "nonnegative supported source" with d = 1 is all of that part of the cone.
- When the d coordinate is ≤ 0, the source is a nonnegative combination of minor generators. Then Phi ≥ 0 holds trivially.

The proof is in C, step S1.

## B. Dependencies

| ref | kind | used for | status |
|---|---|---|---|
| chain (a,e): Phi(z f5 + x f0 + y f2 + t f40) ≥ 0 if z ≥ t²/x + t²/y | common lemma (bad70 chain orbit) | pointwise allocation | CITED; numeric pass |
| chain (q,j): the same with leaves f7, f29 | common lemma | allocation | CITED; numeric pass |
| chain (j,b): the same with leaves f29, f41 | common lemma | allocation | CITED; numeric pass |
| old-H: Phi(a f29 + b f5 + d f40) ≥ 0 if a(b − d) ≥ d² | common lemma | allocation | CITED; numeric pass |
| X ⊂ {grad C = 0} | core (paper App. A) | realizable ⇒ G = 0 | EXACT in lib (only this direction of HMSV is needed) |
| zero face q=0: contacts T ∪ {16} ⊇ g·T_39 ("T38") | paper case | boundary | CITED; containment EXACT |
| p=0: T ∪ {11} ⊇ g·T_1 ("T0") | paper case | boundary | CITED; containment EXACT |
| j=0: T ∪ {2} ⊇ g·T_3 ("T3") | paper case | boundary | CITED; containment EXACT |
| e=0: T ∪ {13} ⊇ g·T_24 ("T23") | paper case | boundary | CITED; containment EXACT |
| b=0: T ∪ {7} ⊇ g·T_43 ("T42") | paper case | boundary | CITED; containment EXACT |
| a=0: T ∪ {23} ⊇ g·T_44 ("T43") | paper case | boundary | CITED; containment EXACT |

The supplement labels are historical labels. In `cases/index.tsv` they map to paper cases as T0 = T_1, T3 = T_3, T23 = T_24, T38 = T_39, T42 = T_43 and T43 = T_44. Each zero-face contact set is the set of permutations where the other seven role generators all vanish. Each contains a group image of the listed paper case (lib `find_g`). Each zero-face set also contains T_49 itself.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| S1 | the 8 roles are generators of N_{T_49}; N_T ∩ {d>0} is the nonneg orthant in role coordinates; d ≤ 0 is minor-only | identity + argument | EXACT / PROOF-READ (argument supplied by auditor) | rank 8, f1 and f9 relations. Proof: starting from a representation by all 10 generators, move weight y from f1 or f9 onto (q,c,a) or (q,e,b) and lower d by y. |
| S2 | the four capacity atoms | citation | CITED | numeric (D) |
| S3 | the j resources sum to j, the non-minor currents sum to 1, and the atom requirements sum to p* | identity | EXACT | sympy |
| S4 | u(1/a+1/e) − v/j = 1/2 and p* = 1 − u/2 + v/j − Q/4 | identity | EXACT | sympy |
| S5 | p* − (1/2 + v/j) = (2v − Q)/4 ≥ 0 | identity | EXACT | sympy |
| S6 | p ≥ p* ⇒ Phi ≥ 0 (sum of atoms, plus (p − p*) f5 and c f39, which are minors) | inequality | PROOF-READ, modulo S2 | zero currents are dropped correctly |
| S7 | realizable ⇒ C(y) = 0 and det H(y) = 0, where y = Re x = c + rκ and r = Re F_{π*} | analytic | PROOF-READ (argument supplied by auditor; the document only says "already verified") | see note R below |
| S8 | C(c + rκ) = C(c) + r·Q0, with Q0 = (c+1)(1+e+j) − p(e+j) | identity | EXACT | flint |
| S9 | Q0 = 0 is impossible for realizable sources: (e+j)C − ejQ0 has 42 positive terms, including e and j | identity + sign | EXACT | flint |
| S10 | det H(c − (C/Q0)κ) = −4G²/Q0², with G primitive and 1772 terms | identity | EXACT | 14x14 Bareiss determinant over Z[q,p,j,c,e,b,a,r]; the determinant is cubic in r; the result is −4 times an exact square |
| S11 | the chart is onto the cell ∩ {p < p*, Q < 2v}, the e formula is forced, and all denominators are > 0 | identity | EXACT | inverse maps recomputed |
| S12 | D·G is a polynomial with 149933 terms | identity | EXACT | substitution and exact division |
| S13 | −DG has negative coefficients only in 18 C1-groups, each equal to (C1−1)²R(C1) with R ≥ 0; the term 432B⁶S² is present; hence G < 0 | sign certificate | EXACT | flint (54 negative coefficients in all) |
| S14 | walls: A = 0 and c = 0 are covered by the same signs; Q = 2v is the B⁶ slice (degree in B is 6, the 18 groups lie in that slice, the marker survives); p = p* is S6 | wall | EXACT / PROOF-READ (I supplied the continuity in W_B for B > 0) | |
| S15 | zero faces go to the cases T_39, T_1, T_3, T_24, T_43, T_44 | citation | CITED (containment EXACT) | find_g |
| S16 | d = 0 is minor; homogeneity removes the normalization | argument | PROOF-READ | f40 is the only non-minor generator |

**Note R (S7).** Write the cubic C as T(u,u,u), with T symmetric and trilinear. If x = y + iw is in X, then grad C(x) = 0. Taking real and imaginary parts gives T(y,w,·) = 0, that is H(y)w = 0, and T(y,y,·) = T(w,w,·). Pair the second equation with y: 3C(y) = T(y,y,y) = T(w,w,y) = T(y,w,·)·w = 0.

- If w ≠ 0, then H(y)w = 0, so det H(y) = 0.
- If w = 0, then H(y)y = 6T(y,y,·) = 0 with y ≠ 0, so again det H(y) = 0.

Since y = c + rκ and C(c + rκ) = C + rQ0 (S8), we get r = −C/Q0 when Q0 ≠ 0 (S9), and then G = 0 (S10). The argument uses only the closed condition grad C = 0, so it also covers collision limits and degenerate spectra. Numeric sanity check: on 200 random points of X, max |C(Re x)|/|Re x|³ = 4e-17 and min singular-value ratio of H(Re x) = 7e-17.

## D. Falsification (NUMERIC; lib `min_phi_over_SO4`, multistart BFGS)

- The atoms at equality capacity, 25 random leaf values each (numpy seed 2026), normalized worst min Phi:
  - chain_ae: 1e-17
  - chain_qj: 9e-18
  - chain_jb: −9e-17
  - old-H: 5e-18
- Power check: the optimizer does find negative values when the centre coefficient p is too small. f40 alone gives −0.0674. chain_ae at p = 0.2 × capacity gives −2.5e-4. old-H at p = 0.8 × capacity gives −5.3e-5. The chains show no negative value already at 0.5 × capacity, so their stated capacities are not sharp. That is harmless.
- The full allocation at p = p*, on 30 random cell points (seed 7): worst min Phi = −1.3e-15.
- Realizable sources with d = 1, found by least squares over spectra (seeds 11 and 23; 3 points):
  - |G|/Σ|terms| ≤ 3e-17, consistent with S7–S10;
  - min Phi ≥ 2.7e-14;
  - none of the 3 points lies in the cell.
- The library sampler found 396 points on the T_49 stratum. All had d ≤ 0, where the source is minor-only (S1).

I found no violation.

## E. Verdict: confirmed-modulo-citations

Citations: the four capacity atoms (two-leaf chain and old-H) and the paper cases T_1, T_3, T_24, T_39, T_43 and T_44 for the zero faces.

Presentation-level gaps, each repaired above:
1. p.2, "provides the primitive 1772-term polynomial G … All spectrally realizable sources satisfy G = 0". The document does not say why realizability forces the radial equation C(c + rκ) = 0 or det H = 0; it says only that the input was "already verified". The claim is true and the derivation is short (Note R).
2. p.1, the roles. The document drops generators f1 and f9 of N_{T_49} without comment. This is harmless by S1.
3. p.2, "The face c = 0 is not covered and was stated throughout." The wording is garbled. c = 0 is in fact covered by the coefficient certificate, since c ≥ 0 is allowed there.
