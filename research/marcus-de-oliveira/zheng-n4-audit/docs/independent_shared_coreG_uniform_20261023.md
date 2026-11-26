# Audit: independent_shared_coreG_uniform_20261023 (seven-contact set G)

Evidence level: every algebraic identity and every positivity certificate EXACT (own code, replayed from the definitions); the case logic PROOF-READ; the three kinds of "atoms" (sharp SOC, old-H) CITED and only NUMERIC-tested. Verdict: **confirmed-modulo-citations**.

PDF: `proofs/manin_transfer_20261021/independent_shared_coreG_uniform_20261023.pdf` (5 pages). Replay: `replay/independent_shared_coreG_uniform_20261023__check.py` (output in `__check_output.txt`; about 4 minutes; exits 0).

## Deviations and limits

- The supplement's companion scripts and JSON (`independent_shared_coreG_radial_20261023.py/json`) were neither run nor read. R9 was re-derived from scratch: the lift w(z) = sum of coefficients times f_j + z*kappa is built from the shared library's f_j and kappa. The 14x14 Hessian of the cubic (A.2) is formed and its determinant computed by my own fraction-free Bareiss elimination (819 exact divisions). Then z = -C/Q is substituted and the result factored.
- The bounded-sector triangulation is my own exact pulling triangulation of the polytope (10). The document describes its triangulation only in words. Mine also has 13 simplices, and its 13 certificate sizes coincide with the document's list as a multiset (total 352583).
- I did not audit the cited atoms (sharp SOC from the reference family (f1, f9; f22), old-H from (f3, f15; f26)). I checked exactly that a group element maps the reference profiles onto each target, and I tested every target atom numerically at its boundary.
- The numeric falsification of the theorem is weak. Every floating point I found on the realizable stratum with contacts containing G has d <= 0 (up to 1e-9), so it falls in the trivial branch.

## A. Statements

Notation: g = P c(x) is the gap vector, with permutation indices 0..23 as in the paper. The document's f_i is the paper's f_i (Appendix D). This was checked: generators(G) from the library is [2,3,4,11,12,15,21,26], exactly the document's "active source list".

**Theorem 1 (p. 1).** Let G = {0,1,2,7,9,10,23}. This is literally the paper's Table 5 seven-set no. 10 (EXACT). Take x in X with P c(x) >= 0 and (P c(x))_pi = 0 for every pi in G (non-strict contacts). Then Phi_{c(x)}(O) >= 0 for every O in SO(4), hence for every O in O(4).

- Spectra: the proof uses only that Re x is a sum of two points of Sing C. So it covers every x in X, including repeated eigenvalues, a_i + b_j = 0 and limit points.
- The unitary extension via the "two-real-response theorem" is not needed for the paper's (*) and was not audited.
- Coordinates: g = a f2 + t f3 + u f11 + b f12 + c f15 + v f21 + d f26, with (a,t,u,b,c,v) = (g15,g14,g17,g8,g4,g16) >= 0 and d = g3 - g16 of either sign (EXACT).

## B. Dependencies

| dependency | kind | status | use |
|---|---|---|---|
| Sharp SOC atom: x f_i + y f_j + delta f26 has Phi >= 0 for all O when xy >= delta^2, for {i,j} among {2,12,15}; reference family (f1,f9;f22), source document not named precisely (several core documents state it) | other-doc | CITED; group maps (f1,f9;f22) -> (f_i,f_j;f26) EXACT for all three pairs; NUMERIC at xy = delta^2: 3 x 20 samples, worst min Phi 2.6e-14 | (4) |
| Old-H atom: X f_aux + m f_main + delta f26 has Phi >= 0 when X(m - delta) >= delta^2; reference (f3,f15;f26) | other-doc | CITED; ordered group maps to (f3,f15;f26) and (f11,f12;f26) EXACT; NUMERIC 2 x 20 samples, worst 1.9e-16 | (5),(6) |
| Phi of a minor vector is a squared minor, so it is >= 0 (paper (4.2)) | other-doc | EXACT in shared lib, with a PROOF-READ general proof | d <= 0 branch, leftover minors |
| X is contained in Sing C (paper A.1), and conj(X) = X | other-doc | EXACT in shared lib; conjugation closure PROOF-READ (B_j have integer coefficients) | F(w) = 0 and det H(w) = 0 |

Both atoms are used as unconditional inequalities on V_R: the split pieces are not themselves realizable. That is how they are stated in the document. The numeric tests are consistent with that reading. No other contact case is used, and the document says so.

## C. Steps

| # | claim (location) | kind | status | how checked |
|---|---|---|---|---|
| 1 | contact rows have rank 7; [f2,f3,f11,f12,f15,f21,f26] is a basis of the contact kernel; coefficient formulas (p. 1) | identity | EXACT | own linear algebra on library P and f_j |
| 2 | (3): g21 = a+d, g20 = t+d, g11 = u+d, g3 = v+d; f4 + f26 = f2 + f3 + f11 + f21 | identity | EXACT | all 24 coordinates |
| 3 | d <= 0: g is a nonnegative combination of the minor profiles f2,f3,f4,f11,f12,f15,f21, so Phi >= 0 | identity + (4.2) | EXACT | rewrite checked in all 24 coordinates; minor flags |
| 4 | (4): SOC atoms sum to give sufficiency for d <= H3(a,b,c), via the balanced and dominant allocations | inequality | PROOF-READ (atoms CITED) | balanced capacities sum to s/2; dominant capacities sum to sqrt(M(s-M)); scaling down is allowed because xy >= delta^2 is monotone in delta |
| 5 | (6): l + n >= d is sufficient, with l = cap(u,b) and n = cap(t,c) | inequality | PROOF-READ (atoms CITED) | cap solves X = lambda^2/(m - lambda) (EXACT, sympy); the two resource pairs are disjoint |
| 6 | F(w(z)) = C + Qz with C, Q as printed; C > 0, so Q != 0 | identity | EXACT | own polynomial expansion of (A.2) |
| 7 | a realizable w = Re x satisfies F(w) = 0 and det H(w) = 0 | analytic | PROOF-READ | w = (x + conj x)/2 with both in Sing C. Cubic polarization gives C(x+y) = 0 and H(x+y)(x-y) = 0. If x = conj x, then H(w)w = 2 grad C(w) = 0 with w != 0 (because d > 0). |
| 8 | (7): Q^4 det H(w(-C/Q)) = -4 Q^2 R9^2; R9 has 536 terms, degree 9 | identity | EXACT | own Bareiss, 819 exact divisions; deg_z det H = 3; exact sqrt |
| 9 | (8): R9 is invariant under (t,u,b,c) -> (u,t,c,b), and so are the failure conditions | identity | EXACT / PROOF-READ | polynomial equality; H3 is symmetric, and {l,n} are swapped |
| 10 | c = 0 face: R9 = (u+d)(a+d)(v+d)^2(t+d)^2 K with K > 0 when ab < d^2; b = 0 by symmetry | identity + inequality | EXACT / PROOF-READ | exact factorization; (u+b+d)(d^2 - ab) > 0 and the other terms are >= 0 |
| 11 | bounded sector: a,b,c <= 1 and H3 < 1 imply a + b + c < 2; the point lies in polytope (10) with the four forms strictly positive | inequality | PROOF-READ | if s >= 2 the balanced case applies and H3 = s/2 >= 1 |
| 12 | polytope (10) has 13 integer vertices from 252 hyperplane 5-sets | finite | EXACT | own Fraction enumeration |
| 13 | 13-simplex closed coverage | finite | EXACT | own pulling triangulation; volume 17/120 equals that of a second pulling triangulation (17 simplices, other apex) |
| 14 | (11): each (B-L)^2(C-N)^2 R9(list) has positive coefficients, and every barycentric support meeting the four strict forms contains a positive V^0 monomial | certificate | EXACT | sizes match (multiset; total 352583); support test over all 63 subsets per simplex |
| 15 | a > 1 implies a(b+c) < 1; (12) covers the hard domain; X^2Y^2 R9(list12) has 56416 positive coefficients with marker 2X^2Y^2 | certificate + inverse | EXACT / PROOF-READ | inverse map verified with sympy; X, Y > 0 because b > l and c > n |
| 16 | b > 1 (and c > 1 by symmetry); (13); Y^2 V_D^2 R9(list13) has 141628 positive coefficients with marker 2Y^2; V_n and V_D have positive coefficients | certificate + inverse | EXACT / PROOF-READ | inverse map verified with sympy; K >= 0 because l + n < 1 |
| 17 | the regions are exhaustive after d = 1 | case split | PROOF-READ | at most one main exceeds 1 (pairwise product < 1); c > 1 maps to b > 1 by symmetry (8) |

## D. Falsification

- Atoms, at the boundary of each atom's own condition: 5 atoms x 20 random coefficient samples, multistart BFGS over SO(4) with 8 restarts, seed 20261125. The worst minimum is 2.6e-14. All minima are about 0, so the atoms are attained but not violated.
- Theorem: 12 floating points with contacts containing G (library stratum search, seeds 0..11), 10 restarts each. The worst min Phi / max gap is -1.0e-10, which is roundoff at the stratum tolerance.
- Sanity check of the necessary condition at these points: C(Re x) is about 0 and H(Re x) is singular to 1e-17 relative. Where R9 != 0, Q is also about 0, as expected (scratch run, not in the replay).

## E. Verdict

**confirmed-modulo-citations.** Every computation the proof relies on was reproduced exactly from the definitions: the radial polynomial R9, its symmetry, the zero-face factor, the 13 vertices, the triangulation, the 15 certificates (with matching term counts) and the inverse maps. The case logic is complete. The citations are the sharp-SOC and old-H atoms.

Presentation notes, none affecting correctness:
1. Two displays are both labelled "(1)" on p. 1.
2. The triangulation is described only in words ("standard recursive convex-pyramid pulling").
3. The statement refers to "the stated support inequality" and "the stated simultaneous two-real-response reduction" without naming the source documents.
