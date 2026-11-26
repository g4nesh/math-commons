# Finite classification (Section 4): independent audit

Evidence level: EXACT for every finite assertion listed below (own integer/rational code,
`finite_classification.py`, 51 claim lines, all PASS); PROOF-READ for the completeness
argument of Proposition 4.1 and for the proof of Lemma 4.2; NUMERIC only for the
comparison with the float data `minors/moc/facets4.npy` (task 7).

Paper: D. Zheng, "The Marcus-de Oliveira conjecture for n = 4" (24 Oct 2026), Section 4
(pp. 21-25), Appendix A (p. 35), Appendix D (p. 59), Appendix E (p. 60), Table 3 (p. 31),
Table 5 (pp. 62-63).

## Deviations and limits

- Coordinates. The work is done in the indicator embedding e_sigma in {0,1}^52
  (coordinates chi_{I,J}(sigma) = [sigma(I) = J], |I| = |J| in {1,2}). To make the
  identification with the paper's Q exact rather than asserted, the script also re-derives
  the matrix P from the matching basis (A.1) by expanding prod_i (z_{4+sigma(i)} - z_i) in
  B_0..B_13, and checks that P and [1 | E] have the same column space. So the polytope,
  facets and orbits below are the paper's, not an analogue.
- Transcriptions. Appendix D, Appendix E, Table 3 and Table 5 were transcribed from the
  text extract and compared with the page images (pp. 59, 60, 62, 63). One entry of
  Appendix D (f_79 = h_{1243,-}) was dropped in my first transcription; the script caught
  this, because it requires the 82 entries to be a bijection onto the 82 facets.
- Supplement. Only two JSON files were read, with `json.load`
  (`data/root_permutation_polytope_complete_20261022.json`, `data/facet-definitions.json`).
  No supplement code was executed or imported, and `moc4lib.py` was not used.
- Scope. Out of scope here: the citation [19] (Kovacec-Ribeiro), the G-equivariance of X
  and of the Leibniz sum (eq. (4.7)), and the ninety case proofs. Section 5's claim that
  the gaps Gamma_j, j in {2,3,8,10,13,16,19,20,21,22,23}, are strictly positive on X85+
  was not audited. Only its bookkeeping (Table 3) was checked.
- Runtime. About 90 s through RUNJOB. Most of it (about 55 s) is the exhaustive search over
  13-subsets.

Reproduce (from this directory):
`MOC4_SUPP=<supplement root> python -I -B finite_classification.py --json finite_classification.json --lane-npy ../../minors/moc/facets4.npy`

## 1. Facets of Q (Proposition 4.1): EXACT

- The 52 functions chi_{I,J} on S_4 span a space W of dimension 14, and W contains the
  constants. The affine hull of the 24 points e_sigma has dimension 13.
- P re-derived from (A.1) is integral. It has rank 14 and last column 1 (P kappa = 1). The
  14x14 minor on rows {0,1,2,3,4,6,7,8,9,10,12,13,16,18} is 1, as stated in Appendix A.
  P equals the supplement's `vertices` array entry by entry.
- The candidate list has 16 n_{I,J} (|I| = 1), 18 n_{I,J} (|I| = 2, up to complement) and
  48 h_{pi,eps}. Each lies in W and is >= 0 at all vertices.
  - Value patterns: (18 zeros, 6 ones), (20, 4) and (13 zeros, 10 ones, 1 two).
  - Each zero set has rank 13, and the 82 zero sets are distinct.
  - G permutes these 82 vectors in three orbits (sizes 16, 18, 48) and preserves the
    minor/non-minor type.
- My own algorithm finds them. Exact double description on the cone {c : v_sigma(c) >= 0}
  in R^14 gives exactly 82 extreme rays, and they are exactly the candidates: 16 + 18 + 48,
  with no difference from the paper.

The completeness of my list does not depend on the double-description code. Two separate
exact certificates support it:

1. Exhaustive search. Every facet of the 13-polytope Q contains 13 affinely independent
   vertices, and these determine the facet's hyperplane. The script runs over all
   C(24,13) = 2,496,144 vertex 13-subsets. 1,632,360 of them are affinely independent. For
   each of those it computes the unique hyperplane through the subset and tests whether all
   24 vertices lie on one side. The supporting hyperplanes found are exactly the 82
   candidates.
2. Ridge closure. The ridges of each listed facet are computed exactly.
   - For the 48 simplicial facets, the ridges are the thirteen 12-subsets.
   - For one representative of each of the two minor-facet orbits, the ridges come from an
     exhaustive search over the 12-subsets of its vertex set. There are 39 ridges for the
     n_{i,j} type and 57 for the n_{ij,kl} type. They are transported to the rest of each
     orbit by G.
   - Every ridge lies in exactly two listed facets, and the facet-ridge graph is connected.
   - Since the dual graph of a polytope is connected and a ridge lies in exactly two
     facets, a set of genuine facets that is closed under crossing ridges is the set of
     all facets.

## 2. Replay of the paper's certificate (pp. 23-24): EXACT

All the following hold:
- There are 6846 distinct, sorted 13-element vertex sets.
- Each lies in exactly one of the 82 hyperplanes.
- Every listed facet contains at least one of the simplices (between 1 and 306 per facet).
- det[24(v_i - b)] (13x13, vertices in sorted order) is nonzero for all of them (3423
  positive, 3423 negative).
- The boundary of z = sum_k sgn(det_k) [sorted simplex_k] is zero. It involves 44,499
  distinct 12-faces, and each face is shared by exactly two simplices with opposite induced
  signs.

Two consistency checks also pass. The JSON facet list is in Appendix D order, and
`facet-definitions.json` agrees with my transcription of Appendix D. As a control, the
unsigned chain (all coefficients +1 in sorted order) is not a cycle: 24,839 faces fail to
cancel. So the orientation is essential, and the supplement's data supports the oriented
claim.

Orientation and the ordered-versus-oriented question. The data stores unsigned sorted
13-tuples. The orientation is not data; it is the sign of det[v_{i_0}-b, ..., v_{i_12}-b]
in that order. That sign is the orientation of the 13-simplex [b, v_{i_0}, ..., v_{i_12}].
The boundary formula gives that 13-simplex the face [v_{i_0}, ..., v_{i_12}] with
coefficient +1. So sgn(det) times the sorted simplex is the 12-simplex oriented as the
boundary of a positively oriented 13-simplex with apex b. Since b is interior (b is the
barycentre of all vertices) and the face lies in a supporting hyperplane strictly
separating it from b (det != 0), this is the outward boundary orientation of Q. Every
term of z therefore has local multiplicity +1 with respect to one fixed orientation of the
sphere bd Q. "Consistently oriented outward" in the paper means exactly this, and it is
what the replay checks: the 12-faces are compared as sorted tuples with coefficient
sgn(det)(-1)^m.

Validity of the topological argument (PROOF-READ): valid. The steps are:
- (a) Each simplex lies in bd Q, because its vertices are zeros of an f_j that is >= 0 on
  Q (checked).
- (b) z is a singular 12-cycle in bd Q, which is homeomorphic to S^12, so
  [z] = d [bd Q] in H_12(bd Q) = Z.
- (c) Take a point p outside the image of all 11-faces. The image of [z] in
  H_12(bd Q, bd Q - p) = Z is the sum of the local degrees of the affine simplices
  containing p. Each is +1 by the orientation argument above. Since bd Q is a connected
  oriented manifold, the map H_12(bd Q) -> H_12(bd Q, bd Q - p) is an isomorphism for
  every p. Hence the number of simplices containing p equals d for every such p.
- (d) Points that lie in exactly one simplex exist, because a finite union of 11-simplices
  has empty interior. Hence d >= 1.
- (e) Every point outside the 11-skeleton is covered. That set is dense and the support of
  z is closed, so the support of z is all of bd Q.
- (f) The relative interior of a facet F meets no other facet, and it is a nonempty open
  subset of bd Q. If F were unlisted, its relative interior would avoid the support of z,
  which lies in the union of the listed facets. Contradiction.

The argument does not need the simplices to have disjoint interiors, as the paper says.
Here each 12-face happens to be shared by exactly two simplices, so the data is in fact a
pseudomanifold.

Minor presentational points, not gaps. The paper does not say that the orientation sign
is taken with respect to the sorted vertex order. It also does not say why the relative
interior of an unlisted facet avoids the listed ones (two distinct facets meet only in
their relative boundaries).

Verdict on Proposition 4.1: correct. The facet list is complete, by three independent
exact routes and by an exact replay of the paper's own certificate.

## 3. Numbering (Appendix D): EXACT

f_0..f_81 are stored in `finite_classification.json` (`facets`) with spec, type, the 24
vertex values and the zero set. The Appendix D table is a bijection onto the 82 facets and
matches the supplement's ordering.

## 4. Orbits (Section 4.2, Appendix E): EXACT

- G was generated as a permutation group on the 24 vertices. It has order 1152 and equals
  {s -> r s r'} together with {s -> r s^-1 r'}.
- The 134,596 six-subsets fall into **214 orbits**. Orbit sizes and counts:
  16:1, 24:2, 36:1, 72:6, 96:7, 144:9, 192:1, 288:40, 576:85, 1152:62.
- **124 orbits are minor-only and 90 are not.** Among the 90 non-minor orbits the sizes
  are 24:1, 72:3, 96:2, 144:2, 192:1, 288:16, 576:31, 1152:34.
- The paper's lexicographically smallest representatives are reproduced (`orbits` in the
  JSON).
- **Appendix E matches.** T_1..T_90 lie in 90 distinct orbits, and these are exactly the 90
  non-minor orbits. Only T_1 and T_2 are the lexicographic minima of their orbits; the paper says the listed sets are "equivalent representatives", which is
  consistent.
- Non-minor generators per case: one in 51 cases, two in 35, three in 3 (T_13, T_33, T_54)
  and four in 1 (T_22). Per-case generator lists are in `cases[*].generators`.
- For every non-minor orbit, dim N_T = 8.

## 5. Larger contact sets and Table 5: EXACT, with logic PROOF-READ

The logical point holds. Let S be a contact set with |S| >= 6 and T any 6-subset. Then
N_S is contained in N_T (more equalities), and T lies in one of the 214 orbits, because the
orbits partition all 6-subsets (checked: 134,596 sets assigned). So a proof for every
N_{T_nu} in full covers every S. The contact conditions are non-strict, as the paper says
on p. 21. This also needs the G-transfer (4.7) and the invariance of X, which were not
audited here.

Of the 90 x 18 = 1620 one-element extensions T_nu + {j}:
- 797 have a minor-only cone N_{T+j}, so they are trivially fine.
- 823 still have a non-minor generator.
- Per case, the count still needing an argument is 7 for 51 cases, 11 for 16, 12 for 15,
  13 for 4, 14 for 2 and 15 for 2.
- The full lists are in `cases[*].sup7_minor_only` and `cases[*].sup7_needs_argument`. The
  second gives the surviving non-minor f_j and the T_nu of each six-subset T - x + j.

Example, T_85: j in {2,10,13,16,21,23} is minor-only, and j in
{1,3,4,6,8,11,12,15,17,19,20,22} is not. Section 5 handles {1,4,6,11,12,15,17} through
Table 3, and every row of Table 3 checks exactly: A_{rho,rho',eps}(T_nu) equals the stated
image, which is a six-subset of T_85 + {j} containing j. It excludes {3,8,19,20,22}
together with the six minor-only ones by the claim that those gaps are strictly positive
on X85+. That claim is a Section 5 argument and is outside this lane.

Table 5. All 27 auxiliary seven-sets contain at least one six-subset from a non-minor
orbit, and each has a non-minor generator itself. Contained T_nu orbits:

| j | set | T_nu orbits among its 6-subsets | non-minor f |
|---|---|---|---|
| 1 | {0,1,2,3,6,15,17} | 8,9,15,36,45,61,67 | 19 |
| 2 | {0,1,2,3,6,15,22} | 8,14,38,61,68 | 19 |
| 3 | {0,1,2,3,6,17,19} | 9,11,17,38,51,65,76 | 19 |
| 4 | {0,1,2,3,6,17,22} | 9,14,40,64,71 | 19 |
| 5 | {0,1,2,6,9,11,20} | 24,26,37,44,65,67,73 | 18 |
| 6 | {0,1,2,6,9,11,23} | 24,28,39,42,47,49,66 | 18 |
| 7 | {0,1,2,7,9,10,14} | 50,72,80 | 26 |
| 8 | {0,1,2,7,9,10,16} | 50,51,55,61,62,63,73 | 26 |
| 9 | {0,1,2,7,9,10,17} | 50,52,56,62,72,73,74 | 26 |
| 10 | {0,1,2,7,9,10,23} | 50,52,54,55,61 | 26 |
| 11 | {0,1,2,7,9,16,17} | 51,52,63,73,87 | 26 |
| 12 | {0,1,2,7,9,16,23} | 51,54,64 | 26 |
| 13 | {0,1,2,7,10,16,17} | 55,56,63,80 | 26 |
| 14 | {0,1,2,7,10,16,23} | 55,63,64,85 | 26 |
| 15 | {0,1,2,7,14,16,23} | 61,64,87 | 20,26 |
| 16 | {0,1,2,9,10,14,17} | 72,74,88 | 26 |
| 17 | {0,1,2,9,10,16,17} | 73,74,80 | 26 |
| 18 | {0,1,3,10,11,12,16} | 76,77,78,79 | 30 |
| 19 | {0,1,4,6,9,12,17} | 53,54,60 | 35 |
| 20 | {0,1,5,7,9,12,17} | 57,58,59,78,86 | 35 |
| 21 | {0,1,6,9,16,17,22} | 81,83,85 | 31 |
| 22 | {0,2,3,7,11,12,23} | 32,33,89 | 18,40 |
| 23 | {0,3,5,6,9,12,17} | 48,49,59,84 | 40 |
| 24 | {0,3,7,11,12,16,20} | 89,90 | 30 |
| 25 | {1,3,4,8,11,12,15} | 30,32,79,82,89 | 56 |
| 26 | {1,3,5,6,7,8,14} | 47,48,67,77 | 28 |
| 27 | {1,4,5,7,9,12,14} | 67,68,71,76,87 | 35 |

## 6. Lemma 4.2 (p. 24): PROOF-READ correct; EXACT for all 214 representatives

The proof is correct given Proposition 4.1:
- Elements c of V_R correspond bijectively to affine functions on aff Q (P is injective of
  rank 14, and aff Q has dimension 13). Under this correspondence kappa is the constant 1.
- Since the 82 facets give an irredundant H-description of Q, the affine Farkas lemma
  writes any affine g >= 0 on Q as g = lambda_0 + sum_j lambda_j f_j with all lambda >= 0.
- Evaluate at v^T (the average of the vertices in T). Every term is >= 0 and the total is
  0, so lambda_0 = 0, and lambda_j = 0 unless f_j(v^T) = 0. Since f_j >= 0 on T,
  f_j(v^T) = 0 is the same as j in F_T.
- The paper's kappa remark is exactly the step lambda_0 = 0 (v^T(kappa) = 1). It is right,
  although "lifting back to V_R may leave a multiple of kappa" would be clearer stated as
  the constant term of affine Farkas.
- The alternative argument is also valid: the cone {c : Pc >= 0} is pointed with extreme
  rays f_0..f_81 (Prop. 4.1), and complementary slackness at T follows.
- Nothing uses |T| = 6.

Exact check: for each of the 214 lexicographic representatives, the script runs a fresh
double description of N_T. It restricts to the subspace {c : v_pi(c) = 0, pi in T} and
uses the 24 inequalities, with no use of the facet list. The extreme rays equal
{f_j : j in F_T} in every case. The dimensions of N_T over the 214 orbits are
1:1, 2:2, 3:5, 4:13, 5:21, 6:33, 7:22, 8:116, 9:1. The facet lattice also makes this
expected: N_T is a face of the cone, so its extreme rays are those f_j that vanish on T.

## 7. Cross-check with `minors/moc/facets4.npy`: NUMERIC (data are floats)

I recomputed the lane's 13 reduced coordinates exactly (M1[0:3,0:3] and the pairing matrix
N[0:2,0:2]). Exactly, [VERTS | 1] has rank 14 and the same column space as [1 | E], so the
lane's Q is the same polytope. The 82 rows of `facets4.npy` evaluated at the 24 vertices
have zero sets, at tolerance 1e-9, exactly equal to my 82 facets. Their value vectors are
proportional to mine, with a maximum relative deviation of 4.4e-16. **They are the same 82
hyperplanes.** That lane's facets came from scipy ConvexHull (floating point), so they are
numeric evidence only.

## Discrepancies

None in the finite claims. All counts agree with the paper:
- 82 = 16 + 18 + 48 facets
- 6846 simplices
- 44,499 faces
- 214 orbits = 124 minor-only + 90 non-minor
- Appendix E consistent with the 90 non-minor orbits
- Table 5 consistent

Points for the rest of the audit:
- (i) The orientation of the certificate is implicit, not stored. It is the sign of the
  determinant in sorted order relative to b.
- (ii) For T_85, five of the twelve non-minor one-element extensions ({3,8,19,20,22}) are
  handled not by smaller cases but by a strict-positivity claim in Section 5.2. That needs
  checking by the Section 5 lane.
- (iii) Proposition 4.4's "any six zero coordinates" is logically sound as far as the
  finite part goes. It depends on (4.7) and the G-invariance of X for the transfer from
  T_nu to its orbit.

Files: `finite_classification.py`, `finite_classification.json`,
`finite_classification_output.txt`, this report.
