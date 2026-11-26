# Second-pass audit: independent_moc4_six_contact_reduction_20261022

Evidence level: confirmed-modulo-citations. Algebra and finite combinatorics are EXACT (own code). The curve, subharmonic and induction arguments are PROOF-READ. The matching-cone realization, the one-circle theorem (paper Thm C.1), the two-orthogonal theorem (paper Thm 1.2), facet completeness and the minor-square identity (4.2) are CITED. Falsification tests are NUMERIC.

Document: `proofs/manin_transfer_20261021/independent_moc4_six_contact_reduction_20261022.pdf`, 8 pages, dated 22 Oct 2026. It is the primary proof of paper case T_1 = {0,1,2,3,6,7} (its Prop. 10) and of paper case T_2 = {0,1,2,3,6,8} (its Prop. 13, which the document labels "T1").

Replay: `replay2/independent_moc4_six_contact_reduction_20261022__check.py`, with output in `..._check_output.txt`. All 57 EXACT checks pass and the script exits 0.

## Deviations and limits

- This is an independent second pass. I did not open the first auditor's report.
- The only supplement file I used is the facet JSON `root_permutation_polytope_complete_20261022.json`, loaded as data. I used it only to translate the document's facet numbering f_j into 24-entry profiles. Those profiles were then matched exactly to the shared library's own facet list, which was derived from the paper. None of the supplement's scripts named in Remarks 1, 9, 11 and 14 were run or relied on.
- The library's exact facet-completeness certificate (Prop. 4.1) and its proof of (4.2) were not redone here. They are CITED from `lib/README.md`.
- The group equivariance (4.7) is also CITED from the library tests.
- Standard facts are recorded as standard rather than re-proved: Birkhoff's theorem, removability of bounded subharmonic singularities, the maximum principle on compact Riemann surfaces, the affine dimension theorem, and the existence of an irreducible curve through two points of an irreducible variety.
- The numeric samplers mostly land on degenerate configurations (see D). Direct sampling of the nondegenerate T_1 stratum was not achieved. This is consistent with the document, which proves that this stratum forces six more contacts.

## A. Statements (paper notation; F_sigma = lambda prod_i (a_i + b_sigma(i)); gaps are real parts; permutation indices are 0-based lexicographic)

- **(M1)-(M4) (p. 1-2)** are the inputs. They say:
  - P is 24x14 of rank 14, and X is the 6-dimensional irreducible matching cone.
  - Every nonzero point of X has finite spectra with lambda != 0.
  - X lies in Sing(C), grad C(kappa) = 0, and (2) holds: grad C(x + t kappa) = grad C(x) + t H_kappa x.
  - ker H_kappa = span{e_ij} and rank H_kappa = 4.
  - For each O in O(4) there is a real linear functional Phi_O with Phi_O(x) = lambda det(D_a + O D_b O^t) on matching points, Phi_O(kappa) = 1 and Phi_O(e_ij) = O_ij^2.
- **Lemma 2 (lem:curve, p. 2)** concerns an irreducible affine curve C carrying regular functions f_1..f_m and h. Assume:
  - not all f_j are constant, and sum_j Re f_j = 1 on C;
  - L = Re h is bounded on {all g_j = Re f_j >= 0};
  - every finite point with all g_j >= 0 and some g_j = 0 has L >= 0.

  Then every point with all g_j > 0 has L >= 0. The same holds after dropping coordinates whose real part is identically 0.
- **Lemma 4 (lem:safe, p. 3)**: if x + C kappa lies in X, then Phi_O(Re x) >= min_sigma (P Re x)_sigma for all O in O(4).
- **Theorem 5 (thm:six, p. 4)**: fix O in O(4). Suppose the support inequality Re Phi_O(x) >= min_sigma Re F_sigma(x) holds for every nonzero x in X with at least 6 contacts. Then it holds for every x in X. This is paper Thm 3.6 in contrapositive form.
- **Corollary 7 (p. 5)**: given the six-contact hypothesis for every O, MOC(4) follows for all normal A and B.
- **Proposition 8 (prop:inventory, p. 5)**: the 6-subsets fall into 214 G-orbits, and exactly 90 of them have a non-minor generator. Positivity on those 90 families, for all O, suffices. The listed counts are 134596, 214, 113, 90, 9 and 81, and the distribution of non-minor generators per orbit is 51/35/3/1.
- **Proposition 10 (prop:one-safe, p. 6)**: every x in X whose contact set contains T = {0,1,2,3,6,7} satisfies (*) for all O (and for U(4)). This covers additional contacts and every spectrum, including a_i + b_j = 0 and repeated eigenvalues (the latter via the one-circle theorem). Checked: T equals paper T_1 exactly (identity group element).
- **Proposition 13 (prop:T1-safe, p. 7)**: the same conclusion for T = {0,1,2,3,6,8}. Outside the cited degenerations, c(x) is a nonnegative combination of minor vectors. Checked: T equals paper T_2 exactly.

Both case statements are at least as strong as (*). Contacts are non-strict, the result holds for all O in O(4), and all spectra are covered. Repeated within-spectrum eigenvalues rest on paper Prop. 1.4(i) / Thm C.1. A zero contact product, and any other a_i + b_j = 0, are handled inside the document.

## B. Dependencies

| ref | kind | used for | status |
|---|---|---|---|
| Paper Sect. 3 / App. A, (M1): GIT realization of nonzero points of X by finite spectra with lambda != 0 | core-doc | Thm 5 applies the hypothesis at curve points; Prop 10 changes coordinates | CITED (standard GIT surjectivity plus a projective coordinate change; PROOF-READ as a sketch) |
| dim X = 6, irreducible | core-doc | fiber dimension in Thm 5 | CITED (lib Jacobian rank 6, PROOF-READ upper bound) |
| Paper Thm 1.2 (two-orthogonal) | core-doc | Cor 7; "every unitary response" in Props 10 and 13 | CITED |
| Paper Prop 1.4(i) / Thm C.1 (one spectrum on a generalized circle) | common-lemma | repeated within-spectrum points in Props 10 and 13 | CITED. Statement read on p. 42 of the paper. Its hypothesis (at most 3 distinct points lie on a generalized circle) matches. |
| Paper Prop 4.1 (82 facets complete) + Lemma 4.2 | core-doc | Prop 8; the generator sets in Props 10 and 13 | CITED (lib exact certificate); cone generation PROOF-READ |
| Paper (4.2): Phi_{n_IJ}(O) = det O[I,J]^2 | core-doc | minor-only cones are nonnegative | CITED (lib, PROOF-READ there) |
| Group equivariance (4.7) | core-doc | Prop 8 orbit reduction | CITED (lib tests) |
| Birkhoff-von Neumann | external | Lemma 4 | standard |
| Removable singularities of bounded subharmonic functions; maximum principle | external | Lemma 2 | standard |
| An irreducible curve through two points of an irreducible variety | external | Thm 5, curve choice | standard (auditor-supplied justification) |

No other contact case and no seven-element set is cited.

## C. Steps

| # | claim | kind | status | how |
|---|---|---|---|---|
| 1 | (M1): P has rank 14, P kappa = 1, and (P m)_sigma = prod (a_i + b_sigma(i)) symbolically | identity | EXACT | own sympy, 8 symbols |
| 2 | (M2): grad C(kappa) = 0; H_kappa is (3.10); (2) holds; grad C vanishes on the symbolic matching vector | identity | EXACT | sympy |
| 3 | (M3): the e_ij span rank 10, rank H_kappa = 4, H_kappa e_ij = 0, and ker H_kappa = span e_ij | identity | EXACT | Fraction linear algebra |
| 4 | (M4): Phi_O(kappa) = 1; Phi_O(e_ij) = O_ij^2; lambda det(D_a + O D_b O^t) = Leibniz functional. For det O = -1, use O diag(1,1,1,-1), which gives the same determinant. | identity | EXACT | 8 exact rational SO(4) matrices with Gaussian-rational spectra |
| 5 | Lemma 2: finite boundary of D has L = 0; u = max-extension is subharmonic and bounded; punctures are removable; u is constant, so D is everything; the g_j lie in [0,1], so there are no poles and the f_j are constant | analytic | PROOF-READ | every inference checked; the sub-mean-value property at boundary points uses u >= 0 |
| 6 | Lemma 4: a line in X forces H_kappa x = 0 by (2). Then Re x = sum c_ij e_ij, and Birkhoff gives Phi_O(q) = sum_sigma t_sigma (Pq)_sigma >= min. | identity + analytic | PROOF-READ (identities EXACT in steps 2 and 3) | |
| 7 | The common kernel of the gap forms is C kappa, and h_{O,tau} is in their real span. So L is bounded on the simplex. | identity | EXACT | tau = 0, 7, 23; random rational O |
| 8 | Thm 5 normalization: the all-equal case, and the positive real scaling (4) | elementary | PROOF-READ | |
| 9 | Fiber of r <= 5 affine cuts: each component has dim >= 6 - r >= 1 | dimension | PROOF-READ | affine dimension theorem; X irreducible of dim 6 |
| 10 | A constant gap map on Y forces Y = x + C kappa, and Lemma 4 then gives a contradiction | algebraic | PROOF-READ | |
| 11 | Choice of an irreducible curve C in Y through x with a nonconstant gap | algebraic | PROOF-READ (argument supplied by auditor) | Pick y in Y with gap(y) != gap(x) and join x to y by an irreducible curve in Y. The document's "generic hyperplane sections" also works. |
| 12 | Apply Lemma 2. On C the contact gaps are identically 0 and their sum is 1. A finite point with a new zero has >= r+1 contacts, is nonzero (s_tau != 0) and is realizable. Induction gives L >= 0 there. | analytic | PROOF-READ (relies on CITED (M1)) | |
| 13 | Cor 7: two-orthogonal convex combination plus separation in R^2 | citation | CITED (Thm 1.2) + PROOF-READ | |
| 14 | Prop 8: cone generation by facet profiles; minor-only orbits are nonnegative | polyhedral | PROOF-READ + CITED (Prop 4.1, (4.2)) | |
| 15 | Prop 8 counts: 214 orbits (sizes sum to 134596); 90 non-minor; distribution {1:51, 2:35, 3:3, 4:1}; 113 exact six-vertex-face orbits; all 90 are six-vertex faces; 9 retained orbits contain an image of {2,7,19,20}; T_1..T_90 lie in distinct retained orbits; non-minor profiles are 13/10/1 | finite | EXACT | own enumeration with the lib group |
| 16 | Prop 10/13, degenerate case: a repeated within-spectrum point is covered by Thm C.1 | citation | CITED | |
| 17 | Zero contact product: the whole entry face {sigma(i) = j} becomes contacts. All 16 entry faces have exactly 15 generators, all minor (6 singleton + 9 two-by-two). | finite | EXACT + PROOF-READ | |
| 18 | Factorizations of D06, D17, D02, D13, D01 and D23, the three displayed ratios, and nonvanishing of each factor under the nondegeneracy assumptions | identity | EXACT (+ PROOF-READ for nonvanishing: each factor is a within-spectrum difference or a factor of a contact product F_0, F_1, F_2) | sympy |
| 19 | A real cross-ratio of 4 distinct points means they are concyclic. The Mobius normalization keeps all F_sigma (lambda absorbs the common factor), so x is the same point of X and Phi_O is unchanged (P is injective). | geometric | PROOF-READ | |
| 20 | Imaginary numerators xy(x-y)v and -tu(u-t)q | identity | EXACT | sympy with real symbols |
| 21 | Collision u = 0: D01 = -L beta t(x-y) and D02 = -L xy(beta - t), so Re L = 0 and then q = 0 | identity + elementary | EXACT + PROOF-READ | |
| 22 | All six points real and Re L = 0, so the products with sigma(0) = 0 are L times a real number. Hence {0..5} are contacts, an entry face. | identity | EXACT + PROOF-READ | |
| 23 | Prop 13: identity (6) D79 D02 = D68 D13; D02 != 0, so D79 is imaginary and g7 = g9 | identity | EXACT + PROOF-READ | Only D02 != 0 is needed. The document's further claim that D68 and D13 are nonzero is true (their factors come from F6, F1 and within-spectrum differences) but unused. |
| 24 | J = {0,1,2,3,4,5,6,9,10,15,16,17} in the JSON numbering. It equals lib generators(T_2) by profile, and only f17 is non-minor. | finite | EXACT | |
| 25 | f17 + f9 = f4 + f10 + f16; l values -1, 1, -1, 0; l(g) = 0 iff a9 = a6 + a17; the formal coefficient transfer | identity | EXACT | sympy with formal a_j |
| 26 | The resulting combination of minors is nonnegative for every O | citation | CITED ((4.2)) | |

## D. Falsification (NUMERIC)

The library samplers were used, with seeds as printed in the output file and tolerance 1e-9.

- **T_1, 12 stratum points.** Every sample is degenerate: within-spectrum distance 0, as the document predicts. Its Prop 10 argument forces {4,5} to be contacts in the nondegenerate case. The worst min Phi over SO(4) is -2.0e-10, which is 0 to working precision (12 BFGS restarts).
- **T_1 with contact set forced to be exactly T_1 (margin 1e-2), 6 points.** All six have a repeated within-spectrum point, again as predicted.
- **T_2, 12 points.** The worst min Phi is -3.1e-15. Points with g7 != g9 all have a zero product a_i + b_j and extra contacts that contain an entry face. That is the degenerate branch, so it is consistent with the document.
- **T_2, 2 nondegenerate points (seeds 33 and 46).** |g7 - g9| is at most 9e-16, as predicted by (7).

No violation was found.

## E. Verdict

**confirmed-modulo-citations.** I found no gap.

Theorem 5 (six-contact reduction) is correct as written, given the CITED realization (M1). The one step that is only sketched, the curve choice (step 11), is closed by the standard two-point curve lemma.

Props 10 and 13 cover paper T_1 and T_2 completely, for all O and all spectra. The only external case they need is Thm C.1, for repeated eigenvalues within a spectrum.

Presentation remarks, none of which counts as a gap:
- The document labels paper T_2 as "T1".
- Its f_j numbering is the supplement's JSON order, not paper App. D.
- The nonvanishing of D68 and D13 is asserted without reason, but it is true and not needed.
- Remark 6's statement about strictly semistable configurations is not used.
- p. 3 and p. 4 contain the typo "configurationsl".
