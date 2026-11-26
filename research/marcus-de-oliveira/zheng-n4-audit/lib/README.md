# moc4lib: shared exact library for the Zheng n = 4 audit

Evidence level: EXACT for every test whose name does not start with `numeric_`; NUMERIC for the `numeric_` tests (seeds and sample sizes are printed in `test_output.txt`).

All objects are rebuilt from their definitions in the paper (Sections 3.1-3.3, 4.1-4.2, Appendices A, D, E, G.1). Nothing from the supplement was executed, imported or read as data. The tables (A.1), (A.2), A.2, Appendix D, Appendix E, Table 3 and Table 5 were transcribed from the page images of the paper PDF (pp. 31, 35-36, 59, 60, 62-63) and cross-checked against the text extract.

Run: `python -I -B test_moc4lib.py` from this directory (Python 3.14, sympy 1.14, numpy, scipy). It runs in about 13 s and exits 0 iff all checks pass. The last run is in `test_output.txt`: 77/77 checks passed.

## Deviations and limits

- `P` is computed by exact linear algebra on polynomial coefficients rather than by bracket straightening. Every row is then checked coefficient by coefficient as a polynomial identity in the sixteen variables of eight symbolic 2-vectors. Eq. (4.5) is checked again, independently, with sympy in the dehomogenized variables.
- Prop. 4.1 (exactly 82 facets) is covered by an exact check of a certificate. qhull, in floating point, proposes the boundary simplices (`scipy.spatial.ConvexHull`), and every property the argument uses is then re-checked exactly. The paper's own 6846-simplex tables were not used. See the claims table.
- The test `eq4.1_minor_vectors_exist_all_24_rows` passes by construction once the solves return, because `solve_P` raises on any inconsistent row and a crash counts as a failure.
- The numeric helpers are best effort. `min_phi_over_SO4` returns an upper bound for the minimum, not a certified minimum. If `random_point_on_contact_stratum` returns None, that says nothing about whether the stratum is empty.
- Not done here: the 90 case inequalities, Theorem 3.6, Appendices B and C, and HMSV's ideal-generation statement (X = {grad C = 0} is CITED [15]; only X ⊂ Sing C is replayed).

## Conventions

- **Permutations.** `PERMS[j]` is a 0-based tuple with `sigma[i]` the image of `i`. The list is in lexicographic order of one-line notation, so `PERMS[0]` is the identity and `PERMS[j]` is the paper's `pi_j`, shifted down by one. A matrix `O` enters a sum through `prod_i O[i][sigma[i]]`, which matches the paper's `P_sigma` with a 1 in row i, column sigma(i). `compose(p, q)` is p o q, so `compose(p, q)[i] = p[q[i]]`. `perm_from_oneline("2143") = (1,0,3,2)`.
- **Brackets.** Labels run 1..8 and `[ij] = x_i y_j - x_j y_i`. An element of R_1 is stored as `{mask: int}`, where bit l-1 set means the monomial contains x_l rather than y_l. This is the exact 16-variable polynomial, and setting y = 1 gives the paper's `[ij] = z_i - z_j`.
- **V_R.** A vector c is the list of its 14 coordinates m_0..m_13, with `m_j = c(B_j)`. `values(c)` = P c, which equals (v_sigma(c)) in `PERMS` order. `KAPPA = e_13`.
- **Index sets.** I and J are 1-based tuples. For 2-sets the canonical key has 1 in I, as in Appendix D.
- **Group G.** An element is `g = (rho, rho2, eps)`. It acts on permutations by `g.sigma = rho o sigma o rho2` when eps = 0, and by `rho o sigma^{-1} o rho2` when eps = 1. This is the paper's A_{rho,rho',eps} (Sect. 5.2). It is the only one of the three conventions tried that reproduces all 7 rows of Table 3: the reversed order matches 2 rows and the inverted one matches 6.
  - On matrices: `(g.O)[i][j] = O[rho2(i)][rho^{-1}(j)]`, i.e. R_{rho'}^{-1} O R_rho^{-1}, and eps = 1 first transposes O.
  - On V_R: `g.c` is defined by `v_{g.sigma}(g.c) = v_sigma(c)`.
  - On spectra: `g_act_spectra` gives a' = a o rho2 and b' = b o rho^{-1}, with a and b swapped first when eps = 1.

## Functions

| name | returns |
|---|---|
| `PERMS, SGN, perm_index, sgn, compose, inverse, perm_from_oneline, perm_to_oneline` | permutation helpers (above) |
| `rref, rank, det, inv, nullspace, matmul, matvec, transpose, eye, det_leibniz4` | exact Fraction linear algebra; `det_leibniz4` works over any ring |
| `QI(re, im)`, `re_part(x)` | exact Gaussian rationals, real part |
| `MATCHINGS, B_POLYS` | the 14 matchings of (A.1) and their 16-variable polynomials |
| `bracket, matching_poly, poly_add, poly_relabel` | R_1 polynomial helpers |
| `express_in_basis(poly)` | exact coefficients in B_0..B_13. It raises unless the identity holds on every monomial, so a return is a proof |
| `p_sigma_poly(sigma)`, `P` | p_sigma, and the 24 x 14 integer matrix (entries 0, 1, 2) |
| `P_PIVOT_ROWS`, `KAPPA`, `values(c)`, `solve_P(vals)` | Appendix A rows; kappa; P c; the unique preimage, checked on all 24 rows (raises if the vector is not in im P) |
| `cubic_C(m)`, `cubic_sympy()`, `C_PARTIAL_MONOMIALS`, `transposition_matrix(i)`, `hess_kappa()` | (A.2); T_i with m(t_i v) = T_i m(v); the claimed H_kappa |
| `E_coords(Q)` | coordinates of E(Q), Appendix A.2 (Q is 0-based 4x4) |
| `n_values, n_vec(I,J), MINOR_KEYS, canonical_minor_key, t_vec(pi), d_vec(pi), h_vec(pi,eps)` | Section 4.1 vectors (Fractions) |
| `F_TABLE, F_VECS, F_VALUES, F_ZERO, F_IS_MINOR, f_name(j)` | Appendix D: f_j specs, vectors, 24 integer values, zero sets, minor flag |
| `leibniz_weights(O), Phi(c,O), Phi_values(vals,O), Phi_O4(c,O)` | (3.4); the O(4) extension of Sect. 4.2 |
| `submatrix, minor_det, O_STAR, cayley(S), plane_rotation, random_rational_SO4(rng), is_special_orthogonal` | exact rational SO(4) tools |
| `matching_vector(a,b,alpha)`, `matching_vector_from_points(pts)`, `eval_poly`, `eval_poly_hom` | m with (P m)_sigma = alpha prod_i (a_i + b_sigma(i)), using v_i = (a_i,1) and v_{4+j} = (-b_j,1) |
| `sigma_points(a,b)`, `gap_vector(m, pi_star)`, `contact_set(c)`, `det_Da_plus_ODbOt(a,b,O)` | (3.1), (3.3), the contact set (raises if P c has a negative entry), the determinant |
| `G_ELEMENTS, g_act_perm, g_index_map, g_act_values, g_act_c, g_act_O, g_act_spectra, g_act_set, orbit_rep(T), find_g(T,T2)` | the group of order 1152 |
| `CONTACT_SETS[nu]`, `AUX7[j]`, `TABLE3`, `generators(T)`, `enumerate_six_set_orbits()` | Appendix E, Table 5, Table 3; F_T; the 214 orbits with their lexicographically first representatives |
| NUMERIC: `P_numpy, leibniz_weights_np, Phi_np, so4_from_params(theta), min_phi_over_SO4(c, restarts, seed), F_values_np, matching_vector_np, random_point_on_contact_stratum(T, seed, tries, tol, pi_star, margin)` | floating point. SO(4) = expm(skew(theta)). The minimizer is multistart BFGS. The stratum search uses least squares in (a, b, phase of alpha) with mean off-T gap = 1, and `margin > 0` forces the contact set to be exactly T |

## Paper claims replayed (all EXACT unless marked)

| claim (paper location) | result |
|---|---|
| (A.1) B_0..B_13 are exactly the 14 noncrossing matchings and are independent; all 105 matching products lie in their span | pass |
| P: integer, rank 14, p_sigma = sum P B as a polynomial identity in 8 symbolic 2-vectors; (4.5) re-checked with sympy | pass |
| App. A: the minor on rows 0,1,2,3,4,6,7,8,9,10,12,13,16,18 has det **+1** | pass |
| (3.2) P e_13 = 1; kappa in X (columns (0,1)^4, (1,1)^4) | pass |
| A.1: grad C vanishes identically on X (16 variables); T_i^2 = I and C(T_i m) = -C(m) for i = 1..7; the 14x14 coefficient minor of the partials on the listed monomials has det **+1** (rows are dC/dm_0..dC/dm_13) | pass |
| (3.10) Hess C(kappa) = e0e6^t + e6e0^t - e1e5^t - e5e1^t, rank 4; (3.11) its kernel = {m0=m1=m5=m6=0} = E; (3.12) grad C(kappa) = 0 and the polarization identity | pass |
| dim X = 6: the Jacobian of (lambda, z) -> lambda B(z) has rank 6 at three rational points (an upper bound of 6 follows from the PGL_2 invariance, PROOF-READ) | pass |
| (3.7)/A.2: (P E(Q))_sigma = sum_i q_{i,sigma(i)} symbolically; dim E = 10 | pass |
| (4.1): n_{I,J} exists for all 16 + 36 (I,J); n_{I,J} = n_{I^c,J^c}; 16 + 18 = 34 distinct vectors; 18 and 20 zeros respectively | pass |
| Sect. 4.1: t_pi counts agreements; d_pi is nonzero exactly at the 8 sigma that agree with pi in exactly one place, four +1 and four -1; h has profile 13/10/1 with its 2 at pi; the 48 h are distinct | pass |
| App. D: the 82 f_j are exactly the 34 minor vectors and the 48 h; profiles correct; f_j >= 0 at all vertices; zero sets pairwise distinct; each zero set has affine rank 12, so each f_j defines a facet | pass |
| Prop. 4.1 completeness, via an exact check of a qhull certificate: 9420 boundary 12-simplices, each inside a listed zero set; their mod-2 boundary is 0; a rational point inside facet f_17 lies in the interior of exactly 1 simplex and on the boundary of none. The class is then nonzero in H_12(dQ; Z/2), so the support is all of dQ and no facet is unlisted | pass (EXACT check; qhull only proposes the simplices) |
| (4.2) Phi_{n_{I,J}}(O) = det O[I,J]^2 for all 52 (I,J), on 24 random exact rational SO(4) matrices | pass (see the proof below for general O) |
| Phi_kappa = 1; (3.8) Phi_{E(Q)} = sum q_ij O_ij^2; (4.4); Phi at a signed permutation matrix with det 1 equals (Pc)_sigma | pass |
| (4.3) at O*: Phi_t = 4/49, Phi_d = 512/2401, Phi_{h_id,-} = -158/2401; O* is in SO(4); Example 1.3 multiplicities 13, 10, 1 | pass |
| (3.1) (P m)_sigma = alpha prod (a_i + b_sigma(i)), and (3.3), on Gaussian-rational spectra; (3.5) Phi_{c(x)}(O) = Re(alpha det(D_a + O D_b O^t) - F_pi*) on 12 Gaussian-rational spectra; grad C = 0 at a Gaussian-rational matching vector | pass |
| G: 1152 distinct permutations of the 24 indices; all 7 rows of Table 3, each image contained in T_85 ∪ {j}; (4.7) w_{g.sigma}(g.O) = w_sigma(O) and Phi_{g.c}(g.O) = Phi_c(O) for all 1152 g at one rational point, plus 20 random; the column flip for det -1; G maps matching vectors to matching vectors of relabelled spectra; G permutes {f_j} and permutes the minor vectors | pass |
| Sect. 4.2: the 6-subsets form **214** orbits (sizes sum to 134596); in **124** all generators are minor and in **90** they are not; T_1..T_90 of App. E lie in exactly those 90 orbits, one each | pass |
| Table 5: 27 distinct 7-sets, all with a non-minor generator (INFO only, not a paper claim) | info |
| NUMERIC: min Phi_{n_12,34} >= -1e-10 (10 restarts); min Phi_{h_id,-} = -0.067442 < Phi(O*) = -0.065806 (20 restarts, seed 1); a floating-point point with contact set exactly T_85 (seed 3), and one with contact set exactly T_1 (margin 1e-2) | pass |

No discrepancy with the paper was found in any item above.

## Proof of (4.2) for general O in SO(4)

Fix I and J with |I| = |J| = k. By definition, (P n_{I,J})_sigma = 1 when sigma(I) = J and 0 otherwise, so

Phi_{n_{I,J}}(O) = sum over sigma with sigma(I) = J of sgn(sigma) prod_i O_{i,sigma(i)}.

Each such sigma is a pair: a bijection I -> J and a bijection I^c -> J^c. Its sign is (-1)^{s(I)+s(J)} times the signs of the two bijections, where s(I) is the sum of the elements of I and each bijection's sign is computed in the increasing orders. This is the sign rule of the generalized Laplace expansion. So the sum factors as

Phi_{n_{I,J}}(O) = (-1)^{s(I)+s(J)} det O[I,J] det O[I^c,J^c].

Jacobi's complementary-minor identity says that for invertible O,

det(O^{-1})[J,I] = (-1)^{s(I)+s(J)} det O[I^c,J^c] / det O.

For O in SO(4), O^{-1} = O^t and det O = 1, and det(O^t)[J,I] = det O[I,J]. Hence

det O[I^c,J^c] = (-1)^{s(I)+s(J)} det O[I,J],

and substituting gives Phi_{n_{I,J}}(O) = det O[I,J]^2 >= 0. For k = 1 this is the statement that the cofactor matrix of O equals O. The same computation, applied entry by entry, proves (3.8):

Phi_{E(Q)}(O) = sum_{i,j} q_ij Phi_{n_{i,j}}(O) = sum_{i,j} q_ij O_ij^2.

This works because P E(Q) = sum_{i,j} q_ij P n_{i,j} and P is injective. (PROOF-READ; the tests replay it on 24 exact rational matrices.)

## Facts used by the exact Prop. 4.1 certificate

Let z be the mod-2 sum of the affine 12-simplices on the 13-vertex sets returned by qhull. The test checks exactly that each simplex lies in the zero set of a listed f_j, which is a face of Q because f_j >= 0 on all vertices. So supp z is contained in the union of the listed facets, which lies in dQ, a 12-sphere. The test also checks that every 12-vertex face occurs an even number of times, so dz = 0 mod 2. If supp z missed a point q of dQ, then [z] would lie in the image of H_12(dQ \ {q}; Z/2) = 0. The local degree mod 2 at every point off supp(dz) would then be 0. But at the chosen generic rational point p, the local degree is the number of simplices whose relative interior contains p, which is 1. Hence supp z = dQ. A facet outside the list would have a relative interior point that lies on no listed facet, a contradiction.
