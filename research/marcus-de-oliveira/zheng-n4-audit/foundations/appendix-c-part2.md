Evidence level: EXACT replays of every displayed identity and rational witness in Sections C.4-C.7 and of (C.3), (C.5), (C.6), (C.8) on samples; PROOF-READ for Lemmas C.8-C.11 and the final assembly; CITED for the inputs listed in Section 2; NUMERIC for the attacks on Lemma C.10 and Theorem C.1.

# Appendix C, part 2: Lemmas C.8-C.11, the proof of Theorem C.1, Section C.7

Paper: Zheng, "The Marcus-de Oliveira conjecture for n = 4", Appendix C, pages 50-58
(Sections C.4-C.7 and "Proof of Theorem C.1"). Parallel supplement text: `one_circle.pdf`
(22 Oct 2026, 12 pages), which is the same argument with different names (F for d, kappa for
alpha, Lemma 4-8 for Lemma C.7-C.11, r(U) for rho(U), e for kappa).

Own code: `appC2_checks.py` (this directory); captured run: `appC2_checks_output.txt`.
Nothing from the supplement was executed or imported. The code rebuilds the matching basis
(A.1), the matrix P, kappa, the cubic C (A.2), the minor vectors n_{I,J} and rho(U) from the
printed definitions.

## Deviations and limits

- (C.4) (Segre-polar identity) is in Section C.1, so it is part 1's; I did not replay it.
- Lemma C.7, Proposition C.4, Corollary C.5, Proposition C.6 and the statements (C.2) quoted from
  Howard-Millson-Snowden-Vakil (HMSV) are part 1's. I treat them as inputs (Section 2). I did
  however test (C.8), the identity C|K_O = alpha(O) det S_O on which Lemma C.9 rests, exactly at
  four orthogonal matrices, because part 1's text only argues it (see Section 4.2).
- Lemma C.10's positivity rests on Lam's Theorems 2.1 and 3.1 (Postnikov networks, double
  dimers). I read the argument against the theorem as quoted; the source was not re-read. I
  replaced the citation by an exact LP replay on samples (EXACT per sample, not a proof).
- The irreducibility of the Pluecker divisors inside the stable Gale-fixed locus (proof of
  Lemma C.9, p. 53) is only sketched in the paper. I checked the sketch for coherence but did not
  write out the proof.
- All "for all O" statements tested numerically are NUMERIC; sizes, seeds and tolerances are in
  the output file.
- The supplement's "alternate_onecircle ... elementary" note (3 pages, scalar and two-rotation
  inequalities) is not used by Appendix C. It belongs to the separate principal-minor route
  (`principal_minors.pdf`, Corollary 11), which claims Theorem C.1 by a different proof. I did
  not audit that route; I only tested its main claim numerically (Section 6.3).


## 1. Verdicts

| Item | Verdict |
|---|---|
| Lemma C.8 (p. 51-52) | Proved as written, modulo (C.9). (C.9) replayed EXACT as a polynomial identity. |
| Lemma C.9 (p. 52-53) | Proved modulo part-1 inputs ((C.8), the line-bundle section alpha, smoothness of the stable fixed locus) and one sketched step (Pluecker divisors irreducible, one S8 orbit). All its computations replayed EXACT. No error found. |
| Lemma C.10 (p. 53-55) | Proved modulo the citation of Lam [21, Thms 2.1, 3.1] / Postnikov [27]. Minor expository gap in the density step (Section 5.2). Conclusion replayed EXACT on 60 Cauchy samples (unique strictly positive weights on the 14 Lam permutations), 25 zero-entry and 20 decomposable samples; NUMERIC on 500 root-found divisor points and a 399-point walk on d = 0: no failure. |
| Lemma C.11 (p. 56) | Proved modulo (C.2) (HMSV, cited). All printed numbers replayed EXACT. |
| Final assembly, Theorem C.1 (p. 56-58) | Proved as written, modulo the cited inputs of Section 2. No missing case found. NUMERIC attacks found no violation. |
| Section C.7 (p. 58) | Every identity in my part replayed EXACT with own code (list in Section 7). |

No error was found. The appendix is not self-contained: its force rests on part 1
(Lemma C.7, (C.8), Proposition C.4) and on HMSV (C.2). The strongest independent support I can
give for those is the exact test of (C.8) in Section 4.2. The critical-point search of Section 6.1
found only decomposable critical points and indecomposable ones lying exactly on the Cauchy divisor,
so the delicate branch (stable, d != 0) of the chain was never exercised numerically; it was not
contradicted either.

## 2. What part 2 takes from part 1 and from outside

- A1 (C.2), CITED to HMSV [15, Thm 3.4, proof of Prop. 3.13, Prop. 3.19]: Sing C = M8 (cone X);
  Sigma2 = Sec(M8) is contained in {C = 0}; det H restricted to {C = 0} equals gamma0 s^2 with
  gamma0 != 0, where the septic s cuts out Sigma2 as a reduced Cartier divisor on {C = 0}. Used in
  Lemma C.11 and the uniqueness step. My exact side checks (Section 5.4) are consistent with it.
- A2 (C.3) rho(U) and the determinant pairing: replayed EXACT here (4 rational orthogonal and
  one Gaussian-rational unitary matrix, all 70 minor pairings).
- A3 (C.7), main-text Theorem 1.2 (two orthogonal matrices represent every unitary response):
  CITED to the main text. Used only to place a maximiser on O(4).
- A4 Section C.3 (text after Lemma C.7): an indecomposable real or complex orthogonal O gives a
  stable octad [I4, O]; the stable Gale-fixed locus is smooth of dimension 6, M8-bar is smooth of
  dimension 9 there with tangent eigenspaces of dimensions 6 (+1) and 3 (-1); the 3 anti-invariant
  directions are the unitary phase directions modulo row/column phases. PROOF-READ by part 1;
  I confirmed EXACT at four matrices that the gauge space {D1 + O D2 O^t} has dimension 7 and
  that it is null for the phase Hessian (Section 4.2).
- A5 Lemma C.7 (conormal limit) and Proposition C.4 ((M8')* = Sigma2): part 1. Used to get
  [c] in Sigma2.
- A6 (C.8): C|K_O = alpha(O) det S_O, dim K_O = 7, S_O : K_O -> Sym^2 surjective, alpha a regular
  section of a line bundle on the stable fixed locus with chart-independent zero set. Part 1
  argues this; I tested it EXACT at four matrices (Section 4.2).
- A7 Lemma C.3 (Drury-Cload zero-determinant case): part 1. Used in the last line.
- A8 Lemma C.2 / Birkhoff for blocks of size <= 3: used for decomposable maximisers and in
  Lemma C.10. Replayed EXACT on 20 decomposable samples.
- A9 Projective biduality (reflexivity in characteristic 0): standard, CITED.
- A10 Lam, "Dimers, webs, and positroids", Theorems 2.1 and 3.1, and Postnikov's network
  theorem: external citation, hypotheses checked against the statement as quoted (a totally
  nonnegative point of Gr(4,8); complementary boundary sets), source not re-read.

## 3. The final argument for Theorem C.1, step by step

Notation: z = (a1..a4, -b1..-b4), m = m(z) in X, theta in C, theta m = y + i w with y, w real,
f(U) = <rho(U), y> = Re(theta det(D_a + U D_b U*)) by (C.3).

| # | Step | Status |
|---|---|---|
| 0 | Reduce to real spectrum a: a Moebius map f with pole on the circle (off the 8 parameters) sends the circle to R; m(f(z)) = (ad-bc)^4 / prod(gamma z_j + delta) m(z), so the determinant value (via (C.3), for every U) and all 24 sigma-points are multiplied by one nonzero scalar and hull membership transfers. Swap A and B if b is the concyclic spectrum. | Identity EXACT (symbolic); argument PROOF-READ. Correct. |
| 1 | Fix (a real, b, theta) in the dense open set of Lemma C.11. | Lemma C.11, see 5.4. |
| 2 | max over U(4) of f(U) exists (U(4) compact, f continuous) and by (C.7) is attained on O(4): f(U) is a convex combination of f(O-), f(O+). A row sign change keeps rho, so O in SO(4). | Existence PROOF-READ; reduction CITED (A3). The extremum is over U(4); SO(4) suffices because squared minors ignore row signs. |
| 3 | Case O decomposable: blocks of size <= 3 have Birkhoff weights representing all squared minors; the product weights give rho(O) = sum p_sigma (sigma-point functional), hence h = f(O) <= max_sigma Re(theta z_sigma). | PROOF-READ; EXACT on samples. |
| 4 | Case d(O) = 0: Lemma C.10 gives the same conclusion. | Lemma C.10, see 5.1-5.3. |
| 5 | Otherwise O is indecomposable, hence [I4, O] stable (A4), and d(O) != 0. Put h = f(O), c = y - h kappa. Real criticality on SO(4) extends C-linearly to SO(4, C); <rho(O), c> = 0; so c in K_O. | PROOF-READ. Only first-order criticality is used from here on. |
| 6 | If c = 0: y = h kappa, so Re(theta z_sigma) = (P y)_sigma = h for all sigma; the inequality holds with equality. | PROOF-READ (P kappa = 1 EXACT). |
| 7 | If c != 0: C(c) = C(y) - h q(y) = 0 because y in Sigma2-hat subset {C = 0} (A1) and q(y) = 0 by (C.6). | (C.5), (C.6) EXACT (symbolic in all 17 real parameters); A1 CITED. |
| 8 | Lemma C.9 gives alpha(O) != 0, so (C.8) gives det S_O(c) = 0. | Lemma C.9 see 4.2; (C.8) is A6 (tested EXACT at 4 points). |
| 9 | Lemma C.7: ([rho(O)], [c]) in Con(M8'), so [c] in (M8')* = Sigma2 (Prop. C.4). | Part 1 (A5). Hypotheses of Lemma C.7 match K_O and S_O (checked against the text). |
| 10 | c = y - h kappa and y both lie on Sigma2-hat and on the line y + C kappa; Lemma C.11 says y is the only such point, so h = 0, c = y. | PROOF-READ given Lemma C.11. |
| 11 | Biduality: ([y], [rho(O)]) in Con(Sigma2); y is a smooth point of Sigma2-hat, so rho(O) annihilates T_y Sigma2-hat. | PROOF-READ; A9 CITED. Over a smooth point the fibre of the closed conormal variety is exactly the tangent annihilator, so the closure in the definition causes no trouble. |
| 12 | y + i w = theta m and y - i w = conj(theta m) lie in X (X is defined over R), so their span lies in Sigma2-hat; y + t w lies on it for real t, so w in T_y Sigma2-hat and <rho(O), w> = 0. With h = 0: <rho(O), theta m> = 0, i.e. theta det(D_a + O D_b O^t) = 0, theta != 0. | PROOF-READ. |
| 13 | Lemma C.3: 0 is in the hull, so max_sigma Re(theta z_sigma) >= 0 = h = max_U f. | A7 (part 1). |
| 14 | The support inequality holds on a dense open set of (a, b, theta); both sides are continuous (compactness of U(4); finite max), so it holds everywhere; all directions theta give hull membership. | PROOF-READ. The conclusion is a closed condition, so the limit is legitimate. |

Points the task asked me to watch:

- Existence and regularity of the extremum: U(4) is compact; O is a genuine maximiser on U(4)
  (step 2). Nothing about non-degeneracy is assumed: steps 5-12 use only that O is a critical
  point. Degenerate or non-isolated critical points therefore cause no difficulty.
- Local versus global: the argument is local (first order at O) but O is the global maximiser,
  and the conclusion h = 0 is a statement about the global maximum. Correct. Remark: the paper's
  overview (p. 43) speaks of second derivatives "at a maximum", but maximality enters only through
  step 2; the phase Hessian is a device to reach the conormal variety. The chain therefore implies
  a much stronger, checkable statement: on the dense open set, every real critical point of f on
  SO(4) that is indecomposable with d(O) != 0 has f(O) = 0 and det(D_a + O D_b O^t) = 0. I tested
  this numerically (Section 6.1).
- Boundary strata: zero entries are handled (decomposable in step 3; indecomposable with zeros in
  Lemma C.9 via the Pluecker-divisor argument and in Lemma C.10 via (C.13)-(C.14) and density);
  unstable octads of the form [I4, O] are exactly the decomposable ones (A4); repeated eigenvalues
  and a_i + b_j = 0 are excluded from the dense open set and recovered in step 14 (closed
  condition). Nothing is missing.
- Real spectrum to circle or line: step 0. Correct, including the line case (affine map).
- Theorem C.1 does not need the generic set to be dense in any space other than (a real,
  b complex, theta complex); the transfer in step 0 is exact for every configuration.

## 4. Lemmas C.8 and C.9

### 4.1 Lemma C.8 (p. 51-52)

Statement: for real orthogonal O with all entries nonzero, d(O) = 0 iff O_ij = u_i v_j/(a_i - c_j)
with nonzero real u, v and distinct real nodes; equivalently [I4, O] lies on a smooth real twisted
cubic.

- Sign normalisation: d is invariant under every row and column sign change (EXACT, symbolic), so
  "change the sign of one row" is harmless.
- Forward: (C.9) makes all sixteen 3x3 minors of (1/O_ij) vanish (first identity, using that the
  cleared r equals det(1/R) prod R when entries are nonzero), so rank <= 2; rank 1 would make O
  singular. Pairwise independence of row and column factors: two proportional factors give two
  proportional rows or columns of O. A common GL2 change of chart makes every row factor
  lambda_i (a_i, 1) and every column factor mu_j (1, -c_j) simultaneously (finitely many bad
  charts), so 1/O_ij = lambda_i mu_j (a_i - c_j) with a_i distinct, c_j distinct and a_i != c_j.
  PROOF-READ. (The text says "second identity"; it is the first identity of (C.9) that gives
  rank <= 2. The text is slightly loose here but cites (C.9) as a whole.)
- Converse and twisted cubic: Lagrange parametrisation through the four coordinate points;
  conversely the coordinate cubics of a twisted cubic through e1..e4 vanish at three of the a_i,
  which forces the Lagrange form; the O-columns then have Cauchy form. PROOF-READ.
- Replays: 40 exact interlaced rational node sets all give a real orthogonal Cauchy matrix with
  d = 0 (EXACT; d is computed as (u1 u2 u3)^2 (v0 v1 v2)^2 det[C_rk C_sk], all rational); one
  instance with radicals checked directly (EXACT); 40 non-alternating node sets give no real
  orthogonal Cauchy matrix (EXACT), as the alternation step of Lemma C.10 predicts; 200 real O on
  d = 0 found by root-finding independently of the Cauchy form all have sigma_3/sigma_1(1/O) of
  order 1e-15 and alternating nodes (NUMERIC).
- Status: PROOF-READ, EXACT on samples. Cosmetic: the nodes are called a_i, c_j, clashing with
  the spectrum a.

### 4.2 Lemma C.9 (p. 52-53) and the identity (C.8) it rests on

Statement: for stable orthogonal [I4, O] with d(O) != 0, alpha(O) != 0.

First case (all Pluecker coordinates nonzero). Choice of xi in the 2-plane {xi_0 = 0, eta_0 = 0}:
a coordinate form xi_i vanishes on the plane only if two entries of column 0 vanish, and eta_j
only if columns 0 and j are proportional in rows 1..3, i.e. a 2x2 minor vanishes. W o O = xi eta^t,
(W o O) O^t = xi (O eta)^t = xi xi^t is symmetric, which is exactly criticality of
sum w_ij O_ij^2 under left skew perturbations; the value is xi^t xi = h. Then
C(c_W - h kappa) = C(c_W) = det(w_ij)_{i,j>=1} = prod xi_i prod eta_j det(1/O_ij) = d(O) prod xi eta
/ prod O_ij by (C.5), (C.10) and (C.9). Every step PROOF-READ; all of it replayed EXACT on 12 random
rational O in SO(4) (rational xi, exact K_O membership via the six derivatives).

Boundary (Pluecker divisors): the zero set of alpha in {d != 0} is a pure divisor (A6) contained
in the union of Pluecker divisors; those are irreducible and form one S8 orbit; alpha = 0 is
S8-invariant; the example point lies on such a divisor with alpha != 0. Logic PROOF-READ. Notes:
(a) The example matrix is a signed Householder matrix; besides O_00 it has six vanishing 2x2
minors (EXACT), so it lies on seven Pluecker divisors, not on a general point of {O_00 = 0}. This
does no harm: nonvanishing at one point excludes every irreducible divisor through it.
(b) The irreducibility and single-orbit claim is a sketch (isotropic line in a nondegenerate
coordinate 4-space, then a maximal isotropic extension: a P3-bundle over a quadric surface,
dimension 5). It is coherent and the dimensions match; I did not write out the details. Status:
PROOF-READ (sketch).
(c) Pure codimension needs alpha not identically zero on each component of the stable fixed
locus. The text does not say the locus is irreducible. It is: the chart {[I4, O]} is the image of
the irreducible SO(4, C) minus a closed set, every other chart is an S8 translate, and a generic
[I4, O] lies in all charts. Minor expository gap, easily filled.
(d) The S8 action leaves the chart [I4, O]; the argument therefore needs alpha on the whole stable
fixed locus, which part 1 provides (end of Section C.3, gluing across charts asserted briefly).

Direct test of (C.8) (EXACT, my own construction from the definitions): at a rational O I computed
K_O as the nullspace of the 7 linear conditions <rho(O), c> = 0 and its six first derivatives
(dimension 7 found), and the phase Hessian S_O(c) as the t^2 coefficient of
<rho((1 + i t Xi - t^2 Xi^2/2) O), c> on three phase directions Xi (complement of the gauge space,
whose dimension was 7, and which was checked to be null). Then C restricted to K_O, a cubic in 7
variables, was compared with det S_O:

- two random rational O in SO(4) (all minors nonzero): C|K_O = alpha det S_O with
  alpha = -27982932961/599628800 and 41537523/99008;
- the Lemma C.9 boundary example: alpha = 729/16 (= -1/d(O));
- an indecomposable point with d = 0 (row (3/5, 4/5, 0, 0)): C|K_O vanishes identically while
  det S_O does not, i.e. alpha = 0 there.

So (C.8) holds as a polynomial identity at these points, alpha is nonzero off the Cauchy divisor
there, and it vanishes at a Cauchy point, exactly as Lemma C.9 and the case split in the proof
require. Status of Lemma C.9: proved modulo A4-A6 and the sketched step (b).

## 5. Lemmas C.10 and C.11

### 5.1 Lemma C.10, all entries nonzero (p. 54)

- Lemma C.8 gives a real twisted cubic; columns A v(t_k) c_k; orthogonality I - O O^t = 0 becomes
  sum eps_k c_k^2 v(t_k) v(t_k)^t = 0 (A invertible), a relation among the 7 power rows; the
  kernel of the 7x8 Vandermonde is (1/chi'(t_k)), which alternates in sign, so coordinate and O
  labels alternate on RP^1. PROOF-READ, and EXACT/NUMERIC consistent (Section 4.1).
- After positive rescaling and sign changes the point is the row space of a Vandermonde matrix
  with increasing nodes: all maximal minors are products of positive differences, so it is totally
  positive. PROOF-READ.
- (C.12) from Lam: CITED (A10). The compatibility "each pair joins S and S^c" is the
  complementary case of the double-dimer theorem; it needs the boundary vertices of the network
  to have one colour (Postnikov's normalisation), which the text does not state. Minor.
- Alternating order: every noncrossing matching joins odd to even positions, so it is a bijection
  rows -> columns; with S = (S0 \ I) u J the compatible matchings are those with sigma(I) = J (I
  rederived this: a pair (i, j') is cut by S iff (i in I) <=> (j in J)). Normalising by
  Delta_S0 Delta_S0c > 0 gives probabilities; the Jacobi identity identifies the normalised
  products with det O[I,J]^2 (in the totally positive order all products are positive and their
  moduli agree with the original order). PROOF-READ.
- Replay: for 60 exact Cauchy orthogonal matrices (20 with two nodes 1/1000 apart), restricting p
  to the 14 permutations given by the 14 noncrossing matchings in the alternating order (the
  support the Lam argument predicts), the system "sum over sigma(I) = J of p_sigma =
  det O[I,J]^2 for all I, J; sum p = 1" has rank 14 and its unique rational solution is strictly
  positive (smallest weight 5.6e-11) and reproduces all 70 squared minors (EXACT per sample).
  Status: proved modulo A10.

### 5.2 Lemma C.10, zero entries and decomposable O (p. 54-55)

- (C.13): replayed EXACT on 12 rational O with O_00 = 0 and det O = 1 (ratio exactly 1); for
  det O = -1 the sign flips (ratio -1). The paper normalises det O = 1 first, so this is fine.
- Reduction to a row (a, b, 0, 0) with [x y z] in SO(3): PROOF-READ (a row with one nonzero entry
  is excluded by indecomposability; transposition and permutations change d only by a sign, EXACT).
- (C.14): with E_i = e_0 e_i^t - e_i e_0^t acting on the left, dd(E_i O) = - a b x_j x_k y_i z_i for
  i = 1, 2, 3 on 15 rational samples: one common sign, as stated (EXACT).
- The case analysis showing one of the three numbers is nonzero: PROOF-READ, correct.
- Independence of the differentials of two zero entries (bivectors e_i ^ v_j): PROOF-READ, correct
  in all three subcases.
- Density step (p. 55): the text says the smooth real hypersurface d = 0 "cannot locally be
  contained in any zero-entry divisor". What is needed is that it is not locally contained in the
  union of them. Fix: if it were, by Baire it would coincide with one divisor on an open piece
  near a nearby point O', which is still indecomposable and on d = 0, and the same (C.13) argument
  at O' gives the contradiction. Minor expository gap.
- Limit of weights in the compact simplex: correct.
- Decomposable O: block weights from Birkhoff represent all squared minors in sizes <= 3; product
  measures work (minors of block-diagonal matrices factor). EXACT on 20 samples.

### 5.3 Lemma C.10, numeric attack

- 500 real O with d(O) = 0 found by Brent root-finding on random geodesics (not via the Cauchy
  form): LP residual at most 5e-16 (24 permutations) and 1e-16 (reported 0.0) (14 Lam permutations).
- A random walk along d = 0: worst LP residual 0.0 (399 divisor points).
- Sanity: 57 of 300 Haar-random O in SO(4) admit no common weights, so the LP test has
  teeth.
No point on the divisor where the weights fail was found.

### 5.4 Lemma C.11 (p. 56)

- H(kappa) = e0 e6^t + e6 e0^t - e1 e5^t - e5 e1^t has rank 4 (EXACT); H(u) kappa = H(kappa) u
  (symmetry of the cubic's third derivatives); K = ker H(kappa) contains kappa and H(u)|K has kappa
  in its radical, so the h^4 coefficient of det(H(u) + h H(kappa)) vanishes. PROOF-READ.
- On {C = 0, q = 0} the whole line u + C kappa lies in {C = 0} ((C.5), EXACT), and (C.2) turns the
  degree <= 3 determinant into gamma0 s^2, so s has degree <= 1 on the line. This uses that the
  lift of s is well defined modulo C, which is fine on the line. PROOF-READ given A1.
- Example: Re m(a, -b) equals the printed y; q(y) = C(y) = 0; G(y) != 0;
  det H(y + kappa) = -114105877589169227887017984 (all EXACT). Moreover
  det(H(y) + h H(kappa)) = -114105877589169227887017984 h^2 identically (EXACT), as s(y + h kappa)
  = beta h predicts.
- Extra consistency checks of A1 (EXACT): for 6 random rational (distinct real a, nonreal b, theta)
  det(H(y) + h H(kappa)) = c h^2 with c != 0; for 6 random rational u on {C = 0, q = 0} (not of the
  form Re(theta m)) det(H(u) + h H(kappa)) has degree 2 and zero discriminant, i.e. it is a
  constant times a square of a linear polynomial, as det H|C = gamma0 s^2 requires.
- Openness and density: G(y) != 0 and det H(y + kappa) != 0 are polynomial in the real
  parameters and the example shows they are not identically zero; they imply y in Sigma2-hat (always
  true), beta != 0 (uniqueness on the line) and smoothness of {C = s = 0} at y. PROOF-READ.

## 6. Numeric attack on Theorem C.1

### 6.1 Critical points (tests the whole chain, steps 5-12)

For 30 random (a real, b complex, theta on the unit circle) and 40 random starts each
(seed 2026), I solved grad f = 0 on SO(4) (Levenberg-Marquardt in exponential coordinates,
accepted when max |grad| / max |z_sigma| < 1e-10) and classified the critical points:
1166 decomposable, 27 indecomposable with |d(O)| < 1e-9 (largest 2.3e-15), 0 indecomposable with d(O) != 0, 7 runs not converged. The chain predicts that every indecomposable critical point with d(O) != 0 has
f(O) = 0 and det(D_a + O D_b O^t) = 0. Result: the prediction is never contradicted, but also never exercised: no indecomposable critical point with d(O) != 0 was found at all. This is consistent with the chain (it allows such points only at value 0) and is NUMERIC evidence only. A striking by-product: the interior
(indecomposable) critical points found all lie on the Cauchy divisor to machine precision, which
is exactly the case Lemma C.10 handles.

### 6.2 Support-function margins

margin(theta) = max_sigma Re(theta z_sigma) - max_U Re(theta det(D_a + U D_b U*)), normalised by
max |z_sigma|, for theta along every outward edge normal of the hull of the 24 sigma-points plus
random directions, maximising over SO(4) (BFGS, 2 starts, one seeded at the best permutation
matrix) and, for two families, over U(4) (16 parameters). Families: real a with Gaussian b;
clustered a (spread 1e-3); nearly collinear b (imaginary parts 1e-3); modulus ratios up to 1e4;
a real in two tight pairs with b on the unit circle; a on the unit circle searched over U(4); real
a searched over U(4). Total 1148 optimisations (seed 4). Worst relative margin -1.6e-13 (rounding level; the maxima sat at permutation matrices).
No violation found (negative values of order 1e-16 are rounding).

### 6.3 A stronger form: one weight vector for all b

For real a (and for a on a circle) the supplement's principal-minor route claims one probability
vector p, depending on U and a, with det(D_a + U D_b U*) = sum_sigma p_sigma prod_i (a_i + b_sigma(i))
for every b. This implies Theorem C.1. Comparing coefficients of prod_{j in J} b_j gives a 16-row
LP. Result: residual at most 9.8e-10 on 300 Haar U(4) with real a (one third clustered), at most
1e-16 (reported 0.0) on 150 with a on a circle (NUMERIC). Sanity: for the explicit orthogonal matrix O0 of
`principal_minors.pdf` and the non-concyclic a = (0, 1, i, 2) the LP is infeasible (residual
9.87e-02), while for real a it is feasible; the LP therefore has teeth.

## 7. Section C.7 (p. 58): computer-verified identities

| Identity | Owner | My replay | Supplement script (name only, not run) |
|---|---|---|---|
| (C.4) Segre-polar identity | part 1 | not replayed here | `root_moc_universal_segre_polar_20261022.py` |
| (C.9) both identities | part 2 | EXACT polynomial identity in the six quaternion parameters (flint): all 16 signed cofactors equal d (degree 12, 256 terms); det(1/T) prod T = d^2 (degree 24); extended to O(4, C) by sign-flip invariance | `root_arithmetic_moc_reciprocal_universal_20261022.py` |
| (C.10) | part 2 | EXACT in 16 symbolic weights | `root_arithmetic_moc_entry_cost_witness_20261022.py` |
| Lemma C.9 example (critical point, h = 18, d = -16/729, C = -531441/200) | part 2 | EXACT | same script |
| (C.13) | part 2 | EXACT on 12 rational samples (det 1); sign -1 on det -1 | `fresh_moc_cauchy_real_boundary_20261022.py` |
| (C.14) | part 2 | EXACT on 15 rational samples, common sign -1 | same script |
| Lemma C.11 example (y, q = C = 0, G != 0, det H(y + kappa)) | part 2 | EXACT, plus det(H(y) + h H(kappa)) = const h^2 | `root_one_circle_radial_witness_20261022.py` |
| (C.8) | part 1 (no script; the supplement says a full elimination was attempted and dropped) | EXACT at 4 matrices (Section 4.2) | none |
| Lemma C.10 weights | part 2 (no computation claimed) | EXACT on 105 samples, NUMERIC on about 900 divisor points | none |

Also replayed (EXACT): (A.2) equals the 28-term expansion printed in `one_circle.pdf`; grad C
vanishes identically on m(z); (C.5); (C.6) symbolically in all 17 real parameters; the Moebius
scaling of matching vectors. Observation: the docstring of `root_one_circle_radial_witness_20261022.py`
calls the MOC proof "still incomplete"; this looks like a leftover from an earlier stage and
does not point to a specific gap.

## 8. List of issues (none is an error)

1. p. 51-52, proof of Lemma C.8: rank <= 2 of (1/O_ij) comes from the first identity of (C.9)
   (vanishing 3x3 reciprocal minors). Cosmetic: the nodes a_i clash with the spectrum.
2. p. 53, Lemma C.9: irreducibility and single S8 orbit of the Pluecker divisors are sketched; the
   irreducibility of the stable fixed locus, needed for purity, is not stated (easy, Section 4.2 c).
   The example lies on seven Pluecker divisors at once (harmless).
3. p. 54, (C.12): the one-colour boundary hypothesis behind "each pair joins S and S^c" is
   implicit in the citation of Lam.
4. p. 55, density step: "not contained in any zero-entry divisor" should read "not contained in
   their union"; a Baire argument fills it (Section 5.2).
5. p. 43 overview vs p. 56-57 proof: the second-order (maximum) condition is never used; only
   first-order criticality. Not an error; it makes the argument prove more (tested in 6.1).
6. p. 55, (C.13) holds for det O = +1 only (sign -1 on the other component); consistent with the
   normalisation in the proof.
