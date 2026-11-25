# Theorem 1.2, second independent referee (thm12b)

Evidence level: PROOF-READ (paper's Section 2, pp. 7-14, read line by line against the PDF page
images of pp. 11-13) plus an own proof (PROOF-READ, written here in full, with the external
inputs listed) plus NUMERIC tests of the statement in the forward direction (own code,
`thm12b_numeric.py`, output `thm12b_numeric_output.txt`). Nothing is EXACT except the
elementary identities (2.2) and Gale invariance of q, which are proved on paper below (and
also tested numerically).

## Deviations and limits

- Written independently: the first referee's report and code were not opened.
- No elliptic-function or Pic-theoretic computation was done on actual unitaries (mapping a
  given U to its lattice and Abel-Jacobi points was not implemented). The elliptic formulas of
  my proof are checked numerically only on synthetic data (random lattice, random points).
- External inputs that I did not re-prove: the classical facts that a smooth complete
  intersection of two quadrics in P^3 is an elliptic normal quartic, projectively normal, and
  that the base locus of a pencil is smooth iff det(sQ1+tQ2) has four distinct roots; the
  theory of theta functions (quasi-periodicity, argument principle). Standard, not audited.
- Compute: every run through the job limiter; total about 8 minutes of CPU. The multistart
  seed of a matrix is 5000 + its index in the list selected by the run, so the `trace` and
  `mp1` runs restart from other seeds than the main `ls` run and find other members of the
  (non-unique) solution family; all such runs are in the output file.
- The forward test cannot tell which of the many decompositions is the paper's (the one coming
  from the elliptic curve of U); it tests the statement, not the construction. The
  construction itself is tested only on synthetic elliptic data (PART theta).
- Uniqueness and extension claims are numerical with few starts (6-12 and 8): "not found",
  not "impossible".

## 1. Statement and the reduction to q-coordinates (EXACT, paper eq. (2.2))

For V = [I_4, Z] in C^{4x8} and a 4-subset A of [8] put Delta_A = det V[:,A],
q_A = Delta_A Delta_{A^c}. For I, J in [4] with |I|=|J|=k, A = I^c u (4+J) is a bijection between
the 70 pairs (I,J) and the 70 four-subsets of [8]. Expanding along the identity columns,
Delta_A = +-det Z[I,J], Delta_{A^c} = +-det Z[I^c,J^c]; for Z = U unitary Jacobi's identity gives
det U[I^c,J^c] = (-1)^{sum I + sum J} det U conj(det U[I,J]), hence with q_0 = det U

    q_A / q_0 = eps(I,J) |det U[I,J]|^2,   eps(I,J) = (-1)^{k + sum I + sum J}.         (2.2)

The sign depends only on (I,J), not on U. The cases k = 0 and k = 4 (A = {1234}, {5678}) give
q_A/q_0 = 1 identically. I checked the sign rule on 20 Haar unitaries (max error 2.0e-15).

Consequence used below: if a point P of C^70 satisfies P_A/P_0 = t (q_A/q_0)(V1) + (1-t)(q_A/q_0)(V2)
for all A, with V1, V2, U unitary and P = q(U), then (1.3) holds for every (I,J), including
the empty and the full minor, with lambda = t. So point (v) of the brief is not an issue: one
affine relation among the q-vectors is the same statement for all sizes at once.

Gale invariance (replaces the citation of HMSV Prop. 8.2 on p. 8). If V (4x8, rank 4) and
W (4x8) satisfy V W^t = 0 with rank W = 4, the complementary-minor identity
Delta_A(W) = c * sgn(A, A^c) Delta_{A^c}(V) (c independent of A; sgn the sign of the shuffle)
gives q_A(W) = c^2 sgn(A,A^c) sgn(A^c,A) q_A(V) = c^2 q_A(V), because sgn(A,A^c)sgn(A^c,A) =
(-1)^{4*4} = 1. So q is Gale invariant up to one common scalar. EXACT (classical, proved by the
Laplace expansion of det [V; W]); numerically 5.0e-16.

For U unitary: [I, U] diag(I_4, -I_4) [I, U]^* = I - U U^* = 0, so W = conj([I,U]) diag(1,1,1,1,-1,-1,-1,-1)
satisfies V W^t = 0: the Gale dual of [I,U] is its complex conjugate (paper (2.3)), and the
diagonal matrix of the "diagonal relation" at U is diag(+,+,+,+,-,-,-,-).

## 2. Own proof of Theorem 1.2

### 2.0 General position

Call U in U(4) *generic* if the eight columns of V_0 = [I,U] impose 8 independent conditions on
quadrics (so they lie on a pencil), and the base locus E of that pencil is a smooth curve
(det(sQ_1+tQ_2) has four distinct roots on P^1), and the eight points p_1..p_8 are distinct.
Then E is a smooth elliptic normal quartic, the p_i lie on E, and L_0 := O_E(1) is the
embedding bundle. Restriction H^0(O_P3(2)) -> H^0(E, L_0^2) is onto (10 -> 8, kernel the
pencil), so the quadrics through the eight points are 2 + h^0(L_0^2(-Delta)), Delta = p_1+..+p_8;
"only a pencil" is therefore equivalent to L_0^{(x)2} not= O(Delta).

Real orthogonal U are never generic in this sense: for them L_0 turns out to be a ramification
point (L_0^2 = O(Delta)), and indeed the numerical rank test finds a 3-dimensional space of
quadrics for every real orthogonal matrix tried (sigma_8/sigma_1 ~ 1e-16 in the PENCIL table).
This is consistent with the paper (which excludes them by its general-position hypothesis
and recovers them by the limit argument), not a defect.

### 2.1 Uniformization and the explicit Pluecker coordinates

Write E = C/Lambda, Lambda = Z + tau Z, and let theta(z) = theta_1(pi z | tau) (odd, simple zeros
exactly at Lambda, theta(z+1) = -theta(z), theta(z+tau) = -e^{-i pi tau - 2 pi i z} theta(z)).
Let u_1..u_8 in C be lifts of p_1..p_8. For s in C let Th_4(s) be the space of entire f with
f(z+1) = f(z), f(z+tau) = e^{-4 pi i tau - 8 pi i z + 2 pi i s} f(z). Facts (standard): dim Th_4(s) = 4;
every nonzero f in Th_4(s) has 4 zeros mod Lambda, summing to s mod Lambda; products
prod_{m=1}^4 theta(z - a_m) with sum a_m = s lie in Th_4(s) and span it. Th_4(s) is H^0(E, L_s)
for the degree-4 bundle L_s = O(3*0 + s) (Abel), and s -> L_s identifies E = C/Lambda with Pic^4(E).

For a basis f_1..f_4 of Th_4(s) put x(s) = [f_k(u_i)]_{k,i} in C^{4x8} (a matrix representative
of the configuration x(L_s); other lifts u_i + lambda multiply column i by a nonzero factor).

**Lemma A (Frobenius-Stickelberger form).** There is C(s) not= 0 with, for all z_1..z_4,

    det[f_k(z_j)]_{k,j=1..4} = C(s) * theta(z_1+z_2+z_3+z_4 - s) * prod_{j<l} theta(z_j - z_l).

Proof. Fix z_2,z_3,z_4 generic and view both sides as functions of z_1. The left side is a
linear combination of the f_k, so lies in Th_4(s); it vanishes at z_2,z_3,z_4 (repeated
columns), so (if nonzero) its fourth zero is s - z_2 - z_3 - z_4. The right side has the same
multipliers: under z_1 -> z_1+1 four theta factors each change sign; under z_1 -> z_1+tau the
four factors give e^{-4 pi i tau - 8 pi i z_1 + 2 pi i (z_2+z_3+z_4) - 2 pi i (z_2+z_3+z_4 - s)} =
e^{-4 pi i tau - 8 pi i z_1 + 2 pi i s}. Its zeros (z_1 = z_j, z_1 = s - z_2-z_3-z_4) are simple and are
zeros of the left side. So LHS/RHS is an entire elliptic function of z_1, a constant in z_1.
Both sides are alternating in (z_1..z_4) (theta is odd), so the quotient is a symmetric
function constant in each variable separately on a dense open set, hence a constant C(s).
C(s) not= 0 because f_1..f_4 are independent, so the 4x4 evaluation determinant is not
identically zero. []

Hence for A = {a_1<..<a_4}, writing u_A = sum_{i in A} u_i and K_A = prod_{i<j in A} theta(u_i-u_j)
(natural order),

    Delta_A(x(s)) = C(s) K_A theta(u_A - s),
    q_A(x(s))     = C(s)^2 K_A K_{A^c} theta(s - u_A) theta(s - u_{A^c}).                     (*)

(theta odd: two sign flips cancel.) Numerically (PART theta): the ratio
Delta_A / (K_A theta(u_A - s)) is the same for all 70 A to 5.9e-28 relative (30 digits, 3 values of s).

**(i) of the brief: all q_A are sections of one bundle; the ratios depend on L only.**
By (*), the common factor C(s)^2 drops out of every ratio, and
g_A(s) := theta(s - u_A) theta(s - u_{A^c}) satisfies g_A(s+1) = g_A(s),
g_A(s+tau) = e^{-2 pi i tau - 2 pi i (2s - S)} g_A(s), S := u_1+..+u_8 = u_A + u_{A^c}. The multiplier does
not depend on A. The space Th_2^S of entire functions with these multipliers has dimension 2
(standard: a level-2 theta space). So there are a, b in C^70 and a basis h_1, h_2 of Th_2^S with
q(x(s)) = C(s)^2 (h_1(s) a + h_2(s) b): **the projective point f(s) = [q(x(s))] moves on the
fixed line l = P(span(a,b)), and s -> [h_1(s):h_2(s)] is a degree-2 map E -> P^1.** (Degree 2: a
nonzero element of Th_2^S has two zeros mod Lambda.) Numerically: six random s give a 6x70
matrix of normalised q-vectors whose singular values are 2.37, 0.62, then <3e-29.
This is exactly the paper's (2.6)-(2.7), with W = the bundle whose sections are Th_2^S. Also
f is independent of the choice of basis f_k (change of basis multiplies all q_A by det^2) and
of lifts (column factors, each column occurs in every q_A exactly once). Answer to (i):
**agrees; correct.** The paper's route (Poincare bundle, GRR, deg V = -1) is also correct
(I re-did ch(V) = pr_2*(exp(4w + xi)) = 4 - beta: the cubic terms vanish for degree reasons); my
route avoids the Poincare-bundle normalisation and gives the zero divisor div q_A =
{u_A} + {u_{A^c}} directly.

**The involution iota.** g_A(S - s) = theta(u_{A^c} - s) theta(u_A - s) = g_A(s) exactly (two sign
changes). So f(S - s) = f(s): iota(s) = S - s is iota(L) = O(Delta) (x) L^{-1} in the paper. A degree-2
map has exactly the fibres {s, S - s}, its deck involution is iota and its ramification points are
the four solutions of 2s = S mod Lambda (paper (2.8)). Numerically
|q(S-s) - c q(s)|/|q| = 1.2e-29.

**Goppa/Gale duality, own proof of (2.4).** For f in Th_4(s), h in Th_4(S - s) the product f h has
the multipliers of prod_i theta(z - u_i) =: Theta(z), so f h / Theta is an elliptic function with
at most simple poles at the u_i. Residue theorem: sum_i f(u_i) h(u_i) / Theta'(u_i) = 0. Hence
x(s) diag(w) x(S-s)^t = 0 with w_i = 1/Theta'(u_i) not= 0, and x(S-s) has rank 4 (an element of
Th_4(S-s) vanishing at 8 points is 0). So x(S-s) diag(w) is a Gale dual of x(s):
g(x(L)) = x(iota L) up to column scaling. This is the residue argument the paper sketches on
p. 10 (with the paper citing Eisenbud-Popescu); **agrees, correct**.

### 2.2 The real structure (point (ii))

Generic U. By section 1 the Gale dual of x(L_0) = [I,U] is conj([I,U]). By 2.1 it is also x(iota L_0).
The quadrics through conj(p_i) are the conjugates of those through p_i, a pencil with smooth
base curve conj(E); so x(iota L_0) spans a pencil and its base curve is the image of E under
|iota L_0| (equivalently (iota L_0)^2 not= O(Delta), which is the same condition as for L_0).
The projective map T with T(conj p_i) = x(iota L_0)_i carries one pencil onto the other, so
conj(E) onto E_{iota L_0}. Composing E --L_0--> E_{L_0} --conj--> conj(E) --T--> E_{iota L_0} --(|iota L_0|)^{-1}--> E
gives an antiholomorphic bijection sigma_E of E with sigma_E(p_i) = p_i. sigma_E^2 is a
holomorphic automorphism fixing 8 points; automorphisms z -> zeta z + c of C/Lambda other than the
identity fix 0 or |1-zeta|^2 <= 4 points; so sigma_E^2 = id. **Agrees with the paper (p. 11-12);
correct.** The paper's justification that the pencil of x(iota L_0) still has base curve E is
correct; the shorter reason is that x(iota L_0) is projectively equivalent to conj(x(L_0)).

Coordinates. Choose the origin of C/Lambda at p_1 (a fixed point). An antiholomorphic
involution of C/Lambda fixing 0 is z -> c conj(z) with |c| = 1 and c conj(Lambda) = Lambda; rotating
the coordinate by c^{1/2} makes it z -> conj(z), with Lambda = conj(Lambda) (rectangular or rhombic
lattice). Each u_i can be taken with conj(u_i) = u_i mod Lambda. On Pic^4, sigma_J(L_s) = L_{conj s}, and
x(conj s) = conj(x(s)) up to column factors (use the basis conj(f_k(conj z)) of Th_4(conj s), which is
the image of the basis of Th_4(s); conj(u_i) = u_i + lambda_i only rescales column i).
The relation sigma_J(L_0) = iota(L_0) says conj(s_0) = S - s_0.

The circle. tau_J(s) := iota(sigma_J(s)) = S - conj(s) is an antiholomorphic involution of E = Pic^4
fixing s_0. Its fixed set {s : s + conj(s) = S mod Lambda} = s_0 + {w : 2 Re w in Lambda cap R} is a union of
one or two parallel closed geodesics (two for rectangular, one for rhombic Lambda); B is the one
through s_0. **Agrees (p. 12), correct.**

### 2.3 Unitarity along B and the sign pattern (point (iii))

For s in B, conj(x(s)) = x(conj s) = x(S - s) = Gale dual of x(s), all up to projective
equivalence and column factors. So there are M in GL_4 and an invertible diagonal D(s) with
x(s) D(s) x(s)^* = 0 (rows of M conj(x(s)) D span ker x(s)). Uniqueness of D up to scale: if D_0, D_1
are two invertible solutions, D_1 x^* = D_0 x^* N for some N in GL_4, so N^* fixes every marked line,
so N is scalar (Lemma 2.1 of the paper: the stabiliser of the 8 points in PGL_4 is trivial;
its proof is correct: three points never collinear because h^0(L(-p-q-r)) = 1, at most four on a
plane, so strict Hilbert-Mumford inequalities, and a diagonalisable T with all eight points as
eigenvectors and a proper eigenspace would violate them). For a non-invertible solution D,
D + t D_0 is invertible for all but finitely many t. So the solution space is a complex line,
spanned by an invertible matrix; every nonzero solution has all entries nonzero; conj(D) is a
solution, so D can be normalised real with d_1 = 1.

Constancy of the sign pattern. Locally on B choose a real-analytic frame, i.e. a matrix X(s)
depending continuously on s, representing x(s) with nonzero columns. The real linear map
D -> X(s) D X(s)^* on real diagonal matrices has a one-dimensional kernel at every s, so its
rank is constant and the normalised kernel vector d(s) (d_1 = 1) is continuous; its entries never
vanish, so their signs are locally constant. Changing representatives multiplies d_i by
|c_i|^{-2} > 0 and left multiplication does nothing, so the sign vector is a well-defined
function on B, locally constant, hence constant on the connected B. At s_0 it is
(+,+,+,+,-,-,-,-) by section 1. **The paper's statement (p. 12, "the resulting solution varies
locally continuously") is correct; it leaves out the constant-rank argument, which I supply. Not a
real gap.**

From the signs: sum_{i<=4} d_i v_i v_i^* = sum_{i>4} |d_i| v_i v_i^*; both sides are PSD with ranges
span(v_1..v_4) and span(v_5..v_8); they are equal, so both spans are the span of all eight, C^4.
Hence x(s) = [I, Z] after a left change of basis, D_I = Z D_J Z^* with D_I, D_J positive diagonal, and
U(s) := D_I^{-1/2} Z D_J^{1/2} is unitary and represents x(s) after column scalings. In particular
q_0 = Delta_{1234} Delta_{5678} not= 0 on B. **Agrees ((2.10)-(2.11)), correct.** Numerical
illustration (PART theta): on both real circles of a rectangular lattice with eight real points
the sign vector of D is constant over 24 sample points per circle (pattern +-+-+-+- in the
cyclic order of the real points, i.e. a 4/4 split, on both circles), the kernel is exactly
one-dimensional throughout (sigma_8/sigma_1 <= 1.3e-16, sigma_7/sigma_1 >= 1.4e-4), and the
configuration is real (cross-ratios of Z real to 5e-14) exactly at the two ramification points
y = 0 and y = tau/2i of each circle, and visibly not real between them (Im 0.16).

### 2.4 Endpoints (point (iv))

f(tau_J s) = conj f(s) (from x(conj s) = conj x(s)) and f(iota s) = f(s), so on B the affine
coordinates q_A/q_0 are real: f(B) lies in the real affine line l(R) cap {q_0 not= 0}, which is
an honest real affine line because l contains more than one real point. f|B is real analytic
and nonconstant (otherwise f would be constant on E by the identity theorem), so f(B) = [P_-, P_+]
is a nondegenerate compact interval. If s in B has f'(s) not= 0 then f is a local biholomorphism
there and maps B onto a real arc through f(s) contained in l(R), so f(s) is interior. Hence
endpoints are images of ramification points: iota(s) = s, i.e. 2s = S. On B this also gives
conj(s) = s (since tau_J s = s). Conversely B meets Fix(iota) in exactly two points in the
rectangular case (each vertical geodesic meets each of the two real horizontal ones once) and,
since f(B) has two distinct endpoints, in at least two points in any case; so iota|B is a
reflection of the circle with two fixed points, f|B is injective on each of the two arcs, and
f(B) is the image of one arc. **Agrees with p. 13 (my version adds the injectivity, not needed).**

At a ramification point s on B: x(s) is projectively equivalent to conj(x(s)). Let T be the
antilinear map v -> M conj(v) realising this; T preserves each marked line, T^2 is linear and fixes
the eight lines, so T^2 = nu I (Lemma 2.1); T v_i = mu_i v_i gives nu = |mu_i|^2 > 0; T/sqrt(nu) is an
antilinear involution, C^4 = R (+) iR with R its fixed real 4-space, and each marked line contains a
nonzero vector of R. So x(s) has a real representative X. For it D is real with the same signs
(+,+,+,+,-,-,-,-); after a real left change of basis X = [I, Z] with Z real, and
D_I^{-1/2} Z D_J^{1/2} is real and unitary, i.e. in O(4). **Agrees, correct.** The positivity needed
for an orthogonal (not O(p,q)) matrix is exactly the constancy of the 4/4 sign split along B,
which holds.

Conclusion for generic U: with O_-, O_+ the orthogonal matrices at the two ends,
q(U)/q_0 = lambda q(O_-)/q_0 + (1-lambda) q(O_+)/q_0 in R^70 with lambda in [0,1] (f(s_0) in [P_-,P_+]),
and by (2.2) this is (1.3) for all 70 pairs (I,J) at once.

### 2.5 Non-generic U (point (vi))

Genericity is a Zariski-open condition on Z in GL_4(C): 8 independent conditions (a minor of the
8x10 evaluation matrix is nonzero), then the discriminant of the binary quartic
det(sQ_1+tQ_2) (Q_i rational in Z on that open set) is nonzero, and the points are distinct.
Nonempty: the paper's example on p. 14 (eight general points on a quartic with
O(Delta) not= L^2 and both quadruples independent) is correct, and numerically all 30 Haar samples
are generic. U(4) is Zariski dense in GL_4(C) (a polynomial vanishing on U(4) vanishes on a
totally real submanifold of maximal dimension 16, hence on GL_4). So the defining polynomial is
not identically zero on the connected real-analytic manifold U(4), its zero set there has empty
interior, and generic unitaries are dense. For any U take generic U_n -> U, the triples
(O_{-,n}, O_{+,n}, lambda_n) have a convergent subsequence in the compact O(4)^2 x [0,1], and r is
continuous. **Agrees (pp. 13-14), correct.** Data need not converge along the whole sequence
(and need not be unique for degenerate U), but subsequence convergence suffices.

## 3. Step-by-step comparison with the paper

| Paper step (page, eq.) | Verdict |
|---|---|
| (2.1) minor expansion of det(D_a + U D_b U^*) | agrees (Cauchy-Binet; not needed for Thm 1.2) |
| (2.2) q_A/q_0 = sign * |minor|^2 | agrees; EXACT on paper, sign rule (-1)^{k+sum I+sum J} numerically 2e-15 |
| p. 8 Gale dual g([I,Z_0]) = [I, Z_0^{-t}], (2.3) | agrees |
| p. 8 "q_A invariant under Gale duality" via HMSV Prop. 8.2 | correct; citation replaceable by the complementary-minor identity (section 1) |
| p. 9 general-position hypotheses | agrees; note real orthogonal U are never in general position (extra quadric), harmless |
| Lemma 2.1 (stability, trivial stabiliser) | agrees, correct |
| (2.4) Goppa duality g(x(L)) = x(iota L) | agrees (residue argument; own proof in 2.1) |
| (2.5)-(2.6) deg V = -1, deg W = 2, all q_A sections of W | differs but correct (I use theta functions, Lemma A) |
| p. 11 spanning, base-point freeness, (2.7) degree 2, (2.8) ramification | agrees |
| p. 11-12 sigma_E exists, is an involution | agrees, correct |
| p. 12 sigma_J, tau, circle B | agrees |
| (2.9) diagonal relation, 1-dim solution space | agrees, correct |
| p. 12 continuity of normalised D along B | correct; the constant-rank kernel argument is implicit (expository, not a gap) |
| (2.10) constant signs, (2.11) unitary representative, q_0 not= 0 on B | agrees, correct |
| p. 13 f(B) closed interval, endpoints are ramification points | agrees, correct |
| p. 13 real representative at endpoints, orthogonal | agrees, correct |
| p. 13 "lambda in [0,1] ... by (2.2) exactly (1.3)" | agrees; covers all sizes incl. empty and full minor |
| 2.5 density and compactness | agrees, correct |

No step is wrong. No step needs a repair beyond one-line expansions.

## 4. Numerical evidence (forward direction, NUMERIC)

All from `thm12b_numeric.py`; output in `thm12b_numeric_output.txt`. Test set: 30 Haar
unitaries (seeds 1000-1029, QR of complex Gaussian with phase fix), the 4x4 Fourier matrix, a
random real orthogonal matrix (seed 2000), the matrix O* of Example 1.3, block-diagonal 2+2
(seed 2001) and 3+1 (seed 2002), a unitary with two zero entries (seed 2003), one with exactly
one zero entry (seed 2004), P exp(i eps H) near a permutation P (eps = 0.3, 0.1, 0.03, 0.01,
seed 2005) and O exp(i eps H) near a random orthogonal O (same eps, seed 2006).

(a) Pencil of quadrics. All 30 Haar samples, `twozeros`, and all near-permutation and
near-orthogonal samples are in general position: 8 independent conditions (sigma_8/sigma_1 of
the 8x10 evaluation matrix >= 2.1e-3), smooth base quartic (min chordal root separation of
det(sQ_1+tQ_2) >= 2.1e-2 for Haar), distinct points. Not in general position: every real
orthogonal matrix (orth, O*, block 2+2 which is phase-equivalent to a real one): a third
quadric (sigma_8/sigma_1 ~ 1e-16), as the theory predicts (ramification point); the Fourier
matrix and `onezero` (pencil with a singular base curve, root separation ~1e-8); block 3+1
(coincident points e_4 = u_4). These are covered only by the limit argument of Sec. 2.5.

(b) Decomposition by multistart Levenberg-Marquardt on SO(4) x SO(4) x [0,1] (Lie-algebra
retraction, lam = sin^2 mu; 6 starts per Haar matrix, 12 per structured one). **Every one of
the 47 matrices admits a decomposition at double precision**: best max-residual over the 70
coordinates <= 5.7e-14 for every Haar sample except haar[1022] (3.6e-13 in the main run,
9.5e-15 on a rerun), <= 1.5e-14 for every structured one; lam in [0,1] by construction; most
starts converge (typically 5/6 or 6/6; worst 1/6 for haar[1022] and 1/12 for the eps = 0.01
near-permutation, where the problem is badly conditioned).

Tracing (9 matrices: haar 1000-1002, Fourier, block 3+1, twozeros, onezero, near-perm 0.1,
near-orth 0.1): continuation from V = U along the segment, solving r(V) = P(t) for V in U(4)
at 11 points towards each end. Every point of the segment is r of a unitary (worst residual along a path 1.3e-6, attained
at the last step, which is the fold point where convergence is only of square-root type, and
far smaller in the interior); the realness defect max |Im V_ij V_kl conj(V_il V_kj)| decreases
from its value at U (2.6e-2 to 7.2e-2) to 1e-8..2e-4 at the endpoints: the path
of unitaries ends at matrices phase-equivalent to real orthogonal ones, as in Sec. 2.4.
Extension test: at t = 1 + 0.02, 1 + 0.1 (and symmetrically below 0) the best residual over 8
starts (the continuation point and 7 Haar starts) is 6e-5..2e-2 resp. 3e-4..9e-2, of the order of
the offset (0.05-1 times it): the segment is not extendable within r(U(4)) (NUMERIC, 8 starts
only, so "not found", not "proved impossible").

(c) High precision. The three worst Haar cases of the main run (haar[1022], [1004], [1016])
plus haar[1000] were re-solved at 40 working digits after re-unitarising U and the endpoint
matrices by the Newton polar iteration: damped Gauss-Newton reaches max-residual 7.2e-40,
6.5e-38, 6.9e-41 (and 2.3e-41) in two steps, with |O^tO - I| ~ 1e-41 and lam in [0,1]. Nothing
stays above rounding level: **no counterexample candidate.**

(d) Uniqueness: **the decomposition is NOT unique, even up to the obvious symmetries**
(O -> S_1 O S_2 with sign diagonals, swapping the endpoints). Different starts give different
pairs {r(O-), r(O+)} with different lam (e.g. haar[1000]: lam = 0.742, 0.787, 0.785, 0.793,
0.842, 0.579, with different segment directions). The Jacobian of
(O-, O+, lam) -> lam r(O-) + (1-lam) r(O+) at a solution has rank 12 of 13 for all generic U
(13th singular value ~1e-10 in double, ~5e-19 relative at 40 digits with step 1e-23), so the
solutions form (at least) a one-parameter family; the observed classes are samples of it.
Every segment found is tangent at r(U) to the 9-dimensional image r(U(4)) (residual of the
direction after projection onto the tangent space ~1e-10 in all cases), and the traced ones are
covered by unitaries. The paper's construction produces one particular member (the line of the
elliptic curve of U); nothing in Theorem 1.2 or its use (Corollary 2.2) needs uniqueness, so this
is an observation, not an objection. The affine span of r(U(4)) and of r(O(4)) is 13-dimensional
(numerical rank of 60 samples), consistent with the 14 independent q_A of HMSV.

Limits. Near a permutation matrix (eps -> 0) the decompositions found have lam -> 1 (lam up
to 0.99998 at eps = 0.01, all classes >= 0.94 at eps = 0.03), i.e. r(U) approaches the endpoint
r(O-), which approaches r(P); the opposite endpoint does not converge to anything specific
(many classes, segment lengths 0.1-1). Near a real orthogonal O the behaviour is different:
lam stays spread over (0.5, 1) but the segments shrink (|r(O-) - r(O+)| = 0.23-0.95 at
eps = 0.3, 0.04-0.14 at eps = 0.01), i.e. both endpoints approach r(O).
For U itself real orthogonal the problem degenerates: lam = 1 with O- = U (up to signs) and
arbitrary O+ (Jacobian corank 7), plus genuine short segments. Fourier: lam = 1/2 exactly in all
eleven classes found (symmetry). No continuity of (O-, O+, lam) in U is claimed or observed;
the paper's limit argument only uses subsequences, which is all that is needed.

## 5. Verdict

**Theorem 1.2: proved as written** (modulo the standard external facts listed under
"Deviations and limits", and two citations, HMSV Prop. 8.2 for Gale invariance of the q_A and
Eisenbud-Popescu for Goppa duality, both of which I replace by elementary proofs above).
My own proof (Sec. 2) follows the same architecture with explicit theta functions, is complete
at the PROOF-READ level, and agrees with the paper at every step on points (i)-(vi) of the
brief. Expository points, none a gap: (1) p. 12, continuity of the normalised diagonal solution
D_0(L) should cite the constant-rank kernel argument; (2) p. 8-9, real orthogonal U are never in
general position (they have a third quadric), so they really need Sec. 2.5; (3) the
representation is not unique (one-parameter family numerically), which the paper never claims.
Numerically: 47 matrices, every one decomposes (max residual <= 3.6e-13 in double, <= 7e-38 at
40 digits for the four refined cases); no counterexample candidate.

