Evidence level: EXACT for the ranks, Catalan numbers, P_5, kappa, dim X_5, minor facets of Q_5 and the Burnside orbit counts (integer/rational arithmetic, python-flint); EXACT modulo p = 2^31 - 1 at random points (randomized, Schwartz-Zippel type) for the quadric ideal of X_5 and the Lemma 3.4 analogue K = E; NUMERIC (multistart local optimization, LP column generation) for every unitary-versus-orthogonal statement; PROOF-READ for the dimension counts and the architectural map; the rest is labelled JUDGMENT.

# Can the n = 4 route reach n = 5? A map of what extends

Auditor label `n5`. Paper: D. Zheng, "The Marcus-de Oliveira conjecture for n = 4", 24 Oct 2026 (Sections 1.2, 2, 3, 4, 7, Appendices A, B, C). Script: `generalization_n5.py` (parts `alg`, `burnside`, `facets`, `uo1`-`uo4`); output: `generalization_n5_output.txt`.

## Deviations and limits

- The facet enumeration of the 41-dimensional polytope Q_5 (120 vertices) was **not done**. I checked exactly that all 125 minor vectors are facets and found one non-minor facet (Section 3). A random-ray facet search was tried as a cheap probe; it is biased towards large facets and found only minor facets even for n = 4 (where 48 of the 82 facets are non-minor), so it says nothing about how many non-minor facets Q_5 has.
- The quadric ideal of X_5 was computed modulo the prime p = 2^31 - 1 from 963 random points, and each of the 300 quadrics found was re-checked at 40 fresh random points. Dimensions over Q are inferred from the computation mod p (they agree for all but finitely many primes). No computation over Q was done for this part.
- The unitary-versus-orthogonal test is numerical only. A minimum over a group is found by local search from random starts, so each reported minimum is an upper bound. The comparison can expose a gap. It cannot prove that two minima are equal.
- I did not read the 90 case proofs of the supplement. The section on case-proof tools relies on the paper's own description of the proof groups (Table 1, p. 26, and the text after it) and on the statements of Appendices B and C. I did not audit those appendices.
- The cited classical facts (First Fundamental Theorem / Kempe generation, GIT surjectivity, the dimension of degree-one invariants as a count of standard tableaux, Birkhoff-von Neumann, Riemann-Roch) are CITED, not re-proved.

## 0. Summary table

| ingredient (paper location) | n = 5 status | level |
|---|---|---|
| (2.1) Cauchy-Binet minor expansion, (2.2) complementary products, Gale self-duality [I,U] <-> [I,U bar] | extends verbatim (all n) | PROOF-READ |
| Theorem 1.2 (two-orthogonal representation) via elliptic normal quartic | **proof does not extend**: for n >= 5 a general [I,U] lies on no elliptic normal curve (codimension n(n-4) = 5) | PROOF-READ (dimension count) |
| Statement r_5(U(5)) ⊆ conv r_5(O(5)) | **open, numerically supported**. Same linear span (dim 42) proved. No separating functional found among 17 n = 5 test functionals; 3/3 Haar U(5) lie inside the hull of sampled orthogonal points | EXACT (span), NUMERIC (rest) |
| Matching cone X_n: irreducible, dim 2n - 2, P injective, kappa, closedness | extends: dim X_5 = 8, rank P_5 = 42 = C_5 | EXACT + CITED |
| Lemma 3.4 (line x + C kappa ⊂ X forces c(x) ∈ E) | extends **with a new proof**. There is no cubic, but the 300 quadrics through X_5 give K = E_5 (dim 17) | EXACT mod p |
| Lemma 3.5, Theorem 3.6 (contact reduction) | extend verbatim with "six" replaced by 2n - 2 = 8. They also work for a fixed unitary U, with no orthogonal reduction | PROOF-READ |
| Prop. 4.1-type facet list of Q_5, Lemma 4.2 | Lemma 4.2 is general. The facet list is unknown: 125 minor facets plus at least 1 non-minor facet found | EXACT (what is listed) |
| Classification: orbits of 8-subsets of S_5 under the group of order 28800 | 29,342,815 orbits (n = 4: 214 orbits of 6-subsets) | EXACT |
| Case proofs, groups IV and V (matching cubic, Hessian) | **do not extend as written**: X_5 is cut out by 300 quadrics, which is more than the 42 partials of one cubic | EXACT mod p |
| Appendix C (one circle) | does not extend as written: it uses Theorem 1.2, the n = 4 Gale quotient M_8 ⊂ P^13 and the cubic | PROOF-READ (statements only) |
| Minor decompositions, Gram arguments, Lemma 3.5 maximum principle, Birkhoff | general in n | PROOF-READ / JUDGMENT |
| Appendix B (coincidences) | proved only for 4x4; the method (common weights) is not obviously n-specific | JUDGMENT |

## 1. Unitary versus orthogonal (Theorem 1.2)

### 1.1 What the n = 4 proof uses (PROOF-READ, pp. 7-14)

The 2n columns of [I_n, U] are viewed as points of P^{n-1}. Each |det U[I,J]|^2 equals ± q_A/q_0, where q_A = Δ_A Δ_{A^c} is a product of complementary maximal minors (2.2). The steps are:

1. The 8 points lie on the base curve E of the pencil of quadrics through them, an elliptic normal quartic.
2. Re-embedding E by every L ∈ Pic^4(E) gives a family of configurations. The q_A are sections of a degree-2 line bundle W on Pic^4(E), so they map Pic^4(E) two-to-one onto a line ℓ.
3. Gale duality is L ↦ O(Δ) ⊗ L^{-1}. Together with complex conjugation it gives a real circle B through the initial L_0, and every point of B is again unitary.
4. f(B) is a closed interval whose endpoints are ramification points, and these points are real, so they are orthogonal.

### 1.2 Dimension count for general n (PROOF-READ)

- Configurations of 2n points in P^{n-1} modulo PGL_n: 2n(n-1) - (n^2 - 1) = (n-1)^2. This matches U(n) modulo left and right diagonal phases, of real dimension n^2 - (2n - 1) = (n-1)^2. Gale self-duality [I,U]^Gale = [I, U^{-T}] = [I, conj U] holds for every n.
- Quadrics in P^{n-1}: n(n+1)/2. Quadrics through 2n general points: n(n+1)/2 - 2n = n(n-3)/2. That is a pencil for n = 4, and 5 quadrics for n = 5.
- Configurations of 2n marked points on an elliptic normal curve of degree n, modulo PGL_n, have dimension 1 (the j-invariant) + 1 (L ∈ Pic^n) + 2n (the points) - 1 (translations) = 2n + 1. For n = 4 this is 9 = (n-1)^2, so the elliptic locus is dense. That is the paper's statement that eight general points lie on an elliptic quartic. For n = 5 it is 11 < 16. In general the elliptic locus has codimension (n-1)^2 - (2n+1) = n(n-4), which is 5 for n = 5. A second count agrees: elliptic normal curves of degree n in P^{n-1} form an n^2-dimensional family, and each point on the curve imposes n - 2 conditions, giving n^2 - 2n(n-2) = -n(n-4).
- U(n) is Zariski dense in GL_n (paper p. 14), so for n >= 5 the configuration [I,U] of a **general** unitary U lies on no elliptic normal curve. The first step of the proof has nothing to act on. Nothing in steps 2-4 replaces it.
- JUDGMENT, not checked: on the codimension-n(n-4) elliptic locus the remaining steps look dimension-free. The Grothendieck-Riemann-Roch computation gives deg V = -1 and deg W = 2 for every n, so the q_A still map Pic^n(E) two-to-one onto a line. For unitary U on that thin locus a two-point orthogonal representation is therefore plausible. It says nothing about general U.

### 1.3 What remains true for every n

- **(EXACT + CITED) The linear spans agree.** For O ∈ SO(n), the generalized Laplace expansion plus Jacobi's identity give Σ_{σ(I)=J} sgn σ Π O_{iσ(i)} = det O[I,J]^2. The proof in `lib/README.md` is written for n = 4 but uses nothing about n. So every functional ℓ_w(r) = Σ w_{IJ} r_{IJ} satisfies ℓ_w(r(O)) = Σ_σ sgn σ v_σ(w) Π O_{iσ(i)}, where v_σ(w) = ℓ_w(r(P_σ)), and span r_n(SO(n)) has dimension rank{χ_{I,J}} = C_n (EXACT: 5, 14, 42 for n = 3, 4, 5). On the unitary side, the q_A are degree-(1,...,1) SL_n-invariants of 2n vectors in C^n, a space of dimension C_n (CITED: number of standard tableaux of shape 2^n). So span r_n(U(n)) ⊆ a C_n-dimensional space that contains span r_n(SO(n)), and the two spans are equal. NUMERIC confirmation: the sample matrices have rank 42 for U(5), for SO(5) and for both together. Functionals in the 84-dimensional space of k ≤ 2 functionals that vanish at all 120 permutation matrices are at most 2.9e-15 at random U(5).
- **Consequence for the architecture (PROOF-READ).** For unitary U, Φ_c(U) := ℓ_w(r(U)) with v(w) = Pc is well defined and linear in c. Φ_{n_{I,J}}(U) = |det U[I,J]|^2 ≥ 0, and by (2.1), Φ_{c(x)}(U) = Re(α det(D_a + U D_b U*) - F_{π*}(x)). Proposition 3.3, Lemma 3.5, Theorem 3.6 and every minor decomposition therefore hold verbatim with SO(n) replaced by U(n). The orthogonal reduction is not needed for the contact reduction. It is needed only where the case proofs use real entries of O.

### 1.4 Numerical test of r_5(U(5)) ⊆ conv r_5(O(5)) (NUMERIC)

Method. Every functional is written in the coordinates |U_ij|^2 and |det U[R,C]|^2 with |R| = |C| = 2. Minors with k ≥ 3 are folded onto complements, using |det U[I,J]| = |det U[I^c,J^c]| on U(n). The folding and the Leibniz identity were checked on SO(5) to 2.6e-15. Minimization uses multistart L-BFGS-B in the chart X_0 · exp(·) around Haar-random X_0: U = X_0 exp(iH) with H Hermitian (25 parameters), O = X_0 exp(S) with S skew (10 parameters). The best point is then polished. Seeds are in the output.

Results (min over SO(n) versus min over U(n)):

| functional | min over perms | min SO(n) | min U(n) |
|---|---|---|---|
| n = 3, two random functionals (control) | -2.8671299884, -0.9437297927 | same | same |
| n = 4, h_{id,-1} and h_{id,+1} (control) | 0 | -0.0674422489 | -0.0674422489 |
| n = 5, h_{id,±1} and h_{2143,-1} on the block rows 1-4 x cols 1-4, plus c(1 - r_55) to make it ≥ 0 at all permutations | 0 | -0.0674422489 | -0.0674422489 |
| n = 5, h_{id,-1} with index 5 appended to every I, J (plus correction) | 0 | -0.0674422489 | -0.0674422489 |
| n = 5, h_{id,-1} on rows 1-4 x cols 2-5 (plus correction) | 0 | -0.0674422489 | -0.0674422489 |
| n = 5, E4 = block embedding + appended embedding (a **non-minor facet** of Q_5, Section 3) | 0 | -0.1348844977 | -0.1348844977 |
| n = 5, four positive combinations of three random G-images of E4 | 0 | 0, 0, 0, -1.990e-6 | 0, 0, 0, -1.990e-6 |
| n = 5, six realizable functionals Re(θ det(D_a + U D_b U*)), random Gaussian a, b, θ | σ-point minimum | = perms | = perms |

A seventh embedding (E5: block plus shifted block) could not be lifted to be nonnegative at all permutations: its minimum at a permutation is -0.5, and SO(5) and U(5) give the same -0.5. In every row the two minima agree to within 1e-7. Most agree to 1e-12. One run first showed min U(5) = -1.99e-6 below min SO(5) = 0 (combination #3, 40 starts). A follow-up (part `uo4`) with 600 SO(5) starts found -1.990e-6 on SO(5) as well, with a difference of 1.6e-13. That was a search artifact, not a gap. The E-type functionals have block-diagonal minimizers (value exactly -0.0674422489, the n = 4 constant, and twice that for the sum), so in effect they re-test n = 4. The realizable ones are minimized at σ-points. Only the facet-image combinations probe genuinely five-dimensional directions.

Membership test (column generation, part `uo3`). For a Haar U, the test alternates two steps. First it solves the LP "separate r(U) from the convex hull of the sampled r(O)" (3000 Haar O plus the 120 permutation matrices, in the 41 affine coordinates). Then it minimizes the separating functional over SO(n) and adds the minimizers to the sample. Results:

- n = 3 (1 U) and n = 4 (2 U): r(U) was in the initial sample hull. For n = 4 this is consistent with Theorem 1.2.
- n = 5, 3 Haar U(5): r(U) **inside** the sample hull. For one U this held at the start. For the other two it held after 2 and 7 rounds. In those rounds the minimized separating functionals undercut r(U) by 0.06-0.22, and the LP margin then dropped to ≤ 1e-9.
- So these three random unitary points are numerically convex combinations of orthogonal points (a floating-point LP certificate). Haar points tend to lie deep inside, so this does not probe the boundary.

**Verdict on item 1.** No functional with min over U(5) clearly below min over SO(5) was found, and three Haar U(5) points were found inside the convex hull of sampled r(O). This is evidence only: it does not show that the inclusion holds. The proof of Theorem 1.2 certainly does not extend, because a general unitary configuration lies on no elliptic normal quintic. Whether the statement holds for n = 5 is open. It is the first of the decisive questions in Section 5.

## 2. The matching cone X_n

All EXACT (part `alg`):

| n | C_n | rank of noncrossing products | rank of the n! bipartite p_σ | rank of all (2n-1)!! matchings | rank of χ_{I,J} on S_n | dim E_n = rank of χ_{i,j} | dim X_n |
|---|---|---|---|---|---|---|---|
| 3 | 5 | 5 | 5 | 5 (15) | 5 | 5 | |
| 4 | 14 | 14 | 14 (24) | 14 (105) | 14 | 10 | 6 |
| 5 | 42 | 42 | **42 (120)** | 42 (945) | 42 | 17 | **8** |
| 6 | 132 | 132 | 132 (720) | not computed | | | |

- The n! σ-point functionals span the whole degree-one space for n = 3, 4, 5, 6, so P_n is injective. P_5 is an integer 120 x 42 matrix with entries {0, 1, 2, 3}, and im P_5 = span{χ_{I,J}}. In the noncrossing basis, kappa = e_41 for n = 5 (e_13 for n = 4), and P kappa = 1.
- dim X_5 = 8. The Jacobian of the 42 noncrossing coordinates has rank 8 at an integer point. The 3 sl_2 fields and the 9 torus fields with Π t = 1 are independent there and lie in its kernel. In general dim X_n = 4n - (2n + 2) = 2n - 2.
- Irreducibility (closure of the image of an irreducible variety), closedness of the image and the finite spectral representation (3.1) (GIT surjectivity, CITED) do not depend on n.

**Lemma 3.4 for n = 5.** The paper's proof uses the cubic C, whose 14 partials generate the ideal of X_4, and its Hessian at kappa. For n = 5 the ideal is different:

- dim I_2(X_5) = 300 quadrics and dim R_2 = 603 (n = 4 control: 14 and 91, matching HMSV). 0 failures at 40 fresh points. EXACT mod p.
- 300 > 42, so these quadrics are not the partials of one cubic in 42 variables. There is no analogue of C.

The polarization argument does not need a cubic:

- For any quadric f vanishing on X, f(x + tκ) = f(x) + t ∇f(κ)·x + t^2 f(κ), and κ ∈ X.
- So if the line x + Cκ lies in X, then x ∈ K := {y : ∇f(κ)·y = 0 for all f ∈ I_2}.
- Computed: rank{∇f(κ)} = 25, dim K = 17, and **K = E_5** (rank[K; E_5] = 17). n = 4 control: rank 4 (= rank H_κ), K = E_4, dim 10.
- K is defined over Q and κ ∈ K, so c(x) = Re x - Re F_{π*} κ ∈ E_5 ∩ V_R. Then Φ_{c}(O) = Σ q_ij O_ij^2 ≥ 0 by the k = 1 cofactor identity and Birkhoff, for any n. For unitary U it equals Σ q_ij |U_ij|^2, again Birkhoff.

So **the Lemma 3.4 analogue holds for n = 5** (EXACT mod p, randomized). The argument uses only quadrics that vanish on X, not HMSV's generation theorem.

Tangent test (mod p). At 3 random points of X_5, dim T_x X = 8 and κ ∉ T_x X, so a generic point has no κ-line. At a point whose first spectrum has a_2 = ... = a_5, κ ∈ T_x X and x ∈ E_5. There F_σ depends only on σ(1), and adding a constant to it can be realized again by a spectrum of the same shape: values A + B/b'_j, solved explicitly (PROOF-READ). So {x : x + Cκ ⊂ X_5} = X_5 ∩ E_5. It contains the locus where one spectrum has n - 1 equal entries. I did not determine whether it contains more.

**Theorem 3.6 for n = 5 (PROOF-READ).**

- With k contacts the proof fixes k complex affine conditions (k - 1 differences and the sum). Every component through x has dimension ≥ 2n - 2 - k, which is ≥ 1 for k ≤ 2n - 3.
- The exceptional case, where all g_σ are constant, is a κ-line inside X. It is handled by the Lemma 3.4 analogue above, using injectivity of P (common kernel of the g_σ = Cκ, EXACT rank 42).
- Lemma 3.5 is general.
- The argument goes through verbatim and gives **at least 2n - 2 = 8 contact permutations** for n = 5. The paper's "six = dim X" becomes 8 = dim X_5. With the remark in 1.3 it works for fixed U ∈ U(5) as well as for O ∈ SO(5).

## 3. The classification step

- Q_5 = conv{v_σ} has 120 vertices and linear rank 42, so dim Q_5 = 41 (EXACT). For comparison, Q_4: 24 vertices, dim 13, 82 facets in 3 G-classes (16 + 18 minor, 48 non-minor h_{π,ε}; Prop. 4.1, replayed in `lib/`).
- All 125 distinct minor vectors of Q_5 (25 with |I| = 1, 100 with |I| = 2, complements identified) define facets: each zero set has linear rank 41 (EXACT).
- **A non-minor facet exists.** E4 = h_{id,-1} on the block [4] x [4] plus h_{id,-1} with 5 appended to every I, J. Its values are 2 x {0^65, 1^40, 2^14, 4^1}. Its zero set has 65 vertices spanning a 40-dimensional face (EXACT), and its Leibniz minimum on SO(5) is -0.1349 < 0 (NUMERIC). So Q_5 has non-minor facets, as Q_4 does. The lifted single embeddings (zero sets of 17 or 22 vertices) are faces of dimension 16 or 21, not facets.
- Full facet enumeration: **not done**. qhull in dimension 41 with 120 points was judged infeasible within the budget and was not attempted. The random-ray probe found 53 distinct facets, all minor (|zero set| = 108), and for n = 4 it also found only minor facets. It is not informative about the non-minor facets.
- Symmetry group G_5 = (S_5 x S_5) ⋊ C_2 of order 2 · 120^2 = 28800, acting faithfully on S_5 (EXACT).
- Burnside orbit counts of k-subsets (EXACT; the control k = 6, n = 4 gives **214** of C(24,6) = 134596):

| k | n = 4 orbits | n = 5 orbits |
|---|---|---|
| 6 | 214 | 133,816 |
| 7 | 447 | 2,099,444 |
| 8 | 891 | **29,342,815** |
| 9 | 1451 | 363,781,693 |

The reduction for n = 5 lands on contact sets of size 8, which is 29.3 million classes against 214. Classes whose generators are all minor vectors are disposed of automatically, as the 124 were for n = 4. Their number cannot be computed without the facet list. Even if only a small fraction need a real argument, a hand-written case-by-case proof like the paper's 90 cases (89 in the supplement) is out of reach. Any n = 5 proof along this route needs either an automated certificate per class, or a uniform argument of the kind the paper asks for in Section 7 (the curve method of Section 5 applied to all cases).

## 4. The case-proof tools (JUDGMENT, from Table 1 and the appendix statements)

- **General in n.**
  - Minor decompositions (Lemma 4.2 plus Φ_{n_IJ} = det^2, all n).
  - Gram/PSD arguments, which are linear algebra in the entries of O or U.
  - The maximum principle on curves (Lemma 3.5, Group VI) and the curve construction of Section 5. The construction needs an analogue of the "fix one coefficient" bound for each case.
  - Birkhoff on the entry subspace.
  - Multiplicative relations among σ-points (Group IV's F_{π3}F_{π11} = F_{π5}F_{π9} type identities exist for every n, but which ones are useful is case-specific).
- **Specific to n = 4.**
  - The matching cubic C and its Hessian (Groups IV and V: "C(u) = 0 and det Hess C(u) = 0" for the real part u of x ∈ X). For n = 5 the replacement is the 300-dimensional space of quadrics. The analogues of "C(u) = 0" and the Hessian equation would be the real and imaginary parts of f(u + iv) = 0 for all f ∈ I_2, with no single polynomial to exploit.
  - Appendix C (one spectrum on a circle). Its intro says it uses Theorem 1.2, the n = 4 Gale quotient M_8 ⊂ P^13 with Sing C = M_8, the secant variety and the degree-7 Hessian factor (C.2). All of these are n = 4 objects.
  - The zero-determinant Lemma C.3 is the exception. Its argument (at most 5 positive weights, then reduction to size n - 1) carries over to n = 5 if MOC(4) is available, because then each spectrum has support at most 4 (PROOF-READ of the paper's argument, with n - 1 = 4 in place of 3).
  - Appendix B (coincidences) is stated and proved for 4x4. The common-weight method is not obviously n-specific, but I did not check it.
- **Real orthogonality.** The paper notes that the group II-V arguments work with real O (p. 26-27). With no Theorem 1.2 analogue for n = 5 they would have to be redone over U(5) (16 real parameters instead of 10), or the orthogonal reduction would have to be proved first.

## 5. Bottom line

- **Extends:**
  - (2.1), (2.2), Gale self-duality.
  - Equality of the linear spans of r_n(U(n)) and r_n(O(n)).
  - The matching cone (dim 2n - 2, C_n = 42 coordinates, P injective, κ).
  - Lemma 3.5.
  - Theorem 3.6 with 8 contacts, for SO(5) or directly for U(5).
  - Lemma 4.2.
  - Minor decompositions and Birkhoff.
- **Extends with a new idea:**
  - Lemma 3.4. There is no cubic. Quadric polarization gives K = E_5, checked mod p.
  - The facet list of Q_5. A new enumeration is needed: 125 minor facets plus at least one non-minor facet are known.
  - The case classification. 29,342,815 orbits of 8-sets need automation or a uniform argument.
- **Does not extend as written:**
  - The proof of Theorem 1.2: no elliptic normal curve exists for general n ≥ 5.
  - Groups IV and V: the cubic and its Hessian.
  - Appendix C, apart from Lemma C.3.
  - The 90 hand-proved cases, by scale.

Questions that would decide whether an n = 5 proof by this route is conceivable:

1. Is r_5(U(5)) ⊆ conv r_5(O(5)), or at least is min_U Φ_{c(x)}(U) = min_O Φ_{c(x)}(O) for the realizable gap vectors c(x), x ∈ X_5? The numerical tests found no gap, but they do not prove equality. If it fails, every real-orthogonal case argument must be redone over U(5). If it holds, a new proof is needed, since the elliptic-curve one has nothing to act on.
2. Can the contact cones N_T of the 8-subsets that are not covered by minor vectors be handled uniformly? That means one curve/maximum-principle argument, the paper's own Section 7 question, using the quadrics of X_5 in place of the cubic. A 29-million-class case split rules out the n = 4 style of proof. Before that, what are the facets of Q_5 (an exact enumeration in dimension 41)?
3. Is K_n = E_n for all n, which would give Lemma 3.4 uniformly? It holds for n = 4 and n = 5. The generic statement is checked here only for those two values.
