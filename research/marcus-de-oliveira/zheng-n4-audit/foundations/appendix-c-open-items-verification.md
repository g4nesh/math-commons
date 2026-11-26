Evidence level: EXACT for every PASS line in `appC_fill_verify_output.txt` except D4, E2 and E3. Those were recomputed by my own code from (A.1), (A.2), (4.1) and the definitions of Appendix C.1, not by rerunning the audited script. PROOF-READ for the geometric steps of the six claims, each re-argued below in my own words. NUMERIC for D4, E2 and E3 (seeds and sizes are in the output). CITED for (H3) at every stable fixed point (HMSV) and for standard facts. No ERROR was found. Two claims need small repairs, written out in full below.

# Verification of `appendix-c-open-items.md` (label appC-verify)

## Deviations and limits

- **Inputs.** The paper's Appendix C (text pp. 42-58), `appendix-c-part1.md`, `appendix-c-part2.md` and `appendix-c-open-items.md`. I read the audited script `appC_fill_checks.py` to see what it claims. I did not run or import it.
- **The library.** From `../lib/moc4lib.py` I used only the transcribed tables (MATCHINGS, cubic_C, n_vec) and small Fraction helpers.
- **rho~ is built by a different route.** It is defined as the unique functional with <rho~(V), m(z)> = det((z_i - z_{4+j}) V_ij). It is solved at 14 points and certified on the grid {0,1}^8; both sides are multi-affine in z, so the grid check proves the identity. Then it is compared with the complementary-minor definition on all 70 pairs (A1, A2).
- **(H3) is not proved for all stable fixed points.** It is proved generically, and EXACT at six special points. See Section 3.
- **External sources were not re-read.** This covers HMSV, Drury-Cload, Lam and Postnikov; there was no network.
- **Compute.** The final full run took about 2 to 3 minutes, all through RUNJOB. One extra probe, about 3 minutes, is folded into E1.

## 1. Verdicts

| # | Claim (audit section) | Verdict |
|---|---|---|
| 1 | Prop. C.4: (M8')^* = Sigma2, dim Sigma2 = 11; Cor. C.5 (Sec. 1) | CERTIFIED. One gap in a side remark is filled (Repair R1). |
| 2 | G4: F irreducible; alpha's zero locus Z intrinsic; Krull; S8-invariance (Secs. 3.1-3.2), given (H3) | CERTIFIED WITH REPAIRS (R2, R3) |
| 3 | Pluecker divisors irreducible, D_S = D_{S^c}, one S8 orbit (Sec. 3.3) | CERTIFIED (PROOF-READ) |
| 4 | Lemma C.11 without (C.2) (Sec. 4.1) | CERTIFIED (structure PROOF-READ, constants and witnesses EXACT) |
| 5 | Lam replaced by 14 Vandermonde weights (Sec. 5.1) | CERTIFIED (EXACT, plus a 50-digit end-to-end test) |
| 6 | Density step of Lemma C.10; (C.13), (C.14) (Sec. 5.2) | CERTIFIED WITH REPAIR (R4) |

## 2. Claim by claim

### Claim 1: Proposition C.4 and Corollary C.5

**What is needed downstream:**
- the inclusion (M8')^* ⊂ Sigma2, for (C.8) and assembly step 9;
- the equality, for biduality in step 11;
- dim Sigma2-hat = 12, for Corollary C.5 and Lemma C.11.

**Line-by-line check:**
- **Lemma 1.1 (rank of d rho~).** The torus identity is immediate from A1: D_a D1 V D2 + D1 V D2 D_b = D1 (D_a V + V D_b) D2.
  - At V with all entries nonzero, (D1, D2) -> D1 V + V D2 has a 7-dimensional image, which d rho~ sends into C rho~(V). So the kernel has dimension at least 6.
  - Rank <= 10 is a closed condition, so it holds everywhere. Checked at a V with zero entries (A7).
  - Rank exactly 10 at V0 (A3).
- **Lemma 1.2 (conormal criterion).** The function V' -> det(L(V')) with L linear vanishes to second order where L(V) has rank <= 2. Correct.
- **Lemma 1.3 (Segre configurations).** Rederived: A0 D_z = S A0 + E0 with S = e1 e2^t + e3 e4^t, and E0 supported on rows 2 and 4. Correct. Also checked at 3 random points (A6).
- **Prop. 1.4.**
  - The rank <= 2 locus of 4x4 matrices has codimension 4 and is smooth at rank 2, with normal space detected by l^t (.) r. Under (iv), Inc_V is a 4-dimensional submanifold near z.
  - Phi = s m(z') + t m(w') maps into N(V) ∩ Sigma2-hat. If Phi has rank 4, its image contains a Euclidean-open subset of the irreducible 4-dimensional space N(V), so the closed set Sigma2-hat ∩ N(V) is all of N(V).
  - Correct. Inc_V contains the 2-dimensional family s z + t(1,..,1), on which m only scales. So m(Inc_V) is a 3-dimensional cone, the elliptic quartic. Rank 4 is consistent with that.
  - Conditions (i)-(v) hold at the audit's point (A3) and at a second random rational point (A5).
- **Thm. 1.5, inclusion.**
  - Image of Omega x C^* dense (A4, rank 10).
  - At p = lambda rho~(V) in D0: im d rho~_V (dimension 10) ⊂ T_p M-hat (dimension 10, p smooth; the cone has the same tangent at p and p/lambda). So the conormal fibre is P(N(V)).
  - Con(M8') is the closure of a bundle over the irreducible smooth locus, so its part over the dense open P(D0) is dense. pr2 of a closure lies in the closure of pr2.
  - Correct.
- **Thm. 1.5, dimension.**
  - E is a rank-4 bundle over U0, of dimension 20. Its part over the dense open set where rho~(V) is smooth maps into the cone over (M8')^*, so pr2(E) does too.
  - The projection has rank 12 at (V0, xi0), recomputed with my own xi0 (A10). Generic rank >= 12.
  - So 11 <= dim (M8')^* <= dim Sigma2 <= 11. Both are irreducible, hence equal.
  - Correct.
- **A shorter proof of dim Sigma2-hat = 12 (A9).** The differential of (lambda, z, mu, w) -> lambda m(z) + mu m(w) has rank 12 at one rational point. Corollary C.5 needs only this, so it does not depend on the dual-variety argument at all.
- **Corollary C.5.** G maps a dense subset of Sigma2-hat into M-hat (by (C.4): G(m(z) + m(w)) = det(A0)^2 rho~(V)). So the generic rank of H(c) restricted to T_c Sigma2-hat is at most 10, which gives dim ker H(c) >= 2 at generic points, and then everywhere by closedness. Correct.

**Consequence tests.**
- C vanishes identically on N(V0) (a cubic in 4 variables).
- rank H = 12 at 3 random points of N(V0) (A8).
- Every intermediate claim holds at both points.

**R1 (gap in the remark on Lemma C.7, audit Sec. 1, last bullet).** The audit asserts that smooth points of M-hat are dense on the immersed 10-dimensional germ through rho(0, epsilon v), without a reason.
- **Reason.** Sing(M-hat) has dimension <= 9, and the germ is a 10-dimensional complex manifold. So the germ is not contained in the singular locus, and its intersection with the smooth locus is open and dense in the germ.
- **The rest follows.** At such points T M-hat = span[rho, d_x rho, d_z rho]. Left kernels converge because the rank is constant at 10 nearby, and Con is closed.
- So Lemma C.7 holds for the image variety M8' under (H3) alone; injectivity of the quotient map is not needed.

### Claim 2: G4, given (H3)

**F irreducible: correct.**
- For stable configurations, equality in M8-bar means one SL4 x (C^*)^8 orbit. So Lambda^⊥ = Lambda D with D an invertible diagonal matrix, and Lambda is isotropic for diag(d).
- Rescaling columns by sqrt(d_k) gives Lambda in OG(4,8) for q0 = sum x_k^2. A sign change of one column swaps the two families.
- So F is the image of the open, nonempty, stable part of the irreducible OG+(4,8).

**R2 (cosmetic).**
- Step 3.1(2) ("all d_k != 0") is redundant, because D is a torus element.
- Its bound "a plane contains at most six points of a stable configuration" should read at most five. Hilbert-Mumford for 8 points in P^3 gives fewer than 2(r+1) points on an r-plane; six is the semistable bound.
- The conclusion is unaffected.

**Z intrinsic, local equation, Krull: correct.**
- (H3) is equivalent to: rank[rho0, A] = 7 and S_O maps K_O onto Sym^2. In particular K_p is a rank-7 bundle on F°.
- Near p, alpha = (a coefficient of C|K) / (a coefficient of det S) in a regular frame, with the denominator nonzero at p.
- F° is irreducible as an open subset of F. alpha is nonzero at one point, so its zero set is not all of F°. Krull's principal ideal theorem applies in any Noetherian local domain, so smoothness of F is not needed.

**R3 (missing step in S8-invariance).**
- Z is defined only on F°, so its S8-invariance also needs F° to be S8-invariant. The audit does not say this.
- **Proof.** Relabelling commutes with Gale duality (permuting columns commutes with taking Lambda^⊥). The complementary-product map is S8-equivariant up to a linear action on V_C^*. The cotangent-generation condition (H3) is intrinsic. Hence sigma(F°) = F°.

**Conclusion (alpha(O) != 0 for d(O) != 0).** This holds only for O in F°, as the audit says.

**Exact consistency (E1).** (C.8) holds as an identity of cubics on the 7-dimensional K_O at five real points, and dim K_O = 7 at each:
- the Lemma C.9 example;
- the Cayley matrix of part 1;
- a random rational O;
- a random O with O00 = 0;
- the indecomposable zero-entry point with d = 0, where alpha = 0.

**New observation (EXACT at 4 points, not proved).**
- In my phase frame F01, F02, F03, the 10x10 frame determinant equals -d(O). So that frame degenerates exactly on the Cauchy divisor.
- alpha_frame equals 1/d(O) at all four points with d != 0 (E1 output).
- Normalising by the frame, alpha * framedet^2 = d(O) at all four.
- This suggests the identity C|K_O = d(O) det S^_O, with S^ the phase Hessian in a gauge-normalised frame. If proved, it would give Lemma C.9 directly, with no purity or Pluecker-divisor argument. Recorded as a conjecture only.

### Claim 3: Pluecker divisors

All six steps are correct:
1. Delta_S(Lambda) = 0 iff Lambda ∩ span(e_{S^c}) != 0.
2. (Lambda + W0^⊥)^⊥ = Lambda ∩ W0. So Lambda ∩ W0 != 0 iff Lambda ∩ W0^⊥ != 0, hence D_S = D_{S^c}. Consistent with Jacobi in the chart.
3. The incidence variety I fibres over the smooth quadric surface of isotropic lines in W0, with fibres OG+(3,6) ≅ P^3. A surjective closed map onto an irreducible base with irreducible fibres of constant dimension has an irreducible source (standard). So I is irreducible of dimension 5.
4. Its image is irreducible, and it is cut from the 6-dimensional OG+ by one equation, so it has pure dimension 5.
5. The stable part is open in it, and the map to F has finite fibres (the self-dual relation is unique up to scale and signs).
6. S8 is transitive, after a column sign change.

PROOF-READ.

**Completion of Lemma C.9.** The complement of {chart, all minors nonzero, d != 0} is the union of the D_S and the closure of {d = 0}. Each component of Z is therefore one D_S ∩ F° or lies in {d = 0}. The example point is in F°: (H3) holds there (EXACT, E1), alpha != 0 there, and it lies on D_S for O00 = 0. That excludes every D_S. Correct.

### Claim 4: Lemma C.11 without (C.2)

**Structure (PROOF-READ, with B1 and B2 checked exactly):**
- C is linear in u13 with coefficient q, so H(kappa) = Hess q, supported on W × W, and column 13 of H(u) is supported on W.
- **h^4 term.** It is det Hess q[W,W] · det H[K,K], and H[K,K] has a zero column 13.
- **h^3 term.** Only partner pairs (s, t) contribute. The double expansion along column and row 13 gives ±d_s q d_t q det H[K',K'].
- So p3 = det H[K',K'] (alpha u0 u6 + beta u1 u5).

**Constants.** alpha = -2 and beta = 2 hold at 6 random rational u with my own seed (B3). Given the proved form, that is a proof.

**Prop. 4.2.**
- q(y) = 0; rank H(y) <= 12 kills p0, and adj H(y) = 0 kills p1.
- So det H(y + h kappa) = p2 h^2, and y + h kappa is not in Sigma2-hat for h != 0. Correct.

**Prop. 4.3.**
- The Zariski tangent space of Z_red has dimension <= 12, and dim O_{Z,y} >= 12. So O_{Z,y} is regular, hence a domain. So Sigma2-hat = Z near y, and y is smooth. Correct.
- At the example, dC_y together with the differentials of all 196 rank-12 minors spans rank exactly 2 (B5). The audit's three minors suffice.

**Openness.** Both conditions are polynomial in (a, Re b, Im b, theta) and hold at the example. They also hold at 3 random parameter points (B6). They fail, as expected, at non-generic points: a1 = a2, b real, a1 + b1 = 0 (B7). Those are recovered by the paper's continuity step.

### Claim 5: Lam replaced

**Exact checks (C1-C4):**
- 14 noncrossing matchings, each joining odd to even positions.
- The 70 x 14 compatibility matrix has rank 14, so the weights are unique.
- All 70 identities hold as polynomials in the gaps. I solved with rows chosen in reverse order.
- All weights are nonzero with nonnegative coefficients. The term counts are identical to the audit's.

**The use in Lemma C.10 (PROOF-READ):**
- |Delta_S Delta_{S^c}([I|O])| = det O[I,J]^2, by Jacobi.
- Columns equal A v(t_k) c_k exactly, so the ratio of Vandermonde products is exact.
- Compatibility of a matching with S = (S0\I) ∪ J is equivalent to sigma_pi(I) = J.
- Distinct matchings give distinct bijections.

**50-digit end-to-end test (C5).** 12 Cauchy orthogonal matrices built from interlaced nodes, with both parities and one node gap of 1e-6. In all 70 identities, sum over sigma(I) = J of p_sigma equals det O[I,J]^2 to within 3e-50.

### Claim 6: density step, (C.13), (C.14)

**(C.13).** An exact polynomial identity on a bijective rational parametrization of SO(4) ∩ {O00 = 0}, namely O = diag(1, -Q) H(u) with u stereographic and Q Cayley (D2).
- This replaces the audit's informal remark that "the parametrization omits only a 4-dimensional set".

**(C.14).** An exact identity in my own parametrization, with common sign -1 (D3).

**R4 (missing step).** Step 4 of the corrected density argument uses (C.13) "transported to (i,j)" without proof.
- **Proof.** Take the first identity of (C.9) at (i,j) = (1,1). With O00 = 0, the reciprocal-minor expansion keeps only the permutations fixing (0,0), so r(O-hat_{11}) = O02 O03 O20 O30 (O23 O32 - O22 O33).
- For O in SO(4), Jacobi gives O22 O33 - O23 O32 = O00 O11 - O01 O10 = -O01 O10. So d = O01 O02 O03 O10 O20 O30.
- A signed permutation P1 O P2 in SO(4) moves (i,j) to (0,0), and d is preserved there with sign +1. Checked EXACT at all 16 positions and 10 random matrices (D1).

**The Baire argument.** It is correct as written:
- O' near O is still indecomposable.
- H and E_ij coincide near O'.
- Germs of real-analytic functions form a domain.
- Two zero-entry divisors cannot coincide at an indecomposable point.

**NUMERIC observation (D4).** At the zero-entry point with row (3/5, 4/5, 0, 0):
- 20 tangent steps projected back to {d = 0} give all entries nonzero, but O02 grows only quadratically. So H is tangent to {O02 = 0} there.
- This does not affect the argument, which needs only density, but it is why naive first-order tests look marginal.

## 3. (H3)

**Reduction.** By S8-equivariance and R3, (H3) holds on all of F if and only if it holds at every stable [I4, O] with O in SO(4, C). Every point of F has a nonzero Pluecker coordinate, and a relabelling moves it into the chart.

**Generic (H3) is proved.** F is irreducible (claim 2), the failure locus is Zariski closed, and rank[rho0, A, B] = 13 holds exactly at a point.

**EXACT points (E1).** (H3) holds at six points, all with gauge dimension 7:
- the Lemma C.9 example, which lies on 7 Pluecker divisors;
- the part-1 Cayley matrix;
- a random rational O;
- a random O with O00 = 0;
- the zero-entry point with d = 0;
- a complex Gaussian-rational O in SO(4, C).

**NUMERIC evidence.**
- Over 1000 Haar real O (E2), the minimum smallest singular value of the column-normalised 14x13 matrix is 5.5e-3, attained only where min |O_ij| < 0.02, close to decomposable matrices where stability degenerates. Over the 120 samples with all |O_ij| >= 0.1 the minimum is 6.5e-2.
- Local minimisation (E3) with all |O_ij| >= 0.1: minima at least 0.117.
- No real near-failure was seen. Complex failures cannot be excluded this way.

**Proof for all points: not achieved.**
- **Route 1.** Show the 14 maximal minors of the 14x13 matrix have no common zero on stable SO(4, C). This is an elimination problem in degree about 100 in 6 quaternion parameters, out of budget here.
- **Route 2.** Rely on the closed embedding of M8-bar/<g> (HMSV Props 8.1-8.2, degree-one generation of R^g). That implies cotangent generation at every point, which is exactly (H3). The audit's derivation of (H3) from that citation is correct.

**Is generic (H3) enough?**
- **Yes for Lemma C.9's purity argument,** which lives on the open set F°.
- **No for the assembly as written.** (H3) must hold at the actual maximiser O for (a, b, theta) in a dense set. A repair would have to show that, for a dense set of parameters, no real critical point that is stable, has d != 0 and fails (H3) exists.
- The naive count is favourable: the failure locus has real dimension <= 5; for fixed O the critical y form an 8-dimensional affine space; the parameter set Y = {Re(theta m(a, -b))} has real dimension 11; and 5 + 8 + 11 - 14 = 10 < 11.
- But where rank A drops, K_O grows, and the count fails. Not done. So (H3) at all stable fixed points stays CITED.

## 4. What Theorem C.1 rests on after this check

| Input | Status |
|---|---|
| (C.3)-(C.7) identities, Moebius rule | EXACT (parts 1 and 2; A1 and A2 here); (C.7) modulo Theorem 1.2 |
| Lemma C.2, Lemma C.3, Lemma C.8, Lemma C.7 (abstract), stability and tangent decomposition, assembly | proved in the paper (PROOF-READ by parts 1 and 2) |
| Prop. C.4 (equality), dim Sigma2 = 11, Cor. C.5 | proved by the audit, checked here (CERTIFIED; R1) |
| (C.8) | proved given (H3), Lemma C.7 and Prop. C.4; EXACT at 5 points |
| F irreducible, Z a divisor, Pluecker divisors, Lemma C.9 | proved by the audit given (H3), checked here (R2, R3) |
| Lemma C.10 (weights, density) | proved by the audit, checked here (R4); no Lam or Postnikov needed |
| Lemma C.11 | proved by the audit without (C.2), checked here |
| Theorem 1.2 (main text) | CITED (audited elsewhere) |
| (H3) at every stable Gale-fixed point | CITED: HMSV Props. 8.1-8.2, as quoted on paper p. 44, "the Gale-invariant subring is generated by the degree-one invariants". Proved here only generically; OPEN without HMSV |
| Standard facts | CITED: projective biduality in characteristic 0; Chevalley; Krull; GIT smoothness at stable points with trivial stabiliser and Luna-slice equivariant coordinates; irreducibility of a closed surjection with irreducible equidimensional fibres; Baire |

No part of (C.2) beyond X ⊂ Sing C (EXACT) is used by Theorem C.1. The four supplement case documents that use HMSV's set equality (audit Sec. 4.2) are outside this check.

## 5. Files

- `appC_fill_verify.py`: sections A-E, deterministic.
- `appC_fill_verify_output.txt`: the captured run.
