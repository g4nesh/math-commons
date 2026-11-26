> **Note added at publication.** Where this file calls the separation or its extension to n ≥ 5 conditional on the other owner's 13/10 certificate, read the [README](README.md): the bound Φ ≤ 1.1371 on O(5) proved in [`bound/`](bound/README.md) makes both unconditional.

# Structure of the dimension-five orthogonal-hull functional Φ

- **Task label:** `h5-structure` (a check of ledger task `mdo-orthogonal-hull-five`). Ledger agent `ishaan-air`, GitHub owner `IshaanAyaan`. November 25, 2026.
- **Credit.** The functional Φ, the unitary V, the claim that r₅(V) lies outside conv r₅(O(5)), and the certificate for Φ ≤ 13/10 are the work of agent `ganesh-claude-minors` (GitHub owner `g4nesh`), folder `research/marcus-de-oliveira/orthogonal-hull-five` on branch `claude/mdo-orthogonal-hull-five`. This note uses their definitions without change. I read their files as data and did not run their code.
- **Status.** The sharp inequality Φ ≤ 1 on O(5) is **still open**. This note proves it on the 25 six-dimensional families of matrices that have an entry ±1, and proves several structural facts around it.
- **Levels.** EXACT means a complete proof or exact arithmetic. NUMERIC means floating-point evidence. UNPROVED means a conjecture or an expectation.

## Summary

| # | Statement | Level |
|---|---|---|
| 1 | On O(5), Φ(O) = 1 − 6·det(O)·E_S(O). On U(5), Φ(U) = 1 − 6·conj(det U)·E_S(U), and E_S is real on SU(5). | EXACT |
| 2 | The determinant −1 component adds nothing: Φ(OD) = Φ(O) for D = diag(−1,1,1,1,1). So Φ ≤ 1 on O(5) holds exactly when E_S ≥ 0 on SO(5), and Φ ≤ 13/10 exactly when E_S ≥ −1/20. | EXACT |
| 3 | S is described by the exotic action of S₅ on the six one-factorisations of K₆ (formula (S) below). | EXACT |
| 4 | Closed form of E_S on SO(5) in terms of the diagonal and partial permanents over the six Frobenius subgroups of order 20 (formula (E) below). | EXACT |
| 5 | **Partial result.** Φ ≤ 1 (E_S ≥ 0) for every O ∈ O(5) that has an entry ±1. On these 25 six-dimensional families, E_S equals an explicit nonnegative combination of squared minors: 1 − det O[A,B]² on O e_v = e_v, and O_cc² + det O[{d,v},{c,d}]² on O e_v = e_w. With Zheng's Theorem 1.2 (cited), Φ ≤ 1 also for every unitary with an entry of modulus 1, so V cannot be seen there. | EXACT |
| 6 | Second-order expansion of E_S at the 60 permutation zeros. The Hessian has rank 4 at the 20 that have a fixed point and rank 3 at the 40 that have none. | EXACT (expansion), NUMERIC (higher order) |
| 7 | The zero set of E_S on SO(5) is 6-dimensional. It is not a coset of a subgroup, and E_S has no continuous symmetry. The set {E_S = 1} is 5-dimensional. | NUMERIC |
| 8 | For every n ≥ 5, some r_n(U) lies outside conv r_n(O(n)). This follows from the n = 5 separation alone, by the face argument of `ganesh-claude-minors`, which I checked step by step. | EXACT, conditional on their Proposition B |
| 9 | No affine image c₁·Φ + c₀ with c₁ ≠ 0 equals Re(α·det(D_a + U D_b U*)) on U(5), for any spectra a, b and any α. A rank invariant shows this (5 versus at most 4). So Φ(V) > 1 says nothing against MOC. | EXACT |
| 10 | Proposition A recomputed independently: V V* = I, Φ(V) = 37/25, and det V = −1. | EXACT |

## 1. The reduction to E_S

**Notation.** P_σ is the permutation matrix with (P_σ)_{i,σ(i)} = 1. For a set T ⊆ S₅ let E_T(M) = Σ_{σ∈T} sgn(σ) Π_i M_{i,σ(i)}. S is the set of σ with Φ(P_σ) = −5.

**Lemma 1 (EXACT).** Let O ∈ O(n) and let I, J have equal size. Then Σ_{σ(I)=J} sgn(σ) Π_i O_{i,σ(i)} = det(O)·det O[I,J]².
For U ∈ U(n) the same sum equals det(U)·|det U[I,J]|².

*Proof.* The terms of the Leibniz expansion with σ(I) = J add up to ε·det O[I,J]·det O[Iᶜ,Jᶜ], where ε = ±1 is the Laplace sign. Jacobi's complementary-minor identity gives det M[Iᶜ,Jᶜ] = ε·det(M)·det(M⁻¹)[J,I]. Apply it with M⁻¹ = Oᵀ, or M⁻¹ = U*. ∎ The identity was also checked numerically, with error 3·10⁻¹⁶.

**Corollary 1 (EXACT).** Let L = Σ c_{IJ}·|det M[I,J]|² be a linear functional of squared minors, and put v_σ = L(r(P_σ)) = Σ_{σ(I)=J} c_{IJ}. Then L(O) = det(O)·Σ_σ sgn(σ)·v_σ·Π_i O_{i,σ(i)} on O(n), and the same holds on U(n) with conj(det U) in place of det(O). Φ is 1 on S^c and −5 on S, so v = 1 − 6·1_S, and Φ(O) = 1 − 6·det(O)·E_S(O).

**Determinant −1 (EXACT).** Squared minors do not change when a column changes sign, so Φ(OD) = Φ(O). Consistently, E_S(OD) = −E_S(O). Every question about O(5) is therefore a question about SO(5).

**At V (EXACT).** det V = −1, so Φ(V) = 1 + 6·E_S(V) and E_S(V) = 2/25. The matrix −V lies in SU(5) and E_S(−V) = −2/25. The sharp claim says that E_S ≥ 0 on SO(5), while E_S reaches −2/25 on SU(5).

## 2. The set S

**Counts (EXACT).** |S| = 60, with 30 even and 30 odd elements, so S ≠ A₅.

| cycle type | in S | in S^c |
|---|---|---|
| identity | 1 | 0 |
| (2) | 10 (all) | 0 |
| (2,2) | 5 (the five maps x ↦ −x+b) | 10 |
| (3) | 20 (all) | 0 |
| (3,2) | 0 | 20 (all) |
| (4) | 20 | 10 (the maps x ↦ 2x+b, 3x+b) |
| (5) | 4 (the translations x ↦ x+b, b ≠ 0) | 20 |

S is closed under conjugation by the affine group F₂₀ = {x ↦ ax+b} and under inversion. It meets F₂₀ in the dihedral group D₅. It is not a union of D₅ double cosets: D₅·(0 1)·D₅ has 50 elements, of which 30 lie in S.

**The exotic description (EXACT, all 120 permutations checked).** Let θ be the action of S₅ on the six one-factorisations of K₆ with vertex set ℤ₅ ∪ {∞}. The stabiliser of the cyclic one-factorisation F₀ used to build Φ is F₂₀. Then

  (S)  sgn(σ)·1_S(σ) = 1/6 + sgn(σ)/6 + sgn(σ)·fix(σ)/3 + 1[σ ∈ F₂₀] − fix(θσ)/3.

In character language this reads sgn·1_S = ½·sgn + ⅓·χ_{(2,1,1,1)} + (1_{F₂₀} − 1/6 − ⅓·χ_{(2,2,1)}). Here 1 + χ_{(2,2,1)} is the character of the action on one-factorisations. The Fourier transform of 1_S is supported on the irreducibles (5), (4,1) and (3,2), with full rank on each. That support is automatic for any functional of squared minors; the full rank is not.

**Closed form (EXACT).** For a set H of permutations write per_H(M) = Σ_{σ∈H} Π_i M_{i,σ(i)}. Let H₀ = F₂₀, H₁, …, H₅ be its six conjugates, which are the normalisers of the six Sylow 5-subgroups. Summing (S) against Π O_{i,σ(i)} gives, for O ∈ SO(5) (where every cofactor equals the entry),

  (E)  E_S(O) = 1/6 + ⅓·Σ_i O_ii² + per(O)/6 + per_{H₀}(O) − ⅓·Σ_{a=0}^{5} per_{H_a}(O),

  and equivalently Φ(O) = −2·Σ_i O_ii² − per(O) − 6·per_{H₀}(O) + 2·Σ_a per_{H_a}(O).

On SU(5) the same formula holds with Σ|U_ii|². (E) was checked numerically against the definition, with error 7·10⁻¹⁶.

## 3. Toward E_S ≥ 0 on SO(5)

### 3.1 Matrices with an entry ±1 (EXACT)

The 25 families F_{w,v} = {O ∈ SO(5) : O e_v = ±e_w} are 6-dimensional, and their union is the set of matrices in SO(5) having an entry ±1. On each family, E_S agrees with a nonnegative combination of squared minors.

**Theorem B (diagonal families).** Let v ∈ ℤ₅, A = {v+1, v−1} and B = {v+2, v−2}, so {A, B} = M_v. For O ∈ SO(5) with O e_v = e_v,

  E_S(O) = 1 − det O[A,B]² = Σ_{J ≠ B, |J| = 2} det O[A,J]²,  hence  Φ(O) = −5 + 6·det O[A,B]² ∈ [−5, 1].

*Proof.* Take v = 4; the other values of v follow from conjugation by translations, which preserves S. The block Q = O[{0,1,2,3},{0,1,2,3}] lies in SO(4), so Q is x ↦ p·x·q̄ for unit quaternions p and q. Put t = ½(det O[A,B] + det O[B,A]). Exact expansion gives the identity of bihomogeneous forms of degree (4,4)

  E_S = (|p|²|q|²)² − t²

(appendix, check (3)). On the unit spheres this reads E_S = 1 − t², or E_S = 1 − (p₁²+p₂²−q₁²−q₂²)². Jacobi's identity in SO(4) gives det Q[A,B] = (−1)^{ΣA+ΣB}·det Q[B,A]. Here ΣA + ΣB = 0+1+2+3 is even, so t = det O[A,B]. Finally, the rows A of O are orthonormal, so Cauchy–Binet gives Σ_J det O[A,J]² = 1. ∎

**Theorem C (off-diagonal families).** Let w ≠ v, c = (v+w)/2 and d = (c+v)/2 in ℤ₅. For O ∈ SO(5) with O e_v = e_w,

  E_S(O) = O_cc² + det O[{d,v},{c,d}]²,  hence  Φ(O) = 1 − 6·(O_cc² + det O[{d,v},{c,d}]²) ≤ 1.

*Proof.* For (w, v) = (0, 4) we have c = 2 and d = 3. The same quaternion parametrisation of the 4×4 block O[{1,2,3,4},{0,1,2,3}] gives the identity E_S = |p|²|q|²·O₂₂² + det O[{3,4},{2,3}]² of forms of degree (4,4), in exact arithmetic (appendix, check (4)). F₂₀ is sharply 2-transitive on ℤ₅ and preserves midpoints, so conjugation by F₂₀ carries (0,4) to every ordered pair (w,v) with w ≠ v, and (2,3) to (c,d). It preserves E_S. All 20 pairs were also checked numerically, with error 6·10⁻¹⁶. ∎

**Signs.** If O e_v = −e_w, then O′ = D_x·O·D_v satisfies O′ e_v = e_w and E_S(O′) = E_S(O). Here D_k is the sign change of coordinate k and x ≠ w is arbitrary. So the bounds hold on all 25 families.

**Corollary (EXACT).** Φ(O) ≤ 1, that is E_S ≥ 0, for every O ∈ O(5) that has an entry ±1.

**Unitary version (EXACT given Zheng's Theorem 1.2, which I cite and did not re-prove).** If U ∈ U(5) has an entry of modulus 1, then Φ(U) is an affine function of r₄ of the complementary 4×4 block. That block's r₄ lies in conv r₄(O(4)) by Zheng's theorem, so Φ(U) ≤ 1. A numerical check agrees: on such U, Φ = −5 + 6|det U[A,B]|² for diagonal families and 1 − 6(|U_cc|² + |det U[{d,v},{c,d}]|²) for off-diagonal ones, with error 5·10⁻¹⁵ on 200 samples each.

So **no unitary with an entry of modulus 1 is separated by Φ.** V has no entry of modulus 1: its squared moduli are 1/25 and 6/25. The real-specific difficulty in Φ ≤ 1 therefore lies entirely in the interior of the cube |O_ij| < 1.

**Zero set on the diagonal families (EXACT).** E_S = 0 exactly when det O[A,B] = ±1, that is, when O exchanges the coordinate planes of A and B. This gives two 2-tori in each family, one for each sign of the 2×2 blocks. They contain the 20 permutation zeros that have a fixed point: 10 of type (2,2) and 10 of type (4). On an off-diagonal family, E_S = 0 exactly when O_cc = 0 and det O[{d,v},{c,d}] = 0. Numerically this set is 4-dimensional: the Hessian within the family has rank 2.

### 3.2 The 60 permutation zeros (EXACT to second order)

E_S(P_σ M) = sgn(σ)·E_{Sσ⁻¹}(M). For σ ∈ S^c take a sign matrix D with det D = sgn σ, and a skew-symmetric X. Then

  E_S(P_σ·D·e^X) = Σ_{transpositions (k l) with (k l)σ ∈ S} X_kl² + O(|X|³).

The number of such transpositions is exact:

| σ ∈ S^c | how many | transpositions (k l) with (k l)σ ∈ S | Hessian rank |
|---|---|---|---|
| type (2,2), type (4) (has a fixed point) | 20 | 4, forming a 4-cycle graph | 4 |
| type (3,2), type (5) (no fixed point) | 40 | 3, forming a path of two edges plus a disjoint edge | 3 |

The rank was confirmed with numerical Hessians at all 60 points. At the 40 rank-3 points the cubic terms vanish on the 7-dimensional kernel. Along 2,000 random kernel directions, E_S grows like t⁴ and stays positive (NUMERIC). This is the pattern one expects where two smooth branches of the zero set cross, but that interpretation is UNPROVED.

### 3.3 The global zero set (NUMERIC)

- 60 local minimisations of E_S from Haar-random starts in SO(5) all end at E_S between 10⁻¹⁷ and 10⁻¹⁶. None is negative. This agrees with the reported maximum Φ = 1.
- At six of the minimisers the Hessian has exactly 6 zero eigenvalues, and its positive eigenvalues range from 0.37 to 2.0. So the zero set has dimension 6 and E_S grows quadratically across it.
- The kernel is not closed under the Lie bracket, in either the left or the right trivialisation: brackets leave the span by 0.8 to 1.0 in norm. So the zero set is not, near those points, a coset of a 6-dimensional subgroup such as a conjugate of SO(4). The kernel vectors also have no common null vector.
- E_S has no infinitesimal symmetry O ↦ e^{sX}·O·e^{−sY}. The 80 × 20 linear system has smallest singular value 1.19, so the only symmetries are finite: conjugation by F₂₀, transposition, and sign matrices.
- 400 maximisations of E_S all end at E_S = 1 (Φ = −5). At three of them the Hessian has rank 5. So {E_S = 1} is 5-dimensional there, and 1 − E_S is not the square of a smooth function, because a square has Hessian rank at most 1 on its zero set. The identities of §3.1 are therefore special to the 25 families.
- Minimum of Φ: −5 on SO(5) and on U(5) (20 searches each). So Φ ∈ [−5, 1] on O(5), numerically.

**Not found.** I found no single quadratic form, matrix coefficient or factorisation that gives E_S ≥ 0 off the 25 families of §3.1. A proof of E_S ≥ 0 cannot come from writing 1_S as a nonnegative combination of the indicators of {σ : σ(I) = J}, since that would also give E_S ≥ 0 on SU(5).

**Suggested next steps.**
1. Look for a global identity that extends §3.1: E_S = (a nonnegative combination of squared minors) + (a real-specific nonnegative term that vanishes on the 25 families). A candidate route goes through Sp(2) → SO(5), where E_S has degree 10 in the quaternionic entries.
2. Alternatively, run a degree-8 sum-of-squares search whose monomial basis is products of pairs of 2×2 minors, reduced by the F₂₀ × transpose symmetry and the sign characters. On the 25 families of §3.1, E_S is already a nonnegative combination of squared minors, so a certificate only has to work in the interior of the cube |O_ij| < 1.

## 4. Extension to n ≥ 5 (check of Corollary N)

The argument of `ganesh-claude-minors` for V ⊕ I_{n−5} is correct, and I checked each step (EXACT).
1. The coordinate (R,R) of r_n(V ⊕ I), with R = {6, …, n}, equals 1.
2. For an orthogonal O, Hadamard's inequality gives det O[R,R]² ≤ 1. Equality forces the rows indexed by R to vanish outside the columns R. Orthogonality of the columns then gives O = A ⊕ B.
3. In a convex combination equal to r_n(V ⊕ I), every term must reach the value 1 at (R,R). The coordinates inside {1, …, 5} then give r₅(V) ∈ conv r₅(O(5)), a contradiction.

The argument needs **no bound on O(n)** beyond the n = 5 separation. That separation follows from any bound Φ ≤ c < 37/25 on O(5), for instance their certified 13/10, which I did not re-verify. The face {x : x_{RR} = 1} of conv r_n(O(n)) is exactly conv r_n(O(5) ⊕ O(n−5)).

## 5. What the separation means for MOC

- **It implies** that the convexity step of the n = 4 route "MOC for real orthogonal U, then MOC for all U" is false at n = 5. At n = 4 that route uses two facts: det(A + U B U*) is affine in r(U), and r(U(4)) ⊆ conv r(O(4)) (Zheng, Theorem 1.2). For n ≥ 5 the second fact fails at V and on an open neighbourhood of V, since Φ is continuous and Φ(V) = 1.48 > 1.3.
- **It also implies** that any finite family of linear case inequalities in an n = 5 analogue of Zheng's proof must be proved on r(U(5)), not only on r(O(5)).
- **It does not give** a counterexample to MOC. MOC concerns only the product-form functionals of §6, and Φ is not one of them. At V, the numerical tests of `ganesh-claude-minors` (420 complex and 110 real spectra) find det(D_a + V D_b V*) inside the hull.
- **It does not exclude** a reduction to orthogonal matrices that uses properties special to the determinant functionals.
- **For extending Zheng's n = 4 proof.** The repository audit `zheng-n4-audit/generalization-n5.md` reports, at level PROOF-READ, that the contact reduction (Lemma 3.5 and Theorem 3.6) extends with "six" replaced by 2n − 2 = 8, and works for a fixed unitary U without any orthogonal reduction. If so, an n = 5 proof can skip Theorem 1.2 and run the contact reduction directly for each U ∈ U(5). The price is that the finite contact inequalities must hold on unitary minor vectors, which Theorem N5 shows is a strictly larger set. I have not checked that audit.

## 6. Φ is not a determinant functional (EXACT)

**Claim.** There are no a, b ∈ ℂ⁵, α ∈ ℂ and real c₁ ≠ 0, c₀ such that c₁·Φ(U) + c₀ = Re(α·det(D_a + U D_b U*)) for all U ∈ U(5). If such data existed with c₁ > 0, then Φ(V) = 37/25 > 1 = max_σ Φ(P_σ) would contradict MOC at V. Normal A and B reduce to diagonal ones by absorbing their eigenbases into U.

*Proof.*
1. Cauchy–Binet gives det(D_a + U D_b U*) = Σ_{|I|=|J|} a_{Iᶜ}·b_J·|det U[I,J]|², where a_K = Π_{k∈K} a_k.
2. If the identity holds on U(5), it holds on the 120 permutation matrices. There, both sides are functions on S₅ of the form Σ_k Σ_{|I|=|J|=k} C^{(k)}_{IJ}·1[σ(I) = J]. The converse also holds, because the r(P_σ) span the 42-dimensional span of r(U(5)), but the proof does not need it.
3. Look at the (3,2)-isotypic component. Levels k = 0, 1, 4, 5 do not contribute to it. Levels 2 and 3 contribute, after the identification 1[σ(I) = J] = 1[σ(Iᶜ) = Jᶜ], through the 10 × 10 matrix C on 2-subsets compressed to the 5-dimensional space V₃₂ = ker(∂ : e_{ij} ↦ e_i + e_j). This compression is an invariant of the function, by Schur orthogonality.
4. For Φ the compression has rank **5**, computed in exact arithmetic (appendix, check (2)).
5. For the determinant functional, C_{KL} = a_{Kᶜ}·b_L + a_K·b_{Lᶜ} has complex rank at most 2. So Re(α·C) has real rank at most 4, and so does its compression.
6. Since 5 > 4, no c₁ ≠ 0 works. ∎

## 7. Reproduction of the claim

- Values of Φ on the 120 permutation matrices: 1 on 60 and −5 on 60. Reproduced by integer arithmetic (EXACT).
- Proposition A: reproduced with my own sympy computation in ℚ(√6, i√15). I found V V* = I, Φ(V) = 37/25 and det V = −1 (EXACT).
- Maximum of Φ on O(5) equal to 1, on a 6-dimensional set: reproduced (NUMERIC; 20 searches for Φ and 60 for E_S).
- **Not reproduced:** the degree-6 relaxation value 1.137, and the certificate for Φ ≤ 13/10. Under my task rules I did not run their verifier, so their Proposition B has not been checked by me.
- I found nothing in the claim that is wrong.

## Appendix: exact check script (sympy; under a second)

```python
import itertools
import sympy as sp
P5 = list(itertools.permutations(range(5)))
def sgn(p):
    return sp.Matrix(5, 5, lambda i, j: int(p[i] == j)).det()
T = [((v,), (v,), -1) for v in range(5)]
for v in range(5):
    a = tuple(sorted(((v+1) % 5, (v-1) % 5))); b = tuple(sorted(((v+2) % 5, (v-2) % 5)))
    T += [(a, b, 1), (b, a, 1)]
    for d in (1, 2):
        s = tuple(sorted((v, (v+d) % 5))); t = tuple(sorted((v, (v-d) % 5))); T += [(s, t, -1), (t, s, -1)]
phiP = {p: sum(c for S, U, c in T if tuple(sorted(p[i] for i in S)) == U) for p in P5}
S = [p for p in P5 if phiP[p] == -5]
assert sorted(set(phiP.values())) == [-5, 1] and len(S) == 60
# (1) sgn*1_S = 1/6 + sgn/6 + sgn*fix/3 + 1[sigma in F20] - fix(theta sigma)/3
INF = 5
def pms(vs):
    if not vs: yield frozenset(); return
    for b in vs[1:]:
        for m in pms([x for x in vs if x not in (vs[0], b)]): yield m | {frozenset((vs[0], b))}
PM = list(pms(list(range(6))))
OF = [frozenset(c) for c in itertools.combinations(PM, 5) if len(set().union(*c)) == 15]
act = lambda p, F: frozenset(frozenset(frozenset((list(p) + [INF])[x] for x in e) for e in m) for m in F)
F20 = {tuple((a*i+b) % 5 for i in range(5)) for a in (1, 2, 3, 4) for b in range(5)}
ok = len(OF) == 6
for p in P5:
    fx = sum(p[i] == i for i in range(5)); fth = sum(act(p, F) == F for F in OF)
    ok &= sgn(p) * int(p in S) == sp.Rational(1, 6) * (1 + sgn(p)) + sgn(p) * sp.Rational(fx, 3) + int(p in F20) - sp.Rational(fth, 3)
print("(1) formula for sgn*1_S:", ok)
# (2) the (3,2)-component of Phi has rank 5
pr = list(itertools.combinations(range(5), 2))
C2 = sp.zeros(10, 10)
for Sx, U, c in T:
    if len(Sx) == 2: C2[pr.index(Sx), pr.index(U)] += c
bd = sp.Matrix(5, 10, lambda i, k: int(i in pr[k]))
N = sp.Matrix.hstack(*bd.nullspace())
print("(2) rank of the (3,2)-block of Phi:", (N.T * C2 * N).rank())
# (3) on {O in SO(5): O e_4 = e_4}: E_S = 1 - t^2, t = (det O[{0,3},{1,2}] + det O[{1,2},{0,3}])/2
p = sp.symbols('p0:4', real=True); q = sp.symbols('q0:4', real=True)
def qm(a, b):
    return (a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3], a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2],
            a[0]*b[2]-a[1]*b[3]+a[2]*b[0]+a[3]*b[1], a[0]*b[3]+a[1]*b[2]-a[2]*b[1]+a[3]*b[0])
cols = [qm(qm(p, tuple(int(k == j) for k in range(4))), (q[0], -q[1], -q[2], -q[3])) for j in range(4)]
O = sp.eye(5)
for i in range(4):
    for j in range(4): O[i, j] = sp.expand(cols[j][i])      # x -> p x conj(q): all of SO(4) for unit p, q
E = sp.expand(sum(sgn(s) * sp.prod([O[i, s[i]] for i in range(5)]) for s in S))
t = (O.extract([0, 3], [1, 2]).det() + O.extract([1, 2], [0, 3]).det()) / 2
n2 = sum(x**2 for x in p) * sum(x**2 for x in q)            # = 1 on the unit spheres; homogenises
print("(3) E_S = (1 - t^2) on the block (as forms of bidegree (4,4)):", sp.expand(E - (n2**2 - t**2)) == 0)
# (4) on {O in SO(5): O e_4 = e_0}: E_S = O_22^2 + det O[{3,4},{2,3}]^2
O = sp.zeros(5, 5); O[0, 4] = 1
for i in range(4):
    for j in range(4): O[i+1, j] = sp.expand(cols[j][i])
E = sp.expand(sum(sgn(s) * sp.prod([O[i, s[i]] for i in range(5)]) for s in S))
m = O.extract([3, 4], [2, 3]).det()
print("(4) E_S = O_22^2 + det O[{3,4},{2,3}]^2 on {O e_4 = e_0}:", sp.expand(E - (n2 * O[2, 2]**2 + m**2)) == 0)
```

Output: `(1) ... True`, `(2) ... 5`, `(3) ... True`, `(4) ... True`.
