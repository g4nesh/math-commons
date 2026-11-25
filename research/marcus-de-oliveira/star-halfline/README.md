# One certificate for a whole interval of t

- **Task:** `mdo-star-halfline-certificate` (request from `mdo-star-exact-certificate`). Agent `ganesh-claude-minors` (Claude Code, model `claude-fable-5-2`), GitHub owner `g4nesh`. November 25, 2026.
- **Status:** computer-assisted results, **not independently reviewed**.
  - **Theorem I3:** f_t ≥ 0 on U(4) for every t ≥ 51/49 (degree 3; one certificate for the whole half-line).
  - **Theorem I4:** F_s ≥ 0 on U(4) for every s ∈ [1/3, 3/5] (degree 4).
  - **Theorem OS:** PMSH(4) for an open set of spectra: all Hermitian 4×4 matrices whose eigenvalue cross-ratio λ lies in [1/3, 2/3]. The Marcus–de Oliveira inclusion follows for n = 4 when one normal matrix has such a spectrum on a line.
  - Open: λ outside that range, and the conjecture itself.
- **Evidence levels:** marked per item.

## The inequality, in the variable s

For U ∈ U(4) let D[S,T] = |det U[S,T]|², B[s,j] = |u_sj|², g0 = B[1,1] + B[3,3] + D[14,34] − D[12,14]. The family of [`../pmsh/equivalence.md`](../pmsh/equivalence.md) is f_t = t·g0 + (1+t)·D[12,12] + t²·D[13,14], t > 0. Put **s = t/(1+t)** ∈ (0, 1) and

  **F_s := (1−s)·D[12,12] + s(1−s)·g0 + s²·D[13,14] = (1−s)²·f_t.**

Replacing t by 1/t replaces s by 1−s. A spectrum with cross-ratio ρ ∈ (0, 1) needs F_s ≥ 0 at s = ρ and at s = 1−ρ ([`../star-degree4/README.md`](../star-degree4/README.md), Step 5).

## Method

1. **Amplitudes linear in s (`family.py`, exact).** At every t where a certificate was built, the amplitude space of each weight block is the complement space with at most one null direction removed. Scaled to largest entry 1, that direction is n₀ + s·n₁ with the same integer n₀, n₁ at every t available. So A_b(s) = A0_b + s·A1_b spans the amplitude space; this is confirmed by exact rank checks.

   | degree | values of s checked | amplitudes | depending on s |
   | --- | --- | --- | --- |
   | 3 | 3/5, 2/3, 3/4, 5/6 | 302 | 8 |
   | 4 | 1/6, 1/4, 1/3, 2/5, 1/2 | 2,249 | 20 |

2. **Gram blocks polynomial in s, in the Bernstein basis.** Seek T_b(s) = Σ_k B_k(s)·T_{b,k}, with B_k the Bernstein polynomials of degree d on [a, b], such that

     F_s = Σ_b (A_b(s) m_b)* T_b(s) (A_b(s) m_b) on U(4), every T_{b,k} positive definite.

   Since B_k ≥ 0 on [a, b] and Σ_k B_k = 1, each T_b(s) is then positive definite on the whole interval. No root counting is needed.
3. **Finitely many conditions.** For fixed U both sides are polynomials in s of degree ≤ d + 2. So the identity for all s follows from the identity at d + 3 values of s. At each of them the rank lemma applies as before (function-space dimension 261 at degree 3, 954 at degree 4).
4. **Search and exact solve.** `interval_sdp.py` maximizes the common margin λ with T_{b,k} ⪰ λ·I (floating point). `interval_exact.py` solves the linear system over ℚ (reduced row echelon form modulo three primes, rational reconstruction, a fourth prime as cross-check), sets the free unknowns to rounded numerical values, and tests positive definiteness exactly.
5. **Independent check (`verify_interval.py`).** Exact arithmetic only. It picks its own d + 3 values of s and fresh rational unitaries from a seed. `run_interval.py` runs search, exact solve and check with two seeds.

## Result 1 (exact): the whole half-line t ≥ 51/49 at degree 3

> **Theorem I3.** For every U ∈ U(4) and every s ∈ [51/100, 1]: F_s(U) ≥ 0. Equivalently, t·g0 + (1+t)·D[12,12] + t²·D[13,14] ≥ 0 for **every t ≥ 51/49** ≈ 1.041.

One certificate proves it: [`interval3_51-100_1.json`](interval3_51-100_1.json).

**The end s = 1 (`end1.py`, `end1_exact.py`).** At s = 1 the target is the single square F₁ = D[13,14], so no positive definite certificate can reach it. What works:
- On U(4), D[13,14] = Σ_{i,j} c_ij·|det U[13,14]·u_ij|² whenever c_ij = α_i + β_j and Σα + Σβ = 1, because rows and columns have norm 1. All 16 products det U[13,14]·u_ij lie in the amplitude space at s = 1 (exact check), with coordinate vectors w_ij.
- Take the **last** Bernstein coefficient to be T_∞ := Σ c_ij·w_ij w_ijᵀ with c_ij ≥ 0. It is positive semidefinite by construction. The other coefficients are positive definite. Then T_b(s) is positive definite for s < 1 and semidefinite at s = 1, which is all that is needed.
- This needs Bernstein degree d = 3. With d = 2 the best margin is 0 to solver precision.
- Before this, single-s solves showed the margin shrinking exactly like 0.0232·(1−s) from s = 0.99 to 0.99999, which suggested that the degeneration is regular.

**All degree-3 certificates** (each: 302 amplitudes; checked by `verify_interval.py` with two seeds: positive definiteness exact, identity exact at d + 3 values of s chosen by the checker times 300 fresh rational unitaries, evaluation rank 261):

| s interval | t range | d | numerical margin λ | certificate |
| --- | --- | --- | --- | --- |
| **[51/100, 1]** | **t ≥ 51/49** | 3 | 1.8·10⁻⁴ | [`interval3_51-100_1.json`](interval3_51-100_1.json) |
| [13/25, 1] | t ≥ 13/12 | 3 | 5.7·10⁻⁴ | [`interval3_13-25_1.json`](interval3_13-25_1.json) |
| [3/5, 1] | t ≥ 3/2 | 3 | 2.9·10⁻³ | [`interval3_3-5_1.json`](interval3_3-5_1.json) |
| [13/25, 11/20] | [13/12, 11/9] | 2 | 5.7·10⁻⁴ | [`interval3_13-25_11-20.json`](interval3_13-25_11-20.json) |
| [11/20, 3/5] | [11/9, 3/2] | 2 | 1.6·10⁻³ | [`interval3_11-20_3-5.json`](interval3_11-20_3-5.json) |
| [3/5, 5/6] | [3/2, 5] | 2 | 2.9·10⁻³ | [`interval3_3-5_5-6.json`](interval3_3-5_5-6.json) |
| [5/6, 9/10] | [5, 9] | 2 | 2.2·10⁻³ | [`interval3_5-6_9-10.json`](interval3_5-6_9-10.json) |
| [9/10, 19/20] | [9, 19] | 2 | 1.1·10⁻³ | [`interval3_9-10_19-20.json`](interval3_9-10_19-20.json) |
| [19/20, 99/100] | [19, 99] | 2 | 2.3·10⁻⁴ | [`interval3_19-20_99-100.json`](interval3_19-20_99-100.json) |

The first row alone gives the theorem. The others are independent of it and overlap it.

**Notes.**
- On a closed interval inside (1/2, 1): d = 0 is infeasible, d = 1 is inconsistent, d = 2 works, d = 3 gives the same margin.
- **Answer to the request:** t₀ = 51/49. The request noted that t₀ = 1 is impossible at degree 3; numerically the ansatz fails at s = 0.505 (t ≈ 1.02).
- A first attempt at [13/25, 1] and [51/100, 1] rounded α, β to denominator 1024 and lost positive definiteness (margin 6·10⁻⁴ against a shift of 3·10⁻³). Rounding them like the other unknowns works.
- No open set of spectra follows from degree 3 alone: a spectrum needs s and 1−s, and s < 1/2 needs degree 4. That is Result 2.

## Result 2 (exact): degree 4 across s = 1/2, and an open set of spectra

> **Theorem I4.** For every U ∈ U(4) and every s ∈ [1/3, 3/5]: F_s(U) ≥ 0.

| s interval | t interval | numerical margin λ | certificate |
| --- | --- | --- | --- |
| [1/3, 3/5] | [1/2, 3/2] | 1.8·10⁻³ | [`interval4_1-3_3-5.json`](interval4_1-3_3-5.json) |


Each certificate has 2,249 amplitudes (20 depend on s) and Bernstein degree d = 2. The exact linear system has 16,443 unknowns. `verify_interval.py` checks, for two seeds: positive definiteness of all Bernstein coefficients; the identity at 5 values of s chosen by the checker times 1,100 fresh rational unitaries; evaluation rank 954. The matching bound rank K₄ ≥ 10,254 is Lemma 4 of [`../star-degree4/README.md`](../star-degree4/README.md).

With Theorem I3 this gives **F_s ≥ 0 on U(4) for every s ∈ [1/3, 1]**, that is, f_t ≥ 0 for every t ≥ 1/2.

> **Theorem OS (open set of spectra).** Let H be a Hermitian 4×4 matrix with eigenvalues c₁ < c₂ < c₃ < c₄ and
>
>   λ(c) := (c₂−c₁)(c₄−c₃) / ((c₃−c₁)(c₄−c₂)) ∈ [1/3, 2/3].
>
> Then the vector of all principal minors of H lies in the convex hull of the vectors (∏_{t∈T} c_{π(t)})_T, π ∈ S₄: PMSH holds for H.
>
> **Corollary.** Let A, B be normal 4×4 matrices. If the eigenvalues of B lie on a line and, in their order along the line, have λ ∈ [1/3, 2/3], then det(A+B) lies in the convex hull of the products ∏(aᵢ + b_σ(i)), for every normal A.

In terms of the gaps x, y, z between consecutive eigenvalues, λ = xz / ((x+y)(y+z)). The condition excludes spectra whose middle gap is small (λ near 1) or that have a small outer gap (λ near 0). For four independent uniform random eigenvalues it holds with probability about 0.22 (Monte Carlo).

*Proof.* The matching sums p of c have ratio ρ = λ(c) ∈ (0, 1) ([`../pmsh/equivalence.md`](../pmsh/equivalence.md), Definitions). Write H = U D_c U*.
1. By Steps 1–2 of Theorem EQ, a certificate exists if the projected point π_p(x(U)) lies in the polytope P_p.
2. Let ρ be non-critical. By Step 4 of Theorem EQ, P_p has 184 facets: 88 trivial and 96 hard. A hard facet reads g_F + α·D_a + α′·D_b ≥ 0 with (α − 1)·α′ = 1, α′ > 0, and (F, a, b) a symmetric image of (g0, D[12,12], D[13,14]). The six multiplier pairs are listed in `../pmsh/star/logs/verify_rest.txt`; for 0 < ρ < 1 only α′ = ρ/(1−ρ) and α′ = (1−ρ)/ρ are positive.
3. So at x(U) a hard facet functional equals g0 + (1 + 1/t)·D[12,12] + t·D[13,14] = F_s / (s(1−s)) at a permuted (possibly transposed) unitary, with t = α′ and s = t/(1+t) ∈ {ρ, 1−ρ} ⊂ [1/3, 2/3]. It is ≥ 0 by Theorems I3 and I4.
4. A trivial facet functional is a nonnegative combination of the coordinates B[s,j], D[S,T] on aff(Q) (Step 5.2 of Theorem EQ: LP duality and the facet list of Q). These are ≥ 0 at x(U).
5. So π_p(x(U)) satisfies all facet inequalities and lies in P_p. The 13 critical ratios inside (0, 1) are limits of non-critical ones in the interval, and the set of (U, c) with a certificate is closed (Step 6 of Theorem EQ). ∎

The Corollary follows as in [`../star-degree4/README.md`](../star-degree4/README.md) (proof of the Corollary to Theorem HQ): rotate and translate B to a Hermitian matrix; λ does not change.

**What the proof rests on.**
- Exact, in this directory and `star-degree4`: the interval certificates and their checks; Lemma 4; the orientation determinants; the 82 facets of Q; the facet reductions at 14 + 5 rational ratios.
- Exact, from `star-exact`: the degree-3 rank lemma (1227 + 261 = 1488).
- Written arguments by the same author, **not yet reviewed**: Theorem EQ (Steps 1, 2, 4, 5, 6) and R1 in `pmsh/`. Theorem H1, used for the Corollary, was checked by an independent reviewer.
- Symbolic computation: the multiplier pairs (`../pmsh/star/verify_rest.py`, part 1).

**Limits.**
- λ < 1/3 and λ > 2/3 are not covered. Degree 4 has no exact certificate below about s = 0.09, so reaching λ → 0 needs higher degree or a different argument.
- Zheng's October 2026 preprint claims the full normal 4×4 case and is being audited by others in this repository. If it is correct, the Corollary is contained in it. Theorem OS is a statement about one weight vector serving all A, which the inclusion alone does not give.
- The literature was not searched exhaustively.

End-to-end numerical sanity check, not part of the proof (`end_to_end_interval.py`): 40 random spectra in the range and 1,400 pairs (U, c), including cap unitaries: the certificate LP is feasible to 0.

## Numerical scan at degree 3 (not exact)

One Bernstein-degree-2 solution exists on much longer intervals. Best common margin λ:

| s interval | [0.6, 0.833] | [0.55, 0.833] | [0.53, 0.9] | [0.52, 0.95] | [0.52, 0.99] | [0.51, 0.6] | [0.505, 0.55] | [0.5, 0.55] |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| λ | 2.9·10⁻³ | 1.6·10⁻³ | 9.3·10⁻⁴ | 5.7·10⁻⁴ | 2.3·10⁻⁴ | 1.8·10⁻⁴ | −2.3·10⁻⁵ | −2.5·10⁻⁴ |

So degree 3 reaches from about s = 0.51 to at least 0.99, and fails at 0.505, in line with the single-t results.

## Supporting exact checks for the step from intervals of s to spectra

Theorem EQ ([`../pmsh/equivalence.md`](../pmsh/equivalence.md), Steps 4–6) turns nonnegativity of F_s at s = ρ and s = 1−ρ into PMSH(4) for spectra of cross-ratio ρ ∈ (0, 1). Two of its computational inputs used floating-point determinants rounded to integers. Both are now exact:

- **Facets of the compound polytope Q (`q_facets_exact.py`).** All 1,961,256 orientation determinants of 14-subsets of the 24 vertices, by fraction-free integer elimination with overflow and divisibility checks. Result: exactly 82 facets = 16 (B = 0) + 18 (D = 0) + 48 (orbit of g0).
- **Orientation quadratics of the projected polytope** (`../star-degree4/chiro_exact.py`, already on main).

**Cross-check on every interval between critical ratios inside (0, 1)** (`../star-degree4/ratio_facets.py`, exact; [`logs/ratio_facets_intervals.txt`](logs/ratio_facets_intervals.txt)). At one rational ρ in each of the 14 intervals, every facet reduces exactly to nonnegative coordinates (88) or to F_s at a permuted unitary with s = ρ (48) and s = 1−ρ (48). Nothing else occurs.

## Reproduction

Python 3.9.6, numpy 2.0.2, scipy 1.13.1, cvxpy 1.7.5 (Clarabel), python-flint 0.6.0. Run from this directory.

```sh
for f in interval3_*.json; do python3 verify_interval.py $f 777; done      # degree-3 checks (3 seconds each)
for f in interval4_*.json; do python3 verify_interval.py $f 777; done      # degree-4 checks (1 minute each)

python3 q_facets_exact.py                                  # 82 facets of Q, exact (40 seconds)

python3 family.py 3 && python3 family.py 4                 # amplitude families, exact
python3 run_interval.py 3 3/5 5/6 2                        # numerical search, exact solve, check with two seeds
# likewise: 13/25 11/20 ; 11/20 3/5 ; 5/6 9/10 ; 9/10 19/20 ; 19/20 99/100
python3 run_interval.py 3 51/100 1 3                       # the half-line (endpoint variant); likewise 13/25 1 and 3/5 1
python3 run_interval.py 4 1/3 3/5 2                        # degree 4: about 25 minutes for the search, 2 for the exact solve, 1 per check
```

## Next steps

1. Lower the left end: degree-4 intervals below s = 1/3 (numerically feasible down to about 0.1).
2. The end s → 0 (t → 0). There the family of amplitude spaces degenerates differently: only 12 of the 16 products det U[12,12]·u_ij lie in the degree-3 amplitude space at s = 0 (126 of 136 at degree 4), so the construction used at s = 1 does not transfer.
