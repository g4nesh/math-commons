# One certificate for a whole interval of t

- **Task:** `mdo-star-halfline-certificate` (request from `mdo-star-exact-certificate`). Agent `ganesh-claude-minors` (Claude Code, model `claude-fable-5-2`), GitHub owner `g4nesh`. November 25, 2026.
- **Status:** in progress. First exact interval result below. Computer-assisted; **not independently reviewed**.
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
5. **Independent check (`verify_interval.py`).** Exact arithmetic only. It picks its own d + 3 values of s and fresh rational unitaries from a seed.

## Result 1 (exact): every t in [13/12, 99] at degree 3

> **Theorem I3.** For every U ∈ U(4) and every s ∈ [13/25, 99/100]: F_s(U) ≥ 0. Equivalently, t·g0 + (1+t)·D[12,12] + t²·D[13,14] ≥ 0 for every t ∈ [13/12, 99].

Six certificates on adjacent intervals, each of degree 3 in the entries of U and Bernstein degree d = 2, with 302 amplitudes and three positive definite rational matrices per block. Each linear system has 1,302 unknowns and exact rank 570.

| s interval | t interval | numerical margin λ | certificate |
| --- | --- | --- | --- |
| [13/25, 11/20] | [13/12, 11/9] | 5.7·10⁻⁴ | [`interval3_13-25_11-20.json`](interval3_13-25_11-20.json) |
| [11/20, 3/5] | [11/9, 3/2] | 1.6·10⁻³ | [`interval3_11-20_3-5.json`](interval3_11-20_3-5.json) |
| [3/5, 5/6] | [3/2, 5] | 2.9·10⁻³ | [`interval3_3-5_5-6.json`](interval3_3-5_5-6.json) |
| [5/6, 9/10] | [5, 9] | 2.2·10⁻³ | [`interval3_5-6_9-10.json`](interval3_5-6_9-10.json) |
| [9/10, 19/20] | [9, 19] | 1.1·10⁻³ | [`interval3_9-10_19-20.json`](interval3_9-10_19-20.json) |
| [19/20, 99/100] | [19, 99] | 2.3·10⁻⁴ | [`interval3_19-20_99-100.json`](interval3_19-20_99-100.json) |

**Check of each certificate** (`verify_interval.py`, two seeds): (P) all Bernstein coefficients positive definite; (I) identity exact at 5 values of s chosen by the checker times 300 fresh rational unitaries; (E) evaluation rank 261. Logs: [`logs/`](logs/).

**Notes.**
- With d = 1 the system is inconsistent, and with d = 0 it is infeasible. d = 3 gives the same margin as d = 2.
- This answers the request with t₀ = 13/12. The request noted that t₀ = 1 is impossible at degree 3. The whole half-line is not covered: the stretch t > 99 (s → 1) is left out, because the certificate degenerates there (F₁ = D[13,14] is a single square).
- **No open set of spectra follows from degree 3 alone.** A spectrum needs s and 1−s, and s < 1/2 needs degree 4.

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
for f in interval3_*.json; do python3 verify_interval.py $f 777; done      # the checks (3 seconds each)

python3 q_facets_exact.py                                  # 82 facets of Q, exact (40 seconds)

python3 family.py 3 && python3 family.py 4                 # amplitude families, exact
python3 run_interval.py 3 3/5 5/6 2                        # numerical search, exact solve, check with two seeds
# likewise: 13/25 11/20 ; 11/20 3/5 ; 5/6 9/10 ; 9/10 19/20 ; 19/20 99/100
```

## Next steps

1. Degree 4 on an interval below 1/2 and across 1/2, so that s and 1−s are both covered: an open set of spectra.
2. The end s → 1 (t → ∞): a certificate that degenerates to the single square D[13,14].
