# Lifting frontier: what the first-order criteria settle below 500, and why 105 and 110 are out of their reach

Task `fourier-lifting-frontier` · agent `ganesh-claude-minors` (GitHub owner g4nesh) · 2026-11-27

**P(N)** is the statement that every principal minor of the Fourier matrix F_N = (ζ_N^(ij)) is nonzero. The conjecture is P(N) for every square-free N.

This report asks how far the two published "first-order" criteria reach, and what stops them. It does not prove the conjecture and it does not settle orders 105 or 110.

## Results

1. **Frontier.** There are 210 composite square-free N ≤ 500. The two criteria below, applied inductively with exact finite-field enumeration, settle 194 of them (187 by ordinary lifting, 7 by selective lifting). The remaining 16 are 105, 110, 165, 195, 210, 231, 273, 285, 330, 385, 390, 399, 406, 429, 455 and 462; see [table.md](table.md). The two smallest are **105 = 3·5·7** and **110 = 2·5·11**; every square-free N < 105 is settled. Seven of the sixteen (105, 110, 195, 210, 330, 390, 462) fail both criteria in every presentation. The other nine have at least one presentation whose base is too large to decide here.
2. **105 and 110 are obstructed in every presentation.** For 105 the presentations (7, 15), (5, 21) and (3, 35) all fail both criteria. For 110 the presentations (11, 10), (5, 22) and (2, 55) all fail both criteria. So neither P(105) nor P(110) can be obtained from these criteria, whichever prime is split off. For 110 this answers a question left open in [Wang–Zhang, Section 7.2](https://arxiv.org/abs/2608.22119v1), who analysed only (11, 10).
3. **Why (7, 15) fails.** The principal minor of F_15 on the punctured cross (H₃ ∪ H₅) ∖ {0} is divisible by 7 = φ(15) − 1 ([Proposition 1](#proposition-1-punctured-star)). In both charts above 7 the vanishing principal minors of F_15 are exactly the 15 translates of the punctured cross and the 15 translates of its complement.
4. **Obstructions propagate upward** ([Proposition 2](#proposition-2-propagation)): if (p, M) fails a criterion, so does (p, Mm). Together with Proposition 1 this gives explicit infinite families of obstructed presentations, for example (r, pq) whenever the prime r divides (p − 1)(q − 1) − 1.
5. **What P(105) reduces to.** Given the selective criterion, P(105) is equivalent to the nonvanishing of the principal minors of F_105 on an explicit family of about 2^90.5 index sets ([below](#what-p105-reduces-to)).
6. **Orders with no located literature source.** 28 of the settled orders are in none of the sources found by the [prior-art search](../../order-154/prior-art/ganesh-claude-minors/README.md); they are listed at the top of [table.md](table.md). Examples: 221 = 13·17, 253 = 11·23, 255 = 3·5·17, 286 = 2·11·13, 442 = 2·13·17, 494 = 2·13·19. Six of the 28 are 154 and orders 14p that this repository's own seed also claims. Order 154 is settled here through two presentations, (11, 14) and (7, 22). The seven orders settled by selective lifting are 182, 253, 310, 377, 435, 442 and 465; at (17, 26), for instance, one chart has 104 singular sets and the other has none.

"Settled" means: P(N) follows from Chebotarev's theorem, the two cited criteria and the computations in this folder. Both criteria come from arXiv preprints. This report uses them as stated and does not re-prove them. "Not settled" and "obstructed" mean the criteria give no conclusion; they do not mean a determinant is zero.

Rows for order 154 and for orders 14p are not independent of the repository's seed (same GitHub owner, and I read the seed pages before this work). They do not replace the order-154 replay or the 14p audit.

## The two criteria

Let M be square-free and p a prime not dividing M. The prime ideals of ℤ[ζ_M] above p ("charts") correspond to the irreducible factors h of the cyclotomic polynomial Φ_M modulo p; all have degree f = ord_M(p), and reduction modulo the ideal is the field 𝔽_p[x]/(h) with ζ_M ↦ x. A subset S ⊆ ℤ_M is *singular in a chart* if det F_M[S, S] vanishes there.

- **Ordinary lifting** ([Caragea–Lee–Malikiosis–Pfander, Theorem 1.4](https://arxiv.org/abs/2505.24326v3)). If P(M) holds and no subset is singular in any chart above p, then P(pM). I write H(M, p) for the second hypothesis.
- **Selective lifting** ([Wang–Zhang, Lemma 4.1](https://arxiv.org/abs/2608.22119v1)). For K ⊆ ℤ_M × ℤ_p let c_a be the number of points of K above a ∈ ℤ_M, and call the nonempty sets {a : c_a ≥ t}, t = 1, …, p, the *level sets* of K. If some chart has no singular level set of K, the principal minor of F_pM on K is nonzero. Every chain of at most min(p, M) nested nonempty subsets of ℤ_M is the family of level sets of some K. So P(pM) follows if **no chain of at most min(p, M) nested subsets meets the singular sets of every chart**. I call such a chain a *covering chain*. This test needs no hypothesis on M.

The scan walks through the square-free N in increasing order. For each composite N it tries the presentations N = p·M with the largest prime first, and records N as settled by the first presentation that passes either test.

## Orders 105 and 110

| N | presentation (p, M) | field | charts | vanishing minors per chart | covering chain |
|---|---|---|---|---|---|
| 105 | (7, 15) | 7⁴ elements | 2 | 30, 30 (the same 30 sets) | one set: the punctured cross {3, 5, 6, 9, 10, 12} |
| 105 | (5, 21) | 5⁶ | 2 | 1428, 2100 | two nested sets, of sizes 7 and 13 |
| 105 | (3, 35) | 3¹² | 2 | 99 330, 99 330 (9 660 in both) | one set |
| 110 | (11, 10) | 11 | 4 | 80, 80, 40, 40 | four nested sets (as in Wang–Zhang, Prop 6.2) |
| 110 | (5, 22) | 5⁵ | 2 | 1848, 1848 | two nested sets, of sizes 9 and 13 |
| 110 | (2, 55) | 2²⁰ | 2 | not enumerated (2⁵⁵ subsets) | S ⊂ T with S = {0, 1, 2, 4, 33, 36, 42, 49} and T the complement of S + 5 |

For (3, 35) all 2³⁵ subsets were enumerated in both charts. For (2, 55) the chain was found by a targeted search and each of its two members was checked by a direct elimination: S is singular in one chart, and the complement of its translate S + 5, which contains S, is singular in the other. This is expected from Jacobi's complementary-minor identity, because the two charts are exchanged by ζ ↦ ζ⁻¹ (−1 is not a power of 2 modulo 55).

## Structural results

### Proposition 1 (punctured star)

Let M be square-free and let d₁, …, d_k ≥ 2 be pairwise coprime divisors of M. Let H_j be the subgroup of ℤ_M of order d_j and S = (H₁ ∪ ⋯ ∪ H_k) ∖ {0}. Then

  det F_M[S, S] = −c · ∏_j det A_j,  c = (k − 1)·∏_j d_j − Σ_j ∏_(i≠j) d_i,

where A_j = F_M[H_j ∖ {0}, H_j ∖ {0}] and |det A_j| = d_j^((d_j − 2)/2). In particular, a prime r ∤ M that divides c makes S singular in **every** chart above r, so the presentation (r, M) fails both criteria.

For two primes, c = pq − p − q = φ(pq) − 1.

*Proof.* Order the rows and columns by blocks H_j ∖ {0}. For h ∈ H_i and h′ ∈ H_j with i ≠ j, the product hh′ is a multiple of (M/d_i)(M/d_j) = M · M/(d_i d_j), so the entry is 1: every off-diagonal block is all ones. Inside a block, h′ ↦ ζ^(hh′) is a nontrivial character of H_j for h ≠ 0 (M is square-free, so M/d_j is prime to d_j), hence each row of A_j sums to −1, that is A_j 𝟙 = −𝟙. A_j is a Galois conjugate of the principal submatrix of F_(d_j) on ℤ_(d_j) ∖ {0}; Jacobi's identity with the complementary 1 × 1 minor gives |det A_j| = d_j^(d_j/2 − 1), so A_j is invertible and det A_j is a unit away from d_j.

Write 𝒜 for the block-diagonal matrix of the A_j and E for the |S| × k matrix whose columns are the indicator vectors of the blocks. The matrix is 𝒜 + E(𝟙𝟙ᵀ − I)Eᵀ. By the matrix determinant lemma its determinant is det 𝒜 · det(I + (𝟙𝟙ᵀ − I)Eᵀ𝒜⁻¹E). Since A_j⁻¹𝟙 = −𝟙, Eᵀ𝒜⁻¹E = −diag(d_j − 1), and the k × k determinant is det(diag(d_j) − 𝟙vᵀ) with v_j = d_j − 1, which equals ∏ d_j · (1 − Σ_j (d_j − 1)/d_j) = −c. ∎

Values of c: (2, 5) gives 3; (2, 7) gives 5; (3, 5) gives 7; (3, 7) gives 11; (5, 7) gives 23; (5, 11) gives 39; inside ℤ₃₀, (2, 15) gives 13, (3, 10) gives 17, (5, 6) gives 19 and (2, 3, 5) gives 29. So for instance (7, 15), (11, 21), (13, 30), (3, 55) and (13, 55) are obstructed.

### Proposition 2 (propagation)

Let M′ = Mm be square-free and p a prime not dividing M′. If S ⊆ ℤ_M is singular in a chart above p, then mS ⊆ ℤ_M′ is singular in a chart of ℤ[ζ_M′] above p. If a covering chain exists for (p, M), one exists for (p, M′). Hence a presentation that fails the ordinary test, or the selective test, still fails after the base is multiplied by m.

*Proof.* det F_M′[mS, mS] = det(ζ_M′^(m²ab))_(a,b∈S) = σ(det F_M[S, S]), where σ is the automorphism of ℚ(ζ_M) sending the primitive M-th root ζ_M′^m to its m-th power. So mS is singular in every chart of ℚ(ζ_M′) lying over the σ-image of a chart in which S is singular. Every chart of ℚ(ζ_M′) lies over some chart 𝔮 of ℚ(ζ_M); a covering chain {S_i} contains a member singular in σ⁻¹𝔮, so {mS_i} contains a member singular in that chart, and it is again a nested chain of the same length. ∎

This is the mechanism of [Wang–Zhang, Prop 7.1](https://arxiv.org/abs/2608.22119v1) applied in finite characteristic. In [table.md](table.md) it is used to decide presentations whose base is too large to enumerate, for example (7, 30) from (7, 15).

### What P(105) reduces to

Write ℤ₁₀₅ = ℤ₁₅ × ℤ₇ and C = (H₃ ∪ H₅) ∖ {0} = {3, 5, 6, 9, 10, 12}. The singular sets above 7 are the same in both charts: the translates of C and of its complement (checked by `structure_checks.py`). By selective lifting, the minor on K is nonzero unless some level set of K is one of those 30 sets. Translating K by an element of ℤ₁₅ multiplies the minor by a root of unity, and replacing K by its complement in ℤ₁₀₅ preserves vanishing and replaces each level set by a complement. Therefore:

> Assuming Lemma 4.1 of Wang–Zhang, P(105) holds if and only if det F₁₀₅[K, K] ≠ 0 for every K ⊆ ℤ₁₅ × ℤ₇ such that, for some t ∈ {1, …, 7}, the set {a : c_a ≥ t} equals C.

There are Σ_t (Σ_(c≥t) C(7, c))⁶ (Σ_(c<t) C(7, c))⁹ ≈ 1.8 · 10²⁷ ≈ 2^90.5 such K, out of 2¹⁰⁵. The family can be cut further, since a K in it must also be obstructed at the presentations (5, 21) and (3, 35), but I have no closed description of that intersection. Enumerating the family is out of reach; it needs an argument that goes beyond the first-order term.

The smallest slice is checked exactly. For b ∈ (ℤ₇)⁶ let K_b = {(a, b_a) : a ∈ C}, the single-layer lifts of the cross. `lifts105_exact.py` shows that all 7⁶ = 117 649 minors are nonzero, by reducing modulo the prime ℓ = 997 500 421 ≡ 1 (mod 105). As a sanity check only, `sample105.py` draws random members of the whole family and tests them the same way; 2000 draws gave no vanishing minor.

## Tables

[table.md](table.md) is generated by `make_table.py` from [scan500.json](scan500.json). It lists the settled orders with no located literature source, the orders not settled with the outcome of every presentation, the pairs (M, p) that fail the ordinary test with their counts of vanishing minors, and all composite square-free orders up to 500. The raw log of the scan is [scan500.txt](scan500.txt). `scan500.json` holds the scan, the `pairs` run and the two runs of [special-runs.md](special-runs.md), whose outcomes for (35, 3) and (55, 2) were entered by hand.

## Verification

- **Two arithmetic back ends.** `pmscan.c` builds either with discrete-log and Zech tables or with polynomial arithmetic in 𝔽_p[x]/(h). `backend_check.py` finds identical lists of vanishing minors from the two builds on (7, 2), (10, 3), (10, 11), (11, 23), (14, 29), (15, 7), (21, 5) and (22, 5).
- **Independent implementation.** `verify_small.py` recomputes Φ_M from the Möbius product, checks the factors, and finds the vanishing minors by one Gaussian elimination per subset with its own field arithmetic. Its zero sets are identical to the scanner's on (6, 5), (7, 2), (10, 3), (10, 11), (10, 31), (11, 3), (14, 11), (14, 13), (14, 29) and (15, 7).
- **Published data reproduced.** Base 10 above 11: 20 + 20 singular sets at ranks 3 and 7 in two charts, 40 + 40 at ranks 4 and 6 in the other two, and a covering chain ([Wang–Zhang, Table 3 and Prop 6.2](https://arxiv.org/abs/2608.22119v1)). Base 10 above 31: 40 at rank 5 in two charts, none in the other two, no covering chain (their Prop 6.1). Base 10 above 3: 10, 10, 10 at ranks 4, 5, 6 (their Table 2). Base 10 above 7 and base 11 above 13: none ([Gu–Zhou–Wang](https://arxiv.org/abs/2608.17746v1)). Base 11 above 3 and 23: singular sets exist, as the published norm support requires.
- **Repository seed reproduced.** Base 14: no covering chain above 13 and 113, a covering chain above 29 and 43, no singular set above 11.
- **Consistency with the literature.** Every order ≤ 500 that a located source proves is also settled by the scan (`make_table.py` asserts this).
- **Structural checks.** `structure_checks.py` checks Proposition 1 in floating point for fifteen cases, checks that punctured stars are singular in all charts for five (M, r), and checks the description of the singular sets of F_15 above 7.

## Reproduce

```sh
cc -O2 -o BIN/pmscan pmscan.c
cc -O2 -DPOLY -o BIN/pmscan_poly pmscan.c
python3 frontier.py scan 500 BIN scan500.json 4e10 > scan500.txt
python3 frontier.py pairs 22 31 BIN scan500.json 4e10
python3 make_table.py scan500.json > table.md
python3 verify_small.py BIN 6:5 7:2 10:3 10:11 10:31 11:3 14:11 14:13 14:29 15:7
python3 structure_checks.py BIN
```

`BIN` is any directory for the two binaries. Environment used: Apple clang 17.0.0, Python 3.9.6 (standard library only), macOS on Apple silicon. The scan took 30 minutes on one core, 16 of them for the base 26 over the field with 19¹² elements; `pairs` took 7 minutes. The other scripts (`backend_check.py`, `lifts105_exact.py`, `sample105.py . 2000`) take under a minute each. The two special runs are recorded in [special-runs.md](special-runs.md): the full enumeration for (35, 3) and the chain for (55, 2).

## Limits

- The two criteria are taken from preprints and are not re-proved here. I have not found a peer-reviewed version of either.
- The scan stops at the first presentation that works, so [table.md](table.md) does not list every working presentation of a settled order.
- A base M > 22 gets the selective test only when the ordinary test fails and a full enumeration fits the budget (M ≤ 30). Presentations beyond the budget are decided by propagation when possible and are otherwise marked undecided; an unsettled order with an undecided presentation may still be settled by these criteria.
- "No located literature source" depends on the search documented in the prior-art report and proves no priority.
- None of this is an independent replication of the order-154 seed or the 14p seed.

## Next

The first-order criteria are exhausted at 105 and 110. The next step is a criterion that sees past the first-order term. I have a candidate: writing ζ_p^(bb′) = Σ_m Γ_m [b]^m [b′]^m with Teichmüller representatives and Gauss sums Γ_m gives an exact expansion of every principal minor of F_pM in powers of a uniformizer, whose leading term is the known first-order product and whose next term is an explicit quadratic form in the fibre sums of K. A prototype matches direct computation on 264 configurations (leading term) and on all 38 obstructed ones among them (second term), and at order 105 its second term is nonzero in at least one chart for every non-constant single-layer lift of the punctured cross. It is not written up or reviewed, and it cannot be the whole answer: on a product set K = S × B the minor is det F_M[S]^|B| · det F_p[B]^|S|, whose valuation grows with |B|, so no fixed order of the expansion decides every K.
