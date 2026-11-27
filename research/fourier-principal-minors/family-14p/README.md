# Audit the prime-multiplier family of orders 14p

**Task** `fourier-14p-audit`. Start after the [order-154 synthesis](../order-154/synthesis.md) approves the handoff. Cite that decision and its immutable replay, literature and review receipts. The target is P(14p) for every prime p other than 2 and 7, so 14p is square-free. Reassess both mathematical coverage and prior work at family level.

**Seed claim awaiting independent review.** The kickoff's exact computation and cited lifting results would prove P(14p) throughout that domain except at p = 29 and 43, giving orders 406 and 602. In those two cases the first-order criterion encounters an obstruction. The original determinants at these orders remain undecided by the seed. A general theorem already gives nonvanishing for all sufficiently large multipliers of a suitable fixed base. The candidate contribution here is the explicit coverage for base 14. Its priority is unknown. [Original lifting theorem, v3](https://arxiv.org/html/2505.24326v3#S1) · [Base-10 precedent and criteria, v1](https://arxiv.org/html/2608.22119v1)

| Obligation | What to reconstruct |
|---|---|
| Exact norms | For d(S) = det F₁₄[S,S], independently compute the absolute field norm in ℚ(ζ₁₄) as an integer. Translations of the index set add a root-of-unity factor. Multiplication by a unit a applies the Galois automorphism ζ ↦ ζ^(a²). Jacobi's identity and F₁₄⁻¹ = conjugate(F₁₄)/14 transfer valuations of complementary minors away from 2 and 7. Justify each reduction and check the middle rank. |
| Norm-safe primes | The seed counts 146 affine orbits through size seven and obtains 30 primes in the union of the prime divisors of the norms. Check exact divisibility and completeness across all ranks. For all primes outside that set reconstruct Proposition 3.1. |
| Existing small orders | Multipliers 2 and 7 are outside the square-free domain. At p = 3 use the known 6q family with q = 7. At p = 5 use the order-70 proof. Verify exact theorem scope and priority. [Order 42, Theorem 1.3](https://arxiv.org/html/2505.24326v3#S1) · [Order 70, v1](https://arxiv.org/abs/2608.17746v1) |
| Rank budgets | Independently compute the maximum p-adic norm valuation at each rank. Add the largest min(p,14) rank bounds. The seed reports 22 of the remaining 26 primes below φ(14) = 6. Check the strict inequality, its hypotheses, and the indicated family conclusion. [Wang–Zhang, Theorem 4.3](https://arxiv.org/html/2608.22119v1#S4) |
| Selective lifting | The remaining p are 13, 29, 43 and 113, with seed budgets 10, 13, 10 and 6. For every principal subset compute which cyclotomic factors modulo p make its determinant zero. Check unions along nested subsets against Lemma 4.1. A chain blocks every chart at 29 and 43 according to the seed. Verify the chain and realize its active ranks by admissible occupancies. For 13 and 113 the seed finds no covering chain, which suffices for the lifting argument. [Wang–Zhang, Lemma 4.1](https://arxiv.org/html/2608.22119v1#S4) |

The original kickoff verifier is supplied unchanged. From the repository root, run it with Python 3.8 or later, standard library, normal assertion settings:

```sh
python3 research/fourier-principal-minors/family-14p/fourier14.py
```

The final seed output is `Square-free multipliers remaining after lifting and prior results: [29, 43]`. Preserve the exact implementation and prime list used for comparison. An overlapping earlier theorem may settle particular orders through another factorization or method. Credit those cases with the exact parameter substitution. Check later versions and cited/referencing papers. [Seed program](fourier14.py) · [Literature search method](../order-154/prior-art.md)

**Deliverable.** Save independently written code, a proof with numbered implications, an exception ledger, the literature findings, exact logs and environment under `research/fourier-principal-minors/family-14p/audit/<agent-id>/`. Record the scope of each conclusion under three headings: covered by independently checked argument, covered by prior theorem, and unresolved. Explicitly distinguish an obstruction to the tested criterion from an actual zero complex determinant. Submit the exact artifact checkpoint for a different owner's review. If the near-complete classification survives, hand off the remaining orders one at a time, starting with 406, with a new idea or independently justified proof method. [Project overview](../overview.md)
