# Marcus–de Oliveira · research overview

**Objective `marcus-de-oliveira` · open.** For all complex normal matrices $A,B$ of size $n$ with eigenvalues $a_i,b_i$, does $\det(A+B)$ belong to the closed planar convex hull of $\prod_i(a_i+b_{\sigma(i)})$ over all permutations $\sigma$? This repository has neither a general proof nor a certified outside pair. [Exact target and boundary cases](README.md#exact-target)

| Choose a question | Current finding and limit |
| --- | --- |
| [MDO-RED · Which reductions hold?](context/reductions.md) | Independently reviewed reductions give a rational witness if any counterexample exists. A claimed dimension-four proof still lacks a complete audit here. |
| [MDO-FOREST · When do gate weights work?](context/forest.md) | A written proof gives weights when each forest edge is used once. A triangle defeats those weights. Lean checks a cut-history lemma. |
| [MDO-MINORS · Which spectra are covered?](context/minors.md) | Theorem D gives a written proof for three distinct eigenvalues on one side with one simple. Independent review and the unrestricted three-value case remain open. |
| [MDO-TEST · What was searched and certified?](context/witnesses.md) | All 924 exact stress cases in dimensions 4–7 are inside. An earlier 200-sample floating-point search has a separate review. |
| [MDO-WORK · What are others pursuing?](context/work.md) | PRs track corrected cycle weights, a Lean subset-trajectory bridge and a dimension-four Hermitian route. Each has a distinct deliverable. |

Based on [main at `c183060`](https://github.com/g4nesh/math-commons/tree/c183060438e58cce0c29e922a7d48fd2443ef3e7). For a synthesis, consult the [proof obligations](proof-obligations.md). For ownership or a new result, check the [live board](https://math-commons-dun.vercel.app/#work) and [current PRs](https://github.com/g4nesh/math-commons/pulls) before acting.
