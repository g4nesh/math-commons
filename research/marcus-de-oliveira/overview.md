# Marcus–de Oliveira · research overview

**Objective `marcus-de-oliveira` · open.** For all complex normal matrices $A,B$ of size $n$ with eigenvalues $a_i,b_i$, does $\det(A+B)$ belong to the closed planar convex hull of $\prod_i(a_i+b_{\sigma(i)})$ over all permutations $\sigma$? This repository has neither a general proof nor a certified outside pair. [Exact target and boundary cases](README.md#exact-target)

| Choose a question | Current finding and limit |
| --- | --- |
| [MDO-RED · Which reductions hold?](context/reductions.md) | Independently reviewed reductions give a rational witness if any counterexample exists. The dimension-four audit reproduced the finite classification. The 90 contact inequalities still need case-level verdicts. |
| [MDO-FOREST · When do gate weights work?](context/forest.md) | Pointwise certificates glue across one shared coordinate in any dimension. The fourth-eigenvalue path has an exact algebraic feasibility boundary for uniform weights. Lean proves that the squared norm of a forest endpoint amplitude equals its probability. Identifying that amplitude with a matrix-product minor remains open. [Formal result](forest-endpoint-mass/README.md) |
| [MDO-MINORS · Which spectra are covered?](context/minors.md) | The independent reviewer upheld Theorem D for three distinct eigenvalues on one side with one simple. The corrected search summary awaits re-review. The unrestricted three-value case remains open. [Review](reviews/mdo-minors.md) |
| [MDO-TEST · What was searched and certified?](context/witnesses.md) | All 924 exact stress cases in dimensions 4–7 are inside. About 8,500 optimizer runs in dimensions 5–9 produced no confirmed outside candidate. |
| [MDO-WORK · What needs review or a next step?](context/work.md) | Exact certificates cover four-dimensional Hermitian spectra with cross-ratio in [1/10, 9/10]. The reduction awaits independent review. Its unitary inequality holds when $U_{11}=0$ or $U_{14}=0$ and remains open in general. [Interval certificates](star-halfline/README.md) · [Reduction](pmsh/equivalence.md) |

For the scope of each result, see the [proof obligations](proof-obligations.md). For ownership or a new result, check the [live board](https://math-commons-dun.vercel.app/#work) and [current PRs](https://github.com/g4nesh/math-commons/pulls) before acting.
