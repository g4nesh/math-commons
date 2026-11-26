# Marcus–de Oliveira · research overview

**Objective `marcus-de-oliveira` · open.** For all complex normal matrices $A,B$ of size $n$ with eigenvalues $a_i,b_i$, does $\det(A+B)$ belong to the closed planar convex hull of $\prod_i(a_i+b_{\sigma(i)})$ over all permutations $\sigma$? This repository has neither a general proof nor a certified outside pair. [Exact target and boundary cases](README.md#exact-target)

| Choose a question | Current finding and limit |
| --- | --- |
| [MDO-RED · Which reductions hold?](context/reductions.md) | Independently reviewed reductions give a rational witness if any counterexample exists. The dimension-four audit reproduced the finite classification. The 90 contact inequalities still need case-level verdicts. |
| [MDO-FOREST · When do gate weights work?](context/forest.md) | Pointwise certificates glue across one shared coordinate in any dimension. The fourth-eigenvalue path has an exact algebraic feasibility boundary for uniform weights. Lean now checks the whole ordered product of the local probabilities. |
| [MDO-MINORS · Which spectra are covered?](context/minors.md) | Theorem D gives a written proof for three distinct eigenvalues on one side with one simple. Independent review and the unrestricted three-value case remain open. |
| [MDO-TEST · What was searched and certified?](context/witnesses.md) | All 924 exact stress cases in dimensions 4–7 are inside. An earlier 200-sample floating-point search has a separate review. |
| [MDO-WORK · What needs review or a next step?](context/work.md) | Review the spectral certificates and the Lean gate proof. Claude’s merged [interval certificates](star-halfline/README.md) cover four-dimensional Hermitian spectra with cross-ratio in [1/6, 5/6]. The remaining parameter range and independent review are open. |

For the scope of each result, see the [proof obligations](proof-obligations.md). For ownership or a new result, check the [live board](https://math-commons-dun.vercel.app/#work) and [current PRs](https://github.com/g4nesh/math-commons/pulls) before acting.
