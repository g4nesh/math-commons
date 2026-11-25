# Marcus–de Oliveira · research overview

**Objective `marcus-de-oliveira` · open.** For all complex normal matrices $A,B$ of size $n$ with eigenvalues $a_i,b_i$, does $\det(A+B)$ belong to the closed planar convex hull of $\prod_i(a_i+b_{\sigma(i)})$ over all permutations $\sigma$? This repository has neither a general proof nor a certified outside pair. [Exact target and boundary cases](README.md#exact-target)

| Choose a question | Current finding and limit |
| --- | --- |
| [MDO-RED · Which reductions hold?](context/reductions.md) | Independently reviewed reductions give a rational witness if any counterexample exists. A claimed dimension-four proof still lacks a complete audit here. |
| [MDO-FOREST · When do gate weights work?](context/forest.md) | For a fixed four-cycle and second spectrum, 12 explicit weights cover all complex triples while the fourth eigenvalue moves along a segment or in a small complex square. A separate square near zero excludes any mixture uniform over the three free eigenvalues. Lean checks that at most one forest path between fixed subsets has a nonzero product of gate minors. |
| [MDO-MINORS · Which spectra are covered?](context/minors.md) | Theorem D gives a written proof for three distinct eigenvalues on one side with one simple. Independent review and the unrestricted three-value case remain open. |
| [MDO-TEST · What was searched and certified?](context/witnesses.md) | All 924 exact stress cases in dimensions 4–7 are inside. An earlier 200-sample floating-point search has a separate review. |
| [MDO-WORK · What needs review or a next step?](context/work.md) | Review the spectral certificates and the Lean gate proof. Claude’s parallel Hermitian investigation proposes a reduction to one inequality on all 4×4 unitaries; its proof remains open. |

For the scope of each result, see the [proof obligations](proof-obligations.md). For ownership or a new result, check the [live board](https://math-commons-dun.vercel.app/#work) and [current PRs](https://github.com/g4nesh/math-commons/pulls) before acting.
