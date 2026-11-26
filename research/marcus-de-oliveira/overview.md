# Marcus–de Oliveira · research overview

**Objective `marcus-de-oliveira` · open.** For all complex normal matrices $A,B$ of size $n$ with eigenvalues $a_i,b_i$, does $\det(A+B)$ belong to the closed planar convex hull of $\prod_i(a_i+b_{\sigma(i)})$ over all permutations $\sigma$? This repository has neither a general proof nor a certified outside pair. [Exact target and boundary cases](README.md#exact-target)

| Choose a question | Current finding and limit |
| --- | --- |
| [MDO-RED · Which reductions hold?](context/reductions.md) | Independently reviewed reductions give a rational witness if any counterexample exists. The dimension-four audit reproduced the finite classification. The 90 contact inequalities still need case-level verdicts. |
| [MDO-FOREST · When do gate weights work?](context/forest.md) | Lean identifies actual order-two forest product minors with endpoint amplitudes and proves their squared magnitudes form normalized probabilities. Higher-order work now focuses on the one-gate predecessor recurrence. [Lean proof](forest-pair-history/README.md) |
| [MDO-MINORS · Which spectra are covered?](context/minors.md) | The independent reviewer upheld Theorem D for three distinct eigenvalues on one side with one simple. The corrected search summary awaits re-review. The unrestricted three-value case remains open. [Review](reviews/mdo-minors.md) |
| [MDO-TEST · What was searched and certified?](context/witnesses.md) | All 924 exact stress cases in dimensions 4–7 are inside. About 8,500 optimizer runs in dimensions 5–9 produced no confirmed outside candidate. |
| [MDO-WORK · What needs review or a next step?](context/work.md) | The archived four-cycle permits unbounded motion of three eigenvalues on two slices. Ordinary Fourier bases through order eight admit common weights for both arbitrary spectra. Order nine has an explicit obstruction to the stronger common-minor method. [Four-cycle weights](full-four-cycle/README.md) · [Fourier results](fourier-composite-laws/README.md) |

For the scope of each result, see the [proof obligations](proof-obligations.md). For ownership or a new result, check the [live board](https://math-commons-dun.vercel.app/#work) and [current PRs](https://github.com/g4nesh/math-commons/pulls) before acting.
