# MDO-ZHENG4 · Does the claimed proof for $n=4$ hold up?

**Status.** The dimension-four proof remains under audit. The supplement has been located and pinned. The finite classification has an exact replay of 214 contact-set orbits, of which 124 reduce to nonnegative squared minors and 90 require separate inequalities. [Finite classification and limits](reduction/finite-classification.md) · [pinned inputs](METHOD.md#pinned-inputs)

**Current coverage.** The manifest lists 208 proof documents. This checkpoint contains 38 first-pass document reports with replay scripts. A separate numerical search found 246 configurations across 82 of the 90 contact sets. The eight remaining searches produced no examples. The reports and sampled configurations do not establish the ninety inequalities. [Document inventory](manifest.json) · [first-pass reports](docs/) · [numerical coverage and tolerances](contact-configs.md)

**Next question.** Which case can be settled once its shared dependencies and remaining proof steps receive an independent second pass? For a specific case, look up its document in the manifest, then open the matching report and replay, if present. Follow the stated evidence levels and verdict rules when recording a result. [Verdict rules](METHOD.md#verdict-per-case)

## What is being audited

D. Zheng, *The Marcus–de Oliveira conjecture for n = 4* (24 October 2026, [doi:10.13140/RG.2.2.36790.25923](https://doi.org/10.13140/RG.2.2.36790.25923)), and its public supplement at commit `53ec16567d2307266ce52a8e9678410b2457cc27`. Exact identities of both are in the [specification](METHOD.md#pinned-inputs).

The paper's argument has five parts:

1. **Unitary to orthogonal (Theorem 1.2).** The squared moduli of all minors of a $4\times4$ unitary matrix are one convex combination of those of two real orthogonal matrices.
2. **Contact reduction (Theorem 3.6).** A counterexample at a fixed orthogonal matrix can be moved, along a curve in the six-dimensional cone of eight-point matching invariants, to one with at least six permutation products on the supporting line.
3. **Finite classification (Section 4).** The 82 facets of the minor polytope and the action of a group of order 1152 sort the six-element contact sets into 214 orbits. In 124 of them the inequality follows from nonnegativity of squared minors.
4. **Ninety contact inequalities (Proposition 4.3).** The paper gives its argument for one case and directs readers to the supplement for the other 89, which use shared auxiliary results.
5. **Special spectra (Appendices B and C).** Coincident eigenvalues, and one spectrum on a circle or line.

## Where to look

| Question | Evidence |
| --- | --- |
| How are contact sets obtained and classified? | [Contact reduction](reduction/contact-reduction.md) · [finite classification and exact checks](reduction/finite-classification.md) |
| What supports the reduction from unitary to orthogonal matrices? | [First review of Theorem 1.2](foundations/theorem-1-2.md) · [second independent review](foundations/theorem-1-2-second.md) |
| What has been checked for coincident eigenvalues? | [Appendix B review and limitations](foundations/appendix-b.md) |
| How can I reproduce a computation? | [Library, dependencies and 77 checks](lib/README.md) · [per-document replay scripts](replay/) |
| What was assigned to other workers? | [Worker instructions and slice assignments](CLOUD-WORKER.md) · [document inventory](manifest.json) |
| Which arguments depend on dimension four? | [Dimension-five analysis](generalization-n5.md) |

[Research overview](../overview.md)
