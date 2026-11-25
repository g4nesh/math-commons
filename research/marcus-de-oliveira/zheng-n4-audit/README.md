# MDO-ZHENG4 · Does the claimed proof for $n=4$ hold up?

**Status: audit in progress. Nothing below is a verdict on the proof yet.**

**Finding.** None yet. The inputs are pinned and the audit plan is fixed in the [specification](METHOD.md).

**Limit.** No step of the claimed proof has been confirmed or refuted here so far.

**Next question.** Do the paper's Sections 2–4 reduce the conjecture for $n=4$ to the ninety contact inequalities of its Appendix E, and does each of the ninety supplementary proofs hold?

**Evidence.** [Specification, evidence levels and verdict rules](METHOD.md). Task `mdo-zheng-n4-audit` in the live ledger.

## What is being audited

D. Zheng, *The Marcus–de Oliveira conjecture for n = 4* (24 October 2026, [doi:10.13140/RG.2.2.36790.25923](https://doi.org/10.13140/RG.2.2.36790.25923)), and its public supplement at commit `53ec16567d2307266ce52a8e9678410b2457cc27`. Exact identities of both are in the [specification](METHOD.md#pinned-inputs).

The paper's argument has five parts:

1. **Unitary to orthogonal (Theorem 1.2).** The squared moduli of all minors of a $4\times4$ unitary matrix are one convex combination of those of two real orthogonal matrices.
2. **Contact reduction (Theorem 3.6).** A counterexample at a fixed orthogonal matrix can be moved, along a curve in the six-dimensional cone of eight-point matching invariants, to one with at least six permutation products on the supporting line.
3. **Finite classification (Section 4).** The 82 facets of the minor polytope and the action of a group of order 1152 sort the six-element contact sets into 214 orbits. In 124 of them the inequality follows from nonnegativity of squared minors.
4. **Ninety contact inequalities (Proposition 4.3).** One is proved in the paper. The other 89 are proved in the supplement, with shared auxiliary results.
5. **Special spectra (Appendices B and C).** Coincident eigenvalues, and one spectrum on a circle or line.

## Layout

| Path | Content |
| --- | --- |
| [METHOD.md](METHOD.md) | Specification, evidence levels, verdict rules, deviations |
| `lib/` | Own exact library for the objects of the paper, with 77 tests |
| `reduction/` | Re-derivation of parts 2 and 3 and their exact replays |
| `foundations/` | Audits of parts 1 and 5 |
| [manifest.json](manifest.json) | The 208 proof documents, their roles, and who audits each batch |
| `docs/`, `replay/` | First-pass report and own-code replay for each proof document |
| [CLOUD-WORKER.md](CLOUD-WORKER.md), `cloud/` | Instructions for workers auditing two slices of documents, and their results |
| [generalization-n5.md](generalization-n5.md) | Which ingredients are specific to dimension four |

First-pass document reports are not verdicts on cases. Case verdicts are assembled after a second pass.

[Research overview](../overview.md)
