# Lean: no interference in a forest endpoint fiber

[`ForestEndpointMass.lean`](../../../formal/ForestEndpointMass.lean) proves that, for a fixed starting and ending occupied subset in a forest schedule, the squared norm of the sum of signed path amplitudes equals the sum of their path probabilities. It includes an empty endpoint fiber. **The sum is an explicitly defined finite path sum; its identification with a minor of a matrix product remains open.**

## Checked statement

The file imports `ForestPathNormalization` and reuses its canonical, duplicate-free list `histories edges start n`. It defines `endpointHistories` by filtering this list for the ending subset `finish`. The endpoint test checks all vertices in the existing exhaustive finite enumeration, so it is equivalent to equality of the Boolean subset functions.

Writing $H_{S,T}$ for this finite list, the main theorem is

$$
\left|\sum_{h\in H_{S,T}} \operatorname{pathAmplitude}(h)\right|^2
=\sum_{h\in H_{S,T}} \operatorname{pathProbability}(h).
$$

Here amplitudes and their squared norms are the existing Gaussian-rational definitions. The hypotheses are `EveryEdgeIsBridge edges n` and a `UnitaryBlock` at each time below `n`. Both endpoints and the schedule length are arbitrary within the existing finite-vertex interface.

All declarations below are in namespace `MathCommons.ForestEndpointMass`:

| Declaration | Kernel-checked conclusion |
|---|---|
| `sumWithZero` | A finite rational sum of zero terms is zero. |
| `normSumUnique` | For a duplicate-free list of Gaussian values with at most one nonzero value, the squared norm of their sum equals the sum of squared norms. |
| `endpointMember` | Membership in the filtered list is exactly canonical-history membership and equality to the requested endpoint. |
| `endpointHistoriesNodup` | Filtering does not duplicate histories. |
| `endpointNonzeroUnique` | Two histories in the same forest endpoint fiber cannot both have nonzero amplitude unless they are equal. This support statement requires no unitarity. |
| `endpointMagnitude` | The displayed endpoint amplitude/probability identity. |
| `emptyEndpointFiber` | If the endpoint list is empty, its amplitude sum is Gaussian zero and its probability mass is rational zero. This assertion needs neither the forest nor unitary hypotheses. |

## Proof and attribution

The forest support theorem `ForestGateSupport.nonzeroAmplitudesUnique` shows that two nonzero paths with the same endpoints agree through time `n`. The canonical tail supplied by `ForestPathNormalization.historyDeterminedByPrefix` turns this prefix equality into equality of the represented histories. Their enumeration is duplicate-free, so there is at most one nonzero summand, rather than potentially repeated copies of one path.

The elementary `normSumUnique` lemma then removes all cross terms without expanding a quadratic sum. Finally, `ForestPathProbability.pathMagnitude` replaces each squared amplitude by its existing probability. This last result is Raghav's merged [PR 29](https://github.com/g4nesh/math-commons/pull/29); this file does not modify or re-prove that module. The canonical enumeration and its normalization come from [PR 36](https://github.com/g4nesh/math-commons/pull/36), pinned at [`d1e4672380b38dc10408c332875a9eef7e65c551`](https://github.com/g4nesh/math-commons/blob/d1e4672380b38dc10408c332875a9eef7e65c551/formal/ForestPathNormalization.lean).

## Verification and scope

The module was compiled with Lean **4.24.0** and `-DwarningAsError=true`, using rebuilt dependency modules. All seven printed transitive axiom reports contain only the permitted standard axioms `propext`, `Classical.choice`, and `Quot.sound`. There is no `sorry`, `native_decide`, or new axiom. The repository checker registers this module after its imports:

```sh
python3 formal/check.py
```

The checker's optional first argument is a path to the pinned standalone Lean binary. Source checking and the dependency audit establish these formal statements; independent-owner mathematical review remains a separate ledger action.

This result resolves the finite-path no-interference step. It does not prove iterated Cauchy–Binet, identify the list determinant with the usual complex matrix determinant, or identify `endpointAmplitude` with the corresponding minor of a gate product. A common permutation-mixture theorem and the resulting MDO conclusion are also outside this file. The scalar domain remains Gaussian rationals, and the inherited total edge-family interface still requires a supplied exhaustive vertex enumeration and distinct edge endpoints even when `n = 0`.
