# Lean: finite legal paths have total probability one

[`ForestPathNormalization.lean`](../../../formal/ForestPathNormalization.lean) extends the merged path-probability module. It proves that, from any starting occupied subset, the probabilities of all legal histories of a finite gate schedule sum to one. A second theorem transfers this identity to the squared Gaussian-rational path amplitudes. **This is normalization of paths, not a Cauchy–Binet identity or a proof of MDO.**

## Exact checked statement

Use the existing finite vertex enumeration, `Edge vertices`, Boolean-function `Subset V`, `stepProb`, and `pathProbability`. `successors e s` is `[s]` when the two endpoints have equal occupancy, and `[s, swapSubset e s]` otherwise. Thus a collision does not count a stay and an invisible swap twice.

`histories edges start n` builds the finite list of histories by appending one legal successor at each time. Each history is represented by a function on natural-number times; its values after time `n` are fixed to `start`. This tail is only a canonical representation, not an extra path constraint.

The main declarations are in namespace `MathCommons.ForestPathNormalization`:

| Declaration | Kernel-checked conclusion |
|---|---|
| `successorMass` | The existing `stepProb` values over the collision-deduplicated successors sum to one. |
| `successorLegal`, `everyLegalSuccessor`, `successorsNodup` | The successor list contains exactly the legal next subsets, without duplicates. |
| `historyLegal` | Every enumerated history begins at `start` and has legal steps through the specified finite length. |
| `everyLegalPrefix` | Every legal history prefix has a representative in the enumeration. No condition is imposed on the input history's later times. |
| `historyFuture`, `historyDeterminedByPrefix`, `historiesNodup` | The canonical tail makes finite-prefix equality unambiguous, and no history is counted twice. |
| `pathMassOne` | The sum of the **existing** `pathProbability edges bs history n` over this enumeration equals one. |
| `stepProbabilityNonnegative`, `pathProbabilityNonnegative` | A schedule of unitary Gaussian-rational blocks gives nonnegative step and path weights. |
| `pathSquaredMagnitudeMassOne` | Under the same unitary hypotheses, the sum of `normSq (pathAmplitude ...)` over all these histories equals one. |

Five elementary finite-sum lemmas and two prefix/extension lemmas are also checked and audited. There are twenty printed transitive axiom reports in total. The allowed dependencies are only Lean's standard `propext`, `Classical.choice`, and `Quot.sound`; no new axiom, `sorry`, or `native_decide` is used.

## Why normalization holds

For equal endpoint occupancies, the only successor has weight one. Otherwise the two distinct successors have weights `1−q` and `q`, where `q=normSq b`. This proves the single-step mass identity. The prefix lemma shows that appending a new subset leaves the earlier probability product unchanged. Summing the appended step over its successors multiplies each old history's mass by one. Induction on the schedule length proves `pathMassOne`, including length zero.

Unitarity is needed for nonnegativity and for the squared-amplitude identification. The algebraic mass identity itself remains true for arbitrary rational `q`. The result needs **no forest hypothesis**: a gate schedule with cycles also defines a normalized path distribution. Forest uniqueness becomes relevant when comparing paths with fixed endpoints; interference and the matrix-product identity remain separate issues.

The existing edge type requires a supplied exhaustive, duplicate-free enumeration and distinct endpoints. The zero-length theorem uses no gate, but the current total schedule interface still takes an edge family. This does not introduce an encoding of every zero- or one-vertex graph.

## Dependencies and reproduction

Upstream status was checked before the claim:

* [PR 22](https://github.com/g4nesh/math-commons/pull/22), `ForestGateMagnitude`, merged at `cd6ad6e0ff50576015f8b72401a4ece9209ce5fa`.
* [PR 29](https://github.com/g4nesh/math-commons/pull/29), Raghav's `ForestPathProbability`, merged at `cdda55919a2cb514f1bf95791ccc24cf5237a9f8`.

The new file imports `ForestPathProbability` without modifying it. The dependency chain through `ForestCut`, `ForestTrajectory`, `ForestGateSupport`, and `ForestGateMagnitude` was rebuilt using Lean **4.24.0**, with `-DwarningAsError=true`, before compiling this module. The repository checker integrates the new file after these dependencies; from the repository root run:

```sh
python3 formal/check.py
```

The checker supports an explicit path to a standalone pinned Lean binary as its first argument. It rebuilds imports in a fresh temporary directory and rejects unexpected or missing axiom reports. Source compilation and axiom inspection establish the formal statements; independent-owner mathematical review remains a separate ledger step.

## Remaining bridge

The file does not identify the sum of signed path amplitudes with the compound matrix of the gate product, translate the list determinant into a standard complex matrix determinant, or derive the final permutation-mixture theorem. In particular, summing squared path amplitudes is different from squaring the sum of amplitudes when several paths interfere. Cauchy–Binet and the standard-complex bridge remain explicit next tasks.
