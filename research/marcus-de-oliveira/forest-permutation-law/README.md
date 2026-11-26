# One common permutation law for all forest product minors

`MathCommons.ForestPermutationLaw.commonSquaredMinorLaw` constructs **one finite list of genuine permutations with nonnegative weights summing to one** that simultaneously represents every squared minor of an actual unitary forest gate product. The list and its weights are independent of both selected subsets. This is the common-law step after the actual arbitrary-order minor identity in [PR73](https://github.com/g4nesh/math-commons/pull/73).

[Lean source](../../../formal/ForestPermutationLaw.lean), [declarations](declarations.json), [audit manifest](audit.json). Scalars remain Gaussian rationals, using the repository's actual matrix product and ordered Laplace determinant. This is not a formal arbitrary-complex matrix-library bridge or a proof of MDO. The author's kernel check is not independent-owner mathematical review.

## Exact common law

Let $G_k$ have the unitary two-coordinate block $\begin{pmatrix}a_k&b_k\\c_k&d_k\end{pmatrix}$ and put $q_k=|b_k|^2$. Let $\tau_k$ be the actual vertex transposition exchanging the two endpoints. Begin with the identity permutation of weight one. At step $k$, replace each current weighted permutation $(\sigma,w)$ by

$$
(\sigma,w(1-q_k)),\qquad (\tau_k\circ\sigma,wq_k).
$$

This is `coinLaw edges bs n`, a finite list depending only on the gate schedule, the off-diagonal gate entries, and the product length. Zero weights and repeated permutation images are retained. They cause no problem for a convex combination: repeated occurrences are separate nonnegative contributions and can be aggregated if desired. The theorem does not claim a duplicate-free enumeration of the symmetric group.

`VertexPermutation` stores `forward`, `inverse`, and **proved two-sided inverse identities**. `permutationInjective` and `permutationSurjective` establish bijectivity. Subset transport is the direct image, represented through the inverse:

$$
\operatorname{subsetImage}(\sigma,S)(x)=S(\sigma^{-1}(x)).
$$

`subsetImageForward` checks the convention explicitly: membership at $\sigma(x)$ after transport equals membership at $x$ before transport. `subsetImageSwap` proves that left-composing the permutation by a transposition gives exactly the Boolean swap used by the existing canonical subset-history formalization.

For $P_n=G_{n-1}\cdots G_0$, the final declaration has the following logical order:

1. Every entry of the single `coinLaw` has a nonnegative rational weight.
2. Its weights sum to one.
3. **For all** Boolean subsets $S,T$,

   $$
   |\det(P_n[T,S])|^2
   =\sum_{(\sigma,w)\in\operatorname{coinLaw}}w\,\mathbf 1_{\sigma(S)=T}.
   $$

Both index lists use the order of `vertices.filter`. The permutation law is fixed before the universal subset quantifiers, so this is simultaneous common weighting across all minor orders, rather than one separately chosen law per minor.

## Why collision sectors do not change the law

The global coin law always contains both stay and swap choices. If a chosen subset occupies both endpoints or neither endpoint, these choices give the same subset. The canonical subset-history list counts that stay only once. `successorExpectation` proves exactly that the two contributions aggregate:

$$
\sum_{T\in\operatorname{successors}(S)}p(S,T)f(T)
=(1-q)f(S)+qf(\tau S),
$$

including $\tau S=S$. No division, positive-weight assumption or conditioning is used. `lawExpectationStep` and `historyExpectationStep` therefore obey the same recurrence. `commonPushforward` proves equality for **every rational-valued test function**, every starting subset and every finite prefix. Taking an endpoint indicator gives `commonEndpointLaw`; taking the constant function one gives `coinLawMassOne`.

These push-forward and normalization statements need neither a forest nor unitarity. Without probability bounds, their weights can be signed. `coinWeightsNonnegative` proves positivity from $0\le q_k\le1$, and `unitaryCoinWeightsNonnegative` derives those bounds from the actual unitary block equations.

Finally, `commonSquaredMinorLaw` combines the internally proved push-forward with `ForestSubsetExchange.productMagnitude`. That imported theorem identifies the **actual** forest gate-product minor with the canonical endpoint probability. No matrix/history identity, common-law assumption, or Cauchy–Binet theorem is added as a premise here.

## Domain and exact assumptions

- Lean **4.24.0**, `Std`, no Mathlib; Gaussian scalars are exact `Rat × Rat`.
- `vertices` is finite, duplicate-free and complete through the existing `Edge vertices` certificates; each gate has two distinct endpoints. Since `edges : Nat → Edge vertices` exists even for an empty used prefix, this API does not instantiate on zero- or one-vertex ambient types. Empty and full **selected subsets** are included.
- The squared-minor conclusion assumes `EveryEdgeIsBridge edges n`: deleting each indexed edge disconnects its endpoints within the used prefix. This excludes repeated parallel edges and cycles. No equivalence with a standard graph-library acyclicity predicate is asserted.
- Every used block satisfies the explicit `UnitaryBlock` equations. Vanishing transition weights, collision sectors, singular selected minors, unequal selected cardinalities and zero-length prefixes are included.
- Each listed permutation is an actual bijection of the ambient vertex type, with its inverse supplied by construction. A permutation may occur more than once; the common law remains a normalized finite weighted list.

## Replay and review

The source has 15 named theorem reports. A fresh build checks these and all 286 earlier declarations: **301 total**. Only `propext`, `Classical.choice`, and `Quot.sound` occur in the transitive axiom reports. No `sorry`, `admit`, `native_decide`, or new axiom is used. The manifest records the checked source/dependency hashes and exact declaration inventory.

After the separate checker registration lands:

```sh
python3 formal/check.py /absolute/path/to/lean
```

The source-only commit can be replayed before that integration with:

```sh
LEAN=/absolute/path/to/lean
BUILD_DIR=$(mktemp -d)
for MODULE in ForestCut ForestTrajectory ForestGateSupport ForestGateMagnitude ForestPathProbability ForestPathNormalization ForestEndpointMass ForestProductMinor ForestPairHistory ForestMinorAlternation ForestSubsetExchange ForestPermutationLaw; do
  LEAN_PATH="$BUILD_DIR" "$LEAN" -DwarningAsError=true -o "$BUILD_DIR/$MODULE.olean" "formal/$MODULE.lean" || exit 1
done
```

Independent-owner review should inspect the actual inverse convention, composition order, law independence from selected subsets, collision aggregation, positivity premises, and the exact imported actual-minor theorem, then cite the submitted immutable checkpoint. Source merge alone is not a mathematical review verdict.

## What remains

The common permutation law is now formalized for this Gaussian-rational forest gate family. A formal determinant-expansion identity connecting this common-minor law to the determinant of a sum of normal matrices and its permutation-product convex hull remains separate. So do an arbitrary-complex library interface and any continuity/density transfer. No unrestricted MDO statement or general arbitrary-matrix all-order Cauchy–Binet theorem follows merely from naming this law.
