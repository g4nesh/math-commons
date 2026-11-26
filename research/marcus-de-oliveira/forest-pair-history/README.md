# Actual order-two forest minors and canonical probabilities

`MathCommons.ForestPairHistory.productPairMagnitude` proves that the squared norm of an **actual 2×2 minor of an arbitrary-length gate product** equals the probability of its endpoint in the canonical legal-subset history enumeration. `productPairMassOne` proves that these actual squared minors sum to one across all canonical output pairs. This closes the order-two matrix-to-history gap left by [PR47](https://github.com/g4nesh/math-commons/pull/47).

The theorem is over Gaussian rationals, with arbitrary finite ambient dimension and arbitrary finite schedule length. It does not establish all minor orders, a standard complex-matrix interface, or the general Marcus–de Oliveira conjecture. The formal artifact awaits review by a different GitHub owner; the author's kernel run is not an independent-owner review.

## Exact statement

Let `vertices` be the exhaustive, duplicate-free enumeration supplied by `Edge vertices`. Write

$$
G_k=\texttt{twoCoordinateGate}(e_k,a_k,b_k,c_k,d_k),\qquad
P_0=I,\qquad P_{n+1}=G_nP_n.
$$

A subset is a Boolean function on vertices. Filtering `vertices` gives its canonical row or column order. Suppose the initial and final filters each have length two. Then:

$$
\det(P_n[\mathrm{finish},\mathrm{start}])
=\sum_{h\in\mathcal H_n(\mathrm{start},\mathrm{finish})}
 \prod_{k=0}^{n-1}\det(G_k[h_{k+1},h_k]).
$$

This is `productPairEndpointAmplitude`. It needs **neither unitarity nor the forest hypothesis**. Every determinant and matrix product here uses the existing explicit Gaussian-rational arithmetic and signed Laplace determinant; the equality is proved, not supplied as a transport premise.

If each block is `UnitaryBlock` and `EveryEdgeIsBridge edges n` holds, then:

$$
\left|\det(P_n[\mathrm{finish},\mathrm{start}])\right|^2
=\texttt{endpointProbability}(\mathrm{start},\mathrm{finish},n).
$$

This is `productPairMagnitude`. The right side is the sum of the existing `pathProbability` over canonical histories ending at `finish`. Each step stays or swaps the gate endpoints. For single occupancy, the exchange probability is $|b_k|^2$ and the stay probability is $1-|b_k|^2$; equal occupancies have one stay state with probability one. The existing unitary theorem supplies the nonnegative step probabilities.

For a fixed two-element start, `productPairMassOne` states:

$$
\sum_{(i,j)\in\texttt{orderedPairs}(\mathrm{vertices})}
 \left|\det(P_n[[i,j],\mathrm{start}])\right|^2=1.
$$

The sum includes every ordered-by-enumeration pair exactly once. Nonnegativity also follows from `ForestGateMagnitude.normSqNonnegative`.

## How the missing bridge is proved

1. `pairMemberSublist`, `pairFilter`, `twoSubsetRepresentative`, and `pairSubsetInjective` prove the canonical pair/subset correspondence. `pairSubsetsNodup` prevents multiplicities.
2. `identitySubsetMinor` proves the actual identity-matrix minor is one for equal subsets and zero otherwise. It applies to all subset sizes. The nonzero-matching support theorem rules out unequal subsets; identity spectator elimination proves the diagonal value.
3. `endpointAmplitudeRecurrence` derives the recurrence from the existing finite history enumeration. `successorMinorSum` uses actual gate-minor support to show that omitted illegal transitions vanish. Collision stays are deduplicated by the existing `successors` list.
4. `historySizeTwo` proves every enumerated intermediate subset still has size two. PR47's actual order-two Cauchy–Binet recurrence is contracted against the finite history sums. A unique canonical pair represents each history endpoint.
5. `productPairMagnitude` applies the existing forest endpoint no-interference theorem to the now-identified actual minor. `endpointPairMassOne` partitions all canonical histories by their unique endpoint pair; existing path normalization yields `productPairMassOne`.

No generic determinant alternation premise, desired product-minor equality, or probability normalization is assumed.

## Boundaries and degeneracies

- Schedule length zero is included: the actual product is the identity and the history enumeration contains the constant history.
- Gate endpoints may be nonadjacent in the supplied enumeration; signed minors retain that enumeration's order.
- Repeated occupancy, empty endpoint fibers, and zero amplitudes are allowed. `productPairEmptyFiber` separately proves the actual minor is zero for an empty endpoint fiber, without unitarity or a forest assumption.
- `EveryEdgeIsBridge` uses explicit reachability after deleting each indexed edge in the finite schedule. This module does not replace it by an external graph library's forest predicate.
- Only order-two product minors are connected here. The arbitrary-order tuple expansion exists in PR47, but its full alternating-subset Cauchy–Binet reduction remains a separate task.
- Gaussian rationals are pairs of `Rat`. No completeness/continuity argument to arbitrary complex entries, library matrix identification, determinant convex-hull conclusion, or universal MDO proof is claimed.

## Source and replay

- [Lean source](../../../../formal/ForestPairHistory.lean): 28 declarations with printed transitive axiom reports.
- [Declaration manifest](declarations.json): exact audited names.
- [Audit manifest](audit.json): source hashes, dependency hashes, compiler version, and axiom dependencies.

Lean 4.24.0, `Std`, no Mathlib. With the compiler installed and this PR's checker registration present, run from the repository root:

```sh
python3 formal/check.py /absolute/path/to/lean
```

For the source-only commit before CI registration, the following builds the same dependencies and the new module in a fresh directory, with warnings fatal:

```sh
LEAN=/absolute/path/to/lean
BUILD_DIR=$(mktemp -d)
for MODULE in ForestCut ForestTrajectory ForestGateSupport ForestGateMagnitude ForestPathProbability ForestPathNormalization ForestEndpointMass ForestProductMinor ForestPairHistory; do
  LEAN_PATH="$BUILD_DIR" "$LEAN" -DwarningAsError=true -o "$BUILD_DIR/$MODULE.olean" "formal/$MODULE.lean" || exit 1
done
```

The full validation also checks the existing `MathCommons.lean` artifact. All 220 declaration reports in that run are restricted to `propext`, `Classical.choice`, and `Quot.sound`. There are no new axioms, `sorry`, `admit`, or `native_decide` in this module.

Review request: replay the pinned source and dependency hashes, inspect the canonical enumeration and successor contraction, and confirm the exact theorem assumptions above. Please record any independent verdict against the immutable submitted checkpoint under your own GitHub owner.
