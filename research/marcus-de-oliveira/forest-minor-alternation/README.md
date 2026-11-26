# Arbitrary-order row alternation for the actual ordered minor

`MathCommons.ForestMinorAlternation.gateMinorBoth` proves the actual both-occupied gate update without the three unproved alternation premises retained by [PR47](https://github.com/g4nesh/math-commons/pull/47). The proof establishes arbitrary-position row exchange, duplicate-row cancellation, and the exact sign of moving a row past a block of spectators for the existing signed Laplace determinant.

This is an arbitrary-minor-order algebraic result over Gaussian rationals. It is **not** a proof of all-order Cauchy–Binet, the all-order forest product probability law, or MDO. The full tuple-to-canonical-subset grouping remains separate. The artifact is submitted for review by a different GitHub owner; the author's compiler run is not independent-owner review.

## Exact actual gate formula

Let $G$ be the two-coordinate gate with block $\begin{pmatrix}a&b\\c&d\end{pmatrix}$ on the distinct endpoints `e.left`, `e.right`. Let `rows` be a duplicate-free list containing both endpoints, and let `cols` be any list. Then:

$$
\det((GM)[\mathrm{rows},\mathrm{cols}])
=(ad-bc)\det(M[\mathrm{rows},\mathrm{cols}]).
$$

`gateMinorBoth` uses the actual `matrixMul`, `twoCoordinateGate`, and `minorDet` definitions. It assumes the existing `Edge` certificate, including a complete duplicate-free vertex enumeration. It assumes neither unitarity, a forest schedule, square minors, nonzero determinants, nor any division condition. When the row and column lengths differ, the same identity holds for their defined rectangular determinant values.

The only premises beyond those concrete definitions are `rows.Nodup` and membership of the two gate endpoints. The earlier `gateMinorBothOfAlternation` remains a correctly stated conditional theorem in its original module; this new result supplies each of its three premises by proved lemmas.

## Reusable statements

| Declaration | Checked result |
| --- | --- |
| `doubleChoiceSwap` | The double signed column-selection sum changes sign when its two selected arguments are exchanged. No distinct-column assumption. |
| `rowDetHeadSwap`, `rowDetAdjacentSwap` | An adjacent row-vector exchange negates the actual determinant, at the head or after any prefix. |
| `rowDetExchange` | Exchange of any two row positions negates the determinant while preserving all intervening positions. |
| `rowDetDuplicate` | Two equal row vectors, separated by an arbitrary block, give zero. |
| `minorEqualRowsZero` | Two distinct selected row indices with equal row vectors force the actual matrix minor to vanish. The row list need not be duplicate-free. |
| `minorDuplicateRow` | Replacing one selected matrix row by another selected row gives zero. |
| `minorExchangeRows` | The explicit two-row `replaceRow` operation has negative determinant sign on a duplicate-free selected row list. |
| `rowDetMove`, `minorMoveRow` | Moving one row across a block of selected spectators has exactly the parity of that block. |
| `gateMinorBoth` | The unconditional actual both-occupied gate formula above. |

For example, `minorMoveRow` says, for arbitrary ordered lists $B,L,A$ and any selected index $r$:

$$
\det M[B\mathbin{+\!+}L\mathbin{+\!+}[r]\mathbin{+\!+}A,C]
=(-1)^{|L|}\det M[B\mathbin{+\!+}[r]\mathbin{+\!+}L\mathbin{+\!+}A,C].
$$

No adjacency of vertex labels is required. This theorem expresses a list-order identity; it does not silently identify a row replacement with a Boolean subset's sorted filter.

## Proof structure

The first two Laplace selections form a double sum. Splitting the column list at its head separates three terms: the head selected first, the head selected second, and neither selection at the head. The first two contributions have opposite signs. Induction handles the remaining double sum after retaining the unselected head. This establishes `doubleChoiceSwap` for arbitrary remaining-row expressions and arbitrary column lists.

Multiplication commutativity converts that finite-sum identity into adjacent row exchange. A common row prefix preserves the negative relation. Three exchanges extend the result across any intervening block. Setting the two row vectors equal yields $D=-D$, hence $D=0$ in Gaussian-rational arithmetic. No cancellation by a potentially zero determinant is used.

A two-endpoint list decomposition translates the row-vector identities to actual matrix row replacements. Duplicate-free selected rows ensure that exchanging the two matrix rows changes exactly those two selected positions. Finally the four-term actual gate expansion from PR47 collapses to the block determinant factor.

## Scope still open

- Group the arbitrary-order ordered-tuple product expansion into canonical occupied subsets, including its multiplicity and permutation signs, to obtain full Cauchy–Binet.
- Derive the single-occupied gate row-replacement formula directly in `vertices.filter` order using the proved movement parity and explicit filter decomposition. The list-order movement sign is proved here; this Boolean-subset bridge is not assumed or claimed.
- Connect the resulting arbitrary-order product formula to the canonical history law. The actual order-two bridge is already proved by [PR53](https://github.com/g4nesh/math-commons/pull/53).
- Connect these Gaussian-rational definitions to a standard arbitrary-complex matrix/determinant library and any resulting MDO convex-hull theorem.

## Replay and audit

[Source](../../../formal/ForestMinorAlternation.lean), [declaration manifest](declarations.json), and [audit manifest](audit.json). The source has 23 named theorem reports. The source-commit fresh rebuild checks 243 declarations including all 220 existing reports. Transitive axioms are restricted to `propext`, `Classical.choice`, and `Quot.sound`; no `sorry`, `admit`, `native_decide`, or new axiom is used. Source and dependency hashes identify exactly what was checked.

With this PR's checker registration present, from the repository root:

```sh
python3 formal/check.py /absolute/path/to/lean
```

To replay the source-only commit before its separately reserved checker integration:

```sh
LEAN=/absolute/path/to/lean
BUILD_DIR=$(mktemp -d)
for MODULE in ForestCut ForestTrajectory ForestGateSupport ForestGateMagnitude ForestPathProbability ForestPathNormalization ForestEndpointMass ForestProductMinor ForestMinorAlternation; do
  LEAN_PATH="$BUILD_DIR" "$LEAN" -DwarningAsError=true -o "$BUILD_DIR/$MODULE.olean" "formal/$MODULE.lean" || exit 1
done
```

Use Lean 4.24.0 and `Std`, without Mathlib. Review should inspect the actual definitions, row-list quantifiers, sign orientation and transitive axiom reports, then cite the immutable submitted coordination checkpoint. Protected source merge alone is not an independent mathematical review.
