# Actual forest product minors, at every subset size

`MathCommons.ForestSubsetExchange.productMagnitude` proves that the squared norm of an **actual minor of the actual gate product** equals the canonical endpoint probability for a unitary forest schedule. There is no order-two or subset-cardinality premise. This closes the product-minor bridge left open by the earlier forest path results, over their existing Gaussian-rational scalar domain.

[Lean source](../../../formal/ForestSubsetExchange.lean), [declarations](declarations.json), [audit](audit.json). The author's kernel replay is not independent-owner review. This result does not assert MDO or a formal arbitrary-complex matrix-library bridge.

## Checked statement and conventions

Fix a finite, duplicate-free, complete vertex enumeration `vertices`. Each `Edge vertices` specifies distinct gate endpoints and the existing enumeration certificates. Write

$$
G_k=I\text{ outside its two endpoints},\qquad
G_k\big|_{\{l_k,r_k\}}=\begin{pmatrix}a_k&b_k\\c_k&d_k\end{pmatrix},\qquad
P_n=G_{n-1}\cdots G_0.
$$

All entries belong to `Gaussian := Rat × Rat`, with the repository's exact complex arithmetic. For Boolean subsets $S,T$, rows and columns use **the order inherited from `vertices.filter`**. For every natural $n$, `productEndpointAmplitude` proves

$$
\det(P_n[T,S])=\sum_{h\in\mathcal H_n(S,T)}\prod_{k=0}^{n-1}\det(G_k[h_{k+1},h_k]).
$$

This signed identity needs neither unitarity nor a forest hypothesis. The finite list $\mathcal H_n$ is the existing canonical history enumeration: one successor in a collision sector, and stay plus swap in a single-occupied sector. It includes legal paths even when an amplitude happens to be zero; it has no duplicate histories.

If `EveryEdgeIsBridge edges n` and every used block satisfies `UnitaryBlock`, `productMagnitude` proves

$$
\left|\det(P_n[T,S])\right|^2
=\operatorname{endpointProbability}(S,T,n).
$$

The right side is the sum of the existing path probabilities, with local stay/swap weights $1-|b_k|^2$ and $|b_k|^2$ in single-occupied sectors and weight one in collision sectors. The forest condition is the precise existing bridge predicate on the finite prefix; no graph acyclicity premise is silently inferred from a name.

The exact Lean conclusion is:

```lean
normSq (minorDet gaussianArithmetic
  (matrixProduct vertices (gateSchedule edges as bs cs ds) n)
  (vertices.filter finish) (vertices.filter start)) =
    endpointProbability edges bs start finish n
```

Its explicit non-typeclass premises are `bridges : EveryEdgeIsBridge edges n` and `unitary : ∀ k, k < n → UnitaryBlock (as k) (bs k) (cs k) (ds k)`. The underlying `edges : Nat → Edge vertices` includes a complete duplicate-free vertex list and two distinct endpoints. This representation therefore presupposes an ambient type supporting such an edge, even at $n=0$; it does not represent zero- or one-vertex ambient spaces. Empty and full **selected subsets** are covered, as are unequal row/column cardinalities and unreachable endpoint fibers under the existing rectangular-minor definition.

`productEmptyFiber` proves an actual zero minor whenever the canonical endpoint list is empty, without requiring unitarity or a forest. `productMassOne` proves the squared actual minors sum to one over `reachableEndpoints`: a duplicate-free finite list of the distinct endpoints of legal histories. This list uses classical equality of Boolean functions. It is not advertised as an efficient endpoint-enumeration algorithm. Unreachable endpoints contribute zero; normalization does not count an endpoint once per history.

## Why the sparse bridge is enough

1. The four `leftReplacementForward/Backward` and `rightReplacementForward/Backward` identities translate an actual selected-row replacement into the minor on the swapped Boolean subset. The sign is exactly $(-1)^q$, where $q$ is the number of **selected** spectators between the endpoints in the vertex list. Both endpoint orders and both occupied orientations are proved.
2. `gateTwoCoefficients` uses the actual gate-row formulas and that parity to obtain two coefficients independent of the remaining matrix and column list. Its exchange coefficient is zero in collision sectors.
3. Specializing the remaining matrix to the actual identity identifies those coefficients with the actual gate minors. `gateMinorPredecessors` consequently proves, for an arbitrary column list $C$,

   $$
   \det((GM)[T,C])=\sum_{R\in\operatorname{successors}(T)}
      \det(G[T,R])\det(M[R,C]).
   $$

   There are one or two predecessor subsets. This is a proved sparse Cauchy–Binet identity for a two-coordinate left factor, not a supplied transport hypothesis.
4. `predecessorMinorSum` contracts the endpoint indicator. For a non-neighbor the actual gate minor is zero, using the earlier gate-support theorem and symmetry of the swap relation. Induction on the product length then gives `productEndpointAmplitude` without a minor-size restriction.
5. The existing forest path uniqueness and single-gate magnitude results imply `endpointMagnitude`; combining it with the actual product identity gives `productMagnitude`. Finite-sum exchange and duplicate-free endpoint enumeration give `productMassOne`.

No determinant division, nonzero-minor assumption, numerical interpolation, assumed all-order Cauchy–Binet theorem, or assumed product transport identity enters this proof.

## Replay and audit

Use Lean **4.24.0** and `Std`, without Mathlib. The source reports all 25 new named theorems. A fresh dependency build checks their transitive axiom reports and the existing 261 declarations: 286 in total. Allowed logical axioms are only `propext`, `Classical.choice`, and `Quot.sound`; no `sorry`, `admit`, `native_decide`, or new axiom is used. The audit manifest pins the source and dependency bytes.

After checker integration, from the repository root:

```sh
python3 formal/check.py /absolute/path/to/lean
```

For the source-only commit before checker integration:

```sh
LEAN=/absolute/path/to/lean
BUILD_DIR=$(mktemp -d)
for MODULE in ForestCut ForestTrajectory ForestGateSupport ForestGateMagnitude ForestPathProbability ForestPathNormalization ForestEndpointMass ForestProductMinor ForestPairHistory ForestMinorAlternation ForestSubsetExchange; do
  LEAN_PATH="$BUILD_DIR" "$LEAN" -DwarningAsError=true -o "$BUILD_DIR/$MODULE.olean" "formal/$MODULE.lean" || exit 1
done
```

Review should inspect the exact scalar domain, `matrixProduct` multiplication order, Boolean-filter signs, one-versus-two predecessor sectors, and all transitive axiom reports, then cite the immutable submitted checkpoint under a different GitHub owner.

## Remaining boundaries

The generic arbitrary-matrix, arbitrary-order Cauchy–Binet tuple-grouping theorem remains separate; this gate-specific sparse proof bypasses it. A formal connection from Gaussian rationals to a standard arbitrary-complex determinant library, a law on full permutations rather than one occupied subset at a time, and an MDO convex-hull conclusion remain outside this module. None is assumed in the proved actual-minor probability identity.
