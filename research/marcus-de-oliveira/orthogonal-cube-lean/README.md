# Lean formalization of the unit-cube quadratic lemma

Work in progress on Claude's Lemma L in [PR48](https://github.com/g4nesh/math-commons/pull/48), source commit `ec52b46c2e2a9ed21947fcf6d91a27b1c56d6658`.

[`formal/OrthogonalCube.lean`](../../../formal/OrthogonalCube.lean) currently proves the two-dimensional positive-semidefinite criterion, its Lagrange identity, both bilinear squared inequalities, the specialized F1/F2 consequences, complement polarization, and the balanced value at one. All ten named declarations compile under Lean4.24.0 with warnings fatal; their only transitive axioms are `propext`, `Classical.choice`, and `Quot.sound`.

The domain is an arbitrary `Lean.Grind.Field R` with `LE R`, `LT R`, `Std.LawfulOrderLT R`, `Std.IsLinearOrder R`, and `Lean.Grind.OrderedRing R`. These are the standard field and compatible linear-order laws in Lean's core hierarchy. This is a faithful general ordered-field formulation, not a rational sample check or a rational-only theorem. No completeness or square-root operation is assumed. A concrete Mathlib `Real` instance is not imported.

**Still in progress:** derive F4 from nonnegative cube coordinates, combine the strict inequalities, and remove the product-order case assumption by the symmetry. The final Lemma L is not claimed by this initial draft. The SO(3)/SO(4), Hodge, matrix-minor, and PMSH reductions are outside this task.

```sh
lean -DwarningAsError=true formal/OrthogonalCube.lean
```

Same-owner work; independent-owner statement review is required.
