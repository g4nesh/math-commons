# Lean verification

Lean 4.24.0 checks an exact three-dimensional instance of the Marcus–de Oliveira conjecture. The proof starts from the spectra and relative eigenbasis, computes the determinant and the six permutation products, and constructs a convex combination for the determinant. The tree on the website draws these data. The [Lean source](MathCommons.lean) uses Lean's standard library.

## Reproduce it

From the repository root, with [elan](https://lean-lang.org/install/) installed, run:

```sh
python3 formal/check.py
```

The script uses the toolchain pinned in this directory and treats warnings as errors. It checks the transitive axiom report for each named theorem. The allowlist contains Lean's standard propositional extensionality, classical choice and quotient soundness. A missing report or an extra axiom fails the check. The same command runs in the required GitHub check. The executable path can be supplied as the first argument when using a standalone installation of the pinned version.

## What is checked

Write the spectra as $a=(0,1,i)$ and $b=(-1-i,1,-1+i)$. The relative eigenbasis is the real symmetric matrix with diagonal entries −1/3 and off-diagonal entries 2/3. The formalization represents complex rationals as pairs of rationals with explicit complex multiplication.

| Lean declaration | Result |
| --- | --- |
| `allPairings` | Exhausts all triples of distinct indices from a three-element set. Each triple defines one pairing of the two spectra. |
| `orthogonal` | Checks all nine entries of the row Gram matrix against the identity. |
| `determinantValue` | Evaluates the three-by-three Leibniz determinant of the diagonal matrix plus the conjugated second diagonal matrix. Its coordinates are (2/3, 2/9). |
| `witnessInHull` | Constructs two binary convex combinations of three computed endpoints. The endpoints (2,1), (0,−1), (0,2) have weights 9/27, 13/27 and 5/27. |
| `inEveryConvexSet` | Proves by induction that every constructed certificate belongs to every set containing the vertices and closed under rational convex combinations. |

The rational coordinates embed in the real plane. Addition, multiplication and nonnegative rational interpolation have their usual meaning under this embedding. The code formalizes the rational operations and interpolation predicate directly. The bridge to Mathlib's complex matrices, spectral theorem, determinant and real convex hull is still a separate formalization task. The [independent arithmetic implementation](../research/marcus-de-oliveira/literature/check_source_examples.py) reproduces the example.

## An unbounded cut-history lemma

[ForestCut.lean](ForestCut.lean) also proves the arithmetic step behind the forest rotation proof. For any finite number of steps, if the integer count on one side of a cut can change at only one step, its signed change is fixed by the endpoint counts. Two histories with the same endpoints have identical counts throughout that cut. If a crossing decision is encoded as stay, enter or leave, the signed change also determines the decision. The checked declarations include `cutHistoryUnique` and `crossingDecisionUnique`.

## From edge cuts to occupied subsets

In [ForestTrajectory.lean](ForestTrajectory.lean), an occupied subset is a Boolean function on a finite vertex set. At each indexed edge a step either keeps the current subset or swaps the membership of its two endpoints. The source proves cardinality preservation and handles both endpoints occupied, neither occupied, and empty and full subsets.

`MathCommons.ForestTrajectory.forestTrajectoryUnique` assumes that deleting each indexed edge disconnects its endpoints. For any two legal paths with the same initial and final subsets, their subsets agree at every intermediate step. The proof constructs a cut from reachability after deleting the edge, proves each surviving edge stays on one side, then applies the arithmetic result in `ForestCut`. This establishes uniqueness of occupied subsets. Swaps can still describe the same subset when both endpoints have the same occupancy.

`python3 formal/check.py` rebuilds the imports in a fresh directory and audits all 17 declarations in this module. A Mathlib interface could also supply the vertex enumeration and indexed-edge bridge condition from its standard finite graph definitions.

## Which minors of a two-coordinate gate can be nonzero?

[`ForestGateSupport.lean`](ForestGateSupport.lean) defines a signed Laplace expansion on ordered lists of row and column indices. The proof extracts a matching of nonzero entries from any nonzero expansion. In a matrix supported on the gate's two coordinates and the diagonal, a matching must pair every occupied coordinate outside the gate with itself. The matching also gives equal cardinalities. These facts imply that the output subset either equals the input subset or exchanges its two gate coordinates. Equal occupancy of the two endpoints forces the subset to stay fixed.

`MathCommons.ForestGateSupport.twoCoordinateNonzeroMinorStep` accepts four arbitrary block entries. Its scalar interface specifies addition, multiplication, negation, zero and one, with the zero identities used in the proof. The support implication applies to singular blocks. Rational and Gaussian-rational instances instantiate the interface.

`MathCommons.ForestGateSupport.nonzeroAmplitudesUnique` combines this support result and the forest theorem. For fixed endpoint subsets, at most one history has a nonzero ordered product of its gate minors. Empty and full subsets, repeated occupancy and vanishing factors are covered. The checker audits all 23 named results in the module. The local unitary computation below supplies the transition magnitudes. Identifying the ordered expansion with the standard complex determinant and composing the minors by Cauchy–Binet remain open. [Written forest argument](../research/marcus-de-oliveira/baselines/README.md)

## Squared probabilities for a single gate

[ForestGateMagnitude.lean](ForestGateMagnitude.lean) uses Gaussian-rational entries and the ordered Laplace expansion in `ForestGateSupport`. For the gate with block entries $a,b,c,d$ in row order, let $p=|b|^2$. Unitarity gives $|c|^2=p$, $|a|^2=|d|^2=1-p$, and $|ad-bc|^2=1$. The following table describes occupied **subsets**, so a transposition acting on two occupied endpoints produces the same subset.

| Occupied endpoints in the input | Output subset | Squared minor |
| --- | --- | --- |
| Neither or both | Same subset | 1 |
| Exactly one | Same subset | $1-p$ |
| Exactly one | Occupied endpoint exchanged | $p$ |
| Any | A change outside the gate, or another illegal transition | 0 by the [support theorem](ForestGateSupport.lean) |

The checked declaration `MathCommons.ForestGateMagnitude.orderedGateMagnitudes` takes three lists of occupied spectators, covering those before the first endpoint, between the endpoints and after the second. `decomposeTwoEndpoints` and `filterTwoEndpoints` express the corresponding ordered-list decompositions; `reverseGate` and `reverseUnitary` handle the opposite endpoint order. The subsequent path module assembles the cases for arbitrary Boolean subset histories. The module also checks $0≤p≤1$.

Both the empty and the full subset occur in the first row. Occupied spectator coordinates may precede, separate or follow the two gate endpoints. An exchange picks up the parity of the occupied spectators between the endpoints in the signed minor. Squaring the norm removes that sign. The local result and the forest path-uniqueness theorem have different hypotheses. To combine them into a matrix-product identity still requires a Cauchy–Binet expansion in the same minor convention.

## Contributing a theorem

1. Claim the specific lemma or translation in the live ledger. State its quantifiers, all hypotheses, the Lean declaration you intend to prove and the source file you own. Separate mathematical hypotheses from any unproved bridge from the paper's statement.
2. Pin dependencies. For Mathlib, commit its revision and Lake manifest with a matching Lean toolchain. Agree on matrix indices, normality, conjugate transpose, the scalar field for convexity and the direction of each permutation before building dependent results.
3. Cover repeated eigenvalues, singular sums, point or segment hulls and the intended dimension range. If a proof divides by a quantity, prove it nonzero or give the zero case. A strict-separation argument needs its exact sign and quantifiers.
4. Add a reproducible checker to the required CI job. For this project add a named `#print axioms` report, list the declaration in `check.py`, and keep warnings fatal. Inspect the statement and imports as well as its axiom report. Drafts with placeholders belong in clearly marked notes until they pass. Avoid `native_decide` or introducing an axiom to discharge the claim. The current finite computations use kernel reduction.
5. Checkpoint the declaration names, exact commit, local checker result, and unresolved dependencies. Submit the evidence to a different GitHub owner, who should rebuild it and compare the formal statement against the mathematical claim. Merging the source makes it available to collaborators. A ledger review records the independent verdict.

A natural next lemma is the relation between the rational interpolation constructor and convexity of subsets of the complex plane over the reals. A separate, larger direction is to formalize the coefficient identity in the [forest-rotation argument](../research/marcus-de-oliveira/baselines/README.md). Neither the finite witness nor that restricted family establishes the unrestricted conjecture.

[ForestPathProbability.lean](ForestPathProbability.lean) proves that the squared amplitude of each ordered path equals the product of its gate transition probabilities. The required check compiles this module after its dependencies and audits the transitive axioms of all fourteen printed declarations, including `pathMagnitude` and `nonzeroProbabilityUnique`. Connecting these finite path products to standard complex matrix minors by Cauchy–Binet remains open.

## Total probability over finite histories

[ForestPathNormalization.lean](ForestPathNormalization.lean) enumerates legal occupied-subset histories for any finite gate schedule and any starting subset. It counts a stay once when a swap leaves the same subset. `pathMassOne` gives total mass one. For unitary Gaussian-rational gates, `pathProbabilityNonnegative` proves each mass is nonnegative and `pathSquaredMagnitudeMassOne` gives total squared path amplitude one. Cyclic schedules are allowed. Connecting amplitudes that have the same endpoint to the minor of a matrix product still requires an interference argument and Cauchy–Binet. [Statements, assumptions and replay](../research/marcus-de-oliveira/forest-path-normalization/README.md)
