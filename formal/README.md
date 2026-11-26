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

`MathCommons.ForestGateSupport.nonzeroAmplitudesUnique` combines this support result and the forest theorem. For fixed endpoint subsets, at most one history has a nonzero ordered product of its gate minors. Empty and full subsets, repeated occupancy and vanishing factors are covered. The checker audits all 23 named results in the module. The local unitary computation below supplies the transition magnitudes. Actual gate-product composition and the canonical history bridge at every subset size are now proved in `ForestSubsetExchange` below; the standard arbitrary-complex determinant interface remains open. [Written forest argument](../research/marcus-de-oliveira/baselines/README.md)

## Squared probabilities for a single gate

[ForestGateMagnitude.lean](ForestGateMagnitude.lean) uses Gaussian-rational entries and the ordered Laplace expansion in `ForestGateSupport`. For the gate with block entries $a,b,c,d$ in row order, let $p=|b|^2$. Unitarity gives $|c|^2=p$, $|a|^2=|d|^2=1-p$, and $|ad-bc|^2=1$. The following table describes occupied **subsets**, so a transposition acting on two occupied endpoints produces the same subset.

| Occupied endpoints in the input | Output subset | Squared minor |
| --- | --- | --- |
| Neither or both | Same subset | 1 |
| Exactly one | Same subset | $1-p$ |
| Exactly one | Occupied endpoint exchanged | $p$ |
| Any | A change outside the gate, or another illegal transition | 0 by the [support theorem](ForestGateSupport.lean) |

The checked declaration `MathCommons.ForestGateMagnitude.orderedGateMagnitudes` takes three lists of occupied spectators, covering those before the first endpoint, between the endpoints and after the second. `decomposeTwoEndpoints` and `filterTwoEndpoints` express the corresponding ordered-list decompositions; `reverseGate` and `reverseUnitary` handle the opposite endpoint order. The subsequent path module assembles the cases for arbitrary Boolean subset histories. The module also checks $0≤p≤1$.

Both the empty and the full subset occur in the first row. Occupied spectator coordinates may precede, separate or follow the two gate endpoints. An exchange picks up the parity of the occupied spectators between the endpoints in the signed minor. Squaring the norm removes that sign. The local result and the forest path-uniqueness theorem have different hypotheses. The order-two identity is proved below by Cauchy–Binet, and `ForestSubsetExchange` proves arbitrary-order gate-product composition by a sparse predecessor recurrence in this same minor convention.

## An unconditional quadratic lemma on the closed unit cube

[`OrthogonalCube.lean`](OrthogonalCube.lean) proves Claude's Lemma L over an arbitrary field with a compatible lawful linear order. For nonnegative $p,q,r,u$ satisfying $p+q=r+u$, put

$$Q(a,b,c,d)=pr(a^2+c^2)+qu(b^2+d^2)-2pqab-2rucd.$$

`MathCommons.OrthogonalCube.lemmaL` proves that, for every $(a,b,c,d)\in[0,1]^4$, at least one of $Q(a,b,c,d)$ and $Q(1-a,1-b,1-c,1-d)$ is nonnegative. `lemmaLMax` gives the equivalent maximum statement. All four auxiliary facts F1–F4 and the symmetry removing the product-order case assumption are proved internally. Zero parameters and all cube boundaries are included.

The standard `Lean.Grind.Field` and compatible ordered-ring hierarchy supply the algebra and order laws. This is a generic ordered-field proof, not a rational-only check; it does not import a concrete Mathlib `Real` instance. The required checker audits all 18 declarations, pins the source hash, checks the exact generic theorem type, and compiles a rational instance. Only `propext`, `Classical.choice`, and `Quot.sound` are allowed. The matrix, SO(3)/SO(4), Hodge, and PMSH reductions are **not formalized by this result**. [Proof organization, precise domain and reproducible audit](../research/marcus-de-oliveira/orthogonal-cube-lean/README.md)

## Contributing a theorem

1. Claim the specific lemma or translation in the live ledger. State its quantifiers, all hypotheses, the Lean declaration you intend to prove and the source file you own. Separate mathematical hypotheses from any unproved bridge from the paper's statement.
2. Pin dependencies. For Mathlib, commit its revision and Lake manifest with a matching Lean toolchain. Agree on matrix indices, normality, conjugate transpose, the scalar field for convexity and the direction of each permutation before building dependent results.
3. Cover repeated eigenvalues, singular sums, point or segment hulls and the intended dimension range. If a proof divides by a quantity, prove it nonzero or give the zero case. A strict-separation argument needs its exact sign and quantifiers.
4. Add a reproducible checker to the required CI job. For this project add a named `#print axioms` report, list the declaration in `check.py`, and keep warnings fatal. Inspect the statement and imports as well as its axiom report. Drafts with placeholders belong in clearly marked notes until they pass. Avoid `native_decide` or introducing an axiom to discharge the claim. The current finite computations use kernel reduction.
5. Checkpoint the declaration names, exact commit, local checker result, and unresolved dependencies. Submit the evidence to a different GitHub owner, who should rebuild it and compare the formal statement against the mathematical claim. Merging the source makes it available to collaborators. A ledger review records the independent verdict.

A natural next lemma is the relation between the rational interpolation constructor and convexity of subsets of the complex plane over the reals. A separate, larger direction is to formalize the coefficient identity in the [forest-rotation argument](../research/marcus-de-oliveira/baselines/README.md). Neither the finite witness nor that restricted family establishes the unrestricted conjecture.

[ForestPathProbability.lean](ForestPathProbability.lean) proves that the squared amplitude of each ordered path equals the product of its gate transition probabilities. The required check compiles this module after its dependencies and audits the transitive axioms of all fourteen printed declarations, including `pathMagnitude` and `nonzeroProbabilityUnique`. Connecting these finite path products to standard complex matrix minors by Cauchy–Binet remains open.

## Total probability over finite histories

[ForestPathNormalization.lean](ForestPathNormalization.lean) enumerates legal occupied-subset histories for any finite gate schedule and any starting subset. It counts a stay once when a swap leaves the same subset. `pathMassOne` gives total mass one. For unitary Gaussian-rational gates, `pathProbabilityNonnegative` proves each mass is nonnegative and `pathSquaredMagnitudeMassOne` gives total squared path amplitude one. Cyclic schedules are allowed. The endpoint and composition modules below now supply no-interference and actual gate-product matrix composition at every subset size. [Statements, assumptions and replay](../research/marcus-de-oliveira/forest-path-normalization/README.md)

## Squared amplitude at a fixed forest endpoint

[ForestEndpointMass.lean](ForestEndpointMass.lean) restricts the canonical histories to one final subset. For a forest gate schedule, at most one of those histories has a nonzero amplitude. `endpointMagnitude` proves that the squared norm of their signed sum equals their total probability, including zero for an unreachable endpoint. The signed sum is identified with the actual gate-product minor at every subset size by `ForestSubsetExchange` below. [Exact hypotheses and proof](../research/marcus-de-oliveira/forest-endpoint-mass/README.md)

## Actual matrix products and order-two Cauchy–Binet

[ForestProductMinor.lean](ForestProductMinor.lean) defines finite Gaussian-rational matrix multiplication and chronological gate products. `productMinorTupleExpansion` proves the general expansion of an actual product minor over ordered intermediate row tuples. It retains repeated indices explicitly. `cauchyBinetTwo` proves the complete Cauchy–Binet identity for two-by-two minors in **arbitrary ambient dimension**, and `iteratedCauchyBinetTwo` identifies their actual arbitrary-length product with the order-two compound-kernel evolution. Neither theorem assumes Cauchy–Binet, unitarity or a forest schedule. Orders zero and one and actual single-gate squared probabilities are also checked.

The module audits 52 declarations, including finite-sum identities and row-linearity at arbitrary spectator positions. The original higher-order collision theorem `gateMinorBothOfAlternation` retains three explicit determinant-alternation premises. The later `ForestMinorAlternation.gateMinorBoth` now discharges all three with actual proved determinant lemmas. The subsequent pair-history module identifies this actual order-two evolution with canonical Boolean-subset histories. Arbitrary-order tuple grouping and the standard complex matrix interface remain open; row alternation is proved below. [Precise statements, manifest and remaining obligations](../research/marcus-de-oliveira/forest-product-minor/README.md)


## Actual squared order-two minors form the endpoint law

[ForestPairHistory.lean](ForestPairHistory.lean) closes the order-two gap between actual gate-product minors and the canonical legal-history enumeration. The complete duplicate-free vertex list supplies the canonical order of every two-element subset. The module proves its pair/subset correspondence, the actual identity-minor delta, the successor recurrence, and vanishing of illegal transitions. It then contracts the actual order-two Cauchy–Binet recurrence against the finite history sums.

| Lean declaration | Exact scope |
| --- | --- |
| `productPairEndpointAmplitude` | For any finite gate sequence, the actual signed order-two product minor equals its canonical endpoint amplitude. No unitarity or forest hypothesis. |
| `productPairMagnitude` | With unitary Gaussian-rational blocks and `EveryEdgeIsBridge`, the actual squared order-two minor equals the canonical endpoint probability. |
| `productPairEmptyFiber` | An empty canonical endpoint fiber gives an actual zero order-two minor, without forest or unitarity assumptions. |
| `productPairMassOne` | For a two-element start and a forest unitary schedule, the actual squared minors sum to one over all canonical output pairs. |

Zero-length schedules, nonadjacent endpoints, collision stays and vanishing factors are included. No matrix-to-history equality, generic alternation premise or probability normalization is assumed. The required checker rebuilds this module after `ForestProductMinor` and audits all 28 reports in this module. The subsequent `ForestSubsetExchange` module removes the subset-size restriction using sparse predecessors. The standard arbitrary-complex matrix interface remains separate; this is not a general MDO proof. [Precise statement, assumptions and hashes](../research/marcus-de-oliveira/forest-pair-history/README.md)


## Arbitrary-order actual row alternation and gate updates

[ForestMinorAlternation.lean](ForestMinorAlternation.lean) proves that exchanging two row positions negates the actual ordered Gaussian-rational determinant, and two equal row vectors force it to vanish. The proof pairs two successive signed column selections, then extends adjacent exchange across an arbitrary spectator block. `minorMoveRow` gives the exact parity when a row moves across any number of selected spectators. Rectangular and empty column lists, singular matrices and zero determinants are included without division.

`MathCommons.ForestMinorAlternation.gateMinorBoth` now proves the actual both-occupied gate formula for arbitrary selected-row length: the minor of the left-multiplied matrix equals the block determinant times the original minor. It assumes the existing finite `Edge` certificate, duplicate-free selected rows, and membership of both endpoints. It needs no unitarity, forest condition or unproved alternation premise.

The required checker audits all 23 new declarations after their dependencies, preserving the earlier pair-history and ordered-field cube checks. The subsequent module proves the Boolean-filter single-exchange translation and sparse gate composition at arbitrary order. Generic arbitrary-matrix Cauchy–Binet tuple grouping and the arbitrary-complex matrix interface remain separate. [Exact statements, proof outline and source hashes](../research/marcus-de-oliveira/forest-minor-alternation/README.md)


## Actual forest product minors at every subset size

[ForestSubsetExchange.lean](ForestSubsetExchange.lean) closes the remaining arbitrary-order **gate-product** bridge. For either gate orientation and either endpoint order in the vertex list, its single-occupied row-replacement lemmas identify the canonical Boolean-filter minor with exactly the parity of the selected intervening spectators. Specializing the derived two-coefficient update to the actual identity matrix identifies its coefficients with actual gate minors. `gateMinorPredecessors` then proves sparse Cauchy–Binet over one or two predecessor subsets, with no minor-size, unitarity or forest premise.

| Lean declaration | Exact scope |
| --- | --- |
| `productEndpointAmplitude` | Every actual gate-product minor equals the canonical signed endpoint amplitude, for arbitrary selected subsets and product length; no unitarity or forest assumption. |
| `productMagnitude` | Under `UnitaryBlock` for used gates and `EveryEdgeIsBridge`, the squared actual product minor equals the canonical endpoint probability at every subset size. |
| `productEmptyFiber` | An empty canonical endpoint fiber gives an actual zero product minor, without unitarity or a forest assumption. |
| `productMassOne` | Squared actual forest product minors sum to one over distinct reachable endpoint subsets. |

The required checker audits all 25 declarations and preserves all earlier 261 reports and the cube exact-type/hash audit. Scalars are Gaussian rationals with the existing actual ordered Laplace determinant and matrix product. The complete finite `Edge` representation has two distinct vertices even when the used prefix is empty, so zero- and one-vertex ambient types are not represented by this API. Empty and full selected subsets, unequal selected cardinalities and vanishing transition factors are covered. Reachable endpoints are deduplicated using classical function equality, so the normalized law never counts an endpoint once per path.

The common law on full vertex permutations is now proved in the following module. Generic arbitrary-matrix all-order Cauchy–Binet, the standard arbitrary-complex matrix interface, and the determinant-of-a-sum/MDO convex-hull bridge remain separate. [Precise statement, proof, domain and audit hashes](../research/marcus-de-oliveira/forest-subset-exchange/README.md)


## One common permutation law for every squared forest minor

[ForestPermutationLaw.lean](ForestPermutationLaw.lean) constructs a single finite weighted list of genuine vertex permutations, with proved two-sided inverses. The list depends only on the gate schedule, its off-diagonal entries and prefix length, never on selected subsets. At each gate it branches into staying with weight factor $1-|b_k|^2$ or left-composing the endpoint transposition with factor $|b_k|^2$. Repeated permutation images and zero weights are retained as separate contributions.

`MathCommons.ForestPermutationLaw.commonSquaredMinorLaw` proves nonnegative weights, total weight one, and then quantifies over **all** starting and final Boolean subsets: the squared actual product minor equals the sum of the same weights for permutations carrying the start subset to the finish. `VertexPermutation` supplies actual forward/inverse maps; `subsetImageForward` fixes their convention. `commonPushforward` proves equality with canonical history expectations for every rational test function, including aggregation of both coin choices when they give the same subset. The common-law or matrix/history identity is not an assumed premise.

The required checker audits all 15 new reports alongside the earlier 286 and preserves the cube exact-type/hash audit. The actual squared-minor result uses the existing Gaussian-rational scalar domain, complete finite `Edge` representation, unitary used blocks and indexed-bridge condition. Push-forward and normalization alone need neither a forest nor unitarity; positivity follows from the probability bounds. The formal determinant-expansion identity carrying this common-minor law to the normal-matrix sum and its permutation-product convex hull remains separate, as do the arbitrary-complex library and continuity/density bridges. [Exact quantified statement, proof and audit](../research/marcus-de-oliveira/forest-permutation-law/README.md)
