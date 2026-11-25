# MDO-FOREST · When do rotation histories produce a certificate?

**Written family proof.** With arbitrary complex spectra, suppose the relative unitary is a product of two-coordinate gates on distinct edges of a forest, each used once, in any order. Independent random transpositions supply a single distribution matching every squared minor. The argument allows left and right monomial factors. Eleven exact forest fixtures in dimensions 1–6 checked 3,332 compound entries. The written proof awaits independent-owner review. [Full proof, conventions and reproduction](../baselines/README.md) · [fixture output](../baselines/results.json)

**Where the argument breaks.** Different paths interfere on the triangle fixture. The product weights predict 6121/15625 for a coefficient whose exact value is 11881/15625. Both values lie in the original segment. A cyclic extension must account for the interfering paths. Reusing an edge also invalidates the unique-crossing argument. [Triangle derivation](../baselines/README.md#4-a-cycle-breaks-these-weights-not-the-conjecture)

**Lean coverage.** At the indexed revision, `cutHistoryUnique` and `crossingDecisionUnique` concern abstract integer histories on one cut. A separate three-dimensional rational example constructs convex weights for one pair. Relating subset histories to actual exterior-power amplitudes and the matrix convex-hull statement remains a separate obligation. [Lean statements](../../../formal/ForestCut.lean) · [checker and exact example](../../../formal/README.md) · [follow-up work](work.md)

[Research overview](../overview.md)
