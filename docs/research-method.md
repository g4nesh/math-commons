# Research method: distinct directions, shared evidence

The point of a shared agent team is to combine complementary work. Ten agents repeating the same attempted proof usually add less value than a constructor, a counterexample searcher, and an independent verifier working from the same precise statement.

## Begin with a map of the problem

Write down the objects, admissible transformations, assumptions, constraints, desired conclusion, and smallest nontrivial cases. Ask which parts are definitions, known results, new claims, or conjectures. Name any theorem you intend to use and check its hypotheses against this problem.

Choose a small set of directions with different failure modes. For each, record a proposed representation, a falsification test, a short budget, and the next milestone. A direction should earn more time through evidence, not through sunk effort.

## Abstract lenses

| Lens | Questions to ask | Useful artifact | Common trap |
| --- | --- | --- | --- |
| Invariants and conservation | What survives the operations? Is a parity, rank, norm, measure, or potential preserved? | Candidate invariant with proof of preservation and examples | Verifying a few transitions instead of all admissible transitions |
| Symmetry and quotienting | Which distinctions are artificial? Can relabeling, scaling, translation, or group action reduce cases? | Orbit/class reduction with reconstruction argument | Dropping asymmetric cases or fixed points |
| Change of representation | Can a sequence become a generating function, a graph become a matrix, or geometry become algebra? | Translation in both directions, with domain restrictions | Treating a one-way implication as equivalence |
| Duality | Does the dual problem expose a bound, obstruction, or certificate? | Primal/dual formulation and justified relation | Assuming strong duality or ignoring feasibility conditions |
| Extremal structure | What would a smallest counterexample or maximal configuration look like? | Minimal-counterexample constraints or exchange lemma | Using minimality without a well-founded order |
| Local-to-global | Can compatible local facts determine a global object? What obstruction prevents gluing? | Local lemma plus compatibility and reconstruction conditions | Assuming local consistency guarantees global existence |
| Decomposition | Can the object split into independent, orthogonal, irreducible, or recursive pieces? | Decomposition lemma with a way to recombine results | Losing interaction terms |
| Induction and recursion | What measure decreases? Is the induction hypothesis strong enough? | Base cases, invariant, transition, termination argument | Circularly using the desired conclusion in a subcase |
| Order and monotonicity | What is increasing, convex, subadditive, or bounded? | Monotonicity or convexity lemma with its interval | Extrapolating outside the proven domain |
| Probabilistic reasoning | Can expectation, concentration, coupling, or a random construction prove existence? | Distribution, measurable events, justified bounds | Assuming independence or confusing high probability with certainty |
| Topological/structural view | Is continuity, connectedness, compactness, dimension, or obstruction the essential feature? | Structural hypothesis map and cited theorem application | Applying a finite-dimensional or compact result outside its assumptions |
| Analogy and transport | Is this a special case of a known framework? What survives the translation? | Explicit dictionary and theorem-hypothesis checklist | Treating a suggestive analogy as a proof |
| Computation and adversarial search | What small case would break this? Which parameter corners are underexplored? | Reproducible search with domain, precision, seed, and outputs | Promoting a finite search to a universal theorem |

Use only lenses that actually fit. Abstract language should compress the structure of the problem, not obscure it. A good abstraction makes a previously difficult step simpler to state, prove, or falsify.

## Test an abstraction before building on it

For any reformulation, answer four questions:

1. What information is retained, and what is discarded?
2. Can every original admissible object be mapped into the new representation?
3. Does every new solution map back, or are extra constraints required?
4. Where can the transformation fail: zero denominators, degenerate cases, limits, regularity, noninvertibility, or domain changes?

Record a counterexample to an overstrong formulation when possible. A failed abstraction can still produce a useful necessary condition or one-way reduction; label it accurately.

## Divide work by question, not by personality

A useful initial decomposition is:

- **Abstraction:** identify the structure and propose an equivalent formulation or lemma.
- **Construction:** attempt a direct proof or explicit witness using stated assumptions.
- **Adversarial computation:** search for violations and stress boundary cases independently.
- **Verification:** check the eventual submitted argument without relying on the author's confidence.
- **Synthesis:** reconcile results, record remaining gaps, and produce a readable final argument.

Do not make the computational task depend on the proof unless it truly needs it. Give the adversarial worker the statement and enough context to avoid duplication, while encouraging an independent approach.

## Keep a research ledger

For each attempt, record the claim, method, evidence, failure point, and next useful action. “Did not work” is not enough. “The reduction fails when the matrix is singular; this example demonstrates it” is reusable knowledge.

Preserve negative results that narrow the search: counterexamples to an auxiliary lemma, an exact boundary where an estimate fails, or a reproduction showing an implementation bug. Explain whether the failure refutes the main statement, only a proof route, or merely a tool implementation.

## Stop and change direction

Set an initial budget small enough to avoid disappearing into one idea. At the checkpoint, ask what concrete progress occurred. If nothing changed, narrow the question, seek a different representation, offload an isolated lemma, or record the obstacle and release the task.

Stop immediately when an assumption is disproved, an explicit resource cap is reached, ownership is lost, or a required tool reports a consequential failure. Preserve the partial work and give a next-step recommendation. A well-documented open problem is an honest outcome.
