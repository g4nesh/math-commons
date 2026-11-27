---
role: map
schema: directions
budget: 9000
effort: high
maxTurns: 60
writes: ["MAP.md", "experiments/**"]
inputs: [{"use": "file:STATEMENT.md", "priority": 1}, {"use": "digest:SCREEN.md", "priority": 2}, {"use": "section:SCREEN.md#attacks", "priority": 2}]
---
# Step: map the problem and choose directions

Build the map that every later worker starts from, then choose a small portfolio of directions with **different failure modes**. Ten workers repeating one proof idea add less than a constructor, a counterexample hunter and a structural reformulator working from the same precise statement.

## MAP.md

Write these sections, each short and concrete:

1. `## Objects and conventions`: every object, its domain, normalisations, and the admissible transformations (symmetries, scalings, relabellings) that preserve the statement.
2. `## Smallest cases`: the smallest nontrivial instances. Work out at least one by hand or by exact computation (scripts in `experiments/`, exact arithmetic, fixed seeds) and record what it shows. If the statement is false, a counterexample is often small; if true, the small cases show which inequalities are tight.
3. `## Known tools`: theorems that apply, with their exact hypotheses checked against this problem, and the closest known partial results.
4. `## Lenses`: which of these fit and what each would produce: invariants/conservation, symmetry/quotient, change of representation, duality/certificates (LP/SDP/SOS), extremal/minimal counterexample, local-to-global, decomposition, induction/recursion, monotonicity/convexity, probabilistic method, topology/structure, transport from an analogous solved problem, adversarial computation. For each reformulation say what is retained, what is discarded, and where the translation can fail (degenerate cases, division by zero, limits, non-invertibility).
5. `## Where it should break`: the parameter corners and structures most likely to harbour a counterexample, and why.

## Directions

Choose 2–4 directions. At least one must be `adversarial` (search for counterexamples or stress the strongest special case) and at least one `constructive` or `abstraction`. For each give: a hypothesis it tests, a **falsification test** that is cheap and decisive, a **milestone** that would justify more budget, and a token budget (typically 200k–1.5M). Prefer directions whose failure would still teach something reusable (an obstruction, a tight example).
