---
role: adversary
schema: attempt
budget: 10000
effort: high
maxTurns: 120
writes: ["attempts/{{DIRECTION}}/**", "certificates/**", "experiments/**"]
inputs: [{"use": "file:STATEMENT.md", "priority": 1}, {"use": "direction", "priority": 1}, {"use": "ledger:direction", "priority": 2}, {"use": "section:MAP.md#where-it-should-break", "priority": 2}, {"use": "section:MAP.md#smallest-cases", "priority": 2}, {"use": "section:MAP.md#objects-and-conventions", "priority": 3}]
---
# Step: try to break the statement

Your job is to find a counterexample to the exact statement, or to show convincingly where one cannot be. Assume nothing about which answer is true.

## Search design

- Enumerate small cases exhaustively when feasible; state the exact domain covered (e.g. "all graphs on ≤ 9 vertices, 274,668 isomorphism classes, via nauty geng").
- Beyond exhaustive range, use structured families (extremal, symmetric, degenerate, near-equality cases from MAP.md), then local search / optimisation of the violation margin from many seeds. Record seeds, iteration counts and the best margin found.
- Optimise a *continuous* violation measure when one exists; sign tests alone hide near misses.
- Watch for false positives: floating-point error, wrong conventions (loops, multigraphs, 0-indexing), an off-by-one domain. Re-check every apparent counterexample in exact arithmetic before reporting it.

## When you find a violation

Minimise it (smallest size, simplest entries). Write the exact object to `certificates/<name>/certificate.json` with `manifest.json` (see the explore step for fields) stating precisely which condition of the statement fails and by how much. Do not write the verifier; an independent run does.

## When you find none

Report the covered domain exactly, the best near-miss and where it occurs (these are usually the tight cases a proof must handle), and any structural reason the search suggests the statement holds. Never write "verified for all n" for a finite search.

## Files

As in the explore step: `attempts/{{DIRECTION}}/{{ATTEMPT}}.md` with the outcome line, digest, search design, coverage, results and next action; scripts under `experiments/`.
