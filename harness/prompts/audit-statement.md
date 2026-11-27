---
role: audit-statement
schema: statement-audit
budget: 9000
effort: high
maxTurns: 80
writes: ["STATEMENT_REVIEW.md", "formal/.harness-statement-tests/**"]
inputs: [{"use": "file:STATEMENT.md", "priority": 1}, {"use": "lean:statement", "priority": 1}, {"use": "section:STATEMENT_REVIEW.md#formalization-notes", "priority": 2}]
---
# Step: audit the Lean statement against the source

Decide whether `formal/**/Statement.lean` says exactly what STATEMENT.md#exact-statement says. You did not write it. Its fingerprint will be pinned on your verdict, and every later proof will be checked against it, so an error here silently invalidates the whole formal result.

## Procedure

1. Build a table in `STATEMENT_REVIEW.md` under `## Audit`: `| Source condition or quantity | Literal Lean representation | Match | Note |` with one row for every quantifier, domain, hypothesis, constant, convention and conclusion. Match is `exact`, `stronger`, `weaker` or `different`.
2. Check every project definition the statement uses (the context pack lists them): does it mean what its name suggests on all inputs, including 0, empty sets, degenerate cases?
3. Vacuity tests: argue (or check with small `example`s or `#eval` in `formal/.harness-statement-tests/`, which the gate never imports) that the hypotheses are satisfiable and that the conclusion is not trivially true or false for a wrong reason. For a known small case of the source, confirm the Lean statement gives the same answer.
4. Polarity: if the campaign disproves the statement, `¬ Claim` must be the negation of the source statement, not of a variant.

## Verdict

`faithful` only if every row is `exact` (or a documented, harmless notational difference). Otherwise `needs-change` with the precise fix, or `unfaithful` if the formal statement targets a different problem. First line of `## Audit`: `**Verdict:** <verdict>`.
