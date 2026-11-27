---
role: explore
schema: attempt
budget: 12000
effort: high
maxTurns: 120
writes: ["attempts/{{DIRECTION}}/**", "certificates/**", "experiments/**"]
inputs: [{"use": "file:STATEMENT.md", "priority": 1}, {"use": "direction", "priority": 1}, {"use": "ledger:direction", "priority": 2}, {"use": "section:MAP.md#objects-and-conventions", "priority": 2}, {"use": "section:MAP.md#smallest-cases", "priority": 3}, {"use": "section:MAP.md#known-tools", "priority": 3}, {"use": "outline:MAP.md", "priority": 4}]
---
# Step: push one direction as far as it honestly goes

You own exactly one direction (below). Earlier attempts on it are in the context pack; failed approaches across the campaign are under "Do not retry" in the capsule. Do not repeat them; if you must revisit one, say precisely what is different.

## Method

1. Start from the direction's falsification test. Run it first, before building theory on top of it.
2. Work in small, checkable claims. After each one, try to break it: plug in the smallest cases, degenerate cases and the corners MAP.md names.
3. Use computation to discover, then convert discoveries into proofs or exact certificates. Exploration may use floats and optimisation; anything a claim relies on must be re-derived exactly.
4. When a reformulation helps, prove the translation in both directions or label it one-way.
5. Stop when you reach the milestone, when the falsification test kills the approach, or when your budget is nearly spent. A precise obstruction ("the reduction fails when the matrix is singular; example: …") is a valuable outcome.

## Computed claims

If a claim rests on computation, create `certificates/<name>/` containing the exact data (`certificate.json`, rationals as "p/q" strings, integers as strings when large), the generator you used, and `manifest.json`:
`{"name", "claim", "description", "notClaimed", "generator": {"path", "command"}, "verifier": {"path": null, "command": null}}`. Do **not** write the verifier: a separate run writes it independently so that a bug in your generator cannot also hide in the check.

## Files

`attempts/{{DIRECTION}}/{{ATTEMPT}}.md` with a first-line status (`**Outcome:** … · **Status:** …`), a `<!-- digest -->` block, then: claim, method, what was established (with evidence level per item), the failure point or the remaining gap, and the next concrete action. Scripts go next to it or under `experiments/`.

## Outcome vocabulary

`proof-sketch` (complete idea, details to write), `counterexample` (explicit object violating the exact statement), `reduction` (target follows from a cleaner statement), `partial` (precise weaker theorem), `failed` (approach refuted; say why), `inconclusive`. Status: `promising` (ready to draft or to continue with a clear next step), `open`, `killed`, `resolved`.
