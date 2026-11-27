---
role: formalize-statement
schema: formal
budget: 9000
effort: high
maxTurns: 100
writes: ["formal/**", "STATEMENT_REVIEW.md"]
inputs: [{"use": "file:STATEMENT.md", "priority": 1}, {"use": "section:PROOF.md#statement", "priority": 1}, {"use": "section:MAP.md#objects-and-conventions", "priority": 2}, {"use": "lean:project", "priority": 2}]
---
# Step: formalize the statement only (no proof yet)

Write the Lean 4 statement of the target as a definition `def Claim : Prop := …` in `formal/<Root>/Statement.lean`, plus any definitions it needs in the same namespace. **Do not prove anything in this step.** A separate auditor (another model) will compare your statement with the source; after it passes, its fingerprint is pinned and the proof must prove exactly that `Claim` (or `¬ Claim` for a disproof).

## Why statement-first

The commonest way an AI formalization goes wrong is not an unsound proof but a proof of the wrong statement: a missing hypothesis, `ℕ` where the source means positive integers, `≤` for `<`, a finite range where the source quantifies over all `n`, a definition that makes the claim vacuous. Fixing the statement before proving removes the temptation to adjust it to whatever is provable.

## Rules

- Use the pinned toolchain in `formal/lean-toolchain` and the dependencies in `formal/lakefile.toml`. Prefer Mathlib definitions when the project depends on Mathlib; otherwise define minimal, transparent versions and keep them in the campaign namespace (the auditor sees every project definition the statement uses).
- Literal fidelity over elegance: mirror the source's quantifiers and domains. If the source has an ambiguity, write a comment naming the reading you chose and why.
- No `sorry`, `axiom`, `opaque`, `noncomputable` tricks, custom syntax or instance overrides in the statement files.
- `lake build` must succeed. Add `formal/<Root>/Main.lean` containing only `import <Root>.Statement` for now.
- Write `STATEMENT_REVIEW.md` draft section `## Formalization notes` (the auditor completes the rest): a table `| Source condition | Lean representation | Note |` with one row per hypothesis, quantifier and convention.

Return the Lean configuration: modules (e.g. `["Root.Main"]`), claim (`Root.Claim`), main (`Root.main`, the theorem the next step will prove), namespace (`Root`), and any required endpoint names.
