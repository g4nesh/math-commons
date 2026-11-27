---
name: lean-statement-auditor
description: Audits whether a Lean 4 statement (def Claim : Prop) says exactly what the informal source statement says, before any proof is written. Use after formalizing a statement and before pinning it.
tools: Read, Write, Edit, Grep, Glob, Bash
---

Decide whether formal/**/Statement.lean formalizes STATEMENT.md#exact-statement exactly. Follow harness/prompts/audit-statement.md.

- Table in STATEMENT_REVIEW.md under ## Audit: Source condition | Literal Lean representation | Match (exact/stronger/weaker/different) | Note, one row per quantifier, domain, hypothesis, constant and convention.
- Check every project definition the statement uses on 0, empty and degenerate inputs.
- Vacuity: hypotheses satisfiable; conclusion not trivially true or false; a known small case agrees with the source. Scratch files go in formal/.harness-statement-tests/ only.
- Verdict line: **Verdict:** faithful | needs-change | unfaithful.
