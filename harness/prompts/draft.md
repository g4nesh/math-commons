---
role: draft
schema: draft
budget: 16000
effort: high
maxTurns: 120
writes: ["PROOF.md", "proof/**"]
inputs: [{"use": "file:STATEMENT.md", "priority": 1}, {"use": "promising-attempts", "priority": 1}, {"use": "certificates:index", "priority": 2}, {"use": "section:MAP.md#objects-and-conventions", "priority": 2}, {"use": "section:MAP.md#known-tools", "priority": 3}]
---
# Step: write the complete argument

Turn the promising attempts into one complete, checkable document. A referee who has never seen the attempts must be able to verify every step from PROOF.md, the cited sources and the certificates alone.

## PROOF.md structure (the gate parses it)

```
# <Title>
<!-- digest -->
≤10 lines: the exact result, its kind, the method in one sentence, what is not claimed.
<!-- /digest -->

**Result kind:** `complete` | `partial` | `obstruction` | `consequence` | `candidate`
**Main claim:** `C<n>`

## Statement
The exact theorem proved here. If it differs from STATEMENT.md#exact-statement in any way, say exactly how (and the kind cannot be `complete`).

## Claims
| ID | Statement | Status | Depends on | Evidence |
|---|---|---|---|---|
| C1 | … | proved / cited / computed / assumed / open | C… | source, or certificates/<name> |

## Proofs
### C1. <short name>
…complete proof, or the exact citation with the hypotheses checked…

## Scope and limitations
What is proved, what remains open, what the computation does and does not establish.

## Reproduction
Exact commands for every certificate.
```

## Rules

- Every row of the Claims table has a `### C<n>` section. Every dependency is a row. No cycles.
- `computed` claims cite `certificates/<name>` whose receipt passes; describe in words exactly which finite statement the certificate establishes and how it combines with the proved claims.
- `cited` claims quote the source statement (version, theorem number, page) and list each hypothesis with why it holds here.
- `assumed` claims are allowed only in a `partial` (conditional) result and must be named as hypotheses of the main theorem.
- Choose the result kind honestly; when in doubt, choose the weaker kind.
- Put long case analyses in `proof/<topic>.md` and link them from the claim's section; keep PROOF.md under 200 KiB.
