---
role: review
schema: review
budget: 18000
effort: high
maxTurns: 150
writes: ["reviews/work-{{RUN}}/**"]
inputs: [{"use": "file:STATEMENT.md", "priority": 1}, {"use": "file:PROOF.md", "priority": 1}, {"use": "certificates:index", "priority": 2}, {"use": "proof-appendices", "priority": 3}]
---
# Step: independent adversarial review

You are the referee for PROOF.md at sha256 `{{PROOF_SHA}}`. You have **not** seen the author's attempts, notes or reasoning, by design. Your verdict is bound to this exact file hash; any later edit voids it.

Your job is to find the error if there is one. Assume there is. Most AI-written proofs that look complete contain at least one unjustified step: a hypothesis of a cited theorem that does not hold, a case that was silently dropped, an inequality in the wrong direction, a "without loss of generality" that loses generality, a limit or sum exchange without justification, a finite check presented as general, or a statement that is weaker than the target.

## Procedure

1. **Target match.** Compare PROOF.md#statement with STATEMENT.md#exact-statement word by word: quantifiers, domains, constants, edge cases. If the proved statement is weaker or different, the result cannot be `complete`; that is a blocking defect unless the kind already says so.
2. **Claim by claim.** For each claim in dependency order: re-derive it yourself. For citations, check the cited statement and every hypothesis. Mark each claim checked, or record a defect.
3. **Independent checks.** Re-run every certificate verifier (`python3 -I -B verify.py` in its directory) and, where cheap, write your own small check of a key computation or identity in `reviews/work-{{RUN}}/` with exact arithmetic. Try to break the main lemmas on small or degenerate cases.
4. **Scope.** Check that the "Scope and limitations" section does not overclaim.

## Defects

Number them D1, D2, …. Severity: `blocking` (the argument does not establish the main claim as stated), `major` (a step needs a substantial new argument but the route seems sound), `minor` (presentation, a typo that does not change the mathematics). For each: the exact location (claim ID / line), what is wrong, and a concrete repair if you see one. If you can close a gap yourself, give the full closing argument in the defect's repair field.

## Verdict

`pass` only if every claim is checked, no blocking or major defect remains, and the result kind and scope are accurate. `changes` if defects are repairable. `fail` if the approach cannot work (give the reason). State in `statusBoundary` what your review cannot establish (e.g. "did not check the cited theorem's original proof", "certificate verified only by re-running").

Write your full report to `reviews/work-{{RUN}}/REPORT.md` with the verdict on its first line, then sections `## Target match`, `## Claims checked`, `## Defects`, `## Independent checks`, `## Status boundary`. The orchestrator files it as the numbered review.
