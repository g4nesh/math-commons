---
role: formalize
schema: formal
budget: 16000
effort: high
maxTurns: 200
writes: ["formal/**"]
inputs: [{"use": "lean:statement", "priority": 1}, {"use": "lean:receipt", "priority": 1}, {"use": "file:PROOF.md", "priority": 2}, {"use": "lean:project", "priority": 3}]
---
# Step: prove the pinned statement in Lean

Prove `{{MAIN}} : {{GOAL}}` in `formal/` following PROOF.md. The statement files are frozen: their fingerprint is pinned, and changing `Statement.lean` or any definition it uses fails the gate.

## The gate you must pass (exactly what `discover lean-gate` checks)

1. No `sorry`, `admit`, `axiom`, `native_decide`, `decide +native`, `implemented_by`, `extern`, `unsafe`, `#exit` or `debug.skipKernelTC` in the audited modules or the project modules they import (comments are ignored).
2. `lake build` succeeds with no "declaration uses 'sorry'".
3. `{{MAIN}}` is a `theorem` whose type is literally `{{GOAL}}`.
4. Its transitive axioms are within {propext, Classical.choice, Quot.sound}.
5. The statement fingerprint equals the pinned one.

## Method

- Mirror the claim structure of PROOF.md: one lemma per claim, named after it (`c3_bound`), proved bottom-up. Build after each lemma.
- For finite computations, prefer a checker function proved sound once plus `decide` (kernel reduction) on bounded data, joined by lemmas, over one giant `decide`. Never `native_decide`.
- If a step of PROOF.md is wrong or underspecified, do not force it: record it in `formal/FORMALIZATION.md` under `## Gaps found` with the exact claim ID; this is valuable review evidence.
- Keep files under 200 KiB; split large generated data into modules.
- Run `node harness/bin/discover.mjs lean-gate {{SLUG}}` yourself before finishing and include its first failures in your notes if it fails.

The context pack shows the last gate failures, if any. Fix those first.
