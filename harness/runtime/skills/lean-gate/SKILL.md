---
name: lean-gate
description: Formalize a campaign result in Lean 4 so it passes the discovery harness's Lean gate - statement-first with an audited, fingerprint-pinned Claim, no sorry/axioms/native_decide, exact theorem type. Use when writing or fixing Lean for a campaign, auditing a Lean statement, or interpreting a lean-gate receipt.
---

# Lean gate

The gate (`node harness/bin/discover.mjs lean-gate SLUG`) passes only if all hold:

1. No `sorry`, `admit`, `axiom`, `native_decide`/`decide +native`, `implemented_by`, `extern`, `unsafe`, `#exit`, `debug.skipKernelTC` in the audited modules and their project imports (comments ignored).
2. `lake build` succeeds, with no "declaration uses 'sorry'".
3. `Claim` is a `def … : Prop`; its transitive project definitions are fingerprinted.
4. `main : Claim` (prove) or `main : ¬ Claim` (disprove), literally.
5. Transitive axioms of `main`, `Claim` and required endpoints ⊆ {propext, Classical.choice, Quot.sound}.
6. The fingerprint equals the one pinned after the statement audit.
7. Optional kernel replay (`formal.replay: "auto"` uses `leanchecker` when the toolchain has it).

## Order of work

1. **Statement only**: `formal/<Root>/Statement.lean` with `def Claim : Prop`, minimal transparent definitions, comments naming any reading of an ambiguity. `Main.lean` just imports it. Build.
2. **Audit** by a different identity: table `Source condition | Lean representation | Match` in `STATEMENT_REVIEW.md`; vacuity checks; every helper definition on 0 and degenerate inputs.
3. **Pin**: the orchestrator pins automatically after a faithful audit (`discover pin-statement SLUG` by hand).
4. **Prove** bottom-up, one lemma per PROOF.md claim. Do not touch the statement files; any change to them or to a definition they use breaks the pin and sends the statement back to audit.

## Finite computation in Lean

Prove a checker sound once, then `decide` (kernel reduction) on bounded chunks joined by lemmas. Never `native_decide`. Split generated data into modules under 200 KiB.

## Reading a receipt

`campaigns/SLUG/receipts/lean-latest.json`: `failures` (what to fix), `mainTheorem.type`/`axioms`, `statement.closure` (definitions the claim depends on), `findings` (warnings: custom syntax, opaque, heartbeats), `scope.notClaimed`. Full logs: `receipts/lean-build.txt`, `receipts/lean-audit.txt`.
