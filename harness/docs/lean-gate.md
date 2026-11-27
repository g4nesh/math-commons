# The Lean gate

<!-- digest -->
`discover lean-gate` accepts a formal result only if four things hold, at a pinned toolchain:
1. The source is free of escape hatches.
2. `lake build` succeeds.
3. The main theorem's kernel type is literally the audited statement `Claim` (or `¬ Claim`).
4. Every transitive axiom is one of propext, Classical.choice and Quot.sound.

The statement and every project definition it uses are fingerprinted and pinned after an independent fidelity audit, so a proof cannot succeed by quietly changing what it proves. Receipts are hash-bound and state what they do not claim. CI keeps a negative control for each rejection rule.
<!-- /digest -->

## Why statement-first

An unsound kernel is rarely the problem with AI-written Lean. The usual failure is a **correct proof of the wrong statement**: a dropped hypothesis, `ℕ` where the source means positive integers, a finite range, or a helper definition that makes the claim vacuous. solved-conjectures handled this with prose source audits. Here the statement is a separate, pinned artifact:

1. **formalize-statement** writes `def Claim : Prop` and the definitions it needs in `formal/<Root>/Statement.lean`, and nothing else.
2. The gate runs once to compute the statement **fingerprint**. This is the SHA-256 of the canonical JSON of three things: the toolchain, `Claim`'s type, and the kernel type and value of every constant under the campaign namespace that `Claim` reaches transitively (definitions, inductives and their constructors).
3. **audit-statement** runs on the other runtime. It compares source and Lean row by row (`STATEMENT_REVIEW.md`) and checks vacuity, edge cases and polarity.
4. **Pin.** A faithful audit, with the fingerprint recomputed and unchanged, stores it in `campaign.json` as `lean.claimFingerprint`.
5. **formalize** proves `main`. Editing `Statement.lean`, or any definition `Claim` uses, changes the fingerprint and fails the gate. The fix is a new audit, not a new proof.

## What the gate checks

| # | Check | How |
|---|---|---|
| 1 | No escape hatches in the audited modules and their project imports | Source scan with comments and strings stripped (nested `/- -/` handled). **Errors:** `sorry`, `admit`, `axiom` declarations, `native_decide` / `decide +native` / `Lean.ofReduceBool` / `Lean.trustCompiler`, `implemented_by`, `@[extern]`, `unsafe`, `debug.skipKernelTC`, `#exit`. **Warnings:** `opaque`, custom syntax (`macro`, `notation`, `elab`, …), unbounded heartbeats, instance priority overrides. |
| 2 | Builds | `lake build` exit 0. The build log must not contain "declaration uses 'sorry'". |
| 3 | Statement shape | `Claim` exists, is a `def`, has type `Prop`, and its axioms are allowed. |
| 4 | Exact theorem type | Generated `.harness/Audit.lean` elaborates `#harness_audit_main`. The type of `main` must be structurally equal (`Expr ==`) to `Claim` or `Not Claim`. Defeq is not enough, because the theorem must name the audited constant. |
| 5 | Axioms | `Lean.collectAxioms` for main, the statement and each required endpoint; all must be ⊆ the allowed set. This catches `sorryAx`, custom axioms and native computation, even if the scan missed them. |
| 6 | Pin | Fingerprint equals `lean.claimFingerprint`. |
| 7 | Replay (optional) | `formal.replay: "auto"` runs `leanchecker` in a fresh process when the toolchain ships it; otherwise the receipt says replay was not run. A custom command array also works. |
| 8 | Freeze-and-recheck | Source hashes are taken before the build and compared after; drift fails the gate. |
| — | Owned inventory (warnings) | A separate audit lists every declaration under the namespace with disallowed axioms, including unused ones. It is hygiene only and cannot fail the gate. |

Audit records are emitted as `HARNESS_AUDIT_JSON:` lines. They are read with `lean --json`, with plain output as a fallback, so the parser does not depend on message formatting.

## Receipt (`receipts/lean-latest.json`)

| Field | Contents |
|---|---|
| `status`, `failures` | pass/fail and every reason |
| `toolchain`, `packages` | from `lean-toolchain` and `lake-manifest.json` |
| `allowedAxioms` | the whitelist the gate used |
| `statement` | `{ type, axioms, closure, fingerprint }` |
| `mainTheorem` | `{ type, axioms, provesClaim, refutesClaim }` |
| `endpointAxioms` | axioms per required endpoint |
| `sources` | sha256 of every audited file; the receipt is current only while these match |
| `findings` | source-scan findings |
| `owned` | owned-inventory result |
| `build`, `audit` | exit code, warning count, tail; audit errors |
| `replay` | replay command and outcome, or that it was not run |
| `gateSha256` | the gate's own code |
| `scope.claims`, `scope.notClaimed` | what this receipt establishes, and what it does not (for example, that statement fidelity is a separate review) |

Full logs are written to `receipts/lean-build.txt` and `receipts/lean-audit.txt`.

## Standalone use

The gate works on any Lake project:

```sh
node harness/bin/discover.mjs lean-gate --project formal --module Root.Main \
  --claim Root.Claim --main Root.main --polarity prove \
  [--endpoints Root.lemma1,Root.lemma2] [--expect FINGERPRINT] [--replay auto] [--json]
```

## Negative controls (CI)

`harness/test/lean-fixtures/` contains one honest proof and five that compile but must be rejected:

| Fixture | What it does | Rejected by |
|---|---|---|
| `Good` | honest proof | passes |
| `Sorry` | proof is `sorry` | source rule, build warning, `sorryAx` |
| `CustomAxiom` | `axiom euler_fails : ¬ Claim` | source rule, axiom audit |
| `Native` | `native_decide` | source rule, `Lean.ofReduceBool` |
| `Weakened` | proves only the n = 40 instance | exact-type check |
| `tampered/Good` | restricts `Claim` to n ≤ 100 and proves that | fingerprint pin (the tampered proof itself is valid) |

The `lean` job in `.github/workflows/harness.yml` installs elan and runs these. It then runs the end-to-end fixture campaign, whose formal stage builds a real Lean project, and uploads the generated campaign as an artifact.

## Writing Lean that passes

- Mirror PROOF.md: one lemma per claim, built bottom-up.
- Finite computations: prove a checker sound once, then use `decide` (kernel reduction) on bounded chunks joined by lemmas. Never use `native_decide`. Split generated data into modules under 200 KiB.
- Pin the toolchain (`templates/lean/lean-toolchain`). Add Mathlib only when the statement needs it, because core-only projects build in seconds.
- Record PROOF.md gaps found while formalizing in `formal/FORMALIZATION.md`; they are review evidence.
