# Changing the discovery harness

These rules are for editing `harness/` itself. To run campaigns, use the `discovery-loop` skill and [README.md](README.md).

- **Gates are the contract.** A new failure mode becomes a check in `lib/gates.mjs` (or `lib/lean.mjs`) with a test in `test/`, in the same change. Never weaken a gate to make a campaign pass.
- **Facts come from files, not from model answers.** Hashes, run IDs, identities, provenance and labels are computed by the orchestrator. Do not add schema fields that let a model assert them.
- **Schemas stay strict.** Every object is closed and every property is required (use `["string", "null"]` for optional fields). Codex rejects anything else. `test/units.test.mjs` enforces this.
- **Prompts keep their cache order.** Preamble, then role, then campaign constants, then capsule, then pack, then step. Keep volatile text out of `_preamble.md` and the role bodies.
- **Runtime files are generated.** Edit `runtime/`, then run `node harness/bin/discover.mjs install-runtime`. CI fails if `.claude/`, `.agents/` or `.codex/` drift from it.
- **Before pushing:**
  - Run `node --test harness/test/*.test.mjs`.
  - For any change to `lib/lean.mjs`, the Lean job in `.github/workflows/harness.yml` must pass on the real toolchain. A skipped Lean test is not a pass.
- **No live model runs in tests.** Fake executables (`test/fixtures/fake-runtime.mjs`, `fake-lake.mjs`) speak the real CLI protocols. Ask the operator before any live smoke run.
