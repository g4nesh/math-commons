# Runtimes, routing and operations

<!-- digest -->
Roles run as fresh Claude Code (`claude -p`) or Codex (`codex exec`) sessions. Each run has:
- a strict JSON output schema;
- the operator's own permissions (no bypass flags);
- a write scope that the orchestrator enforces afterwards;
- usage receipts from the existing Math Commons collectors.

A `manual` runtime turns any step into a prompt file, so an interactive Claude or Codex session (or a person) can do it and record the answer. Routing sends review, statement audit and certificate verification to the runtime that did not author the work. Skills and subagents extend both runtimes' own harnesses for interactive use.
<!-- /digest -->

## Drivers (`lib/runtimes.mjs`)

**Claude Code.** Invoked as:

```
claude -p --output-format stream-json --verbose --json-schema <schema>
  --session-id <uuid> --no-session-persistence --max-turns N
  --permission-mode <config> [--model M] [--effort E] [--allowedTools …] [--max-budget-usd X]
```

- The prompt goes on stdin.
- The structured answer is `result.structured_output`, parsed by `adapters/claude-usage.mjs`.
- Usage comes from `modelUsage`, which includes subagents.

**Codex.** Invoked as:

```
codex exec --sandbox workspace-write --color never --json --ephemeral
  --output-schema <schema> --output-last-message <file> [-m M] [-p profile]
  [-c model_reasoning_effort="E"] [extraArgs…] -
```

- The answer is read from the `-o` file.
- Usage comes from `turn.completed`, parsed by `adapters/codex-usage.mjs`.

**Manual.** Writes `.runs/<id>/prompt.md` and `schema.json`, then marks the campaign `pending`. Do the step, write `.runs/<id>/answer.json`, and run `discover record SLUG RUN-ID`. Pass `--as <identity>` so independence can be checked, for example `claude-interactive` or `human:g4nesh`.

**Schemas.** All are strict, which Codex requires: every object is closed and every property is required (nullable where optional). A unit test enforces this.

### After every run

The orchestrator does four things:
1. Diffs the campaign tree. Changes outside the role's `writes:` patterns (other than orchestrator-managed files) mark the run `scope-violation`, and the run is not ingested.
2. Ingests the answer. Facts it can establish itself, such as hashes, run IDs, identities and generator/verifier provenance, are never taken from the model.
3. Appends to `runs.jsonl`.
4. Rebuilds the capsule.

Three consecutive failed runs stop `run`.

### Permissions

Neither driver bypasses permissions.
- **Claude:** `permissionMode` defaults to `default`, and `allowedTools` lists exactly what roles may use. Read, Write and Edit are allowed; Web is allowed for screening; Bash is limited to python3, lake, read-only git and `discover`. Edit the list in `discovery.config.json`.
- **Codex:** runs in `workspace-write`. Web search for the screen role needs your Codex configuration, for example a profile passed through `runtimes.codex.profile` or `extraArgs`.

## Routing (`config.routing`)

| Role | Default | Why |
|---|---|---|
| scout, screen | claude | web research and source reading |
| map, draft, repair, formalize-statement, formalize | author (codex) | one consistent author identity per campaign |
| explore | codex, claude (alternating by direction) | diversity of approach |
| adversary | codex | search and computation |
| verify-certificate | reviewer, never the generator's runtime | independent check |
| review | reviewer, then author (alternating) | at least one review from another identity, plus a fresh-context second review |
| audit-statement | reviewer | the formalizer must not audit its own statement |

Swap `author` and `reviewer` to make Claude the author. Set any role to `manual` to do it interactively. One-off override: `discover step SLUG --runtime claude`.

## Configuration

Defaults live in `harness/config.json`. Put overrides in `discovery.config.json` at the repository root, or pass `--config FILE`; it is deep-merged. Environment overrides:

| Variable | Overrides |
|---|---|
| `DISCOVER_CLAUDE_BIN`, `DISCOVER_CODEX_BIN` | runtime binaries |
| `DISCOVER_LAKE_BIN` | lake |
| `DISCOVER_ROOT` | repository root |

Main settings:

| Setting | Controls |
|---|---|
| `runtimes.*.model`, `effort`, `maxTurns`, `timeoutMinutes` | how each runtime is invoked |
| `review.passes` (2), `review.requireDifferentRuntime`, `review.maxRepairRounds` | review |
| `explore.parallel`, `explore.maxRounds`, `explore.minDirections` | exploration |
| `formal.required`, `formal.allowedAxioms`, `formal.replay` (`null`, `"auto"` or a command), `formal.maxRounds` | the Lean stage |
| `budgets.campaignTokens`, `context.capsuleTokens`, `context.maxArtifactBytes` | budgets and context limits |

Both CLIs on this machine ship inside their desktop apps. If they are not on `PATH`, point the config at them. `discover doctor` shows what it found.

## Extending Claude Code and Codex for interactive work

`discover install-runtime` copies `harness/runtime/` into the locations each tool reads. `--check` (run in CI) fails on drift.

| Kind | Claude Code | Codex |
|---|---|---|
| Skills (`discovery-loop`, `lean-gate`, `context-economy`, `math-review`, `target-screening`) | `.claude/skills/*/SKILL.md` | `.agents/skills/*/SKILL.md` |
| Subagents (`math-referee`, `counterexample-hunter`, `certificate-verifier`, `lean-statement-auditor`, `literature-screener`) | `.claude/agents/*.md` | `.codex/agents/*.toml` |
| Campaign-local rules | `campaigns/*/CLAUDE.md` (`@AGENTS.md`) | `campaigns/*/AGENTS.md` |

Skills load on demand by description. Subagents run when invoked, for example "use the math-referee agent on campaigns/x/PROOF.md". A useful interactive pattern is to do the author roles yourself through `--runtime manual`, and to give reviews to the `math-referee` subagent or the other CLI. The subagent records its review through `discover record` with its own identity.

## Scheduled continuation

solved-conjectures ran an hourly continuation pointed at NEXT.md. The equivalent here is a scheduled `discover run SLUG --max-steps N`. Every step starts from the capsule, so the schedule needs no prompt of its own:

```sh
# cron, launchd, a Claude Code /loop or a Codex automation
node harness/bin/discover.mjs run SLUG --max-steps 6 && node harness/bin/discover.mjs reconcile --write
```

Bound each run with `--max-steps` and the campaign token budget. Run `reconcile` at least every four hours: it flags manual runs pending for more than 4 h, reviewed-but-unpublished campaigns, stale reviews, changed statements, non-current receipts, over-budget campaigns and oversized files.

## Math Commons coordination

The harness never writes the ledger. On transitions it appends to `outbox.jsonl` (stage changes, statement pins). To publish a checkpoint for a campaign tied to a ledger task, use the normal CLI from an up-to-date `main` checkout, keeping details short because the ledger has a size cap:

```sh
node bin/commons.mjs checkpoint --agent ID --task TASK --summary "…" --details "…" --evidence <commit URL>
```

Bounded ledger workers (`commons work … -- node adapters/claude.mjs`) and the harness can coexist. Use the worker for a ledger task with fixed paths, and the harness for a campaign that needs the gated pipeline. A campaign's directory is its owned path.
