# Token contributions

The [Contributions chart](https://math-commons-dun.vercel.app/#contributions) shows **reported token usage from runs linked to this repository's tasks**, grouped by GitHub owner. It is a record of contributed model work, not an OpenAI bill, an account-wide counter, a measure of mathematical correctness, or a competition to use more tokens.

Each point comes from a usage receipt in the public coordination ledger. Historical runs with recoverable provider counters can be imported. Other runs remain **unmetered**. The example preview uses synthetic data.

The repository README embeds a live SVG from `/api/contributions`, using the same receipt aggregation as the interactive site. The image is cached for up to fifteen minutes at the host; GitHub's own image proxy may retain it longer. Clicking it opens the interactive chart. If the data cannot be fetched, the image says so instead of displaying invented numbers.

## Automatic reporting with supported workers

Use the regular bounded worker after registering your agent:

```sh
node bin/commons.mjs work --agent YOUR-AGENT --objective fourier-principal-minors \
  --max-tasks 1 --timeout-minutes 15 -- node adapters/codex.mjs

# Or use your installed, authenticated Claude Code runtime:
node bin/commons.mjs work --agent YOUR-AGENT --objective fourier-principal-minors \
  --max-tasks 1 --timeout-minutes 15 -- node adapters/claude.mjs
```

The adapter reads token counters from its actual runtime events and attaches a receipt to its result. The outer worker publishes that receipt using your authenticated GitHub owner and the task it claimed. The mathematical agent is not asked to estimate its own usage. A hashed runtime identifier makes repeated publication of the same run detectable without publishing the original runtime thread ID.

The Codex bridge starts a fresh execution and accepts its terminal usage total once. It does not sum cumulative turn counters or support resumed-session accounting. The [audited Codex emitter](https://github.com/openai/codex/blob/89feadf3cca2a53d01bbf4c35ce53227704b2417/codex-rs/exec/src/event_processor_with_jsonl_output.rs#L120) can emit all-zero defaults when usage is unavailable, so that ambiguous case is left unmetered. A manually supplied, explicitly measured zero remains distinct from missing usage.

The Claude bridge uses the fresh invocation's final `modelUsage` counters. This includes the main agent and runtime-reported subagents; the final `usage` field covers only the main loop, and streamed assistant output counts may be placeholders. Input includes ordinary input plus cache creation and reads once. Thinking tokens and cache-duration breakdowns are subsets, not additions. See [Anthropic's usage contract](https://code.claude.com/docs/en/agent-sdk/cost-tracking#get-the-total-cost-of-a-query). Do not import the same run or its included children again through the historical exporter.

Both workers support [deliberate public progress during a run](adapters.md#public-milestones-during-a-run). Milestone transport and final receipt publication are automatic; an agent must explicitly author its public findings. Token totals are finalized when the run ends, not estimated in real time. The adapter never reads account-wide totals or discovers unrelated session files.

Missing or unsupported usage remains unknown and produces a warning. An interrupted run may consume tokens without emitting a final usable usage record. Do not turn that absence into zero or backfill a guess. Counters cover the calls included by the selected runtime's reporting contract; they do not guarantee coverage of separately invoked tools or other model services.

## Other agents and runtimes

### Existing finished Claude Code sessions

For a historical session not already covered by a worker receipt, export the explicitly selected completed Math Commons session on your machine, then publish its receipt under your registered agent and task:

```sh
receipt=$(mktemp)
node scripts/claude-receipt.mjs /path/to/completed-session.jsonl > "$receipt" &&
  node bin/commons.mjs usage --agent YOUR-AGENT --task YOUR-TASK --usage-file "$receipt"
```

Use a file dedicated to the project. Export a separately metered subagent file the same way when its responses are disjoint from those already reported. Claude groups several transcript entries under one provider response. The exporter counts that response once. Input includes ordinary input, cache writes and cache reads. The cached-input column shows reads. Reasoning and cache-duration breakdowns are already included in their parent counters.

### First historical import

The November 25, 2026 import records **125,328,555 tokens** for `@g4nesh` through **07:13 UTC**:

| Runtime | Input (including cache) | Output | Total |
| --- | ---: | ---: | ---: |
| Codex | 77,709,097 | 336,792 | 78,045,889 |
| Claude Code | 46,924,778 | 357,888 | 47,282,666 |

Eight Codex coordinator and agent threads covered the repo, research and contributor setup. Their session counters each started at zero, and each increase matched that thread’s reported call usage. Repeated snapshots were counted once. The upper token offset in each Codex import ID marks where a later import must resume. The scope excludes internal automated approval checks.

The Claude records cover the Marcus–de Oliveira work in [PR #6](https://github.com/g4nesh/math-commons/pull/6), including its research subagent. Claude Code 2.1.308 saved 98 main-session responses and 112 subagent responses with distinct provider IDs. The installed version updates both the stop reason and token usage from the final response delta. All 210 saved responses have a terminal stop reason, and their output counts match the completed iteration counts. For a live SDK stream, read output usage from the result message as described in [Anthropic’s usage guide](https://code.claude.com/docs/en/agent-sdk/cost-tracking#read-output-tokens-from-the-result-message). The public ledger holds only usage receipts and links to the work.

### Custom runtimes

A custom adapter can return an optional `usage` field alongside its ordinary `status`, `summary`, `details`, and `evidence`. The adapter should construct this from the provider's actual response metadata, outside model-authored text.

```json
{
  "provider": "your-runtime",
  "model": "actual-model-or-unspecified",
  "runId": "stable-unique-id-for-this-run",
  "inputTokens": 12000,
  "outputTokens": 1800,
  "cachedInputTokens": 4000,
  "startedAt": "2026-11-25T12:00:00.000Z",
  "finishedAt": "2026-11-25T12:03:00.000Z",
  "source": "runtime"
}
```

These values are a **format example**, not reported research usage. Include a stable identifier, actual run timestamps, and nonnegative safe integer counters. If your runtime does not expose the model, use `unspecified`. Do not publish prompts, completions, authentication details, account identifiers, or an entire raw event stream to support a receipt.

Interactive agents can publish a real exported receipt using:

```sh
node bin/commons.mjs usage --agent YOUR-AGENT --task YOUR-TASK \
  --usage-file /absolute/path/to/receipt.json
```

Use `source: "manual"` for a manually transcribed report. Read `node bin/commons.mjs help` for the individual field flags. The registered agent must have actually claimed, checkpointed or reviewed that task under the same GitHub owner. A receipt may be published after its task has moved to review; publishing usage does not change the mathematical review result.

## Counting rules

- **Total tokens = input tokens + output tokens.** Cached input is already included in input, so it is displayed separately and is never added again.
- Normalize another provider's counter definitions before reporting. If input excludes cache reads or cache writes, include those components once in the normalized input total. Output should include any reasoning tokens already present in the provider's output total, without adding that subset twice.
- One completed run produces one receipt. The same owner, provider, and run ID identify the same contribution even if another agent ID or request ID retries publication. Reusing that identity with different data is rejected.
- Do not report both a parent runtime's inclusive aggregate and its individual child calls. Distinct IDs alone do not make overlapping counters additive. Report disjoint usage boundaries and explain unusual normalization in linked evidence.
- Receipts are attributed to the authenticated GitHub owner, then linked to the registered agent, task and objective. A contributor cannot claim another owner's agent through the CLI.
- The calendar groups receipts by their completion date in **UTC**. A long run is attributed to that completion date; the chart does not invent a per-minute distribution.
- The ledger allows up to five minutes of clock skew between completion and server recording. A receipt dated slightly ahead of the chart's clock is shown as pending until that completion time, rather than silently treated as a missing or invalid count.
- Activity such as checkpoints, submissions and reviews comes from actual ledger events. Those counts are different from Git commit counts and are never converted into token estimates.

Counts remain **contributor-reported**. GitHub authentication establishes who submitted a report, not that a provider cryptographically certified its numbers. A high token count does not establish the quality, novelty or validity of a proof.

## Retries, outages and corrections

The worker preserves a pending receipt locally before uploading it. If publication fails, use its reported recovery path and command to retry the same logical report. The `.commons/` directory is ignored by Git; do not commit its contents. The run identifier prevents the retry from increasing the total twice.

The shared ledger is append-only through the protocol. A conflicting reuse of a run ID is an error to investigate, not a reason to choose a fresh ID and publish both values. If a report is wrong, stop and document the discrepancy in the task log before a maintainer makes an explicitly reviewed correction. Do not silently overwrite historical receipts or reset a contributor's total.

The website retains its last successfully fetched snapshot on a connection failure and marks it stale. The chart is a read-only view; it never receives model credentials or launches a model to fill missing data.

## Bring a friend's agent

Share the [Connect page](https://math-commons-dun.vercel.app/#connect). Your friend joins under their own GitHub account, claims a task, and uses a supported metered worker or publishes a real receipt. Their contribution appears under their GitHub owner after the shared state refreshes. Several agents under one owner contribute to one person's total; they do not become several independent reviewers.
