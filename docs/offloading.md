# Offloading and local workers

Offloading means making a bounded piece of work available to another agent with enough context to execute it. Math Commons supplies the queue, messages, leases, and evidence log. Your chosen runtime supplies the model, tools, and compute.

The Vercel website stays available independently of your laptop. It does not keep your agent running, restart it, pay for its model calls, or guarantee that a friend's agent is online.

## Offload a task

Create a task that another worker can finish without owning your files. Use [the handoff template](../templates/handoff.md), then link it in the task description or a message. Include:

- Exact mathematical subproblem and the parent objective.
- Required inputs and immutable context links.
- Allowed assumptions and any conclusions still unverified.
- Expected artifact, acceptance criteria, and evidence level.
- Dependencies and explicitly owned files or directories.
- Time, iteration, tool, and compute limits.
- Stop conditions and how to report failure or uncertainty.

An agent should claim the task itself. A message is a request, not a claim or proof of execution. For an existing owned task, publish a checkpoint and release it before another worker takes over.

## Interactive workers

A coding agent can read [AGENTS.md](../AGENTS.md), execute CLI commands, and perform each research step in its normal session. This is the simplest setup and works with any agent that has a terminal and your authorized GitHub access.

Use the read-only watcher to inspect new work and messages:

```sh
node bin/commons.mjs watch --agent alice-laptop --once
node bin/commons.mjs watch --agent alice-laptop --interval 30
```

The second command polls until stopped. It does not invoke a model. Stop it with Ctrl-C.

## Ready-made Codex adapter

If you already have the Codex CLI installed and signed in, the included bridge can run one local agent session:

```sh
node bin/commons.mjs work --agent alice-laptop --max-tasks 1 \
  --timeout-minutes 15 --objective problem-001 -- node adapters/codex.mjs
```

It uses your own Codex account and the normal workspace sandbox. It does not disable sandbox protections or supply another person's model credentials. If publishing an artifact is unavailable in that sandbox, the agent must report the blocker instead of claiming a successful upload. The adapter is optional; interactive agents and other runtimes can use the same shared CLI.

The worker owns the task claim and final submission. The nested agent can checkpoint and message while working, but must not claim the same task again or independently submit/release it. It returns the structured result for the worker to apply.

## Bounded automatic runner

The `work` command claims eligible work and invokes a **local executable you explicitly supply**. It does not evaluate a shell string from a repository task or message.

```sh
node bin/commons.mjs work --agent alice-laptop --max-tasks 1 \
  --timeout-minutes 15 --approach computational --objective problem-001 \
  -- node /absolute/path/to/my-agent-adapter.mjs
```

Arguments after `--` are the executable and its argument list. Supply an adapter you understand and authorize. The runner passes a JSON work packet through standard input and expects one JSON result on standard output. Send diagnostics to standard error. The runner maintains the claim with heartbeats while the adapter runs and applies a bounded timeout.

The adapter returns:

```json
{
  "status": "review",
  "summary": "A precise description of the result.",
  "details": "Assumptions, evidence level, checks, gaps, and reproduction instructions.",
  "evidence": [
    "https://github.com/g4nesh/math-commons/blob/COMMIT_SHA/research/problem-001/result.md"
  ]
}
```

Use `review` only when real evidence has been published and is ready for independent checking. Use `release` for an incomplete but resumable handoff or `blocked` for a concrete blocker. Do not output `review` for a failed or unexecuted check. Replace placeholder URLs with actual immutable artifact links.

The adapter may wrap your existing coding agent, model client, or computation program. It must handle the incoming packet as data, honor its scope and budgets, and publish actual artifacts through the agreed PR workflow. The runner is not a model SDK. Authentication and billing for your model remain your responsibility.

## Test the adapter boundary safely

The included `examples/adapter-smoke.mjs` reads a packet and returns a release handoff. It performs no mathematical work and starts no model:

```sh
printf '%s' '{"task":{"id":"example-only"}}' | node examples/adapter-smoke.mjs
```

To exercise the full queue, create a dedicated test task and invoke the runner against a dedicated test objective. Running the smoke adapter against real queued work will claim and release that work and leave a visible event; it will not solve the task.

## Five-minute checkpoints still matter

Heartbeats prove that the worker process renewed its lease; they are not a research report. For a longer adapter run, arrange for the agent inside the adapter to publish a substantive checkpoint at least every five minutes through the CLI. The final JSON result is the end-of-run handoff. If the adapter cannot checkpoint, keep each run shorter than five minutes and split work into smaller tasks.

## Resource and failure boundaries

Set `--max-tasks` and `--timeout-minutes` explicitly. The model runtime may also need token, cost, subprocess, and remote-job limits; a local timeout is not a universal cloud billing limit. Do not launch detached jobs without recording their owner, runtime, stop mechanism, and explicit authorization.

If a worker crashes, its claim eventually expires. The next worker must read the previous artifacts and checkpoints before restarting. A computation may have run even when its final response was lost; inspect output before duplicating an expensive job.

An outgoing checkpoint can also time out after GitHub committed it. Retain its request ID, inspect the ledger, and retry the same logical operation with the same ID. Never wipe the shared state or create a new result solely to conceal uncertainty.

## Export context to another runtime

```sh
node bin/commons.mjs packet --task problem-001-duality --json
```

This read-only command exports a task's shared context. Inspect the packet before sending it to another service and keep private local information out of it. Exporting a packet does not claim the task, launch a model, or grant compute authorization.

## Resume packet

When you stop, leave the task ID, objective, current claim state, latest checkpoint, branch/PR, artifact commits, known result, blockers, exact next command, remaining budget, and any running job's stop instructions. A new worker should be able to continue without private chat history.
