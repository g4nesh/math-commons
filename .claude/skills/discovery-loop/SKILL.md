---
name: discovery-loop
description: Drive a Math Commons discovery campaign (screen → map → explore → draft → review → formalize → publish) with harness/bin/discover.mjs, either autonomously with claude/codex runtimes or step by step inside this interactive session. Use when asked to start, continue, resume or check a campaign on an open conjecture, or to run "the harness".
---

# Discovery loop

The harness owns state, gates and routing; you supply mathematics one step at a time. Never hand-edit `campaign.json`, `CAPSULE.md`, ledgers, receipts or review records.

## Orient (cheap)

```sh
node harness/bin/discover.mjs status            # all campaigns
node harness/bin/discover.mjs status SLUG       # gate checks + the next action
cat campaigns/SLUG/CAPSULE.md                    # the only state you need to start
```

Expand a locator only when the next decision needs it: `discover show campaigns/SLUG/PROOF.md#c3`, `discover outline FILE`.

## Run

- Autonomous, bounded: `node harness/bin/discover.mjs run SLUG --until review --max-steps 20`. It stops at a stage, a parked/rejected campaign, a pending manual run, three failed runs or the token budget.
- One step: `discover step SLUG`. See what it would do first: `discover next SLUG`.
- Inside this session (no subprocess model): `discover step SLUG --runtime manual --as claude-interactive`, read `campaigns/SLUG/.runs/RUN/prompt.md`, do exactly that step (write only the listed paths), write the JSON answer to `.runs/RUN/answer.json` matching `.runs/RUN/schema.json`, then `discover record SLUG RUN`.
  - Do not review your own draft in the same session: independence is checked by identity. Use a subagent (math-referee) or the other runtime for reviews.

## New campaign

`discover new SLUG --title "..." --source URL` then fill `STATEMENT.md` (exact statement, source with statement number, conventions, scope) before the screen step. Use `discover scout --area "..."` to find candidates; the screen gate rejects prior resolutions and low-impact targets.

## When something stops

- `parked`: read the last history line in the capsule; fix the cause, then `discover set SLUG stage STAGE --reason "..."`.
- Gate failures name the file and fix: `discover gate SLUG`.
- Before a new target: `discover reconcile` and finish any reviewed-but-unpublished campaign first.
