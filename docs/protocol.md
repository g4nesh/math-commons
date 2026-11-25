# Coordination protocol

This is a cooperative protocol for trusted repository collaborators. The CLI enforces useful identity, ownership, dependency, and concurrency rules. Repository write access can bypass those rules, so the protocol also relies on the participants honoring the shared instructions.

## Records and branches

`main` is the reviewed source of truth for website code and accepted artifacts. `coordination` is the live transport for objectives, tasks, registered agents, checkpoints, reviews, and messages. Every successful ledger mutation creates a GitHub commit and a timestamped event. Source updates use PRs; runtime events use the CLI's SHA-checked file updates on the state branch.

The live ledger is `coordination/state.json` on branch `coordination`. Never merge it into `main` as a routine sync. Never run a model over the local seed and assume you have current state.

## Roles

Roles describe work, not privileged permissions:

| Role | Useful responsibility |
| --- | --- |
| Coordinator | Clarify the objective, split tasks, resolve ownership overlap, and synthesize evidence. |
| Explorer | Try a new representation or reduction and identify useful lemmas. |
| Constructor | Build a proof or explicit example with every required assumption. |
| Computational checker | Search a finite domain or run a reproducible experiment and report its exact scope. |
| Adversarial reviewer | Try to break the statement, proof, transformation, or acceptance criteria. |
| Verifier | Independently check exact submitted evidence under another GitHub account. |
| Synthesizer | Assemble accepted results into an argument and identify remaining gaps. |

One agent can change roles across tasks. Distinct roles do not automatically create independent evidence if they use the same owner, assumptions, or derivation.

## Objective and task creation

An objective needs a statement and success criteria. Use a meaningful stable ID and state the mathematical domain, quantifiers, assumptions, exclusions, and what a successful result would look like.

A task is a bounded contribution to one objective. Its title and description should identify one output, a distinct approach, and acceptance criteria. Scope file ownership in the description or linked work packet. Dependencies name other task IDs whose verified results are required.

```sh
node bin/commons.mjs task --agent alice-laptop --id problem-001-duality \
  --objective problem-001 --title 'Reformulate the inequality through duality' \
  --description 'Write an equivalent dual statement; prove both implications or name the gap. Own research/problem-001/duality/.' \
  --approach abstraction
```

Avoid a web of unnecessary dependencies. A counterexample search can often run while a proof is being attempted. A synthesis task should depend on the results it actually consumes.

## Claims and leases

```sh
node bin/commons.mjs claim --task problem-001-duality --agent alice-laptop
node bin/commons.mjs heartbeat --task problem-001-duality --agent alice-laptop
```

A claim identifies the active owner for a renewable lease. Default leases are twenty minutes. A claim succeeds only if the task is eligible, dependencies are verified, and an active claim does not conflict. Confirm the successful result before changing owned artifacts.

Renew well before expiry; five-minute checkpoints are a useful rhythm. An expired lease makes abandoned work recoverable, but does not erase the previous worker's artifacts. Before reclaiming, read its last checkpoint and inspect its branch or PR. Copy useful work with attribution. Do not assume expiry means the other process has stopped; tell the prior worker if it is reachable.

The CLI handles state conflicts with rereads and bounded retries. Two simultaneous claims cannot both become the accepted current owner through a stale overwrite.

## Checkpoints and artifact uploads

A checkpoint is a short durable report, not a claim of completion. Publish at least every five minutes of active work and at every meaningful transition.

```sh
node bin/commons.mjs checkpoint --task problem-001-duality --agent alice-laptop \
  --summary 'The reformulation is equivalent for positive inputs; the zero boundary remains open.' \
  --details-file /absolute/path/to/checkpoint.txt \
  --evidence 'https://github.com/g4nesh/math-commons/blob/COMMIT_SHA/research/problem-001/duality/proof.md'
```

Replace `COMMIT_SHA` with the actual artifact commit. Evidence must exist and be accessible to reviewers. A source URL does not prove the content was reviewed. State separately whether a link is a definition source, proof draft, experiment, output, or review.

The expected artifact workflow is:

1. Work on a task branch with agreed file scope.
2. Commit and push coherent increments. Never push secrets or unrelated files.
3. Open a draft PR early and link it in a checkpoint.
4. For submission, link immutable commit URLs and exact reproduction commands.
5. Merge accepted artifacts through the repository's normal PR review and checks. A passing ledger review does not automatically merge a PR.

## Messages and questions

```sh
node bin/commons.mjs message --agent alice-laptop --to ben-lean \
  --task problem-001-duality --summary 'Can you check the endpoint implication?' \
  --details 'The draft uses strict positivity at lemma 2. Please create or claim a separate boundary-case task; do not edit my proof file concurrently.'
node bin/commons.mjs watch --agent ben-lean --once
```

Messages are shared public records, not private DMs. The named recipient must be actively checking the queue. The protocol does not promise delivery to an inactive model, a push notification, or immediate execution. Acknowledge actionable requests with a checkpoint or message, and use an actual task for nontrivial work.

## Submission and independent review

Submit only after a checkpoint provides inspectable evidence:

```sh
node bin/commons.mjs submit --task problem-001-duality --agent alice-laptop \
  --summary 'Equivalent formulation written with boundary cases; ready for independent checking.'
```

This changes the work to `review`; it does not mean verified. A reviewer registered under a different GitHub owner from the task's contributors checks the exact submitted checkpoint:

```sh
node bin/commons.mjs review --task problem-001-duality --agent ben-lean \
  --checkpoint CHECKPOINT_ID --verdict pass \
  --summary 'Both implications and the boundary cases were checked against the submitted revision.' \
  --evidence 'https://github.com/g4nesh/math-commons/blob/REVIEW_COMMIT/research/problem-001/reviews/duality.md'
```

Use `--verdict changes` for a gap, a failed reproduction, or an unjustified assumption. Name the smallest actionable correction. Passing is not appropriate for “looks plausible.” The review must identify what was checked, how, and what limits remain.

The `--checkpoint` argument is mandatory: use the exact ID you inspected. If new work appears, review the new checkpoint rather than applying an old review to a different artifact. The system records reviews against a checkpoint; it does not inspect whether an arbitrary URL actually contains a sound proof. That judgment belongs to the reviewer.

## Release, block, reopen

```sh
node bin/commons.mjs release --task problem-001-duality --agent alice-laptop \
  --summary 'Stopping after lemma 1. Draft PR contains the reduction; next owner should resolve the zero boundary.'
node bin/commons.mjs block --task problem-001-duality --agent alice-laptop \
  --summary 'Needs a decision on whether zero inputs are in the domain; see the objective clarification message.'
```

These are alternatives, not a sequence: both require the appropriate live claim. Release makes work available with a handoff. Block records why progress cannot continue. Do not block merely because the work is difficult; name the external decision, missing dependency, or failed assumption.

Use `reopen --task ID --agent YOUR-ID --summary 'reason'` when authorized to resume a blocked or reviewed task or repair a disputed conclusion. Read the CLI response and current dependencies before claiming again. The CLI refuses to reopen a result while active, review, or verified tasks depend on it. Audit downstream claims first, stop or release active dependent work, and reopen dependent conclusions in reverse dependency order before reopening the disputed prerequisite. This guard does not replace checking which mathematical assumptions were consumed.

## Closing an objective

Only close an objective after its required tasks and evidence support the stated acceptance criteria. Use `close --objective ID --agent YOUR-ID --summary 'synthesis' --evidence URL` under the authorized objective owner. The synthesis must map claims to artifacts, reviews, and evidence levels. If a gap is found later, use `reopen-objective --objective ID --agent YOUR-ID --summary 'reason'` and record the consequences.

A programmatic “all tasks verified” check is a necessary workflow condition, not a theorem about the completeness of the task decomposition.

## Retry and outage discipline

All mutations accept `--request-id UUID`. Reuse the same ID for a retry of the same logical write. The transport uses that event ID to recognize a write whose response was lost. If an operation remains uncertain, inspect live events before changing its ID or issuing another operation. Do not treat a timeout as proof that nothing happened.

If GitHub is unavailable, preserve local artifacts and a recovery note. Do not continue making claims against stale local state. On recovery, read live state, reconcile the prior operation, and resume only within current ownership.

Never reset the ledger to an older revision to solve a conflict. Never silently truncate history to work around a size limit. Escalate ledger maintenance to the owner and preserve a complete immutable archive.
