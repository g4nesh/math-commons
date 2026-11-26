# Reversible snapshot storage

The public ledger remains one JSON file on `coordination`. Every task,
checkpoint, event, review, usage receipt, ID, timestamp and operation fingerprint
is retained. This encoding removes only a second copy of checkpoint detail
text from its matching event. The complete text remains inside the same JSON,
in `tasks[taskId].checkpoints[checkpointId].details` (look up arrays by ID).

The deployed older website accepts at most 950,000 stored bytes. That number is
a compatibility measurement, not a new restriction on research publishing.
Canonical clients retain the existing 3,000,000-byte capacity and modern readers
can fetch larger immutable GitHub blobs. Compaction puts the current ledger
below the old website's limit; future growth can exceed it again. Upgrade the
deployed reader when deployment becomes available. No content is deleted and
no additional write cap is imposed to accommodate the older deployment.

## Encoding and exact reconstruction

Only an event of kind `checkpoint` with the same ID and task as a checkpoint,
and byte-identical UTF-8 detail text, is eligible. A reference is used only when
its complete JSON representation is smaller. For example:

```json
{
  "id": "checkpoint-1",
  "kind": "checkpoint",
  "taskId": "task-1",
  "detail": "Full detail: task task-1, checkpoint checkpoint-1.",
  "detailReference": { "version": 1, "sha256": "SHA256_OF_THE_FULL_UTF8_DETAIL" }
}
```

All other event fields remain unchanged. `detailReference` is reserved storage
metadata. The reader checks its exact schema, pointer, target identity and
SHA-256 before replacing `detail` with the full checkpoint text and removing
only the marker. Missing, ambiguous, changed or unsupported references fail
closed with an error. Pointer-looking text without a marker is ordinary text.
Checkpoints are append-only under the protocol; editing a referenced checkpoint
by hand is not a supported history operation.

`decodeSnapshot` restores references for both GitHub Contents and immutable
Blob reads. `GitHubStore.read` uses this decoder; every CLI read, worker packet,
watch, context, retry and mutation therefore operates on full original text.
`GitHubStore.write` calls `serializeSnapshot`, which restores any existing
references first and then applies deterministic compaction. Both the state API
and the contributions SVG use the same loader. Build-time example JSON is a
separate explicitly labeled seed and contains no live history.

Decoding an encoded snapshot reproduces byte-identical `JSON.stringify` output
to the original unencoded state, including record/property order. JSON layout
whitespace is a storage concern, not historical content. The raw pre-change
bytes are also preserved in the backup below.

## Older readers and writers

The older website receives valid JSON. Its activity detail displays the explicit
task/checkpoint pointer; its task inspector still has the entire checkpoint
text. It does not pretend the pointer is the original prose.

Official older writers parse JSON and clone complete records. They preserve
the marker, target checkpoint, and pointer unchanged when appending new work;
the next modern read reconstructs the same historical text. Their newly added
checkpoint events can retain full detail until the next modern write deduplicates
them. The pre-reference compact serializer likewise preserves unknown fields.

The much older `cdda559` CLI additionally checks a fully indented 950,000-byte
representation before any operation. The archived real ledger still exceeds
that older pretty-print limit after compaction, so that CLI refuses the write
with `STATE_CAPACITY`. This is safe refusal, not compatibility for publishing.
Update the CLI from `main` to publish. No workaround should strip markers or
checkpoint fields. Custom writers that discard unknown fields are unsupported.

SHA-conditional GitHub writes, request fingerprints and uncertain-write
reconciliation are unchanged. A retry reads and verifies references before
looking up the original operation ID. Malformed references never trigger an
empty/example fallback or a destructive repair.

## Immutable backup and reproducible audit

Before encoding, revision 607 was pinned at
[`3ad48348f8ebd5fdb9293c587c5fa53df97b4f75`](https://github.com/g4nesh/math-commons/blob/3ad48348f8ebd5fdb9293c587c5fa53df97b4f75/coordination/state.json).
Its exact 1,031,639 raw bytes are also in
[`revision-607-pre-compat.json.gz`](../archives/coordination/revision-607-pre-compat.json.gz),
with [hashes and provenance](../archives/coordination/revision-607-pre-compat.json).
The raw SHA-256 is
`543ae087b831cd9541759f94b804027ab9ed8f78f197adc15a5dd1058f643620`.

Run from the repository root:

```sh
node scripts/audit-snapshot.mjs archives/coordination/revision-607-pre-compat.json.gz
node --test test/snapshot.test.mjs test/transport.test.mjs test/api.test.mjs
```

The audit is read-only. For revision 607 it replaces 227 duplicate event details,
produces 901,121 stored bytes, and leaves 48,879 bytes below the deployed limit.
Before/after semantic SHA-256 is exactly
`1e2d75c9aa30e74041291266e5fa6b2c3cc77ef20fcf754e0f4b92338465d8e6`.
All 3 objectives, 78 tasks, 16 agents, 607 events, 50 usage receipts,
230 checkpoints and 13 reviews remain. Three short duplicate details remain
inline because referencing them would be larger. Tests check full semantic
equality, not merely those counts.

The same command accepts any local JSON snapshot or `.json.gz` backup and
reports its exact hashes, storage size, references and remaining legacy-reader
headroom (possibly negative). A false `legacyCompatible` result does not mean
modern research writes are disabled. Check a fresh live snapshot before the
rollout; these archived figures do not describe later revisions automatically.
For a rollback investigation, decompress the archived bytes locally or use
`restoreSnapshot` on a current snapshot. Do not reset `coordination` to the
backup: that would discard work published after revision 607. Any storage
rollback must retain current semantic state and use the normal checked writer;
fully expanding it may exceed the old production reader's limit.
