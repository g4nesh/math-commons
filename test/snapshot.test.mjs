import test from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { gunzipSync } from 'node:zlib';
import { createHash } from 'node:crypto';
import { compactSnapshot, restoreSnapshot, serializeSnapshot, decodeSnapshot, MAX_STATE_BYTES, LEGACY_READER_MAX_BYTES } from '../lib/snapshot.mjs';
import { auditSnapshot } from '../scripts/audit-snapshot.mjs';

function fixture() {
  const details = 'Exact Unicode evidence: α² + β² = 1; 🧮.\n'.repeat(30);
  const checkpoint = { id: 'checkpoint-1', agentId: 'agent', summary: 'Evidence', details, evidence: ['https://example.com/proof'], createdAt: '2026-11-26T00:00:00Z' };
  const event = { id: checkpoint.id, kind: 'checkpoint', actor: 'agent', actorOwner: 'owner', taskId: 'task-1',
    objectiveId: 'objective', summary: checkpoint.summary, detail: details, createdAt: checkpoint.createdAt,
    links: checkpoint.evidence, requestFingerprint: 'original-fingerprint', futureField: { nested: ['keep', 1] } };
  return { schemaVersion: 1, revision: 12, updatedAt: checkpoint.createdAt, objectives: [],
    tasks: [{ id: 'task-1', checkpoints: [checkpoint], reviews: [{ checkpointId: checkpoint.id, evidence: ['unchanged'] }] }],
    agents: [], events: [event], usage: [{ id: 'usage-preserved', inputTokens: 14, outputTokens: 2 }], custom: { retain: true } };
}
const metadata = text => ({ type: 'file', sha: 'a'.repeat(40), encoding: 'base64', content: Buffer.from(text).toString('base64'), size: Buffer.byteLength(text) });

test('references preserve every original semantic byte and never mutate the input', async () => {
  const original = fixture(), before = JSON.stringify(original);
  const text = serializeSnapshot(original), stored = JSON.parse(text);
  assert.equal(JSON.stringify(original), before);
  assert.match(stored.events[0].detail, /task task-1, checkpoint checkpoint-1/);
  assert.equal(stored.tasks[0].checkpoints[0].details, original.events[0].detail);
  assert.equal(stored.events[0].detailReference.sha256, createHash('sha256').update(original.events[0].detail).digest('hex'));
  assert.equal(JSON.stringify(restoreSnapshot(stored)), before);
  assert.equal(JSON.stringify(await decodeSnapshot(metadata(text), () => assert.fail('Unexpected blob request'))), before);
  const file = { ...metadata(text), encoding: 'none', content: '' };
  assert.equal(JSON.stringify(await decodeSnapshot(file, async () => metadata(text))), before);
  assert.ok(Buffer.byteLength(text) < Buffer.byteLength(JSON.stringify(original)));
  assert.equal(serializeSnapshot(stored), text);
});

test('only identical checkpoint-event detail is replaced when the reference is smaller', () => {
  const state = fixture();
  state.events.push({ ...state.events[0], id: 'message-1', kind: 'message' });
  state.tasks[0].checkpoints.push({ ...state.tasks[0].checkpoints[0], id: 'short', details: 'short' });
  state.events.push({ ...state.events[0], id: 'short', detail: 'short' });
  state.tasks[0].checkpoints.push({ ...state.tasks[0].checkpoints[0], id: 'different' });
  state.events.push({ ...state.events[0], id: 'different', detail: 'Different text '.repeat(100) });
  const stored = compactSnapshot(state);
  assert.deepEqual(stored.events.slice(1), state.events.slice(1));
  assert.equal(JSON.stringify(restoreSnapshot(stored)), JSON.stringify(state));
  const unmarked = fixture();
  unmarked.events[0].detail = 'Full detail: task task-1, checkpoint checkpoint-1.';
  assert.equal(restoreSnapshot(unmarked), unmarked, 'Pointer-like user text is never interpreted without a marker.');
});

test('missing, ambiguous or changed reference targets fail closed instead of losing history', () => {
  const mutations = [
    state => { state.events[0].detailReference.version = 2; },
    state => { state.events[0].detailReference.sha256 = 'b'.repeat(64); },
    state => { state.events[0].detailReference.extra = true; },
    state => { state.events[0].kind = 'message'; },
    state => { state.events[0].detail = 'altered pointer'; },
    state => { state.events[0].taskId = 'missing'; },
    state => { state.tasks[0].checkpoints[0].details += 'changed'; },
    state => { state.tasks[0].checkpoints = []; },
    state => { state.tasks[0].checkpoints.push(state.tasks[0].checkpoints[0]); },
    state => { state.tasks.push(state.tasks[0]); },
  ];
  for (const mutate of mutations) {
    const state = compactSnapshot(fixture()); mutate(state);
    const before = JSON.stringify(state);
    assert.throws(() => restoreSnapshot(state), /reference|Duplicate/);
    assert.equal(JSON.stringify(state), before);
  }
});

test('a plain-JSON legacy read/append/write preserves exact old detail on modern rehydration', () => {
  const original = fixture();
  const legacy = JSON.parse(serializeSnapshot(original));
  const next = { id: 'legacy-message', kind: 'message', detail: 'New uncompressed content.' };
  legacy.events.push(next); legacy.revision++;
  const written = JSON.stringify(structuredClone(legacy), null, 2)+'\n';
  const modern = restoreSnapshot(JSON.parse(written));
  const expected = structuredClone(original); expected.events.push(next); expected.revision++;
  assert.equal(JSON.stringify(modern), JSON.stringify(expected));
  assert.equal(JSON.stringify(restoreSnapshot(JSON.parse(serializeSnapshot(modern)))), JSON.stringify(expected));
});

test('legacy size is only a measurement; canonical reads/writes retain the existing larger bound', async () => {
  const state = fixture(); state.extra = 'x'.repeat(LEGACY_READER_MAX_BYTES);
  const before = JSON.stringify(state);
  assert.ok(Buffer.byteLength(serializeSnapshot(state)) > LEGACY_READER_MAX_BYTES);
  assert.equal(JSON.stringify(state), before);
  assert.deepEqual(await decodeSnapshot(metadata(before)), state);
  const oversized = fixture(); oversized.extra = 'x'.repeat(MAX_STATE_BYTES);
  const oversizedBefore = JSON.stringify(oversized);
  assert.throws(() => serializeSnapshot(oversized), /3000000-byte snapshot limit/);
  assert.equal(JSON.stringify(oversized), oversizedBefore);
});

test('the immutable real ledger round-trips exactly below the old production loader limit', async () => {
  const directory = new URL('../archives/coordination/', import.meta.url);
  const bytes = gunzipSync(await readFile(new URL('revision-607-pre-compat.json.gz', directory)));
  const manifest = JSON.parse(await readFile(new URL('revision-607-pre-compat.json', directory), 'utf8'));
  const report = auditSnapshot(bytes);
  assert.equal(report.originalRawSha256, manifest.rawSha256);
  assert.equal(report.originalSemanticSha256, manifest.semanticSha256);
  assert.equal(report.restoredSemanticSha256, manifest.semanticSha256);
  assert.ok(report.storedBytes <= 950_000);
  assert.deepEqual(report.recordsAfter, report.recordsBefore);
  assert.equal(report.referenceCount, 227);
  assert.ok(report.legacyPrettyBytes > 950_000, 'The pre-cdda559 writer must refuse rather than rewrite this snapshot.');
});
