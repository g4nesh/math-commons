import test from 'node:test';
import assert from 'node:assert/strict';
import { aggregateUsage, normalizeReceipts, ownerActivity, exampleUsage } from '../public/usage.js';
const now = Date.parse('2026-11-25T12:00:00Z');
const receipt = (overrides = {}) => ({ id: 'r1', owner: 'Alice', agentId: 'a1', runId: 'run1', provider: 'openai', model: 'model', inputTokens: 100, outputTokens: 20, cachedInputTokens: 50, startedAt: '2026-11-25T01:00:00Z', finishedAt: '2026-11-25T02:00:00Z', recordedAt: '2026-11-25T03:00:00Z', source: 'runtime', ...overrides });
test('missing receipts stay unknown, but a recorded zero remains a measured zero', () => {
  assert.equal(aggregateUsage(undefined, { now }).total, null);
  assert.equal(aggregateUsage([], { now }).days.at(-1).total, null);
  const data = aggregateUsage([receipt({ inputTokens: 0, outputTokens: 0, cachedInputTokens: 0 })], { now });
  assert.equal(data.total, 0); assert.equal(data.receiptCount, 1); assert.equal(data.days.at(-1).total, 0); assert.equal(data.days.at(-2).total, null);
});
test('cached input is a subset and never added twice', () => {
  const data = aggregateUsage([receipt()], { now });
  assert.equal(data.total, 120); assert.equal(data.input, 100); assert.equal(data.cachedInput, 50);
});
test('duplicate IDs and replayed owner/provider run identities count once', () => {
  const data = normalizeReceipts([receipt(), receipt(), receipt({ id: 'other-id', owner: 'alice' }), receipt({ id: 'r2', runId: 'run2' })], now);
  assert.equal(data.receipts.length, 2); assert.equal(data.duplicateCount, 2);
});
test('separate owners or providers may use the same run ID', () => {
  const data = aggregateUsage([receipt(), receipt({ id: 'r2', owner: 'Bob' }), receipt({ id: 'r3', provider: 'other' })], { now });
  assert.equal(data.total, 360); assert.deepEqual(data.availableOwners, ['alice', 'bob']);
});
test('UTC bucket follows finishedAt including offset and 7d includes today plus six days', () => {
  const data = aggregateUsage([receipt({ finishedAt: '2026-11-24T20:00:00-07:00', startedAt: '2026-11-24T19:00:00-07:00' }), receipt({ id: 'r2', runId: 'run2', startedAt: '2026-11-18T00:00:00Z', finishedAt: '2026-11-18T00:10:00Z' })], { period: '7d', now });
  assert.equal(data.startDay, '2026-11-19'); assert.equal(data.endDay, '2026-11-25'); assert.equal(data.days.length, 7); assert.equal(data.total, 120); assert.equal(data.days.at(-1).total, 120);
});
test('owner filter and cumulative values count only the selected window', () => {
  const data = aggregateUsage([receipt(), receipt({ id: 'r2', owner: 'bob' }), receipt({ id: 'r3', runId: 'old', startedAt: '2026-10-01T00:00:00Z', finishedAt: '2026-10-01T01:00:00Z' })], { owners: ['ALICE'], period: '7d', now });
  assert.equal(data.total, 120); assert.equal(data.days[0].cumulative, 0); assert.equal(data.days.at(-1).cumulative, 120); assert.equal(data.series.length, 1);
  assert.equal(aggregateUsage([receipt()], { owners: [], now }).total, null);
});
test('all time starts at first valid receipt and provenance counts are separate', () => {
  const data = aggregateUsage([receipt(), receipt({ id: 'r2', runId: 'run2', source: 'manual', startedAt: '2026-11-01T00:00:00Z', finishedAt: '2026-11-01T01:00:00Z' })], { period: 'all', now });
  assert.equal(data.startDay, '2026-11-01'); assert.equal(data.days.length, 25); assert.equal(data.runtimeCount, 1); assert.equal(data.manualCount, 1);
});
test('reject malformed dates, future receipts, bad counts and invalid cached subsets', () => {
  const bad = [receipt({ id: 'bad1', inputTokens: -1 }), receipt({ id: 'bad2', inputTokens: '100' }), receipt({ id: 'bad3', cachedInputTokens: 101 }), receipt({ id: 'bad4', finishedAt: '2026-02-30T00:00:00Z' }), receipt({ id: 'bad5', finishedAt: '2026-11-26T00:00:00Z' }), receipt({ id: 'bad6', outputTokens: Number.MAX_SAFE_INTEGER }), receipt({ id: 'bad7', recordedAt: 'bad' })];
  const data = normalizeReceipts(bad, now); assert.equal(data.receipts.length, 0); assert.equal(data.invalidCount, 7);
});
test('ledger activity stays actual counts and is never converted to tokens', () => {
  const rows = ownerActivity({ agents: [{ owner: 'Alice' }, { owner: 'alice' }], events: [{ actorOwner: 'ALICE', kind: 'checkpoint' }, { actorOwner: 'alice', kind: 'message' }] });
  assert.deepEqual(rows, [{ owner: 'alice', agents: 2, checkpoints: 1, messages: 1, reviews: 0, events: 2 }]); assert.equal('total' in rows[0], false);
});
test('synthetic example is deterministic, distinctly named, and has valid receipts', () => {
  const example = exampleUsage(now); assert.deepEqual(example, exampleUsage(now)); assert.ok(example.length > 90);
  assert.equal(normalizeReceipts(example, now).invalidCount, 0); assert.ok(example.every((item) => item.id.startsWith('synthetic-') && item.owner.endsWith('-example')));
});
test('accepted clock-skew receipts remain pending until their completion time', () => {
  const pending = receipt({ startedAt: '2026-11-25T11:59:00Z', finishedAt: '2026-11-25T12:04:00Z', recordedAt: '2026-11-25T12:00:00Z' });
  const before = aggregateUsage([pending], { now });
  assert.equal(before.pendingCount, 1); assert.equal(before.invalidCount, 0); assert.equal(before.total, null);
  const after = aggregateUsage([pending], { now: now + 240_000 });
  assert.equal(after.pendingCount, 0); assert.equal(after.total, 120);
  assert.equal(normalizeReceipts([receipt({ ...pending, finishedAt: '2026-11-25T12:06:00Z' })], now).invalidCount, 1);
});
test('hours 24 and malformed offsets are rejected rather than normalized', () => {
  for (const finishedAt of ['2026-11-24T24:00:00Z', '2026-11-24T23:00:00+25:00', '2026-11-24T23:00:00+00:99']) {
    assert.equal(normalizeReceipts([receipt({ startedAt: '2026-11-24T01:00:00Z', finishedAt })], now).invalidCount, 1);
  }
});
test('conflicting duplicates invalidate both records and do not choose first', () => {
  const result = aggregateUsage([receipt(), receipt({ id: 'retry', outputTokens: 30 }), receipt({ id: 'r2', runId: 'run2', owner: 'Bob' })], { now });
  assert.equal(result.invalidCount, 2); assert.equal(result.total, 120); assert.deepEqual(result.availableOwners, ['bob']);
});
test('aggregate precision overflow fails explicitly instead of rounding a total', () => {
  const large = receipt({ inputTokens: Number.MAX_SAFE_INTEGER - 20, outputTokens: 20, cachedInputTokens: 0 });
  assert.throws(() => aggregateUsage([large, receipt({ id: 'r2', runId: 'run2' })], { now }), /safe integer/);
});
