import test from 'node:test';
import assert from 'node:assert/strict';
import { mkdtemp, writeFile, copyFile, readFile, rm } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { spawnSync } from 'node:child_process';
import { ClaudeReceiptCollector, exportClaudeReceipt } from '../scripts/claude-receipt.mjs';

const response = (change = {}) => ({
  type: 'assistant', sessionId: 'synthetic-private-session', requestId: 'req_synthetic1',
  timestamp: '2026-11-24T01:00:00.000Z',
  message: { type: 'message', role: 'assistant', id: 'msg_synthetic1', model: 'claude-fixture',
    stop_reason: 'end_turn', content: [{ type: 'text', text: 'SECRET transcript 🧮' }],
    usage: { input_tokens: 11, cache_creation_input_tokens: 22, cache_read_input_tokens: 33,
      output_tokens: 44, cache_creation: { ephemeral_5m_input_tokens: 22 },
      iterations: [{ input_tokens: 11, cache_creation_input_tokens: 22, cache_read_input_tokens: 33, output_tokens: 44 }] } },
  ...change
});
function collect(records) {
  const collector = new ClaudeReceiptCollector();
  for (const record of records) collector.add(record);
  return collector.finish();
}
async function temporary(run) {
  const directory = await mkdtemp(join(tmpdir(), 'claude-receipt-test-'));
  try { return await run(directory); } finally { await rm(directory, { recursive: true, force: true }); }
}
const jsonl = records => records.map(record => JSON.stringify(record)).join('\n') + '\n';

test('normalizes top-level input/cache components once and deduplicates complete responses', () => {
  const first = response();
  const second = response({ requestId: 'req_synthetic2', timestamp: '2026-11-24T01:05:00Z',
    message: { ...first.message, id: 'msg_synthetic2', stop_reason: 'tool_use',
      usage: { input_tokens: 20, cache_creation_input_tokens: 10, cache_read_input_tokens: 5, output_tokens: 7 } } });
  const receipt = collect([
    { type: 'user', content: 'SECRET user text', usage: { input_tokens: 99999 } },
    first, first, { ...first, timestamp: '2026-11-24T01:00:01Z' }, second,
    { type: 'result', usage: { input_tokens: 99999, output_tokens: 99999 } }
  ]);
  assert.equal(receipt.provider, 'claude');
  assert.equal(receipt.source, 'runtime');
  assert.equal(receipt.model, 'claude-fixture');
  assert.equal(receipt.inputTokens, 101);
  assert.equal(receipt.outputTokens, 51);
  assert.equal(receipt.cachedInputTokens, 38);
  assert.equal(receipt.totalTokens, 152);
  assert.equal(receipt.startedAt, '2026-11-24T01:00:00.000Z');
  assert.equal(receipt.finishedAt, '2026-11-24T01:05:00.000Z');
  assert.match(receipt.runId, /^claude:[a-f0-9]{64}$/);
  assert.doesNotMatch(JSON.stringify(receipt), /SECRET|synthetic-private-session|req_synthetic|msg_synthetic/);
});

test('stable identity distinguishes separately selected subagent responses sharing a session', () => {
  const parent = response();
  const child = response({ requestId: 'req_child', message: { ...parent.message, id: 'msg_child' } });
  assert.notEqual(collect([parent]).runId, collect([child]).runId);
  assert.equal(collect([parent]).runId, collect([parent, parent]).runId);
  const later = response({ requestId: 'req_later', timestamp: '2026-11-24T02:00:00Z',
    message: { ...parent.message, id: 'msg_later' } });
  assert.equal(collect([parent]).runId, collect([parent, later]).runId);
});

test('inconsistent duplicates, reused IDs, or multiple sessions fail closed', () => {
  const first = response();
  for (const next of [
    response({ message: { ...first.message, usage: { ...first.message.usage, output_tokens: 45 } } }),
    response({ message: { ...first.message, model: 'different-model' } }),
    response({ message: { ...first.message, id: 'msg_different' } }),
    response({ requestId: 'req_different' }),
    response({ sessionId: 'different-session' })
  ]) assert.throws(() => collect([first, next]), /inconsistent|Multiple sessions/);
});

test('incomplete provider records, SDK envelopes, and malformed counters are rejected', () => {
  const first = response();
  for (const next of [
    response({ requestId: undefined }),
    response({ sessionId: undefined, session_id: 'sdk-only-session' }),
    response({ message: { ...first.message, stop_reason: null } }),
    response({ message: { ...first.message, usage: undefined } }),
    response({ message: { ...first.message, usage: { ...first.message.usage, output_tokens: '44' } } }),
    response({ message: { ...first.message, usage: { ...first.message.usage, cache_read_input_tokens: undefined } } }),
    response({ message: { ...first.message, usage: { ...first.message.usage, input_tokens: -1 } } }),
    response({ message: { ...first.message, usage: { ...first.message.usage, input_tokens: 1.5 } } }),
    response({ message: { ...first.message, usage: { ...first.message.usage, input_tokens: Number.MAX_SAFE_INTEGER } } })
  ]) assert.throws(() => collect([next]));
});

test('timestamps are validated; no usable provider response does not become a zero receipt', () => {
  for (const timestamp of ['not-a-date', '2026-02-31T00:00:00Z', '2026-11-24T25:00:00Z', '2026-11-24T01:00:00']) {
    assert.throws(() => collect([response({ timestamp })]), /timestamp/);
  }
  assert.throws(() => collect([{ type: 'user' }]), /No completed provider/);
  assert.throws(() => collect([response({ message: { model: '<synthetic>' } })]), /No completed provider/);
});

test('explicit measured zero is preserved, and multiple actual models are labeled mixed', () => {
  const first = response();
  const zero = response({ message: { ...first.message, usage: {
    input_tokens: 0, output_tokens: 0, cache_creation_input_tokens: 0, cache_read_input_tokens: 0 } } });
  assert.equal(collect([zero]).totalTokens, 0);
  const different = response({ requestId: 'req_other', message: { ...first.message, id: 'msg_other', model: 'other-fixture-model' } });
  assert.equal(collect([first, different]).model, 'mixed');
});

test('explicit file export is path-independent and leaves source bytes unchanged', async () => temporary(async directory => {
  const path = join(directory, 'one.jsonl'), copied = join(directory, 'copied.jsonl');
  const bytes = jsonl([response()]);
  await writeFile(path, bytes); await copyFile(path, copied);
  assert.deepEqual(await exportClaudeReceipt(path), await exportClaudeReceipt(copied));
  assert.equal(await readFile(path, 'utf8'), bytes);
}));

test('file parser handles fragmented Unicode and CRLF but rejects broken JSON and size limits', async () => temporary(async directory => {
  const path = join(directory, 'fixture.jsonl');
  const large = response();
  large.message.content = [{ type: 'text', text: '🧮'.repeat(20_000) }];
  await writeFile(path, jsonl([large]).replaceAll('\n', '\r\n'));
  assert.equal((await exportClaudeReceipt(path)).totalTokens, 110);
  await assert.rejects(exportClaudeReceipt(path, { maxLineBytes: 100 }), /size limit/);
  await assert.rejects(exportClaudeReceipt(path, { maxFileBytes: 100 }), /size limit/);
  await assert.rejects(exportClaudeReceipt(directory), /regular JSONL/);
  await writeFile(path, '{"SECRET":"unfinished');
  await assert.rejects(exportClaudeReceipt(path), /invalid or incomplete JSONL/);
}));

test('CLI prints only receipt JSON and sanitizes malformed input errors', async () => temporary(async directory => {
  const path = join(directory, 'fixture.jsonl');
  await writeFile(path, jsonl([response()]));
  const good = spawnSync(process.execPath, ['scripts/claude-receipt.mjs', path], { encoding: 'utf8' });
  assert.equal(good.status, 0, good.stderr);
  assert.equal(JSON.parse(good.stdout).totalTokens, 110);
  assert.equal(good.stderr, '');
  assert.doesNotMatch(good.stdout, /SECRET|synthetic-private-session|req_synthetic|msg_synthetic|fixture\.jsonl/);
  await writeFile(path, '{"SECRET": not-valid-json}\n');
  const bad = spawnSync(process.execPath, ['scripts/claude-receipt.mjs', path], { encoding: 'utf8' });
  assert.equal(bad.status, 1);
  assert.equal(bad.stdout, '');
  assert.doesNotMatch(bad.stderr, /SECRET|fixture\.jsonl/);
  const missing = spawnSync(process.execPath, ['scripts/claude-receipt.mjs'], { encoding: 'utf8' });
  assert.equal(missing.status, 1);
  assert.equal(missing.stdout, '');
}));
