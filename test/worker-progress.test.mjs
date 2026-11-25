import test from 'node:test';
import assert from 'node:assert/strict';
import { mkdtemp, writeFile, readFile, appendFile, rm } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { createWorkerProgress, validateMilestone } from '../lib/worker-progress.mjs';
import { publishMilestone, replayMilestone } from '../adapters/progress.mjs';

const frame = (id, summary = 'Checked a bounded case') => ({ version: 1, type: 'checkpoint', id, summary, details: 'Public derivation summary', evidence: ['https://example.test/proof'] });
async function fixture(run) { const directory = await mkdtemp(join(tmpdir(), 'commons-progress-test-')); try { return await run(directory); } finally { await rm(directory, { recursive: true, force: true }); } }
const pause = ms => new Promise(resolve => setTimeout(resolve, ms));
async function until(test) { for (let i = 0; i < 100; i++) { if (test()) return; await pause(10); } throw new Error('Fixture condition timed out'); }

test('milestones require deliberate bounded public fields and safe URLs', () => {
  assert.deepEqual(validateMilestone(frame('one')), frame('one'));
  for (const changed of [{ type: 'assistant' }, { transcript: 'private' }, { summary: '' }, { details: 'x'.repeat(4001) },
    { evidence: ['file:///secret'] }, { evidence: ['https://user:secret@example.test'] }, { id: '../unsafe' }]) assert.throws(() => validateMilestone({ ...frame('one'), ...changed }), /Invalid public milestone/);
});

test('progress publishes while a run is active; duplicates, bursts and final flush are bounded', async () => fixture(async directory => {
  const operations = [];
  const progress = await createWorkerProgress({ repository: 'owner/repo', agent: 'alice', task: 't', runId: 'run', directory, pollMs: 5, throttleMs: 10000, mutate: async op => { operations.push(op); } });
  await publishMilestone(frame('one'), progress.path); await until(() => operations.length === 1);
  await publishMilestone(frame('one'), progress.path);
  await publishMilestone(frame('two'), progress.path); await publishMilestone(frame('three'), progress.path);
  await pause(30); assert.equal(operations.length, 1);
  await progress.finish(); assert.equal(operations.length, 2); assert.equal(operations[1].summary, frame('three').summary);
  assert.match(operations[0].id, /^progress-[a-f0-9]{64}$/); assert.notEqual(operations[0].id, operations[1].id);
  await assert.rejects(publishMilestone(frame('late'), progress.path));
}));

test('fragmented UTF-8 frames are decoded and conflicting identity aborts without publishing raw input', async () => fixture(async directory => {
  const operations = [], warnings = [];
  const progress = await createWorkerProgress({ repository: 'owner/repo', agent: 'alice', task: 't', runId: 'run', directory, pollMs: 5, mutate: async op => { operations.push(op); }, onProgress: x => warnings.push(x) });
  const bytes = Buffer.from(JSON.stringify(frame('one', 'A public 🧮 check')) + '\n');
  for (let i = 0; i < bytes.length; i += 3) await appendFile(progress.path, bytes.subarray(i, i + 3));
  await until(() => operations.length === 1);
  await publishMilestone(frame('one', 'Different content'), progress.path);
  await until(() => warnings.some(x => x.kind === 'progress-warning'));
  await assert.rejects(progress.finish(), /publication stopped/); assert.equal(operations.length, 1);
  assert.doesNotMatch(JSON.stringify(warnings), /Different content/);
}));

test('uncertain checkpoint write retains its exact operation for deduplicated recovery', async () => fixture(async directory => {
  const operations = [];
  const progress = await createWorkerProgress({ repository: 'owner/repo', agent: 'alice', task: 't', runId: 'run', directory, pollMs: 5, mutate: async op => { operations.push(op); throw new Error('PRIVATE transport detail'); } });
  await publishMilestone(frame('one'), progress.path);
  await assert.rejects(progress.finish(), error => error.code === 'PROGRESS_FAILED' && !error.message.includes('PRIVATE'));
  const saved = JSON.parse(await readFile(progress.recoveryFile, 'utf8'));
  assert.deepEqual(saved.operation, operations[0]);
  await replayMilestone(progress.recoveryFile, { config: { repository: 'owner/repo' }, mutate: async op => { assert.deepEqual(op, operations[0]); return { requestId: op.id }; } });
  await assert.rejects(replayMilestone(progress.recoveryFile, { config: { repository: 'other/repo' } }), /Invalid milestone recovery/);
}));

test('malformed and oversized artifacts stop safely without any ledger write', async () => fixture(async directory => {
  for (const text of ['PRIVATE malformed\n', JSON.stringify({ ...frame('one'), privateReasoning: 'PRIVATE' }) + '\n', 'x'.repeat(32001)]) {
    let writes = 0;
    const progress = await createWorkerProgress({ repository: 'owner/repo', agent: 'alice', task: 't', runId: 'run', directory, mutate: async () => { writes++; } });
    await writeFile(progress.path, text);
    await assert.rejects(progress.finish(), /publication stopped/); assert.equal(writes, 0);
  }
}));

test('finish rereads frames appended while the previous EOF poll awaits a slow mutation', async () => fixture(async directory => {
  const operations = []; let release;
  const firstMutation = new Promise(resolve => { release = resolve; });
  const progress = await createWorkerProgress({ repository: 'owner/repo', agent: 'alice', task: 't', runId: 'run', directory, pollMs: 5, throttleMs: 10000,
    mutate: async op => { operations.push(op); if (operations.length === 1) await firstMutation; } });
  await publishMilestone(frame('one', 'First milestone'), progress.path);
  await until(() => operations.length === 1);
  await publishMilestone(frame('last', 'Final milestone'), progress.path);
  const finishing = progress.finish(); release(); await finishing;
  assert.deepEqual(operations.map(op => op.summary), ['First milestone', 'Final milestone']);
}));
