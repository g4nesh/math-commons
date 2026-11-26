import test from 'node:test';
import assert from 'node:assert/strict';
import { GitHubStore, GitHubError, validateConfig } from '../lib/github.mjs';
import { serializeSnapshot } from '../lib/snapshot.mjs';
import { applyOperation } from '../lib/core.mjs';

const now = '2026-11-25T00:00:00Z';
const config = { repository: 'owner/math-commons', stateBranch: 'coordination', statePath: 'coordination/state.json', leaseMinutes: 20 };
function initial() {
  let state = { schemaVersion: 1, revision: 0, updatedAt: now, objectives: [], tasks: [], agents: [], events: [] };
  const ops = [
    { id: 'j1', kind: 'join', agent: 'a', name: 'A' },
    { id: 'j2', kind: 'join', agent: 'b', name: 'B' },
    { id: 'o', kind: 'objective', agent: 'a', objectiveId: 'p', title: 'P', statement: 'Prove P', success: 'Reviewed' },
    { id: 't', kind: 'task', agent: 'a', objectiveId: 'p', taskId: 't', title: 'T', description: 'Prove T' },
  ];
  for (const op of ops) state = applyOperation(state, op, { owner: 'owner', now }).state;
  return state;
}
function server({ loseWriteResponse = false, failReadsAfterWrite = false, initialState = initial() } = {}) {
  let state = initialState, sha = 'sha0', writes = 0, conflicts = 0;
  const reply = (body, status = 200) => new Response(JSON.stringify(body), { status, headers: { date: now, 'content-type': 'application/json' } });
  return {
    get state() { return state; }, get writes() { return writes; }, get conflicts() { return conflicts; },
    fetch: async (url, options) => {
      assert.equal(options.headers.Authorization, 'Bearer fake-token');
      if (url.endsWith('/user')) return reply({ login: 'owner' });
      if (options.method === 'GET') {
        if (failReadsAfterWrite && writes) throw new Error('disconnected');
        return reply({ type: 'file', encoding: 'base64', sha, size: Buffer.byteLength(JSON.stringify(state, null, 2)), content: Buffer.from(JSON.stringify(state)).toString('base64') });
      }
      const body = JSON.parse(options.body);
      assert.equal(body.branch, 'coordination');
      if (body.sha !== sha) { conflicts++; return reply({ message: 'conflict' }, 409); }
      state = JSON.parse(Buffer.from(body.content, 'base64').toString()); writes++; sha = `sha${writes}`;
      if (loseWriteResponse) throw new Error('response lost');
      return reply({ commit: { sha: `commit-${writes}` } });
    },
  };
}
const storeFor = fake => new GitHubStore(config, { token: 'fake-token', fetchImpl: fake.fetch, sleep: async () => {}, maxRetries: 2 });

test('concurrent claims produce exactly one winner, never overwrite newer ownership', async () => {
  const fake = server(), a = storeFor(fake), b = storeFor(fake);
  const result = await Promise.allSettled([
    a.mutate({ id: 'claim-a', kind: 'claim', agent: 'a', task: 't' }),
    b.mutate({ id: 'claim-b', kind: 'claim', agent: 'b', task: 't' }),
  ]);
  assert.equal(result.filter(item => item.status === 'fulfilled').length, 1);
  assert.equal(result.find(item => item.status === 'rejected').reason.code, 'TASK_BUSY');
  assert.equal(fake.writes, 1); assert.equal(fake.conflicts, 1);
  assert.equal(fake.state.events.filter(event => event.kind === 'claim').length, 1);
});

test('response lost after successful commit reconciles by the original operation ID', async () => {
  const fake = server({ loseWriteResponse: true }), store = storeFor(fake);
  const op = { id: 'stable-request', kind: 'claim', agent: 'a', task: 't' };
  const result = await store.mutate(op);
  assert.equal(result.reconciled, true); assert.equal(result.duplicate, true); assert.equal(fake.writes, 1);
  const retry = await store.mutate(op);
  assert.equal(retry.duplicate, true); assert.equal(fake.writes, 1);
});

test('unreconciled ambiguous writes stop and include reusable request ID', async () => {
  const fake = server({ loseWriteResponse: true, failReadsAfterWrite: true });
  await assert.rejects(storeFor(fake).mutate({ id: 'retry-this', kind: 'claim', agent: 'a', task: 't' }), error => error.code === 'WRITE_UNCERTAIN' && error.requestId === 'retry-this' && /identical command/.test(error.message));
  assert.equal(fake.writes, 1);
});

test('authorization errors are not retried and never expose token or raw response', async () => {
  let calls = 0;
  const store = new GitHubStore(config, { token: 'super-secret', fetchImpl: async () => { calls++; return new Response(JSON.stringify({ message: 'super-secret' }), { status: 403 }); } });
  await assert.rejects(store.mutate({ kind: 'claim', agent: 'a', task: 't' }), error => error instanceof GitHubError && error.status === 403 && !error.message.includes('super-secret'));
  assert.equal(calls, 1);
});

test('read detects malformed or oversized state and config refuses source branch writes', async () => {
  assert.throws(() => validateConfig({ ...config, stateBranch: 'main' }), /separate branch/);
  assert.throws(() => validateConfig({ ...config, repository: '../evil' }), /owner\/repository/);
  const store = new GitHubStore(config, { fetchImpl: async () => new Response(JSON.stringify({ type: 'file', encoding: 'base64', sha: 'sha', size: 50, content: Buffer.from('{broken').toString('base64') }), { status: 200 }) });
  await assert.rejects(store.read(), /not valid JSON/);
});


test('optional configured default objective must be a nonempty objective slug', () => {
  assert.equal(validateConfig({ ...config, defaultObjectiveId: 'marcus-de-oliveira' }).defaultObjectiveId, 'marcus-de-oliveira');
  assert.equal(validateConfig(config).defaultObjectiveId, undefined);
  for (const defaultObjectiveId of ['', null, 12, '../other', 'Uppercase', 'a'.repeat(81)]) {
    assert.throws(() => validateConfig({ ...config, defaultObjectiveId }), /defaultObjectiveId/);
  }
});


test('ledger encoding retains every nested field and keeps individual records diffable', () => {
  const state = initial();
  const text = serializeSnapshot(state);
  assert.deepEqual(JSON.parse(text), state);
  assert.ok(text.length < JSON.stringify(state, null, 2).length);
  for (const event of state.events) assert.ok(text.split('\n').some(line => line.trim().replace(/,$/, '') === JSON.stringify(event)));
});

test('official read/mutate/write restores references and retains all prior content on retry', async () => {
  let state = initial();
  state = applyOperation(state, { id: 'claim', kind: 'claim', agent: 'a', task: 't' }, { owner: 'owner', now }).state;
  state = applyOperation(state, { id: 'evidence', kind: 'checkpoint', agent: 'a', task: 't',
    summary: 'An exact result', details: 'Full preserved evidence. '.repeat(100), evidence: ['https://example.com/proof'] }, { owner: 'owner', now }).state;
  const stored = JSON.parse(serializeSnapshot(state));
  assert.ok(stored.events.at(-1).detailReference);
  const fake = server({ initialState: stored, loseWriteResponse: true }), store = storeFor(fake);
  assert.deepEqual((await store.read()).state, state);
  const operation = { id: 'after-encoding', kind: 'message', agent: 'a', task: 't', to: 'all', summary: 'Still working', details: 'New content.' };
  const expected = applyOperation(state, operation, { owner: 'owner', now }).state;
  const result = await store.mutate(operation);
  assert.equal(result.reconciled, true);
  assert.equal(fake.writes, 1);
  assert.deepEqual((await store.read()).state, expected);
  assert.deepEqual(fake.state.events.find(event => event.id === 'evidence'), stored.events.find(event => event.id === 'evidence'));
  // The pre-reference protocol also clones whole records. Plain-JSON operations
  // retain the marker rather than converting a pointer into historical content.
  const legacy = applyOperation(stored, { ...operation, id: 'legacy-append' }, { owner: 'owner', now }).state;
  assert.deepEqual(legacy.events.find(event => event.id === 'evidence'), stored.events.find(event => event.id === 'evidence'));
});
