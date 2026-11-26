import test from 'node:test';
import assert from 'node:assert/strict';
import { mkdtemp, readFile, readdir, rm, stat, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { applyOperation, validateState } from '../lib/core.mjs';
import { normalizeUsage, parseTokenFlag, publishWorkerUsage, usageOperationId } from '../lib/usage.mjs';
import { usageOperationFromFlags } from '../bin/commons.mjs';
import { runAdapter, work } from '../lib/worker.mjs';

const NOW = '2026-11-25T06:00:00.000Z';
const receipt = (overrides = {}) => ({ provider: 'codex', model: 'unspecified', runId: 'codex:run-one', inputTokens: 120, outputTokens: 30, cachedInputTokens: 80, startedAt: '2026-11-25T05:58:00Z', finishedAt: '2026-11-25T05:59:00Z', source: 'runtime', ...overrides });
let sequence = 0;
function apply(state, op, owner = 'alice') {
  return applyOperation(state, { id: `op-${++sequence}`, agent: 'a', ...op }, { owner, now: NOW });
}
function fixture({ claim = true } = {}) {
  let state = { schemaVersion: 1, revision: 0, updatedAt: NOW, agents: [], tasks: [], objectives: [], events: [] };
  state = apply(state, { kind: 'join', name: 'A' }).state;
  state = apply(state, { kind: 'join', agent: 'a-alt', name: 'A alternative' }).state;
  state = apply(state, { kind: 'join', agent: 'b', name: 'B' }, 'bob').state;
  state = apply(state, { kind: 'objective', objectiveId: 'p', title: 'P', statement: 'Proof', success: 'Review' }).state;
  state = apply(state, { kind: 'task', taskId: 't', objectiveId: 'p', title: 'T', description: 'One contribution' }).state;
  if (claim) state = apply(state, { kind: 'claim', task: 't' }).state;
  return state;
}
function mockStore() {
  let state = fixture({ claim: false });
  return { config: { repository: 'owner/repo', leaseMinutes: 20 }, get state() { return state; }, read: async () => ({ state, serverTime: NOW }), mutate: async op => { const result = apply(state, op); state = result.state; return result; } };
}
async function temporary(run) {
  const directory = await mkdtemp(join(tmpdir(), 'commons-usage-test-'));
  try { await run(directory); } finally { await rm(directory, { recursive: true, force: true }); }
}

test('usage is backward compatible and total excludes cached input', () => {
  const before = fixture();
  assert.equal(before.usage, undefined); validateState(before);
  const result = apply(before, { kind: 'usage', task: 't', ...receipt(), owner: 'mallory' });
  assert.equal(result.state.usage.length, 1);
  assert.equal(result.usage.totalTokens, 150);
  assert.equal(result.usage.cachedInputTokens, 80);
  assert.equal(result.usage.owner, 'alice');
  assert.equal(result.usage.agentId, 'a');
  assert.equal(result.usage.objectiveId, 'p');
  assert.equal(before.usage, undefined);
});

test('usage accepts post-submit and expired claims but rejects noncontributors and impersonation', () => {
  let state = fixture();
  state = apply(state, { kind: 'checkpoint', task: 't', summary: 'Proof', evidence: ['https://example.com/proof'] }).state;
  state = apply(state, { kind: 'submit', task: 't', summary: 'Ready' }).state;
  state = apply(state, { kind: 'usage', task: 't', ...receipt() }).state;
  assert.equal(state.tasks[0].status, 'review'); assert.equal(state.tasks[0].leaseUntil, null);
  assert.throws(() => apply(state, { kind: 'usage', task: 't', agent: 'b', ...receipt({runId:'different'}) }, 'bob'), { code: 'USAGE_CONTRIBUTION_REQUIRED' });
  assert.throws(() => apply(state, { kind: 'usage', task: 't', ...receipt() }, 'mallory'), { code: 'FORBIDDEN' });
  const expired = fixture(); expired.tasks[0].leaseUntil = '2026-11-25T01:00:00Z';
  assert.equal(apply(expired, { kind: 'usage', task: 't', ...receipt() }).state.usage.length, 1);
});

test('usage accepts an independent reviewer of the task under the same owner only', () => {
  let state = fixture();
  state = apply(state, { kind: 'checkpoint', task: 't', summary: 'Proof', evidence: ['https://example.com/proof'] }).state;
  state = apply(state, { kind: 'submit', task: 't', summary: 'Ready' }).state;
  state = apply(state, { kind: 'join', agent: 'c', name: 'C' }, 'carol').state;
  assert.throws(() => apply(state, { kind: 'usage', task: 't', agent: 'c', ...receipt({runId:'carol-run'}) }, 'carol'), { code: 'USAGE_CONTRIBUTION_REQUIRED' });
  state = apply(state, { kind: 'review', task: 't', agent: 'b', checkpointId: state.tasks[0].submittedCheckpointId, verdict: 'pass', summary: 'Checked', evidence: ['https://example.com/review'] }, 'bob').state;
  assert.equal(apply(state, { kind: 'usage', task: 't', agent: 'b', ...receipt({runId:'bob-run'}) }, 'bob').usage.agentId, 'b');
  const archived = structuredClone(state); archived.tasks[0].reviews = [];
  assert.equal(apply(archived, { kind: 'usage', task: 't', agent: 'b', ...receipt({runId:'bob-run'}) }, 'bob').usage.agentId, 'b');
  for (const event of archived.events) if (event.kind === 'review') event.actorOwner = 'mallory';
  assert.throws(() => apply(archived, { kind: 'usage', task: 't', agent: 'b', ...receipt({runId:'bob-run'}) }, 'bob'), { code: 'USAGE_CONTRIBUTION_REQUIRED' });
});

test('same run deduplicates across request IDs and same-owner agent aliases', () => {
  let state = fixture();
  const first = apply(state, { id: 'first-usage', kind: 'usage', task: 't', ...receipt() }); state = first.state;
  const duplicate = apply(state, { id: 'retry-new-id', kind: 'usage', task: 't', ...receipt() });
  assert.equal(duplicate.duplicate, true); assert.equal(duplicate.state.revision, state.revision);
  state = apply(state, { kind: 'release', task: 't', summary: 'Handoff' }).state;
  state = apply(state, { kind: 'claim', task: 't', agent: 'a-alt' }).state;
  const alias = apply(state, { id: 'alias-retry', kind: 'usage', agent: 'a-alt', task: 't', ...receipt({provider:'CODEX'}) });
  assert.equal(alias.duplicate, true); assert.equal(alias.state.usage.length, 1);
  assert.equal(alias.usage.agentId, 'a'); assert.equal(alias.usage.id, 'first-usage');
});

test('run identity conflicts reject changed totals, cache, task, timing, model, and source', () => {
  let state = apply(fixture(), { kind: 'usage', task: 't', ...receipt() }).state;
  for (const change of [{inputTokens:121}, {outputTokens:31}, {cachedInputTokens:79}, {finishedAt:'2026-11-25T06:00:00Z'}, {model:'another-model'}, {source:'manual'}]) {
    assert.throws(() => apply(state, { kind: 'usage', task: 't', ...receipt(change) }), { code: 'USAGE_RUN_CONFLICT' });
  }
  state = apply(state, { kind: 'release', task: 't', summary: 'Done' }).state;
  state = apply(state, { kind: 'task', taskId: 'other', objectiveId: 'p', title: 'Other', description: 'Other work' }).state;
  state = apply(state, { kind: 'claim', task: 'other' }).state;
  assert.throws(() => apply(state, { kind: 'usage', task: 'other', ...receipt() }), { code: 'USAGE_RUN_CONFLICT' });
});

test('the same run string under different owners or providers is independently namespaced', () => {
  let state = apply(fixture(), { kind: 'usage', task: 't', ...receipt() }).state;
  state = apply(state, { kind: 'usage', task: 't', ...receipt({provider:'other-runtime'}) }).state;
  state = apply(state, { kind: 'release', task: 't', summary: 'Handoff' }).state;
  state = apply(state, { kind: 'claim', task: 't', agent: 'b' }, 'bob').state;
  state = apply(state, { kind: 'usage', task: 't', agent: 'b', ...receipt() }, 'bob').state;
  assert.equal(state.usage.length, 3);
});

test('token validation rejects negatives, fractions, strings, overflow, and additive cache totals', () => {
  for (const field of ['inputTokens', 'outputTokens', 'cachedInputTokens']) {
    for (const value of [-1, 0.1, '2', true, null, NaN, Infinity, Number.MAX_SAFE_INTEGER + 1]) {
      assert.throws(() => normalizeUsage(receipt({[field]:value})), /safe integer/);
    }
  }
  assert.throws(() => normalizeUsage(receipt({cachedInputTokens:121})), /subset/);
  assert.throws(() => normalizeUsage(receipt({inputTokens:Number.MAX_SAFE_INTEGER,outputTokens:1})), /Total tokens/);
  assert.throws(() => normalizeUsage(receipt({totalTokens:230})), /not additive/);
  assert.equal(normalizeUsage(receipt({inputTokens:Number.MAX_SAFE_INTEGER,outputTokens:0})).totalTokens, Number.MAX_SAFE_INTEGER);
  for (const value of ['1e3','-1','1.5','', 'Infinity']) assert.throws(() => parseTokenFlag(value,'Tokens'));
});

test('receipt dates and metadata are validated and normalized', () => {
  assert.equal(normalizeUsage(receipt({startedAt:'2026-11-24T22:58:00-07:00'})).startedAt, '2026-11-25T05:58:00.000Z');
  for (const change of [{startedAt:'tomorrow'}, {startedAt:'2026-02-31T00:00:00Z'}, {startedAt:'2026-11-25T07:00:00Z'}, {finishedAt:'2026-11-25T24:00:00Z'}, {source:'estimated'}, {provider:''}, {runId:'has spaces'}, {model:'x'.repeat(161)}, {evidence:['javascript:bad']}]) assert.throws(() => normalizeUsage(receipt(change)));
});

test('state validation rejects corrupt aggregate receipts and duplicate runtime runs', () => {
  const state = apply(fixture(), { kind: 'usage', task: 't', ...receipt() }).state;
  const duplicate = structuredClone(state); duplicate.usage.push({...duplicate.usage[0],id:'other'});
  assert.throws(() => validateState(duplicate), /more than once/);
  const changed = structuredClone(state); changed.usage[0].totalTokens++;
  assert.throws(() => validateState(changed), /Total tokens/);
  assert.throws(() => validateState({...state,usage:{}}), /array/);
});

test('missing or invalid worker capture does not write zero or any fabricated receipt', async () => {
  const store = mockStore(), messages = [];
  const missing = await publishWorkerUsage(store, {agent:'a',owner:'alice',task:'t',onProgress:m=>messages.push(m)});
  const invalid = await publishWorkerUsage(store, {agent:'a',owner:'alice',task:'t',receipt:{},onProgress:m=>messages.push(m)});
  assert.equal(missing.status, 'missing'); assert.equal(invalid.status, 'invalid');
  assert.equal(store.state.usage, undefined); assert.equal(messages.length, 2);
});

test('failed upload preserves a private recovery envelope and exact replay deduplicates', async () => temporary(async directory => {
  const store = mockStore(); await store.mutate({kind:'claim',agent:'a',task:'t'});
  const mutate = store.mutate;
  store.mutate = async op => { if (op.kind === 'usage') throw new Error('offline'); return mutate(op); };
  const pending = await publishWorkerUsage(store, {agent:'a',owner:'alice',task:'t',receipt:receipt(),recoveryDirectory:directory});
  assert.equal(pending.status, 'pending'); assert.ok(pending.recoveryPath);
  const envelope = JSON.parse(await readFile(pending.recoveryPath,'utf8'));
  assert.equal(envelope.receipt.totalTokens,150);
  assert.equal((await stat(pending.recoveryPath)).mode & 0o777, 0o600);
  const op = await usageOperationFromFlags({'usage-file':pending.recoveryPath},store.config);
  assert.equal(op.id,pending.requestId); assert.equal(op.task,'t');
  store.mutate = mutate;
  await store.mutate(op);
  const replay = await publishWorkerUsage(store, {agent:'a',owner:'alice',task:'t',receipt:receipt(),recoveryDirectory:directory});
  assert.equal(replay.status,'already-recorded'); assert.equal(store.state.usage.length,1);
  assert.deepEqual(await readdir(directory),[]);
}));

test('response lost after recording can be replayed from the retained receipt without double counting', async () => temporary(async directory => {
  const store=mockStore(); await store.mutate({kind:'claim',agent:'a',task:'t'});
  const mutate=store.mutate; let lost=true;
  store.mutate=async op=>{ const result=await mutate(op); if(op.kind==='usage'&&lost){lost=false;const error=new Error('Lost response');error.code='WRITE_UNCERTAIN';throw error;} return result; };
  const pending=await publishWorkerUsage(store,{agent:'a',owner:'alice',task:'t',receipt:receipt(),recoveryDirectory:directory});
  assert.equal(pending.status,'pending');assert.equal(store.state.usage.length,1);
  const retry=await publishWorkerUsage(store,{agent:'a',owner:'alice',task:'t',receipt:receipt(),recoveryDirectory:directory});
  assert.equal(retry.status,'already-recorded');assert.equal(store.state.usage.length,1);
}));

test('CLI flags and recovery-file identity remain explicit', async () => temporary(async directory => {
  const flags={agent:'a',task:'t',provider:'codex',model:'unspecified','run-id':'manual-run','input-tokens':'120','output-tokens':'30','cached-input-tokens':'80','started-at':'2026-11-25T05:58:00Z','finished-at':'2026-11-25T05:59:00Z'};
  const op=await usageOperationFromFlags(flags,{repository:'owner/repo'});
  assert.equal(op.source,'manual');assert.equal(op.totalTokens,150);
  const missingCache={...flags};delete missingCache['cached-input-tokens'];
  await assert.rejects(usageOperationFromFlags(missingCache,{repository:'owner/repo'}),/Cached input tokens/);
  const path=join(directory,'recovery.json');await writeFile(path,JSON.stringify({schemaVersion:1,repository:'owner/repo',agent:'a',task:'t',requestId:'stable-id',receipt:receipt()}));
  await assert.rejects(usageOperationFromFlags({'usage-file':path,agent:'other'},{repository:'owner/repo'}),/differs/);
  await assert.rejects(usageOperationFromFlags({'usage-file':path},{repository:'another/repo'}),/another repository/);
  await assert.rejects(usageOperationFromFlags({'usage-file':path,provider:'codex'},{repository:'owner/repo'}),/mixed/);
  assert.equal(usageOperationId('Alice','a',normalizeUsage(receipt())),usageOperationId('alice','a',normalizeUsage(receipt())));
}));

test('worker records actual optional usage and surfaces runtime warnings without replacing research outcome', async () => temporary(async directory => {
  const store=mockStore(),progress=[];
  const outcomes=await work(store,{agent:'a',command:'fake',recoveryDirectory:directory,onProgress:event=>progress.push(event),run:async()=>({status:'release',summary:'Bounded run complete',usage:receipt(),usageWarning:'Contributor-reported runtime counts.'})});
  assert.equal(outcomes[0].status,'release');assert.equal(outcomes[0].usage.status,'recorded');
  assert.equal(store.state.usage[0].totalTokens,150);
  assert.ok(progress.some(event=>event.kind==='usage-warning'&&event.summary==='Contributor-reported runtime counts.'));
  const missingStore=mockStore();
  const missing=await work(missingStore,{agent:'a',command:'fake',recoveryDirectory:directory,run:async()=>({status:'release',summary:'No counter available'})});
  assert.equal(missing[0].usage.status,'missing');assert.equal(missingStore.state.usage,undefined);
}));

test('usage upload failure preserves research handoff and returns an explicit pending receipt', async () => temporary(async directory => {
  const store=mockStore(),mutate=store.mutate;
  store.mutate=async op=>{if(op.kind==='usage')throw new Error('offline');return mutate(op);};
  const outcomes=await work(store,{agent:'a',command:'fake',recoveryDirectory:directory,run:async()=>({status:'release',summary:'Actual work finished',usage:receipt()})});
  assert.equal(outcomes[0].status,'release');assert.equal(outcomes[0].usage.status,'pending');
  assert.equal(store.state.tasks[0].status,'open');assert.equal(store.state.usage,undefined);
  assert.equal((await readdir(directory)).length,1);
}));


test('future usage beyond five-minute server-clock tolerance is rejected', () => {
  const state=fixture();
  assert.throws(()=>apply(state,{kind:'usage',task:'t',...receipt({finishedAt:'2026-11-25T06:05:01Z'})}),{code:'USAGE_FUTURE_TIMESTAMP'});
  assert.equal(apply(state,{kind:'usage',task:'t',...receipt({finishedAt:'2026-11-25T06:05:00Z'})}).state.usage.length,1);
});

test('dual local-save and upload failure still returns the captured recovery receipt', async () => temporary(async directory => {
  const store=mockStore();await store.mutate({kind:'claim',agent:'a',task:'t'});
  const file=join(directory,'not-a-directory');await writeFile(file,'block mkdir');
  store.mutate=async()=>{throw new Error('offline');};
  const pending=await publishWorkerUsage(store,{agent:'a',owner:'alice',task:'t',receipt:receipt(),recoveryDirectory:file});
  assert.equal(pending.status,'pending');assert.equal(pending.recoveryPath,null);
  assert.equal(pending.recovery.receipt.totalTokens,150);
  assert.equal(pending.recovery.receipt.runId,'codex:run-one');
}));


test('runtime and manual receipts both require explicit cache counts', () => {
  const missing=receipt();delete missing.cachedInputTokens;
  assert.throws(()=>normalizeUsage(missing),/missing capture is not zero/);
  assert.throws(()=>normalizeUsage({...missing,source:'manual'}),/missing capture is not zero/);
});

test('subprocess result errors retain captured usage instead of losing a real completed run', async () => temporary(async directory => {
  const packet={agent:{id:'a'},repository:'owner/repo'};
  for(const code of [0,2]){
    const value={status:'review',summary:'Missing evidence',usage:receipt(),usageWarning:'Captured counts survived failure.'};
    const script=`process.stdout.write(${JSON.stringify(JSON.stringify(value))});process.exitCode=${code};`;
    await assert.rejects(runAdapter(process.execPath,['-e',script],packet),error=>error.usage?.runId==='codex:run-one'&&error.usageWarning==='Captured counts survived failure.');
  }
  const store=mockStore();
  const result=await work(store,{agent:'a',command:process.execPath,args:['-e',`process.stdout.write(${JSON.stringify(JSON.stringify({status:'review',summary:'Incomplete result',usage:receipt()}))});`],recoveryDirectory:directory});
  assert.equal(result[0].failed,true);assert.equal(result[0].usage.status,'recorded');
  assert.equal(store.state.usage[0].totalTokens,150);assert.equal(store.state.tasks[0].status,'blocked');
}));
