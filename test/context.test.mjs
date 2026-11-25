import test from 'node:test';
import assert from 'node:assert/strict';
import { liveContext, modelContext, taskSignals } from '../lib/context.mjs';
import { mkdtemp, writeFile, rm } from 'node:fs/promises';
import { join } from 'node:path';
import { tmpdir } from 'node:os';
import { runAdapter } from '../lib/worker.mjs';
import { main } from '../bin/commons.mjs';
import { readFileSync } from 'node:fs';
import { requestSourceContext } from '../lib/worker-context.mjs';

const now = Date.parse('2026-11-25T18:00Z');
const objective = { id: 'o', status: 'open', statement: 'x'.repeat(600), successCriteria: 'Exact criteria', example: false };
const agent = { id: 'a', owner: 'alice' };
const task = { id: 't', objectiveId: 'o', status: 'active', scopePaths: ['owned'], assignee: 'a', description: 'Prove exactly this.',
  leaseUntil: new Date(now - 1000).toISOString(), checkpoints: [], reviews: [], dependsOn: ['source'], request: { fromTaskId: 'source', to: 'a' } };
const checkpoint = { id: 'cp', summary: 'Bounded case', details: 'full evidence '.repeat(600), evidence: Array.from({length: 7}, (_, i) => 'https://example.net/' + i) };
const source = { id: 'source', objectiveId: 'o', status: 'verified', reviews: [], dependsOn: [], checkpoints: [checkpoint] };
const event = (id, kind, taskId, extra = {}) => ({ id, kind, taskId, objectiveId: 'o', summary: id, detail: 'full evidence', ...extra });
const state = { revision: 10, objectives: [objective], agents: [agent], tasks: [task, source], events: [] };

test('the objective card pages tasks, flags abridged contracts and retrieves the exact objective alone', () => {
  const tasks = Array.from({length: 19}, (_, i) => ({...source, id: 't' + i}));
  const s = {...state, tasks};
  const first = liveContext(s, {objective:'o'}, now), second = liveContext(s, {objective:'o',page:'2'}, now);
  assert.equal(first.tasks.total, 19); assert.equal(first.tasks.pages, 3); assert.equal(first.statusCounts.verified, 19);
  assert.equal(first.tasks.items.length, 8); assert.equal(new Set([...first.tasks.items,...second.tasks.items].map(x=>x.id)).size, 16);
  assert.match(first.objective.exactContract, /context --objective o --exact/);
  assert.equal(JSON.stringify(first).includes('full evidence'), false);
  assert.deepEqual(liveContext(s, {objective:'o',exact:true}, now), {revision:10, objective});
  for(const opts of [{objective:'o',page:'4'},{objective:'o',page:'-1'},{objective:'o',agent:'unknown'},{task:'t',objective:'o'},{task:'t',exact:true},{event:'cp',page:'2'},{event:'absent'}])
    assert.throws(() => liveContext(s, opts, now));
});

test('a closed synthesis stays attributable after reopening and the exact leaf retains its evidence', () => {
  const closed = {...objective, statement: 'Short', status: 'closed', conclusion: 'Prior conclusion '.repeat(50), closedAt: '2026-11-20T00:00Z', evidence: ['https://example.net/report']};
  const closedState = {...state, objectives: [closed], events: [event('end', 'close', null)]};
  const first = liveContext(closedState, {objective:'o'}, now).objective;
  assert.equal(first.synthesis.evidenceCount, 1); assert.equal(first.latestDecision.id,'end');
  assert.match(first.exactContract, /--exact/); assert.ok(first.synthesis.summary.length < closed.conclusion.length);
  const reopened = {...closedState, objectives:[{...closed, status:'open'}], events:[...closedState.events,event('again','reopen-objective',null)]};
  const view = liveContext(reopened,{task:'t'},now).objective;
  assert.equal(view.synthesis,undefined); assert.equal(view.previousSynthesis.evidenceCount,1); assert.equal(view.latestDecision.id,'again');
  assert.equal(liveContext(reopened,{objective:'o',exact:true},now).objective.conclusion,closed.conclusion);
});

test('active reservations and newly claimable requests surface ahead of older completed branches', () => {
  const closed = Array.from({ length: 12 }, (_, i) => ({ ...source, id: 'old-' + i, updatedAt: '2026-11-26T01:00Z' }));
  const pending = { ...source, id: 'requested', status: 'open', updatedAt: '2026-11-25T09:00Z' };
  const active = { ...task, dependsOn: [], request: undefined, leaseUntil: new Date(now + 1000).toISOString() };
  const view = liveContext({ ...state, tasks: [...closed, pending, active] }, { objective: 'o' }, now);
  assert.deepEqual(view.tasks.items.slice(0, 2).map(item => item.id), ['t', 'requested']);
  assert.equal(view.claimableTotal, 1); assert.equal(view.statusCounts.verified, 12);
});

test('task scope, lease expiry and submitted checkpoint stay exact; full legacy records can be addressed by id', () => {
  const review = {id:'r', verdict:'changes', checkpointId:'cp', summary:'Check the boundary.'};
  const selected = {...task, checkpoints:[checkpoint], reviews:[review], submittedCheckpointId:'cp'};
  const s = {...state,tasks:[selected,{...source,checkpoints:[]}]};
  const card = liveContext(s,{task:'t'},now);
  assert.equal(card.task.leaseExpired,true); assert.equal(card.task.submittedCheckpointId,'cp');
  assert.deepEqual(card.task.scopePaths,['owned']); assert.equal(card.task.description,task.description);
  assert.equal(card.task.latestCheckpoint.evidenceCount,7); assert.equal(card.task.latestCheckpoint.evidence.length,1);
  assert.equal(card.task.latestReview.checkpointId,'cp'); assert.doesNotMatch(JSON.stringify(card),/full evidence/);
  assert.deepEqual(liveContext(s,{event:'cp'},now).checkpoint,checkpoint);
  const r = event('r','review','t'); assert.deepEqual(liveContext({...s,events:[r]},{event:'r'},now), {revision:10,event:r,verdict:'changes',checkpointId:'cp'});
  assert.deepEqual(liveContext(s,{event:'r'},now).review,review);
  assert.throws(() => liveContext({...s,tasks:[selected,source]},{event:'cp'},now), /Ambiguous legacy id/);
});

test('a busy feed retains a direct message and upstream decision; other objectives and accounting history stay outside the card', () => {
  const events = [event('old-msg','message',null,{to:'a',objectiveId:null}),event('source-block','block','source'),
    event('foreign-msg','message','foreign-task',{to:'a',objectiveId:'foreign'}),event('noise','usage','t'),
    ...Array.from({length:14},(_,i)=>event('cp'+i,'checkpoint','t'))];
  const signals = taskSignals(events,task,'a');
  assert.equal(signals.total,16); assert.deepEqual(signals.earlierSignals.map(e=>e.id), ['old-msg','source-block']);
  assert.ok(!JSON.stringify(signals).includes('full evidence'));
  const older = taskSignals(events,task,'a',2);
  assert.equal(older.items.length,8); assert.equal(older.earlierSignals,undefined);
  assert.equal(older.items.some(e=>e.id==='foreign-msg'||e.id==='noise'),false);
  const card = modelContext({repository:'alice/repo',agent,objective,task,dependencies:[source],sourceTask:requestSourceContext(state,task),recentEvents:events});
  assert.equal(card.task.leaseExpired,undefined);
  const claimed = modelContext({agent,task,snapshotAt: new Date(now - 2000).toISOString()});
  assert.equal(claimed.task.leaseExpired,false);
  assert.equal(card.sourceTask.latestCheckpoint.id,'cp'); assert.equal(card.sourceTask.latestCheckpoint.evidenceCount,7);
  assert.deepEqual(card.signals,signals); assert.doesNotMatch(JSON.stringify(card),/full evidence/);
});

test('CLI context can read a legacy objective without writing or loading other objectives into the result', async () => {
  const seed = readFileSync(new URL('../coordination/state.json',import.meta.url));
  const fetch = globalThis.fetch, log = console.log, token = process.env.GH_TOKEN;
  let output;
  try {
    process.env.GH_TOKEN='fixture-token';
    globalThis.fetch=async (_url, options) => {
      assert.equal(options.method,'GET');
      return new Response(JSON.stringify({type:'file',encoding:'base64',sha:'test',content:seed.toString('base64')}));
    };
    console.log = value => { output = value; };
    await main(['context','--objective','example-sum-of-cubes','--json']);
    assert.equal(JSON.parse(output).objective.example,true);
    await main(['context','--task','example-proof','--json']);
    assert.equal(JSON.parse(output).task.id,'example-proof');
  } finally { globalThis.fetch=fetch; console.log=log; if(token===undefined) delete process.env.GH_TOKEN; else process.env.GH_TOKEN=token; }
});

test('an alternate worker configuration follows a child lookup to the exact branch and state file', async () => {
  const dir = await mkdtemp(join(tmpdir(), 'commons alternate config '));
  try {
    const configPath = join(dir, 'custom config.json');
    await writeFile(configPath, JSON.stringify({repository:'owner/repo',stateBranch:'alternate',statePath:'alternate/state.json'}));
    const moduleURL = new URL('../bin/commons.mjs', import.meta.url).href;
    const seedURL = new URL('../coordination/state.json', import.meta.url).href;
    const code = `
      const { main } = await import(${JSON.stringify(moduleURL)});
      const { readFileSync } = await import('node:fs');
      const seed = readFileSync(new URL(${JSON.stringify(seedURL)}));
      globalThis.fetch = async (url, options) => {
        if (url !== 'https://api.github.com/repos/owner/repo/contents/alternate/state.json?ref=alternate' || options.method !== 'GET')
          throw new Error('Lookup used the wrong ledger');
        return new Response(JSON.stringify({type:'file',encoding:'base64',sha:'test',content:seed.toString('base64')}));
      };
      let out; const log = console.log; console.log = value => { out = JSON.parse(value) };
      await main(['context','--task','example-proof']);
      log(JSON.stringify({status:'release',summary:out.configPath,details:out.task.id,evidence:[]}));
    `;
    const result = await runAdapter(process.execPath,['--input-type=module','-e',code],{agent,repository:'owner/repo'},{configPath});
    assert.equal(result.summary,configPath); assert.equal(result.details,'example-proof');
  } finally { await rm(dir,{recursive:true,force:true}); }
});
