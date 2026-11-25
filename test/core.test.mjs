import test from 'node:test';
import assert from 'node:assert/strict';
import { applyOperation, eligibleTasks, requestDirectory, validateState, MAX_STATE_BYTES } from '../lib/core.mjs';
import { parseArgs, main } from '../bin/commons.mjs';

const START = '2026-11-25T00:00:00.000Z';
let sequence = 0;
const empty = () => ({ schemaVersion: 1, revision: 0, updatedAt: START, objectives: [], tasks: [], agents: [], events: [] });
function apply(state, op, owner = 'alice', now = START) {
  return applyOperation(state, { id: `request-${++sequence}`, agent: 'alice-agent', ...op }, { owner, now });
}
function fixture() {
  let state = apply(empty(), { kind: 'join', name: 'Alice' }).state;
  state = apply(state, { kind: 'join', agent: 'bob-agent', name: 'Bob' }, 'bob').state;
  state = apply(state, { kind: 'objective', objectiveId: 'problem', title: 'A problem', statement: 'Prove a theorem', success: 'Independent review' }).state;
  state = apply(state, { kind: 'task', taskId: 'lemma', objectiveId: 'problem', title: 'Prove the lemma', description: 'All cases' }).state;
  return state;
}
function submitted() {
  let state = fixture();
  state = apply(state, { kind: 'claim', task: 'lemma' }).state;
  state = apply(state, { kind: 'checkpoint', task: 'lemma', summary: 'Proof', details: 'Evidence', evidence: ['https://github.com/org/repo/blob/abc/proof.md'] }).state;
  return apply(state, { kind: 'submit', task: 'lemma', summary: 'Ready' }).state;
}

test('operations preserve input, identity, audit trail, and exact idempotency', () => {
  const input = empty(), op = { id: 'stable', kind: 'join', agent: 'a', name: 'A' };
  const result = applyOperation(input, op, { owner: 'alice', now: START });
  assert.equal(input.agents.length, 0);
  assert.equal(result.state.events[0].actorOwner, 'alice');
  const retry = applyOperation(result.state, op, { owner: 'ALICE', now: START });
  assert.equal(retry.duplicate, true); assert.equal(retry.state.revision, 1);
  assert.throws(() => applyOperation(result.state, { ...op, name: 'Changed' }, { owner: 'alice', now: START }), { code: 'IDEMPOTENCY_CONFLICT' });
  assert.throws(() => applyOperation(result.state, op, { owner: 'mallory', now: START }), { code: 'IDEMPOTENCY_CONFLICT' });
});

test('agent IDs cannot be taken over by other accounts', () => {
  const state = fixture();
  assert.throws(() => apply(state, { kind: 'join', name: 'New owner' }, 'mallory'), { code: 'FORBIDDEN' });
  assert.throws(() => apply(state, { kind: 'claim', task: 'lemma' }, 'mallory'), { code: 'FORBIDDEN' });
});

test('claim is exclusive and expiration allows another registered agent', () => {
  let state = apply(fixture(), { kind: 'claim', task: 'lemma' }).state;
  assert.throws(() => apply(state, { kind: 'claim', task: 'lemma', agent: 'bob-agent' }, 'bob'), { code: 'TASK_BUSY' });
  assert.throws(() => apply(state, { kind: 'checkpoint', task: 'lemma', agent: 'bob-agent', summary: 'No' }, 'bob'), { code: 'LEASE_REQUIRED' });
  assert.throws(() => apply(state, { kind: 'checkpoint', task: 'lemma', summary: 'Expired' }, 'alice', '2026-11-25T00:21:00Z'), { code: 'LEASE_REQUIRED' });
  state = apply(state, { kind: 'claim', task: 'lemma', agent: 'bob-agent' }, 'bob', '2026-11-25T00:21:00Z').state;
  assert.equal(state.tasks[0].assignee, 'bob-agent');
  assert.match(state.events.at(-1).summary, /Reclaimed/);
});

test('heartbeat renews only live own leases and agent cannot hoard two active tasks', () => {
  let state = fixture();
  state = apply(state, { kind: 'task', taskId: 'other', objectiveId: 'problem', title: 'Other task', description: 'Another direction' }).state;
  state = apply(state, { kind: 'claim', task: 'lemma' }).state;
  assert.throws(() => apply(state, { kind: 'claim', task: 'other' }), { code: 'AGENT_BUSY' });
  state = apply(state, { kind: 'heartbeat', task: 'lemma' }, 'alice', '2026-11-25T00:10:00Z').state;
  assert.equal(state.tasks[0].leaseUntil, '2026-11-25T00:30:00.000Z');
  assert.throws(() => apply(state, { kind: 'heartbeat', task: 'lemma' }, 'alice', '2026-11-25T00:31:00Z'), { code: 'LEASE_REQUIRED' });
});

test('dependencies require verified results and reject invalid references', () => {
  let state = fixture();
  state = apply(state, { kind: 'task', taskId: 'synthesis', objectiveId: 'problem', title: 'Synthesis', description: 'Use lemma', dependsOn: ['lemma'] }).state;
  assert.deepEqual(eligibleTasks(state).map(task => task.id), ['lemma']);
  assert.throws(() => apply(state, { kind: 'claim', task: 'synthesis' }), { code: 'NOT_ELIGIBLE' });
  assert.throws(() => apply(state, { kind: 'task', taskId: 'bad', objectiveId: 'problem', title: 'Bad', description: 'Unknown', dependsOn: ['missing'] }), /Unknown task/);
  const corrupt = structuredClone(state); corrupt.tasks[0].dependsOn = ['synthesis'];
  assert.throws(() => validateState(corrupt), /cycle/);
});

test('submission requires own latest checkpoint with reproducible evidence', () => {
  let state = apply(fixture(), { kind: 'claim', task: 'lemma' }).state;
  assert.throws(() => apply(state, { kind: 'submit', task: 'lemma', summary: 'Ready' }), /checkpoint/);
  state = apply(state, { kind: 'checkpoint', task: 'lemma', summary: 'Progress', details: 'Still exploring' }).state;
  assert.throws(() => apply(state, { kind: 'submit', task: 'lemma', summary: 'Ready' }), /evidence/);
  for (const url of ['javascript:alert(1)', 'http://example.com', 'https://user:pass@example.com']) {
    assert.throws(() => apply(state, { kind: 'checkpoint', task: 'lemma', summary: 'Unsafe URL', evidence: [url] }), /https URL/);
  }
});

test('independent review rejects same-owner sockpuppets and stale checkpoint IDs', () => {
  let state = submitted();
  state = apply(state, { kind: 'join', agent: 'alice-reviewer', name: 'Alice reviewer' }).state;
  const review = { kind: 'review', task: 'lemma', checkpointId: state.tasks[0].submittedCheckpointId, verdict: 'pass', summary: 'Checked all cases', evidence: ['https://example.com/review'] };
  assert.throws(() => apply(state, { ...review, agent: 'alice-reviewer' }), { code: 'INDEPENDENT_REVIEW_REQUIRED' });
  assert.throws(() => apply(state, { ...review, checkpointId: 'outdated', agent: 'bob-agent' }, 'bob'), { code: 'STALE_REVIEW' });
  state = apply(state, { ...review, agent: 'bob-agent' }, 'bob').state;
  assert.equal(state.tasks[0].status, 'verified');
  assert.equal(state.tasks[0].reviews[0].checkpointId, state.tasks[0].checkpoints[0].id);
});

test('reclaimed task reviewers cannot be any earlier checkpoint contributor', () => {
  let state = apply(fixture(), { kind: 'claim', task: 'lemma' }).state;
  state = apply(state, { kind: 'checkpoint', task: 'lemma', summary: 'First contribution', evidence: ['https://example.com/a'] }).state;
  state = apply(state, { kind: 'release', task: 'lemma', summary: 'Continue from this result' }).state;
  state = apply(state, { kind: 'claim', task: 'lemma', agent: 'bob-agent' }, 'bob').state;
  state = apply(state, { kind: 'checkpoint', task: 'lemma', agent: 'bob-agent', summary: 'Second contribution', evidence: ['https://example.com/b'] }, 'bob').state;
  state = apply(state, { kind: 'submit', task: 'lemma', agent: 'bob-agent', summary: 'Ready' }, 'bob').state;
  assert.throws(() => apply(state, { kind: 'review', task: 'lemma', checkpointId: state.tasks[0].submittedCheckpointId, verdict: 'pass', summary: 'Looks right', evidence: ['https://example.com/review'] }), { code: 'INDEPENDENT_REVIEW_REQUIRED' });
});

test('changes requested reopen work; new submission cannot use an old review', () => {
  let state = submitted();
  const oldCheckpoint = state.tasks[0].submittedCheckpointId;
  state = apply(state, { kind: 'review', task: 'lemma', agent: 'bob-agent', checkpointId: oldCheckpoint, verdict: 'changes', summary: 'Missing n=0', evidence: ['https://example.com/review'] }, 'bob').state;
  assert.equal(state.tasks[0].status, 'open');
  state = apply(state, { kind: 'claim', task: 'lemma' }).state;
  state = apply(state, { kind: 'checkpoint', task: 'lemma', summary: 'Fixed', evidence: ['https://example.com/fix'] }).state;
  state = apply(state, { kind: 'submit', task: 'lemma', summary: 'Ready again' }).state;
  assert.throws(() => apply(state, { kind: 'review', task: 'lemma', agent: 'bob-agent', checkpointId: oldCheckpoint, verdict: 'pass', summary: 'Good', evidence: ['https://example.com/review'] }, 'bob'), { code: 'STALE_REVIEW' });
});

test('handoffs and blockers preserve work and reopening requires a reason', () => {
  let state = apply(fixture(), { kind: 'claim', task: 'lemma' }).state;
  state = apply(state, { kind: 'checkpoint', task: 'lemma', summary: 'Found an obstruction' }).state;
  state = apply(state, { kind: 'block', task: 'lemma', summary: 'Need an assumption clarified' }).state;
  assert.equal(state.tasks[0].assignee, null); assert.equal(state.tasks[0].checkpoints.length, 1);
  assert.equal(eligibleTasks(state).length, 0);
  assert.throws(() => apply(state, { kind: 'reopen', task: 'lemma' }), /Reason/);
  state = apply(state, { kind: 'reopen', task: 'lemma', summary: 'Assumption clarified' }).state;
  assert.equal(eligibleTasks(state).length, 1);
});

test('reopening a verified dependency cannot silently invalidate active downstream work', () => {
  let state = submitted();
  state = apply(state, { kind: 'review', task: 'lemma', agent: 'bob-agent', checkpointId: state.tasks[0].submittedCheckpointId, verdict: 'pass', summary: 'Checked', evidence: ['https://example.com/review'] }, 'bob').state;
  state = apply(state, { kind: 'task', taskId: 'synthesis', objectiveId: 'problem', title: 'Combine', description: 'Use lemma', dependsOn: ['lemma'] }).state;
  state = apply(state, { kind: 'claim', task: 'synthesis' }).state;
  assert.throws(() => apply(state, { kind: 'reopen', task: 'lemma', summary: 'Found a gap' }), { code: 'DEPENDENTS_ACTIVE' });
});

test('only creator owner can close independently reviewed objective with synthesis', () => {
  let state = submitted();
  assert.throws(() => apply(state, { kind: 'close', objectiveId: 'problem', summary: 'Done', evidence: ['https://example.com/final'] }), /Every task/);
  state = apply(state, { kind: 'review', task: 'lemma', agent: 'bob-agent', checkpointId: state.tasks[0].submittedCheckpointId, verdict: 'pass', summary: 'Checked', evidence: ['https://example.com/review'] }, 'bob').state;
  assert.throws(() => apply(state, { kind: 'close', objectiveId: 'problem', agent: 'bob-agent', summary: 'Done', evidence: ['https://example.com/final'] }, 'bob'), { code: 'FORBIDDEN' });
  state = apply(state, { kind: 'close', objectiveId: 'problem', summary: 'Proved under stated assumptions', evidence: ['https://example.com/final'] }).state;
  assert.equal(state.objectives[0].status, 'closed');
  assert.equal(state.objectives[0].conclusion, 'Proved under stated assumptions');
});

test('state capacity checks the exact pretty-printed bytes written to GitHub', () => {
  const state = fixture();
  state.events[0].detail = 'x'.repeat(MAX_STATE_BYTES - Buffer.byteLength(JSON.stringify(state)) - 100);
  assert.ok(Buffer.byteLength(JSON.stringify(state)) < MAX_STATE_BYTES);
  assert.throws(() => validateState(state), { code: 'STATE_CAPACITY' });
});

test('CLI argument parser preserves shell-sensitive data and separates adapter args', () => {
  const parsed = parseArgs(['checkpoint', '--summary', '$(do-not-run) `literal`', '--agent', 'a']);
  assert.equal(parsed.flags.summary, '$(do-not-run) `literal`');
  assert.deepEqual(parseArgs(['work', '--agent', 'a', '--', 'node', 'adapter.mjs', '--model', 'test']).trailing, ['node', 'adapter.mjs', '--model', 'test']);
  assert.throws(() => parseArgs(['join', '--agent', 'a', '--agent', 'b']), /Repeated/);
});

test('scoped leases protect ancestors across objectives, allow siblings and release on expiry', () => {
  let state = fixture();
  state = apply(state, { kind: 'objective', objectiveId: 'second', title: 'Other problem', statement: 'Another theorem', success: 'Evidence' }).state;
  const propose = (id, paths) => { state = apply(state, { kind: 'task', objectiveId: 'second', taskId: id,
    title: id, description: 'Exact artifact', scopePaths: paths }).state; };
  propose('overlap', ['Research/Part/proof.md']);
  propose('sibling', ['research/partition']);
  state = apply(state, { kind: 'claim', task: 'lemma', scopePaths: ['research/part/'] }).state;
  assert.deepEqual(state.tasks[0].scopePaths, ['research/part']);
  assert.deepEqual(eligibleTasks(state, { now: Date.parse(START) }).map(t => t.id), ['sibling']);
  const claim = { kind: 'claim', agent: 'bob-agent', task: 'overlap' };
  assert.throws(() => apply(state, claim, 'bob'), { code: 'TASK_BUSY' });
  assert.throws(() => apply(state, { ...claim, scopePaths: ['elsewhere'] }, 'bob'), /declared scope/);
  const after = apply(state, claim, 'bob', '2026-11-25T00:21:00Z').state;
  assert.throws(() => apply(after, { kind: 'claim', task: 'lemma' }, 'alice', '2026-11-25T00:22:00Z'), { code: 'TASK_BUSY' });
  state = apply(state, { kind: 'release', task: 'lemma', summary: 'Handoff' }).state;
  assert.equal(apply(state, claim, 'bob').state.tasks[1].assignee, 'bob-agent');
});

test('scope declarations reject ambiguous and nonliteral paths', () => {
  for (const paths of [[], [''], ['../lib'], ['.git/config'], ['/src'], ['lib//src'], ['lib/*.js'], ['a\\b'], ['a/./b']]) {
    assert.throws(() => apply(fixture(), { kind: 'claim', task: 'lemma', scopePaths: paths }), /path/i);
  }
});

const researchRequest = changes => ({ kind: 'request', taskId: 'boundary-request', objectiveId: 'problem',
  title: 'Check the boundary independently', description: 'Write a bounded boundary check with exact evidence.',
  scopePaths: ['research/boundary/'], fromTaskId: 'lemma', to: 'bob-agent', ...changes });

test('request creates one ordinary scoped task and one attributable invitation event', () => {
  const input = fixture(), operation = researchRequest({ id: 'stable-research-request' });
  const result = apply(input, operation), task = result.state.tasks.at(-1);
  assert.equal(result.state.tasks.length, input.tasks.length + 1);
  assert.equal(result.state.events.length, input.events.length + 1);
  assert.equal(task.status, 'open'); assert.equal(task.assignee, null); assert.equal(task.leaseUntil, null);
  assert.equal(task.createdBy, 'alice-agent'); assert.equal(task.createdAt, START);
  assert.deepEqual(task.scopePaths, ['research/boundary']);
  assert.deepEqual(task.request, { fromTaskId: 'lemma', to: 'bob-agent' });
  assert.deepEqual(task.dependsOn, []); assert.deepEqual(task.checkpoints, []);
  assert.equal(result.event.kind, 'request'); assert.equal(result.event.to, 'bob-agent');
  assert.equal(result.event.actorOwner, 'alice'); assert.equal(result.event.taskId, task.id);
  assert.deepEqual(result.event.scopePaths, task.scopePaths); assert.equal(result.event.fromTaskId, 'lemma');
  assert.equal(result.event.detail, task.description);
  assert.equal(input.tasks.length, 1); assert.deepEqual(result.state.tasks[0], input.tasks[0]);
  const retry = apply(result.state, operation);
  assert.equal(retry.duplicate, true); assert.equal(retry.state.revision, result.state.revision);
  for (const change of [{ to: 'all' }, { fromTaskId: null }, { scopePaths: ['another/path'] }]) {
    assert.throws(() => apply(result.state, { ...operation, ...change }), { code: 'IDEMPOTENCY_CONFLICT' });
  }
});

test('request invitations never restrict claims and reuse the existing review lifecycle', () => {
  let state = apply(fixture(), researchRequest({ to: 'alice-agent' })).state;
  // Bob was not the invitee and can still claim an eligible task.
  state = apply(state, { kind: 'claim', agent: 'bob-agent', task: 'boundary-request' }, 'bob').state;
  state = apply(state, { kind: 'checkpoint', agent: 'bob-agent', task: 'boundary-request', summary: 'Checked exact boundary', evidence: ['https://example.com/proof'] }, 'bob').state;
  state = apply(state, { kind: 'submit', agent: 'bob-agent', task: 'boundary-request', summary: 'Ready for review' }, 'bob').state;
  const task = state.tasks.at(-1);
  state = apply(state, { kind: 'review', task: task.id, checkpointId: task.submittedCheckpointId, verdict: 'pass', summary: 'Replayed the exact cases', evidence: ['https://example.com/review'] }).state;
  assert.equal(state.tasks.at(-1).status, 'verified');
  assert.deepEqual(state.tasks.at(-1).request, { fromTaskId: 'lemma', to: 'alice-agent' });
});

test('requests require valid paths, a known invitation target and a same-objective origin', () => {
  const input = fixture();
  for (const change of [{ scopePaths: undefined }, { scopePaths: [] }, { scopePaths: ['../outside'] }, { to: 'missing-agent' },
    { fromTaskId: 'missing-task' }, { fromTaskId: 'boundary-request' }, { taskId: 'lemma' }]) {
    assert.throws(() => apply(input, researchRequest(change)));
  }
  let state = apply(input, { kind: 'objective', objectiveId: 'other', title: 'Other', statement: 'Other scope', success: 'Evidence' }).state;
  state = apply(state, { kind: 'task', taskId: 'foreign', objectiveId: 'other', title: 'Foreign', description: 'Other task' }).state;
  assert.throws(() => apply(state, researchRequest({ fromTaskId: 'foreign' })), /same objective/);
  state.objectives[0].status = 'closed';
  assert.throws(() => apply(state, researchRequest()), /closed objective/);
  assert.throws(() => apply(input, researchRequest(), 'mallory'), { code: 'FORBIDDEN' });
  const open = apply(input, researchRequest({ fromTaskId: undefined, to: undefined })).state;
  assert.deepEqual(open.tasks.at(-1).request, { fromTaskId: null, to: 'all' });
});

test('request directory preserves actual status and owner while marking lease expiry and overlap', () => {
  let state = apply(fixture(), researchRequest()).state;
  let directory = requestDirectory(state, { objective: 'problem', now: Date.parse(START) });
  assert.equal(directory.requests.length, 1); assert.equal(directory.requests[0].claimable, true);
  assert.deepEqual(directory.requests[0].requester, { id: 'alice-agent', name: 'Alice', owner: 'alice' });
  assert.deepEqual(directory.requests[0].fromTask, { id: 'lemma', title: 'Prove the lemma' });
  state = apply(state, { kind: 'claim', task: 'boundary-request', agent: 'bob-agent' }, 'bob').state;
  directory = requestDirectory(state, { objective: 'problem', now: Date.parse(START) });
  assert.equal(directory.requests[0].liveLease, true); assert.equal(directory.requests[0].claimable, false);
  assert.deepEqual(directory.requests[0].assigneeAgent, { id: 'bob-agent', name: 'Bob', owner: 'bob' });
  const expired = requestDirectory(state, { objective: 'problem', now: Date.parse('2026-11-25T00:21:00Z') }).requests[0];
  assert.equal(expired.status, 'active'); assert.equal(expired.leaseExpired, true); assert.equal(expired.liveLease, false); assert.equal(expired.claimable, true);
  assert.equal(expired.assigneeAgent.id, 'bob-agent'); // Expiry does not rewrite historical attribution.
  state = apply(state, researchRequest({ taskId: 'overlapping-request', scopePaths: ['research/boundary/other.md'] })).state;
  assert.equal(requestDirectory(state, { objective: 'problem', now: Date.parse(START) }).requests.find(task => task.id === 'overlapping-request').claimable, false);
});

test('legacy tasks/messages are not inferred as requests and malformed metadata fails validation', () => {
  let state = fixture();
  state = apply(state, { kind: 'message', to: 'bob-agent', task: 'lemma', summary: 'Please inspect this' }).state;
  assert.deepEqual(requestDirectory(state, { objective: 'problem', now: Date.parse(START) }).requests, []);
  state = apply(state, researchRequest()).state;
  for (const metadata of [null, [], {}, { fromTaskId: null }, { fromTaskId: 'boundary-request', to: 'all' }, { fromTaskId: null, to: 'missing' }, { fromTaskId: null, to: 'all', hidden: true }]) {
    const changed = structuredClone(state); changed.tasks.at(-1).request = metadata;
    assert.throws(() => validateState(changed));
  }
  const unscoped = structuredClone(state); delete unscoped.tasks.at(-1).scopePaths;
  assert.throws(() => validateState(unscoped), /scope paths/);
  assert.throws(() => requestDirectory(state, { objective: 'missing' }), /Unknown objective/);
  assert.throws(() => requestDirectory(state, { objective: 'problem', now: NaN }), /current time/);
  const copy = requestDirectory(state, { objective: 'problem' }); copy.requests[0].scopePaths.push('changed');
  assert.deepEqual(state.tasks.at(-1).scopePaths, ['research/boundary']);
});

test('request CLI serializes explicit metadata and requests directory performs only a read', async () => {
  let state = fixture();
  const oldFetch = globalThis.fetch, oldLog = console.log, oldError = console.error, oldToken = process.env.GH_TOKEN;
  const calls = [], output = [];
  process.env.GH_TOKEN = 'fixture-token-not-a-credential';
  console.log = value => output.push(String(value)); console.error = () => {};
  globalThis.fetch = async (url, options = {}) => {
    calls.push({ url, method: options.method });
    if (url === 'https://api.github.com/user') return new Response(JSON.stringify({ login: 'alice' }), { status: 200 });
    if (options.method === 'PUT') {
      state = JSON.parse(Buffer.from(JSON.parse(options.body).content, 'base64').toString());
      return new Response(JSON.stringify({ commit: { sha: 'fixture-commit' } }), { status: 200 });
    }
    return new Response(JSON.stringify({ type: 'file', encoding: 'base64', sha: 'fixture-sha', size: Buffer.byteLength(JSON.stringify(state)), content: Buffer.from(JSON.stringify(state)).toString('base64') }), { status: 200, headers: { date: new Date(START).toUTCString() } });
  };
  try {
    await main(['request', '--agent', 'alice-agent', '--id', 'cli-request', '--objective', 'problem', '--title', 'Check endpoints', '--description', 'Exact endpoint calculation', '--paths', 'research/endpoint/', '--from-task', 'lemma', '--to', 'bob-agent', '--request-id', 'stable-cli-request']);
    assert.deepEqual(state.tasks.at(-1).request, { fromTaskId: 'lemma', to: 'bob-agent' });
    assert.deepEqual(state.tasks.at(-1).scopePaths, ['research/endpoint']);
    assert.equal(state.events.at(-1).kind, 'request');
    calls.length = 0; output.length = 0;
    await main(['requests', '--objective', 'problem', '--json']);
    assert.equal(calls.length, 1); assert.equal(calls[0].method, 'GET');
    const directory = JSON.parse(output[0]); assert.equal(directory.requests.length, 1); assert.equal(directory.requests[0].requester.owner, 'alice');
    output.length = 0; await main(['requests', '--objective', 'problem']);
    assert.match(output.join('\n'), /invitation: bob-agent/); assert.match(output.join('\n'), /assigned: unassigned/);
    const count = calls.length;
    await assert.rejects(main(['requests']), /requires --objective/);
    await assert.rejects(main(['request', '--agent', 'alice-agent']), /requires --paths/);
    await assert.rejects(main(['requests', '--objective', 'problem', '--agent', 'alice-agent']), /Unknown option/);
    assert.equal(calls.length, count);
  } finally {
    globalThis.fetch = oldFetch; console.log = oldLog; console.error = oldError;
    if (oldToken === undefined) delete process.env.GH_TOKEN; else process.env.GH_TOKEN = oldToken;
  }
});
