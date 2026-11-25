import test from 'node:test';
import assert from 'node:assert/strict';
import { buildResearchTree, taskStatus } from '../public/tree-model.js';

const now = Date.parse('2026-11-25T12:00:00Z');
const task = (id, change = {}) => ({ id, objectiveId: 'math', title: id, approach: 'constructive', status: 'open',
  createdAt: '2026-11-25T10:00:00Z', dependsOn: [], checkpoints: [], reviews: [], ...change });
const agent = (id, name = id) => ({ id, name, owner: 'owner' });
const event = (id, change = {}) => ({ id, objectiveId: 'math', taskId: null, actor: 'author', createdAt: '2026-11-25T11:00:00Z', ...change });
const state = (change = {}) => ({ objectives: [{ id: 'math', title: 'A conjecture', createdBy: 'creator' }, { id: 'platform' }],
  tasks: [], agents: [], events: [], ...change });
const ids = entries => entries.map(entry => entry.id);

test('tasks and events stay within one objective, including legacy task-only events', () => {
  const input = state({ tasks: [task('proof'), task('site', { objectiveId: 'platform' })], events: [
    event('objective'), event('task', { taskId: 'proof' }), event('legacy', { objectiveId: null, taskId: 'proof' }),
    event('global', { objectiveId: null }), event('other', { objectiveId: 'platform' }),
    event('conflicting-objective', { objectiveId: 'platform', taskId: 'proof' }),
    event('conflicting-task', { taskId: 'site' }), event('unknown-task', { taskId: 'missing' })
  ] });
  const tree = buildResearchTree(input, 'math', { now });
  assert.equal(tree.objective, input.objectives[0]);
  assert.deepEqual(ids(tree.tasks), ['proof']);
  assert.deepEqual(ids(tree.events), ['legacy', 'task', 'objective']);
  assert.equal(tree.events[0], input.events[2]);
});

test('branches use approach strings and predictable natural order; tasks use stable creation order', () => {
  const input = state({ tasks: [task('t10', { title: 'Lemma 10', approach: 'direction-10' }),
    task('t2', { title: 'Lemma 2', approach: 'direction-2' }),
    task('later', { approach: 'direction-2', createdAt: '2026-11-25T10:30:00Z' }),
    task('earlier', { approach: 'direction-2', createdAt: '2026-11-25T09:00:00Z', updatedAt: '2026-11-25T11:59:59Z' })] });
  const tree = buildResearchTree(input, 'math', { now });
  assert.deepEqual(tree.branches.map(branch => [branch.id, branch.label]), [['direction-2', 'Direction 2'], ['direction-10', 'Direction 10']]);
  assert.deepEqual(ids(tree.tasks), ['earlier', 't2', 't10', 'later']);
  assert.deepEqual(ids(tree.branches[0].tasks), ['earlier', 't2', 'later']);
  assert.deepEqual(buildResearchTree({ ...input, tasks: [...input.tasks].reverse() }, 'math', { now }), tree);
});

test('lease display expires at its exact boundary and never promotes another status to active', () => {
  const live = task('live', { status: 'active', assignee: 'worker', leaseUntil: '2026-11-25T12:00:00.001Z' });
  assert.equal(taskStatus(live, now), 'active');
  assert.equal(taskStatus(live, now + 1), 'expired');
  for (const leaseUntil of [null, 'invalid', '2026-11-25T11:59:59Z']) assert.equal(taskStatus({ ...live, leaseUntil }, now), 'expired');
  assert.equal(taskStatus({ ...live, assignee: null }, now), 'expired');
  for (const status of ['open', 'review', 'verified', 'blocked', 'closed']) assert.equal(taskStatus({ ...live, status }, now), status);
  assert.equal(taskStatus(null, now), 'unknown');
  const tree = buildResearchTree(state({ tasks: [live, { ...live, id: 'expired', leaseUntil: '2026-11-25T12:00:00Z' }, task('review', { status: 'review' })] }), 'math', { now });
  assert.equal(tree.branches[0].activeCount, 1);
});

test('registered participants include creators, assignments, contribution authors and explicit message recipients only', () => {
  const registered = ['creator', 'proposer', 'worker', 'checkpoint-author', 'reviewer', 'author', 'recipient', 'unrelated', 'global-actor'].map(id => agent(id));
  const input = state({ agents: registered, tasks: [task('proof', { assignee: 'worker', createdBy: 'proposer',
    checkpoints: [{ agentId: 'checkpoint-author' }], reviews: [{ agentId: 'reviewer' }] })], events: [
    event('message', { to: 'recipient' }), event('broadcast', { to: 'all' }), event('unknown-agent', { actor: 'not-registered' }),
    event('global', { objectiveId: null, actor: 'global-actor', to: 'unrelated' })] });
  const tree = buildResearchTree(input, 'math', { now });
  assert.deepEqual(ids(tree.agents), ['author', 'checkpoint-author', 'creator', 'proposer', 'recipient', 'reviewer', 'worker']);
  assert.ok(tree.agents.every(row => registered.includes(row)));
  assert.deepEqual(buildResearchTree(input, 'platform', { now }).agents, []);
});

test('event ordering follows timestamp then reverse ledger order; malformed dates are last', () => {
  const input = state({ events: [event('new', { createdAt: '2026-11-25T11:30:00Z' }),
    event('old'), event('tie-first', { createdAt: '2026-11-25T11:30:00Z' }), event('tie-last', { createdAt: '2026-11-25T11:30:00Z' }),
    event('missing-date', { createdAt: undefined })] });
  assert.deepEqual(ids(buildResearchTree(input, 'math', { now }).events), ['tie-last', 'tie-first', 'new', 'old', 'missing-date']);
});

test('unknown objective and missing collections stay empty without leaking a default objective', () => {
  const empty = { objective: null, branches: [], tasks: [], agents: [], events: [] };
  assert.deepEqual(buildResearchTree(state(), 'missing', { now }), empty);
  assert.deepEqual(buildResearchTree(undefined, 'math', { now }), empty);
  assert.deepEqual(buildResearchTree(state(), undefined, { now }), empty);
  const objective = { id: 'math', status: 'closed' };
  assert.deepEqual(buildResearchTree({ objectives: [objective] }, 'math', { now }), { ...empty, objective });
  assert.throws(() => buildResearchTree(state(), 'math', { now: NaN }), /finite current time/);
  assert.throws(() => taskStatus(task('proof'), Infinity), /finite current time/);
});

test('raw records, dependency arrays and source ordering are preserved without fabricated edges', () => {
  const first = task('first'), second = task('second', { dependsOn: ['first'], approach: 'verification' });
  const input = state({ tasks: [second, first], agents: [agent('creator')], events: [event('created', { actor: 'creator' })] });
  function freeze(value) { if (value && typeof value === 'object') { Object.freeze(value); Object.values(value).forEach(freeze); } return value; }
  freeze(input);
  const before = JSON.stringify(input), tree = buildResearchTree(input, 'math', { now });
  assert.equal(JSON.stringify(input), before);
  assert.equal(tree.tasks[0], first);
  assert.equal(tree.tasks[1], second);
  assert.equal(tree.tasks[1].dependsOn, second.dependsOn);
  assert.deepEqual(tree.tasks[1].dependsOn, ['first']);
  assert.equal('edges' in tree, false);
  assert.equal(tree.branches.find(branch => branch.id === 'verification').tasks[0], second);
});
