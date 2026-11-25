import { createHash } from 'node:crypto';
import { normalizeUsage, sameUsage, usageRunKey, USAGE_CLOCK_SKEW_MS } from './usage.mjs';

export const MAX_STATE_BYTES = 850_000;
export const TASK_STATUSES = ['open', 'active', 'review', 'verified', 'blocked'];
const slug = /^[a-z0-9][a-z0-9-]{0,79}$/;
const ownerKey = value => String(value).toLowerCase();

export class ProtocolError extends Error {
  constructor(message, code = 'PROTOCOL_ERROR') { super(message); this.name = 'ProtocolError'; this.code = code; }
}
const fail = (message, code) => { throw new ProtocolError(message, code); };
const required = (value, label, max = 12_000) => {
  if (typeof value !== 'string' || !value.trim() || value.length > max) fail(`${label} must be nonempty text, at most ${max} characters.`);
  return value.trim();
};
const identifier = (value, label = 'ID') => {
  if (typeof value !== 'string' || !slug.test(value)) fail(`${label} must use 1–80 lowercase letters, digits, or hyphens, starting with a letter or digit.`);
  return value;
};
const textList = (value, label, max = 30) => {
  if (!Array.isArray(value) || value.length > max) fail(`${label} must be an array with at most ${max} entries.`);
  return [...new Set(value.map(item => required(item, label, 2_000)))];
};
const evidenceList = value => textList(value ?? [], 'Evidence URL').map(value => {
  let url;
  try { url = new URL(value); } catch { fail('Evidence must be an absolute https URL.'); }
  if (url.protocol !== 'https:' || url.username || url.password) fail('Evidence must be an https URL without credentials.');
  return url.href;
});
const fingerprint = op => createHash('sha256').update(JSON.stringify(op, Object.keys(op).sort())).digest('hex');
const isLive = (task, now) => task.status === 'active' && Date.parse(task.leaseUntil) > now;
function normalizePaths(paths) {
  if (!Array.isArray(paths) || !paths.length || paths.length > 30) fail('Declare between 1 and 30 repository-relative paths.');
  return [...new Set(paths.map(value => {
    const path = required(value, 'Scope path', 350).replace(/\/$/, '');
    if (path.split('/').some(part => !/^[\w.-]+$/.test(part) || ['.', '..', '.git'].includes(part.toLowerCase())))
      fail('Scope paths must be literal repository-relative files or directories, without wildcards or traversal.');
    return path;
  }))].sort();
}
function scopeConflict(task, active) {
  const paths = task.scopePaths || [];
  return active.find(other => other.id !== task.id && (other.scopePaths || []).some(scope =>
    paths.some(path => {
      const a = scope.toLowerCase(), b = path.toLowerCase();
      return a === b || a.startsWith(b + '/') || b.startsWith(a + '/');
    })));
}
const findTask = (state, id) => state.tasks.find(task => task.id === id) ?? fail(`Unknown task: ${id}.`);
const findObjective = (state, id) => state.objectives.find(item => item.id === id) ?? fail(`Unknown objective: ${id}.`);
const findAgent = (state, id) => state.agents.find(agent => agent.id === id) ?? fail(`Unknown agent: ${id}. Join first.`);

export function validateState(state) {
  if (!state || state.schemaVersion !== 1 || !Number.isSafeInteger(state.revision) || state.revision < 0) fail('Unsupported or invalid coordination state.');
  if (!Number.isFinite(Date.parse(state.updatedAt))) fail('Invalid state update timestamp.');
  for (const key of ['objectives', 'tasks', 'agents', 'events']) {
    if (!Array.isArray(state[key])) fail(`State ${key} must be an array.`);
    const ids = new Set();
    for (const item of state[key]) {
      if (!item || typeof item.id !== 'string' || ids.has(item.id)) fail(`Invalid or duplicate ID in ${key}.`);
      ids.add(item.id);
    }
  }
  for (const agent of state.agents) {
    identifier(agent.id, 'Agent ID'); required(agent.owner, 'Agent owner', 100);
  }
  for (const objective of state.objectives) {
    identifier(objective.id, 'Objective ID');
    if (!['open', 'closed'].includes(objective.status)) fail('Invalid objective status.');
  }
  for (const task of state.tasks) {
    identifier(task.id, 'Task ID'); findObjective(state, task.objectiveId);
    if (task.scopePaths !== undefined) {
      const normalized = normalizePaths(task.scopePaths);
      if (JSON.stringify(normalized) !== JSON.stringify(task.scopePaths)) fail('Stored scope paths must be canonical and unique.');
    }
    if (task.request !== undefined) {
      const request = task.request;
      if (!request || typeof request !== 'object' || Array.isArray(request) ||
          Object.keys(request).some(key => !['fromTaskId', 'to'].includes(key)) ||
          !Object.hasOwn(request, 'fromTaskId') || !Object.hasOwn(request, 'to')) fail('Invalid research request metadata.');
      if (!task.scopePaths) fail('Research requests require declared scope paths.');
      findAgent(state, task.createdBy);
      if (!Number.isFinite(Date.parse(task.createdAt))) fail('Research requests require a creation timestamp.');
      if (request.to !== 'all') findAgent(state, identifier(request.to, 'Request target'));
      if (request.fromTaskId !== null) {
        const origin = findTask(state, identifier(request.fromTaskId, 'Origin task ID'));
        if (origin.objectiveId !== task.objectiveId || origin.id === task.id) fail('A request origin must be another task in the same objective.');
      }
    }
    if (!TASK_STATUSES.includes(task.status) || !Array.isArray(task.dependsOn) || !Array.isArray(task.checkpoints) || !Array.isArray(task.reviews)) fail(`Invalid task: ${task.id}.`);
    if (task.assignee != null) findAgent(state, task.assignee);
    if (task.status === 'active' && (!task.assignee || !Number.isFinite(Date.parse(task.leaseUntil)))) fail(`Active task ${task.id} requires an assignee and valid lease.`);
    for (const dependency of task.dependsOn) {
      const dep = findTask(state, dependency);
      if (dep.objectiveId !== task.objectiveId || dep.id === task.id) fail(`Invalid dependency on ${dependency}.`);
    }
    const checkpointIds = new Set();
    for (const checkpoint of task.checkpoints) {
      if (!checkpoint?.id || checkpointIds.has(checkpoint.id)) fail(`Invalid checkpoint in ${task.id}.`);
      checkpointIds.add(checkpoint.id); findAgent(state, checkpoint.agentId);
      evidenceList(checkpoint.evidence);
    }
    for (const review of task.reviews) {
      findAgent(state, review.agentId);
      if (!checkpointIds.has(review.checkpointId) || !['pass', 'changes'].includes(review.verdict)) fail(`Invalid review in ${task.id}.`);
    }
  }
  if (state.usage !== undefined) {
    if (!Array.isArray(state.usage)) fail('State usage must be an array.');
    const ids = new Set(), runs = new Set();
    for (const receipt of state.usage) {
      if (!receipt || typeof receipt.id !== 'string' || ids.has(receipt.id)) fail('Invalid or duplicate usage receipt ID.');
      const normalized = normalizeUsage(receipt);
      const agent = findAgent(state, receipt.agentId), task = findTask(state, receipt.taskId);
      if (ownerKey(agent.owner) !== ownerKey(receipt.owner) || receipt.objectiveId !== task.objectiveId || !Number.isFinite(Date.parse(receipt.recordedAt))) fail('Usage receipt identity, objective, or timestamp is invalid.');
      if (receipt.totalTokens !== normalized.totalTokens) fail('Usage receipt total is invalid.');
      if (Date.parse(normalized.finishedAt) > Date.parse(receipt.recordedAt) + USAGE_CLOCK_SKEW_MS) fail('Usage finishedAt exceeds recordedAt by more than the five-minute clock-skew allowance.');
      const key = usageRunKey(receipt.owner, normalized);
      if (runs.has(key)) fail('A runtime run was counted more than once.');
      ids.add(receipt.id); runs.add(key);
    }
  }
  const visited = new Set(), visiting = new Set();
  function visit(id) {
    if (visiting.has(id)) fail('Task dependency cycle detected.');
    if (visited.has(id)) return;
    visiting.add(id);
    findTask(state, id).dependsOn.forEach(visit);
    visiting.delete(id); visited.add(id);
  }
  state.tasks.forEach(task => visit(task.id));
  if (Buffer.byteLength(JSON.stringify(state, null, 2) + '\n', 'utf8') > MAX_STATE_BYTES) fail('Coordination state is near the GitHub Contents API limit. Archive reviewed history in a maintainer PR before accepting more work.', 'STATE_CAPACITY');
  return state;
}

export function eligibleTasks(state, { now = Date.now(), approach, objective } = {}) {
  const active = state.tasks.filter(task => isLive(task, now));
  return state.tasks.filter(task =>
    (!objective || task.objectiveId === objective) &&
    (!approach || task.approach === approach) &&
    findObjective(state, task.objectiveId).status === 'open' &&
    (task.status === 'open' || (task.status === 'active' && !isLive(task, now))) &&
    task.dependsOn.every(id => findTask(state, id).status === 'verified') &&
    !scopeConflict(task, active)
  );
}

/** Explicit requests remain ordinary tasks; invitation targets never restrict claims. */
export function requestDirectory(state, { objective, now = Date.now() } = {}) {
  validateState(state);
  findObjective(state, identifier(objective, 'Objective ID'));
  if (!Number.isFinite(now)) fail('Request directory requires a valid current time.');
  const eligible = new Set(eligibleTasks(state, { objective, now }).map(task => task.id));
  const identity = id => {
    if (!id) return null;
    const agent = findAgent(state, id);
    return { id: agent.id, name: agent.name || agent.id, owner: agent.owner };
  };
  return { revision: state.revision, updatedAt: state.updatedAt, objectiveId: objective,
    requests: state.tasks.filter(task => task.objectiveId === objective && task.request !== undefined).map(task => {
      const origin = task.request.fromTaskId ? findTask(state, task.request.fromTaskId) : null;
      return { ...structuredClone(task), requester: identity(task.createdBy), assigneeAgent: identity(task.assignee),
        fromTask: origin ? { id: origin.id, title: origin.title } : null,
        liveLease: isLive(task, now), leaseExpired: task.status === 'active' && !isLive(task, now), claimable: eligible.has(task.id) };
    }) };
}

/** Pure state transition. All transport retries replay the same operation ID. */
export function applyOperation(input, operation, { owner, now = new Date().toISOString(), leaseMinutes = 20 } = {}) {
  validateState(input);
  required(owner, 'Authenticated GitHub owner', 100);
  required(operation?.id, 'Request ID', 120);
  const at = new Date(now).toISOString(), millis = Date.parse(at);
  if (!Number.isInteger(leaseMinutes) || leaseMinutes < 1 || leaseMinutes > 120) fail('Lease duration must be 1–120 minutes.');
  const digest = fingerprint(operation);
  const existing = input.events.find(event => event.id === operation.id);
  if (existing) {
    if (ownerKey(existing.actorOwner) !== ownerKey(owner) || existing.requestFingerprint !== digest) fail('Request ID was already used for a different operation.', 'IDEMPOTENCY_CONFLICT');
    return { state: input, event: existing, duplicate: true };
  }
  const state = structuredClone(input), op = operation;
  const agentId = identifier(op.agent, 'Agent ID');
  let actor = state.agents.find(agent => agent.id === agentId);
  if (actor && ownerKey(actor.owner) !== ownerKey(owner)) fail(`Agent ${agentId} belongs to another GitHub account.`, 'FORBIDDEN');
  if (op.kind !== 'join' && !actor) fail('Join with this agent ID before publishing work.');
  let task, objectiveId = null, taskId = null, summary = '', detail = '', links = [], to = null;
  const leaseUntil = () => new Date(millis + leaseMinutes * 60_000).toISOString();
  const ownedLease = () => {
    if (task.assignee !== agentId || !isLive(task, millis)) fail('This operation requires your active, unexpired task lease. Claim the task again if its lease expired.', 'LEASE_REQUIRED');
  };
  if (op.task) { task = findTask(state, op.task); taskId = task.id; objectiveId = task.objectiveId; }
  switch (op.kind) {
    case 'join': {
      const name = required(op.name ?? agentId, 'Agent name', 160);
      const role = required(op.role ?? 'explorer', 'Agent role', 100);
      const capabilities = textList(op.capabilities ?? [], 'Capability', 30);
      if (!actor) {
        actor = { id: agentId, owner, name, role, capabilities, registeredAt: at, lastSeenAt: at };
        state.agents.push(actor);
      } else Object.assign(actor, { name, role, capabilities });
      summary = `${name} joined as ${role}.`;
      break;
    }
    case 'objective': {
      const id = identifier(op.objectiveId, 'Objective ID');
      if (state.objectives.some(item => item.id === id)) fail(`Objective ${id} already exists.`);
      state.objectives.push({ id, title: required(op.title, 'Title', 240), statement: required(op.statement, 'Problem statement'), successCriteria: required(op.success, 'Success criteria'), status: 'open', tags: textList(op.tags ?? [], 'Tag', 12), example: false, createdAt: at, createdBy: agentId, createdByOwner: owner });
      objectiveId = id; summary = `Opened objective: ${op.title}`;
      break;
    }
    case 'task':
    case 'request': {
      const id = identifier(op.taskId, 'Task ID'), objective = findObjective(state, op.objectiveId);
      if (objective.status !== 'open') fail('Cannot add work to a closed objective.');
      if (state.tasks.some(item => item.id === id)) fail(`Task ${id} already exists.`);
      const dependsOn = textList(op.dependsOn ?? [], 'Dependency', 30);
      for (const id of dependsOn) {
        const dependency = findTask(state, id);
        if (dependency.objectiveId !== objective.id) fail('Dependencies must belong to the same objective.');
      }
      let request;
      if (op.kind === 'request') {
        const fromTaskId = op.fromTaskId == null ? null : identifier(op.fromTaskId, 'Origin task ID');
        if (fromTaskId && findTask(state, fromTaskId).objectiveId !== objective.id) fail('A request origin must belong to the same objective.');
        to = op.to ?? 'all';
        if (to !== 'all') findAgent(state, identifier(to, 'Request target'));
        request = { fromTaskId, to };
        normalizePaths(op.scopePaths);
      }
      state.tasks.push({ id, objectiveId: objective.id, title: required(op.title, 'Task title', 240), description: required(op.description, 'Task description'), approach: required(op.approach ?? 'constructive', 'Approach', 100), status: 'open', dependsOn, assignee: null, leaseUntil: null, createdAt: at, updatedAt: at, createdBy: agentId, checkpoints: [], reviews: [] });
      if (op.scopePaths !== undefined) state.tasks.at(-1).scopePaths = normalizePaths(op.scopePaths);
      if (request) state.tasks.at(-1).request = request;
      taskId = id; objectiveId = objective.id; summary = `Proposed ${op.approach ?? 'constructive'} task: ${op.title}`;
      if (request) { summary = `Requested research: ${state.tasks.at(-1).title}`; detail = state.tasks.at(-1).description; }
      break;
    }
    case 'claim': {
      if (!task) fail('A task is required.');
      if (isLive(task, millis)) fail(`Task is leased to ${task.assignee} until ${task.leaseUntil}.`, 'TASK_BUSY');
      if (op.scopePaths !== undefined) {
        const paths = normalizePaths(op.scopePaths);
        if (task.scopePaths && JSON.stringify(paths) !== JSON.stringify(task.scopePaths)) fail('Claim must retain the declared scope paths.');
        task.scopePaths = paths;
      }
      const conflict = scopeConflict(task, state.tasks.filter(item => isLive(item, millis)));
      if (conflict) fail(`Overlapping work is leased on ${conflict.id} by ${conflict.assignee} until ${conflict.leaseUntil}.`, 'TASK_BUSY');
      if (!eligibleTasks(state, { now: millis }).some(item => item.id === task.id)) fail('Task is not claimable: check its status, objective, and verified dependencies.', 'NOT_ELIGIBLE');
      const active = state.tasks.find(item => item.assignee === agentId && isLive(item, millis));
      if (active) fail(`Finish or release ${active.id} before claiming another task.`, 'AGENT_BUSY');
      const previous = task.assignee;
      Object.assign(task, { status: 'active', assignee: agentId, leaseUntil: leaseUntil() });
      summary = previous ? `Reclaimed task after ${previous}'s lease expired.` : 'Claimed task.';
      break;
    }
    case 'heartbeat':
      if (!task) fail('A task is required.'); ownedLease();
      task.leaseUntil = leaseUntil(); summary = 'Renewed task lease.'; break;
    case 'checkpoint': {
      if (!task) fail('A task is required.'); ownedLease();
      summary = required(op.summary, 'Checkpoint summary', 2_000); detail = required(op.details ?? op.summary, 'Checkpoint details', 24_000); links = evidenceList(op.evidence);
      task.checkpoints.push({ id: op.id, agentId, summary, details: detail, evidence: links, createdAt: at });
      task.leaseUntil = leaseUntil(); break;
    }
    case 'submit': {
      if (!task) fail('A task is required.'); ownedLease();
      const latest = task.checkpoints.at(-1);
      if (!latest || latest.agentId !== agentId || !latest.evidence.length) fail('Publish your checkpoint with evidence links before submitting.');
      summary = required(op.summary, 'Submission summary', 2_000);
      Object.assign(task, { status: 'review', leaseUntil: null, submittedCheckpointId: latest.id });
      links = latest.evidence; break;
    }
    case 'review': {
      if (!task || task.status !== 'review') fail('Only a task awaiting review can be reviewed.');
      const latest = task.checkpoints.at(-1);
      if (!latest || op.checkpointId !== latest.id || task.submittedCheckpointId !== latest.id) fail('Review must name the exact currently submitted checkpoint ID.', 'STALE_REVIEW');
      const author = findAgent(state, latest.agentId);
      if (task.checkpoints.some(checkpoint => ownerKey(findAgent(state, checkpoint.agentId).owner) === ownerKey(owner)) || ownerKey(author.owner) === ownerKey(owner)) fail('Independent review requires a different GitHub owner from every checkpoint contributor.', 'INDEPENDENT_REVIEW_REQUIRED');
      if (!['pass', 'changes'].includes(op.verdict)) fail('Review verdict must be pass or changes.');
      summary = required(op.summary, 'Review summary', 6_000); links = evidenceList(op.evidence);
      if (!links.length) fail('Review requires reproducible evidence URLs.');
      task.reviews.push({ id: op.id, agentId, verdict: op.verdict, summary, evidence: links, checkpointId: latest.id, createdAt: at });
      Object.assign(task, { status: op.verdict === 'pass' ? 'verified' : 'open', leaseUntil: null, assignee: op.verdict === 'pass' ? task.assignee : null });
      summary = `${op.verdict === 'pass' ? 'Independent review passed' : 'Changes requested'}: ${summary}`;
      break;
    }
    case 'release':
    case 'block': {
      if (!task) fail('A task is required.'); ownedLease();
      summary = required(op.summary, op.kind === 'release' ? 'Handoff summary' : 'Blocker summary', 6_000);
      Object.assign(task, { status: op.kind === 'release' ? 'open' : 'blocked', assignee: null, leaseUntil: null });
      break;
    }
    case 'reopen': {
      if (!task || !['blocked', 'review', 'verified'].includes(task.status)) fail('Only blocked, review, or verified tasks can be reopened.');
      summary = required(op.summary, 'Reason to reopen', 6_000);
      const downstream = state.tasks.filter(item => item.dependsOn.includes(task.id) && ['active', 'review', 'verified'].includes(item.status));
      if (downstream.length) fail(`Reopen downstream tasks first: ${downstream.map(item => item.id).join(', ')}.`, 'DEPENDENTS_ACTIVE');
      if (findObjective(state, task.objectiveId).status !== 'open') fail('Reopen the objective before changing its tasks.');
      Object.assign(task, { status: 'open', assignee: null, leaseUntil: null, submittedCheckpointId: null });
      break;
    }
    case 'usage': {
      if (!task) fail('Usage requires a task.');
      const contributed = task.checkpoints.some(checkpoint => checkpoint.agentId === agentId) || state.events.some(event => event.kind === 'claim' && event.actor === agentId && event.taskId === task.id && ownerKey(event.actorOwner) === ownerKey(owner));
      if (!contributed) fail('Usage requires this agent to have claimed or checkpointed the task under its authenticated owner.', 'USAGE_CONTRIBUTION_REQUIRED');
      const receipt = normalizeUsage(op);
      if (Date.parse(receipt.finishedAt) > millis + USAGE_CLOCK_SKEW_MS) fail('Usage finishedAt is more than five minutes ahead of GitHub server time.', 'USAGE_FUTURE_TIMESTAMP');
      const key = usageRunKey(owner, receipt);
      const previous = (state.usage ?? []).find(item => usageRunKey(item.owner, item) === key);
      if (previous) {
        if (previous.taskId !== task.id || !sameUsage(previous, receipt)) fail('This owner/provider/run ID already has a different immutable usage receipt.', 'USAGE_RUN_CONFLICT');
        return { state: input, event: input.events.find(event => event.id === previous.id), usage: previous, duplicate: true };
      }
      state.usage ??= [];
      state.usage.push({ id: op.id, agentId, owner, taskId: task.id, objectiveId: task.objectiveId, ...receipt, recordedAt: at });
      summary = `Reported ${receipt.totalTokens} tokens for ${receipt.provider}/${receipt.model} (${receipt.source}; contributor-reported).`;
      detail = `Input: ${receipt.inputTokens}; output: ${receipt.outputTokens}; cached input subset: ${receipt.cachedInputTokens}. Run: ${receipt.runId}.`;
      links = receipt.evidence ?? [];
      break;
    }
    case 'message': {
      summary = required(op.summary, 'Message summary', 2_000); detail = required(op.details ?? op.summary, 'Message details', 12_000);
      to = op.to ?? 'all'; if (to !== 'all') findAgent(state, to);
      links = evidenceList(op.evidence); break;
    }
    case 'close': {
      const objective = findObjective(state, op.objectiveId);
      const creatorOwner = objective.createdByOwner ?? findAgent(state, objective.createdBy).owner;
      if (ownerKey(creatorOwner) !== ownerKey(owner)) fail('Only the objective creator’s GitHub owner may close it.', 'FORBIDDEN');
      const tasks = state.tasks.filter(item => item.objectiveId === objective.id);
      if (!tasks.length || tasks.some(item => item.status !== 'verified')) fail('Every task must be independently reviewed before closing an objective.');
      summary = required(op.summary, 'Synthesis conclusion', 12_000); links = evidenceList(op.evidence);
      if (!links.length) fail('Closing an objective requires a synthesis artifact URL.');
      objective.status = 'closed'; objective.conclusion = summary; objective.evidence = links; objective.closedAt = at; objectiveId = objective.id; break;
    }
    case 'reopen-objective': {
      const objective = findObjective(state, op.objectiveId);
      summary = required(op.summary, 'Reason to reopen', 6_000);
      objective.status = 'open'; objectiveId = objective.id; break;
    }
    default: fail(`Unknown operation: ${op.kind}.`);
  }
  actor.lastSeenAt = at;
  if (task && op.kind !== 'usage') task.updatedAt = at;
  const event = { id: op.id, kind: op.kind, actor: agentId, actorOwner: owner, objectiveId, taskId, summary, detail, to, createdAt: at, links, requestFingerprint: digest };
  if (op.kind === 'request') {
    const requested = findTask(state, taskId);
    event.scopePaths = [...requested.scopePaths]; event.fromTaskId = requested.request.fromTaskId;
  }
  state.events.push(event); state.revision += 1; state.updatedAt = at;
  validateState(state);
  return { state, event, duplicate: false, ...(op.kind === 'usage' ? { usage: state.usage.at(-1) } : {}) };
}
