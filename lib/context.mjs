import { eligibleTasks } from './core.mjs';

const PAGE_SIZE = 8;
const noise = new Set(['heartbeat', 'usage', 'join']);
const decisions = new Set(['block', 'reopen', 'release', 'review', 'close', 'reopen-objective']);
const cli = 'node bin/commons.mjs context';

function excerpt(text, limit) {
  const marker = '… [abridged]';
  return typeof text === 'string' && text.length > limit ? text.slice(0, limit - marker.length) + marker : text;
}

export function objectiveDecision(events, objectiveId) {
  const event = events.findLast(item => item.objectiveId === objectiveId && ['close', 'reopen-objective'].includes(item.kind));
  return event && eventCard(event);
}

function objectiveCard(objective, decision) {
  if (!objective) return null;
  const { id, title, status, example, statement, successCriteria } = objective;
  const card = { id, title, status, example, statement: excerpt(statement, 480), successCriteria: excerpt(successCriteria, 320) };
  if (objective.conclusion) {
    const field = status === 'closed' ? 'synthesis' : 'previousSynthesis';
    card[field] = { summary: excerpt(objective.conclusion, 240), closedAt: objective.closedAt, ...evidencePreview(objective.evidence) };
  }
  if (decision) card.latestDecision = eventCard(decision);
  if (objective.conclusion || status === 'closed' || statement?.length > 480 || successCriteria?.length > 320) {
    card.exactContract = `${cli} --objective ${id} --exact`;
  }
  return card;
}

function checkpointCard(checkpoint) {
  if (!checkpoint) return null;
  return { id: checkpoint.id, summary: excerpt(checkpoint.summary, 240), ...evidencePreview(checkpoint.evidence, checkpoint.evidenceCount) };
}

function evidencePreview(links = [], count = links.length) {
  return { evidenceCount: count, evidence: links.filter(url => typeof url === 'string' && url.length <= 320).slice(0, 1) };
}

function taskCard(task, now) {
  const { id, title, status, assignee, leaseUntil } = task;
  const card = { id, title, status, assignee };
  if (status === 'active') {
    card.leaseUntil = leaseUntil;
    if (Number.isFinite(now)) card.leaseExpired = Date.parse(leaseUntil) <= now;
  }
  return card;
}

function taskContext(task, now) {
  const { objectiveId, description, approach, scopePaths, dependsOn, request, submittedCheckpointId } = task;
  const review = task.reviews?.at(-1);
  return {
    ...taskCard(task, now), objectiveId, description, approach, scopePaths, dependsOn, request, submittedCheckpointId,
    latestCheckpoint: checkpointCard(task.checkpoints?.at(-1)),
    latestReview: review ? { id: review.id, verdict: review.verdict, checkpointId: review.checkpointId, summary: excerpt(review.summary, 240) } : null
  };
}

function eventCard({ id, kind, taskId, actor, to, summary, createdAt }) {
  return { id, kind, taskId, actor, ...(to ? { to } : {}), summary: excerpt(summary, 180), createdAt };
}

function pageOf(items, number) {
  const pages = Math.max(1, Math.ceil(items.length / PAGE_SIZE));
  if (!Number.isSafeInteger(number) || number < 1 || number > pages) throw new Error(`Page must be an integer from 1 to ${pages}.`);
  return { page: number, pages, total: items.length, items: items.slice((number - 1) * PAGE_SIZE, number * PAGE_SIZE) };
}

// Show lagging decisions and direct messages even if routine checkpoints have filled the recent page.
export function taskSignals(events, task, agent, page = 1) {
  const objectiveId = task.objectiveId;
  const related = new Set([task.id, ...(task.dependsOn || [])]);
  if (task.request?.fromTaskId) related.add(task.request.fromTaskId);
  const relevant = events.filter(event => {
    if (noise.has(event.kind) || (event.objectiveId && event.objectiveId !== objectiveId)) return false;
    return related.has(event.taskId) || (event.kind === 'message' && ((agent && event.to === agent) || event.to === 'all')) ||
      (!event.taskId && event.objectiveId === objectiveId && decisions.has(event.kind));
  }).reverse();
  const window = pageOf(relevant, page);
  const direct = agent && relevant.find(event => event.kind === 'message' && event.to === agent);
  const decision = relevant.find(event => decisions.has(event.kind));
  const attention = [direct, decision].filter(event => event && !window.items.includes(event));
  return { ...window, items: window.items.map(eventCard), ...(attention.length ? { earlierSignals: attention.map(eventCard) } : {}) };
}

function sourceCard(source) {
  if (!source) return null;
  return { id: source.id, title: source.title, latestCheckpoint: checkpointCard(source.latestCheckpoint) };
}

function dependencyCards(dependencies, now) {
  return { dependencyCount: dependencies.length, dependencies: dependencies.slice(0, PAGE_SIZE).map(dep => taskCard(dep, now)) };
}

export function modelContext(packet) {
  const dependencies = packet.dependencies || [];
  const now = Date.parse(packet.snapshotAt);
  const decision = packet.objectiveDecision || objectiveDecision(packet.recentEvents || [], packet.objective?.id);
  const signals = packet.contextSignals || taskSignals(packet.recentEvents || [], packet.task, packet.agent.id);
  // Project optional v1 extensions too: an upstream adapter may include additional historical fields.
  const { page, pages, total } = signals;
  return {
    repository: packet.repository, agent: { id: packet.agent.id, owner: packet.agent.owner, role: packet.agent.role },
    revision: packet.revision, snapshotAt: packet.snapshotAt,
    objective: objectiveCard(packet.objective, decision), task: taskContext(packet.task, now),
    ...dependencyCards(dependencies, now),
    sourceTask: sourceCard(packet.sourceTask),
    signals: { page, pages, total, items: signals.items.slice(0, PAGE_SIZE).map(eventCard),
      ...(signals.earlierSignals?.length ? { earlierSignals: signals.earlierSignals.slice(0, 2).map(eventCard) } : {}) }
  };
}

// The CLI fetches the validated GitHub snapshot once. Only the selected projection reaches stdout.
export function liveContext(state, flags, now, { replaySeedOwners } = {}) {
  const { objective, task, event, agent, exact } = flags;
  if ([objective, task, event].filter(Boolean).length !== 1) throw new Error('Select exactly one of --objective ID, --task ID, or --event ID.');
  if (exact && (!objective || flags.page || agent)) throw new Error('--exact selects the objective contract alone.');
  if (event && (flags.page || agent)) throw new Error('--event selects one exact record; omit --page and --agent.');
  const page = Number(flags.page ?? 1);
  const revision = state.revision;
  if (event) {
    const record = state.events.find(item => item.id === event);
    // Older seeds can contain a checkpoint without a matching event; reviews store the exact checkpoint only on the task.
    if (record) {
      const selected = state.tasks.find(item => item.id === record.taskId);
      const review = selected?.reviews.find(item => item.id === event);
      return { revision, event: record, ...(review ? { verdict: review.verdict, checkpointId: review.checkpointId } : {}) };
    }
    let match;
    for (const item of state.tasks) {
      const checkpoint = item.checkpoints.find(entry => entry.id === event);
      const review = item.reviews.find(entry => entry.id === event);
      if (!checkpoint && !review) continue;
      if (match || (checkpoint && review)) throw new Error('Ambiguous legacy id. Use packet --task ID and inspect that task’s matching entry.');
      match = { revision, taskId: item.id, ...(checkpoint ? { checkpoint } : { review }) };
    }
    if (match) return match;
    throw new Error(`Unknown event or checkpoint: ${event}.`);
  }
  const owner = agent && state.agents.find(item => item.id === agent)?.owner;
  if (agent && !owner) throw new Error(`Unknown agent: ${agent}. Join first or omit --agent.`);
  if (task) {
    const selected = state.tasks.find(item => item.id === task);
    if (!selected) throw new Error(`Unknown task: ${task}.`);
    const obj = state.objectives.find(item => item.id === selected.objectiveId);
    return {
      revision, objective: objectiveCard(obj, objectiveDecision(state.events, obj.id)), task: taskContext(selected, now),
      claimable: eligibleTasks(state, { now, objective: obj.id, owner, replaySeedOwners }).some(item => item.id === task),
      ...(replaySeedOwners && Object.hasOwn(replaySeedOwners, task) ? { replaySeedOwners: replaySeedOwners[task] } : {}),
      ...dependencyCards(state.tasks.filter(item => selected.dependsOn.includes(item.id)), now),
      signals: taskSignals(state.events, selected, agent, page),
      expand: 'Read one event with context --event ID. Use --page N for additional signals; packet --task ID exports the complete task history.'
    };
  }
  const selected = state.objectives.find(item => item.id === objective);
  if (!selected) throw new Error(`Unknown objective: ${objective}.`);
  if (exact) return { revision, objective: selected };
  const tasks = state.tasks.filter(item => item.objectiveId === objective);
  const statusCounts = {};
  for (const item of tasks) statusCounts[item.status] = (statusCounts[item.status] || 0) + 1;
  const eligible = new Set(eligibleTasks(state, { now, objective, owner, replaySeedOwners }).map(item => item.id));
  function priority(item) {
    if (item.status === 'active' && Date.parse(item.leaseUntil) > now) return 0;
    if (eligible.has(item.id)) return 1;
    if (item.status === 'review') return 2;
    if (item.status === 'verified') return 4;
    return 3;
  }
  tasks.sort((a, b) => priority(a) - priority(b) || (b.updatedAt || '').localeCompare(a.updatedAt || '') || a.id.localeCompare(b.id));
  const window = pageOf(tasks, page);
  const inbox = state.events.filter(item => item.kind === 'message' && (!item.objectiveId || item.objectiveId === objective) &&
    (item.to === 'all' || (agent && item.to === agent)));
  const direct = agent && inbox.findLast(item => item.to === agent);
  const latest = inbox.slice(-3);
  const cards = window.items.map(item => ({
    ...taskCard(item, now), claimable: eligible.has(item.id), latestFinding: excerpt(item.checkpoints.at(-1)?.summary, 130)
  }));
  return {
    revision, updatedAt: state.updatedAt, objective: objectiveCard(selected, objectiveDecision(state.events, objective)),
    statusCounts, claimableTotal: eligible.size, tasks: { ...window, items: cards },
    inbox: { total: inbox.length, latest: latest.map(eventCard), ...(direct && !latest.includes(direct) ? { latestDirect: eventCard(direct) } : {}) },
    expand: 'Use context --task ID for scope and findings; --page N for more tasks. context --event ID opens one message or finding.'
  };
}
