/** Read-only projection of ledger records. Grouping is not a dependency graph. */
const array = value => Array.isArray(value) ? value : [];
const timestamp = value => typeof value === 'string' ? Date.parse(value) : NaN;
const natural = new Intl.Collator('en', { numeric: true, sensitivity: 'base' });
const compare = (left, right) => {
  const a = String(left ?? ''), b = String(right ?? '');
  return natural.compare(a, b) || (a < b ? -1 : a > b ? 1 : 0);
};
const label = value => value.replace(/[_-]+/g, ' ').replace(/^./, character => character.toUpperCase());
function checkNow(now) {
  if (!Number.isFinite(now)) throw new TypeError('Research tree requires a finite current time.');
}

/** A live lease is a reservation, not evidence that an agent process is running. */
export function taskStatus(task, now = Date.now()) {
  checkNow(now);
  if (task?.status === 'active') {
    return typeof task.assignee === 'string' && task.assignee.length > 0 && timestamp(task.leaseUntil) > now ? 'active' : 'expired';
  }
  return typeof task?.status === 'string' && task.status ? task.status : 'unknown';
}

export function pendingDependencies(task, tasks) {
  return array(task.dependsOn).filter(id => tasks.find(item => item.id === id)?.status !== 'verified');
}

/**
 * Keep raw records by reference, but never sort or modify the input arrays.
 * Branches sort naturally by approach; tasks stay in creation order (title/ID
 * break ties). Events sort newest first, with later ledger entries first on a
 * timestamp tie. Missing dates sort last. Unknown objectives do not fall back.
 *
 * Only explicit objective events and events attached to this objective's tasks
 * belong here. Global activity and contradictory objective/task links do not.
 * All dependency information remains exactly as stored in task.dependsOn.
 */
export function buildResearchTree(state, objectiveId, { now = Date.now() } = {}) {
  checkNow(now);
  const objective = typeof objectiveId === 'string' && objectiveId
    ? array(state?.objectives).find(item => item?.id === objectiveId) ?? null : null;
  if (!objective) return { objective: null, branches: [], tasks: [], agents: [], events: [] };

  const tasks = array(state?.tasks).filter(task => task?.objectiveId === objectiveId).slice().sort((a, b) => {
    const first = timestamp(a.createdAt), second = timestamp(b.createdAt);
    const aTime = Number.isFinite(first) ? first : Infinity, bTime = Number.isFinite(second) ? second : Infinity;
    return (aTime < bTime ? -1 : aTime > bTime ? 1 : 0) || compare(a.title, b.title) || compare(a.id, b.id);
  });
  const taskIds = new Set(tasks.map(task => task.id));
  const branchesById = new Map();
  for (const task of tasks) {
    const id = typeof task.approach === 'string' && task.approach ? task.approach : 'unspecified';
    if (!branchesById.has(id)) branchesById.set(id, { id, label: label(id), tasks: [], activeCount: 0 });
    const branch = branchesById.get(id);
    branch.tasks.push(task);
    if (taskStatus(task, now) === 'active') branch.activeCount++;
  }
  const branches = [...branchesById.values()].sort((a, b) => compare(a.label, b.label) || compare(a.id, b.id));
  const events = array(state?.events).map((event, index) => ({ event, index })).filter(({ event }) => {
    if (!event || (event.objectiveId != null && event.objectiveId !== objectiveId)) return false;
    return event.taskId ? taskIds.has(event.taskId) : event.objectiveId === objectiveId;
  }).sort((a, b) => {
    const first = timestamp(a.event.createdAt), second = timestamp(b.event.createdAt);
    const aTime = Number.isFinite(first) ? first : -Infinity, bTime = Number.isFinite(second) ? second : -Infinity;
    return (aTime > bTime ? -1 : aTime < bTime ? 1 : 0) || b.index - a.index;
  }).map(({ event }) => event);

  const relevantAgents = new Set([objective.createdBy]);
  for (const task of tasks) {
    relevantAgents.add(task.assignee);
    relevantAgents.add(task.createdBy);
    for (const contribution of [...array(task.checkpoints), ...array(task.reviews)]) relevantAgents.add(contribution?.agentId);
  }
  for (const event of events) {
    relevantAgents.add(event.actor);
    if (event.to && event.to !== 'all') relevantAgents.add(event.to);
  }
  const agents = array(state?.agents).filter(agent => agent && relevantAgents.has(agent.id))
    .sort((a, b) => compare(a.name || a.id, b.name || b.id) || compare(a.id, b.id));
  return { objective, branches, tasks, agents, events };
}
