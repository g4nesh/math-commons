const list = (value) => Array.isArray(value) ? value.filter(Boolean) : [];
const label = (value = '') => String(value).replace(/[_-]/g, ' ').replace(/^./, (character) => character.toUpperCase());
const timestamp = (value) => value == null || value === '' ? NaN : new Date(value).getTime();
const short = (value, limit = 320) => String(value || '').length > limit ? `${String(value).slice(0, limit).trimEnd()}…` : String(value || '');

function element(tag, attributes = {}, children = []) {
  const node = document.createElement(tag);
  for (const [key, value] of Object.entries(attributes)) {
    if (value === null || value === undefined || value === false) continue;
    if (key === 'onClick') node.addEventListener('click', value);
    else node.setAttribute(key, String(value));
  }
  for (const child of [children].flat(Infinity)) {
    if (child !== null && child !== undefined && child !== false) node.append(typeof child === 'object' ? child : document.createTextNode(String(child)));
  }
  return node;
}

function safeUrl(value) {
  try {
    const url = new URL(value);
    return url.protocol === 'https:' && !url.username && !url.password ? url : null;
  } catch { return null; }
}

/** A read-only, objective-scoped inspector for nodes in the research tree. */
export function renderTreeInspector({ tree, selection, onSelect = () => {}, state, now = Date.now(), openTask, openObjective, repoUrl } = {}) {
  const objective = tree?.objective;
  const parsedNow = typeof now === 'number' ? now : timestamp(now);
  const currentTime = Number.isFinite(parsedNow) ? parsedNow : Date.now();
  const tasks = objective?.id ? list(tree?.tasks).filter((task) => task.objectiveId === objective.id) : [];
  const taskMap = new Map(tasks.map((task) => [task.id, task]));
  // Explicitly mismatched objective/task references never enter a selected view.
  const events = objective?.id ? list(tree?.events).filter((event) => {
    if (event.objectiveId && event.objectiveId !== objective.id) return false;
    if (event.taskId) return taskMap.has(event.taskId);
    return event.objectiveId === objective.id;
  }) : [];
  const branches = list(tree?.branches).map((branch) => ({ ...branch, tasks: list(branch.tasks).map((task) => taskMap.get(task.id)).filter(Boolean) }));
  const agentIds = new Set([objective?.createdBy]);
  for (const task of tasks) {
    agentIds.add(task.assignee);
    agentIds.add(task.createdBy);
    for (const record of [...list(task.checkpoints), ...list(task.reviews)]) agentIds.add(record.agentId);
  }
  for (const event of events) { agentIds.add(event.actor); if (event.to !== 'all') agentIds.add(event.to); }
  agentIds.delete(undefined); agentIds.delete(null); agentIds.delete('');
  const agents = new Map([...list(state?.agents), ...list(tree?.agents)].filter((agent) => agentIds.has(agent.id)).map((agent) => [agent.id, agent]));
  const agentName = (id) => agents.get(id)?.name || id || 'Unassigned';
  const taskStatus = (task) => ledgerTaskStatus(task, currentTime);
  const live = (task) => taskStatus(task) === 'active';
  const statusText = (status) => ({ active: 'Live lease', review: 'In review', expired: 'Lease expired', unknown: 'Lease unavailable' })[status] || label(status || 'Unknown');
  const badge = (status, text) => element('span', { class: `tree-inspector-status tree-inspector-status-${['active', 'review', 'verified', 'blocked', 'closed', 'expired'].includes(status) ? status : 'neutral'}` }, text || statusText(status));
  const action = (text, handler, className = 'tree-inspector-link') => element('button', { type: 'button', class: className, onClick: handler }, text);
  const select = (type, id, text, className) => action(text, () => onSelect(type, id), className);
  const agentLink = (id) => agents.has(id) ? select('agent', id, agentName(id)) : element('span', {}, agentName(id));
  const note = (text) => element('p', { class: 'tree-inspector-note' }, text);
  const section = (title, children) => element('section', { class: 'tree-inspector-section' }, [element('h3', {}, title), children]);
  const date = (value) => {
    const millis = timestamp(value);
    if (!Number.isFinite(millis)) return element('span', { class: 'tree-inspector-time' }, 'Time not recorded');
    const delta = currentTime - millis;
    let text;
    if (delta < -60_000) text = new Intl.DateTimeFormat(undefined, { month: 'short', day: 'numeric' }).format(millis);
    else if (delta < 60_000) text = 'Just now';
    else if (delta < 3_600_000) text = `${Math.floor(delta / 60_000)}m ago`;
    else if (delta < 86_400_000) text = `${Math.floor(delta / 3_600_000)}h ago`;
    else text = new Intl.DateTimeFormat(undefined, { month: 'short', day: 'numeric' }).format(millis);
    return element('time', { class: 'tree-inspector-time', datetime: new Date(millis).toISOString(), title: new Date(millis).toUTCString() }, text);
  };
  const external = (text, href) => {
    const url = safeUrl(href);
    return url ? element('a', { class: 'tree-inspector-link', href: url.href, target: '_blank', rel: 'noopener noreferrer' }, text) : null;
  };
  const evidence = (values) => {
    const urls = [...new Set(list(values).filter((value) => typeof value === 'string'))].map(safeUrl).filter(Boolean);
    if (!urls.length) return null;
    const links = urls.map((url) => {
      const parts = url.pathname.split('/').filter(Boolean);
      let title = url.hostname;
      if (url.hostname === 'github.com' && parts[2] === 'pull') title = `Pull request #${parts[3] || ''}`;
      else if (url.hostname === 'github.com' && parts[2] === 'commit') title = `Commit ${(parts[3] || '').slice(0, 7)}`;
      else if (parts.length) {
        let file = parts.at(-1);
        try { file = decodeURIComponent(file); } catch { /* Keep the encoded filename. */ }
        title = `${safeUrl(repoUrl)?.origin === url.origin ? 'Artifact' : url.hostname} · ${short(file, 48)}`;
      }
      const link = external(title, url.href);
      link.setAttribute('title', url.href);
      return element('li', {}, link);
    });
    return element('div', { class: 'tree-inspector-evidence' }, [
      element('ul', {}, links.slice(0, 3)),
      links.length > 3 ? element('details', {}, [element('summary', {}, `${links.length - 3} more evidence links`), element('ul', {}, links.slice(3))]) : null
    ]);
  };
  const header = (type, title, status) => element('header', { class: 'tree-inspector-header' }, [
    element('div', { class: 'tree-inspector-kicker' }, [element('span', {}, type), status]),
    element('h2', { id: 'tree-inspector-heading', 'data-inspector-heading': '', tabindex: '-1' }, title)
  ]);
  const panel = element('section', { class: 'tree-inspector', 'aria-labelledby': 'tree-inspector-heading' });
  const append = (...children) => children.filter((child) => child !== null && child !== undefined).forEach((child) => panel.append(child));
  const taskRow = (task) => element('li', { class: 'tree-inspector-task-row' }, [
    select('task', task.id, task.title || task.id), badge(taskStatus(task))
  ]);
  const taskList = (items, limit = 5) => element('div', {}, [
    element('ul', { class: 'tree-inspector-task-list' }, items.slice(0, limit).map(taskRow)),
    items.length > limit ? element('details', { class: 'tree-inspector-more' }, [element('summary', {}, `${items.length - limit} more tasks`), element('ul', { class: 'tree-inspector-task-list' }, items.slice(limit).map(taskRow))]) : null
  ]);
  const taskEvents = (items) => { const ids = new Set(items.map((task) => task.id)); return events.filter((event) => ids.has(event.taskId)); };
  const activity = (items, title = 'Connected activity', limit = 3) => section(title, items.length ? element('ol', { class: 'tree-inspector-activity' }, items.slice(0, limit).map((event) => element('li', {}, [
    element('div', { class: 'tree-inspector-record-meta' }, [element('span', {}, label(event.kind)), date(event.createdAt)]),
    select('event', event.id, short(event.summary || label(event.kind), 200), 'tree-inspector-event-link'),
    element('div', { class: 'tree-inspector-byline' }, [agentLink(event.actor)])
  ]))) : note('No activity recorded for this node yet.'));
  const substantive = (items) => items.filter((event) => !['heartbeat', 'usage', 'claim', 'task', 'objective'].includes(event.kind));
  const latestCheckpoint = (items) => items.flatMap((task) => list(task.checkpoints).map((checkpoint, index) => ({ task, checkpoint, index })))
    .sort((a, b) => (timestamp(b.checkpoint.createdAt) || 0) - (timestamp(a.checkpoint.createdAt) || 0) || b.index - a.index)[0];
  const checkpointCard = (record, showTask = false) => {
    if (!record) return null;
    const { checkpoint, task } = record;
    return section('Latest checkpoint', element('div', { class: 'tree-inspector-checkpoint' }, [
      showTask ? select('task', task.id, task.title || task.id) : null,
      element('div', { class: 'tree-inspector-record-meta' }, [agentLink(checkpoint.agentId), date(checkpoint.createdAt)]),
      element('p', { class: 'tree-inspector-summary' }, short(checkpoint.summary, 260)),
      checkpoint.details || String(checkpoint.summary || '').length > 260 ? element('details', { class: 'tree-inspector-notes' }, [
        element('summary', {}, 'Read checkpoint notes'),
        String(checkpoint.summary || '').length > 260 ? element('p', {}, checkpoint.summary) : null,
        checkpoint.details ? element('p', {}, checkpoint.details) : null,
        element('span', { class: 'tree-inspector-record-id' }, checkpoint.id)
      ]) : null,
      evidence(checkpoint.evidence)
    ]));
  };
  const stats = (items) => element('dl', { class: 'tree-inspector-stats' }, [
    [items.length, 'Tasks'], [items.filter(live).length, 'Live leases'], [items.filter((task) => task.status === 'review').length, 'In review'], [items.filter((task) => task.status === 'verified').length, 'Verified']
  ].map(([value, text]) => element('div', {}, [element('dt', {}, text), element('dd', {}, value)])));
  const fullTask = (task) => typeof openTask === 'function' ? action('Open full task ↗', () => openTask(task.id), 'tree-inspector-full') : null;
  const fullObjective = () => typeof openObjective === 'function' ? action('Read objective & success criteria ↗', () => openObjective(objective.id), 'tree-inspector-full') : null;
  const selected = selection || { type: 'objective', id: objective?.id };

  if (objective && selected.type === 'objective' && selected.id === objective.id) {
    append(header('Objective', 'Research notes', badge(objective.status)), fullObjective());
    if (objective.id === 'marcus-de-oliveira') append(section('The diagram', [
      note('Each black point is an eigenvalue pairing. The red point is the determinant, ⅔ + ²⁄₉ i. This example has dimension three.'),
      external('Read the checked calculation ↗', `${repoUrl}/blob/main/formal/README.md`)
    ]));
    append(stats(tasks));
    if (objective.conclusion) append(section('Recorded conclusion', [note(short(objective.conclusion)), evidence(objective.evidence)]));
    const checkpoint = latestCheckpoint(tasks);
    if (checkpoint) append(checkpointCard(checkpoint, true));
    else append(note('No findings published yet.'));
    const liveTasks = tasks.filter(live);
    if (liveTasks.length) append(section('Reserved now', taskList(liveTasks, 3)));
    append(activity(substantive(events).filter((event) => event.id !== checkpoint?.checkpoint.id), 'Latest activity', 2));
  } else if (selected.type === 'branch') {
    const branch = branches.find((item) => item.id === selected.id);
    if (branch) {
      append(header('Approach', branch.label || label(branch.id)), stats(branch.tasks));
      append(section('Tasks in this direction', branch.tasks.length ? taskList(branch.tasks) : note('No tasks have been added to this direction.')));
      const checkpoint = latestCheckpoint(branch.tasks);
      if (checkpoint) append(checkpointCard(checkpoint, true));
      append(activity(substantive(taskEvents(branch.tasks)).filter((event) => event.id !== checkpoint?.checkpoint.id)));
    }
  } else if (selected.type === 'task') {
    const task = taskMap.get(selected.id);
    if (task) {
      append(header('Task', task.title || task.id, badge(taskStatus(task))));
      const branch = branches.find((item) => item.tasks.some((itemTask) => itemTask.id === task.id));
      if (branch) append(select('branch', branch.id, `↳ ${branch.label || label(branch.id)}`, 'tree-inspector-context'));
      append(fullTask(task));
      if (task.description) append(note(short(task.description)));
      const leaseTime = timestamp(task.leaseUntil);
      let leaseNote = 'No live lease on this task.';
      if (live(task)) leaseNote = `Reserved until ${new Date(leaseTime).toUTCString()}. Process status unknown.`;
      else if (task.status === 'active') leaseNote = Number.isFinite(leaseTime) && leaseTime <= currentTime ? `Lease expired ${new Date(leaseTime).toUTCString()}.` : 'No live lease: assignment or expiry is unavailable.';
      append(section(task.assignee ? 'Assigned agent' : 'Assignment', [
        task.assignee ? element('div', { class: 'tree-inspector-agent-row' }, [agentLink(task.assignee), agents.get(task.assignee)?.owner ? element('span', {}, `@${agents.get(task.assignee).owner}`) : null]) : note('Unassigned'), note(leaseNote)
      ]));
      const checkpoint = latestCheckpoint([task]);
      append(checkpoint ? checkpointCard(checkpoint) : section('Latest checkpoint', note('No checkpoint published yet.')));
      if (list(task.dependsOn).length) append(section('Depends on', element('ul', { class: 'tree-inspector-task-list' }, list(task.dependsOn).map((id) => taskMap.has(id) ? taskRow(taskMap.get(id)) : element('li', { class: 'tree-inspector-note' }, 'Dependency outside this objective or unavailable.')))));
      const dependents = tasks.filter((item) => list(item.dependsOn).includes(task.id));
      if (dependents.length) append(section('Supports', taskList(dependents, 3)));
      const review = list(task.reviews).filter((item) => item.checkpointId === task.submittedCheckpointId).at(-1);
      if (review) append(section('Review of submitted checkpoint', [element('div', { class: 'tree-inspector-record-meta' }, [agentLink(review.agentId), date(review.createdAt)]), note(`${review.verdict === 'pass' ? 'Passed' : 'Changes requested'}: ${review.summary || ''}`), evidence(review.evidence)]));
      append(activity(events.filter((event) => event.taskId === task.id && event.id !== checkpoint?.checkpoint.id && event.kind !== 'heartbeat')));
    }
  } else if (selected.type === 'agent') {
    const agent = agents.get(selected.id);
    if (agent) {
      const assigned = tasks.filter((task) => task.assignee === agent.id);
      const relatedEvents = events.filter((event) => event.actor === agent.id || event.to === agent.id);
      const contributed = tasks.filter((task) => task.assignee === agent.id || task.createdBy === agent.id || [...list(task.checkpoints), ...list(task.reviews)].some((record) => record.agentId === agent.id) || relatedEvents.some((event) => event.taskId === task.id));
      const liveTasks = assigned.filter(live);
      append(header('Agent', agent.name || agent.id, badge(liveTasks.length ? 'active' : 'neutral', liveTasks.length ? 'Live lease' : 'No live lease here')));
      append(element('div', { class: 'tree-inspector-agent-row' }, [agent.owner ? external(`@${agent.owner}`, `https://github.com/${encodeURIComponent(agent.owner)}`) : null, element('span', {}, label(agent.role || 'Agent'))]));
      if (list(agent.capabilities).length) append(element('p', { class: 'tree-inspector-capabilities' }, list(agent.capabilities).join(' · ')));
      if (liveTasks.length) append(note('Lease reserved; process status unknown.'));
      if (assigned.length) append(section('Assigned in this objective', taskList(assigned)));
      const otherTasks = contributed.filter((task) => task.assignee !== agent.id);
      if (otherTasks.length) append(section('Contributed or connected tasks', taskList(otherTasks, 3)));
      const authored = contributed.map((task) => ({ ...task, checkpoints: list(task.checkpoints).filter((checkpoint) => checkpoint.agentId === agent.id) }));
      const checkpoint = latestCheckpoint(authored);
      if (checkpoint) append(checkpointCard(checkpoint, true));
      append(activity(relatedEvents.filter((event) => event.kind !== 'heartbeat' && event.id !== checkpoint?.checkpoint.id), 'Activity in this objective'));
    }
  } else if (selected.type === 'event') {
    const event = events.find((item) => item.id === selected.id);
    if (event) {
      append(header('Activity', label(event.kind)), element('div', { class: 'tree-inspector-record-meta' }, [agentLink(event.actor), date(event.createdAt)]));
      if (event.to) append(element('div', { class: 'tree-inspector-byline' }, ['To ', event.to === 'all' ? 'everyone' : agentLink(event.to)]));
      append(element('p', { class: 'tree-inspector-summary' }, event.summary || 'No summary recorded.'));
      if (event.detail && event.detail !== event.summary) append(section('Published notes', element('p', { class: 'tree-inspector-event-detail' }, event.detail)));
      const links = evidence(event.links);
      if (links) append(section('Evidence', links));
      const task = taskMap.get(event.taskId);
      if (task) append(section('Connected task', taskList([task])), fullTask(task));
      const connected = events.filter((item) => item.id !== event.id && (task ? item.taskId === task.id : !item.taskId));
      if (connected.length) append(activity(connected.filter((item) => item.kind !== 'heartbeat'), 'Related activity', 2));
      append(element('p', { class: 'tree-inspector-record-id' }, event.id));
    }
  }
  if (!panel.childNodes.length) {
    append(header('Research tree', 'Select a node'), note('Choose a task, approach, agent, or activity in this objective to inspect its record.'));
    if (objective) append(select('objective', objective.id, 'Return to this objective', 'tree-inspector-full'));
  }
  return panel;
}
import { taskStatus as ledgerTaskStatus } from './tree-model.js';

