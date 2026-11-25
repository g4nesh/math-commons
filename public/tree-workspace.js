import { buildResearchTree, taskStatus } from './tree-model.js';
import { renderResearchMap } from './research-map.js';
import { renderTreeInspector } from './tree-inspector.js';

export const isTreeView = view => ['work', 'activity', 'agents'].includes(view);
const nodeSets = { task: 'tasks', agent: 'agents', branch: 'branches', event: 'events' };

export function treeHref(view, objectiveId, selection) {
  const query = new URLSearchParams({ objective: objectiveId });
  if (selection && selection.type !== 'objective') query.set(selection.type, selection.id);
  return '#' + view + '?' + query.toString();
}

function readSelection(tree, query) {
  for (const [type, collection] of Object.entries(nodeSets)) {
    const id = query.get(type);
    if (id && tree[collection].some(record => record.id === id)) return { type, id };
  }
  return { type: 'objective', id: tree.objective.id };
}

export function renderTreeWorkspace({ state, view, objective, query, el, time, openTask,
  openObjective, newObjective, repoUrl, navigate, copyLink, now = Date.now() }) {
  const tree = buildResearchTree(state, objective.id, { now });
  const selection = readSelection(tree, query);
  const root = tree.objective;
  const onSelect = (type, id) => navigate(treeHref('work', root.id, { type, id }),
    type === 'branch' || type === 'objective' ? 'path' : 'node');
  const count = (n, singular, plural = singular + 's') => n + ' ' + (n === 1 ? singular : plural);
  const plainButton = (label, click, className) => el('button', { type: 'button', class: className, onClick: click }, label);
  const indexLink = (label, href) => el('a', { class: 'tree-index-link', href }, label);
  const panelHeader = el('nav', { class: 'tree-panel-nav', 'aria-label': 'Explore this objective' },
    [['work', 'Details'], ['activity', 'Activity'], ['agents', 'Agents']].map(([target, label]) => {
      const href = treeHref(target, root.id, selection);
      return el('a', { href, 'aria-current': target === view ? 'page' : null,
        onClick: e => {
          if (e.metaKey || e.ctrlKey || e.shiftKey || e.altKey) return;
          e.preventDefault(); navigate(href, 'tab');
        } }, label);
    }));
  const rows = (items, type, row) => items.map(item => plainButton(row(item), () => onSelect(type, item.id), 'tree-index-row'));
  const listHeading = title => el('h2', { class: 'tree-list-heading', 'data-inspector-heading': '', tabindex: '-1' }, title);
  const selectedRecord = tree[nodeSets[selection.type]]?.find(record => record.id === selection.id);
  const selectedTasks = selection.type === 'branch' ? selectedRecord.tasks :
    selection.type === 'task' ? [selectedRecord] : selection.type === 'event' && selectedRecord.taskId ?
      tree.tasks.filter(task => task.id === selectedRecord.taskId) : null;
  const taskIds = selectedTasks && new Set(selectedTasks.map(task => task.id));
  const visibleEvents = taskIds ? tree.events.filter(event => taskIds.has(event.taskId)) :
    selection.type === 'agent' ? tree.events.filter(event => event.actor === selection.id || event.to === selection.id) : tree.events;
  function scopeNote() {
    if (!selectedRecord) return null;
    return el('div', {class:'tree-scope-note'}, [
      el('span', {}, selectedRecord.title || selectedRecord.label || selectedRecord.name || selectedRecord.summary),
      plainButton('Show whole objective', () => navigate(treeHref(view, root.id), 'tab'), 'tree-index-link')]);
  }
  function activity() {
    return el('div', { class: 'tree-index' }, [listHeading('Activity'), scopeNote(),
      visibleEvents.length ? el('div', { class: 'tree-event-list' }, rows(visibleEvents, 'event', event => [
        el('span', { class: 'tree-row-meta' }, [event.kind, time(event.createdAt)]),
        el('span', { class: 'tree-row-title' }, event.summary)
      ])) : el('p', { class: 'tree-empty' }, 'No activity recorded here yet.'),
      indexLink('All workspace activity', '#log')]);
  }
  function agents() {
    const relevantIds = selectedTasks && new Set([
      ...selectedTasks.flatMap(t => [t.assignee, t.createdBy, ...(t.checkpoints || []).map(c => c.agentId), ...(t.reviews || []).map(r => r.agentId)]),
      ...visibleEvents.flatMap(e => [e.actor,e.to])]);
    const participants = tree.agents.filter(a => (!relevantIds || relevantIds.has(a.id)) && (selection.type !== 'agent' || a.id === selection.id));
    const withLease = new Set((selectedTasks || tree.tasks).filter(t => taskStatus(t,now) === 'active').map(t => t.assignee));
    participants.sort((a, b) => Number(withLease.has(b.id)) - Number(withLease.has(a.id)));
    const liveCount = participants.filter(a => withLease.has(a.id)).length;
    return el('div', { class: 'tree-index' }, [listHeading('Agents'), scopeNote(), el('p', { class: 'tree-index-subtitle' },
      liveCount ? count(liveCount, 'agent') + ' with a live task lease' : 'No active task leases'),
      rows(participants, 'agent', agent => {
        const task = (selectedTasks || tree.tasks).find(t => t.assignee === agent.id && taskStatus(t, now) === 'active');
        return [el('span', { class: 'tree-row-title' }, agent.name || agent.id),
          el('span', { class: 'tree-person-context' }, task ? 'Active · ' + task.title : '@' + agent.owner + ' · ' + (agent.role || 'Participant'))];
      }), participants.length ? null : el('p', {class:'tree-empty'}, 'No participants recorded here yet.'), indexLink('Full agent directory', '#directory')]);
  }
  const selector = el('select', { class: 'tree-objective-select', 'aria-label': 'Switch objective',
    onChange: e => navigate(treeHref(view, e.target.value), 'objective') },
    state.objectives.map(o => el('option', { value: o.id }, (o.example ? 'Example · ' : '') + o.title)));
  selector.value = root.id;
  return el('section', { class: 'tree-workspace', 'aria-label': 'Research tree' }, [
    el('header', { class: 'tree-page-heading' }, [
      el('div', { class: 'tree-heading-main' }, [
        el('p', { class: 'tree-kicker' }, root.example ? 'Example research tree' : 'Research tree'),
        el('h1', {}, root.title),
        el('p', { class: 'tree-counts' }, [count(tree.branches.length, 'path'), el('span', {}, ' / '),
          count(tree.tasks.length, 'task'), el('span', {}, ' / '), el('a', {
            class: 'tree-agents-count', href: treeHref('agents', root.id),
            onClick: event => {
              if (event.metaKey || event.ctrlKey || event.shiftKey || event.altKey) return;
              event.preventDefault(); navigate(treeHref('agents', root.id), 'tab');
            }
          }, count(tree.agents.length, 'agent'))])]),
      el('div', { class: 'tree-switcher' }, [
        el('label', { class: 'tree-switch-button' }, [el('span', { 'aria-hidden': 'true' }, 'Switch objective'), selector]),
        plainButton('New objective', newObjective, 'tree-new-objective')])]),
    el('div', { class: 'tree-layout' }, [
      el('div', { id: 'research-canvas', class: 'tree-canvas-wrap', tabindex: '-1' },
        renderResearchMap({ tree, selection, onSelect, now })),
      el('aside', { id: 'research-panel', class: 'tree-panel', 'aria-label': 'Research details' }, [
        panelHeader, el('div', { class: 'tree-panel-scroll' }, view === 'activity' ? activity() : view === 'agents' ? agents() :
          renderTreeInspector({ tree, selection, onSelect, state, now, openTask, openObjective, repoUrl })),
        plainButton('Back to tree', () => { const canvas = document.getElementById('research-canvas'); canvas.focus({ preventScroll:true }); canvas.scrollIntoView({block:'start'}); }, 'tree-back')])]),
    el('div', { class: 'tree-links' }, [
      plainButton('Copy selection link', () => copyLink(new URL(treeHref(view, root.id, selection), location.href).href), 'tree-index-link'),
      indexLink('Search task index', '#tasks'), indexLink('Request contributor access', '/join.html')])
  ]);
}
