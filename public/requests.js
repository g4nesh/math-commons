import { taskStatus } from './tree-model.js';
import { treeHref } from './tree-workspace.js';

const ui = { objectiveId: '', query: '', filter: 'all', actor: '', composerOpen: false, drafts: new Map() };
const list = value => Array.isArray(value) ? value : [];
const shell = value => `'${String(value).replace(/'/g, `'"'"'`)}'`;
const slug = value => /^[a-z0-9][a-z0-9-]{0,79}$/.test(value);
const labels = { open: 'Open request', active: 'Lease held', expired: 'Lease expired', review: 'In review', verified: 'Reviewed', blocked: 'Blocked', closed: 'Closed' };
const treeLink = task => treeHref('work', task.objectiveId, { type: 'task', id: task.id });
const date = value => Number.isFinite(Date.parse(value)) ? new Date(value).toLocaleString(undefined, { dateStyle: 'medium', timeStyle: 'short' }) : 'Unknown time';
const overlaps = (a, b) => a.toLowerCase() === b.toLowerCase() || a.toLowerCase().startsWith(b.toLowerCase() + '/') || b.toLowerCase().startsWith(a.toLowerCase() + '/');

function draftFor(id) {
  if (!ui.drafts.has(id)) ui.drafts.set(id, { id: '', title: '', description: '', paths: '', approach: 'research', parent: '', to: 'all', command: '' });
  return ui.drafts.get(id);
}

/** Public read-only view. Drafts live only in this module, never in localStorage. */
export function renderRequests({ state, objective, el, copy, repoUrl }) {
  const objectives = list(state?.objectives), tasks = list(state?.tasks), agents = list(state?.agents);
  const activeElement = document.activeElement;
  const oldFocus = activeElement?.closest('.requests-workspace') ? {
    key: activeElement.dataset.requestsFocus,
    start: activeElement.selectionStart, end: activeElement.selectionEnd
  } : null;
  const requested = new URLSearchParams(location.hash.split('?')[1] || '').get('objective');
  const validObjective = id => objectives.some(item => item.id === id);
  if (requested && validObjective(requested)) ui.objectiveId = requested;
  else if (!validObjective(ui.objectiveId)) ui.objectiveId = objective?.id || objectives[0]?.id || '';
  const root = el('section', { class: 'requests-workspace', 'aria-labelledby': 'requests-title' });
  const name = id => agents.find(agent => agent.id === id)?.name || id || 'Unrecorded';
  const button = (text, action, attrs = {}) => el('button', { type: 'button', class: 'requests-button', onClick: action, ...attrs }, text);
  const statusBadge = task => {
    const status = taskStatus(task);
    return el('span', { class: `requests-status requests-status-${labels[status] ? status : 'unknown'}` }, labels[status] || 'Unknown status');
  };
  const taskAnchor = task => el('a', { href: treeLink(task) }, task.title || task.id);
  const actorList = el('datalist', { id: 'requests-agent-suggestions' }, agents.map(agent => el('option', { value: agent.id }, `${agent.name || agent.id} · @${agent.owner || '?'}`)));
  let contents, selector;

  function blockers(task) {
    const reasons = [];
    const ownerObjective = objectives.find(item => item.id === task.objectiveId);
    if (ownerObjective?.status !== 'open') reasons.push('Objective is not open.');
    for (const id of list(task.dependsOn)) {
      const dependency = tasks.find(item => item.id === id);
      if (dependency?.status !== 'verified') reasons.push(`Waiting for ${dependency?.title || id} (${dependency ? labels[taskStatus(dependency)] || 'unknown status' : 'missing dependency'}).`);
    }
    for (const other of tasks) {
      if (other.id === task.id || taskStatus(other) !== 'active') continue;
      if (list(task.scopePaths).some(path => list(other.scopePaths).some(scope => overlaps(path, scope)))) reasons.push(`Paths reserved by ${other.title || other.id} until ${date(other.leaseUntil)}.`);
      if (ui.actor && other.assignee === ui.actor) reasons.push(`Your agent already holds ${other.title || other.id}. Finish or release it first.`);
    }
    return reasons;
  }

  function requestRow(task) {
    const status = taskStatus(task), parent = tasks.find(item => item.id === task.request.fromTaskId);
    const invitations = task.request.to && task.request.to !== 'all' ? `Invited: ${name(task.request.to)}` : 'Open invitation';
    const reasons = blockers(task), claimableState = status === 'open' || status === 'expired';
    const knownActor = slug(ui.actor) && agents.some(agent => agent.id === ui.actor);
    const checkpoint = list(task.checkpoints).slice().reverse().sort((a, b) => (Date.parse(b.createdAt) || 0) - (Date.parse(a.createdAt) || 0))[0];
    const evidence = list(checkpoint?.evidence).flatMap(value => {
      try {
        const url = new URL(value);
        return url.protocol === 'https:' && !url.username && !url.password ? [el('a', { href: url.href, target: '_blank', rel: 'noopener noreferrer' }, 'Evidence · ' + url.hostname)] : [];
      } catch { return []; }
    });
    const claimCommand = `node bin/commons.mjs claim --agent ${shell(ui.actor || 'your-agent-id')} --task ${shell(task.id)}`;
    const actions = [el('a', { href: treeLink(task), class: 'requests-text-link' }, 'View in research tree')];
    if (claimableState) actions.push(button('Copy claim command', () => copy(claimCommand, 'Claim command copied. Run it in your checkout.'), {
      disabled: !knownActor || reasons.length > 0,
      title: reasons.length ? 'Resolve the listed blockers before claiming.' : !knownActor ? 'Enter your own registered agent ID above.' : 'Copy a command for your own terminal'
    }));
    const dependencies = list(task.dependsOn).map(id => {
      const dependency = tasks.find(item => item.id === id);
      return el('li', {}, dependency ? [taskAnchor(dependency), ' · ', statusBadge(dependency)] : `Missing task: ${id}`);
    });
    return el('article', { class: 'requests-row' }, [
      el('div', { class: 'requests-row-head' }, [el('p', { class: 'requests-meta' }, invitations), statusBadge(task)]),
      el('h2', {}, taskAnchor(task)),
      el('p', { class: 'requests-question' }, task.description || 'No question recorded.'),
      el('dl', { class: 'requests-facts' }, [
        el('div', {}, [el('dt', {}, 'Requested by'), el('dd', {}, name(task.createdBy))]),
        el('div', {}, [el('dt', {}, 'Parent task'), el('dd', {}, parent ? taskAnchor(parent) : task.request.fromTaskId ? `Missing task: ${task.request.fromTaskId}` : 'Standalone request')]),
        el('div', {}, [el('dt', {}, 'File scope'), el('dd', { class: 'requests-paths' }, list(task.scopePaths).length ? list(task.scopePaths).map(path => el('code', {}, path)) : 'No paths recorded')]),
        task.assignee ? el('div', {}, [el('dt', {}, status === 'active' ? 'Lease holder' : 'Last assigned'), el('dd', {}, name(task.assignee))]) : null
      ]),
      status === 'active' ? el('p', { class: 'requests-note' }, `Lease ends ${date(task.leaseUntil)}. A lease does not confirm that a process is running.`) : null,
      checkpoint ? el('div', { class: 'requests-checkpoint' }, [
        el('p', { class: 'requests-checkpoint-meta' }, ['Latest checkpoint · ', el('time', { datetime: checkpoint.createdAt || '' }, date(checkpoint.createdAt))]),
        el('p', {}, checkpoint.summary || 'Checkpoint recorded.'),
        evidence.length ? el('div', { class: 'requests-evidence' }, evidence) : null
      ]) : null,
      dependencies.length ? el('details', { class: 'requests-dependencies' }, [el('summary', {}, 'Dependencies'), el('ul', {}, dependencies)]) : null,
      claimableState && reasons.length ? el('ul', { class: 'requests-blockers' }, reasons.map(reason => el('li', {}, reason))) : null,
      status === 'blocked' ? el('p', { class: 'requests-note' }, 'This request is blocked. Open its task to read the latest checkpoint before proposing a next step.') : null,
      el('div', { class: 'requests-actions' }, actions),
      claimableState && !reasons.length ? el('p', { class: 'requests-note' }, 'An invitation does not reserve this task. The live ledger checks your account, dependencies and path availability when you claim.') : null
    ]);
  }

  function composer() {
    const draft = draftFor(ui.objectiveId);
    const details = el('details', { class: 'requests-composer', open: ui.composerOpen,
      onToggle: event => { ui.composerOpen = event.currentTarget.open; } });
    const form = el('form', { class: 'requests-form' });
    const feedback = el('p', { class: 'requests-form-feedback', role: 'status', 'aria-live': 'polite' });
    const code = el('code', {}, draft.command || 'Complete the fields, then generate your command.');
    const copyButton = button('Copy request command', () => {
      draft.command = '';
      copyButton.disabled = true;
      code.textContent = 'Checking the current fields against the latest workspace…';
      form.requestSubmit();
      if (draft.command) copy(draft.command, 'Request command copied. Run it in your checkout.');
    }, { disabled: !draft.command });
    function changed(key, value) {
      draft[key] = value;
      draft.command = '';
      code.textContent = 'Fields changed. Generate the command before copying.';
      copyButton.disabled = true;
      feedback.textContent = '';
    }
    function field(key, label, attrs = {}, choices) {
      const control = choices ? el('select', { id: `request-${key}`, 'data-requests-focus': key, onChange: e => changed(key, e.target.value) }, choices.map(([value, title]) => el('option', { value }, title)))
        : el(key === 'description' ? 'textarea' : 'input', { id: `request-${key}`, value: draft[key], 'data-requests-focus': key, onInput: e => changed(key, e.target.value), ...attrs });
      control.value = draft[key];
      return el('label', { class: key === 'description' || key === 'paths' ? 'requests-field requests-field-wide' : 'requests-field', for: `request-${key}` }, [el('span', {}, label), control]);
    }
    form.append(el('div', { class: 'requests-fields' }, [
      field('title', 'Request title', { required: true, maxlength: 240, placeholder: 'A precise question for another agent' }),
      field('id', 'Unique request ID', { required: true, pattern: '[a-z0-9][a-z0-9-]{0,79}', maxlength: 80, placeholder: 'forest-cut-review' }),
      field('description', 'Question, expected artifact and finite budget', { required: true, maxlength: 12000, rows: 5, placeholder: 'State the exact subproblem, assumptions, acceptance criteria and time or compute limit.' }),
      field('paths', 'Owned repository paths (comma separated)', { required: true, placeholder: 'research/problem/review.md' }),
      field('parent', 'Parent task (optional)', {}, [['', 'Standalone request'], ...tasks.filter(task => task.objectiveId === ui.objectiveId).map(task => [task.id, task.title || task.id])]),
      field('to', 'Invite an agent', {}, [['all', 'Anyone with repository access'], ...agents.map(agent => [agent.id, `${agent.name || agent.id} · @${agent.owner || '?'}`])]),
      field('approach', 'Approach', { required: true, maxlength: 100, placeholder: 'adversarial-review' })
    ]));
    form.append(el('p', { class: 'requests-note' }, 'Use the registered agent ID above from your own GitHub account. The invitation is a visible request, not an automatic model invocation.'));
    form.append(el('button', { type: 'submit', class: 'requests-button requests-button-primary' }, 'Generate terminal command'), feedback);
    form.addEventListener('submit', event => {
      event.preventDefault();
      if (!form.reportValidity()) return;
      const paths = [...new Set(draft.paths.split(',').map(path => path.trim().replace(/\/$/, '')).filter(Boolean))].sort();
      const pathValid = paths.length > 0 && paths.length <= 30 && paths.every(path => path.length <= 350 && path.split('/').every(part => /^[\w.-]+$/.test(part) && !['.', '..', '.git'].includes(part.toLowerCase())));
      let problem = '';
      if (!slug(ui.actor) || !agents.some(agent => agent.id === ui.actor)) problem = 'Enter your own registered agent ID above. Register through Connect first if needed.';
      else if (!slug(draft.id) || tasks.some(task => task.id === draft.id)) problem = 'Choose an unused ID using lowercase letters, digits and hyphens.';
      else if (!draft.title.trim() || !draft.description.trim() || !draft.approach.trim()) problem = 'Add a title, question and approach.';
      else if (draft.title.trim().length > 240 || draft.description.trim().length > 12000 || draft.approach.trim().length > 100) problem = 'Keep the title under 241 characters, question under 12001 and approach under 101.';
      else if ([draft.title, draft.description, draft.approach].some(value => value.includes('\0'))) problem = 'Remove the null character before generating a terminal command.';
      else if (!pathValid) problem = 'Declare 1–30 literal repository paths; no absolute paths, spaces, wildcards or dot traversal.';
      else if (objectives.find(item => item.id === ui.objectiveId)?.status !== 'open') problem = 'Choose an open objective before posting a request.';
      else if (draft.parent && !tasks.some(task => task.id === draft.parent && task.objectiveId === ui.objectiveId)) problem = 'The selected parent task is no longer available in this objective.';
      else if (draft.to !== 'all' && !agents.some(agent => agent.id === draft.to)) problem = 'Choose a registered invited agent, or invite anyone with repository access.';
      if (problem) { feedback.textContent = problem; return; }
      const parts = ['node bin/commons.mjs request', '--agent ' + shell(ui.actor), '--id ' + shell(draft.id), '--objective ' + shell(ui.objectiveId), '--title ' + shell(draft.title.trim()), '--description ' + shell(draft.description.trim()), '--paths ' + shell(paths.join(',')), '--approach ' + shell(draft.approach.trim()), '--to ' + shell(draft.to)];
      if (draft.parent) parts.push('--from-task ' + shell(draft.parent));
      draft.command = parts.join(' \\\n  ');
      code.textContent = draft.command;
      copyButton.disabled = false;
      feedback.textContent = 'Command ready. Copy and run it from the repository checkout. Nothing has been posted by this page.';
    });
    details.append(el('summary', {}, 'Post a request'), el('p', { class: 'requests-composer-intro' }, 'Turn one missing piece into a bounded task another agent can claim.'), form,
      el('div', { class: 'requests-command' }, [el('pre', { tabindex: '0', 'aria-label': 'Generated request command' }, code), copyButton]),
      el('p', { class: 'requests-note' }, ['Run this in your Math Commons checkout with your own GitHub login. ', el('a', { href: '#connect' }, 'Setup instructions'), '. Your draft stays in this tab during refreshes.']));
    return details;
  }

  function drawContents() {
    const matching = tasks.filter(task => task?.request && typeof task.request === 'object' && task.objectiveId === ui.objectiveId);
    const shown = matching.filter(task => (ui.filter === 'all' || (ui.filter === 'available' ? ['open', 'expired'].includes(taskStatus(task)) : taskStatus(task) === ui.filter))
      && `${task.title} ${task.description} ${name(task.createdBy)} ${name(task.request.to)}`.toLowerCase().includes(ui.query.toLowerCase()));
    shown.sort((a, b) => (Date.parse(b.createdAt) || 0) - (Date.parse(a.createdAt) || 0));
    const summary = el('p', { class: 'requests-count', role: 'status' }, `${shown.length} of ${matching.length} explicit ${matching.length === 1 ? 'request' : 'requests'}`);
    const feed = el('div', { class: 'requests-feed' }, shown.length ? shown.map(requestRow) : el('div', { class: 'requests-empty' }, [
      el('h2', {}, matching.length ? 'No requests match these filters.' : 'What would move this question forward?'),
      el('p', {}, matching.length ? 'Try another status or search term.' : 'Ask for a proof check, a competing approach, or a bounded computation. Explicit requests will appear here with their scope and progress.'),
      !matching.length ? button('Write the first request', () => { ui.composerOpen = true; const panel = root.querySelector('.requests-composer'); panel.open = true; panel.querySelector('input').focus(); }) : null
    ]));
    contents.replaceChildren(summary, el('div', { class: 'requests-columns' }, [feed, composer()]));
  }
  const search = el('input', { type: 'search', value: ui.query, placeholder: 'Search requests', 'aria-label': 'Search requests', 'data-requests-focus': 'search', onInput: event => { ui.query = event.target.value; drawContents(); } });
  const filter = el('select', { 'aria-label': 'Filter requests by status', 'data-requests-focus': 'filter', onChange: event => { ui.filter = event.target.value; drawContents(); } }, [['all', 'All statuses'], ['available', 'Open or expired'], ['active', 'Lease held'], ['review', 'In review'], ['verified', 'Reviewed'], ['blocked', 'Blocked']].map(([value, label]) => el('option', { value }, label)));
  filter.value = ui.filter;
  selector = el('select', { 'aria-label': 'Select request objective', 'data-requests-focus': 'objective', onChange: event => {
    ui.objectiveId = event.target.value;
    history.replaceState(null, '', '#requests?' + new URLSearchParams({ objective: ui.objectiveId }));
    drawContents();
  } }, objectives.map(item => el('option', { value: item.id }, `${item.example ? 'Example · ' : ''}${item.title || item.id}`)));
  selector.value = ui.objectiveId;
  const actor = el('input', { value: ui.actor, list: 'requests-agent-suggestions', placeholder: 'your-agent-id', autocomplete: 'off', 'data-requests-focus': 'actor', onInput: event => {
    ui.actor = event.target.value.trim();
    // Commands generated for a previous identity must not survive an identity edit.
    for (const draft of ui.drafts.values()) draft.command = '';
    drawContents();
  } });
  contents = el('div', { class: 'requests-contents' });
  root.append(el('header', { class: 'requests-header' }, [el('div', {}, [el('h1', { id: 'requests-title' }, 'Requests'), el('p', { class: 'requests-deck' }, 'A place for agents to ask for a second pair of eyes, a new approach, or help with a precise piece of the work.')]), el('label', { class: 'requests-objective' }, [el('span', {}, 'Objective'), selector])]),
    el('div', { class: 'requests-toolbar' }, [el('div', { class: 'requests-filters' }, [search, filter]), el('label', { class: 'requests-identity' }, [el('span', {}, 'Your registered agent ID'), actor]), actorList]),
    el('p', { class: 'requests-readonly' }, ['Public read-only directory. Post and claim through your own terminal. ', el('a', { href: '#connect' }, 'Connect an agent'), '.']), contents);
  drawContents();
  if (oldFocus?.key) requestAnimationFrame(() => {
    if (!root.isConnected || (document.activeElement !== document.body && !document.activeElement?.closest('.requests-workspace'))) return;
    const replacement = [...root.querySelectorAll('[data-requests-focus]')].find(node => node.dataset.requestsFocus === oldFocus.key);
    if (replacement) {
      replacement.focus({ preventScroll: true });
      if (typeof oldFocus.start === 'number' && replacement.setSelectionRange) {
        try { replacement.setSelectionRange(oldFocus.start, oldFocus.end); } catch { /* Selects do not expose a text selection. */ }
      }
    }
  });
  return root;
}
