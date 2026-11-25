import { taskStatus } from './tree-model.js';

const SVG = 'http://www.w3.org/2000/svg';
const viewports = new Map();
let activeInstance = null;
let pendingFocus = null;
let sequence = 0;

function html(tag, attributes = {}, children = []) {
  const node = document.createElement(tag);
  for (const [key, value] of Object.entries(attributes)) {
    if (value === undefined || value === null) continue;
    if (key === 'class') node.className = value;
    else if (key.startsWith('on') && typeof value === 'function') node.addEventListener(key.slice(2).toLowerCase(), value);
    else node.setAttribute(key, String(value));
  }
  append(node, children);
  return node;
}
function svg(tag, attributes = {}, children = []) {
  const node = document.createElementNS(SVG, tag);
  for (const [key, value] of Object.entries(attributes)) node.setAttribute(key, String(value));
  append(node, children);
  return node;
}
function append(node, children) {
  for (const child of [children].flat(Infinity)) {
    if (child !== undefined && child !== null && child !== false) node.append(child instanceof Node ? child : document.createTextNode(String(child)));
  }
}
const keyOf = (type, id) => `${type}:${id}`;
const labels = { open: 'Open', active: 'Lease held', expired: 'Lease expired', review: 'In review', verified: 'Reviewed', blocked: 'Blocked', closed: 'Closed', unknown: 'Unknown' };
const statusOf = (task, now) => Object.hasOwn(labels, taskStatus(task, now)) ? taskStatus(task, now) : 'unknown';
const cubic = (x1, y1, x2, y2) => `M${x1},${y1} C${(x1 + x2) / 2},${y1} ${(x1 + x2) / 2},${y2} ${x2},${y2}`;

// Coordinates are the six permutation products in formal/MathCommons.lean.
function witnessDiagram() {
  const points = [[6,-2], [2,0], [2,1], [0,-1], [0,2], [2,-2]];
  const scale = ([x,y]) => [11 + 14*x, 39 - 14*y];
  const at = point => scale(point).join(',');
  const [cx,cy] = scale([2/3,2/9]);
  return svg('svg', { class:'rmap-witness', viewBox:'0 0 110 79', 'aria-hidden':'true' }, [
    svg('path', { d:'M4,39 H104 M11,5 V73', class:'rmap-axis' }),
    svg('polygon', { points:[3,5,0,2,4].map(i => at(points[i])).join(' '), class:'rmap-polygon' }),
    ...points.map(point => { const [x,y]=scale(point); return svg('circle',{cx:x,cy:y,r:2,class:'rmap-vertex'}); }),
    svg('circle', { cx,cy,r:3.8,class:'rmap-determinant' })
  ]);
}

function relatedSelection(tree, selection, now) {
  const event = selection?.type === 'event' ? tree.events.find((item) => item.id === selection.id) : null;
  const selectedTask = selection?.type === 'task' ? tree.tasks.find((item) => item.id === selection.id) : event?.taskId ? tree.tasks.find((item) => item.id === event.taskId) : null;
  const agentTasks = selection?.type === 'agent' ? tree.tasks.filter((task) => task.assignee === selection.id
    || task.createdBy === selection.id
    || (task.checkpoints || []).some((item) => item.agentId === selection.id)
    || (task.reviews || []).some((item) => item.agentId === selection.id)
    || tree.events.some((item) => item.taskId === task.id && (item.actor === selection.id || item.to === selection.id))) : [];
  const primary = selectedTask || agentTasks.find((task) => task.assignee === selection?.id && statusOf(task, now) === 'active') || agentTasks[0];
  const related = new Set(selectedTask ? [selectedTask.id] : agentTasks.map((task) => task.id));
  const focusBranch = selection?.type === 'branch' ? tree.branches.find((branch) => branch.id === selection.id) : primary ? tree.branches.find((branch) => branch.tasks.some((task) => task.id === primary.id)) : null;
  return { selectedTask, related, primary, focusBranch };
}

function overviewLayout(tree, availableWidth, narrow) {
  const width = narrow ? Math.max(260, availableWidth - 2) : 820;
  const half = Math.ceil(tree.branches.length / 2);
  const overviewHeight = Math.max(420, Math.max(0, half - 1) * 77 + 100);
  const root = narrow ? { x: width / 2 - 106, y: 16, width: 212, height: 96 } : { x: 328, y: overviewHeight / 2 - 55, width: 164, height: 110 };
  const branches = [], tasks = [];
  tree.branches.forEach((branch, index) => {
    const side = index < half ? 'left' : 'right';
    const row = side === 'left' ? index : index - half;
    const sideCount = side === 'left' ? half : tree.branches.length - half;
    const branchWidth = narrow ? Math.min(190, (width - 64) / 2) : 175;
    const branchNode = {
      data: branch, color: index % 7, side,
      x: narrow ? side === 'left' ? 20 : width - 20 - branchWidth : side === 'left' ? 64 : 581,
      y: narrow ? 132 + row * 77 : 35 + (sideCount === 1 ? (overviewHeight - 112) / 2 : row * (overviewHeight - 112) / (sideCount - 1)),
      width: branchWidth, height: narrow ? 38 : 44
    };
    branches.push(branchNode);
    const shown = branch.tasks.slice(0, 3);
    shown.forEach((task, taskIndex) => tasks.push({
      data: task, branchId: branch.id, color: index % 7, dot: true, side,
      x: narrow ? branchNode.x + branchNode.width / 2 - 10 + (taskIndex - (shown.length - 1) / 2) * 25 : side === 'left' ? 28 : 772,
      y: narrow ? branchNode.y + branchNode.height + 11 : branchNode.y + branchNode.height / 2 - 11 + (taskIndex - (shown.length - 1) / 2) * 26,
      width: 22, height: 22
    }));
  });
  return { root, branches, tasks, narrow, overview: true, width, height: narrow ? Math.max(330, 132 + Math.max(half, tree.branches.length - half) * 77 + 4) : overviewHeight };
}

function taskMeta(task, tree, now) {
  const status = statusOf(task, now);
  const agent = status === 'active' ? tree.agents.find((item) => item.id === task.assignee) : null;
  return status === 'active' && task.assignee ? `${labels[status]} · ${agent?.name || task.assignee}` : labels[status];
}

function taskHeight(task, width, tree, now) {
  const context = document.createElement('canvas').getContext('2d');
  if (!context) return Math.max(84, 47 + Math.ceil(String(task.title || task.id).length / Math.max(12, (width - 32) / 9)) * 19);
  context.font = `500 14px ${getComputedStyle(document.documentElement).fontFamily}`;
  const limit = Math.max(80, width - 34);
  let line = '', lines = 1;
  for (const word of String(task.title || task.id).trim().split(/\s+/)) {
    const next = line ? `${line} ${word}` : word;
    if (context.measureText(next).width <= limit) { line = next; continue; }
    if (line) { lines++; line = ''; }
    for (const character of word) {
      if (line && context.measureText(line + character).width > limit) { lines++; line = ''; }
      line += character;
    }
  }
  context.font = `400 10px ${getComputedStyle(document.documentElement).fontFamily}`;
  const metaLines = Math.max(statusOf(task, now) === 'active' ? 2 : 1, Math.ceil(context.measureText(taskMeta(task, tree, now)).width / Math.max(80, width - 50)));
  return Math.max(84, Math.ceil(lines * 18.9) + 48 + (metaLines - 1) * 14);
}

function focusedLayout(branch, tree, availableWidth, narrow, now) {
  const width = narrow ? Math.max(260, availableWidth - 2) : 820;
  const taskWidth = narrow ? width - 62 : 252;
  const heights = branch.tasks.map((task) => taskHeight(task, taskWidth, tree, now));
  const groupHeight = heights.reduce((sum, height) => sum + height, 0) + Math.max(0, heights.length - 1) * 14;
  const height = narrow ? Math.max(360, 261 + groupHeight) : Math.max(350, groupHeight + 70);
  const index = tree.branches.findIndex((item) => item.id === branch.id);
  const root = narrow ? { x: 22, y: 20, width: width - 44, height: 96 } : { x: 26, y: height / 2 - 60, width: 208, height: 120 };
  const approach = { data: branch, color: index % 7, x: narrow ? 38 : 312, y: narrow ? 151 : height / 2 - 23, width: narrow ? width - 76 : 156, height: 46 };
  let top = narrow ? 231 : (height - groupHeight) / 2;
  const tasks = branch.tasks.map((task, i) => {
    const box = { data: task, color: index % 7, branchId: branch.id, x: narrow ? 40 : 538, y: top, width: taskWidth, height: heights[i], dot: false };
    top += heights[i] + 14;
    return box;
  });
  return { root, branches: [approach], tasks, narrow, overview: false, width, height };
}

/** Native node buttons expose data-research-node="type:id". Nodes are drawn
 * in the first animation frame; parent focus restoration can run after two.
 * Viewport is retained per objective and focused approach across remounts. */
export function renderResearchMap({ tree, selection = null, onSelect, now = Date.now() }) {
  activeInstance?.save();
  activeInstance?.disconnect();
  if (!tree?.objective) return html('section', { class: 'research-map rmap-empty' }, [html('h2', {}, 'A question starts the map.'), html('p', {}, 'Choose an objective to explore its research paths.')]);
  tree = { ...tree, branches: tree.branches || [], tasks: tree.tasks || [], agents: tree.agents || [], events: tree.events || [] };
  const context = relatedSelection(tree, selection, now);
  const modeKey = `${tree.objective.id}|${context.focusBranch?.id || 'overview'}`;
  if (!viewports.has(modeKey)) viewports.set(modeKey, { zoom: null, left: 0, top: 0, narrow: null });
  const saved = viewports.get(modeKey);
  const root = html('section', { class: 'research-map', 'aria-label': 'Interactive research map' });
  const stage = html('div', { class: 'rmap-stage', tabindex: '0', role: 'region', 'aria-label': 'Connected research map. Select an approach or task. Tab moves between nodes; scroll or drag to explore.' });
  const canvas = svg('svg', { class: 'rmap-canvas', role: 'group', 'aria-label': context.focusBranch ? `${context.focusBranch.label} research path` : 'Whole research tree' });
  stage.append(canvas);
  const zoomLabel = html('span', { class: 'rmap-zoom', 'aria-label': 'Map zoom', 'aria-live': 'polite' }, '100%');
  const nodeButtons = new Map();
  const marker = `rmap-arrow-${++sequence}`;
  let zoom = saved.zoom || 1, layout, observer, currentWidth = 0, firstDraw = true;
  const select = (type, id) => {
    pendingFocus = { objective: tree.objective.id, key: keyOf(type, id) };
    save();
    onSelect?.(type, id);
  };
  function save() { saved.zoom = zoom; saved.left = stage.scrollLeft; saved.top = stage.scrollTop; }
  function setZoom(next, center = true) {
    if (!layout) return;
    const cx = (stage.scrollLeft + stage.clientWidth / 2) / zoom;
    const cy = (stage.scrollTop + stage.clientHeight / 2) / zoom;
    zoom = Math.max(.5, Math.min(1.75, next));
    canvas.setAttribute('width', Math.ceil(layout.width * zoom));
    canvas.setAttribute('height', Math.ceil(layout.height * zoom));
    zoomLabel.textContent = `${Math.round(zoom * 100)}%`;
    if (center) { stage.scrollLeft = cx * zoom - stage.clientWidth / 2; stage.scrollTop = cy * zoom - stage.clientHeight / 2; }
    save();
  }
  function fit() {
    if (!layout) return;
    setZoom(Math.min(1, (stage.clientWidth - 2) / layout.width, layout.overview ? (stage.clientHeight - 2) / layout.height : 1), false);
    stage.scrollLeft = 0; stage.scrollTop = 0; save();
  }
  const control = (text, label, action, extra = '') => html('button', { type: 'button', class: `rmap-control ${extra}`, 'aria-label': label, title: label, onClick: action }, text);
  const header = html('div', { class: 'rmap-header' }, [
    context.focusBranch ? control('Whole tree', 'Return to the whole research tree', () => select('objective', tree.objective.id), 'rmap-back') : html('span', { class: 'rmap-caption' }, 'Select a path or task dot'),
    html('div', { class: 'rmap-controls', role: 'group', 'aria-label': 'Map viewport' }, [control('−', 'Zoom out', () => setZoom(zoom / 1.2)), zoomLabel, control('+', 'Zoom in', () => setZoom(zoom * 1.2)), control('Fit', 'Fit map', fit, 'rmap-fit')])
  ]);
  const legendStatuses = ['open', 'active', 'review', 'verified', 'blocked', 'expired'];
  for (const extra of ['closed', 'unknown']) if (tree.tasks.some((task) => statusOf(task, now) === extra)) legendStatuses.push(extra);
  const footer = html('div', { class: 'rmap-footer' }, [
    html('div', { class: 'rmap-key', 'aria-label': 'Task status key' }, legendStatuses.map((status) => html('span', {}, [html('i', { class: `rmap-status-dot rmap-status-${status}`, 'aria-hidden': 'true' }), labels[status]]))),
    !context.focusBranch && tree.branches.some((branch) => branch.tasks.length > 3) ? html('span', { class: 'rmap-hint' }, 'First 3 task dots per approach. Select an approach for all tasks.') : null
  ]);
  const dependencyArea = html('div', { class: 'rmap-dependencies' });
  if (context.selectedTask?.dependsOn?.length) {
    dependencyArea.append(html('span', { class: 'rmap-dependency-label' }, 'Prerequisites'), ...context.selectedTask.dependsOn.map((id) => {
      const task = tree.tasks.find((item) => item.id === id);
      return html('button', { type: 'button', class: 'rmap-dependency-button', title: task?.title || id, onClick: () => select('task', id) }, task?.title || id);
    }));
  }
  root.append(header, stage, footer);
  if (dependencyArea.childElementCount) root.append(dependencyArea);

  function addNode(type, id, box, content, classes, label) {
    const selected = selection?.type === type && selection.id === id;
    const button = html('button', { type: 'button', class: `rmap-node rmap-node-${type} ${classes || ''}${selected ? ' is-selected' : ''}`, 'data-research-node': keyOf(type, id), 'aria-label': label, 'aria-pressed': selected ? 'true' : 'false', title: label, onClick: () => select(type, id) }, content);
    nodeButtons.set(keyOf(type, id), button);
    canvas.append(svg('foreignObject', { x: box.x, y: box.y, width: box.width, height: box.height, class: 'rmap-node-host' }, button));
  }
  function draw(width) {
    const focusedKey = canvas.contains(document.activeElement) ? document.activeElement?.getAttribute('data-research-node') : null;
    const rememberedLeft = saved.left, rememberedTop = saved.top;
    const narrow = width < 760;
    const changedKind = saved.narrow !== null && saved.narrow !== narrow;
    layout = context.focusBranch ? focusedLayout(context.focusBranch, tree, width, narrow, now) : overviewLayout(tree, width, narrow);
    root.classList.toggle('rmap-narrow', narrow);
    root.classList.toggle('rmap-overview', layout.overview);
    canvas.setAttribute('viewBox', `0 0 ${layout.width} ${layout.height}`);
    canvas.replaceChildren(svg('defs', {}, svg('marker', { id: marker, markerWidth: 6, markerHeight: 6, refX: 5, refY: 3, orient: 'auto' }, svg('path', { d: 'M0,0 L6,3 L0,6 Z', class: 'rmap-dependency-arrow' }))));
    nodeButtons.clear();
    const r = layout.root;
    for (const branch of layout.branches) {
      const emphasized = branch.data.id === context.focusBranch?.id || branch.data.tasks.some((task) => context.related.has(task.id));
      let path;
      if (layout.overview) {
        const isLeft = branch.side === 'left';
        const fromX = isLeft ? r.x : r.x + r.width;
        const toX = isLeft ? branch.x + (narrow ? 0 : branch.width) : branch.x + (narrow ? branch.width : 0);
        path = narrow ? `M${fromX},${r.y + r.height / 2} C${isLeft ? 7 : layout.width - 7},${r.y + r.height / 2} ${isLeft ? 7 : layout.width - 7},${branch.y + branch.height / 2} ${toX},${branch.y + branch.height / 2}` : cubic(fromX, r.y + r.height / 2, toX, branch.y + branch.height / 2);
      } else path = narrow ? `M${r.x + r.width / 2},${r.y + r.height} L${branch.x + branch.width / 2},${branch.y}` : cubic(r.x + r.width, r.y + r.height / 2, branch.x, branch.y + branch.height / 2);
      canvas.append(svg('path', { d: path, class: `rmap-connector rmap-color-${branch.color}${emphasized ? ' is-emphasized' : ''}` }));
      for (const task of layout.tasks.filter((item) => item.branchId === branch.data.id)) {
        let taskPath;
        if (layout.overview && narrow) taskPath = `M${branch.x + branch.width / 2},${branch.y + branch.height} L${task.x + 11},${task.y + 11}`;
        else if (layout.overview) taskPath = cubic(branch.side === 'left' ? branch.x : branch.x + branch.width, branch.y + branch.height / 2, task.x + 11, task.y + 11);
        else if (narrow) taskPath = `M${branch.x + 12},${branch.y + branch.height} L${branch.x + 12},${task.y + task.height / 2 - 10} Q${branch.x + 12},${task.y + task.height / 2} ${task.x},${task.y + task.height / 2}`;
        else taskPath = cubic(branch.x + branch.width, branch.y + branch.height / 2, task.x, task.y + task.height / 2);
        canvas.append(svg('path', { d: taskPath, class: `rmap-connector rmap-color-${branch.color}${context.related.has(task.data.id) ? ' is-emphasized' : ''}` }));
      }
    }
    const dependencies = new Set(context.selectedTask?.dependsOn || []);
    const target = layout.tasks.find((task) => task.data.id === context.selectedTask?.id);
    if (target) for (const from of layout.tasks.filter((task) => dependencies.has(task.data.id) && task !== target)) {
      canvas.append(svg('path', { d: `M${from.x + from.width},${from.y + from.height / 2} C${layout.width - 8},${from.y + from.height / 2} ${layout.width - 8},${target.y + target.height / 2} ${target.x + target.width + 2},${target.y + target.height / 2}`, class: 'rmap-dependency-link', 'marker-end': `url(#${marker})` }));
    }
    const isDeterminantal = tree.objective.id === 'marcus-de-oliveira';
    const rootContent = isDeterminantal
      ? [witnessDiagram(), html('span', {class:'rmap-figure-caption'}, '3 × 3 witness')]
      : html('span', {class:'rmap-root-title'}, tree.objective.title || tree.objective.id);
    addNode('objective', tree.objective.id, r, rootContent, isDeterminantal ? 'rmap-witness-node' : '',
      isDeterminantal ? 'Marcus–de Oliveira conjecture. Exact three by three example: the determinant lies inside the permutation hull.' : tree.objective.title || tree.objective.id);
    for (const branch of layout.branches) {
      const live = branch.data.activeCount > 0;
      const { label, activeCount } = branch.data;
      const taskCount = branch.data.tasks.length;
      const leaseLabel = `${activeCount} lease${activeCount === 1 ? '' : 's'} held`;
      const content = html('span', { class: 'rmap-branch-main' }, [
        html('i', { class: 'rmap-branch-seed', 'aria-hidden': 'true' }),
        html('span', { class: 'rmap-branch-label' }, label),
        html('span', { class: 'rmap-branch-count' }, taskCount)
      ]);
      addNode('branch', branch.data.id, branch,
        [content, live ? html('span', { class: 'rmap-active-label' }, leaseLabel) : null],
        `rmap-color-${branch.color}${live ? ' has-active-lease' : ''}`,
        `${label}, ${taskCount} task${taskCount === 1 ? '' : 's'}${live ? `, ${leaseLabel}` : ''}. Select to explore.`);
    }
    for (const task of layout.tasks) {
      const status = statusOf(task.data, now);
      const related = context.related.has(task.data.id);
      const label = `${task.data.title || task.data.id}. ${taskMeta(task.data, tree, now)}.`;
      if (task.dot) addNode('task', task.data.id, task, html('span', { class: 'rmap-dot-core', 'aria-hidden': 'true' }), `rmap-task-endpoint rmap-status-${status}${related ? ' is-related' : ''}`, label);
      else addNode('task', task.data.id, task, [html('span', { class: 'rmap-task-title' }, task.data.title || task.data.id), html('span', { class: 'rmap-task-meta' }, [html('i', { class: `rmap-status-dot rmap-status-${status}`, 'aria-hidden': 'true' }), html('span', { class: 'rmap-status-label' }, taskMeta(task.data, tree, now)), dependencies.has(task.data.id) ? html('span', { class: 'rmap-prerequisite' }, 'Prerequisite') : null])], `rmap-color-${task.color}${related ? ' is-related' : ''}${dependencies.has(task.data.id) ? ' is-prerequisite' : ''}`, label);
    }
    if (!tree.branches.length) canvas.append(svg('text', { x: layout.width / 2, y: layout.height - 28, 'text-anchor': 'middle', class: 'rmap-empty-label' }, 'Approaches appear when agents add tasks.'));
    if (saved.zoom === null || changedKind) zoom = Math.max(.9, Math.min(1, (width - 2) / layout.width, layout.overview ? (stage.clientHeight - 2) / layout.height : 1));
    setZoom(zoom, false);
    if (firstDraw && !changedKind) { stage.scrollLeft = rememberedLeft; stage.scrollTop = rememberedTop; }
    if (changedKind) { stage.scrollLeft = 0; stage.scrollTop = 0; }
    firstDraw = false; saved.narrow = narrow; save();
    const active = document.activeElement;
    const mayRestorePending = active === document.body || active === document.documentElement || Boolean(active?.closest?.('#research-canvas')) || canvas.contains(active);
    const requestedFocus = pendingFocus?.objective === tree.objective.id ? pendingFocus.key : null;
    const restore = focusedKey || (mayRestorePending ? requestedFocus : null);
    if (restore) nodeButtons.get(restore)?.focus({ preventScroll: true });
    if (requestedFocus) pendingFocus = null;
  }
  stage.addEventListener('scroll', save, { passive: true });
  let drag = null;
  stage.addEventListener('pointerdown', (event) => {
    if (event.pointerType !== 'mouse' || event.button !== 0 || event.target.closest('button')) return;
    drag = { x: event.clientX, y: event.clientY, left: stage.scrollLeft, top: stage.scrollTop };
    stage.setPointerCapture(event.pointerId); stage.classList.add('is-panning');
  });
  stage.addEventListener('pointermove', (event) => { if (drag) { stage.scrollLeft = drag.left - event.clientX + drag.x; stage.scrollTop = drag.top - event.clientY + drag.y; } });
  const stop = () => { drag = null; stage.classList.remove('is-panning'); };
  for (const event of ['pointerup', 'pointercancel', 'lostpointercapture']) stage.addEventListener(event, stop);
  const disconnect = () => observer?.disconnect();
  activeInstance = { save, disconnect };
  requestAnimationFrame(() => {
    if (!stage.isConnected) return;
    currentWidth = stage.clientWidth; draw(currentWidth);
    observer = new ResizeObserver(() => {
      if (!stage.isConnected) { disconnect(); return; }
      if (Math.abs(stage.clientWidth - currentWidth) > 2) { currentWidth = stage.clientWidth; draw(currentWidth); }
    });
    observer.observe(stage);
  });
  return root;
}
