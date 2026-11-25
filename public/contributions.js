import { aggregateUsage, normalizeReceipts, ownerActivity, exampleUsage, providerLabel } from './usage.js';
const settings = { period: '30d', mode: 'cumulative', owners: null };
let mountedView = null, resizeTimer;
window.addEventListener('resize', () => {
  clearTimeout(resizeTimer);
  resizeTimer = setTimeout(() => { if (mountedView?.root.isConnected) mountedView.draw(); }, 150);
});
const compact = (number) => new Intl.NumberFormat(undefined, { notation: 'compact', maximumFractionDigits: 1 }).format(number);
const whole = (number) => new Intl.NumberFormat().format(number);
const dateLabel = (date) => new Intl.DateTimeFormat(undefined, { month: 'short', day: 'numeric', timeZone: 'UTC' }).format(new Date(`${date}T00:00:00Z`));
const NS = 'http://www.w3.org/2000/svg';
const nodes = (children) => children.flat(Infinity).filter((child) => child !== null && child !== undefined && child !== false);
const replace = (node, ...children) => node.replaceChildren(...nodes(children));
function svg(tag, attributes = {}, children = []) {
  const node = document.createElementNS(NS, tag);
  for (const [key, value] of Object.entries(attributes)) node.setAttribute(key, String(value));
  for (const child of nodes([children])) node.append(child instanceof Node ? child : document.createTextNode(String(child)));
  return node;
}
function doodle(className = 'usage-doodle') {
  return svg('svg', { viewBox: '0 0 64 64', class: className, 'aria-hidden': 'true' }, [svg('path', { d: 'M31 5 L34 23 L50 13 L40 29 L60 33 L40 37 L50 53 L34 43 L29 61 L26 42 L8 51 L19 36 L2 30 L21 26 L12 10 L27 21 Z' })]);
}
function color(owner, allOwners) { return `usage-color-${Math.max(0, allOwners.indexOf(owner)) % 7}`; }
function dataFor(state, example, options = {}) { return aggregateUsage(example ? exampleUsage() : state?.usage, { ...settings, ...options }); }
export function usageTeaser({ state, el, button, onOpen, example, now }) {
  let data;
  try { data = dataFor(state, example, { period: 'all', owners: null, now: now || Date.now() }); }
  catch { return el('section', { class: 'usage-teaser' }, [doodle('usage-teaser-doodle'), el('div', {}, [el('h2', {}, 'Contributions need a closer look.'), el('p', {}, 'Receipt totals could not be safely calculated. The rest of the workspace is available.')]), button('View contributions', onOpen, 'button button-secondary')]); }
  return el('section', { class: 'usage-teaser', 'aria-label': 'Token contributions' }, [
    doodle('usage-teaser-doodle'),
    el('div', {}, [el('p', { class: 'usage-handnote' }, 'a little effort, added up'), el('h2', {}, 'Tokens, transparently.'), el('p', {}, data.total === null ? 'The work is already visible. Token totals begin with the first runtime receipt.' : `${whole(data.total)} ${example ? 'synthetic example' : 'recorded'} tokens from ${data.availableOwners.length} contributors. Follow the work as it adds up.`)]),
    button('View contributions', onOpen, 'button button-secondary')
  ]);
}
export function renderContributions(context) {
  const { state, el, button, externalLink, example, onPreview, onLive, repoUrl } = context;
  const root = el('section', { class: 'contributions-view' });
  function draw() {
    try { drawUnsafe(); }
    catch { const audit = normalizeReceipts(example ? exampleUsage() : state?.usage, context.now || Date.now()); replace(root, el('div', { class: 'usage-data-error', role: 'status' }, [el('h2', {}, 'These receipts need a closer look.'), el('p', {}, 'The contribution chart could not safely aggregate this snapshot. Its totals may exceed the supported integer precision or date range. Other workspace views remain available.'), el('p', {}, `${audit.invalidCount} invalid or conflicting records and ${audit.duplicateCount} duplicate records were excluded.`), externalLink('Inspect the shared ledger', `${repoUrl}/blob/coordination/coordination/state.json`, 'quiet-link inline-link')])); }
  }
  function drawUnsafe() {
    const referenceTime = example ? Date.now() : context.now || Date.now();
    const all = dataFor(state, example, { owners: null, now: referenceTime });
    const data = dataFor(state, example, { now: referenceTime });
    const hasAny = all.allReceiptCount > 0;
    const activity = example ? [] : ownerActivity(state);
    const heading = el('div', { class: 'usage-heading' }, el('div', {}, [
      el('h1', {}, 'Token contributions'),
      el('p', { class: 'usage-intro' }, 'Tokens recorded by each contributor’s agents, grouped by the day a run finished.')
    ]));
    const controls = el('div', { class: 'usage-controls' }, [
      el('div', { class: 'usage-window-label' }, [el('span', { class: 'usage-small-title' }, 'Recorded contributions'), el('span', { class: 'meta' }, `${dateLabel(data.startDay)} – ${dateLabel(data.endDay)} · UTC`)]),
      el('div', { class: 'usage-control-groups' }, [segmented([['7d', '7 days'], ['30d', '30 days'], ['all', 'All time']], settings.period, 'Date window', (value) => { settings.period = value; draw(); }), segmented([['cumulative', 'Cumulative'], ['daily', 'Daily']], settings.mode, 'Chart mode', (value) => { settings.mode = value; draw(); })])
    ]);
    const stats = el('div', { class: 'usage-stats' }, [
      metric(data.total === null ? (hasAny ? 'No receipts' : 'Not metered yet') : compact(data.total), 'Recorded tokens', data.total === null ? 'Historical usage is unknown' : `${whole(data.total)} input + output tokens`, 'usage-total', data.providers),
      metric(data.input === null ? '—' : compact(data.input), 'Input tokens', data.cachedInput === null ? 'No input receipts in this window' : `${whole(data.cachedInput)} cached · included in input`),
      metric(data.output === null ? '—' : compact(data.output), 'Output tokens', data.receiptCount ? `${whole(data.receiptCount)} run receipts in this window` : 'No output receipts in this window'),
      metric(data.series.filter((row) => row.receiptCount).length || '—', 'Metered contributors', example ? 'Synthetic example identities' : 'Grouped by GitHub owner')
    ]);
    const chartPanel = el('section', { class: 'usage-chart-panel', 'aria-labelledby': 'usage-chart-title' }, [
      el('div', { class: 'usage-chart-title-row' }, [el('div', {}, [el('h2', { id: 'usage-chart-title' }, settings.mode === 'cumulative' ? 'Cumulative tokens' : 'Daily tokens'), el('p', {}, settings.mode === 'cumulative' ? 'Cumulative recorded tokens within the selected window.' : 'Recorded tokens by finish date. Blank days have no receipts.')]), example ? el('span', { class: 'usage-stamp' }, 'Synthetic example') : null]),
      hasAny ? ownerControls(all.availableOwners, data.availableOwners) : null,
      data.receiptCount ? chart(data, all.availableOwners) : emptyChart(hasAny),
      el('div', { class: 'usage-chart-foot' }, [el('span', {}, 'Input + output. Cached input is counted once.'), externalLink('Receipt format', `${repoUrl}/blob/main/docs/usage.md`, 'quiet-link inline-link')])
    ]);
    replace(root, heading,
      example ? el('div', { class: 'usage-example-banner', role: 'status' }, [el('div', {}, [el('strong', {}, 'Example chart · synthetic data'), el('p', {}, 'These fictional contributors and receipts illustrate the chart. They are not workspace usage.')]), button('Return to live data', () => { settings.owners = null; onLive(); }, 'button button-secondary button-small')]) : null,
      controls, data.pendingCount ? el('p', { class: 'usage-pending', role: 'status' }, `${data.pendingCount} receipt${data.pendingCount === 1 ? '' : 's'} awaiting the recorded completion time (up to five minutes of clock skew). Counts join the chart after that time.`) : null, data.invalidCount || data.duplicateCount ? el('p', { class: 'usage-pending', role: 'status' }, `${data.invalidCount} invalid or conflicting records excluded; ${data.duplicateCount} duplicate receipts ignored. See provenance below.`) : null, stats, chartPanel, heatmap(data), contributors(all, activity), provenance(data));
  }
  function metric(value, label, detail, className = '', providers = []) {
    return el('div', { class: `usage-stat ${className}` }, [el('span', { class: 'usage-stat-label' }, label), el('strong', {}, value), el('p', {}, detail), ...providers.map(row => el('p', { class: 'usage-provider-total' }, `${providerLabel(row.provider)} · ${whole(row.total)}`))]);
  }
  function segmented(options, value, label, onChange) {
    return el('div', { class: 'usage-segmented', role: 'group', 'aria-label': label }, options.map(([key, text]) => el('button', { type: 'button', class: key === value ? 'is-selected' : '', 'aria-pressed': key === value ? 'true' : 'false', onClick: () => onChange(key) }, text)));
  }
  function ownerControls(owners) {
    return el('div', { class: 'usage-owner-controls', role: 'group', 'aria-label': 'Show contributor series' }, [
      ...owners.map((owner) => {
        const selected = settings.owners === null || settings.owners.includes(owner);
        return el('button', { type: 'button', class: `usage-owner-chip ${color(owner, owners)}${selected ? ' is-selected' : ''}`, 'aria-pressed': selected ? 'true' : 'false', onClick: () => { const next = new Set(settings.owners === null ? owners : settings.owners); if (next.has(owner)) next.delete(owner); else next.add(owner); settings.owners = next.size === owners.length ? null : [...next]; draw(); } }, [el('span', { class: 'usage-color-dot', 'aria-hidden': 'true' }), `@${owner}`]);
      }), settings.owners !== null ? button('Show everyone', () => { settings.owners = null; draw(); }, 'button button-plain usage-reset') : null
    ]);
  }
  function emptyChart(filtered) {
    return el('div', { class: 'usage-empty-chart' }, [
      el('div', { class: 'usage-empty-axis', 'aria-hidden': 'true' }),
      el('div', { class: 'usage-empty-note' }, [el('h3', {}, filtered ? 'No receipts in this selection.' : 'No token receipts yet.'), el('p', {}, filtered ? 'Choose another date window or show everyone. A blank chart means missing receipts, not zero historical usage.' : 'Agents have already been contributing. Earlier runs were not metered, so we’re leaving their token totals unknown.'), el('div', { class: 'usage-empty-actions' }, [filtered ? button('Show all recorded work', () => { settings.period = 'all'; settings.owners = null; draw(); }, 'button button-secondary') : button('Preview example chart', () => { settings.owners = null; settings.period = 'all'; onPreview(); }, 'button usage-coral-button'), externalLink('How receipts work', `${repoUrl}/blob/main/docs/usage.md`, 'quiet-link inline-link')])]),
      doodle('usage-empty-star')
    ]);
  }
  function chart(data, allOwners) {
    const width = Math.max(300, Math.min(1000, window.innerWidth - (window.innerWidth < 740 ? 80 : 140))), height = 310;
    const left = 48, right = width < 500 ? 16 : 30, top = 22, bottom = 40, plotWidth = width - left - right, plotHeight = height - top - bottom;
    const key = settings.mode === 'cumulative' ? 'cumulative' : 'total';
    const values = data.series.flatMap((series) => series.days.map((day) => day[key] ?? 0));
    const maximum = Math.max(4, Math.ceil(Math.max(0, ...values) * 1.12 / 4) * 4);
    const x = (index) => left + plotWidth * (data.days.length === 1 ? 0.5 : index / (data.days.length - 1));
    const y = (value) => top + plotHeight * (1 - value / maximum);
    const plot = svg('svg', { class: 'usage-plot', viewBox: `0 0 ${width} ${height}`, role: 'img', 'aria-labelledby': 'usage-svg-title usage-svg-description' }, [
      svg('title', { id: 'usage-svg-title' }, `${settings.mode === 'cumulative' ? 'Cumulative' : 'Daily'} recorded tokens per GitHub contributor`),
      svg('desc', { id: 'usage-svg-description' }, `${data.receiptCount} receipts totaling ${whole(data.total)} tokens between ${data.startDay} and ${data.endDay} UTC. Use the day inspector or expand the data table for exact values.`)
    ]);
    for (let index = 0; index <= 4; index++) {
      const value = maximum * index / 4;
      plot.append(svg('line', { x1: left, x2: width - right, y1: y(value), y2: y(value), class: 'usage-grid-line' }), svg('text', { x: left - 10, y: y(value) + 4, 'text-anchor': 'end', class: 'usage-axis-label' }, compact(Math.round(value))));
    }
    const tickCount = width < 500 ? 3 : 5;
    const dateIndices = [...new Set(Array.from({ length: tickCount }, (_, index) => Math.round(index * (data.days.length - 1) / (tickCount - 1))))];
    for (const [index, dayIndex] of dateIndices.entries()) {
      plot.append(svg('text', { x: x(dayIndex), y: height - 12, 'text-anchor': dateIndices.length === 1 ? 'middle' : index === 0 ? 'start' : index === dateIndices.length - 1 ? 'end' : 'middle', class: 'usage-axis-label' }, dateLabel(data.days[dayIndex].date)));
    }
    for (const series of data.series) {
      let started = false, path = '';
      series.days.forEach((day, index) => {
        if (settings.mode === 'daily' && day.total === null) { started = false; return; }
        const value = day[key] ?? 0;
        path += `${started ? 'L' : 'M'}${x(index).toFixed(2)},${y(value).toFixed(2)} `; started = true;
      });
      plot.append(svg('path', { d: path, class: `usage-series-line ${color(series.owner, allOwners)}` }));
      if (settings.mode === 'daily') series.days.forEach((day, index) => { if (day.total !== null) plot.append(svg('circle', { cx: x(index), cy: y(day.total), r: 2.5, class: `usage-series-point ${color(series.owner, allOwners)}` })); });
      else if (series.receiptCount) plot.append(svg('circle', { cx: x(series.days.length - 1), cy: y(series.days.at(-1).cumulative), r: 3, class: `usage-series-point ${color(series.owner, allOwners)}` }));
    }
    const marker = svg('g', { class: 'usage-marker', 'aria-hidden': 'true' }); plot.append(marker);
    const readout = el('div', { class: 'usage-day-readout', 'aria-live': 'polite' });
    const range = el('input', { type: 'range', min: 0, max: data.days.length - 1, value: data.days.length - 1, class: 'usage-day-slider', 'aria-label': 'Inspect UTC day on the contribution chart' });
    let lastIndex = -1;
    function inspect(index) {
      if (lastIndex === index) return; lastIndex = index;
      const day = data.days[index]; range.value = index;
      range.setAttribute('aria-valuetext', `${day.date}: ${day.total === null ? 'No receipts; usage unknown' : `${whole(day.total)} recorded tokens`}`);
      marker.replaceChildren(svg('line', { x1: x(index), x2: x(index), y1: top, y2: height - bottom, class: 'usage-guide-line' }));
      readout.replaceChildren(el('strong', {}, `${dateLabel(day.date)} UTC`), ...data.series.map((series) => { const entry = series.days[index]; const value = settings.mode === 'cumulative' ? `${whole(entry.cumulative)} cumulative` : entry.total === null ? 'no receipt' : `${whole(entry.total)} tokens`; if (entry[key] !== null) marker.append(svg('circle', { cx: x(index), cy: y(entry[key]), r: 4, class: `usage-series-point ${color(series.owner, allOwners)}` })); return el('span', { class: `usage-day-value ${color(series.owner, allOwners)}` }, [el('span', { class: 'usage-color-dot', 'aria-hidden': 'true' }), `@${series.owner}: ${value}`]); }));
    }
    range.addEventListener('input', (event) => inspect(Number(event.target.value)));
    plot.addEventListener('pointermove', (event) => { const bounds = plot.getBoundingClientRect(); const coordinate = (event.clientX - bounds.left) / bounds.width * width; inspect(Math.max(0, Math.min(data.days.length - 1, Math.round((coordinate - left) / plotWidth * (data.days.length - 1))))); });
    inspect(data.days.length - 1);
    return el('div', { class: 'usage-chart-body' }, [plot, el('label', { class: 'usage-inspector' }, [el('span', {}, 'Inspect a day'), range]), readout, dataTable(data)]);
  }
  function dataTable(data) {
    return el('details', { class: 'usage-data-details' }, [el('summary', {}, 'View accessible data table'), el('div', { class: 'usage-table-scroll', tabindex: 0 }, el('table', { class: 'usage-data-table' }, [el('caption', {}, 'Recorded tokens only. An em dash means no receipt, not zero usage.'), el('thead', {}, el('tr', {}, [el('th', { scope: 'col' }, 'UTC day'), ...data.series.map((series) => el('th', { scope: 'col' }, `@${series.owner}`)), el('th', { scope: 'col' }, 'Daily total'), el('th', { scope: 'col' }, 'Receipts')])), el('tbody', {}, data.days.map((day, index) => el('tr', {}, [el('th', { scope: 'row' }, day.date), ...data.series.map((series) => el('td', {}, series.days[index].total === null ? '—' : whole(series.days[index].total))), el('td', {}, day.total === null ? '—' : whole(day.total)), el('td', {}, day.receiptCount)])))]))]);
  }
  function heatmap(data) {
    // Daily cells use the selected date window; no-receipt cells are explicitly unknown.
    const limit = 366;
    const days = data.days.slice(-limit);
    const max = Math.max(1, ...days.map((day) => day.total || 0));
    const firstWeekday = new Date(`${days[0].date}T00:00:00Z`).getUTCDay();
    const cells = [];
    for (let index = 0; index < firstWeekday; index++) cells.push(el('span', { class: 'usage-heat-cell usage-heat-padding', 'aria-hidden': 'true' }));
    for (const day of days) {
      const intensity = day.total === null ? 'unknown' : day.total === 0 ? 'zero' : Math.max(1, Math.ceil(day.total / max * 4));
      const label = `${day.date} UTC: ${day.total === null ? 'no receipt; token usage unknown' : `${whole(day.total)} recorded tokens from ${day.receiptCount} receipts`}`;
      cells.push(el('button', { type: 'button', class: `usage-heat-cell usage-heat-${intensity}`, title: label, 'aria-label': label, onClick: () => { detail.textContent = label; } }));
    }
    const detail = el('p', { class: 'usage-heat-detail', role: 'status' }, data.receiptCount ? 'Choose a day to inspect its receipts.' : 'Earlier work is unmetered. Blank days do not mean nobody contributed.');
    return el('section', { class: 'usage-heatmap-panel' }, [el('div', { class: 'usage-section-heading' }, [el('div', {}, [el('h2', {}, 'A little, often.'), el('p', {}, `Daily receipt activity · UTC${data.days.length > limit ? ' · last 366 days of selected range' : ''}`)]), el('span', { class: 'usage-handnote' }, 'small steps count')]), el('div', { class: 'usage-heatmap-scroll', tabindex: 0, 'aria-label': 'Daily receipt heatmap; scroll horizontally if needed' }, el('div', { class: 'usage-heatmap-grid' }, cells)), el('div', { class: 'usage-heat-legend' }, [el('span', { class: 'usage-heat-cell usage-heat-unknown', 'aria-hidden': 'true' }), el('span', {}, 'No receipt'), el('span', { class: 'usage-heat-cell usage-heat-zero', 'aria-hidden': 'true' }), el('span', {}, 'Recorded zero'), el('span', { class: 'usage-heat-less' }, 'Less'), ...[1, 2, 3, 4].map((level) => el('span', { class: `usage-heat-cell usage-heat-${level}`, 'aria-hidden': 'true' })), el('span', {}, 'More')]), detail]);
  }
  function contributors(all, activity) {
    const owners = [...new Set([...all.availableOwners, ...activity.map((entry) => entry.owner)])];
    const rows = owners.map((owner) => ({ owner, usage: all.series.find((entry) => entry.owner === owner), activity: activity.find((entry) => entry.owner === owner) })).sort((a, b) => (b.usage?.total || 0) - (a.usage?.total || 0) || a.owner.localeCompare(b.owner));
    return el('section', { class: 'usage-contributors' }, [el('div', { class: 'usage-section-heading' }, [el('div', {}, [el('h2', {}, 'Everyone brings something.'), el('p', {}, 'Token receipts use the selected date window. Ledger activity below is all time.')]), el('span', { class: 'usage-count-note' }, `${owners.length} contributor${owners.length === 1 ? '' : 's'}`)]),
      rows.length ? el('div', { class: 'usage-contributor-list' }, rows.map(({ owner, usage, activity: ledger }) => el('article', { class: 'usage-contributor-row' }, [el('div', { class: 'usage-contributor-identity' }, [el('span', { class: `usage-person-icon ${color(owner, all.availableOwners)}`, 'aria-hidden': 'true' }, owner.slice(0, 2).toUpperCase()), el('div', {}, [example ? el('strong', {}, `@${owner}`) : externalLink(`@${owner}`, `https://github.com/${owner}`), el('p', {}, example ? 'Fictional example contributor' : `${ledger?.agents || 0} registered agents · ${ledger?.checkpoints || 0} checkpoints · ${ledger?.reviews || 0} reviews`)])]), el('div', { class: 'usage-contributor-total' }, [el('strong', {}, usage?.total === null || !usage ? 'Not metered' : whole(usage.total)), el('span', {}, usage?.receiptCount ? `${usage.receiptCount} receipts · ${usage.manualCount} manual` : 'Historical tokens unknown')])]))) : el('p', { class: 'usage-no-contributors' }, 'Contributors appear here when agents join or a receipt is published.'), el('p', { class: 'usage-note' }, 'Checkpoints and messages show participation. They are never converted into token estimates.')]);
  }
  function provenance(data) {
    const details = el('details', { class: 'usage-provenance' }, [el('summary', {}, ['Follow the numbers ', el('span', { class: 'usage-handnote' }, '(nothing up our sleeve)')]), el('div', { class: 'usage-provenance-body' }, [el('p', {}, `${data.runtimeCount} runtime receipts · ${data.manualCount} manual receipts in this selection. Each receipt records its GitHub owner, provider, model, run ID, task, time, and input/output counts.`), el('p', {}, 'This measures only runs reported to this repository. It is not an account-wide total, invoice, cost calculation, or a comparison of model quality. Missing receipts remain unknown. Cached input is a subset of input.'), data.lastRecordedAt ? el('p', {}, `Latest receipt recorded: ${new Date(data.lastRecordedAt).toISOString()}.`) : el('p', {}, 'No token receipts have been recorded for this selection.'), data.invalidCount || data.duplicateCount ? el('p', {}, `${data.invalidCount} invalid and ${data.duplicateCount} duplicate records were excluded from aggregation.`) : null, externalLink('Read the receipt protocol', `${repoUrl}/blob/main/docs/usage.md`, 'quiet-link inline-link'), data.receipts.length ? el('div', { class: 'usage-receipt-list' }, data.receipts.slice().reverse().slice(0, 20).map((receipt) => el('article', {}, [el('strong', {}, `@${receipt.owner} · ${whole(receipt.totalTokens)} tokens`), el('p', {}, `${receipt.provider} / ${receipt.model || 'model unspecified'} · ${receipt.source} · ${receipt.finishedAt}`), el('p', { class: 'mono' }, `Run ${receipt.runId} · task ${receipt.taskId || 'not specified'}`), typeof receipt.evidence === 'string' ? externalLink('Evidence', receipt.evidence, 'quiet-link inline-link') : Array.isArray(receipt.evidence) ? receipt.evidence.map((link) => externalLink('Evidence', link, 'quiet-link inline-link')) : null]))) : null])]);
    return details;
  }
  mountedView = { root, draw };
  draw(); return root;
}
