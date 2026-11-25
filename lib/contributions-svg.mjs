import { aggregateUsage, normalizeReceipts, ownerActivity, exampleUsage } from '../public/usage.js';

const COLORS = ['#2458a6', '#a44331', '#2c7564', '#815c9d', '#a37420', '#4f788d', '#74746d'];
const PLOT = { x: 72, y: 132, width: 772, height: 164 };
const fmt = value => value.toLocaleString('en-US');
const compact = value => value >= 1e9 ? (value / 1e9).toFixed(1) + 'B' :
  value >= 1e6 ? (value / 1e6).toFixed(1) + 'M' : value >= 1000 ? (value / 1000).toFixed(1) + 'k' : String(Math.round(value));
export const escapeXml = value => String(value ?? '').replace(/[\u0000-\u0008\u000b\u000c\u000e-\u001f\ufffe\uffff]/g, '\ufffd')
  .replaceAll('&', '&amp;').replaceAll('<', '&lt;').replaceAll('>', '&gt;').replaceAll('"', '&quot;').replaceAll("'", '&apos;');
const text = (x, y, value, { size = 11, fill = '#66665d', anchor = 'start', weight = 400 } = {}) =>
  '<text x="' + x + '" y="' + y + '" font-size="' + size + '" fill="' + fill + '" text-anchor="' + anchor + '" font-weight="' + weight + '">' + escapeXml(value) + '</text>';
const line = (x1, y1, x2, y2, color = '#e5e4dd') =>
  '<line x1="' + x1 + '" y1="' + y1 + '" x2="' + x2 + '" y2="' + y2 + '" stroke="' + color + '"/>';
const shell = (description, contents) =>
  '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="420" viewBox="0 0 900 420" role="img" aria-labelledby="chart-title chart-description">' +
  '<title id="chart-title">Math Commons · token contributions</title><desc id="chart-description">' + escapeXml(description) + '</desc>' +
  '<rect x=".5" y=".5" width="899" height="419" rx="14" fill="#fffdf8" stroke="#dcdad2"/>' +
  '<g font-family="ui-monospace, SFMono-Regular, Menlo, Consolas, monospace">' +
  text(32, 36, 'Math Commons', { size: 21, weight: 600, fill: '#252821' }) +
  '<path d="M226 16 L230 24 L239 19 L235 29 L245 34 L234 35 L235 45 L227 38 L219 44 L221 34 L212 29 L223 27 Z" fill="none" stroke="#a66d59" stroke-width="1.5" stroke-linejoin="round" aria-hidden="true"/>' +
  text(32, 57, 'Token contributions', { size: 12 }) + contents + '</g></svg>';

function prepared(receipts, now) {
  const normalized = normalizeReceipts(receipts, now), totals = new Map();
  for (const receipt of normalized.receipts) {
    const sum = (totals.get(receipt.owner) || 0) + receipt.totalTokens;
    if (!Number.isSafeInteger(sum)) throw new RangeError('Owner totals exceed safe integer precision.');
    totals.set(receipt.owner, sum);
  }
  const leaders = [...totals].sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0])).slice(0, 6).map(([owner]) => owner);
  const leading = new Set(leaders);
  let otherKey = 'other-contributors';
  while (leading.has(otherKey)) otherKey += '-x';
  // Group only this in-memory chart projection before allocating day series.
  // Original receipts remain unchanged; at most seven series are constructed.
  const grouped = normalized.receipts.map(receipt => ({
    ...receipt, owner: leading.has(receipt.owner) ? receipt.owner : otherKey,
    runId: JSON.stringify([receipt.owner, receipt.provider, receipt.runId])
  }));
  const metrics = aggregateUsage(grouped, { period: 'all', now });
  const series = leaders.map(owner => metrics.series.find(item => item.owner === owner));
  const rest = metrics.series.find(item => item.owner === otherKey);
  if (rest) series.push({ ...rest, owner: 'Others (' + (totals.size - leaders.length) + ')' });
  return { metrics, series, ...Object.fromEntries(['invalidCount', 'duplicateCount', 'pendingCount'].map(key => [key, normalized[key]])) };
}

function activityText(state) {
  const owners = ownerActivity({
    agents: Array.isArray(state?.agents) ? state.agents : [],
    events: Array.isArray(state?.events) ? state.events : []
  });
  return fmt(owners.length) + ' ledger contributors · ' + fmt(owners.reduce((sum, owner) => sum + owner.checkpoints, 0)) + ' checkpoints';
}

export function renderContributionsSvg(snapshot, { now = Date.now(), demo = false } = {}) {
  if (!demo && snapshot?.source !== 'github') throw new TypeError('Expected a GitHub snapshot.');
  const state = snapshot?.state || {};
  if (!demo && state.usage !== undefined && !Array.isArray(state.usage)) {
    throw new TypeError('Invalid receipt collection.');
  }
  const { metrics, series, invalidCount, pendingCount } = prepared(demo ? exampleUsage(now) : state.usage, now);
  const hasReceipts = metrics.receiptCount > 0;
  const body = [];
  if (demo) {
    body.push('<rect x="660" y="20" width="208" height="30" rx="5" fill="#fff0cd"/>');
    body.push(text(764, 40, 'SYNTHETIC EXAMPLE', { fill: '#8c4e12', weight: 700, anchor: 'middle', size: 12 }));
  } else body.push(text(868, 35, String(snapshot.repository || 'g4nesh/math-commons').slice(0, 90), { anchor: 'end', size: 10 }));
  if (hasReceipts) {
    body.push(text(32, 94, fmt(metrics.total) + ' recorded tokens', { size: 24, weight: 600, fill: '#252821' }));
    body.push(text(32, 115, fmt(metrics.input) + ' input + ' + fmt(metrics.output) + ' output · ' + fmt(metrics.cachedInput) + ' cached input (included)', { size: 11 }));
  } else {
    body.push(text(32, 94, 'Not metered yet', { size: 24, weight: 600, fill: '#252821' }));
    body.push(text(32, 115, 'Historical token usage is unknown. Missing receipts are not zero.', { size: 11 }));
  }

  const largest = Math.max(1, ...series.map(item => item.total || 0));
  for (let index = 0; index <= 4; index++) {
    const y = PLOT.y + PLOT.height - index * PLOT.height / 4;
    body.push(line(PLOT.x, y, PLOT.x + PLOT.width, y));
    if (hasReceipts) body.push(text(PLOT.x - 12, y + 4, compact(largest * index / 4), { anchor: 'end', size: 12 }));
  }
  const lastIndex = metrics.days.length - 1;
  const xAt = index => PLOT.x + (lastIndex ? index / lastIndex : 0.5) * PLOT.width;
  const yAt = value => PLOT.y + PLOT.height - value / largest * PLOT.height;
  const dateIndices = [...new Set([0, Math.floor(lastIndex / 2), lastIndex])];
  for (const index of dateIndices) body.push(text(xAt(index), 315, metrics.days[index].date, { anchor: 'middle', size: 12 }));
  body.push(text(844, 129, 'Cumulative recorded tokens · UTC', { anchor: 'end', size: 9 }));

  if (hasReceipts) {
    for (const [index, owner] of series.entries()) {
      const color = COLORS[index];
      // Keep final SVG bounded even for long histories; preserve every endpoint.
      const stride = Math.max(1, Math.ceil((owner.days.length - 1) / 799));
      const points = owner.days.filter((_, dayIndex) => dayIndex % stride === 0 || dayIndex === lastIndex);
      const path = points.map((day, pointIndex) => {
        const originalIndex = pointIndex === points.length - 1 ? lastIndex : pointIndex * stride;
        return (pointIndex ? 'L' : 'M') + xAt(originalIndex).toFixed(2) + ',' + yAt(day.cumulative).toFixed(2);
      }).join(' ');
      body.push('<path data-owner-line="' + index + '" d="' + path + '" fill="none" stroke="' + color + '" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round"/>');
      body.push('<circle cx="' + xAt(lastIndex).toFixed(2) + '" cy="' + yAt(owner.total).toFixed(2) + '" r="3" fill="' + color + '"/>');
      const x = 32 + (index % 4) * 216, y = 342 + Math.floor(index / 4) * 20;
      body.push(line(x, y - 4, x + 14, y - 4, color));
      const label = owner.owner.length > 18 ? owner.owner.slice(0, 17) + '…' : owner.owner;
      body.push(text(x + 21, y, label + ' ' + compact(owner.total), { size: 12, fill: '#34362f' }));
    }
  } else {
    body.push(text(458, 206, 'The first runtime receipt starts this chart.', { anchor: 'middle', fill: '#484b42', size: 14 }));
    body.push(text(458, 229, 'Work can exist before its tokens are metered.', { anchor: 'middle', size: 11 }));
    body.push(text(32, 345, activityText(state) + ' · recorded work, not token estimates', { size: 11 }));
  }
  const exclusions = [invalidCount ? invalidCount + ' invalid excluded' : '', pendingCount ? pendingCount + ' future receipts pending' : ''].filter(Boolean).join(' · ');
  const source = demo ? 'Synthetic example only · no live receipts' : 'GitHub coordination receipts · contributor-reported';
  body.push(line(32, 375, 868, 375));
  body.push(text(32, 393, source + (hasReceipts ? ' · ' + metrics.runtimeCount + ' runtime / ' + metrics.manualCount + ' manual' : ''), { size: 10 }));
  const last = metrics.lastRecordedAt ? metrics.lastRecordedAt.replace('T', ' ').replace(/\.\d{3}Z$/, ' UTC') : 'none';
  body.push(text(32, 409, 'Last recorded: ' + last + (exclusions ? ' · ' + exclusions : ''), { size: 9 }));
  return shell((demo ? 'SYNTHETIC EXAMPLE. ' : '') + (hasReceipts ?
    fmt(metrics.total) + ' recorded tokens across ' + metrics.receiptCount + ' receipts. Cached input is a subset of input.' :
    'Not metered yet. Historical token usage is unknown; missing receipts are not zero.') + (exclusions ? ' ' + exclusions : ''), body.join(''));
}

export function renderUnavailableSvg({ method = false } = {}) {
  return shell('Usage temporarily unavailable. No token total is inferred.', [
    text(32, 106, method ? 'Read-only chart' : 'Usage temporarily unavailable', { size: 25, weight: 600, fill: '#252821' }),
    text(32, 143, method ? 'Use GET to view recorded contributions.' : 'The shared GitHub ledger could not be read.', { size: 13 }),
    text(32, 169, 'No zero total or example data has been substituted.', { size: 12 }),
    line(32, 375, 868, 375),
    text(32, 398, 'Math Commons · refresh later for recorded usage', { size: 11 })
  ].join(''));
}
