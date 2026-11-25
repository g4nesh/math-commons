import test from 'node:test';
import assert from 'node:assert/strict';

// Minimal DOM surface: exercise the actual page renderer and control handlers
// without adding a browser dependency to the stdlib-only test suite.
class Element {
  constructor(tag, attrs = {}, children = []) { this.tag = tag; this.attrs = attrs; this.children = []; this.handlers = {}; this.append(children); }
  append(...children) { this.children.push(...children.flat(Infinity).filter(value => value !== null && value !== undefined && value !== false)); }
  replaceChildren(...children) { this.children = []; this.append(children); }
  setAttribute(key, value) { this.attrs[key] = String(value); }
  addEventListener(key, handler) { this.handlers[key] = handler; }
  get textContent() { return this.children.map(child => child instanceof Element ? child.textContent : String(child)).join(''); }
}
globalThis.Node = Element;
globalThis.window = { innerWidth: 1000, addEventListener() {} };
globalThis.document = { createElementNS: (_ns, tag) => new Element(tag), createTextNode: text => text };
const { renderContributions } = await import('../public/contributions.js');
const el = (tag, attrs = {}, children = []) => new Element(tag, attrs, children);
function find(root, predicate) {
  return [root, ...root.children.filter(child => child instanceof Element).flatMap(child => find(child, predicate))].filter(predicate);
}
const hasClass = name => node => (node.attrs.class || '').split(' ').includes(name);
const now = Date.parse('2026-11-25T12:00:00Z');
const receipt = (overrides = {}) => ({
  id: 'codex-run', owner: 'alice', agentId: 'alice-agent', provider: 'codex', model: 'unspecified',
  runId: 'run', taskId: 'task', inputTokens: 100, outputTokens: 20, cachedInputTokens: 50,
  startedAt: '2026-11-25T00:00:00Z', finishedAt: '2026-11-25T00:05:00Z',
  recordedAt: '2026-11-25T00:06:00Z', source: 'runtime', ...overrides
});
function render(usage) {
  return renderContributions({ state: { usage }, now, el, repoUrl: 'https://github.com/example/repo',
    button: (text, onClick, className) => el('button', { onClick, class: className }, text),
    externalLink: (text, href, className) => el('a', { href, class: className }, text),
    onPreview() {}, onLive() {} });
}

test('page shows exact provider subtotals beside combined total and centers a single-day point', () => {
  const root = render([receipt(), receipt({ id: 'claude-run', provider: 'claude', inputTokens: 80, outputTokens: 10, cachedInputTokens: 20 })]);
  find(root, node => node.tag === 'button' && node.textContent === 'All time')[0].attrs.onClick();
  const total = find(root, hasClass('usage-total'))[0];
  assert.match(total.textContent, /210 input \+ output tokens/);
  assert.deepEqual(find(total, hasClass('usage-provider-total')).map(node => node.textContent), ['Claude · 90', 'Codex · 120']);
  const plot = find(root, hasClass('usage-plot'))[0];
  assert.equal(find(plot, node => node.tag === 'text' && node.attrs.y === '298').length, 1);
  const point = find(plot, node => node.tag === 'circle' && node.attrs.r === '3')[0];
  assert.equal(point.attrs.cx, '439');
  assert.ok(Number.isFinite(Number(point.attrs.cy)));
  assert.doesNotMatch(root.textContent, /These receipts need a closer look/);

  // Owner selection removes its provider subtotal instead of leaving stale totals.
  find(root, node => node.tag === 'button' && node.textContent === '@alice')[0].attrs.onClick();
  assert.equal(find(root, hasClass('usage-provider-total')).length, 0);
  assert.match(find(root, hasClass('usage-total'))[0].textContent, /No receipts/);
  find(root, node => node.tag === 'button' && node.textContent === 'Show everyone')[0].attrs.onClick();
});

test('unmetered page does not invent a provider breakdown', () => {
  const root = render([]);
  assert.equal(find(root, hasClass('usage-provider-total')).length, 0);
  assert.match(find(root, hasClass('usage-total'))[0].textContent, /Not metered yet/);
});
