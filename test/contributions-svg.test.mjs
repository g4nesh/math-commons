import test from 'node:test';
import assert from 'node:assert/strict';
import { renderContributionsSvg, renderUnavailableSvg, escapeXml } from '../lib/contributions-svg.mjs';
import { createHandler } from '../api/contributions.js';

const now = Date.parse('2026-11-25T12:00:00Z');
const receipt = (change = {}) => ({
  id: 'r-1', owner: 'alice', agentId: 'alice-agent', provider: 'codex', model: 'unspecified',
  runId: 'codex:actual-run', taskId: 't', objectiveId: 'p', inputTokens: 120, outputTokens: 30,
  cachedInputTokens: 80, startedAt: '2026-11-24T00:00:00Z', finishedAt: '2026-11-24T00:05:00Z',
  recordedAt: '2026-11-24T00:06:00Z', source: 'runtime', ...change
});
const snapshot = (usage, change = {}) => ({
  source: 'github', repository: 'g4nesh/math-commons',
  state: { usage, agents: [{ owner: 'alice' }], events: [{ actorOwner: 'alice', kind: 'checkpoint' }] }, ...change
});
function response() {
  return { headers: {}, statusCode: null, body: null,
    setHeader(name, value) { this.headers[name] = value; },
    status(code) { this.statusCode = code; return this; },
    end(value) { this.body = value; return this; } };
}

test('SVG uses exact recorded totals and never adds cached tokens twice', () => {
  const svg = renderContributionsSvg(snapshot([receipt(), receipt({
    id: 'r-2', runId: 'second-run', owner: 'bob', inputTokens: 80, outputTokens: 20,
    cachedInputTokens: 10, source: 'manual'
  })]), { now });
  assert.match(svg, /250 recorded tokens/);
  assert.match(svg, /200 input \+ 50 output · 90 cached input \(included\)/);
  assert.match(svg, /1 runtime \/ 1 manual/);
  assert.match(svg, /Last recorded: 2026-11-24 00:06:00 UTC/);
  assert.match(svg, /2026-11-24/);
  assert.equal((svg.match(/data-owner-line=/g) || []).length, 2);
  assert.doesNotMatch(svg, /340 recorded tokens/);
});

test('unmetered history has no token line or false zero; ledger activity stays distinct', () => {
  const svg = renderContributionsSvg(snapshot(undefined), { now });
  assert.match(svg, /Not metered yet/);
  assert.match(svg, /Historical token usage is unknown/);
  assert.match(svg, /1 ledger contributors · 1 checkpoints/);
  assert.match(svg, /recorded work, not token estimates/);
  assert.doesNotMatch(svg, /data-owner-line=|0 recorded tokens|SYNTHETIC/);
});

test('source labels are XML escaped and untrusted owner markup is rejected', () => {
  const svg = renderContributionsSvg(snapshot([receipt({ owner: '<script>alert(1)</script>' })], {
    repository: 'repo<svg onload="bad">&\''
  }), { now });
  assert.match(svg, /repo&lt;svg onload=&quot;bad&quot;&gt;&amp;&apos;/);
  assert.match(svg, /1 invalid excluded/);
  assert.doesNotMatch(svg, /<script|<foreignObject|onload="bad"|<image|<style|href=/);
  assert.equal(escapeXml('\u0000<&>"\''), '\ufffd&lt;&amp;&gt;&quot;&apos;');
});

test('identical duplicates are counted once and conflicting duplicate runs are excluded', () => {
  const one = receipt();
  const duplicate = renderContributionsSvg(snapshot([one, { ...one }]), { now });
  assert.match(duplicate, /150 recorded tokens/);
  const conflict = renderContributionsSvg(snapshot([one, { ...one, id: 'different', inputTokens: 121 }]), { now });
  assert.match(conflict, /Not metered yet/);
  assert.match(conflict, /2 invalid excluded/);
  assert.doesNotMatch(conflict, /data-owner-line=/);
});

test('top six owners plus Others have deterministic bounded series and exact combined totals', () => {
  const records = Array.from({ length: 9 }, (_, index) => receipt({
    id: 'r-' + index, owner: 'owner-' + index, runId: 'run-' + index,
    inputTokens: (index + 1) * 100, outputTokens: 0, cachedInputTokens: 0
  }));
  const svg = renderContributionsSvg(snapshot(records), { now });
  assert.match(svg, /4,500 recorded tokens/);
  assert.equal((svg.match(/data-owner-line=/g) || []).length, 7);
  assert.match(svg, /Others \(3\) 600/);
  assert.equal(svg, renderContributionsSvg(snapshot([...records].reverse()), { now }));
  for (const color of [...svg.matchAll(/stroke="([^"]+)"/g)].map(match => match[1])) {
    assert.match(color, /^#[0-9a-f]{6}$/);
  }
});

test('unsafe counters and future receipts are not drawn as earned tokens', () => {
  const svg = renderContributionsSvg(snapshot([
    receipt({ inputTokens: '120' }),
    receipt({ id: 'future', runId: 'future', startedAt: '2026-11-25T12:00:00Z',
      finishedAt: '2026-11-25T12:01:00Z', recordedAt: '2026-11-25T12:01:00Z' })
  ]), { now });
  assert.match(svg, /Not metered yet/);
  assert.match(svg, /1 invalid excluded · 1 future receipts pending/);
  assert.doesNotMatch(svg, /data-owner-line=/);
});

test('long histories retain final endpoints with at most 800 plotted points per owner', () => {
  const svg = renderContributionsSvg(snapshot([receipt({
    startedAt: '2020-01-01T00:00:00Z', finishedAt: '2020-01-01T00:05:00Z',
    recordedAt: '2020-01-01T00:06:00Z'
  })]), { now });
  const path = svg.match(/data-owner-line="0" d="([^"]+)"/)[1];
  assert.ok(path.split(' ').length <= 800);
  assert.match(path, /L844\.00,132\.00$/);
});

test('explicit demo is prominently synthetic and live rendering never falls back to it', () => {
  const demo = renderContributionsSvg(null, { now, demo: true });
  assert.match(demo, /SYNTHETIC EXAMPLE/);
  assert.match(demo, /Synthetic example only · no live receipts/);
  assert.match(demo, /data-owner-line=/);
  assert.doesNotMatch(renderContributionsSvg(snapshot([]), { now }), /alex-example|SYNTHETIC/);
});

test('SVG endpoint shares cached live receipts with image headers', async () => {
  let reads = 0;
  const handler = createHandler({ load: async () => { reads++; return snapshot([receipt()]); }, now: () => now });
  const res = response();
  await handler({ method: 'GET', url: '/api/contributions' }, res);
  assert.equal(reads, 1);
  assert.equal(res.statusCode, 200);
  assert.equal(res.headers['Content-Type'], 'image/svg+xml; charset=utf-8');
  assert.equal(res.headers['Cache-Control'], 'public, max-age=0, s-maxage=900');
  assert.equal(res.headers['X-Content-Type-Options'], 'nosniff');
  assert.match(res.headers['Content-Security-Policy'], /default-src 'none'/);
  assert.match(res.body, /150 recorded tokens/);
});

test('failure remains a legible, uncached SVG without fabricated totals or leaked error', async () => {
  const handler = createHandler({ load: async () => { throw new Error('SECRET upstream diagnostic'); }, now: () => now });
  const res = response();
  await handler({ method: 'GET', url: '/api/contributions' }, res);
  assert.equal(res.statusCode, 200);
  assert.equal(res.headers['Cache-Control'], 'no-store');
  assert.match(res.body, /Usage temporarily unavailable/);
  assert.doesNotMatch(res.body, /SECRET|0 recorded tokens|Not metered yet|SYNTHETIC/);
  assert.equal(res.body, renderUnavailableSvg());
});

test('malformed receipt collection fails honestly instead of presenting an empty ledger', async () => {
  const handler = createHandler({ load: async () => snapshot({ unexpected: true }), now: () => now });
  const res = response();
  await handler({ method: 'GET', url: '/api/contributions' }, res);
  assert.match(res.body, /Usage temporarily unavailable/);
  assert.doesNotMatch(res.body, /Not metered yet/);
  assert.equal(res.headers['Cache-Control'], 'no-store');
});

test('demo does not read GitHub; non-GET methods are read-only failures', async () => {
  const handler = createHandler({ load: async () => { throw new Error('Must not fetch'); }, now: () => now });
  const demo = response();
  await handler({ method: 'GET', url: '/api/contributions?demo=1' }, demo);
  assert.equal(demo.statusCode, 200);
  assert.match(demo.body, /SYNTHETIC EXAMPLE/);
  const denied = response();
  await handler({ method: 'POST', url: '/api/contributions' }, denied);
  assert.equal(denied.statusCode, 405);
  assert.equal(denied.headers.Allow, 'GET');
  assert.equal(denied.headers['Cache-Control'], 'no-store');
});
