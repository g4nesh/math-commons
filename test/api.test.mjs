import test from 'node:test';
import assert from 'node:assert/strict';
import { loadSnapshot } from '../api/state.js';
test('snapshot returns GitHub provenance and preserves revision', async () => {
  const state = { schemaVersion: 1, revision: 9, objectives: [], tasks: [], agents: [], events: [] };
  const snapshot = await loadSnapshot({ fetchImpl: async url => {
    assert.match(url, /^https:\/\/raw.githubusercontent.com\/g4nesh\/math-commons\/coordination\//);
    return new Response(JSON.stringify(state));
  }, now: 0 });
  assert.deepEqual(snapshot.state, state);
  assert.equal(snapshot.source, 'github');
  assert.equal(snapshot.fetchedAt, '1970-01-01T00:00:00.000Z');
});
test('GitHub failures never masquerade as example or empty state', async () => {
  await assert.rejects(loadSnapshot({ fetchImpl: async () => new Response('', { status: 404 }) }), /404/);
  await assert.rejects(loadSnapshot({ fetchImpl: async () => new Response('{}') }), /Invalid/);
});
