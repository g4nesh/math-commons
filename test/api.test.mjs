import test from 'node:test';
import assert from 'node:assert/strict';
import { loadSnapshot } from '../api/state.js';
test('snapshot returns GitHub provenance and preserves revision', async () => {
  const state = { schemaVersion: 1, revision: 9, objectives: [], tasks: [], agents: [], events: [] };
  const snapshot = await loadSnapshot({ fetchImpl: async (url, options) => {
    assert.equal(url, 'https://api.github.com/repos/g4nesh/math-commons/contents/coordination/state.json?ref=coordination');
    assert.equal(options.headers['Cache-Control'], 'no-cache');
    return new Response(JSON.stringify({type:'file',encoding:'base64',sha:'current-state-blob',size:100,content:Buffer.from(JSON.stringify(state)).toString('base64')}));
  }, now: 0 });
  assert.deepEqual(snapshot.state, state);
  assert.equal(snapshot.source, 'github');
  assert.equal(snapshot.contentSha, 'current-state-blob');
  assert.equal(snapshot.fetchedAt, '1970-01-01T00:00:00.000Z');
});
test('GitHub failures never masquerade as example or empty state', async () => {
  await assert.rejects(loadSnapshot({ fetchImpl: async () => new Response('', { status: 404 }) }), /404/);
  await assert.rejects(loadSnapshot({ fetchImpl: async () => new Response('{}') }), /Invalid/);
  await assert.rejects(loadSnapshot({ fetchImpl: async () => new Response('{}', {status:403}) }), /403/);
});
