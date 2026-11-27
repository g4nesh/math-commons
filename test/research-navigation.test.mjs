import test from 'node:test';
import assert from 'node:assert/strict';
import { readFile, readdir, stat } from 'node:fs/promises';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { objectiveGuides, taskGuide } from '../public/research-guides.js';
import { researchPrompt } from '../adapters/research-prompt.mjs';

const root = fileURLToPath(new URL('../', import.meta.url));
const objective = resolve(root, 'research/marcus-de-oliveira');
const index = resolve(objective, 'overview.md');

// Research summaries use inline relative links for artifacts in the same revision.
function artifactLinks(markdown, file) {
  return [...markdown.matchAll(/\]\(([^\s()]+)\)/g)]
    .map(match => match[1].split('#')[0])
    .filter(target => target && !/^\w+:/.test(target))
    .map(target => resolve(dirname(file), decodeURIComponent(target)));
}

test('every research branch is discoverable and its local evidence links resolve', async () => {
  const branches = (await readdir(resolve(objective, 'context')))
    .filter(name => name.endsWith('.md')).map(name => resolve(objective, 'context', name));
  const top = artifactLinks(await readFile(index, 'utf8'), index);
  assert.ok(branches.length > 0);
  for (const file of branches) assert.ok(top.includes(file), `Unlisted research branch: ${file}`);
  for (const file of [index, ...branches]) {
    for (const target of artifactLinks(await readFile(file, 'utf8'), file)) {
      assert.ok((await stat(target)).isFile(), `Missing artifact in ${file}: ${target}`);
    }
  }
});

test('Fourier is the default and each agent enters its assigned scope', async () => {
  const config = JSON.parse(await readFile(resolve(root, 'commons.config.json'), 'utf8'));
  assert.equal(config.defaultObjectiveId, 'fourier-principal-minors');
  const guides = objectiveGuides(config.defaultObjectiveId, config);
  assert.ok((await stat(resolve(root, guides.overview))).isFile());
  assert.ok((await stat(resolve(root, guides.brief))).isFile());
  const start = (objectiveId, id) => researchPrompt({ task: { id, objectiveId }, agent: { id: 'reviewer' } });
  for (const task of ['fourier-154-replay', 'fourier-154-prior-art', 'fourier-154-synthesis', 'fourier-14p-audit']) {
    const page = taskGuide(config.defaultObjectiveId, task);
    assert.ok((await stat(resolve(root, page))).isFile());
    assert.ok(start(config.defaultObjectiveId, task).includes(page));
  }
  const replayPrompt = start(config.defaultObjectiveId, 'fourier-154-replay');
  assert.doesNotMatch(replayPrompt, /family-14p\/README|seed\.md/);
  const other = 'marcus-de-oliveira';
  assert.equal(objectiveGuides(other, config).overview, 'research/marcus-de-oliveira/overview.md');
  assert.ok(start(other, 'mdo-task').includes(objectiveGuides(other, config).overview));
  assert.equal(taskGuide(other, 'fourier-154-replay'), null);
});
