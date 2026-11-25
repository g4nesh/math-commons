import test from 'node:test';
import assert from 'node:assert/strict';
import { readFile, readdir, stat } from 'node:fs/promises';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

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
