import { readFile, readdir } from 'node:fs/promises';
import { spawnSync } from 'node:child_process';
const config = JSON.parse(await readFile('commons.config.json', 'utf8'));
if (!/^[\w.-]+\/[\w.-]+$/.test(config.repository)) throw new Error('Invalid repository');
const seed = JSON.parse(await readFile('coordination/state.json', 'utf8'));
if (seed.schemaVersion !== 1 || !seed.objectives.every(item => item.example)) throw new Error('Seed must be labeled as an example');
for (const dir of ['bin', 'lib', 'api', 'scripts', 'public', 'adapters']) {
  for (const file of await readdir(dir)) {
    if (!/\.(mjs|js)$/.test(file)) continue;
    const result = spawnSync(process.execPath, ['--check', `${dir}/${file}`], { stdio: 'inherit' });
    if (result.status !== 0) process.exit(result.status || 1);
  }
}
console.log('Configuration, example labels, and JavaScript syntax passed.');
