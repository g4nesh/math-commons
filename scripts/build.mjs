import { cp, mkdir, readFile, writeFile, rm } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';
import { execFileSync } from 'node:child_process';
const root = fileURLToPath(new URL('../', import.meta.url));
process.chdir(root);
await rm('dist', { recursive: true, force: true });
await mkdir('dist', { recursive: true });
await cp('public', 'dist', { recursive: true });
await cp('commons.config.json', 'dist/config.json');
await cp('coordination/state.json', 'dist/seed-state.json');
await cp('AGENTS.md', 'dist/AGENTS.md');
await cp('docs', 'dist/docs', { recursive: true });
await cp('README.md', 'dist/README.md');
let commit = process.env.GITHUB_SHA || process.env.VERCEL_GIT_COMMIT_SHA || 'local';
if (commit === 'local') { try { commit = execFileSync('git', ['rev-parse', 'HEAD'], { encoding: 'utf8' }).trim(); } catch {} }
await writeFile('dist/deployment.json', JSON.stringify({ commit, builtAt: new Date().toISOString() }) + '\n');
JSON.parse(await readFile('dist/config.json', 'utf8'));
console.log(`Built Math Commons (${commit.slice(0, 12)})`);
