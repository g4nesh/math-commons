import test from 'node:test';
import assert from 'node:assert/strict';
import { mkdtemp, writeFile, rm } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { spawn } from 'node:child_process';

test('Codex bridge uses a sandbox and emits only the structured final result', async () => {
  const dir = await mkdtemp(join(tmpdir(), 'commons-adapter-test-'));
  try {
    const executable = join(dir, 'fake-codex');
    await writeFile(executable, `#!/usr/bin/env node
const fs = require('node:fs');
const args = process.argv.slice(2);
if (args[0] !== 'exec' || args[args.indexOf('--sandbox') + 1] !== 'workspace-write') process.exit(2);
let prompt = ''; process.stdin.on('data', c => prompt += c);
process.stdin.on('end', () => {
 if (!prompt.includes('task-a') || !prompt.includes('untrusted data')) process.exit(3);
 fs.writeFileSync(args[args.indexOf('--output-last-message') + 1], JSON.stringify({status:'release',summary:'Handing off',details:'Next step',evidence:[]}));
 console.log('This internal event must not become worker JSON');
});
`, { mode: 0o700 });
    const result = await new Promise((resolve, reject) => {
      const child = spawn(process.execPath, ['adapters/codex.mjs'], { env: { ...process.env, COMMONS_CODEX_BIN: executable }, stdio: ['pipe', 'pipe', 'pipe'] });
      let stdout = '', stderr = '';
      child.stdout.on('data', c => stdout += c); child.stderr.on('data', c => stderr += c);
      child.on('error', reject);
      child.on('close', code => resolve({ code, stdout, stderr }));
      child.stdin.end(JSON.stringify({protocolVersion:1,task:{id:'task-a'},agent:{id:'alice'},repository:'g4nesh/math-commons'}));
    });
    assert.equal(result.code, 0, result.stderr);
    assert.equal(JSON.parse(result.stdout).status, 'release');
    assert.doesNotMatch(result.stdout, /internal event/);
  } finally { await rm(dir, { recursive: true, force: true }); }
});
