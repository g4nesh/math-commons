#!/usr/bin/env node
// Stands in for `claude -p` and `codex exec` in tests. It reads the harness
// prompt, writes the fixture's files for that role, and speaks each runtime's
// real output protocol so the production drivers and usage collectors are
// exercised unchanged.
import { appendFileSync, existsSync, mkdirSync, readFileSync, readdirSync, writeFileSync } from 'node:fs';
import { dirname, join, resolve } from 'node:path';
import { roleOutput } from './euler.mjs';

const argv = process.argv.slice(2);
const codex = argv[0] === 'exec';
const input = readFileSync(0, 'utf8');
const field = name => (new RegExp(`^${name}: (.*)$`, 'm').exec(input) || [])[1]?.trim();
const role = field('ROLE'), runId = field('RUN_ID');
const dir = resolve(process.cwd(), field('CAMPAIGN_DIR'));
const mode = process.env.FAKE_RUNTIME_MODE || '';
const direction = (/^## Your direction (D\d+)/m.exec(input) || [])[1];
const attempt = direction ? String((existsSync(join(dir, 'ledger.jsonl')) ? readFileSync(join(dir, 'ledger.jsonl'), 'utf8').split('\n').filter(line => line.includes(`"direction":"${direction}"`)).length : 0) + 1).padStart(3, '0') : null;
const certificate = (/^## Certificate ([a-z0-9-]+)/m.exec(input) || [])[1];
const reviewsSoFar = existsSync(join(dir, 'reviews')) ? readdirSync(join(dir, 'reviews')).filter(name => /^review-\d+\.json$/.test(name)).length : 0;
if (process.env.FAKE_RUNTIME_LOG) appendFileSync(process.env.FAKE_RUNTIME_LOG, JSON.stringify({ runtime: codex ? 'codex' : 'claude', role, runId, promptBytes: input.length }) + '\n');

let output;
if (mode === 'fail') process.exit(3);
if (mode === 'logged-out' && !codex) {
  const session = argv[argv.indexOf('--session-id') + 1];
  process.stdout.write(JSON.stringify({ type: 'result', subtype: 'success', is_error: true, session_id: session, result: 'Not logged in · Please run /login', modelUsage: {} }) + '\n');
  process.exit(1);
}
output = roleOutput(role, { runId, direction, attempt, certificate, mode, reviewsSoFar });
for (const [path, text] of Object.entries(output.files || {})) {
  mkdirSync(dirname(join(dir, path)), { recursive: true });
  writeFileSync(join(dir, path), text);
}
for (const [path, text] of Object.entries(output.append || {})) appendFileSync(join(dir, path), text);
if (mode === 'scope-violation' && role === 'map') writeFileSync(join(dir, 'campaign-notes.md'), 'out of scope\n');

const usage = { input: 1200 + input.length / 4 | 0, output: 300 };
if (codex) {
  const out = argv[argv.indexOf('--output-last-message') + 1];
  writeFileSync(out, JSON.stringify(output.answer));
  const events = [
    { type: 'thread.started', thread_id: `fake-${runId}` },
    { type: 'turn.started' },
    { type: 'turn.completed', usage: { input_tokens: usage.input, output_tokens: usage.output, cached_input_tokens: 0 } }
  ];
  process.stdout.write(events.map(event => JSON.stringify(event)).join('\n') + '\n');
} else {
  const session = argv[argv.indexOf('--session-id') + 1];
  const events = [
    { type: 'system', subtype: 'init', session_id: session },
    { type: 'result', subtype: 'success', is_error: false, session_id: session, structured_output: output.answer,
      modelUsage: { 'claude-fake': { inputTokens: usage.input, outputTokens: usage.output, cacheReadInputTokens: 0, cacheCreationInputTokens: 0 } } }
  ];
  process.stdout.write(events.map(event => JSON.stringify(event)).join('\n') + '\n');
}
