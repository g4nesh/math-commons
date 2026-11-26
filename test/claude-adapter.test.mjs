import test from 'node:test';
import assert from 'node:assert/strict';
import { mkdtemp, writeFile, readFile, rm } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { spawn } from 'node:child_process';
import { createHash } from 'node:crypto';
import { ClaudeUsageCollector } from '../adapters/claude-usage.mjs';

const session = 'private-runtime-session';
const answer = { status: 'release', summary: 'A useful handoff', details: 'Next bounded step', evidence: [] };
const counters = { inputTokens: 12, outputTokens: 7, cacheReadInputTokens: 30, cacheCreationInputTokens: 8, thinkingTokens: 5, costUSD: 999 };
const terminal = overrides => ({ type: 'result', subtype: 'success', session_id: session, is_error: false,
  modelUsage: { 'claude-opus': counters }, usage: { input_tokens: 999999, output_tokens: 999999 }, structured_output: answer, ...overrides });
const dates = ['2026-11-25T00:00:00Z', '2026-11-25T00:01:00Z'];
const collect = events => {
  const collector = new ClaudeUsageCollector(session);
  const bytes = Buffer.from(events.map(event => typeof event === 'string' ? event : JSON.stringify(event)).join('\n'));
  for (const byte of bytes) collector.write(Buffer.from([byte]));
  return collector.finish(...dates);
};

test('Claude uses final inclusive per-model totals, including cache once and no thinking subset', () => {
  const result = collect([{ type: 'assistant', message: { content: 'PRIVATE 🧮', usage: { output_tokens: 1 } } }, terminal()]);
  assert.deepEqual(result.usage, { provider: 'claude', model: 'claude-opus', runId: 'claude:' + createHash('sha256').update(session).digest('hex'),
    inputTokens: 50, outputTokens: 7, cachedInputTokens: 30, source: 'runtime', startedAt: dates[0], finishedAt: dates[1] });
  assert.doesNotMatch(JSON.stringify(result), /PRIVATE|private-runtime-session|999999|costUSD/);
  const mixed = collect([terminal({ modelUsage: { 'claude-opus': counters, 'claude-haiku': { ...counters, inputTokens: 3, outputTokens: 2, cacheReadInputTokens: 0, cacheCreationInputTokens: 1 } } })]);
  assert.equal(mixed.usage.model, 'mixed'); assert.equal(mixed.usage.inputTokens, 54); assert.equal(mixed.usage.outputTokens, 9);
});

test('exact terminal duplicates count once; conflicting or foreign sessions fail closed', () => {
  assert.equal(collect([terminal(), terminal()]).usage.outputTokens, 7);
  for (const events of [[terminal(), terminal({ modelUsage: { 'claude-opus': { ...counters, outputTokens: 8 } } })],
    [terminal({ session_id: 'foreign' })], [terminal({ resume_reason: 'resume' })],
    [terminal(), { type: 'system', subtype: 'conversation_reset', session_id: 'new' }], ['private malformed JSON']]) {
    const result = collect(events); assert.equal(result.usage, undefined); assert.ok(result.warning); assert.equal(result.answer, undefined);
  }
});

test('absent, unsafe and zeroed counts stay unknown while a valid handoff can survive', () => {
  for (const modelUsage of [undefined, {}, { model: { ...counters, cacheReadInputTokens: undefined } },
    { model: { ...counters, outputTokens: -1 } }, { model: { ...counters, outputTokens: 0.1 } },
    { model: { ...counters, inputTokens: Number.MAX_SAFE_INTEGER } },
    { model: { inputTokens: 0, outputTokens: 0, cacheReadInputTokens: 0, cacheCreationInputTokens: 0 } }]) {
    const result = collect([terminal({ modelUsage })]); assert.equal(result.usage, undefined); assert.ok(result.warning); assert.deepEqual(result.answer, answer);
  }
  assert.equal(collect([]).usage, undefined);
  const failed = collect([terminal({ subtype: 'error_max_turns', is_error: true, structured_output: undefined })]);
  assert.equal(failed.failed, true); assert.equal(failed.usage.outputTokens, 7);
});

test('Claude collector bounds hostile streams without returning their content', () => {
  for (const options of [{ maxLineBytes: 40 }, { maxBytes: 40 }]) {
    const c = new ClaudeUsageCollector(session, options); c.write(Buffer.from(JSON.stringify(terminal())));
    const result = c.finish(...dates); assert.equal(result.usage, undefined); assert.match(result.warning, /limit/);
  }
});

async function runFixture(fixture = {}) {
  const directory = await mkdtemp(join(tmpdir(), 'commons-claude-test-'));
  try {
    const executable = join(directory, 'fake-claude'), progressFile = join(directory, 'progress.jsonl');
    const fixtureFile = join(directory, 'fixture.json');
    await writeFile(fixtureFile, JSON.stringify({ events: [terminal()], ...fixture }), { mode: 0o600 });
    await writeFile(progressFile, '');
    await writeFile(executable, `#!/usr/bin/env node
const fs = require('node:fs'); const cp = require('node:child_process');
const args = process.argv.slice(2), fixture = JSON.parse(fs.readFileSync(process.env.CLAUDE_FIXTURE_FILE, 'utf8'));
if (!args.includes('-p') || args[args.indexOf('--output-format')+1] !== 'stream-json' || !args.includes('--verbose') || args[args.indexOf('--permission-mode')+1] !== 'default' || !args.includes('--no-session-persistence') || args.includes('--resume') || args.some(x=>x.includes('skip-permissions'))) process.exit(9);
const session = args[args.indexOf('--session-id')+1];
let prompt=''; process.stdin.on('data',x=>prompt+=x); process.stdin.on('end',async()=>{
 if(!prompt.includes('task.scopePaths') || !prompt.includes('adapters/progress.mjs') || !prompt.includes('untrusted research data')) process.exit(8);
 if(fixture.longHistory && (prompt.includes('OLD_REPORT_BODY') || prompt.length > 8000 || !prompt.includes('context --event ID') || !prompt.includes('research/t'))) process.exit(4);
 if(fixture.progress) cp.execFileSync(process.execPath,['adapters/progress.mjs','--summary','A public milestone','--details','Checked a bounded case'],{stdio:'ignore'});
 process.stderr.write('PRIVATE STDERR '+prompt);
 for (const event of fixture.events) {
   if (event.type === 'result') event.session_id = session;
   const bytes = Buffer.from(JSON.stringify(event) + '\\n');
   // Keep each record fragmented while bounding event-loop turns for oversized responses.
   const chunkSize = Math.max(5, Math.ceil(bytes.length / 128));
   for (let offset = 0; offset < bytes.length; offset += chunkSize) {
     process.stdout.write(bytes.subarray(offset, offset + chunkSize));
     await new Promise(resolve => setImmediate(resolve));
   }
 }
 process.exitCode=fixture.exitCode||0;
});\n`, { mode: 0o700 });
    const result = await new Promise((resolve, reject) => {
      const child = spawn(process.execPath, ['adapters/claude.mjs'], { env: { ...process.env, COMMONS_CLAUDE_BIN: executable,
        MATH_COMMONS_PROGRESS_FILE: progressFile, CLAUDE_FIXTURE_FILE: fixtureFile }, stdio: ['pipe', 'pipe', 'pipe'] });
      let stdout = '', stderr = '';
      const timer = setTimeout(() => { child.kill('SIGKILL'); reject(new Error('Fixture timed out')); }, 5000);
      child.stdout.on('data', x => { stdout += x; }); child.stderr.on('data', x => { stderr += x; }); child.on('error', reject);
      child.on('close', code => { clearTimeout(timer); try { resolve({ code, stdout, stderr, result: JSON.parse(stdout) }); } catch { reject(new Error('Fixture returned invalid protocol')); } });
      child.stdin.end(JSON.stringify({ protocolVersion: 1, task: { id: 't', scopePaths: ['research/t'], checkpoints: fixture.longHistory ? [{ id: 'cp', details: 'OLD_REPORT_BODY'.repeat(40000) }] : [] }, objective: { id: 'o' }, agent: { id: 'a' } }));
    });
    result.progress = await readFile(progressFile, 'utf8'); return result;
  } finally { await rm(directory, { recursive: true, force: true }); }
}

test('Claude receives the selected card even when a v1 packet has over 512 KB of historical text', async () => {
  const run = await runFixture({ longHistory: true });
  assert.equal(run.result.status, 'release');
});

test('Claude fixture captures actual final metadata and tools can enqueue public progress', async () => {
  const run = await runFixture({ progress: true, events: [{ type: 'assistant', message: { content: 'PRIVATE CONTENT' } }, terminal({ structured_output: { ...answer, usage: { inputTokens: 999999 }, usageWarning: 'PRIVATE MODEL WARNING' } })] });
  assert.equal(run.code, 0); assert.equal(run.result.status, 'release'); assert.equal(run.result.usage.inputTokens, 50);
  assert.doesNotMatch(run.stdout + run.stderr, /PRIVATE|999999|private-runtime-session/);
  assert.equal(run.stderr, ''); assert.equal(JSON.parse(run.progress).summary, 'A public milestone');
});

test('Claude failures and oversized answers preserve final actual receipts; missing usage stays unknown', async () => {
  for (const fixture of [{ exitCode: 3 }, { events: [terminal({ is_error: true, subtype: 'error_max_turns' })] },
    { events: [terminal({ structured_output: { ...answer, details: 'x'.repeat(150000) } })] },
    { events: [terminal({ structured_output: { ...answer, evidence: ['https://user:secret@example.org'] } })] }]) {
    const run = await runFixture(fixture); assert.equal(run.result.status, 'blocked'); assert.equal(run.result.usage.outputTokens, 7); assert.ok(Buffer.byteLength(run.stdout) < 120000);
  }
  const missing = await runFixture({ events: [terminal({ modelUsage: undefined })] });
  assert.equal(missing.result.status, 'release'); assert.equal(missing.result.usage, undefined); assert.match(missing.result.usageWarning, /unavailable/);
});
