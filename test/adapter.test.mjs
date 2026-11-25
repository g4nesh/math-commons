import test from 'node:test';
import assert from 'node:assert/strict';
import { mkdtemp, writeFile, readFile, rm } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { spawn } from 'node:child_process';
import { createHash } from 'node:crypto';
import { CodexUsageCollector } from '../adapters/codex-usage.mjs';

const answer = { status: 'release', summary: 'Handing off', details: 'Next step', evidence: [] };
const counts = { input_tokens: 123, output_tokens: 17, cached_input_tokens: 80 };
const started = [{ type: 'thread.started', thread_id: 'runtime-thread-secret-id' }, { type: 'turn.started' }];
const completed = { type: 'turn.completed', usage: counts };
const jsonl = events => events.map(event => typeof event === 'string' ? event : JSON.stringify(event)).join('\n');
const dates = ['2026-11-25T06:00:00.000Z', '2026-11-25T06:01:00.000Z'];

async function runAdapter(fixture) {
  const dir = await mkdtemp(join(tmpdir(), 'commons-adapter-test-'));
  try {
    const executable = join(dir, 'fake-codex');
    await writeFile(executable, [
      '#!/usr/bin/env node',
      "const fs = require('node:fs');",
      'const args = process.argv.slice(2);',
      "if (args[0] !== 'exec' || args[args.indexOf('--sandbox')+1] !== 'workspace-write' || !args.includes('--json') || args.includes('resume')) process.exit(2);",
      "const fixture = JSON.parse(process.env.COMMONS_ADAPTER_FIXTURE);",
      "let prompt = ''; process.stdin.on('data', chunk => prompt += chunk);",
      "process.stdin.on('end', async () => {",
      " if (!prompt.includes('task-a') || !prompt.includes('untrusted data')) process.exit(3);",
      " if (fixture.progressFile) require('node:child_process').execFileSync(process.execPath, ['adapters/progress.mjs', '--summary', 'Codex public milestone', '--details', 'A bounded check completed'], { stdio: 'ignore' });",
      " if (fixture.oversizedAnswer) fixture.answer.details = 'x'.repeat(150000);",
      " if (fixture.escapedAnswer) fixture.answer.evidence = Array.from({length:30}, () => 'https://example.test/' + '\\u0000'.repeat(1900));",
      " const final = fixture.largeAnswer ? 'x'.repeat(512001) : fixture.invalidAnswer ? 'invalid secret answer' : JSON.stringify(fixture.answer);",
      " fs.writeFileSync(args[args.indexOf('--output-last-message')+1], final);",
      " process.stderr.write('SECRET child diagnostic: ' + prompt);",
      " const bytes = Buffer.from(fixture.stream || '');",
      " for (let i=0; i<bytes.length; i+=7) { process.stdout.write(bytes.subarray(i,i+7)); await new Promise(resolve=>setImmediate(resolve)); }",
      ' process.exitCode = fixture.exitCode || 0;',
      '});',
      ''
    ].join('\n'), { mode: 0o700 });
    return await new Promise((resolve, reject) => {
      const child = spawn(process.execPath, ['adapters/codex.mjs'], {
        env: { ...process.env, COMMONS_CODEX_BIN: executable,
          MATH_COMMONS_PROGRESS_FILE: fixture.progressFile || '',
          COMMONS_ADAPTER_FIXTURE: JSON.stringify({ answer, ...fixture }) },
        stdio: ['pipe', 'pipe', 'pipe']
      });
      let stdout = '', stderr = '';
      const timeout = setTimeout(() => { child.kill('SIGKILL'); reject(new Error('Fake adapter timed out')); }, 5000);
      child.stdout.on('data', chunk => { stdout += chunk; });
      child.stderr.on('data', chunk => { stderr += chunk; });
      child.on('error', reject);
      child.on('close', code => {
        clearTimeout(timeout);
        try { resolve({ code, stdout, stderr, result: JSON.parse(stdout) }); } catch { reject(new Error('Adapter omitted its result')); }
      });
      child.stdin.end(JSON.stringify({ protocolVersion: 1, task: { id: 'task-a' }, agent: { id: 'alice' }, repository: 'g4nesh/math-commons' }));
    });
  } finally { await rm(dir, { recursive: true, force: true }); }
}

test('Codex adapter captures runtime counts, hashes identity, and publishes no events or diagnostics', async () => {
  const run = await runAdapter({
    stream: jsonl([...started, { type: 'item.completed', item: { text: 'SECRET model transcript 🧮', type: 'agent_message' } }, completed]),
    answer: { ...answer, usage: { inputTokens: 999999 }, usageWarning: 'SECRET forged warning' }
  });
  assert.equal(run.code, 0);
  assert.equal(run.result.status, 'release');
  assert.deepEqual(Object.keys(run.result).sort(), ['details', 'evidence', 'status', 'summary', 'usage']);
  assert.deepEqual({ ...run.result.usage, startedAt: undefined, finishedAt: undefined }, {
    provider: 'codex', model: 'unspecified',
    runId: 'codex:' + createHash('sha256').update(started[0].thread_id).digest('hex'),
    inputTokens: 123, outputTokens: 17, cachedInputTokens: 80, source: 'runtime',
    startedAt: undefined, finishedAt: undefined
  });
  assert.ok(Date.parse(run.result.usage.startedAt) <= Date.parse(run.result.usage.finishedAt));
  assert.doesNotMatch(run.stdout + run.stderr, /SECRET|runtime-thread-secret-id|999999/);
  assert.equal(run.stderr, '');
});

test('Codex tool subprocess can enqueue an explicit public milestone without a special descriptor', async () => {
  const directory = await mkdtemp(join(tmpdir(), 'codex-progress-test-'));
  try {
    const progressFile = join(directory, 'milestones.jsonl'); await writeFile(progressFile, '');
    const run = await runAdapter({ progressFile, stream: jsonl([...started, completed]) });
    assert.equal(run.result.usage.inputTokens, 123);
    const frame = JSON.parse(await readFile(progressFile, 'utf8'));
    assert.equal(frame.type, 'checkpoint'); assert.equal(frame.summary, 'Codex public milestone');
    assert.doesNotMatch(run.stdout + run.stderr, /SECRET/);
  } finally { await rm(directory, { recursive: true, force: true }); }
});

test('missing events and all-zero fallback produce visible missing status, never invented zeros', async () => {
  for (const stream of ['', jsonl([...started, { ...completed, usage: { input_tokens: 0, output_tokens: 0, cached_input_tokens: 0 } }])]) {
    const run = await runAdapter({ stream, answer: { ...answer, usage: { inputTokens: 5 } } });
    assert.equal(run.result.status, 'release');
    assert.equal(run.result.usage, undefined);
    assert.match(run.result.usageWarning, /Usage unavailable/);
    assert.match(run.stderr, /No token receipt/);
    assert.doesNotMatch(run.stdout + run.stderr, /SECRET/);
  }
});

test('duplicate completion is counted exactly once; conflicting or multiple turns are not guessed', async () => {
  const duplicate = await runAdapter({ stream: jsonl([...started, completed, completed]) });
  assert.equal(duplicate.result.usage.inputTokens, 123);
  for (const suffix of [
    [{ ...completed, usage: { ...counts, input_tokens: 124 } }],
    [{ type: 'turn.started' }, completed]
  ]) {
    const run = await runAdapter({ stream: jsonl([...started, completed, ...suffix]) });
    assert.equal(run.result.usage, undefined);
    assert.match(run.result.usageWarning, /ambiguous|multiple/);
  }
});

test('malformed JSONL is not exposed or silently trusted', async () => {
  const run = await runAdapter({ stream: jsonl([...started, '{"SECRET": invalid', completed]) });
  assert.equal(run.result.usage, undefined);
  assert.match(run.result.usageWarning, /JSONL/);
  assert.doesNotMatch(run.stdout + run.stderr, /SECRET/);
});

test('failed child and invalid final answer preserve real usage with blocked status', async () => {
  for (const fixture of [{ exitCode: 4 }, { invalidAnswer: true }, { largeAnswer: true },
    { oversizedAnswer: true }, { escapedAnswer: true }]) {
    const run = await runAdapter({ ...fixture, stream: jsonl([...started, completed]) });
    assert.equal(run.code, 0);
    assert.equal(run.result.status, 'blocked');
    assert.equal(run.result.usage.outputTokens, 17);
    assert.ok(Buffer.byteLength(run.stdout) < 120_000);
    assert.doesNotMatch(run.stdout + run.stderr, /SECRET|invalid secret answer/);
  }
});

test('final-answer field limits preserve captured usage instead of overflowing the worker', async () => {
  for (const changed of [
    { summary: 'x'.repeat(2001) },
    { details: 'x'.repeat(24001) },
    { evidence: Array(31).fill('https://example.test/proof') },
    { evidence: ['https://example.test/' + 'x'.repeat(2000)] },
    { evidence: ['file:///private/secret'] },
    { status: 'review', details: '', evidence: [] }
  ]) {
    const run = await runAdapter({ answer: { ...answer, ...changed }, stream: jsonl([...started, completed]) });
    assert.equal(run.result.status, 'blocked');
    assert.equal(run.result.usage.inputTokens, 123);
    assert.ok(Buffer.byteLength(run.stdout) < 120_000);
  }
});

test('runtime failure without completion stays unmetered; failed event cannot become a success', async () => {
  const run = await runAdapter({ stream: jsonl([...started, { type: 'turn.failed', error: { message: 'SECRET provider error' } }]) });
  assert.equal(run.result.status, 'blocked');
  assert.equal(run.result.usage, undefined);
  assert.match(run.result.usageWarning, /no complete/);
  assert.doesNotMatch(run.stdout + run.stderr, /SECRET/);
});

test('counter validation rejects missing, fractional, negative, unsafe, and inconsistent values', () => {
  const invalid = [
    {}, { ...counts, cached_input_tokens: undefined }, { ...counts, input_tokens: -1 },
    { ...counts, output_tokens: 1.5 }, { ...counts, input_tokens: '123' },
    { ...counts, cached_input_tokens: 124 },
    { ...counts, input_tokens: Number.MAX_SAFE_INTEGER, output_tokens: 1 },
    { ...counts, input_tokens: Number.MAX_SAFE_INTEGER + 1 }
  ];
  for (const usage of invalid) {
    const collector = new CodexUsageCollector();
    collector.write(Buffer.from(jsonl([...started, { ...completed, usage }])));
    const result = collector.finish(...dates);
    assert.equal(result.usage, undefined);
    assert.match(result.warning, /counters/);
  }
});

test('fragmented UTF-8, CRLF, no final newline, and actual zero cache are supported', () => {
  const collector = new CodexUsageCollector();
  const bytes = Buffer.from(jsonl([...started, { type: 'item.updated', item: { text: '🧮 ignored' } },
    { ...completed, usage: { ...counts, cached_input_tokens: 0 } }]).replaceAll('\n', '\r\n'));
  for (const byte of bytes) collector.write(Buffer.from([byte]));
  const result = collector.finish(...dates);
  assert.equal(result.usage.cachedInputTokens, 0);
  assert.equal(result.usage.inputTokens + result.usage.outputTokens, 140);
  assert.equal(result.warning, undefined);
});

test('bounded capture omits receipts on overlong events and total stream limits', () => {
  for (const options of [{ maxLineBytes: 128 }, { maxBytes: 128 }]) {
    const collector = new CodexUsageCollector(options);
    collector.write(Buffer.from(jsonl([...started, { type: 'item.updated', item: { text: 'x'.repeat(256) } }, completed])));
    const result = collector.finish(...dates);
    assert.equal(result.usage, undefined);
    assert.match(result.warning, /capture limit/);
  }
});

test('missing or ambiguous runtime identity and out-of-order lifecycle cannot create a receipt', () => {
  const variants = [
    [{ type: 'turn.started' }, completed],
    [{ type: 'thread.started', thread_id: '' }, { type: 'turn.started' }, completed],
    [...started, started[0], completed],
    [started[0], completed],
    ['[]'],
    [...started, completed, '{"type":']
  ];
  for (const events of variants) {
    const collector = new CodexUsageCollector();
    collector.write(Buffer.from(jsonl(events)));
    assert.equal(collector.finish(...dates).usage, undefined);
  }
});
