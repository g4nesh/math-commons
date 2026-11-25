#!/usr/bin/env node
import { readWorkPacket, researchPrompt } from './research-prompt.mjs';
import { spawn } from 'node:child_process';
import { mkdtemp, open, rm } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { CodexUsageCollector } from './codex-usage.mjs';

// Invoked explicitly by `commons work`; does not start a model at install time.
const packet = await readWorkPacket(process.stdin);
const directory = await mkdtemp(join(tmpdir(), 'math-commons-codex-'));
const output = join(directory, 'result.json');
const schema = fileURLToPath(new URL('./result.schema.json', import.meta.url));
const checkout = fileURLToPath(new URL('../', import.meta.url));
const prompt = researchPrompt(packet);
const collector = new CodexUsageCollector();
const startedAt = new Date().toISOString();
const blocked = summary => ({ status: 'blocked', summary, details:
  'The local Codex run did not return a usable final result. Inspect your local setup and resume the task deliberately. Runtime transcripts were not published.',
  evidence: [] });
let result;
try {
  const exit = await new Promise(resolve => {
    const child = spawn(process.env.COMMONS_CODEX_BIN || 'codex', [
      'exec', '--sandbox', 'workspace-write', '--color', 'never', '--json',
      '--output-schema', schema, '--output-last-message', output, '-'
    ], { cwd: checkout, shell: false, stdio: ['pipe', 'pipe', 'ignore'] });
    child.stdout.on('data', chunk => collector.write(chunk));
    child.once('error', () => resolve({ startFailed: true }));
    child.stdin.on('error', () => {});
    child.once('close', (code, signal) => resolve({ code, signal }));
    child.stdin.end(prompt);
  });
  const finishedAt = new Date().toISOString();
  const metering = collector.finish(startedAt, finishedAt);
  if (exit.startFailed) {
    result = blocked('Codex CLI could not start. Install it and sign in with your own account.');
  } else if (exit.code !== 0 || collector.runtimeFailed) {
    result = blocked('Codex CLI reported a failed or interrupted run.');
  } else {
    try {
      // Read a bounded final answer, independently of the runtime event stream.
      const file = await open(output, 'r');
      let text;
      try {
        const buffer = Buffer.alloc(512_001);
        const { bytesRead } = await file.read(buffer, 0, buffer.length, 0);
        if (bytesRead > 512_000) throw new Error('Final answer exceeds limit.');
        text = buffer.subarray(0, bytesRead).toString('utf8');
      } finally { await file.close(); }
      const answer = JSON.parse(text);
      if (!['review', 'release', 'blocked'].includes(answer.status) ||
          typeof answer.summary !== 'string' || !answer.summary.trim() || answer.summary.length > 2000 ||
          typeof answer.details !== 'string' || answer.details.length > 24000 ||
          !Array.isArray(answer.evidence) || answer.evidence.length > 30 ||
          !answer.evidence.every(value => typeof value === 'string' && value.length <= 2000 &&
            new URL(value).protocol === 'https:')) {
        throw new Error('Invalid final answer.');
      }
      if (answer.status === 'review' && (!answer.details.trim() || !answer.evidence.length)) {
        throw new Error('Review requires a substantive handoff with evidence.');
      }
      // Whitelist answer fields. A model-authored "usage" is never trusted.
      result = { status: answer.status, summary: answer.summary,
        details: answer.details, evidence: answer.evidence };
    } catch {
      result = blocked('Codex CLI did not produce a valid worker result.');
    }
  }
  if (metering.usage) result.usage = metering.usage;
  else {
    result.usageWarning = 'Usage unavailable: ' + metering.warning + '. No token receipt was recorded.';
    process.stderr.write(result.usageWarning + '\n');
  }
  // The outer worker accepts at most 128 KB. Escaped JSON can be much larger
  // than character counts, so preserve captured usage in a compact failure.
  if (Buffer.byteLength(JSON.stringify(result)) > 120_000) {
    result = { ...blocked('Codex CLI final result exceeded the worker protocol limit.'),
      ...(metering.usage ? { usage: metering.usage } : { usageWarning: result.usageWarning }) };
  }
  process.stdout.write(JSON.stringify(result) + '\n');
} finally {
  await rm(directory, { recursive: true, force: true });
}
