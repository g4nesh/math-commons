#!/usr/bin/env node
import { readWorkPacket, researchPrompt } from './research-prompt.mjs';
import { spawn } from 'node:child_process';
import { readFile } from 'node:fs/promises';
import { randomUUID } from 'node:crypto';
import { fileURLToPath } from 'node:url';
import { ClaudeUsageCollector } from './claude-usage.mjs';

const packet = await readWorkPacket(process.stdin);
const schema = await readFile(new URL('./result.schema.json', import.meta.url), 'utf8');
const sessionId = randomUUID(), collector = new ClaudeUsageCollector(sessionId);
const checkout = fileURLToPath(new URL('../', import.meta.url));
const prompt = researchPrompt(packet);
const blocked = summary => ({ status: 'blocked', summary, details: 'Inspect the local Claude Code setup and resume deliberately. Runtime transcripts were not published.', evidence: [] });
const startedAt = new Date().toISOString();
const exit = await new Promise(resolve => {
  // No resume, permission bypass, tool auto-approval, or session log mining.
  const child = spawn(process.env.COMMONS_CLAUDE_BIN || 'claude', [
    '-p', '--output-format', 'stream-json', '--verbose', '--json-schema', schema,
    '--session-id', sessionId, '--no-session-persistence', '--max-turns', '30', '--permission-mode', 'default'
  ], { cwd: checkout, shell: false, stdio: ['pipe', 'pipe', 'ignore'] });
  child.stdout.on('data', chunk => collector.write(chunk));
  child.stdin.on('error', () => {});
  child.once('error', () => resolve({ startFailed: true }));
  child.once('close', code => resolve({ code }));
  child.stdin.end(prompt);
});
const metering = collector.finish(startedAt, new Date().toISOString());
let result;
if (exit.startFailed) result = blocked('Claude Code could not start. Install it and authenticate with your own account.');
else if (exit.code !== 0 || metering.failed) result = blocked('Claude Code reported a failed or interrupted run.');
else {
  try {
    const answer = metering.answer;
    if (!answer || !['review', 'release', 'blocked'].includes(answer.status) ||
        typeof answer.summary !== 'string' || !answer.summary.trim() || answer.summary.length > 2000 ||
        typeof answer.details !== 'string' || answer.details.length > 24000 ||
        !Array.isArray(answer.evidence) || answer.evidence.length > 30 ||
        !answer.evidence.every(value => { try { const url = new URL(value); return typeof value === 'string' && value.length <= 2000 && url.protocol === 'https:' && !url.username && !url.password; } catch { return false; } }) ||
        (answer.status === 'review' && (!answer.details.trim() || !answer.evidence.length))) throw new Error('Invalid result.');
    result = { status: answer.status, summary: answer.summary, details: answer.details, evidence: answer.evidence };
  } catch { result = blocked('Claude Code did not produce a valid structured worker result.'); }
}
if (metering.usage) result.usage = metering.usage;
else result.usageWarning = `Usage unavailable: ${metering.warning}. No token receipt was recorded.`;
if (Buffer.byteLength(JSON.stringify(result)) > 120_000) result = { ...blocked('Claude Code final result exceeded the worker protocol limit.'), ...(metering.usage ? { usage: metering.usage } : { usageWarning: result.usageWarning }) };
process.stdout.write(JSON.stringify(result) + '\n');
