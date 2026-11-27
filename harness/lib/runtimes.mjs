import { spawn } from 'node:child_process';
import { randomUUID } from 'node:crypto';
import { existsSync } from 'node:fs';
import { join } from 'node:path';
import { ClaudeUsageCollector } from '../../adapters/claude-usage.mjs';
import { CodexUsageCollector } from '../../adapters/codex-usage.mjs';
import { HarnessError, readText, writeText } from './util.mjs';

// Each run is a fresh, non-persistent session: no resume, no shared history.
// Continuity lives in files (capsule, ledger, artifacts), which is what keeps
// every prompt small and every run reproducible from the record.

function spawnRun(bin, args, { cwd, input, timeoutMs, onStdout, env }) {
  return new Promise(resolve => {
    let child;
    try { child = spawn(bin, args, { cwd, shell: false, stdio: ['pipe', 'pipe', 'pipe'], env: { ...process.env, ...env } }); }
    catch (error) { resolve({ startFailed: true, error: error.message }); return; }
    let stderr = '', timedOut = false;
    const timer = setTimeout(() => { timedOut = true; child.kill('SIGTERM'); setTimeout(() => child.kill('SIGKILL'), 5000).unref(); }, timeoutMs);
    child.stdout.on('data', chunk => onStdout(chunk));
    child.stderr.on('data', chunk => { if (stderr.length < 20_000) stderr += chunk.toString('utf8'); });
    child.stdin.on('error', () => {});
    child.once('error', error => { clearTimeout(timer); resolve({ startFailed: true, error: error.message }); });
    child.once('close', (code, signal) => { clearTimeout(timer); resolve({ code, signal, stderr, timedOut }); });
    child.stdin.end(input);
  });
}

async function runClaude({ settings, prompt, schemaPath, cwd, role, timeoutMs }) {
  const sessionId = randomUUID();
  const collector = new ClaudeUsageCollector(sessionId);
  const args = [
    '-p', '--output-format', 'stream-json', '--verbose',
    '--json-schema', readText(schemaPath),
    '--session-id', sessionId, '--no-session-persistence',
    '--max-turns', String(role.maxTurns || settings.maxTurns || 60),
    '--permission-mode', settings.permissionMode || 'default'
  ];
  if (settings.model) args.push('--model', settings.model);
  const effort = role.effort || settings.effort;
  if (effort) args.push('--effort', effort);
  if (settings.allowedTools?.length) args.push('--allowedTools', settings.allowedTools.join(','));
  if (settings.maxBudgetUsd) args.push('--max-budget-usd', String(settings.maxBudgetUsd));
  const exit = await spawnRun(settings.bin || 'claude', args, { cwd, input: prompt, timeoutMs, onStdout: chunk => collector.write(chunk) });
  const metering = collector.finish(new Date().toISOString(), new Date().toISOString());
  return { exit, answer: metering.answer, usage: metering.usage || null, warning: metering.warning || null, failed: metering.failed };
}

async function runCodex({ settings, prompt, schemaPath, cwd, role, timeoutMs, runDir }) {
  const collector = new CodexUsageCollector();
  const output = join(runDir, 'answer.json');
  const args = ['exec', '--sandbox', settings.sandbox || 'workspace-write', '--color', 'never', '--json', '--ephemeral',
    '--output-schema', schemaPath, '--output-last-message', output];
  if (settings.model) args.push('-m', settings.model);
  if (settings.profile) args.push('-p', settings.profile);
  const effort = role.effort || settings.effort;
  if (effort) args.push('-c', `model_reasoning_effort="${effort}"`);
  args.push(...(settings.extraArgs || []), '-');
  const exit = await spawnRun(settings.bin || 'codex', args, { cwd, input: prompt, timeoutMs, onStdout: chunk => collector.write(chunk) });
  const metering = collector.finish(new Date().toISOString(), new Date().toISOString());
  let answer;
  if (existsSync(output)) {
    try { answer = JSON.parse(readText(output)); } catch { answer = undefined; }
  }
  const usage = metering.usage ? { ...metering.usage, model: settings.model || metering.usage.model } : null;
  return { exit, answer, usage, warning: metering.warning || null, failed: collector.runtimeFailed };
}

/**
 * Run one role. Returns { status: 'ok'|'failed'|'pending', answer, usage, warning, detail }.
 * `manual` writes the prompt for a human or an interactive agent session and
 * returns pending; `discover record` ingests the answer later.
 */
export async function runRuntime({ config, runtime, prompt, schemaPath, role, runDir, cwd }) {
  writeText(join(runDir, 'prompt.md'), prompt);
  if (runtime === 'manual') {
    writeText(join(runDir, 'schema.json'), readText(schemaPath));
    return { status: 'pending', detail: `Prompt written to ${runDir}/prompt.md. Do the step, write the JSON answer to ${runDir}/answer.json, then run: discover record <slug> <run-id>` };
  }
  const settings = config.runtimes[runtime];
  if (!settings) throw new HarnessError('BAD_RUNTIME', `Unknown runtime ${runtime}`);
  const timeoutMs = (settings.timeoutMinutes || config.budgets.stepTimeoutMinutes) * 60_000;
  const runner = runtime === 'claude' ? runClaude : runtime === 'codex' ? runCodex : null;
  if (!runner) throw new HarnessError('BAD_RUNTIME', `No driver for ${runtime}`);
  const startedAt = new Date().toISOString();
  const result = await runner({ settings, prompt, schemaPath, cwd, role, timeoutMs, runDir });
  const finishedAt = new Date().toISOString();
  if (result.usage) Object.assign(result.usage, { startedAt, finishedAt });
  if (result.exit.startFailed) return { status: 'failed', detail: `${runtime} could not start (${result.exit.error}). Install and authenticate it, or set runtimes.${runtime}.bin.`, usage: null };
  if (result.exit.timedOut) return { status: 'failed', detail: `${runtime} timed out after ${timeoutMs / 60_000} minutes`, usage: result.usage, warning: result.warning };
  if (result.exit.code !== 0 || result.failed) return { status: 'failed', detail: `${runtime} exited ${result.exit.code ?? result.exit.signal}${result.exit.stderr ? `: ${result.exit.stderr.trim().split('\n').slice(-3).join(' | ')}` : ''}`, usage: result.usage, warning: result.warning, answer: result.answer };
  if (!result.answer || typeof result.answer !== 'object') return { status: 'failed', detail: `${runtime} returned no structured answer`, usage: result.usage, warning: result.warning };
  return { status: 'ok', answer: result.answer, usage: result.usage, warning: result.warning };
}

export function pickRuntime(config, role, { index = 0, exclude = null, override = null } = {}) {
  if (override) return override;
  const routing = config.routing;
  let choice = routing.roles[role] ?? 'author';
  if (Array.isArray(choice)) choice = choice[index % choice.length];
  if (choice === 'author') choice = routing.author;
  if (choice === 'reviewer') choice = routing.reviewer;
  if (exclude && choice === exclude) {
    const other = [routing.reviewer, routing.author, 'claude', 'codex'].find(name => name !== exclude);
    if (other) choice = other;
  }
  return choice;
}
