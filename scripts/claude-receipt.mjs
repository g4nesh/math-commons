#!/usr/bin/env node
import { createHash } from 'node:crypto';
import { open } from 'node:fs/promises';
import { resolve } from 'node:path';
import { StringDecoder } from 'node:string_decoder';
import { fileURLToPath } from 'node:url';
import { normalizeUsage } from '../lib/usage.mjs';

// Explicit-file exporter, not discovery or an uploader. Stop the producing
// runtime before exporting. Timestamps bound the observed provider records.
// Claude input = uncached input + cache creation + cache read, exactly once:
// https://platform.claude.com/docs/en/build-with-claude/prompt-caching
// Shared response IDs must be deduplicated; SDK streaming output placeholders
// are not completed provider records:
// https://code.claude.com/docs/en/agent-sdk/cost-tracking
const TERMINAL = new Set(['end_turn', 'tool_use', 'max_tokens', 'stop_sequence',
  'pause_turn', 'refusal', 'model_context_window_exceeded']);
class ReceiptError extends Error {}
const fail = message => { throw new ReceiptError(message); };
const safeCount = value => Number.isSafeInteger(value) && value >= 0;
const identity = value => typeof value === 'string' && value.length > 0 &&
  value.length <= 300 && !/[\x00-\x20\x7f]/.test(value);
function sum(...values) {
  let result = 0;
  for (const value of values) {
    if (!safeCount(value) || !safeCount(result + value)) fail('Missing, invalid, or overflowing provider token counters.');
    result += value;
  }
  return result;
}
function timestamp(value) {
  if (typeof value !== 'string' || !/^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d{1,3})?(?:Z|[+-]\d{2}:\d{2})$/.test(value)) {
    fail('Every provider response needs an explicit timestamp with timezone.');
  }
  const date = new Date(value), calendar = new Date(value.slice(0, 10) + 'T00:00:00Z');
  if (!Number.isFinite(date.getTime()) || calendar.toISOString().slice(0, 10) !== value.slice(0, 10) ||
      Number(value.slice(11, 13)) > 23 || Number(value.slice(14, 16)) > 59 || Number(value.slice(17, 19)) > 59) {
    fail('A provider response timestamp is invalid.');
  }
  return date.getTime();
}

export class ClaudeReceiptCollector {
  constructor() {
    this.records = new Map();
    this.requests = new Map();
    this.messages = new Map();
    this.session = null;
    this.first = null;
    this.started = Infinity;
    this.finished = -Infinity;
  }

  add(record) {
    if (!record || typeof record !== 'object' || Array.isArray(record)) fail('A JSONL record has an unsupported shape.');
    if (record.type !== 'assistant') return;
    const message = record.message;
    if (message?.model === '<synthetic>') return;
    if (!message || message.type !== 'message' || message.role !== 'assistant' ||
        !identity(record.sessionId) || !identity(record.requestId) || !identity(message.id) ||
        !identity(message.model) || message.model.length > 160) {
      fail('Assistant response identity or provider metadata is missing or ambiguous.');
    }
    if (!TERMINAL.has(message.stop_reason)) fail('An assistant response has no supported terminal stop reason.');
    if (this.session !== null && this.session !== record.sessionId) fail('Multiple sessions occur in the selected file.');
    this.session = record.sessionId;
    const at = timestamp(record.timestamp);
    const usage = message.usage;
    if (!usage || typeof usage !== 'object' || Array.isArray(usage)) fail('A provider response has no complete usage counters.');
    const input = usage.input_tokens, writes = usage.cache_creation_input_tokens;
    const reads = usage.cache_read_input_tokens, output = usage.output_tokens;
    const inputTokens = sum(input, writes, reads);
    sum(output, inputTokens);
    // Nested iterations and cache_creation TTL buckets are breakdowns, not
    // additional requests. Never recurse through them or through message content.
    const key = JSON.stringify([record.requestId, message.id]);
    if ((this.requests.has(record.requestId) && this.requests.get(record.requestId) !== message.id) ||
        (this.messages.has(message.id) && this.messages.get(message.id) !== record.requestId)) {
      fail('A request or message identity is reused inconsistently.');
    }
    const fingerprint = JSON.stringify([message.model, message.stop_reason, input, writes, reads, output]);
    if (this.records.has(key) && this.records.get(key).fingerprint !== fingerprint) {
      fail('Duplicate provider responses contain inconsistent counters or metadata.');
    }
    if (!this.records.has(key)) {
      if (this.records.size >= 100_000) fail('The selected file exceeds the supported response limit.');
      this.records.set(key, { fingerprint, model: message.model, inputTokens, outputTokens: output, cachedInputTokens: reads });
      this.requests.set(record.requestId, message.id);
      this.messages.set(message.id, record.requestId);
    }
    if (!this.first || at < this.first.at || (at === this.first.at && key < this.first.key)) this.first = { at, key };
    this.started = Math.min(this.started, at);
    this.finished = Math.max(this.finished, at);
  }

  finish() {
    if (!this.records.size) fail('No completed provider responses with usable counters were found.');
    let inputTokens = 0, outputTokens = 0, cachedInputTokens = 0;
    const models = new Set();
    for (const record of this.records.values()) {
      inputTokens = sum(inputTokens, record.inputTokens);
      outputTokens = sum(outputTokens, record.outputTokens);
      cachedInputTokens = sum(cachedInputTokens, record.cachedInputTokens);
      models.add(record.model);
    }
    const digest = createHash('sha256').update(JSON.stringify(['claude-file-receipt-v1', this.session, this.first.key])).digest('hex');
    return normalizeUsage({
      provider: 'claude', model: models.size === 1 ? [...models][0] : 'mixed',
      runId: 'claude:' + digest, inputTokens, outputTokens, cachedInputTokens,
      startedAt: new Date(this.started).toISOString(), finishedAt: new Date(this.finished).toISOString(),
      source: 'runtime'
    });
  }
}

export async function exportClaudeReceipt(path, { maxFileBytes = 268_435_456, maxLineBytes = 8_388_608 } = {}) {
  if (typeof path !== 'string' || !path) fail('Select exactly one finished JSONL file.');
  let file;
  try { file = await open(path, 'r'); } catch { fail('Could not open the explicitly selected file.'); }
  try {
    const before = await file.stat();
    if (!before.isFile() || before.size > maxFileBytes) fail('Select a regular JSONL file within the supported size limit.');
    const collector = new ClaudeReceiptCollector(), decoder = new StringDecoder('utf8');
    const chunk = Buffer.alloc(65_536);
    let buffer = '', total = 0, lines = 0;
    const addLine = line => {
      if (!line.trim()) return;
      if (Buffer.byteLength(line) > maxLineBytes) fail('A JSONL record exceeds the supported size limit.');
      let record;
      try { record = JSON.parse(line); } catch { fail('The selected file contains invalid or incomplete JSONL.'); }
      collector.add(record);
      lines++;
      if (lines > 1_000_000) fail('The selected file exceeds the supported record limit.');
    };
    while (true) {
      const { bytesRead } = await file.read(chunk, 0, chunk.length, null);
      if (!bytesRead) break;
      total += bytesRead;
      if (total > maxFileBytes) fail('The selected file exceeds the supported size limit.');
      buffer += decoder.write(chunk.subarray(0, bytesRead));
      let newline;
      while ((newline = buffer.indexOf('\n')) !== -1) {
        const line = buffer.slice(0, newline);
        buffer = buffer.slice(newline + 1);
        addLine(line);
      }
      if (Buffer.byteLength(buffer) > maxLineBytes) fail('A JSONL record exceeds the supported size limit.');
    }
    buffer += decoder.end();
    if (buffer) addLine(buffer);
    const after = await file.stat();
    if (before.size !== after.size || before.mtimeMs !== after.mtimeMs || before.ctimeMs !== after.ctimeMs) {
      fail('The selected file changed while reading. Stop the producing runtime before exporting.');
    }
    return collector.finish();
  } finally { await file.close(); }
}

const HELP = 'Usage: node scripts/claude-receipt.mjs /absolute/path/to/finished-session.jsonl\n' +
  'Reads one explicitly selected file; no discovery, network, uploads, or file writes.\n' +
  'Exports observed completed-response counters, not billing totals. Select child files separately.\n' +
  'Stop the producer first. Do not import overlapping copies, partial exports, or a parent aggregate alongside children.\n' +
  'Do not use SDK streaming placeholders as final output counts. Raw identifiers and message content are never emitted.\n';

if (process.argv[1] && resolve(process.argv[1]) === fileURLToPath(import.meta.url)) {
  const args = process.argv.slice(2);
  if (args.length === 1 && ['--help', '-h'].includes(args[0])) process.stdout.write(HELP);
  else if (args.length !== 1 || args[0].startsWith('-')) {
    process.stderr.write('Select exactly one finished JSONL file. Use --help for the supported format.\n');
    process.exitCode = 1;
  } else {
    try { process.stdout.write(JSON.stringify(await exportClaudeReceipt(args[0]), null, 2) + '\n'); }
    catch (error) {
      process.stderr.write((error instanceof ReceiptError ? error.message : 'Receipt export failed; no receipt was emitted.') + '\n');
      process.exitCode = 1;
    }
  }
}
