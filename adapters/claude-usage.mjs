import { createHash } from 'node:crypto';
import { StringDecoder } from 'node:string_decoder';

// Fresh, single-input print invocation only. result.modelUsage includes the
// main agent and its subagents; result.usage does not. Never sum both, resumed
// totals, assistant placeholders, thinking subsets, or cache TTL breakdowns.
// https://code.claude.com/docs/en/agent-sdk/cost-tracking#get-the-total-cost-of-a-query
// https://code.claude.com/docs/en/agent-sdk/typescript#modelusage
export class ClaudeUsageCollector {
  constructor(sessionId, { maxBytes = 33_554_432, maxLineBytes = 1_048_576 } = {}) {
    this.sessionId = sessionId; this.maxBytes = maxBytes; this.maxLineBytes = maxLineBytes;
    this.bytes = 0; this.buffer = ''; this.decoder = new StringDecoder('utf8');
    this.warning = null; this.terminal = null; this.fingerprint = null;
  }
  unavailable(message, discardAnswer = true) { this.warning ||= message; this.buffer = ''; this.terminal = null; if (discardAnswer) this.answer = undefined; }
  write(chunk) {
    if (this.warning) return;
    this.bytes += chunk.length;
    if (this.bytes > this.maxBytes) return this.unavailable('runtime stream exceeded the capture limit');
    this.buffer += this.decoder.write(chunk);
    let newline;
    while (!this.warning && (newline = this.buffer.indexOf('\n')) !== -1) {
      const line = this.buffer.slice(0, newline); this.buffer = this.buffer.slice(newline + 1); this.line(line);
    }
    if (Buffer.byteLength(this.buffer) > this.maxLineBytes) this.unavailable('runtime event exceeded the capture limit');
  }
  line(line) {
    if (!line.trim()) return;
    if (Buffer.byteLength(line) > this.maxLineBytes) return this.unavailable('runtime event exceeded the capture limit');
    let event;
    try { event = JSON.parse(line); } catch { return this.unavailable('runtime stream was not valid JSONL'); }
    if (!event || Array.isArray(event) || typeof event.type !== 'string') return this.unavailable('runtime event shape was unsupported');
    if (event.type === 'system' && ['init', 'conversation_reset'].includes(event.subtype) && event.session_id !== this.sessionId) return this.unavailable('runtime session changed or did not match this fresh invocation');
    if (event.type !== 'result') return; // Never retain conversation/tool events.
    if (event.session_id !== this.sessionId || event.resume_reason || event.local_command || !['success', 'error_max_turns', 'error_during_execution', 'error_max_budget_usd', 'error_max_structured_output_retries'].includes(event.subtype) || typeof event.is_error !== 'boolean') return this.unavailable('terminal runtime identity or result was unsupported');
    this.answer = event.structured_output;
    this.runtimeFailed = event.is_error || event.subtype !== 'success';
    const models = event.modelUsage;
    if (!models || typeof models !== 'object' || Array.isArray(models) || !Object.keys(models).length || Object.keys(models).length > 50) return this.unavailable('final per-model token counters were missing', false);
    let inputTokens = 0, outputTokens = 0, cachedInputTokens = 0;
    for (const [name, counters] of Object.entries(models)) {
      if (!name.trim() || name.length > 160 || !/^[a-zA-Z0-9._:/-]+$/.test(name)) return this.unavailable('runtime model identity was unsupported', false);
      const fields = ['inputTokens', 'outputTokens', 'cacheReadInputTokens', 'cacheCreationInputTokens'];
      if (!fields.every(key => Number.isSafeInteger(counters?.[key]) && counters[key] >= 0)) return this.unavailable('final token counters were missing or unsupported', false);
      inputTokens += counters.inputTokens + counters.cacheReadInputTokens + counters.cacheCreationInputTokens;
      outputTokens += counters.outputTokens; cachedInputTokens += counters.cacheReadInputTokens;
      if (![inputTokens, outputTokens, cachedInputTokens, inputTokens + outputTokens].every(Number.isSafeInteger)) return this.unavailable('final token counters exceeded the safe range', false);
    }
    // A crash can emit zeroed result counters. Do not certify that as zero use.
    if (inputTokens + outputTokens === 0) return this.unavailable('all-zero counters may be a runtime failure fallback', false);
    const names = Object.keys(models).sort();
    const receipt = { provider: 'claude', model: names.length === 1 ? names[0] : 'mixed',
      runId: 'claude:' + createHash('sha256').update(this.sessionId).digest('hex'),
      inputTokens, outputTokens, cachedInputTokens, source: 'runtime' };
    // A single fresh prompt has one terminal result. Exact repeat delivery is
    // harmless; differing results cannot safely be treated as additive calls.
    const terminal = { receipt, failed: event.is_error || event.subtype !== 'success', answer: event.structured_output };
    const fingerprint = JSON.stringify({ ...terminal, models: Object.fromEntries(names.map(name => [name, ['inputTokens', 'outputTokens', 'cacheReadInputTokens', 'cacheCreationInputTokens'].map(key => models[name][key])])) });
    if (this.fingerprint && this.fingerprint !== fingerprint) return this.unavailable('multiple terminal results were ambiguous');
    this.fingerprint = fingerprint; this.terminal = terminal;
  }
  finish(startedAt, finishedAt) {
    if (!this.warning) { this.buffer += this.decoder.end(); if (this.buffer.trim()) this.line(this.buffer); this.buffer = ''; }
    if (this.warning || !this.terminal) return { warning: this.warning || 'no complete terminal usage record was returned', answer: this.answer, failed: this.runtimeFailed };
    return { usage: { ...this.terminal.receipt, startedAt, finishedAt }, answer: this.terminal.answer, failed: this.terminal.failed };
  }
}
