import { createHash } from 'node:crypto';
import { StringDecoder } from 'node:string_decoder';

// Exec emits thread totals, then shuts down after turn.completed. Fresh exec
// has no earlier turns; never use this collector for resume/fork.
// Primary source (lines 120–130, 513–540):
// https://github.com/openai/codex/blob/89feadf3cca2a53d01bbf4c35ce53227704b2417/codex-rs/exec/src/event_processor_with_jsonl_output.rs
// Checked 2026-11-25; local CLI 0.169.1. Model identity is absent from this API.
export class CodexUsageCollector {
  constructor({ maxLineBytes = 1_048_576, maxBytes = 33_554_432 } = {}) {
    this.maxLineBytes = maxLineBytes;
    this.maxBytes = maxBytes;
    this.bytes = 0;
    this.buffer = '';
    this.decoder = new StringDecoder('utf8');
    this.phase = 'thread';
    this.runId = null;
    this.counts = null;
    this.warning = null;
    this.runtimeFailed = false;
  }

  unavailable(reason) {
    this.warning ||= reason;
    this.buffer = '';
  }

  write(chunk) {
    if (this.warning) return;
    this.bytes += chunk.length;
    if (this.bytes > this.maxBytes) {
      this.unavailable('runtime event stream exceeded the capture limit');
      return;
    }
    this.buffer += this.decoder.write(chunk);
    let newline;
    while (!this.warning && (newline = this.buffer.indexOf('\n')) !== -1) {
      const line = this.buffer.slice(0, newline);
      this.buffer = this.buffer.slice(newline + 1);
      this.line(line);
    }
    if (Buffer.byteLength(this.buffer) > this.maxLineBytes) {
      this.unavailable('runtime event exceeded the capture limit');
    }
  }

  line(line) {
    if (!line.trim()) return;
    if (Buffer.byteLength(line) > this.maxLineBytes) {
      this.unavailable('runtime event exceeded the capture limit');
      return;
    }
    let event;
    try { event = JSON.parse(line); } catch {
      this.unavailable('runtime event stream was not valid JSONL');
      return;
    }
    if (!event || Array.isArray(event) || typeof event.type !== 'string') {
      this.unavailable('runtime event shape was unsupported');
      return;
    }
    // Do not retain or forward transcript, command, reasoning, or tool events.
    if (event.type === 'thread.started') {
      if (this.phase !== 'thread' || typeof event.thread_id !== 'string' ||
          !event.thread_id.trim() || event.thread_id.length > 1024) {
        this.unavailable('runtime thread identity was missing or ambiguous');
        return;
      }
      this.runId = 'codex:' + createHash('sha256').update(event.thread_id).digest('hex');
      this.phase = 'turn';
    } else if (event.type === 'turn.started') {
      if (this.phase !== 'turn') {
        this.unavailable('multiple or out-of-order turns cannot be metered safely');
        return;
      }
      this.phase = 'running';
    } else if (event.type === 'turn.completed') {
      const u = event.usage;
      const values = [u?.input_tokens, u?.output_tokens, u?.cached_input_tokens];
      if (!values.every(value => Number.isSafeInteger(value) && value >= 0) ||
          u.cached_input_tokens > u.input_tokens ||
          !Number.isSafeInteger(u.input_tokens + u.output_tokens)) {
        this.unavailable('runtime token counters were missing or unsupported');
        return;
      }
      const counts = { inputTokens: u.input_tokens, outputTokens: u.output_tokens,
        cachedInputTokens: u.cached_input_tokens };
      if (this.phase === 'completed' &&
          JSON.stringify(counts) === JSON.stringify(this.counts)) return;
      if (this.phase !== 'running') {
        this.unavailable('duplicate or out-of-order completion counters were ambiguous');
        return;
      }
      // The official emitter uses Usage::default() when usage is unavailable.
      // All-zero cannot distinguish that fallback from measured zero.
      if (u.input_tokens + u.output_tokens === 0) {
        this.unavailable('all-zero counters may represent unavailable runtime usage');
        return;
      }
      this.counts = counts;
      this.phase = 'completed';
    } else if (event.type === 'turn.failed' || event.type === 'error') {
      this.runtimeFailed = true;
    }
  }

  finish(startedAt, finishedAt) {
    if (!this.warning) {
      this.buffer += this.decoder.end();
      if (this.buffer) this.line(this.buffer);
      this.buffer = '';
    }
    if (!this.warning && (!this.runId || !this.counts || this.phase !== 'completed')) {
      this.unavailable('no complete runtime usage receipt was emitted');
    }
    if (this.warning) return { warning: this.warning };
    return { usage: { provider: 'codex', model: 'unspecified', runId: this.runId,
      ...this.counts, startedAt, finishedAt, source: 'runtime' } };
  }
}
