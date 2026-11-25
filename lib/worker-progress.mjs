import { createHash } from 'node:crypto';
import { mkdir, mkdtemp, open, writeFile, unlink, rm } from 'node:fs/promises';
import { resolve, join } from 'node:path';
import { StringDecoder } from 'node:string_decoder';

export function validateMilestone(value) {
  const allowed = new Set(['version', 'type', 'id', 'summary', 'details', 'evidence']);
  if (!value || Array.isArray(value) || Object.keys(value).some(key => !allowed.has(key)) ||
      value.version !== 1 || value.type !== 'checkpoint' || !/^[a-zA-Z0-9._-]{1,100}$/.test(value.id ?? '') ||
      typeof value.summary !== 'string' || !value.summary.trim() || value.summary.length > 500 ||
      typeof value.details !== 'string' || !value.details.trim() || value.details.length > 4000 ||
      !Array.isArray(value.evidence) || value.evidence.length > 10 || !value.evidence.every(url => {
        try { const parsed = new URL(url); return typeof url === 'string' && url.length <= 2000 && parsed.protocol === 'https:' && !parsed.username && !parsed.password; } catch { return false; }
      }) || Buffer.byteLength(JSON.stringify(value)) > 32_000) {
    throw new Error('Invalid public milestone. Use only the documented bounded checkpoint fields.');
  }
  return { version: 1, type: 'checkpoint', id: value.id, summary: value.summary, details: value.details, evidence: [...value.evidence] };
}

/** Reads only an explicit publication artifact; never watches runtime logs. */
export async function createWorkerProgress({ repository, agent, task, runId, mutate, onProgress = () => {}, onFailure = () => {},
  directory = resolve('.commons/progress'), throttleMs = 30_000, pollMs = 1000, now = Date.now } = {}) {
  await mkdir(directory, { recursive: true, mode: 0o700 });
  const runDirectory = await mkdtemp(join(directory, 'run-'));
  const path = join(runDirectory, 'milestones.jsonl');
  const recoveryFile = join(runDirectory, 'pending.json');
  const file = await open(path, 'wx+', 0o600);
  let buffer = '', bytes = 0, lastWrite = -Infinity, pending = null, failure = null, writing = null, reading = null, closed = false, recoverySaved = false;
  const decoder = new StringDecoder('utf8'), seen = new Map();
  const fail = () => {
    if (failure) return;
    failure = new Error('Public milestone publication stopped. Inspect the local progress files before resuming.');
    failure.code = 'PROGRESS_FAILED'; failure.progressFile = path;
    if (recoverySaved) failure.recoveryFile = recoveryFile;
    onProgress({ kind: 'progress-warning', taskId: task, summary: failure.message, progressFile: path, ...(recoverySaved ? { recoveryFile } : {}) });
    onFailure(failure);
  };
  const flush = async (force = false) => {
    if (writing) { await writing; if (force && pending && !failure) return flush(true); return; }
    if (!pending || failure || (!force && now() - lastWrite < throttleMs)) return;
    const frame = pending; pending = null;
    const operation = { id: 'progress-' + createHash('sha256').update(`${repository}\0${agent}\0${task}\0${runId}\0${frame.id}`).digest('hex'),
      kind: 'checkpoint', agent, task, summary: frame.summary, details: frame.details, evidence: frame.evidence };
    writing = (async () => {
      // Preserve the exact operation before any network write. A lost response
      // must be reconciled with this identity, never replaced with a new one.
      await writeFile(recoveryFile, JSON.stringify({ schemaVersion: 1, repository, operation }, null, 2) + '\n', { mode: 0o600 });
      recoverySaved = true;
      await mutate(operation);
      await unlink(recoveryFile);
      recoverySaved = false;
      lastWrite = now();
      onProgress({ kind: 'checkpoint', taskId: task, requestId: operation.id, summary: frame.summary });
    })().catch(fail).finally(() => { writing = null; });
    await writing;
  };
  const accept = value => {
    const frame = validateMilestone(value), content = JSON.stringify(frame);
    if (seen.has(frame.id)) {
      if (seen.get(frame.id) !== content) throw new Error('A milestone ID was reused with different content.');
      return;
    }
    if (seen.size >= 200) throw new Error('Too many milestones in one worker run.');
    seen.set(frame.id, content);
    // Bursts retain their newest deliberate milestone; the complete local
    // artifact remains available if the run cannot finish publishing.
    pending = frame;
    void flush();
  };
  const line = text => {
    if (!text.trim()) return;
    if (Buffer.byteLength(text) > 32_000) throw new Error('Milestone exceeds frame limit.');
    accept(JSON.parse(text));
  };
  const read = async () => {
    if (closed || failure) return;
    if (reading) return reading;
    reading = (async () => {
      const chunk = Buffer.alloc(16_384);
      while (!failure) {
        const { bytesRead } = await file.read(chunk, 0, chunk.length, null);
        if (!bytesRead) break;
        bytes += bytesRead;
        if (bytes > 1_048_576) throw new Error('Milestone artifact exceeds one megabyte.');
        buffer += decoder.write(chunk.subarray(0, bytesRead));
        let newline;
        while ((newline = buffer.indexOf('\n')) !== -1) {
          const current = buffer.slice(0, newline); buffer = buffer.slice(newline + 1); line(current);
        }
        if (Buffer.byteLength(buffer) > 32_000) throw new Error('Milestone exceeds frame limit.');
      }
      await flush();
    })().catch(fail).finally(() => { reading = null; });
    return reading;
  };
  const timer = setInterval(() => { void read(); }, pollMs); timer.unref?.();
  return {
    path, recoveryFile,
    async finish() {
      clearInterval(timer);
      // A prior poll can be waiting for a slow network write after its EOF.
      // The adapter may append its last frame during that wait. Always read
      // again after the previous poll settles before deleting the artifact.
      if (reading) await reading;
      await read();
      try { buffer += decoder.end(); if (buffer.trim()) line(buffer); } catch { fail(); }
      await flush(true);
      closed = true; await file.close();
      if (failure) throw failure;
      await rm(runDirectory, { recursive: true, force: true });
    }
  };
}
