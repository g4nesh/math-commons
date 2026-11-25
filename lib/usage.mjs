import { createHash } from 'node:crypto';
import { mkdir, readFile, writeFile, unlink } from 'node:fs/promises';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';

export const USAGE_CLOCK_SKEW_MS = 5 * 60_000;

export class UsageError extends Error {
  constructor(message, code = 'INVALID_USAGE') { super(message); this.name = 'UsageError'; this.code = code; }
}
const reject = message => { throw new UsageError(message); };
const boundedText = (value, label, max) => {
  if (typeof value !== 'string' || !value.trim() || value.trim().length > max || /[\x00-\x1f\x7f]/.test(value)) reject(`${label} must be nonempty text, at most ${max} characters, without control characters.`);
  return value.trim();
};
export function tokenCount(value, label) {
  if (!Number.isSafeInteger(value) || value < 0) reject(`${label} must be a nonnegative safe integer.`);
  return value;
}
export function parseTokenFlag(value, label) {
  if (typeof value !== 'string' || !/^\d+$/.test(value)) reject(`${label} must contain only decimal digits.`);
  return tokenCount(Number(value), label);
}
function timestamp(value, label) {
  if (typeof value !== 'string' || value.length > 40 || !/^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d{1,3})?(?:Z|[+-]\d{2}:\d{2})$/.test(value)) reject(`${label} must be an ISO timestamp with a timezone.`);
  const parsed = Date.parse(value);
  const [year, month, day] = value.slice(0, 10).split('-').map(Number);
  const calendar = new Date(`${value.slice(0, 10)}T00:00:00Z`);
  if (!Number.isFinite(parsed) || calendar.getUTCFullYear() !== year || calendar.getUTCMonth() + 1 !== month || calendar.getUTCDate() !== day || Number(value.slice(11, 13)) > 23 || Number(value.slice(14, 16)) > 59 || Number(value.slice(17, 19)) > 59) reject(`${label} is not a valid timestamp.`);
  return new Date(parsed).toISOString();
}

/** Only explicitly reported counts are accepted. Cached input is NOT additive. */
export function normalizeUsage(value, { defaultSource = 'manual' } = {}) {
  if (!value || typeof value !== 'object' || Array.isArray(value)) reject('Usage must be an explicit receipt object.');
  const provider = boundedText(value.provider, 'Provider', 80).toLowerCase();
  const model = boundedText(value.model, 'Model', 160);
  const runId = boundedText(value.runId, 'Run ID', 200);
  if (!/^[A-Za-z0-9][A-Za-z0-9:._/-]*$/.test(runId)) reject('Run ID must use letters, digits, colon, dot, underscore, slash, or hyphen.');
  const source = value.source === undefined ? defaultSource : value.source;
  if (!['runtime', 'manual'].includes(source)) reject('Usage source must be runtime or manual.');
  if (value.cachedInputTokens === undefined) reject('Usage requires an explicit cachedInputTokens count; missing capture is not zero.');
  const inputTokens = tokenCount(value.inputTokens, 'Input tokens');
  const outputTokens = tokenCount(value.outputTokens, 'Output tokens');
  const cachedInputTokens = tokenCount(value.cachedInputTokens, 'Cached input tokens');
  if (cachedInputTokens > inputTokens) reject('Cached input tokens are a subset of input tokens and cannot exceed them.');
  const totalTokens = tokenCount(inputTokens + outputTokens, 'Total tokens');
  if (value.totalTokens !== undefined && value.totalTokens !== totalTokens) reject('Total tokens must equal input plus output; cached input is not additive.');
  const startedAt = timestamp(value.startedAt, 'startedAt'), finishedAt = timestamp(value.finishedAt, 'finishedAt');
  if (Date.parse(startedAt) > Date.parse(finishedAt)) reject('finishedAt must not precede startedAt.');
  if (value.evidence !== undefined && (!Array.isArray(value.evidence) || value.evidence.length > 20)) reject('Usage evidence must contain at most 20 HTTPS URLs.');
  const evidence = [...new Set((value.evidence ?? []).map(item => {
    const text = boundedText(item, 'Evidence URL', 2_000);
    let url;
    try { url = new URL(text); } catch { reject('Usage evidence must be an HTTPS URL.'); }
    if (url.protocol !== 'https:' || url.username || url.password) reject('Usage evidence must be an HTTPS URL without credentials.');
    return url.href;
  }))];
  return { provider, model, runId, inputTokens, outputTokens, cachedInputTokens, totalTokens, startedAt, finishedAt, source, ...(evidence.length ? { evidence } : {}) };
}
export function usageRunKey(owner, receipt) {
  return JSON.stringify([String(owner).toLowerCase(), receipt.provider.toLowerCase(), receipt.runId]);
}
export function sameUsage(left, right) {
  return JSON.stringify(normalizeUsage(left)) === JSON.stringify(normalizeUsage(right));
}
export function usageOperationId(owner, agent, receipt) {
  return `usage-${createHash('sha256').update(JSON.stringify([usageRunKey(owner, receipt), agent])).digest('hex')}`;
}
export const defaultRecoveryDirectory = fileURLToPath(new URL('../.commons/usage/', import.meta.url));

/** Persist only the public receipt/recovery context, never a token or raw runtime log. */
export async function saveUsageRecovery(envelope, directory = defaultRecoveryDirectory) {
  await mkdir(directory, { recursive: true, mode: 0o700 });
  const payload = JSON.stringify(envelope, null, 2) + '\n';
  const digest = createHash('sha256').update(payload).digest('hex').slice(0, 16);
  const path = join(directory, `${envelope.requestId}-${digest}.json`);
  try { await writeFile(path, payload, { flag: 'wx', mode: 0o600 }); }
  catch (error) {
    if (error.code !== 'EEXIST' || await readFile(path, 'utf8') !== payload) throw error;
  }
  return path;
}
export async function clearUsageRecovery(path) {
  if (!path) return;
  try { await unlink(path); } catch (error) { if (error.code !== 'ENOENT') throw error; }
}

/** A failed upload never erases the only captured receipt or pretends zero usage. */
export async function publishWorkerUsage(store, { agent, owner, task, receipt, recoveryDirectory, onProgress = () => {} }) {
  if (receipt === undefined || receipt === null) {
    const result = { status: 'missing', summary: 'Token usage was not captured for this run; no zero or estimated count was recorded.' };
    onProgress({ kind: 'usage-warning', taskId: task, ...result }); return result;
  }
  let usage;
  try { usage = normalizeUsage(receipt, { defaultSource: 'runtime' }); }
  catch (error) {
    const result = { status: 'invalid', summary: `Token usage was rejected: ${error.message}. No count was recorded.` };
    onProgress({ kind: 'usage-warning', taskId: task, ...result }); return result;
  }
  const requestId = usageOperationId(owner, agent, usage);
  const envelope = { schemaVersion: 1, repository: store.config.repository, agent, task, requestId, receipt: usage };
  let recoveryPath = null;
  try { recoveryPath = await saveUsageRecovery(envelope, recoveryDirectory); }
  catch { onProgress({ kind: 'usage-warning', taskId: task, summary: 'Could not save the local usage recovery file. An upload will still be attempted.' }); }
  try {
    const result = await store.mutate({ id: requestId, kind: 'usage', agent, task, ...usage });
    try { await clearUsageRecovery(recoveryPath); }
    catch { onProgress({ kind: 'usage-warning', taskId: task, summary: 'Usage was recorded, but its local recovery file could not be removed. Replaying it is deduplicated.', recoveryPath }); }
    const outcome = { status: result.duplicate ? 'already-recorded' : 'recorded', requestId, receiptId: result.usage?.id ?? result.event?.id, totalTokens: usage.totalTokens };
    onProgress({ kind: 'usage', taskId: task, ...outcome }); return outcome;
  } catch (error) {
    const outcome = { status: 'pending', requestId, summary: 'Token usage upload was not confirmed. Preserve and replay this receipt; do not rerun the model to recover counts.', recoveryPath, recovery: envelope, errorCode: error.code ?? 'UPLOAD_FAILED' };
    onProgress({ kind: 'usage-warning', taskId: task, ...outcome }); return outcome;
  }
}
