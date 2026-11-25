/** Receipt-only contribution metrics. Missing receipts never mean zero usage. */
const DAY = 86_400_000;
export const dayKey = (millis) => new Date(millis).toISOString().slice(0, 10);
export const dayStart = (millis) => Math.floor(millis / DAY) * DAY;
const count = (value) => Number.isSafeInteger(value) && value >= 0;
const safeAdd = (left, right) => {
  const value = left + right;
  if (!Number.isSafeInteger(value)) throw new RangeError('Token totals exceed safe integer precision.');
  return value;
};
function timestamp(value) {
  if (typeof value !== 'string' || !/^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d{1,3})?(?:Z|[+-]\d{2}:\d{2})$/.test(value)) return NaN;
  const [year, month, day] = value.slice(0, 10).split('-').map(Number);
  const check = new Date(Date.UTC(year, month - 1, day));
  if (check.getUTCFullYear() !== year || check.getUTCMonth() !== month - 1 || check.getUTCDate() !== day) return NaN;
  const [hour, minute, second] = value.slice(11, 19).split(':').map(Number);
  if (hour > 23 || minute > 59 || second > 59) return NaN;
  const offset = value.match(/[+-](\d{2}):(\d{2})$/);
  if (offset && (Number(offset[1]) > 23 || Number(offset[2]) > 59)) return NaN;
  return Date.parse(value);
}
export function normalizeReceipts(receipts, now = Date.now()) {
  const candidates = [], ids = new Map(), runs = new Map(), conflicted = new Set();
  let invalidCount = 0, duplicateCount = 0;
  for (const receipt of Array.isArray(receipts) ? receipts : []) {
    const finished = timestamp(receipt?.finishedAt), started = timestamp(receipt?.startedAt), recorded = timestamp(receipt?.recordedAt);
    if (!receipt || typeof receipt.id !== 'string' || !receipt.id || typeof receipt.owner !== 'string' || !/^[a-zA-Z0-9][a-zA-Z0-9-]{0,38}$/.test(receipt.owner)
        || typeof receipt.runId !== 'string' || !receipt.runId || typeof receipt.provider !== 'string' || !receipt.provider
        || !['runtime', 'manual'].includes(receipt.source) || !count(receipt.inputTokens) || !count(receipt.outputTokens) || !count(receipt.cachedInputTokens)
        || receipt.cachedInputTokens > receipt.inputTokens || !count(receipt.inputTokens + receipt.outputTokens)
        || !Number.isFinite(finished) || !Number.isFinite(started) || !Number.isFinite(recorded) || started > finished || finished > recorded + 300_000 || recorded > now + 300_000) {
      invalidCount += 1; continue;
    }
    const owner = receipt.owner.toLowerCase();
    const runKey = JSON.stringify([owner, receipt.provider.toLowerCase(), receipt.runId]);
    const fingerprint = JSON.stringify([runKey, receipt.agentId, receipt.taskId, receipt.objectiveId, receipt.model, receipt.inputTokens, receipt.outputTokens, receipt.cachedInputTokens, started, finished, receipt.source]);
    const previous = [...new Set([ids.get(receipt.id), runs.get(runKey)].filter((value) => value !== undefined))];
    if (previous.length) {
      if (previous.every((index) => !conflicted.has(index) && candidates[index].fingerprint === fingerprint)) duplicateCount++;
      else {
        invalidCount++;
        for (const index of previous) if (!conflicted.has(index)) { conflicted.add(index); invalidCount++; }
      }
      continue;
    }
    ids.set(receipt.id, candidates.length); runs.set(runKey, candidates.length);
    candidates.push({ ...receipt, owner, totalTokens: receipt.inputTokens + receipt.outputTokens, day: dayKey(finished), finished, recorded, fingerprint });
  }
  const valid = candidates.filter((_, index) => !conflicted.has(index));
  valid.sort((a, b) => a.finished - b.finished || a.id.localeCompare(b.id));
  const pendingCount = valid.filter((receipt) => receipt.finished > now).length;
  return { receipts: valid.filter((receipt) => receipt.finished <= now), pendingCount, invalidCount, duplicateCount };
}
export function aggregateUsage(receipts, { period = '30d', owners = null, now = Date.now() } = {}) {
  if (!['7d', '30d', 'all'].includes(period)) throw new TypeError('Unsupported date window.');
  if (!Number.isFinite(now)) throw new TypeError('Invalid current time.');
  const normalized = normalizeReceipts(receipts, now);
  const end = dayStart(now);
  const first = normalized.receipts.length ? dayStart(normalized.receipts[0].finished) : end - 29 * DAY;
  const start = period === 'all' ? first : end - (period === '7d' ? 6 : 29) * DAY;
  if (end - start > 366 * 100 * DAY) throw new RangeError('Receipt history exceeds the supported 100-year window.');
  const selected = owners === null ? null : new Set(owners.map((owner) => owner.toLowerCase()));
  const availableOwners = [...new Set(normalized.receipts.map((receipt) => receipt.owner))].sort();
  const eligible = normalized.receipts.filter((receipt) => receipt.finished >= start && (selected === null || selected.has(receipt.owner)));
  const dates = [];
  for (let date = start; date <= end; date += DAY) dates.push(dayKey(date));
  const totals = { total: null, input: null, output: null, cachedInput: null, receiptCount: 0, runtimeCount: 0, manualCount: 0 };
  function add(target, receipt) {
    target.total = safeAdd(target.total || 0, receipt.totalTokens);
    target.input = safeAdd(target.input || 0, receipt.inputTokens);
    target.output = safeAdd(target.output || 0, receipt.outputTokens);
    target.cachedInput = safeAdd(target.cachedInput || 0, receipt.cachedInputTokens);
    target.receiptCount += 1;
    target[receipt.source === 'runtime' ? 'runtimeCount' : 'manualCount'] += 1;
  }
  const empty = () => ({ total: null, input: null, output: null, cachedInput: null, receiptCount: 0, runtimeCount: 0, manualCount: 0 });
  const days = dates.map((date) => ({ date, ...empty(), cumulative: 0 }));
  const dayMap = new Map(days.map((day) => [day.date, day]));
  const seriesMap = new Map(availableOwners.filter((owner) => selected === null || selected.has(owner)).map((owner) => [owner, { owner, ...empty(), days: dates.map((date) => ({ date, ...empty(), cumulative: 0 })) }]));
  const ownerDays = new Map([...seriesMap].map(([owner, series]) => [owner, new Map(series.days.map((day) => [day.date, day]))]));
  for (const receipt of eligible) {
    add(totals, receipt); add(dayMap.get(receipt.day), receipt);
    add(seriesMap.get(receipt.owner), receipt); add(ownerDays.get(receipt.owner).get(receipt.day), receipt);
  }
  for (const collection of [days, ...[...seriesMap.values()].map((series) => series.days)]) {
    let running = 0;
    for (const day of collection) { running = safeAdd(running, day.total || 0); day.cumulative = running; }
  }
  return { ...totals, period, startDay: dayKey(start), endDay: dayKey(end), days, series: [...seriesMap.values()], receipts: eligible,
    availableOwners, pendingCount: normalized.pendingCount, allReceiptCount: normalized.receipts.length, invalidCount: normalized.invalidCount, duplicateCount: normalized.duplicateCount,
    lastRecordedAt: eligible.length ? new Date(Math.max(...eligible.map((receipt) => receipt.recorded))).toISOString() : null };
}
export function ownerActivity(state) {
  const owners = new Map();
  const get = (owner) => {
    const key = String(owner || '').toLowerCase();
    if (!/^[a-zA-Z0-9][a-zA-Z0-9-]{0,38}$/.test(key)) return null;
    if (!owners.has(key)) owners.set(key, { owner: key, agents: 0, checkpoints: 0, messages: 0, reviews: 0, events: 0 });
    return owners.get(key);
  };
  for (const agent of state?.agents || []) { const row = get(agent.owner); if (row) row.agents++; }
  for (const event of state?.events || []) {
    const row = get(event.actorOwner); if (!row) continue;
    row.events++;
    if (event.kind === 'checkpoint') row.checkpoints++;
    if (event.kind === 'message') row.messages++;
    if (event.kind === 'review') row.reviews++;
  }
  return [...owners.values()].sort((a, b) => b.checkpoints - a.checkpoints || a.owner.localeCompare(b.owner));
}
/** Synthetic receipts only. The caller must explicitly enter labeled example mode. */
export function exampleUsage(now = Date.now()) {
  const today = dayStart(now), receipts = [];
  const profiles = [{ owner: 'alex-example', base: 9200, phase: 0 }, { owner: 'sam-example', base: 7000, phase: 2 }, { owner: 'river-example', base: 5300, phase: 4 }];
  for (let index = 0; index < 56; index++) for (const [profileIndex, profile] of profiles.entries()) {
    if ((index + profile.phase) % 7 === 0 || (index * 3 + profile.phase) % 11 === 0 || index < profile.phase * 2) continue;
    const finished = today - (55 - index) * DAY + 2 * 3600_000;
    if (finished > now) continue;
    const inputTokens = Math.round(profile.base * (1 + index / 45) * (0.65 + ((index * 13 + profileIndex * 7) % 17) / 15));
    receipts.push({ id: `synthetic-${profile.owner}-${index}`, agentId: `${profile.owner}-agent`, owner: profile.owner, taskId: 'example-proof', objectiveId: 'example-sum-of-cubes', provider: 'example-provider', model: 'synthetic-example-model', runId: `example-${profile.owner}-${index}`, inputTokens, outputTokens: Math.round(inputTokens * .2), cachedInputTokens: Math.round(inputTokens * .35), startedAt: new Date(finished - 600_000).toISOString(), finishedAt: new Date(finished).toISOString(), recordedAt: new Date(finished).toISOString(), source: index % 9 === 0 ? 'manual' : 'runtime' });
  }
  return receipts;
}
