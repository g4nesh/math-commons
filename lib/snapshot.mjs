import { createHash } from 'node:crypto';

// Canonical reads/writes retain the existing capacity. The legacy deployment's
// smaller limit is a measurement for rollout checks, not a research-write cap.
export const MAX_STATE_BYTES = 3_000_000;
export const LEGACY_READER_MAX_BYTES = 950_000;

const detailHash = value => createHash('sha256').update(value, 'utf8').digest('hex');
const detailPointer = event => `Full detail: task ${event.taskId}, checkpoint ${event.id}.`;
const requestPointer = event => `Full request: task ${event.taskId}, description.`;
const checkpointSummary = 'Checkpoint recorded; see task checkpoint.';
const referenceKeys = ['detailReference', 'checkpointMetadataReference', 'taskDescriptionReference'];
const validReference = reference => reference && typeof reference === 'object' && !Array.isArray(reference)
  && Object.keys(reference).sort().join(',') === 'sha256,version' && reference.version === 1
  && typeof reference.sha256 === 'string' && /^[0-9a-f]{64}$/.test(reference.sha256);
const smaller = (candidate, original) => Buffer.byteLength(JSON.stringify(candidate), 'utf8') < Buffer.byteLength(JSON.stringify(original), 'utf8');
const requestMatches = (event, task) => event.kind === 'request' && task?.request
  && task.id === event.taskId && task.createdAt === event.createdAt && task.createdBy === event.actor;

function checkpointIndex(state) {
  const result = new Map();
  for (const task of state?.tasks ?? []) {
    if (result.has(task.id)) throw new Error('Duplicate task in coordination snapshot');
    const checkpoints = new Map();
    for (const checkpoint of task.checkpoints ?? []) {
      if (checkpoints.has(checkpoint.id)) throw new Error('Duplicate checkpoint in coordination snapshot');
      checkpoints.set(checkpoint.id, checkpoint);
    }
    result.set(task.id, checkpoints);
  }
  return result;
}

function checkExpandedCapacity(state) {
  if (Buffer.byteLength(JSON.stringify(state), 'utf8') > MAX_STATE_BYTES) {
    throw new Error(`Expanded coordination state exceeds the ${MAX_STATE_BYTES}-byte snapshot limit.`);
  }
  return state;
}

/** Restore storage references without changing IDs, field order, or any historical text. */
export function restoreSnapshot(state) {
  if (!Array.isArray(state?.events) || !state.events.some(event => referenceKeys.some(key => Object.hasOwn(event, key)))) return checkExpandedCapacity(state);
  const checkpoints = checkpointIndex(state);
  const tasks = new Map(state.tasks.map(task => [task.id, task]));
  const restored = { ...state, events: state.events.map(event => {
    if (!referenceKeys.some(key => Object.hasOwn(event, key))) return event;
    const full = { ...event };
    const checkpoint = checkpoints.get(event.taskId)?.get(event.id);
    if (Object.hasOwn(event, 'detailReference')) {
      const reference = event.detailReference;
      if (event.kind !== 'checkpoint' || !validReference(reference)
        || event.detail !== detailPointer(event) || typeof checkpoint?.details !== 'string'
        || detailHash(checkpoint.details) !== reference.sha256) {
        throw new Error(`Invalid checkpoint detail reference for event ${event.id}. Full history was not modified.`);
      }
      full.detail = checkpoint.details;
      delete full.detailReference;
    }
    if (Object.hasOwn(event, 'checkpointMetadataReference')) {
      const reference = event.checkpointMetadataReference;
      if (event.kind !== 'checkpoint' || !validReference(reference)
        || event.summary !== checkpointSummary || !Array.isArray(event.links) || event.links.length !== 0
        || typeof checkpoint?.summary !== 'string' || !Array.isArray(checkpoint.evidence)
        || !checkpoint.evidence.every(link => typeof link === 'string')
        || detailHash(JSON.stringify([checkpoint.summary, checkpoint.evidence])) !== reference.sha256) {
        throw new Error(`Invalid checkpoint metadata reference for event ${event.id}. Full history was not modified.`);
      }
      full.summary = checkpoint.summary;
      full.links = [...checkpoint.evidence];
      delete full.checkpointMetadataReference;
    }
    if (Object.hasOwn(event, 'taskDescriptionReference')) {
      const reference = event.taskDescriptionReference, task = tasks.get(event.taskId);
      if (!validReference(reference) || !requestMatches(event, task) || event.detail !== requestPointer(event)
        || typeof task.description !== 'string' || detailHash(task.description) !== reference.sha256) {
        throw new Error(`Invalid task description reference for event ${event.id}. Full history was not modified.`);
      }
      full.detail = task.description;
      delete full.taskDescriptionReference;
    }
    return full;
  }) };
  return checkExpandedCapacity(restored);
}

/** Encode only exact duplicate text, keeping the full immutable checkpoint in this JSON. */
export function compactSnapshot(input) {
  // An older client can pass already-encoded records through unchanged. Validate
  // and restore them first, so a later modern write never nests or loses a reference.
  const state = restoreSnapshot(input);
  const checkpoints = checkpointIndex(state);
  const tasks = new Map((state.tasks ?? []).map(task => [task.id, task]));
  if (!Array.isArray(state?.events)) return state;
  return { ...state, events: state.events.map(event => {
    const checkpoint = checkpoints.get(event.taskId)?.get(event.id);
    let compact = event;
    if (event.kind === 'checkpoint' && typeof event.detail === 'string' && checkpoint?.details === event.detail) {
      const candidate = { ...compact, detail: detailPointer(event), detailReference: { version: 1, sha256: detailHash(event.detail) } };
      if (smaller(candidate, compact)) compact = candidate;
    }
    if (event.kind === 'checkpoint' && typeof event.summary === 'string' && checkpoint?.summary === event.summary
      && Array.isArray(checkpoint.evidence) && checkpoint.evidence.every(link => typeof link === 'string')
      && JSON.stringify(checkpoint.evidence) === JSON.stringify(event.links)) {
      const candidate = { ...compact, summary: checkpointSummary, links: [],
        checkpointMetadataReference: { version: 1, sha256: detailHash(JSON.stringify([event.summary, event.links])) } };
      if (smaller(candidate, compact)) compact = candidate;
    }
    const task = tasks.get(event.taskId);
    if (requestMatches(event, task) && typeof event.detail === 'string' && event.detail === task.description) {
      const candidate = { ...compact, detail: requestPointer(event), taskDescriptionReference: { version: 1, sha256: detailHash(event.detail) } };
      if (smaller(candidate, compact)) compact = candidate;
    }
    return compact;
  }) };
}

export async function decodeSnapshot(file, getBlob) {
  if (file?.type !== 'file' || !file.sha) throw new Error('Invalid coordination snapshot metadata');
  if (file.size !== undefined && (!Number.isSafeInteger(file.size) || file.size < 0 || file.size > MAX_STATE_BYTES)) {
    throw new Error('Invalid coordination snapshot size');
  }
  let blob = file;
  if (file.encoding === 'none' && file.content === '') {
    if (file.size === undefined || !/^[0-9a-f]{40}$/.test(file.sha)) throw new Error('Invalid coordination blob identity');
    blob = await getBlob(file.sha);
    if (blob.sha !== file.sha || blob.size !== file.size) throw new Error('Coordination blob does not match its snapshot');
  }
  if (blob.encoding !== 'base64' || typeof blob.content !== 'string' || !blob.content) throw new Error('Invalid coordination snapshot content');
  const bytes = Buffer.from(blob.content, 'base64');
  if (bytes.length > MAX_STATE_BYTES || (blob !== file && bytes.length !== file.size)) throw new Error('Invalid coordination snapshot size');
  let state;
  try { state = JSON.parse(bytes.toString('utf8')); }
  catch { throw new Error('Coordination state is not valid JSON.'); }
  return restoreSnapshot(state);
}


// Keep each public record on a separate Git diff line and avoid repeating indentation
// through the full event history. Older readers see plain JSON and explicit
// pointers; the full detail remains in the matching task checkpoint.
export function serializeSnapshot(input) {
  const state = compactSnapshot(input);
  const fields = Object.entries(state).map(([key, value]) => {
    const encoded = Array.isArray(value) && value.length
      ? '[\n' + value.map(record => '    ' + JSON.stringify(record)).join(',\n') + '\n  ]'
      : JSON.stringify(value);
    return '  ' + JSON.stringify(key) + ': ' + encoded;
  });
  const text = '{\n' + fields.join(',\n') + '\n}\n';
  if (Buffer.byteLength(text, 'utf8') > MAX_STATE_BYTES) {
    throw new Error('Coordination state exceeds the 3000000-byte snapshot limit after lossless detail references. No history was removed.');
  }
  return text;
}
