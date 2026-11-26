import { createHash } from 'node:crypto';

// Canonical reads/writes retain the existing capacity. The legacy deployment's
// smaller limit is a measurement for rollout checks, not a research-write cap.
export const MAX_STATE_BYTES = 3_000_000;
export const LEGACY_READER_MAX_BYTES = 950_000;

const detailHash = value => createHash('sha256').update(value, 'utf8').digest('hex');
const detailPointer = event => `Full detail: task ${event.taskId}, checkpoint ${event.id}.`;

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

/** Restore storage references without changing IDs, field order, or any historical text. */
export function restoreSnapshot(state) {
  if (!Array.isArray(state?.events) || !state.events.some(event => Object.hasOwn(event, 'detailReference'))) return state;
  const checkpoints = checkpointIndex(state);
  const restored = { ...state, events: state.events.map(event => {
    if (!Object.hasOwn(event, 'detailReference')) return event;
    const reference = event.detailReference;
    const checkpoint = checkpoints.get(event.taskId)?.get(event.id);
    if (event.kind !== 'checkpoint' || !reference || typeof reference !== 'object' || Array.isArray(reference)
      || Object.keys(reference).sort().join(',') !== 'sha256,version' || reference.version !== 1
      || typeof reference.sha256 !== 'string' || !/^[0-9a-f]{64}$/.test(reference.sha256)
      || event.detail !== detailPointer(event) || typeof checkpoint?.details !== 'string'
      || detailHash(checkpoint.details) !== reference.sha256) {
      throw new Error(`Invalid checkpoint detail reference for event ${event.id}. Full history was not modified.`);
    }
    const full = { ...event, detail: checkpoint.details };
    delete full.detailReference;
    return full;
  }) };
  if (Buffer.byteLength(JSON.stringify(restored), 'utf8') > MAX_STATE_BYTES) throw new Error('Expanded coordination snapshot is too large');
  return restored;
}

/** Encode only exact duplicate text, keeping the full immutable checkpoint in this JSON. */
export function compactSnapshot(input) {
  // An older client can pass already-encoded records through unchanged. Validate
  // and restore them first, so a later modern write never nests or loses a reference.
  const state = restoreSnapshot(input);
  const checkpoints = checkpointIndex(state);
  if (!Array.isArray(state?.events)) return state;
  return { ...state, events: state.events.map(event => {
    const checkpoint = checkpoints.get(event.taskId)?.get(event.id);
    if (event.kind !== 'checkpoint' || typeof event.detail !== 'string' || checkpoint?.details !== event.detail) return event;
    const compact = { ...event, detail: detailPointer(event), detailReference: { version: 1, sha256: detailHash(event.detail) } };
    return Buffer.byteLength(JSON.stringify(compact), 'utf8') < Buffer.byteLength(JSON.stringify(event), 'utf8') ? compact : event;
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
