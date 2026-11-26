import { readFile } from 'node:fs/promises';
import { createHash } from 'node:crypto';
import { gunzipSync } from 'node:zlib';
import { pathToFileURL } from 'node:url';
import { restoreSnapshot, serializeSnapshot, LEGACY_READER_MAX_BYTES } from '../lib/snapshot.mjs';

const hash = text => createHash('sha256').update(text).digest('hex');

/** Read-only audit; compare every semantic byte, not just counts or selected fields. */
export function auditSnapshot(bytes) {
  const original = restoreSnapshot(JSON.parse(bytes.toString('utf8')));
  const encoded = serializeSnapshot(original);
  const stored = JSON.parse(encoded);
  const restored = restoreSnapshot(stored);
  const before = JSON.stringify(original), after = JSON.stringify(restored);
  if (before !== after) throw new Error('Snapshot round trip changed semantic bytes.');
  const counts = state => ({ objectives: state.objectives.length, tasks: state.tasks.length,
    agents: state.agents.length, events: state.events.length, usage: (state.usage ?? []).length,
    checkpoints: state.tasks.reduce((sum, task) => sum + (task.checkpoints ?? []).length, 0),
    reviews: state.tasks.reduce((sum, task) => sum + (task.reviews ?? []).length, 0) });
  return { revision: original.revision, originalBytes: bytes.length, storedBytes: Buffer.byteLength(encoded),
    legacyReaderLimit: LEGACY_READER_MAX_BYTES, legacyCompatible: Buffer.byteLength(encoded) <= LEGACY_READER_MAX_BYTES,
    headroomBytes: LEGACY_READER_MAX_BYTES-Buffer.byteLength(encoded),
    originalRawSha256: hash(bytes), originalSemanticSha256: hash(before), restoredSemanticSha256: hash(after),
    semanticBytesIdentical: before === after, referenceCount: stored.events.filter(event => event.detailReference).length,
    checkpointMetadataReferences: stored.events.filter(event => event.checkpointMetadataReference).length,
    requestDescriptionReferences: stored.events.filter(event => event.taskDescriptionReference).length,
    recordsBefore: counts(original), recordsAfter: counts(restored),
    legacyPrettyBytes: Buffer.byteLength(JSON.stringify(stored, null, 2)+'\n'),
    legacyBehavior: 'Plain JSON readers see explicit event pointers and complete checkpoint details. Writers retaining unknown event fields preserve references; pre-950k pretty-print clients refuse oversized writes.' };
}

if (process.argv[1] && import.meta.url === pathToFileURL(process.argv[1]).href) {
  const path = process.argv[2];
  if (!path) throw new Error('Usage: node scripts/audit-snapshot.mjs PATH.json[.gz]');
  let bytes = await readFile(path);
  if (path.endsWith('.gz')) bytes = gunzipSync(bytes);
  console.log(JSON.stringify(auditSnapshot(bytes), null, 2));
}
