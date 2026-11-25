// Leave room for provenance and other fields in a Vercel function response.
export const MAX_STATE_BYTES = 3_000_000;

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
  try { return JSON.parse(bytes.toString('utf8')); }
  catch { throw new Error('Coordination state is not valid JSON.'); }
}


// Keep each public record on a separate Git diff line and avoid repeating indentation
// through the full event history. Older readers consume the same plain JSON schema.
export function serializeSnapshot(state) {
  const fields = Object.entries(state).map(([key, value]) => {
    const encoded = Array.isArray(value) && value.length
      ? '[\n' + value.map(record => '    ' + JSON.stringify(record)).join(',\n') + '\n  ]'
      : JSON.stringify(value);
    return '  ' + JSON.stringify(key) + ': ' + encoded;
  });
  return '{\n' + fields.join(',\n') + '\n}\n';
}
