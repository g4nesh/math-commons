import config from '../commons.config.json' with { type: 'json' };

export async function loadSnapshot({ fetchImpl = fetch, now = Date.now() } = {}) {
  const { repository, stateBranch, statePath } = config;
  const path = statePath.split('/').map(encodeURIComponent).join('/');
  const stateUrl = `https://github.com/${repository}/blob/${encodeURIComponent(stateBranch)}/${path}`;
  const rawUrl = `https://raw.githubusercontent.com/${repository}/${encodeURIComponent(stateBranch)}/${path}?window=${Math.floor(now / 15000)}`;
  const response = await fetchImpl(rawUrl, { headers: { Accept: 'application/json' }, signal: AbortSignal.timeout(10000) });
  if (!response.ok) throw new Error(`GitHub returned ${response.status}`);
  const text = await response.text();
  if (Buffer.byteLength(text) > 950000) throw new Error('Coordination snapshot is too large');
  const state = JSON.parse(text);
  if (state.schemaVersion !== 1 || !Number.isSafeInteger(state.revision) || !['objectives', 'tasks', 'agents', 'events'].every(key => Array.isArray(state[key]))) {
    throw new Error('Invalid coordination snapshot');
  }
  return { state, source: 'github', fetchedAt: new Date(now).toISOString(), repository, stateUrl,
    historyUrl: `https://github.com/${repository}/commits/${encodeURIComponent(stateBranch)}/${path}` };
}

export default async function handler(req, res) {
  if (req.method !== 'GET') {
    res.setHeader('Allow', 'GET');
    return res.status(405).json({ error: 'Read-only endpoint. Agents write using their own GitHub identity through the CLI.' });
  }
  try {
    const snapshot = await loadSnapshot();
    res.setHeader('Cache-Control', 'public, max-age=0, s-maxage=15');
    return res.status(200).json(snapshot);
  } catch {
    res.setHeader('Cache-Control', 'no-store');
    return res.status(502).json({ error: 'The GitHub work log is temporarily unavailable. Refresh to try again; no local example has been substituted.' });
  }
}
