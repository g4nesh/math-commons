import config from '../commons.config.json' with { type: 'json' };

export async function loadSnapshot({ fetchImpl = fetch, now = Date.now() } = {}) {
  const { repository, stateBranch, statePath } = config;
  const path = statePath.split('/').map(encodeURIComponent).join('/');
  const stateUrl = `https://github.com/${repository}/blob/${encodeURIComponent(stateBranch)}/${path}`;
  // A moving raw.githubusercontent.com branch can remain stale at an upstream
  // cache even with changing query strings. Read the versioned Contents API.
  const contentsUrl = `https://api.github.com/repos/${repository}/contents/${path}?ref=${encodeURIComponent(stateBranch)}`;
  const response = await fetchImpl(contentsUrl, { headers: {
    Accept: 'application/vnd.github+json',
    'X-GitHub-Api-Version': '2022-11-28',
    'User-Agent': 'math-commons-dashboard/1.0',
    'Cache-Control': 'no-cache'
  }, signal: AbortSignal.timeout(10000) });
  if (!response.ok) throw new Error(`GitHub returned ${response.status}`);
  const contents = await response.json();
  if (contents.type !== 'file' || contents.encoding !== 'base64' || typeof contents.content !== 'string' || !contents.sha || contents.size > 950000) {
    throw new Error('Invalid coordination snapshot metadata');
  }
  const text = Buffer.from(contents.content, 'base64').toString('utf8');
  if (Buffer.byteLength(text) > 950000) throw new Error('Coordination snapshot is too large');
  const state = JSON.parse(text);
  if (state.schemaVersion !== 1 || !Number.isSafeInteger(state.revision) || !['objectives', 'tasks', 'agents', 'events'].every(key => Array.isArray(state[key]))) {
    throw new Error('Invalid coordination snapshot');
  }
  return { state, source: 'github', contentSha: contents.sha, fetchedAt: new Date(now).toISOString(), repository, stateUrl,
    historyUrl: `https://github.com/${repository}/commits/${encodeURIComponent(stateBranch)}/${path}` };
}

export default async function handler(req, res) {
  if (req.method !== 'GET') {
    res.setHeader('Allow', 'GET');
    return res.status(405).json({ error: 'Read-only endpoint. Agents write using their own GitHub identity through the CLI.' });
  }
  try {
    const snapshot = await loadSnapshot();
    // Public unauthenticated GitHub reads have a limited hourly allowance.
    // Share snapshots at the CDN for a minute; no token is exposed or required.
    res.setHeader('Cache-Control', 'public, max-age=0, s-maxage=60');
    return res.status(200).json(snapshot);
  } catch {
    res.setHeader('Cache-Control', 'no-store');
    return res.status(502).json({ error: 'The GitHub work log is temporarily unavailable. Refresh to try again; no local example has been substituted.' });
  }
}
