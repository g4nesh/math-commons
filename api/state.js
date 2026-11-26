import { decodeSnapshot } from '../lib/snapshot.mjs';
import config from '../commons.config.json' with { type: 'json' };

export async function loadSnapshot({ fetchImpl = fetch, now = Date.now() } = {}) {
  const { repository, stateBranch, statePath } = config;
  const path = statePath.split('/').map(encodeURIComponent).join('/');
  const stateUrl = `https://github.com/${repository}/blob/${encodeURIComponent(stateBranch)}/${path}`;
  // A moving raw.githubusercontent.com branch can remain stale at an upstream
  // cache even with changing query strings. Read the versioned Contents API.
  const contentsUrl = `https://api.github.com/repos/${repository}/contents/${path}?ref=${encodeURIComponent(stateBranch)}`;
  async function readJson(url, accept = 'application/vnd.github+json') {
    const response = await fetchImpl(url, { headers: {
      Accept: accept, 'X-GitHub-Api-Version': '2022-11-28',
      'User-Agent': 'math-commons-dashboard/1.0', 'Cache-Control': 'no-cache'
    }, signal: AbortSignal.timeout(10000) });
    if (!response.ok) throw new Error(`GitHub returned ${response.status}`);
    return response.json();
  }
  const contents = await readJson(contentsUrl, 'application/vnd.github.object+json');
  const state = await decodeSnapshot(contents, sha => readJson(`https://api.github.com/repos/${repository}/git/blobs/${sha}`));
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
    // Two-minute sharing leaves headroom below the public hourly allowance.
    res.setHeader('Cache-Control', 'public, max-age=0, s-maxage=120');
    return res.status(200).json(snapshot);
  } catch {
    res.setHeader('Cache-Control', 'no-store');
    return res.status(502).json({ error: 'The GitHub work log is temporarily unavailable. Refresh to try again; no local example has been substituted.' });
  }
}
