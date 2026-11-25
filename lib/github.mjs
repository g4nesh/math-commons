import { execFileSync } from 'node:child_process';
import { randomUUID } from 'node:crypto';
import { applyOperation, validateState, MAX_STATE_BYTES } from './core.mjs';

export class GitHubError extends Error {
  constructor(message, status = 0, retryAfter = null) { super(message); this.name = 'GitHubError'; this.status = status; this.retryAfter = retryAfter; }
}
export function resolveToken({ required = true } = {}) {
  const environment = process.env.GH_TOKEN || process.env.GITHUB_TOKEN;
  if (environment) return environment.trim();
  try { return execFileSync('gh', ['auth', 'token', '--hostname', 'github.com'], { encoding: 'utf8', stdio: ['ignore', 'pipe', 'ignore'], timeout: 10_000 }).trim(); }
  catch { if (required) throw new Error('Authenticate with `gh auth login`, or set GH_TOKEN to your own GitHub token with access to this repository.'); }
  return null;
}
const delay = ms => new Promise(resolve => setTimeout(resolve, ms));
export function validateConfig(config) {
  if (!/^[A-Za-z0-9][A-Za-z0-9-]{0,38}\/[A-Za-z0-9_.-]+$/.test(config?.repository ?? '') || ['.', '..'].includes(config.repository.split('/')[1])) throw new Error('config.repository must be owner/repository.');
  if (!/^[A-Za-z0-9_./-]+$/.test(config.stateBranch ?? '') || config.stateBranch.includes('..')) throw new Error('Invalid coordination branch.');
  if (!/^[A-Za-z0-9_./-]+\.json$/.test(config.statePath ?? '') || config.statePath.includes('..')) throw new Error('Invalid coordination state path.');
  if (config.stateBranch === 'main' || config.stateBranch === 'master') throw new Error('Coordination state must use a separate branch from source.');
  if (config.defaultObjectiveId !== undefined && (typeof config.defaultObjectiveId !== 'string' || !/^[a-z0-9][a-z0-9-]{0,79}$/.test(config.defaultObjectiveId))) throw new Error('config.defaultObjectiveId must be a valid objective ID.');
  return config;
}

/** SHA-conditional updates serialize concurrent claims. This is a cooperative protocol for trusted repo collaborators. */
export class GitHubStore {
  constructor(config, { token, fetchImpl = globalThis.fetch, sleep = delay, maxRetries = 5, timeoutMs = 20_000 } = {}) {
    this.config = validateConfig(config); this.token = token; this.fetch = fetchImpl; this.sleep = sleep; this.maxRetries = maxRetries; this.timeoutMs = timeoutMs;
    this.base = `https://api.github.com/repos/${config.repository}`;
    this.path = config.statePath.split('/').map(encodeURIComponent).join('/');
  }
  async request(url, { method = 'GET', body } = {}) {
    const headers = { Accept: 'application/vnd.github+json', 'X-GitHub-Api-Version': '2022-11-28', 'User-Agent': 'math-commons/1.0', 'Cache-Control': 'no-cache' };
    if (this.token) headers.Authorization = `Bearer ${this.token}`;
    if (body !== undefined) headers['Content-Type'] = 'application/json';
    let response;
    try { response = await this.fetch(url, { method, headers, body: body === undefined ? undefined : JSON.stringify(body), signal: AbortSignal.timeout(this.timeoutMs) }); }
    catch { throw new GitHubError(`GitHub ${method} request did not return a response.`, 0); }
    let payload;
    try { payload = await response.json(); }
    catch { throw new GitHubError(`GitHub returned an unreadable response (${response.status}).`, response.status); }
    if (!response.ok) {
      const reasons = { 401: 'GitHub authentication failed. Run gh auth login.', 403: 'GitHub denied access or rate limited this account. Check repository write permissions and rate limits.', 404: 'Repository or coordination state was not found. Check access and initialize the coordination branch.', 409: 'Coordination state changed concurrently.', 422: 'GitHub rejected the state update; the state may have changed concurrently.', 429: 'GitHub rate limited this account.' };
      throw new GitHubError(reasons[response.status] ?? `GitHub request failed (${response.status}).`, response.status, response.headers.get('retry-after'));
    }
    return { payload, serverTime: response.headers.get('date') || new Date().toISOString() };
  }
  async identity() {
    if (!this.token) throw new Error('Publishing requires your own GitHub authentication.');
    const { payload } = await this.request('https://api.github.com/user');
    if (!payload.login) throw new Error('GitHub did not return an authenticated account.');
    return payload.login;
  }
  async read() {
    const { payload, serverTime } = await this.request(`${this.base}/contents/${this.path}?ref=${encodeURIComponent(this.config.stateBranch)}`);
    if (payload.type !== 'file' || payload.encoding !== 'base64' || !payload.sha || !payload.content) throw new Error('Coordination state must be a small JSON file readable through GitHub Contents.');
    if (payload.size > MAX_STATE_BYTES) throw new Error('Coordination state exceeds the supported size limit. Ask a maintainer to archive history.');
    let state;
    try { state = JSON.parse(Buffer.from(payload.content, 'base64').toString('utf8')); }
    catch { throw new Error('Coordination state is not valid JSON.'); }
    validateState(state);
    return { state, sha: payload.sha, serverTime };
  }
  async write(state, sha, message) {
    return this.request(`${this.base}/contents/${this.path}`, { method: 'PUT', body: { message, content: Buffer.from(JSON.stringify(state, null, 2) + '\n').toString('base64'), sha, branch: this.config.stateBranch } });
  }
  async mutate(operation) {
    const op = { ...operation, id: operation.id || randomUUID() };
    const owner = await this.identity();
    let uncertain = false, lastError;
    for (let attempt = 0; attempt <= this.maxRetries; attempt++) {
      let snapshot;
      try { snapshot = await this.read(); }
      catch (error) {
        if (uncertain) {
          const failure = new Error(`The last write may have succeeded, and its state could not be re-read. Do not submit a new request. Retry the identical command with --request-id ${op.id} to reconcile.`);
          failure.code = 'WRITE_UNCERTAIN'; failure.requestId = op.id; throw failure;
        }
        throw error;
      }
      const transition = applyOperation(snapshot.state, op, { owner, now: snapshot.serverTime, leaseMinutes: this.config.leaseMinutes ?? 20 });
      if (transition.duplicate) return { ...transition, requestId: op.id, reconciled: uncertain };
      try {
        const result = await this.write(transition.state, snapshot.sha, `commons(${op.kind}): ${op.agent} [${op.id}]`);
        return { ...transition, requestId: op.id, commit: result.payload.commit?.sha, reconciled: false };
      } catch (error) {
        lastError = error;
        if (!(error instanceof GitHubError) || ![0, 409, 422, 500, 502, 503, 504].includes(error.status)) throw error;
        uncertain ||= error.status === 0 || error.status >= 500;
        if (attempt < this.maxRetries) await this.sleep(Math.min(4_000, 150 * 2 ** attempt + Math.random() * 250));
      }
    }
    // One last read reconciles a successful write with a lost response, even on the final attempt.
    try {
      const snapshot = await this.read();
      const known = snapshot.state.events.find(event => event.id === op.id);
      if (known) {
        const transition = applyOperation(snapshot.state, op, { owner, now: snapshot.serverTime, leaseMinutes: this.config.leaseMinutes ?? 20 });
        return { ...transition, requestId: op.id, reconciled: true };
      }
      if (!uncertain) throw lastError;
    } catch (error) { if (!uncertain) throw error; }
    const error = new Error(`Write outcome is uncertain. Retry the identical command with --request-id ${op.id}; do not create a new operation ID.`);
    error.code = 'WRITE_UNCERTAIN'; error.requestId = op.id; throw error;
  }
}
