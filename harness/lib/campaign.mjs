import { existsSync, readdirSync, mkdirSync, cpSync } from 'node:fs';
import { join, resolve } from 'node:path';
import {
  HARNESS_ROOT, REPO_ROOT, HarnessError, assertSlug, appendJsonl, fileSha256, now, readJson,
  readJsonl, readText, writeJson, writeText
} from './util.mjs';

export const STAGES = ['screen', 'map', 'explore', 'draft', 'review', 'formalize', 'publish', 'done'];
export const TERMINAL = ['rejected', 'parked', 'done'];

function merge(base, override) {
  if (!override || typeof override !== 'object' || Array.isArray(override)) return override ?? base;
  const out = { ...base };
  for (const [key, value] of Object.entries(override)) {
    out[key] = value && typeof value === 'object' && !Array.isArray(value) && base?.[key] && typeof base[key] === 'object'
      ? merge(base[key], value) : value;
  }
  return out;
}

// Defaults live in harness/config.json. An operator overrides them in
// discovery.config.json at the repository root (or --config) without editing the harness.
export function loadConfig({ root = process.env.DISCOVER_ROOT || REPO_ROOT, configPath } = {}) {
  let config = readJson(join(HARNESS_ROOT, 'config.json'));
  const local = configPath ? resolve(configPath) : join(root, 'discovery.config.json');
  if (existsSync(local)) config = merge(config, readJson(local));
  for (const name of ['claude', 'codex']) {
    const variable = process.env[`DISCOVER_${name.toUpperCase()}_BIN`];
    if (variable) config.runtimes[name].bin = variable;
  }
  return { ...config, root };
}

export const campaignsRoot = config => resolve(config.root, config.campaignRoot);
export const campaignDir = (config, slug) => join(campaignsRoot(config), assertSlug(slug));

export function paths(dir) {
  return {
    dir,
    state: join(dir, 'campaign.json'),
    statement: join(dir, 'STATEMENT.md'),
    screen: join(dir, 'SCREEN.md'),
    screenJson: join(dir, 'screen.json'),
    map: join(dir, 'MAP.md'),
    directions: join(dir, 'directions.json'),
    attempts: join(dir, 'attempts'),
    ledger: join(dir, 'ledger.jsonl'),
    certificates: join(dir, 'certificates'),
    proof: join(dir, 'PROOF.md'),
    reviews: join(dir, 'reviews'),
    formal: join(dir, 'formal'),
    receipts: join(dir, 'receipts'),
    capsule: join(dir, 'CAPSULE.md'),
    runs: join(dir, 'runs.jsonl'),
    runDir: join(dir, '.runs'),
    outbox: join(dir, 'outbox.jsonl')
  };
}

export function listCampaigns(config) {
  const root = campaignsRoot(config);
  if (!existsSync(root)) return [];
  return readdirSync(root, { withFileTypes: true })
    .filter(entry => entry.isDirectory() && existsSync(join(root, entry.name, 'campaign.json')))
    .map(entry => entry.name).sort();
}

export function loadCampaign(config, slug) {
  const dir = campaignDir(config, slug);
  const p = paths(dir);
  if (!existsSync(p.state)) throw new HarnessError('NO_CAMPAIGN', `No campaign ${slug} under ${config.campaignRoot}/`);
  const state = readJson(p.state);
  if (state.version !== 1 || state.slug !== slug) throw new HarnessError('BAD_STATE', `campaign.json does not describe ${slug}`);
  return { config, slug, dir, paths: p, state };
}

export function saveCampaign(campaign) {
  campaign.state.updatedAt = now();
  writeJson(campaign.paths.state, campaign.state);
}

export function createCampaign(config, slug, { title, source = '', polarity = 'unknown', fixture = false, statement = '' } = {}) {
  assertSlug(slug);
  const dir = campaignDir(config, slug);
  if (existsSync(join(dir, 'campaign.json'))) throw new HarnessError('EXISTS', `Campaign ${slug} already exists`);
  if (!['prove', 'disprove', 'unknown'].includes(polarity)) throw new HarnessError('BAD_POLARITY', 'Polarity is prove, disprove or unknown');
  mkdirSync(dir, { recursive: true });
  const template = join(HARNESS_ROOT, 'templates', 'campaign');
  for (const name of ['STATEMENT.md', 'AGENTS.md', 'CLAUDE.md']) {
    const target = join(dir, name);
    if (!existsSync(target)) {
      const text = readText(join(template, name))
        .replaceAll('{{TITLE}}', title || slug).replaceAll('{{SLUG}}', slug)
        .replaceAll('{{SOURCE}}', source || 'TODO: primary source (arXiv/DOI/URL and statement number)')
        .replaceAll('{{STATEMENT}}', statement || 'TODO: the exact statement, quoted from the primary source.');
      writeText(target, text);
    }
  }
  const state = {
    version: 1, slug, title: title || slug, fixture: Boolean(fixture), polarity,
    stage: 'screen', status: 'active', createdAt: now(), updatedAt: now(),
    statementSha256: null,
    lean: null,
    budget: { tokens: config.budgets.campaignTokens },
    commons: null,
    counters: { exploreRounds: 0, repairRounds: 0, formalRounds: 0 },
    history: [{ at: now(), from: null, to: 'screen', reason: 'created' }]
  };
  writeJson(join(dir, 'campaign.json'), state);
  return loadCampaign(config, slug);
}

export function transition(campaign, to, reason) {
  if (![...STAGES, ...TERMINAL].includes(to)) throw new HarnessError('BAD_STAGE', `Unknown stage ${to}`);
  const from = campaign.state.stage;
  if (from === to) return;
  campaign.state.history.push({ at: now(), from, to, reason });
  campaign.state.stage = to;
  campaign.state.status = TERMINAL.includes(to) ? to : 'active';
  saveCampaign(campaign);
}

export function recordRun(campaign, run) {
  appendJsonl(campaign.paths.runs, { at: now(), ...run });
}

export function usageTotals(campaign) {
  let tokens = 0, unknown = 0;
  for (const run of readJsonl(campaign.paths.runs)) {
    if (run.usage) tokens += (run.usage.inputTokens || 0) + (run.usage.outputTokens || 0);
    else if (run.runtime !== 'manual') unknown++;
  }
  return { tokens, unknown, runs: readJsonl(campaign.paths.runs).length };
}

export const proofHash = campaign => (existsSync(campaign.paths.proof) ? fileSha256(campaign.paths.proof) : null);
export const statementHash = campaign => (existsSync(campaign.paths.statement) ? fileSha256(campaign.paths.statement) : null);

export function copyTemplateFile(name, target) {
  if (!existsSync(target)) cpSync(join(HARNESS_ROOT, 'templates', 'campaign', name), target);
}
