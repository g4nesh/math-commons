import { spawnSync } from 'node:child_process';
import { mkdirSync, mkdtempSync, readFileSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { createCampaign, loadConfig, saveCampaign } from '../lib/campaign.mjs';
import { STATEMENT } from './fixtures/euler.mjs';

export const HERE = resolve(fileURLToPath(new URL('.', import.meta.url)));
export const FAKE_RUNTIME = join(HERE, 'fixtures', 'fake-runtime.mjs');
export const FAKE_LAKE = join(HERE, 'fixtures', 'fake-lake.mjs');

// Real Lean when the toolchain is installed (CI), the fake otherwise.
export const REAL_LAKE = !process.env.DISCOVER_FORCE_FAKE_LAKE && spawnSync('lake', ['--version'], { encoding: 'utf8' }).status === 0;

export function tempConfig(overrides = {}) {
  // DISCOVER_TEST_KEEP keeps the generated campaigns (CI uploads them as an artifact).
  const parent = process.env.DISCOVER_TEST_KEEP || tmpdir();
  mkdirSync(parent, { recursive: true });
  const root = mkdtempSync(join(parent, 'discover-test-'));
  const config = loadConfig({ root });
  config.runtimes.claude.bin = FAKE_RUNTIME;
  config.runtimes.codex.bin = FAKE_RUNTIME;
  if (!REAL_LAKE) process.env.DISCOVER_LAKE_BIN = FAKE_LAKE;
  return Object.assign(config, overrides);
}

export function eulerCampaign(config, { slug = 'euler-poly', formal = true } = {}) {
  const campaign = createCampaign(config, slug, { title: "Euler's polynomial", source: 'Euler (1772); see https://oeis.org/A005846', fixture: true, statement: STATEMENT });
  const path = campaign.paths.statement;
  writeFileSync(path, readFileSync(path, 'utf8').replace(/TODO: every definition/, 'Every definition').replace(/TODO: what a complete/, 'What a complete'));
  if (formal) { campaign.state.formalRequested = true; saveCampaign(campaign); }
  return campaign;
}

export function withEnv(values, fn) {
  const previous = Object.fromEntries(Object.keys(values).map(key => [key, process.env[key]]));
  Object.assign(process.env, values);
  return Promise.resolve(fn()).finally(() => {
    for (const [key, value] of Object.entries(previous)) {
      if (value === undefined) delete process.env[key];
      else process.env[key] = value;
    }
  });
}
