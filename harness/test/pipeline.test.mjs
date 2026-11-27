import assert from 'node:assert/strict';
import { existsSync, readFileSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { test } from 'node:test';
import { loadCampaign } from '../lib/campaign.mjs';
import { evaluate, reviewSummary, verificationLabel } from '../lib/gates.mjs';
import { nextAction, recordManual, runLoop, step } from '../lib/orchestrator.mjs';
import { reconcile } from '../lib/reconcile.mjs';
import { readJson, readJsonl } from '../lib/util.mjs';
import { REAL_LAKE, eulerCampaign, tempConfig, withEnv } from './helpers.mjs';

const reload = campaign => loadCampaign(campaign.config, campaign.slug);

test('a fixture campaign runs from screening to publication', { timeout: 600_000 }, async () => {
  const config = tempConfig();
  const campaign = eulerCampaign(config);
  const log = await runLoop(campaign, { maxSteps: 40 });
  const done = reload(campaign);
  assert.equal(done.state.stage, 'done', JSON.stringify(log.map(entry => entry.action || entry), null, 1));

  // Statement frozen at screening; the screen record is the model's, the hash is ours.
  assert.match(done.state.statementSha256, /^[0-9a-f]{64}$/);
  assert.equal(readJson(done.paths.screenJson).verdict, 'fixture');

  // Exploration produced ledger entries for both directions and a certificate with a passing receipt.
  const ledger = readJsonl(done.paths.ledger);
  assert.deepEqual([...new Set(ledger.map(entry => entry.direction))].sort(), ['D1', 'D2']);
  const origin = readJson(join(done.paths.certificates, 'euler-forty', 'provenance.json'));
  assert.ok(origin.generator.runId && origin.verifier.runId);
  assert.notEqual(origin.generator.runId, origin.verifier.runId, 'the verifier came from a separate run');
  assert.notEqual(origin.generator.identity, origin.verifier.identity, 'the verifier came from the other runtime');
  const certificateReceipt = readJson(join(done.paths.certificates, 'euler-forty', 'receipt.json'));
  assert.equal(certificateReceipt.status, 'pass');
  assert.deepEqual(certificateReceipt.warnings, []);

  // Reviews are bound to the exact proof hash and include an identity other than the author's.
  const reviews = reviewSummary(done);
  assert.equal(reviews.passes.length, 2);
  assert.ok([...reviews.identities].some(identity => identity !== done.state.author.identity));
  assert.ok(reviews.passes.every(review => review.proofSha256 === reviews.hash));

  // Lean: statement audited and pinned before the proof; the receipt matches the pin.
  assert.equal(done.state.lean.audit.verdict, 'faithful');
  const receipt = readJson(join(done.paths.receipts, 'lean-latest.json'));
  assert.equal(receipt.status, 'pass', receipt.failures.join('; '));
  assert.equal(receipt.statement.fingerprint, done.state.lean.claimFingerprint);
  assert.equal(receipt.mainTheorem.refutesClaim, true);
  if (REAL_LAKE) assert.ok(receipt.mainTheorem.axioms.every(axiom => ['propext', 'Classical.choice', 'Quot.sound'].includes(axiom)));

  // Publication is generated from the record, with a computed label.
  assert.equal(verificationLabel(done).emoji, '✅');
  const page = readFileSync(join(config.root, 'results', 'euler-poly', 'README.md'), 'utf8');
  assert.match(page, /Harness fixture/);
  assert.match(page, /review-001/);
  assert.ok(existsSync(join(config.root, 'results', 'README.md')));

  // Every model run is recorded with runtime usage; the capsule stays small.
  const runs = readJsonl(done.paths.runs);
  assert.ok(runs.length >= 10);
  assert.ok(runs.every(run => run.status === 'ok' && run.usage?.inputTokens > 0));
  assert.ok(readFileSync(done.paths.capsule, 'utf8').length < 6000 * 4);
  assert.deepEqual(reconcile(config).blockers, []);
});

test('a blocking review defect forces a repair and voids the earlier verdict', { timeout: 600_000 }, async () => {
  const config = tempConfig();
  const campaign = eulerCampaign(config, { formal: false });
  await withEnv({ FAKE_RUNTIME_MODE: 'review-changes-once' }, () => runLoop(campaign, { until: 'formalize', maxSteps: 40 }));
  const after = reload(campaign);
  assert.equal(after.state.stage, 'formalize');
  assert.equal(after.state.counters.repairRounds, 1);
  const reviews = reviewSummary(after);
  assert.equal(reviews.all[0].verdict, 'changes');
  assert.ok(reviews.stale >= 1, 'the pre-repair review is stale');
  assert.match(readFileSync(after.paths.proof, 'utf8'), /## Revision notes/);
});

test('a prior resolution is rejected at screening and never explored', async () => {
  const config = tempConfig();
  const campaign = eulerCampaign(config, { formal: false });
  await withEnv({ FAKE_RUNTIME_MODE: 'reject' }, () => runLoop(campaign, { maxSteps: 5 }));
  const after = reload(campaign);
  assert.equal(after.state.stage, 'rejected');
  assert.equal(readJsonl(after.paths.runs).length, 1);
  assert.equal(nextAction(after).kind, 'stop');
});

test('a run that writes outside its scope is not ingested', async () => {
  const config = tempConfig();
  const campaign = eulerCampaign(config, { formal: false });
  await withEnv({ FAKE_RUNTIME_MODE: 'scope-violation' }, async () => {
    await step(campaign);
    await step(campaign);
    const outcome = await step(campaign);
    assert.equal(outcome.results[0].role, 'map');
    assert.equal(outcome.results[0].status, 'scope-violation');
    assert.match(outcome.results[0].detail, /campaign-notes\.md/);
  });
  assert.equal(existsSync(reload(campaign).paths.directions), false);
});

test('an edited statement after screening is caught', async () => {
  const config = tempConfig();
  const campaign = eulerCampaign(config, { formal: false });
  await runLoop(campaign, { until: 'explore', maxSteps: 5 });
  const after = reload(campaign);
  writeFileSync(after.paths.statement, readFileSync(after.paths.statement, 'utf8').replace('is prime', 'is prime or a square'));
  assert.ok(reconcile(config).blockers.some(item => /STATEMENT\.md changed after screening/.test(item)));
});

test('manual runs pause the loop and resume from a recorded answer', async () => {
  const config = tempConfig();
  const campaign = eulerCampaign(config, { formal: false });
  const outcome = await step(campaign, { runtime: 'manual', as: 'claude-interactive' });
  assert.equal(outcome.results[0].status, 'pending');
  const pending = reload(campaign).state.pending;
  const runDir = join(campaign.paths.runDir, pending.runId);
  const prompt = readFileSync(join(runDir, 'prompt.md'), 'utf8');
  assert.match(prompt, /^You are one step of the Math Commons discovery harness/);
  assert.match(prompt, /ROLE: screen/);
  // Play the interactive agent: write the files and the answer the prompt asks for.
  const { roleOutput } = await import('./fixtures/euler.mjs');
  const output = roleOutput('screen', { runId: pending.runId });
  for (const [path, text] of Object.entries(output.files)) writeFileSync(join(campaign.dir, path), text);
  writeFileSync(join(runDir, 'answer.json'), JSON.stringify(output.answer));
  assert.equal(nextAction(reload(campaign)).kind, 'stop');
  const result = recordManual(reload(campaign), pending.runId);
  assert.equal(result.status, 'ok');
  const resumed = reload(campaign);
  assert.equal(evaluate(resumed).ok, true);
  assert.equal(readJsonl(resumed.paths.runs)[0].identity, 'claude-interactive');
});

test('a runtime failure reports the message of the CLI and ingests nothing', async () => {
  const config = tempConfig();
  const campaign = eulerCampaign(config, { formal: false });
  await withEnv({ FAKE_RUNTIME_MODE: 'logged-out' }, async () => {
    const outcome = await step(campaign, { runtime: 'claude' });
    assert.equal(outcome.results[0].status, 'failed');
    assert.match(outcome.results[0].detail, /Not logged in/);
  });
  assert.equal(existsSync(reload(campaign).paths.screenJson), false);
});
