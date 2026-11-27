import { randomBytes } from 'node:crypto';
import { copyFileSync, existsSync, mkdirSync, readdirSync, renameSync } from 'node:fs';
import { join } from 'node:path';
import { proofHash, recordRun, saveCampaign, transition, usageTotals } from './campaign.mjs';
import { writeCapsule } from './capsule.mjs';
import { certificateStatus, certify, provenance } from './certificate.mjs';
import { evaluate, reviewSummary } from './gates.mjs';
import { leanGate } from './lean.mjs';
import { assemblePrompt, loadRole } from './prompt.mjs';
import { pickRuntime, runRuntime } from './runtimes.mjs';
import { HarnessError, appendJsonl, fileSha256, listFiles, now, readJson, readJsonl, rel, writeJson } from './util.mjs';

export const newRunId = role => `${now().replace(/[-:]/g, '').replace(/\..+/, '')}-${role}-${randomBytes(2).toString('hex')}`;

// --- Scope enforcement -----------------------------------------------------------

export function snapshot(dir) {
  return Object.fromEntries(listFiles(dir).map(path => [rel(dir, path), fileSha256(path)]));
}

export function globToRegExp(pattern) {
  let source = '';
  for (let index = 0; index < pattern.length; index++) {
    const char = pattern[index];
    if (char === '*' && pattern[index + 1] === '*') { source += '.*'; index++; if (pattern[index + 1] === '/') index++; }
    else if (char === '*') source += '[^/]*';
    else source += char.replace(/[.+?^${}()|[\]\\]/g, '\\$&');
  }
  return new RegExp(`^${source}$`);
}

// Files the orchestrator itself maintains are never attributed to a model run.
const MANAGED = [/^campaign\.json$/, /^CAPSULE\.md$/, /^runs\.jsonl$/, /^ledger\.jsonl$/, /^outbox\.jsonl$/, /^screen\.json$/, /^directions\.json$/, /^receipts\//, /^reviews\/review-\d+\.(json|md)$/, /^certificates\/[^/]+\/(receipt\.json|output\.txt|provenance\.json)$/];

export function scopeViolations(before, after, writes) {
  const allowed = writes.map(globToRegExp);
  const changed = Object.keys({ ...before, ...after }).filter(path => before[path] !== after[path]);
  return {
    changed,
    violations: changed.filter(path => !MANAGED.some(pattern => pattern.test(path)) && !allowed.some(pattern => pattern.test(path)))
  };
}

// --- Next action ---------------------------------------------------------------------

function directions(campaign) {
  return existsSync(campaign.paths.directions) ? readJson(campaign.paths.directions).directions || [] : [];
}

function certificatesNeedingVerifier(campaign) {
  const root = campaign.paths.certificates;
  if (!existsSync(root)) return [];
  return readdirSync(root, { withFileTypes: true }).filter(entry => entry.isDirectory()).map(entry => entry.name).filter(name => {
    const path = join(root, name, 'manifest.json');
    return existsSync(path) && !readJson(path).verifier?.path;
  });
}

function certificatesNeedingReceipt(campaign) {
  return Object.entries(certificateStatus(campaign))
    .filter(([name, status]) => !status.ok && !/failing/.test(status.reason) && readJson(join(campaign.paths.certificates, name, 'manifest.json'), {}).verifier?.path)
    .map(([name]) => name);
}

const leanConfigured = campaign => Boolean(campaign.state.lean?.claim && campaign.state.lean?.main);

export function nextAction(campaign) {
  const { state, config } = campaign;
  if (['done', 'rejected', 'parked'].includes(state.stage)) return { kind: 'stop', reason: `campaign is ${state.stage}` };
  if (state.pending) return { kind: 'stop', reason: `manual run ${state.pending.runId} is pending: discover record ${campaign.slug} ${state.pending.runId}` };
  const gate = evaluate(campaign);
  switch (state.stage) {
    case 'screen': {
      const screen = existsSync(campaign.paths.screenJson) ? readJson(campaign.paths.screenJson) : null;
      if (screen && !['eligible', 'fixture'].includes(screen.verdict)) return { kind: 'transition', to: 'rejected', reason: `screen verdict ${screen.verdict}: ${screen.priorResolution || screen.summary || ''}`.slice(0, 400) };
      if (gate.ok) return { kind: 'transition', to: 'map', reason: 'screen gate passed; statement frozen', pinStatement: true };
      if (screen && (state.counters.screenRounds || 0) >= 2) return { kind: 'transition', to: 'parked', reason: `screen gate still failing after 2 rounds: ${gate.checks.filter(item => !item.ok).map(item => item.label).join('; ')}` };
      return { kind: 'role', role: 'screen', counter: 'screenRounds', note: screen ? `Previous screen failed these checks: ${gate.checks.filter(item => !item.ok).map(item => item.label).join('; ')}` : null };
    }
    case 'map':
      if (gate.ok) return { kind: 'transition', to: 'explore', reason: 'map gate passed' };
      if ((state.counters.mapRounds || 0) >= 2) return { kind: 'transition', to: 'parked', reason: 'map gate failing after 2 rounds' };
      return { kind: 'role', role: 'map', counter: 'mapRounds', note: existsSync(campaign.paths.directions) ? `Fix: ${gate.checks.filter(item => !item.ok).map(item => item.label).join('; ')}` : null };
    case 'explore': {
      const verify = certificatesNeedingVerifier(campaign);
      if (verify.length) return { kind: 'role', role: 'verify-certificate', tasks: verify.map(certificate => ({ certificate })) };
      const receipts = certificatesNeedingReceipt(campaign);
      if (receipts.length) return { kind: 'certify', names: receipts };
      if (gate.ok) return { kind: 'transition', to: 'draft', reason: 'a promising attempt is ready to draft' };
      const open = directions(campaign).filter(direction => !['killed', 'resolved'].includes(direction.status));
      if (!open.length && (state.counters.remaps || 0) < 1) return { kind: 'remap' };
      if (!open.length || state.counters.exploreRounds >= config.explore.maxRounds) {
        return { kind: 'transition', to: 'parked', reason: open.length ? `no promising attempt after ${state.counters.exploreRounds} explore rounds` : 'every direction is killed, including after one re-map' };
      }
      const batch = open.slice(0, config.explore.parallel);
      return { kind: 'explore', counter: 'exploreRounds', tasks: batch.map((direction, index) => ({ direction, role: direction.kind === 'adversarial' ? 'adversary' : 'explore', index })) };
    }
    case 'draft':
      if (gate.ok) return { kind: 'transition', to: 'review', reason: 'PROOF.md passes the structural gate' };
      if ((state.counters.draftRounds || 0) >= 3) return { kind: 'transition', to: 'parked', reason: 'draft gate failing after 3 rounds' };
      return { kind: 'role', role: 'draft', counter: 'draftRounds', note: existsSync(campaign.paths.proof) ? `PROOF.md fails: ${gate.checks.filter(item => !item.ok).map(item => item.label).slice(0, 8).join('; ')}` : null };
    case 'review': {
      const summary = reviewSummary(campaign);
      if (gate.ok) return { kind: 'transition', to: 'formalize', reason: `${summary.passes.length} independent pass verdicts on PROOF.md ${summary.hash.slice(0, 12)}` };
      const structural = evaluate(campaign, 'draft');
      if (!structural.ok) return { kind: 'transition', to: 'draft', reason: 'PROOF.md no longer passes the draft gate' };
      if (summary.current.some(review => review.verdict === 'fail')) return { kind: 'review-failed' };
      const needsRepair = summary.current.some(review => review.verdict === 'changes') || summary.blocking.length;
      if (needsRepair) {
        if (state.counters.repairRounds >= config.review.maxRepairRounds) return { kind: 'transition', to: 'parked', reason: `defects remain after ${state.counters.repairRounds} repair rounds` };
        return { kind: 'role', role: 'repair', counter: 'repairRounds' };
      }
      const author = state.author?.identity;
      const haveIndependent = summary.current.some(review => (review.identity || review.runtime) !== author);
      return { kind: 'role', role: 'review', runtimeHint: { index: summary.current.length, exclude: haveIndependent ? null : author?.split(':')[0] } };
    }
    case 'formalize': {
      if (!leanConfigured(campaign) && !config.formal.required && !state.formalRequested) return { kind: 'transition', to: 'publish', reason: 'Lean formalization not requested (label stays ⬜ or 🟢 by receipts)' };
      if (state.polarity === 'unknown') return { kind: 'transition', to: 'parked', reason: 'set polarity (prove/disprove) before formalizing: discover set <slug> polarity <value>' };
      const lean = state.lean || {};
      const receipt = existsSync(join(campaign.paths.receipts, 'lean-latest.json')) ? readJson(join(campaign.paths.receipts, 'lean-latest.json')) : null;
      if (!leanConfigured(campaign)) return { kind: 'role', role: 'formalize-statement', counter: 'formalRounds' };
      if (!lean.claimFingerprint) {
        if (!lean.fingerprintForAudit) return { kind: 'lean-gate', purpose: 'fingerprint' };
        if (!lean.audit || lean.audit.fingerprint !== lean.fingerprintForAudit) return { kind: 'role', role: 'audit-statement' };
        if (lean.audit.verdict !== 'faithful') {
          if ((state.counters.formalRounds || 0) >= config.formal.maxRounds) return { kind: 'transition', to: 'parked', reason: `statement audit ${lean.audit.verdict} after ${state.counters.formalRounds} rounds` };
          return { kind: 'role', role: 'formalize-statement', counter: 'formalRounds', note: `The statement auditor returned ${lean.audit.verdict}: ${lean.audit.summary}`, resetFingerprint: true };
        }
        return { kind: 'pin' };
      }
      if (gate.ok) return { kind: 'transition', to: 'publish', reason: 'Lean gate passed against the pinned statement' };
      if (receipt?.status === 'pass' && receipt.expectFingerprint === lean.claimFingerprint && !(lean.lastFormalRun > receipt.finishedAt)) {
        return { kind: 'transition', to: 'parked', reason: `Lean gate passes but the formalize stage gate does not: ${gate.checks.filter(item => !item.ok).map(item => item.label).join('; ')}` };
      }
      if (!receipt || receipt.expectFingerprint !== lean.claimFingerprint || (lean.lastFormalRun && lean.lastFormalRun > receipt.finishedAt)) return { kind: 'lean-gate', purpose: 'verify' };
      if (receipt.status === 'fail' && (state.counters.formalRounds || 0) >= config.formal.maxRounds) {
        return { kind: 'transition', to: 'parked', reason: `Lean gate failing after ${state.counters.formalRounds} rounds: ${(receipt.failures || []).slice(0, 2).join('; ')}` };
      }
      return { kind: 'role', role: 'formalize', counter: 'formalRounds' };
    }
    case 'publish':
      return gate.ok ? { kind: 'transition', to: 'done', reason: 'landing page generated' } : { kind: 'publish' };
    default:
      return { kind: 'stop', reason: `unknown stage ${state.stage}` };
  }
}

// --- Running roles ---------------------------------------------------------------------

function identityOf(runtime, config, as) {
  if (runtime === 'manual') return as || 'manual';
  const model = config.runtimes[runtime]?.model;
  return model ? `${runtime}:${model}` : runtime;
}

function taskFor(campaign, role, base = {}) {
  const task = { ...base };
  if (role === 'explore' || role === 'adversary') {
    const count = readJsonl(campaign.paths.ledger).filter(entry => entry.direction === task.direction.id).length;
    task.attempt = String(count + 1).padStart(3, '0');
  }
  if (role === 'review') task.proofSha256 = proofHash(campaign);
  if (role === 'formalize') {
    const lean = campaign.state.lean;
    task.main = lean.main;
    task.goal = campaign.state.polarity === 'disprove' ? `¬ ${lean.claim}` : lean.claim;
  }
  return task;
}

async function runRole(campaign, role, task, { runtime, as, snapshotBefore = null } = {}) {
  const runId = newRunId(role);
  const runDir = join(campaign.paths.runDir, runId);
  mkdirSync(runDir, { recursive: true });
  const capsule = writeCapsule(campaign);
  const prompt = assemblePrompt({ campaign, role, task, runId, capsule });
  const { schemaPath, meta } = loadRole(role);
  const before = snapshotBefore || snapshot(campaign.dir);
  writeJson(join(runDir, 'meta.json'), { runId, role, runtime, task: { ...task, direction: task.direction?.id }, writes: prompt.writes, promptTokens: prompt.tokens, pack: { used: prompt.pack.used, omitted: prompt.pack.omitted.map(item => item.locator) } });
  const result = await runRuntime({ config: campaign.config, runtime, prompt: prompt.text, schemaPath, role: meta, runDir, cwd: campaign.config.root });
  const identity = identityOf(runtime, campaign.config, as);
  if (result.status === 'pending') {
    writeJson(join(runDir, 'before.json'), before);
    campaign.state.pending = { runId, role, runtime, identity, task: { ...task, direction: task.direction?.id }, writes: prompt.writes, at: now() };
    saveCampaign(campaign);
    return { runId, status: 'pending', detail: result.detail, promptTokens: prompt.tokens };
  }
  return finishRun(campaign, { runId, role, runtime, identity, task, writes: prompt.writes, before, result, promptTokens: prompt.tokens, checkScope: !snapshotBefore });
}

export function finishRun(campaign, { runId, role, runtime, identity, task, writes, before, result, promptTokens, checkScope = true }) {
  const after = snapshot(campaign.dir);
  const scope = checkScope ? scopeViolations(before, after, writes) : { changed: [], violations: [] };
  let status = result.status;
  let detail = result.detail || null;
  if (scope.violations.length) { status = 'scope-violation'; detail = `wrote outside its scope: ${scope.violations.join(', ')}`; }
  let ingest = null;
  if (status === 'ok') {
    try { ingest = ingestAnswer(campaign, role, task, result.answer, { runId, runtime, identity }); }
    catch (error) { status = 'ingest-failed'; detail = error.message; }
  }
  recordRun(campaign, { runId, role, runtime, identity, status, detail, usage: result.usage || null, usageWarning: result.warning || null, promptTokens, changed: scope.changed.length, summary: result.answer?.summary?.slice(0, 500) || null });
  writeCapsule(campaign);
  return { runId, role, runtime, status, detail, ingest, summary: result.answer?.summary };
}

// Turn a validated structured answer into campaign records. Facts the
// orchestrator can establish itself (hashes, run IDs, identities) are never
// taken from the model.
export function ingestAnswer(campaign, role, task, answer, { runId, runtime, identity }) {
  const { paths, state } = campaign;
  switch (role) {
    case 'screen':
      writeJson(paths.screenJson, { ...answer, runId, identity, statementSha256: existsSync(paths.statement) ? fileSha256(paths.statement) : null, at: now() });
      return { verdict: answer.verdict };
    case 'map': {
      const seen = new Set();
      const list = (answer.directions || []).map((direction, index) => {
        let id = /^D\d+$/.test(direction.id) ? direction.id : `D${index + 1}`;
        while (seen.has(id)) id = `D${seen.size + 1}`;
        seen.add(id);
        return { ...direction, id, status: 'open' };
      });
      writeJson(paths.directions, { runId, at: now(), directions: list });
      return { directions: list.length };
    }
    case 'explore':
    case 'adversary': {
      const attempt = answer.attempt || {};
      const id = `${task.direction.id}-${task.attempt}`;
      const file = answer.file || `attempts/${task.direction.id}/${task.attempt}.md`;
      if (!existsSync(join(campaign.dir, file))) throw new HarnessError('MISSING_ARTIFACT', `the attempt report ${file} was not written`);
      appendJsonl(paths.ledger, { id, direction: task.direction.id, role, runtime, identity, runId, at: now(), file, ...attempt });
      const list = readJson(paths.directions);
      const direction = list.directions.find(item => item.id === task.direction.id);
      if (direction) {
        direction.status = attempt.status === 'killed' ? 'killed' : attempt.status === 'resolved' ? 'resolved' : attempt.status === 'promising' ? 'promising' : 'open';
        direction.attempts = (direction.attempts || 0) + 1;
        writeJson(paths.directions, list);
      }
      for (const name of answer.certificates || []) {
        if (existsSync(join(paths.certificates, name, 'manifest.json'))) provenance(campaign, name, { generator: { runId, identity, attempt: id } });
      }
      return { attempt: id, outcome: attempt.outcome, status: attempt.status };
    }
    case 'verify-certificate': {
      const manifest = readJson(join(paths.certificates, task.certificate, 'manifest.json'));
      if (!manifest.verifier?.path || !existsSync(join(paths.certificates, task.certificate, manifest.verifier.path))) throw new HarnessError('MISSING_ARTIFACT', 'the verifier was not written or manifest.verifier.path is unset');
      provenance(campaign, task.certificate, { verifier: { runId, identity } });
      return { certificate: task.certificate, verdict: answer.verdict };
    }
    case 'draft':
      state.author = { identity, runId };
      if (state.polarity === 'unknown' && ['prove', 'disprove'].includes(answer.polarity)) state.polarity = answer.polarity;
      saveCampaign(campaign);
      return { kind: answer.kind, main: answer.main };
    case 'review': {
      mkdirSync(paths.reviews, { recursive: true });
      const number = readdirSync(paths.reviews).filter(name => /^review-\d+\.json$/.test(name)).length + 1;
      const name = `review-${String(number).padStart(3, '0')}`;
      const report = join(paths.reviews, `work-${runId}`, 'REPORT.md');
      if (existsSync(report)) copyFileSync(report, join(paths.reviews, `${name}.md`));
      const summary = reviewSummary(campaign);
      // Bind the verdict to the hashes the reviewer was given, not to the answer.
      const stillCurrent = proofHash(campaign) === task.proofSha256;
      writeJson(join(paths.reviews, `${name}.json`), {
        ...answer, proofSha256: stillCurrent ? task.proofSha256 : `stale:${task.proofSha256}`,
        certificatesKey: summary.certificateHash, runtime, identity, runId, at: now(), report: existsSync(report) ? `${name}.md` : null
      });
      return { review: name, verdict: answer.verdict, defects: (answer.defects || []).length };
    }
    case 'repair':
      return { addressed: (answer.addressed || []).length, kind: answer.kind };
    case 'formalize-statement': {
      const lean = answer.lean || {};
      if (!lean.claim || !lean.main || !lean.modules?.length) throw new HarnessError('BAD_ANSWER', 'formalize-statement must return lean.modules, lean.claim and lean.main');
      state.lean = { project: 'formal', modules: lean.modules, claim: lean.claim, main: lean.main, namespace: lean.namespace || lean.claim.split('.').slice(0, -1).join('.'), endpoints: lean.endpoints || [], claimFingerprint: null, fingerprintForAudit: null, audit: null };
      saveCampaign(campaign);
      return { lean: state.lean.claim };
    }
    case 'audit-statement':
      state.lean.audit = { verdict: answer.verdict, summary: (answer.summary || '').slice(0, 600), fingerprint: state.lean.fingerprintForAudit, identity, runId, at: now() };
      saveCampaign(campaign);
      return { verdict: answer.verdict };
    case 'formalize':
      state.lean.lastFormalRun = now();
      saveCampaign(campaign);
      return { status: answer.status };
    default:
      return null;
  }
}

export async function campaignLeanGate(campaign, { requirePinned = true } = {}) {
  const lean = campaign.state.lean;
  if (!leanConfigured(campaign)) throw new HarnessError('NO_LEAN', 'campaign.lean is not configured (run the formalize-statement step)');
  const { config } = campaign;
  return leanGate({
    project: join(campaign.dir, lean.project), modules: lean.modules, claim: lean.claim, main: lean.main,
    polarity: campaign.state.polarity, namespace: lean.namespace, endpoints: lean.endpoints || [],
    allowedAxioms: config.formal.allowedAxioms, expectFingerprint: lean.claimFingerprint || null, requirePinned,
    buildTimeoutMs: config.formal.buildTimeoutMinutes * 60_000, replay: config.formal.replay,
    receiptsDir: campaign.paths.receipts, errorTailLines: config.context.errorTailLines
  });
}

// --- One step and the loop ------------------------------------------------------------------

export async function step(campaign, { runtime: override = null, as = null, dryRun = false } = {}) {
  const action = nextAction(campaign);
  if (dryRun || action.kind === 'stop') return { action };
  const { state, config } = campaign;
  const bump = name => { if (name) { state.counters[name] = (state.counters[name] || 0) + 1; saveCampaign(campaign); } };
  const usage = usageTotals(campaign);
  if (['role', 'explore'].includes(action.kind) && usage.tokens >= state.budget.tokens) {
    return { action: { kind: 'stop', reason: `token budget spent (${usage.tokens} of ${state.budget.tokens})` } };
  }
  switch (action.kind) {
    case 'transition': {
      if (action.pinStatement) state.statementSha256 = fileSha256(campaign.paths.statement);
      transition(campaign, action.to, action.reason);
      appendJsonl(campaign.paths.outbox, { at: now(), type: 'transition', to: action.to, reason: action.reason });
      writeCapsule(campaign);
      return { action };
    }
    case 'role': {
      bump(action.counter);
      if (action.resetFingerprint && state.lean) { state.lean.fingerprintForAudit = null; saveCampaign(campaign); }
      const tasks = action.tasks || [{}];
      const results = [];
      for (const [index, base] of tasks.entries()) {
        const generator = base.certificate ? provenance(campaign, base.certificate).generator?.identity?.split(':')[0] : null;
        const hint = action.runtimeHint || {};
        const runtime = pickRuntime(config, action.role, { index: hint.index ?? index, exclude: hint.exclude || generator, override });
        const task = taskFor(campaign, action.role, { ...base, note: action.note || null });
        results.push(await runRole(campaign, action.role, task, { runtime, as }));
        if (results.at(-1).status === 'pending') break;
      }
      return { action, results };
    }
    case 'explore': {
      bump(action.counter);
      const before = snapshot(campaign.dir);
      const runs = action.tasks.map(item => {
        const runtime = pickRuntime(config, 'explore', { index: item.index, override });
        const task = taskFor(campaign, item.role, { direction: item.direction });
        return { item, runtime, task };
      });
      if (runs.some(run => run.runtime === 'manual')) {
        return { action, results: [await runRole(campaign, runs[0].item.role, runs[0].task, { runtime: 'manual', as })] };
      }
      const results = await Promise.all(runs.map(run => runRole(campaign, run.item.role, run.task, { runtime: run.runtime, as, snapshotBefore: before })));
      // Concurrent workers: check the union of their scopes once, after all finish.
      const writes = runs.flatMap(run => loadRole(run.item.role).meta.writes.map(pattern => pattern.replaceAll('{{DIRECTION}}', run.task.direction.id)));
      const scope = scopeViolations(before, snapshot(campaign.dir), writes);
      if (scope.violations.length) appendJsonl(campaign.paths.runs, { at: now(), runId: 'explore-batch', status: 'scope-violation', detail: scope.violations.join(', ') });
      return { action, results, scope };
    }
    case 'certify': {
      const receipts = [];
      for (const name of action.names) receipts.push({ name, status: (await certify(campaign, name, null)).status });
      writeCapsule(campaign);
      return { action, receipts };
    }
    case 'lean-gate': {
      const receipt = await campaignLeanGate(campaign, { requirePinned: action.purpose === 'verify' });
      if (action.purpose === 'fingerprint') {
        if (!receipt.statement?.fingerprint) {
          const retry = await retryStatement(campaign, receipt, override, as);
          return { action, receipt: { status: receipt.status, failures: receipt.failures.slice(0, 5) }, ...retry };
        }
        state.lean.fingerprintForAudit = receipt.statement.fingerprint;
        saveCampaign(campaign);
      }
      writeCapsule(campaign);
      return { action, receipt: { status: receipt.status, failures: receipt.failures.slice(0, 5), fingerprint: receipt.statement?.fingerprint } };
    }
    case 'pin': {
      const receipt = await campaignLeanGate(campaign, { requirePinned: false });
      if (receipt.statement?.fingerprint !== state.lean.audit.fingerprint) {
        state.lean.fingerprintForAudit = null;
        saveCampaign(campaign);
        return { action, pinned: false, reason: 'the statement changed after the audit; it will be re-audited' };
      }
      state.lean.claimFingerprint = receipt.statement.fingerprint;
      state.lean.toolchain = receipt.toolchain;
      saveCampaign(campaign);
      appendJsonl(campaign.paths.outbox, { at: now(), type: 'statement-pinned', fingerprint: state.lean.claimFingerprint });
      writeCapsule(campaign);
      return { action, pinned: true, fingerprint: state.lean.claimFingerprint };
    }
    case 'review-failed': {
      const summary = reviewSummary(campaign);
      const failed = summary.current.filter(review => review.verdict === 'fail');
      const ledger = readJsonl(campaign.paths.ledger);
      const routes = ledger.filter(item => item.status === 'promising');
      for (const entry of routes) {
        appendJsonl(campaign.paths.ledger, { ...entry, id: `${entry.id}-x`, supersedes: entry.id, status: 'killed', outcome: 'failed', failurePoint: `review ${failed[0].file}: ${(failed[0].summary || '').slice(0, 240)}`, at: now() });
      }
      const list = readJson(campaign.paths.directions);
      for (const direction of list.directions) if (routes.some(entry => entry.direction === direction.id)) direction.status = 'killed';
      writeJson(campaign.paths.directions, list);
      transition(campaign, 'explore', `review verdict fail (${failed.map(review => review.file).join(', ')}); route marked killed`);
      writeCapsule(campaign);
      return { action };
    }
    case 'remap': {
      state.counters.remaps = (state.counters.remaps || 0) + 1;
      renameSync(campaign.paths.directions, join(campaign.dir, `directions-${state.counters.remaps}.json`));
      transition(campaign, 'map', 'every direction was killed; re-map with the failures in view');
      writeCapsule(campaign);
      return { action };
    }
    case 'publish': {
      const { publish } = await import('./publish.mjs');
      const result = publish(campaign);
      writeCapsule(campaign);
      return { action, result };
    }
    default:
      return { action };
  }
}

async function retryStatement(campaign, receipt, override, as) {
  const { state, config } = campaign;
  if ((state.counters.formalRounds || 0) >= config.formal.maxRounds) {
    transition(campaign, 'parked', `Lean statement does not build after ${state.counters.formalRounds} rounds: ${receipt.failures.slice(0, 2).join('; ')}`);
    return { parked: true };
  }
  state.counters.formalRounds = (state.counters.formalRounds || 0) + 1;
  saveCampaign(campaign);
  const runtime = pickRuntime(config, 'formalize-statement', { override });
  const task = taskFor(campaign, 'formalize-statement', { note: `The statement project failed the gate: ${receipt.failures.slice(0, 3).join('; ')}` });
  return { retry: await runRole(campaign, 'formalize-statement', task, { runtime, as }) };
}

export async function runLoop(campaign, { until = null, maxSteps = 50, runtime = null, as = null, onStep = () => {} } = {}) {
  const log = [];
  for (let count = 0; count < maxSteps; count++) {
    if (until && campaign.state.stage === until) { log.push({ stop: `reached ${until}` }); break; }
    const outcome = await step(campaign, { runtime, as });
    log.push(outcome);
    onStep(outcome);
    if (outcome.action.kind === 'stop') break;
    const failed = (outcome.results || []).filter(result => !['ok', 'pending'].includes(result.status));
    if ((outcome.results || []).some(result => result.status === 'pending')) break;
    if (failed.length && failed.length === (outcome.results || []).length) {
      const recent = readJsonl(campaign.paths.runs).slice(-3);
      if (recent.length === 3 && recent.every(run => run.status !== 'ok')) { log.push({ stop: 'three consecutive failed runs; inspect runs.jsonl' }); break; }
    }
  }
  return log;
}

// Complete a manual (interactive) run from its answer file.
export function recordManual(campaign, runId, answerPath = null) {
  const pending = campaign.state.pending;
  if (!pending || pending.runId !== runId) throw new HarnessError('NO_PENDING', `No pending run ${runId}`);
  const runDir = join(campaign.paths.runDir, runId);
  const path = answerPath || join(runDir, 'answer.json');
  const answer = readJson(path);
  const before = readJson(join(runDir, 'before.json'));
  const task = { ...pending.task };
  if (task.direction) task.direction = directions(campaign).find(direction => direction.id === task.direction);
  delete campaign.state.pending;
  saveCampaign(campaign);
  return finishRun(campaign, { runId, role: pending.role, runtime: 'manual', identity: pending.identity, task, writes: pending.writes, before, result: { status: 'ok', answer, usage: null } });
}

