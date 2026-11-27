import { existsSync, readdirSync } from 'node:fs';
import { join } from 'node:path';
import { proofHash, statementHash, usageTotals } from './campaign.mjs';
import { certificateStatus } from './certificate.mjs';
import { claimProblems, parseProof } from './claims.mjs';
import { leanReceiptCurrent } from './lean.mjs';
import { outline, extract } from './sections.mjs';
import { fileSize, listFiles, readJson, readJsonl, readText, rel } from './util.mjs';

// Each gate is deterministic: it reads files and receipts, never a model's
// opinion of itself. A failing check names the artifact and the fix.
const check = (ok, label, fix = '') => ({ ok: Boolean(ok), label, fix });

export function reviews(campaign) {
  const dir = campaign.paths.reviews;
  if (!existsSync(dir)) return [];
  return readdirSync(dir).filter(name => /^review-\d+\.json$/.test(name)).sort()
    .map(name => ({ file: name, ...readJson(join(dir, name)) }));
}

export function reviewSummary(campaign) {
  const hash = proofHash(campaign);
  const certificates = certificateStatus(campaign);
  const certificateHash = JSON.stringify(Object.fromEntries(Object.entries(certificates).map(([name, value]) => [name, value.receipt?.files || null])));
  const all = reviews(campaign);
  const current = all.filter(review => review.proofSha256 === hash && (review.certificatesKey ?? certificateHash) === certificateHash);
  const passes = current.filter(review => review.verdict === 'pass');
  const blocking = current.flatMap(review => (review.defects || []).filter(defect => defect.severity === 'blocking').map(defect => ({ review: review.file, ...defect })));
  const identities = new Set(passes.map(review => review.identity || review.runtime));
  return { hash, all, current, stale: all.length - current.length, passes, blocking, identities, certificateHash };
}

function statementChecks(campaign) {
  if (!existsSync(campaign.paths.statement)) return [check(false, 'STATEMENT.md exists', 'discover new created it; restore it')];
  const text = readText(campaign.paths.statement);
  const anchors = outline(text).map(entry => entry.anchor);
  const section = anchor => (anchors.includes(anchor) ? extract(text, { anchor }) : '');
  return [
    check(anchors.includes('exact-statement') && !/TODO/.test(section('exact-statement')), 'STATEMENT.md has a filled "## Exact statement"', 'Quote the statement exactly, with its quantifiers'),
    check(anchors.includes('source') && /(https?:\/\/|arXiv:|doi:|\b10\.\d{4,}\/)/i.test(section('source')) && !/TODO/.test(section('source')), 'STATEMENT.md "## Source" cites a primary source (URL, arXiv or DOI) with statement number', 'Add the version-pinned source and statement number'),
    check(anchors.includes('scope'), 'STATEMENT.md has "## Scope" (what counts as full vs partial)', 'Add the scope section')
  ];
}

// Attempts that still count: not superseded by a later record and not on a killed direction.
export function liveLedger(campaign) {
  const entries = readJsonl(campaign.paths.ledger);
  const superseded = new Set(entries.map(entry => entry.supersedes).filter(Boolean));
  const killed = new Set((existsSync(campaign.paths.directions) ? readJson(campaign.paths.directions).directions || [] : []).filter(direction => direction.status === 'killed').map(direction => direction.id));
  return entries.filter(entry => !superseded.has(entry.id) && !entry.supersedes && !killed.has(entry.direction));
}

export const GATES = {
  screen(campaign) {
    const config = campaign.config.screen;
    const out = statementChecks(campaign);
    const screen = existsSync(campaign.paths.screenJson) ? readJson(campaign.paths.screenJson) : null;
    out.push(check(screen, 'screen.json exists (screen step ran)', 'discover step <slug>'));
    if (!screen) return out;
    const allowed = campaign.state.fixture ? ['eligible', 'fixture'] : ['eligible'];
    out.push(check(allowed.includes(screen.verdict), `screen verdict is eligible (got ${screen.verdict})`, 'A prior resolution or low-impact target is rejected, not worked on'));
    out.push(check((screen.priorResolutionSearch?.queries || []).length >= config.minQueries, `≥${config.minQueries} recorded prior-resolution queries`, 'Record every query, including author-name and statement-number searches'));
    out.push(check((screen.priorResolutionSearch?.sources || []).filter(source => /^https?:\/\//.test(source.url || '')).length >= config.minSources, `≥${config.minSources} inspected sources with URLs`, 'List each source read and what it says about the exact claim'));
    out.push(check((screen.significance || '').length >= 80, 'significance case is argued', 'Say what a full resolution would change'));
    out.push(check((screen.attacks || []).length >= 1, 'at least one plausible attack', 'Name a concrete first attack'));
    out.push(check(screen.fullVsPartial, 'full vs partial resolution distinguished'));
    out.push(check(existsSync(campaign.paths.screen), 'SCREEN.md written'));
    return out;
  },
  map(campaign) {
    const directions = existsSync(campaign.paths.directions) ? readJson(campaign.paths.directions).directions || [] : [];
    const lenses = new Set(directions.map(direction => direction.lens));
    return [
      check(existsSync(campaign.paths.map), 'MAP.md written'),
      check(directions.length >= campaign.config.explore.minDirections, `≥${campaign.config.explore.minDirections} directions in directions.json`),
      check(lenses.size >= Math.min(2, directions.length), 'directions use different lenses (different failure modes)'),
      check(directions.every(direction => direction.falsificationTest && direction.milestone && direction.budgetTokens), 'every direction has a falsification test, milestone and token budget'),
      check(directions.some(direction => direction.kind === 'adversarial'), 'one direction is adversarial (counterexample search / stress test)', 'Add a direction of kind "adversarial"')
    ];
  },
  explore(campaign) {
    const ledger = liveLedger(campaign);
    const missing = ledger.filter(entry => !(entry.claim && entry.method && entry.outcome && entry.next));
    const decisive = ledger.filter(entry => ['proof-sketch', 'counterexample', 'reduction', 'partial'].includes(entry.outcome) && ['promising', 'resolved'].includes(entry.status));
    return [
      check(ledger.length > 0, 'at least one attempt recorded in ledger.jsonl'),
      check(missing.length === 0, 'every attempt records claim, method, outcome and next action', `fix: ${missing.map(entry => entry.id).join(', ')}`),
      check(decisive.length > 0, 'an attempt is promising enough to draft (proof sketch, counterexample, reduction or partial theorem)', 'Keep exploring or park the campaign')
    ];
  },
  draft(campaign) {
    const proof = parseProof(campaign.paths.proof);
    const problems = claimProblems(proof, { certificates: certificateStatus(campaign) });
    return [
      check(proof, 'PROOF.md exists'),
      ...problems.map(problem => check(false, problem)),
      check(proof && /<!--\s*digest\s*-->/.test(proof.text), 'PROOF.md opens with a <!-- digest --> abstract', 'Add a ≤10-line digest block'),
      check(fileSize(campaign.paths.proof) <= campaign.config.context.maxArtifactBytes, 'PROOF.md is under the artifact size limit', 'Split appendices into linked files')
    ];
  },
  review(campaign) {
    const config = campaign.config.review;
    const summary = reviewSummary(campaign);
    const author = campaign.state.author?.identity;
    const independent = !config.requireDifferentRuntime || [...summary.identities].some(identity => identity !== author);
    return [
      ...GATES.draft(campaign),
      check(summary.passes.length >= config.passes, `≥${config.passes} pass verdicts bound to the current PROOF.md (have ${summary.passes.length}; ${summary.stale} stale)`),
      check(summary.blocking.length === 0, 'no blocking defect open on the current PROOF.md', summary.blocking.map(defect => `${defect.review}:${defect.id}`).join(', ')),
      check(independent, `a passing review came from an identity other than the author's (${author || 'unknown'})`, 'Route review to the other runtime, or record a manual review with --as <identity>')
    ];
  },
  formalize(campaign) {
    const lean = campaign.state.lean;
    if (!lean) return [check(!campaign.config.formal.required, 'formalization not required by config', 'Configure campaign.lean and run the formal steps')];
    const receipt = existsSync(join(campaign.paths.receipts, 'lean-latest.json')) ? readJson(join(campaign.paths.receipts, 'lean-latest.json')) : null;
    const project = join(campaign.dir, lean.project);
    const statementReview = join(campaign.dir, 'STATEMENT_REVIEW.md');
    return [
      check(existsSync(statementReview), 'STATEMENT_REVIEW.md maps each source condition to its Lean representation'),
      check(lean.claimFingerprint, 'statement fingerprint pinned after audit', 'discover pin-statement <slug>'),
      check(receipt?.status === 'pass', `Lean gate receipt passes${receipt ? ` (${(receipt.failures || []).slice(0, 2).join('; ')})` : ''}`, 'discover lean-gate <slug>'),
      check(receipt && leanReceiptCurrent(project, receipt), 'Lean receipt is current (audited sources unchanged)'),
      check(receipt?.statement?.fingerprint === lean.claimFingerprint, 'receipt fingerprint equals the pinned one')
    ];
  },
  publish(campaign) {
    const target = join(campaign.config.root, campaign.config.resultsRoot, campaign.slug, 'README.md');
    return [check(existsSync(target), `results/${campaign.slug}/README.md generated`, 'discover publish <slug>')];
  }
};

export function evaluate(campaign, stage = campaign.state.stage) {
  const gate = GATES[stage];
  const checks = gate ? gate(campaign) : [];
  return { stage, ok: gate ? checks.every(item => item.ok) : stage === 'done', checks };
}

// The verification label is computed from receipts, never chosen by a model.
export function verificationLabel(campaign) {
  const proof = parseProof(campaign.paths.proof);
  const certificates = Object.values(certificateStatus(campaign));
  const computed = certificates.some(status => status.ok);
  const receiptPath = join(campaign.paths.receipts, 'lean-latest.json');
  const receipt = existsSync(receiptPath) ? readJson(receiptPath) : null;
  const lean = campaign.state.lean;
  const current = receipt && lean && leanReceiptCurrent(join(campaign.dir, lean.project), receipt);
  const full = current && receipt.status === 'pass' && receipt.statement?.fingerprint === lean.claimFingerprint && proof?.kind === 'complete';
  const partial = current && receipt.status === 'pass' && !full;
  if (full) return { emoji: '✅', text: 'Main claim fully Lean-verified against the audited statement' };
  if (partial && computed) return { emoji: '🟢', text: 'Partial Lean verification with a computer-assisted proof' };
  if (partial) return { emoji: '🟡', text: 'Partial Lean verification' };
  return { emoji: '⬜', text: computed ? 'Written proof with exact computer checks; no Lean verification' : 'Written proof; no Lean verification' };
}

export function artifactProblems(campaign) {
  const limit = campaign.config.context.maxArtifactBytes;
  return listFiles(campaign.dir).filter(path => /\.(md|json)$/.test(path) && fileSize(path) > limit)
    .map(path => `${rel(campaign.dir, path)} is ${fileSize(path)} bytes (limit ${limit}); split it and keep old headings as links`);
}

export function fingerprintState(campaign) {
  return {
    statementSha256: statementHash(campaign),
    pinned: campaign.state.statementSha256,
    statementChanged: Boolean(campaign.state.statementSha256 && statementHash(campaign) !== campaign.state.statementSha256),
    proofSha256: proofHash(campaign),
    usage: usageTotals(campaign)
  };
}

