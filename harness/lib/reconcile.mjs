import { existsSync } from 'node:fs';
import { join } from 'node:path';
import { STAGES, listCampaigns, loadCampaign, usageTotals } from './campaign.mjs';
import { certificateStatus } from './certificate.mjs';
import { artifactProblems, evaluate, fingerprintState, reviewSummary, verificationLabel } from './gates.mjs';
import { leanReceiptCurrent } from './lean.mjs';
import { now, readJson, writeText } from './util.mjs';

const STALE_PENDING_HOURS = 4;

// A reconciliation receipt answers: what is each campaign's real state, what
// is finished but unpublished, and what is stale or inconsistent. It reads
// the record, never the capsule's own claims.
export function reconcile(config) {
  const rows = [], blockers = [];
  for (const slug of listCampaigns(config)) {
    const campaign = loadCampaign(config, slug);
    const { state } = campaign;
    const gate = evaluate(campaign);
    const review = reviewSummary(campaign);
    const prints = fingerprintState(campaign);
    const usage = usageTotals(campaign);
    const label = verificationLabel(campaign);
    const issues = [];
    if (prints.statementChanged) issues.push('STATEMENT.md changed after screening: re-screen before any further work');
    if (state.pending) {
      const hours = (Date.now() - Date.parse(state.pending.at)) / 3_600_000;
      if (hours > STALE_PENDING_HOURS) issues.push(`manual run ${state.pending.runId} pending for ${hours.toFixed(1)} h: record it or clear it with discover set ${slug} pending none`);
    }
    if (review.stale && STAGES.indexOf(state.stage) >= STAGES.indexOf('review')) issues.push(`${review.stale} review(s) are bound to an older PROOF.md`);
    const reviewGate = evaluate(campaign, 'review');
    const published = existsSync(join(config.root, config.resultsRoot, slug, 'README.md'));
    if (reviewGate.ok && !published) issues.push('reviewed but unpublished: finish formalize/publish before starting a new target');
    if (published && !reviewGate.ok) issues.push('published page exists but the review gate no longer passes: regenerate or withdraw');
    for (const [name, status] of Object.entries(certificateStatus(campaign))) if (!status.ok) issues.push(`certificate ${name}: ${status.reason}`);
    const receiptPath = join(campaign.paths.receipts, 'lean-latest.json');
    if (state.lean?.claimFingerprint && existsSync(receiptPath) && !leanReceiptCurrent(join(campaign.dir, state.lean.project), readJson(receiptPath))) issues.push('Lean receipt is not current for the sources on disk');
    if (usage.tokens > state.budget.tokens) issues.push(`over budget: ${usage.tokens} of ${state.budget.tokens} tokens`);
    issues.push(...artifactProblems(campaign));
    rows.push({ slug, stage: state.stage, status: state.status, gate: gate.ok, label: label.emoji, passes: review.passes.length, tokens: usage.tokens, issues, fixture: state.fixture });
    for (const issue of issues) blockers.push(`${slug}: ${issue}`);
  }
  return { at: now(), rows, blockers };
}

export function reconcileMarkdown(report) {
  const lines = [
    `# Reconciliation receipt ${report.at}`,
    '',
    `${report.rows.length} campaigns inspected; ${report.blockers.length} open issues.`,
    '',
    '| Campaign | Stage | Gate | Label | Passing reviews | Tokens | Issues |',
    '|---|---|---|---|---|---|---|',
    ...report.rows.map(row => `| ${row.slug}${row.fixture ? ' (fixture)' : ''} | ${row.stage} | ${row.gate ? 'pass' : 'open'} | ${row.label} | ${row.passes} | ${row.tokens.toLocaleString('en-US')} | ${row.issues.length} |`),
    '',
    '## Issues',
    '',
    ...(report.blockers.length ? report.blockers.map(item => `- ${item}`) : ['None.']),
    ''
  ];
  return lines.join('\n');
}

export function writeReconcile(config, report) {
  const path = join(config.root, config.campaignRoot, '_reconcile', `${report.at.slice(0, 10)}.md`);
  writeText(path, reconcileMarkdown(report));
  return path;
}
