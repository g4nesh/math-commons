import { existsSync } from 'node:fs';
import { readText } from './util.mjs';
import { outline, tableRows, extract } from './sections.mjs';

// The claim taxonomy decides what may be counted as a result. Only `complete`
// is a resolution of the target; everything else is published with its label.
export const RESULT_KINDS = {
  complete: 'Complete resolution of the exact target',
  partial: 'Precise partial theorem toward the target',
  candidate: 'Unreviewed candidate (not a result)',
  'prior-resolution': 'Target already resolved in the literature (withdraw)',
  consequence: 'Consequence of known theorems (not a new result)',
  obstruction: 'Obstruction to a proof route (not a counterexample to the target)'
};

// Status of each numbered claim in PROOF.md. `assumed` and `cited` must name
// their source; `computed` must name a certificate with a valid receipt.
export const CLAIM_STATUSES = ['proved', 'cited', 'computed', 'assumed', 'open'];

export function parseProof(path) {
  if (!existsSync(path)) return null;
  const text = readText(path);
  const headings = outline(text);
  const claimsHeading = headings.find(entry => entry.anchor === 'claims');
  const rows = claimsHeading ? tableRows(extract(text, { anchor: 'claims' })) : [];
  const claims = rows.map(row => ({
    id: (row.id || '').replace(/[*`]/g, '').trim(),
    statement: row.statement || row.claim || '',
    status: (row.status || '').toLowerCase().replace(/[*`]/g, '').trim(),
    depends: (row['depends on'] || row.depends || '').split(/[,\s]+/).map(item => item.trim()).filter(item => item && item !== '—' && item !== '-'),
    evidence: row.evidence || row.source || ''
  }));
  const kindLine = /^\*\*Result kind:\*\*\s*`?([a-z-]+)`?/im.exec(text);
  const mainLine = /^\*\*Main claim:\*\*\s*`?([A-Za-z0-9_.-]+)`?/im.exec(text);
  return { text, headings, claims, kind: kindLine?.[1] || null, main: mainLine?.[1] || null };
}

export function claimProblems(proof, { certificates = {} } = {}) {
  const problems = [];
  if (!proof) return ['PROOF.md is missing'];
  if (!proof.claims.length) problems.push('PROOF.md has no "## Claims" table (columns: ID | Statement | Status | Depends on | Evidence)');
  if (!proof.kind || !RESULT_KINDS[proof.kind]) problems.push(`PROOF.md must declare **Result kind:** one of ${Object.keys(RESULT_KINDS).join(', ')}`);
  if (!proof.main) problems.push('PROOF.md must declare **Main claim:** <claim ID>');
  const ids = new Set(proof.claims.map(claim => claim.id));
  if (proof.main && !ids.has(proof.main)) problems.push(`Main claim ${proof.main} is not a row of the Claims table`);
  for (const claim of proof.claims) {
    if (!claim.id) { problems.push('A Claims row has no ID'); continue; }
    if (!CLAIM_STATUSES.includes(claim.status)) problems.push(`${claim.id}: status "${claim.status}" is not one of ${CLAIM_STATUSES.join(', ')}`);
    for (const dependency of claim.depends) if (!ids.has(dependency)) problems.push(`${claim.id}: depends on unknown claim ${dependency}`);
    if (['cited', 'assumed'].includes(claim.status) && !claim.evidence.trim()) problems.push(`${claim.id}: ${claim.status} claims must name their source`);
    if (claim.status === 'computed') {
      const names = [...claim.evidence.matchAll(/certificates\/([a-z0-9-]+)/gi)].map(match => match[1]);
      if (!names.length) problems.push(`${claim.id}: computed claims must cite certificates/<name>`);
      for (const name of names) {
        const status = certificates[name];
        if (!status) problems.push(`${claim.id}: certificate ${name} has no receipt`);
        else if (!status.ok) problems.push(`${claim.id}: certificate ${name} receipt is ${status.reason}`);
      }
    }
    const anchor = claim.id.toLowerCase();
    if (!proof.headings.some(heading => heading.anchor === anchor || heading.anchor.startsWith(`${anchor}-`))) {
      problems.push(`${claim.id}: no "### ${claim.id}" section holding its proof or citation`);
    }
  }
  // A claim cannot be proved from something unproved without saying so.
  const byId = Object.fromEntries(proof.claims.map(claim => [claim.id, claim]));
  const open = new Set();
  const visit = (id, trail = new Set()) => {
    if (trail.has(id)) { problems.push(`Circular dependency through ${id}`); return false; }
    const claim = byId[id];
    if (!claim) return false;
    if (claim.status === 'open') return true;
    const next = new Set(trail).add(id);
    return claim.depends.some(dependency => visit(dependency, next));
  };
  if (proof.main && byId[proof.main] && visit(proof.main)) open.add(proof.main);
  if (open.size && proof.kind === 'complete') problems.push('Result kind is complete but the main claim depends on an open claim');
  const assumed = proof.claims.filter(claim => claim.status === 'assumed').map(claim => claim.id);
  if (assumed.length && proof.kind === 'complete') problems.push(`Result kind is complete but claims ${assumed.join(', ')} are assumed; label it partial (conditional) or prove them`);
  return problems;
}
