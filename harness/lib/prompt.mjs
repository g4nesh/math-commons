import { existsSync, readdirSync } from 'node:fs';
import { join } from 'node:path';
import { certificateStatus } from './certificate.mjs';
import { liveLedger, reviewSummary } from './gates.mjs';
import { digestOf, extract, outlineText, parseLocator } from './sections.mjs';
import { HARNESS_ROOT, HarnessError, estimateTokens, listFiles, readJson, readJsonl, readText, rel, truncate } from './util.mjs';

export const ROLES = ['scout', 'screen', 'map', 'explore', 'adversary', 'verify-certificate', 'draft', 'review', 'repair', 'formalize-statement', 'audit-statement', 'formalize'];

export function loadRole(role) {
  if (!ROLES.includes(role)) throw new HarnessError('BAD_ROLE', `Unknown role ${role}; one of ${ROLES.join(', ')}`);
  const text = readText(join(HARNESS_ROOT, 'prompts', `${role}.md`));
  const match = /^---\n([\s\S]*?)\n---\n/.exec(text);
  if (!match) throw new HarnessError('BAD_PROMPT', `prompts/${role}.md has no frontmatter`);
  const meta = {};
  for (const line of match[1].split('\n')) {
    const at = line.indexOf(':');
    if (at === -1) continue;
    const key = line.slice(0, at).trim(), raw = line.slice(at + 1).trim();
    meta[key] = /^[[{]/.test(raw) ? JSON.parse(raw) : /^\d+$/.test(raw) ? Number(raw) : raw;
  }
  return { meta, body: text.slice(match[0].length).trim(), schemaPath: join(HARNESS_ROOT, 'schemas', `${meta.schema}.schema.json`) };
}

export const fill = (text, vars) => text.replace(/\{\{([A-Z_]+)\}\}/g, (all, key) => (vars[key] ?? all));

// --- Context pack --------------------------------------------------------------

function fileItem(campaign, use, mode) {
  const locator = use.slice(use.indexOf(':') + 1);
  const parsed = parseLocator(locator);
  const path = join(campaign.dir, parsed.path);
  if (!existsSync(path)) return null;
  const text = readText(path);
  if (mode === 'digest') {
    const digest = digestOf(text);
    return { title: `${parsed.path} (digest)`, locator, text: digest || outlineText(text, parsed.path) };
  }
  if (mode === 'outline') return { title: `${parsed.path} (outline)`, locator, text: outlineText(text, parsed.path) };
  try { return { title: locator, locator, text: extract(text, parsed), full: text }; }
  catch { return null; }
}

function special(campaign, use, task) {
  const { paths } = campaign;
  switch (use) {
    case 'direction': {
      const direction = task.direction;
      return direction && { title: `Your direction ${direction.id}`, locator: 'directions.json', text: '```json\n' + JSON.stringify(direction, null, 2) + '\n```' };
    }
    case 'ledger:direction': {
      const entries = readJsonl(paths.ledger).filter(entry => entry.direction === task.direction?.id);
      if (!entries.length) return null;
      return { title: 'Earlier attempts on this direction', locator: 'ledger.jsonl', text: entries.slice(-10).map(entry => `- ${entry.id} [${entry.outcome}/${entry.status}] ${truncate(entry.claim, 120)} — ${truncate(entry.result || entry.failurePoint || '', 200)} → next: ${truncate(entry.next, 120)} (${entry.file || ''})`).join('\n') };
    }
    case 'promising-attempts': {
      const entries = liveLedger(campaign).filter(entry => entry.status === 'promising' || entry.status === 'resolved');
      if (!entries.length) return null;
      const blocks = entries.slice(-4).map(entry => {
        const file = entry.file && join(campaign.dir, entry.file);
        const body = file && existsSync(file) ? readText(file) : '';
        return `### ${entry.id} (${entry.outcome}) — ${entry.file}\n${body}`;
      });
      return { title: 'Promising attempts (full text)', locator: 'ledger.jsonl', text: blocks.join('\n\n'), shrink: () => entries.slice(-4).map(entry => `- ${entry.id} [${entry.outcome}] ${truncate(entry.result, 300)} → ${entry.file}`).join('\n') };
    }
    case 'certificates:index': {
      const status = certificateStatus(campaign);
      const lines = Object.entries(status).map(([name, value]) => {
        const manifestPath = join(paths.certificates, name, 'manifest.json');
        const manifest = existsSync(manifestPath) ? readJson(manifestPath) : {};
        return `- certificates/${name}: receipt ${value.reason}; claim: ${truncate(manifest.claim || '?', 200)}; not claimed: ${truncate(manifest.notClaimed || '?', 160)}; verify: \`${(manifest.verifier?.command || []).join(' ')}\``;
      });
      return lines.length ? { title: 'Certificates', locator: 'certificates/', text: lines.join('\n') } : null;
    }
    case 'certificate': {
      const dir = join(paths.certificates, task.certificate || '');
      if (!task.certificate || !existsSync(dir)) return null;
      const manifest = existsSync(join(dir, 'manifest.json')) ? readText(join(dir, 'manifest.json')) : '{}';
      const files = listFiles(dir).map(path => `- ${rel(dir, path)} (${estimateTokens(readText(path))} tok)`).join('\n');
      return { title: `Certificate ${task.certificate}`, locator: `certificates/${task.certificate}/`, text: `manifest.json:\n\`\`\`json\n${manifest}\n\`\`\`\nFiles:\n${files}` };
    }
    case 'reviews:defects': {
      const summary = reviewSummary(campaign);
      const defects = summary.current.flatMap(review => (review.defects || []).filter(defect => defect.severity !== 'minor').map(defect => `- ${review.file}#${defect.id} [${defect.severity}] at ${defect.location}: ${defect.description}${defect.suggestedRepair ? `\n  Suggested repair: ${defect.suggestedRepair}` : ''}`));
      const minors = summary.current.flatMap(review => (review.defects || []).filter(defect => defect.severity === 'minor').map(defect => `- ${review.file}#${defect.id} [minor] ${truncate(defect.description, 160)}`));
      return { title: 'Defects to repair (current PROOF.md)', locator: 'reviews/', text: [...defects, ...minors].join('\n') || 'No defects recorded.' };
    }
    case 'proof-appendices': {
      const dir = join(campaign.dir, 'proof');
      if (!existsSync(dir)) return null;
      return { title: 'Proof appendices (expand as needed)', locator: 'proof/', text: listFiles(dir).map(path => `- ${rel(campaign.dir, path)} (~${estimateTokens(readText(path))} tok)`).join('\n') };
    }
    case 'lean:project': {
      const dir = join(campaign.dir, 'formal');
      if (!existsSync(dir)) return { title: 'Lean project', locator: 'formal/', text: 'No formal/ project yet. Create formal/lakefile.toml and formal/lean-toolchain (copy harness/templates/lean/).' };
      const files = listFiles(dir).filter(path => /\.(lean|toml|json)$|lean-toolchain$/.test(path)).map(path => `- ${rel(campaign.dir, path)}`);
      const toolchain = existsSync(join(dir, 'lean-toolchain')) ? readText(join(dir, 'lean-toolchain')).trim() : 'unpinned';
      return { title: 'Lean project', locator: 'formal/', text: `Toolchain: ${toolchain}\n${files.join('\n')}` };
    }
    case 'lean:statement': {
      const lean = campaign.state.lean;
      const dir = join(campaign.dir, 'formal');
      const statements = existsSync(dir) ? listFiles(dir).filter(path => /Statement\.lean$/.test(path)) : [];
      if (!statements.length) return null;
      const receiptPath = join(paths.receipts, 'lean-latest.json');
      const closure = existsSync(receiptPath) ? readJson(receiptPath).statement?.closure : null;
      return { title: 'Lean statement (pinned once audited)', locator: rel(campaign.dir, statements[0]), text: statements.map(path => `\`${rel(campaign.dir, path)}\`\n\`\`\`lean\n${readText(path)}\n\`\`\``).join('\n') + (lean ? `\nClaim ${lean.claim}; main ${lean.main}; polarity ${campaign.state.polarity}${lean.claimFingerprint ? `; pinned fingerprint ${lean.claimFingerprint.slice(0, 12)}` : ''}` : '') + (closure ? `\nProject definitions the statement depends on: ${closure.join(', ')}` : '') };
    }
    case 'lean:receipt': {
      const path = join(paths.receipts, 'lean-latest.json');
      if (!existsSync(path)) return null;
      const receipt = readJson(path);
      return { title: `Last Lean gate: ${receipt.status}`, locator: 'receipts/lean-latest.json', text: [...(receipt.failures || []).map(failure => `- ${failure}`), receipt.build?.exitCode ? '\nBuild tail:\n```\n' + receipt.build.tail + '\n```' : ''].join('\n') };
    }
    default: return undefined;
  }
}

export function buildPack(campaign, inputs, task, budget) {
  const items = [];
  for (const [order, input] of inputs.entries()) {
    const use = input.use;
    let item = special(campaign, use, task);
    if (item === undefined) {
      const mode = use.startsWith('digest:') ? 'digest' : use.startsWith('outline:') ? 'outline' : 'full';
      item = fileItem(campaign, use, mode);
    }
    if (item) items.push({ ...item, priority: input.priority ?? 5, order, tokens: estimateTokens(item.text) });
  }
  // Highest priority first; within a priority, the prompt's own order.
  items.sort((a, b) => a.priority - b.priority || a.order - b.order);
  const included = [], omitted = [];
  let used = 0;
  for (const item of items) {
    if (used + item.tokens <= budget) { included.push(item); used += item.tokens; continue; }
    // Degrade gracefully: producer digest, then outline, then just the locator.
    const fallback = item.shrink?.() || (item.full && (digestOf(item.full) || outlineText(item.full, item.locator.split('#')[0])));
    const fallbackTokens = estimateTokens(fallback || '');
    if (fallback && used + fallbackTokens <= budget) {
      included.push({ ...item, title: `${item.title} (shortened: expand with discover show)`, text: fallback, tokens: fallbackTokens });
      used += fallbackTokens;
    } else omitted.push(item);
  }
  return { included, omitted, used, budget };
}

// --- Prompt assembly -------------------------------------------------------------

// Order is deliberate. Static text first (preamble, role) so repeated runs share
// a cacheable prefix; then slow-changing campaign constants; then the capsule;
// then the volatile pack and step. Providers cache by prefix.
export function assemblePrompt({ campaign, role, task = {}, runId, capsule }) {
  const { meta, body } = loadRole(role);
  const vars = {
    RUN: runId, SLUG: campaign?.slug || '', DIRECTION: task.direction?.id || '', ATTEMPT: task.attempt || '',
    CERTIFICATE: task.certificate || '', PROOF_SHA: task.proofSha256 || '', MAIN: task.main || '', GOAL: task.goal || ''
  };
  const preamble = readText(join(HARNESS_ROOT, 'prompts', '_preamble.md')).trim();
  const parts = [preamble, `# Role: ${role}\n\n${fill(body, vars)}`];
  let pack = { included: [], omitted: [], used: 0, budget: meta.budget };
  if (campaign) {
    const statement = existsSync(campaign.paths.statement) ? readText(campaign.paths.statement) : '';
    let exact = '';
    try { exact = extract(statement, { anchor: 'exact-statement' }); } catch { exact = truncate(statement, 4000); }
    parts.push(`# Campaign constants\n\nCampaign \`${campaign.slug}\`: ${campaign.state.title}. Polarity: ${campaign.state.polarity}.${campaign.state.fixture ? ' This is a harness fixture, not research.' : ''}\n\n${exact}`);
    parts.push(capsule.trim());
    pack = buildPack(campaign, meta.inputs || [], task, meta.budget || 8000);
    const packText = pack.included.map(item => `## ${item.title}\n<!-- locator: ${item.locator} · ~${item.tokens} tok -->\n${item.text}`).join('\n\n');
    const omittedText = pack.omitted.length ? `\n\n## Available on demand (not included)\n${pack.omitted.map(item => `- ${item.locator} (~${item.tokens} tok): discover show`).join('\n')}` : '';
    parts.push(`# Context pack (~${pack.used} of ${pack.budget} tokens)\n\n${packText || '(empty)'}${omittedText}`);
  }
  const writes = (meta.writes || []).map(pattern => fill(pattern, vars));
  parts.push([
    '# This step',
    '',
    `RUN_ID: ${runId}`,
    `ROLE: ${role}`,
    `CAMPAIGN_DIR: ${campaign ? rel(campaign.config.root, campaign.dir) : task.root || ''}`,
    `WRITES: ${writes.join(', ') || '(none)'}`,
    task.note ? `TASK: ${task.note}` : null,
    '',
    `Work from the repository root. Write only these paths (relative to CAMPAIGN_DIR): ${writes.join(', ') || 'none'}. Then return the JSON object for schema \`${meta.schema}\`.`
  ].filter(line => line !== null).join('\n'));
  const text = parts.join('\n\n---\n\n') + '\n';
  return { text, meta, writes, pack, tokens: estimateTokens(text) };
}

export function listRolePrompts() {
  return readdirSync(join(HARNESS_ROOT, 'prompts')).filter(name => name.endsWith('.md') && !name.startsWith('_')).map(name => name.replace(/\.md$/, ''));
}
