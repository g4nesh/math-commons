#!/usr/bin/env node
import { spawnSync } from 'node:child_process';
import { existsSync, statSync } from 'node:fs';
import { join, resolve } from 'node:path';
import { STAGES, TERMINAL, campaignsRoot, createCampaign, listCampaigns, loadCampaign, loadConfig, saveCampaign, transition, usageTotals } from '../lib/campaign.mjs';
import { writeCapsule } from '../lib/capsule.mjs';
import { certify } from '../lib/certificate.mjs';
import { evaluate, verificationLabel } from '../lib/gates.mjs';
import { syncRuntime } from '../lib/install.mjs';
import { leanGate } from '../lib/lean.mjs';
import { campaignLeanGate, newRunId, nextAction, recordManual, runLoop, step } from '../lib/orchestrator.mjs';
import { assemblePrompt, loadRole } from '../lib/prompt.mjs';
import { openPullRequest, publish } from '../lib/publish.mjs';
import { reconcile, reconcileMarkdown, writeReconcile } from '../lib/reconcile.mjs';
import { pickRuntime, runRuntime } from '../lib/runtimes.mjs';
import { outlineText, resolveLocator } from '../lib/sections.mjs';
import { HarnessError, listFiles, parseArgs, readJson, readText, rel, slugify } from '../lib/util.mjs';

const HELP = `discover: the Math Commons discovery harness

Campaign lifecycle (stages: ${STAGES.join(' → ')}; terminal: rejected, parked)
  new SLUG --title T [--source URL] [--polarity prove|disprove|unknown] [--fixture]
  status [SLUG] [--json]              Stage, gate, label, budget for one or all campaigns
  next SLUG                           The action the orchestrator would take (no side effects)
  step SLUG [--runtime R] [--as ID]   Take exactly one action
  run SLUG [--until STAGE] [--max-steps N] [--runtime R]
                                      Loop steps until a stage, a stop, a pending manual run or the budget
  record SLUG RUN-ID [--answer FILE]  Ingest a manual run's JSON answer (default .runs/RUN-ID/answer.json)
  set SLUG KEY VALUE [--reason TEXT]  polarity | budget | formal on|off | stage S | pending none | lean JSON
  scout --area TEXT [--count N] [--runtime R]
                                      Propose significant open targets (writes ${'campaigns'}/_scout/)

Context economy
  capsule SLUG                        Rebuild and print CAPSULE.md (deterministic, no model)
  pack SLUG --role R [--direction D] [--certificate C] [--print]
                                      Assemble a role prompt and show its token manifest
  show LOCATOR                        Print one section: path#heading or path:L10-L40
  outline FILE                        Headings with line ranges and token estimates
  changed SLUG --since ISO|GIT-REF    Files changed in the campaign since a time or commit

Verification
  gate SLUG [--stage S] [--json]      Deterministic gate checks for a stage
  certify SLUG NAME [-- CMD ...]      Run a certificate's verifier and write its receipt
  lean-gate SLUG [--json]             Lean gate for the campaign (writes receipts/lean-latest.json)
  lean-gate --project DIR --module M[,M] --claim C --main T --polarity prove|disprove
            [--namespace N] [--endpoints A,B] [--expect FINGERPRINT] [--replay auto] [--json]
  pin-statement SLUG                  Pin the audited Lean statement fingerprint (after a faithful audit)

Publication and upkeep
  publish SLUG [--open-pr] [--base main]
  reconcile [--write] [--json] [--strict]
                                      Inventory receipt across campaigns (--strict: exit 1 on issues)
  install-runtime [--check]           Sync skills/agents into .claude/, .agents/ and .codex/
  doctor                              Check runtimes, Lean, Python and configuration

Runtimes: claude (Claude Code print mode), codex (codex exec), manual (prompt to file,
for interactive sessions). Configure in harness/config.json or discovery.config.json.`;

const out = text => process.stdout.write(text.endsWith('\n') ? text : text + '\n');
const json = value => out(JSON.stringify(value, null, 2));

function describe(outcome) {
  const { action } = outcome;
  const head = action.kind === 'transition' ? `→ ${action.to}: ${action.reason}`
    : action.kind === 'stop' ? `stop: ${action.reason}`
    : action.kind === 'role' || action.kind === 'explore' ? `${action.role || 'explore'}${action.tasks ? ` ×${action.tasks.length}` : ''}`
    : action.kind;
  const lines = [head];
  for (const result of outcome.results || []) lines.push(`  ${result.role || ''} ${result.runId} [${result.runtime || ''}] ${result.status}${result.detail ? `: ${result.detail}` : ''}${result.summary ? ` — ${result.summary}` : ''}`);
  if (outcome.receipt) lines.push(`  lean gate ${outcome.receipt.status}${outcome.receipt.failures?.length ? `: ${outcome.receipt.failures.join('; ')}` : ''}`);
  if (outcome.receipts) lines.push(...outcome.receipts.map(item => `  certificate ${item.name}: ${item.status}`));
  if (outcome.pinned !== undefined) lines.push(`  pinned: ${outcome.pinned}${outcome.fingerprint ? ` ${outcome.fingerprint.slice(0, 16)}` : ''}${outcome.reason ? ` (${outcome.reason})` : ''}`);
  if (outcome.result) lines.push(`  ${JSON.stringify(outcome.result)}`);
  return lines.join('\n');
}

function statusLine(campaign) {
  const gate = evaluate(campaign);
  const usage = usageTotals(campaign);
  return `${campaign.slug.padEnd(28)} ${campaign.state.stage.padEnd(10)} gate ${gate.ok ? 'pass' : 'open'}  ${verificationLabel(campaign).emoji}  ${usage.tokens.toLocaleString('en-US')} tok${campaign.state.fixture ? '  (fixture)' : ''}${campaign.state.pending ? `  pending ${campaign.state.pending.runId}` : ''}`;
}

async function main(argv) {
  const { positional, flags } = parseArgs(argv);
  const [command, ...rest] = positional;
  const config = loadConfig({ configPath: flags.config });
  const campaign = () => {
    if (!rest[0]) throw new HarnessError('USAGE', `${command} needs a campaign slug`);
    return loadCampaign(config, rest[0]);
  };
  switch (command) {
    case undefined: case 'help': return out(HELP);
    case 'new': {
      const slug = rest[0] || slugify(flags.title || '');
      const created = createCampaign(config, slug, { title: flags.title, source: flags.source, polarity: flags.polarity || 'unknown', fixture: Boolean(flags.fixture) });
      writeCapsule(created);
      return out(`Created ${rel(config.root, created.dir)}/. Fill STATEMENT.md (exact statement, source, scope), then: discover step ${slug}`);
    }
    case 'status': {
      if (rest[0]) {
        const c = campaign();
        const gate = evaluate(c);
        if (flags.json) return json({ state: c.state, gate, label: verificationLabel(c), usage: usageTotals(c) });
        return out([statusLine(c), ...gate.checks.map(item => `  ${item.ok ? '✓' : '✗'} ${item.label}${!item.ok && item.fix ? ` → ${item.fix}` : ''}`), `  next: ${JSON.stringify(nextAction(c))}`].join('\n'));
      }
      const all = listCampaigns(config).map(slug => loadCampaign(config, slug));
      if (flags.json) return json(all.map(c => ({ slug: c.slug, stage: c.state.stage, gate: evaluate(c).ok, label: verificationLabel(c).emoji })));
      return out(all.length ? all.map(statusLine).join('\n') : `No campaigns under ${config.campaignRoot}/. Start one: discover new SLUG --title T`);
    }
    case 'capsule': return out(writeCapsule(campaign()));
    case 'next': return json(nextAction(campaign()));
    case 'gate': {
      const c = campaign();
      const result = evaluate(c, flags.stage || c.state.stage);
      if (flags.json) return json(result);
      out([`${result.stage}: ${result.ok ? 'PASS' : 'OPEN'}`, ...result.checks.map(item => `  ${item.ok ? '✓' : '✗'} ${item.label}${!item.ok && item.fix ? ` → ${item.fix}` : ''}`)].join('\n'));
      process.exitCode = result.ok ? 0 : 1;
      return;
    }
    case 'pack': {
      const c = campaign();
      if (!flags.role) throw new HarnessError('USAGE', 'pack needs --role');
      const directions = existsSync(c.paths.directions) ? readJson(c.paths.directions).directions : [];
      const task = { direction: directions.find(item => item.id === flags.direction) || directions[0], certificate: flags.certificate, attempt: 'NNN' };
      const prompt = assemblePrompt({ campaign: c, role: flags.role, task, runId: 'preview', capsule: writeCapsule(c) });
      if (flags.print) return out(prompt.text);
      return out([`${flags.role}: ~${prompt.tokens} tokens total; pack ${prompt.pack.used}/${prompt.pack.budget}`, ...prompt.pack.included.map(item => `  + ${item.title} (~${item.tokens})`), ...prompt.pack.omitted.map(item => `  - omitted ${item.locator} (~${item.tokens})`), `  writes: ${prompt.writes.join(', ')}`].join('\n'));
    }
    case 'step': {
      const outcome = await step(campaign(), { runtime: flags.runtime || null, as: flags.as || null });
      return out(describe(outcome));
    }
    case 'run': {
      const c = campaign();
      if (flags.until && ![...STAGES, ...TERMINAL].includes(flags.until)) throw new HarnessError('USAGE', `--until must be a stage: ${STAGES.join(', ')}`);
      const log = await runLoop(c, { until: flags.until || null, maxSteps: Number(flags['max-steps'] || 50), runtime: flags.runtime || null, as: flags.as || null, onStep: outcome => out(describe(outcome)) });
      const last = log.at(-1);
      if (last?.stop) out(`stop: ${last.stop}`);
      return out(statusLine(c));
    }
    case 'record': {
      const c = campaign();
      const result = recordManual(c, rest[1], flags.answer ? resolve(flags.answer) : null);
      return out(`${result.role} ${result.runId}: ${result.status}${result.detail ? ` (${result.detail})` : ''}\nnext: ${JSON.stringify(nextAction(c))}`);
    }
    case 'set': {
      const c = campaign();
      const [, key, value] = rest;
      const reason = flags.reason || 'set by operator';
      if (key === 'polarity') { if (!['prove', 'disprove', 'unknown'].includes(value)) throw new HarnessError('USAGE', 'polarity is prove, disprove or unknown'); c.state.polarity = value; }
      else if (key === 'budget') c.state.budget.tokens = Number(value);
      else if (key === 'formal') c.state.formalRequested = value === 'on';
      else if (key === 'pending' && value === 'none') delete c.state.pending;
      else if (key === 'lean') c.state.lean = { project: 'formal', claimFingerprint: null, fingerprintForAudit: null, audit: null, endpoints: [], ...JSON.parse(value) };
      else if (key === 'stage') {
        // An operator reset is a fresh start for the retry counters.
        c.state.counters = { exploreRounds: 0, repairRounds: 0, formalRounds: 0, remaps: c.state.counters.remaps || 0 };
        transition(c, value, `operator: ${reason}`);
        writeCapsule(c);
        return out(statusLine(c));
      }
      else throw new HarnessError('USAGE', `Unknown key ${key}`);
      c.state.history.push({ at: new Date().toISOString(), from: c.state.stage, to: c.state.stage, reason: `operator set ${key}=${value}: ${reason}` });
      saveCampaign(c);
      writeCapsule(c);
      return out(statusLine(c));
    }
    case 'show': {
      if (!rest[0]) throw new HarnessError('USAGE', 'show needs a locator');
      return out(resolveLocator(config.root, rest[0]).text);
    }
    case 'outline': {
      if (!rest[0]) throw new HarnessError('USAGE', 'outline needs a file');
      return out(outlineText(readText(resolve(rest[0])), rest[0]));
    }
    case 'changed': {
      const c = campaign();
      const since = flags.since;
      if (!since) throw new HarnessError('USAGE', 'changed needs --since ISO-TIME or a git ref');
      if (/^\d{4}-\d{2}-\d{2}/.test(since)) {
        const cutoff = Date.parse(since);
        return out(listFiles(c.dir).filter(path => statSync(path).mtimeMs > cutoff).map(path => rel(config.root, path)).join('\n') || '(nothing)');
      }
      const result = spawnSync('git', ['diff', '--name-only', since, '--', rel(config.root, c.dir)], { cwd: config.root, encoding: 'utf8' });
      return out(result.stdout.trim() || '(nothing)');
    }
    case 'certify': {
      const c = campaign();
      const receipt = await certify(c, rest[1], flags['--'] || null);
      writeCapsule(c);
      out(`${rest[1]}: ${receipt.status} (exit ${receipt.exitCode})${receipt.warnings.length ? `\n  warnings: ${receipt.warnings.join('; ')}` : ''}\n  ${receipt.output.tail.slice(-5).join('\n  ')}`);
      process.exitCode = receipt.status === 'pass' ? 0 : 1;
      return;
    }
    case 'lean-gate': {
      let receipt;
      if (flags.project) {
        receipt = await leanGate({
          project: flags.project, modules: String(flags.module || '').split(',').filter(Boolean), claim: flags.claim, main: flags.main,
          polarity: flags.polarity, namespace: flags.namespace, endpoints: flags.endpoints ? String(flags.endpoints).split(',') : [],
          expectFingerprint: flags.expect || null, replay: flags.replay === 'auto' ? 'auto' : null,
          receiptsDir: flags.receipts ? resolve(flags.receipts) : null
        });
      } else {
        const c = campaign();
        receipt = await campaignLeanGate(c, { requirePinned: !flags['allow-unpinned'] });
        writeCapsule(c);
      }
      if (flags.json) json(receipt);
      else out([`lean gate: ${receipt.status.toUpperCase()}`, ...receipt.failures.map(item => `  ✗ ${item}`), receipt.statement?.fingerprint ? `  statement fingerprint ${receipt.statement.fingerprint}` : null, receipt.mainTheorem ? `  ${receipt.main}: ${receipt.mainTheorem.type} · axioms [${(receipt.mainTheorem.axioms || []).join(', ')}]` : null, ...receipt.findings.filter(item => item.severity === 'warning').map(item => `  ! ${item.rule} at ${item.file}:${item.line}`)].filter(Boolean).join('\n'));
      process.exitCode = receipt.status === 'pass' ? 0 : 1;
      return;
    }
    case 'pin-statement': {
      const c = campaign();
      if (c.state.lean?.audit?.verdict !== 'faithful') throw new HarnessError('NOT_AUDITED', 'Pin only after a faithful statement audit (STATEMENT_REVIEW.md)');
      const receipt = await campaignLeanGate(c, { requirePinned: false });
      if (receipt.statement?.fingerprint !== c.state.lean.audit.fingerprint) throw new HarnessError('CHANGED', 'The statement changed since the audit; re-audit it');
      c.state.lean.claimFingerprint = receipt.statement.fingerprint;
      saveCampaign(c);
      return out(`pinned ${c.state.lean.claimFingerprint}`);
    }
    case 'publish': {
      const c = campaign();
      const result = publish(c);
      out(`wrote ${result.page} and ${result.index}`);
      if (c.state.stage === 'publish') transition(c, 'done', 'published');
      if (flags['open-pr']) out(`draft PR: ${openPullRequest(c, { base: flags.base || 'main' })}`);
      return;
    }
    case 'reconcile': {
      const report = reconcile(config);
      if (flags.json) return json(report);
      out(reconcileMarkdown(report));
      if (flags.write) out(`written to ${rel(config.root, writeReconcile(config, report))}`);
      if (flags.strict && report.blockers.length) process.exitCode = 1;
      return;
    }
    case 'scout': {
      if (!flags.area) throw new HarnessError('USAGE', 'scout needs --area');
      const runId = newRunId('scout');
      const runDir = join(campaignsRoot(config), '_scout', '.runs', runId);
      const prompt = assemblePrompt({ campaign: null, role: 'scout', runId, task: { root: config.campaignRoot, note: `Area: ${flags.area}. Return at most ${flags.count || 6} candidates. Write ${config.campaignRoot}/_scout/${new Date().toISOString().slice(0, 10)}-${slugify(flags.area)}.md.` } });
      const runtime = pickRuntime(config, 'scout', { override: flags.runtime || null });
      const result = await runRuntime({ config, runtime, prompt: prompt.text, schemaPath: loadRole('scout').schemaPath, role: prompt.meta, runDir, cwd: config.root });
      if (result.status === 'pending') return out(result.detail);
      if (result.status !== 'ok') throw new HarnessError('RUN_FAILED', result.detail);
      return out((result.answer.candidates || []).map(item => `${item.rank}. ${item.title} — ${item.source}\n   ${item.significance}`).join('\n') || result.answer.summary);
    }
    case 'install-runtime': {
      const drift = syncRuntime(config.root, { check: Boolean(flags.check) });
      if (flags.check && drift.length) { out(`out of date:\n${drift.map(path => `  ${rel(config.root, path)}`).join('\n')}\nRun: discover install-runtime`); process.exitCode = 1; return; }
      return out(drift.length ? `updated:\n${drift.map(path => `  ${rel(config.root, path)}`).join('\n')}` : 'runtime files are in sync');
    }
    case 'doctor': {
      const probe = (bin, args) => { const result = spawnSync(bin, args, { encoding: 'utf8' }); return result.status === 0 ? (result.stdout || result.stderr).trim().split('\n')[0] : null; };
      const rows = [
        ['claude', config.runtimes.claude.bin, probe(config.runtimes.claude.bin, ['--version'])],
        ['codex', config.runtimes.codex.bin, probe(config.runtimes.codex.bin, ['--version'])],
        ['lake', process.env.DISCOVER_LAKE_BIN || 'lake', probe(process.env.DISCOVER_LAKE_BIN || 'lake', ['--version'])],
        ['python3', 'python3', probe('python3', ['--version'])],
        ['gh', 'gh', probe('gh', ['--version'])],
        ['node', process.execPath, process.version]
      ];
      out(rows.map(([name, bin, version]) => `${version ? '✓' : '✗'} ${name.padEnd(8)} ${version || `not found (${bin})`}`).join('\n'));
      out(`routing: author=${config.routing.author}, reviewer=${config.routing.reviewer}; campaigns in ${rel(process.cwd(), campaignsRoot(config)) || '.'}/`);
      const drift = syncRuntime(config.root, { check: true });
      out(drift.length ? `✗ runtime files out of date (${drift.length}); run discover install-runtime` : '✓ runtime files in sync');
      return;
    }
    default:
      throw new HarnessError('USAGE', `Unknown command ${command}. Run: discover help`);
  }
}

main(process.argv.slice(2)).catch(error => {
  process.stderr.write(`${error.code ? `${error.code}: ` : ''}${error.message}\n`);
  process.exitCode = 2;
});
