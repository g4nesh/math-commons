import assert from 'node:assert/strict';
import { mkdirSync, readFileSync, readdirSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { test } from 'node:test';
import { certificateStatus, certify } from '../lib/certificate.mjs';
import { claimProblems, parseProof } from '../lib/claims.mjs';
import { syncRuntime } from '../lib/install.mjs';
import { auditSource, claimFingerprint, parseAuditOutput, scanSource, stripComments } from '../lib/lean.mjs';
import { globToRegExp, scopeViolations } from '../lib/orchestrator.mjs';
import { ROLES, assemblePrompt, buildPack, loadRole } from '../lib/prompt.mjs';
import { pickRuntime } from '../lib/runtimes.mjs';
import { anchorOf, digestOf, extract, outline, parseLocator, tableRows } from '../lib/sections.mjs';
import { HARNESS_ROOT, REPO_ROOT, estimateTokens } from '../lib/util.mjs';
import { eulerCampaign, tempConfig } from './helpers.mjs';

test('sections: outline, extraction, locators and digests', () => {
  const text = '# Title\n<!-- digest -->\nshort\n<!-- /digest -->\n## Exact statement\nA\n### Detail\nB\n## Source\nC\n```\n# not a heading\n```\n';
  const headings = outline(text);
  assert.deepEqual(headings.map(entry => entry.anchor), ['title', 'exact-statement', 'detail', 'source']);
  assert.equal(extract(text, { anchor: 'exact-statement' }), '## Exact statement\nA\n### Detail\nB');
  assert.equal(extract(text, { from: 2, to: 4 }), '<!-- digest -->\nshort\n<!-- /digest -->');
  assert.deepEqual(parseLocator('PROOF.md#c3'), { path: 'PROOF.md', anchor: 'c3' });
  assert.deepEqual(parseLocator('a/b.md:L3-L9'), { path: 'a/b.md', from: 3, to: 9 });
  assert.equal(digestOf(text), 'short');
  assert.equal(anchorOf('C1. The value at `forty`'), 'c1-the-value-at-forty');
  assert.deepEqual(tableRows('| ID | Status |\n|---|---|\n| C1 | proved |\n\ntext'), [{ id: 'C1', status: 'proved' }]);
});

test('lean: comments and strings are ignored, cheats are flagged', () => {
  const source = '/- a sorry in a /- nested -/ comment -/\n-- admit here too\ndef s := "sorry"\ntheorem t : True := by\n  native_decide\naxiom bad : False\nset_option debug.skipKernelTC true\n';
  assert.doesNotMatch(stripComments(source), /admit|nested/);
  const rules = scanSource(source, 'X.lean').map(finding => finding.rule);
  assert.ok(!rules.includes('sorry') && !rules.includes('admit'));
  assert.ok(rules.includes('native-decide') && rules.includes('axiom-declaration') && rules.includes('skip-kernel'));
  assert.deepEqual(scanSource('theorem t : 1 = 1 := by sorry\n', 'Y.lean').map(item => item.rule), ['sorry']);
});

test('lean: audit records parse from JSON and plain message formats', () => {
  const record = { role: 'main', name: 'A.main', found: true, axioms: ['propext'] };
  const json = JSON.stringify({ severity: 'information', data: `HARNESS_AUDIT_JSON:${JSON.stringify(record)}` });
  const plain = `Audit.lean:70:0: info: HARNESS_AUDIT_JSON:${JSON.stringify(record)}`;
  assert.deepEqual(parseAuditOutput(json).records, [record]);
  assert.deepEqual(parseAuditOutput(plain).records, [record]);
  assert.equal(parseAuditOutput(JSON.stringify({ severity: 'error', data: 'boom' })).errors[0], 'boom');
  const source = auditSource({ modules: ['A.Main'], claim: 'A.Claim', main: 'A.main', namespace: 'A', endpoints: ['A.lemma'] });
  assert.match(source, /^#harness_audit_claim A\.Claim under A$/m);
  assert.match(source, /^#harness_audit_main A\.main for A\.Claim$/m);
  assert.match(source, /^#harness_audit_endpoint A\.lemma$/m);
});

test('lean: the statement fingerprint ignores order but not content or toolchain', () => {
  const claim = { name: 'A.Claim', type: 'Prop', closure: [{ name: 'A.Claim', value: 'x' }, { name: 'A.P', value: 'y' }] };
  const reordered = { ...claim, closure: [...claim.closure].reverse() };
  assert.equal(claimFingerprint(claim, 'v1'), claimFingerprint(reordered, 'v1'));
  assert.notEqual(claimFingerprint(claim, 'v1'), claimFingerprint({ ...claim, closure: [{ name: 'A.Claim', value: 'x' }, { name: 'A.P', value: 'z' }] }, 'v1'));
  assert.notEqual(claimFingerprint(claim, 'v1'), claimFingerprint(claim, 'v2'));
});

test('claims: the structural checks catch dishonest or broken proofs', () => {
  const config = tempConfig();
  const campaign = eulerCampaign(config, { formal: false });
  const write = text => { writeFileSync(campaign.paths.proof, text); return claimProblems(parseProof(campaign.paths.proof)); };
  const base = (kind, rows, sections = '### C1\nx\n### C2\ny\n') => `# T\n**Result kind:** \`${kind}\`\n**Main claim:** \`C2\`\n## Claims\n| ID | Statement | Status | Depends on | Evidence |\n|---|---|---|---|---|\n${rows}\n## Proofs\n${sections}`;
  assert.deepEqual(write(base('complete', '| C1 | a | proved | — | |\n| C2 | b | proved | C1 | |')), []);
  assert.ok(write(base('complete', '| C1 | a | assumed | — | Smith 2020 |\n| C2 | b | proved | C1 | |')).some(item => /assumed/.test(item)));
  assert.ok(write(base('complete', '| C1 | a | open | — | |\n| C2 | b | proved | C1 | |')).some(item => /open claim/.test(item)));
  assert.ok(write(base('complete', '| C1 | a | computed | — | |\n| C2 | b | proved | C1 | |')).some(item => /must cite certificates/.test(item)));
  assert.ok(write(base('complete', '| C1 | a | proved | C2 | |\n| C2 | b | proved | C1 | |')).some(item => /Circular/.test(item)));
  assert.ok(write(base('solved', '| C1 | a | proved | — | |\n| C2 | b | proved | C9 | |')).some(item => /unknown claim C9/.test(item)));
  assert.ok(write(base('complete', '| C1 | a | proved | — | |\n| C2 | b | proved | C1 | |', '### C1\nx\n')).some(item => /C2: no "### C2"/.test(item)));
});

test('certificates: receipts bind file hashes and lint the verifier', async () => {
  const config = tempConfig();
  const campaign = eulerCampaign(config, { formal: false });
  const dir = join(campaign.paths.certificates, 'demo');
  mkdirSync(dir, { recursive: true });
  writeFileSync(join(dir, 'certificate.json'), '{"x": "3"}');
  writeFileSync(join(dir, 'verify.py'), 'import json\nx = int(json.load(open("certificate.json"))["x"])\nprint("VERIFIED" if x * x == 9 else "FAIL")\nraise SystemExit(0 if x * x == 9 else 1)\n');
  writeFileSync(join(dir, 'manifest.json'), JSON.stringify({ name: 'demo', claim: 'C1', notClaimed: 'nothing else', verifier: { path: 'verify.py', command: ['python3', '-I', '-B', 'verify.py'] } }));
  const receipt = await certify(campaign, 'demo', null);
  assert.equal(receipt.status, 'pass');
  assert.deepEqual(receipt.warnings, ['no recorded verifier run: write the verifier through the verify-certificate step']);
  assert.equal(certificateStatus(campaign).demo.ok, true);
  writeFileSync(join(dir, 'certificate.json'), '{"x": "4"}');
  assert.match(certificateStatus(campaign).demo.reason, /stale/);
  assert.equal((await certify(campaign, 'demo', null)).status, 'fail');
  writeFileSync(join(dir, 'verify.py'), 'import math\nprint(math.sqrt(9))\n');
  assert.ok((await certify(campaign, 'demo', null)).warnings.some(item => /floatJustification/.test(item)));
});

test('prompts: every role loads, has a strict schema, and assembles cache-friendly', () => {
  for (const role of ROLES) {
    const { meta, schemaPath } = loadRole(role);
    assert.equal(meta.role, role);
    assert.ok(Array.isArray(meta.writes) && Array.isArray(meta.inputs) && meta.budget > 0, role);
    const schema = JSON.parse(readFileSync(schemaPath, 'utf8'));
    const strict = node => {
      if (node.type === 'object') {
        assert.equal(node.additionalProperties, false, `${role}: objects must be closed for Codex strict outputs`);
        assert.deepEqual([...node.required].sort(), Object.keys(node.properties).sort(), `${role}: every property required`);
        Object.values(node.properties).forEach(strict);
      }
      if (node.type === 'array') strict(node.items);
    };
    strict(schema);
  }
  const config = tempConfig();
  const campaign = eulerCampaign(config, { formal: false });
  const prompt = assemblePrompt({ campaign, role: 'screen', task: {}, runId: 'R1', capsule: '# Capsule\nX' });
  const order = ['You are one step', '# Role: screen', '# Campaign constants', '# Capsule', '# Context pack', '# This step'].map(marker => prompt.text.indexOf(marker));
  assert.ok(order.every((position, index) => position >= 0 && (index === 0 || position > order[index - 1])), `prompt order ${order}`);
  const other = assemblePrompt({ campaign, role: 'screen', task: {}, runId: 'R2', capsule: '# Capsule\nY' });
  const prefix = prompt.text.slice(0, prompt.text.indexOf('# Capsule'));
  assert.equal(other.text.slice(0, prefix.length), prefix, 'runs share the static prefix');
});

test('context packs degrade to digests and locators under budget', () => {
  const config = tempConfig();
  const campaign = eulerCampaign(config, { formal: false });
  writeFileSync(join(campaign.dir, 'BIG.md'), `# Big\n<!-- digest -->\nthe gist\n<!-- /digest -->\n## Part\n${'word '.repeat(4000)}\n`);
  const pack = buildPack(campaign, [{ use: 'file:STATEMENT.md', priority: 1 }, { use: 'file:BIG.md', priority: 2 }], {}, 600);
  assert.equal(pack.included[0].locator, 'STATEMENT.md');
  assert.match(pack.included[1].title, /shortened/);
  assert.equal(pack.included[1].text, 'the gist');
  assert.ok(pack.used <= 600);
  const tight = buildPack(campaign, [{ use: 'file:STATEMENT.md', priority: 1 }], {}, 10);
  assert.equal(tight.omitted.length + tight.included.length, 1);
  assert.ok(estimateTokens('abcd') === 1);
});

test('routing sends review and verification away from the author', () => {
  const config = tempConfig();
  assert.equal(pickRuntime(config, 'draft'), 'codex');
  assert.equal(pickRuntime(config, 'review', { index: 0 }), 'claude');
  assert.equal(pickRuntime(config, 'review', { index: 0, exclude: 'claude' }), 'codex');
  assert.equal(pickRuntime(config, 'verify-certificate', { exclude: 'claude' }), 'codex');
  assert.equal(pickRuntime(config, 'explore', { index: 1 }), 'claude');
  assert.equal(pickRuntime(config, 'draft', { override: 'manual' }), 'manual');
});

test('scope globs', () => {
  assert.ok(globToRegExp('attempts/D1/**').test('attempts/D1/001.md'));
  assert.ok(globToRegExp('formal/**').test('formal/A/B.lean'));
  assert.ok(!globToRegExp('PROOF.md').test('PROOF.mdx'));
  assert.ok(!globToRegExp('attempts/D1/**').test('attempts/D2/001.md'));
  const result = scopeViolations({ 'a.md': '1' }, { 'a.md': '2', 'campaign.json': 'x', 'b.md': '3' }, ['a.md']);
  assert.deepEqual(result.violations, ['b.md']);
});

test('runtime files in .claude, .agents and .codex match harness/runtime', () => {
  const drift = syncRuntime(REPO_ROOT, { check: true });
  assert.deepEqual(drift, [], 'run: node harness/bin/discover.mjs install-runtime');
  assert.ok(readdirSync(join(HARNESS_ROOT, 'runtime', 'skills')).length >= 3);
});
