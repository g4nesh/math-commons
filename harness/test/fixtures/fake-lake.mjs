#!/usr/bin/env node
// Stands in for `lake` where Lean is not installed. It never claims to check a
// proof: it emits audit records derived from the source text so the gate's
// orchestration, parsing and fingerprint logic can be tested. CI runs the
// same tests against the real toolchain.
import { createHash } from 'node:crypto';
import { existsSync, readFileSync, readdirSync } from 'node:fs';
import { join } from 'node:path';

const argv = process.argv.slice(2);
const sha = text => createHash('sha256').update(text).digest('hex');
const sources = dir => readdirSync(dir, { withFileTypes: true }).flatMap(entry => entry.name.startsWith('.') ? [] : entry.isDirectory() ? sources(join(dir, entry.name)) : entry.name.endsWith('.lean') ? [join(dir, entry.name)] : []);
if (argv[0] === '--version') { console.log('Lake (fake, harness tests)'); process.exit(0); }
if (argv[0] === 'build') { console.log('Build completed successfully.'); process.exit(0); }
if (argv[0] === 'env' && argv[1] === 'which') process.exit(1);
if (argv[0] === 'env' && argv[1] === 'lean') {
  const file = argv.at(-1);
  const audit = readFileSync(file, 'utf8');
  const all = sources(process.cwd()).map(path => readFileSync(path, 'utf8')).join('\n');
  const statementText = sources(process.cwd()).filter(path => path.endsWith('Statement.lean')).map(path => readFileSync(path, 'utf8')).join('\n');
  const emit = record => console.log(JSON.stringify({ severity: 'information', data: 'HARNESS_AUDIT_JSON:' + JSON.stringify(record) }));
  const claim = /^#harness_audit_claim (\S+) under (\S+)/m.exec(audit);
  const main = /^#harness_audit_main (\S+) for (\S+)/m.exec(audit);
  if (claim) emit({ role: 'claim', name: claim[1], found: /def Claim/.test(all), kind: 'def', type: 'Prop', isProp: true, axioms: [],
    closure: [{ name: claim[1], kind: 'def', type: 'Prop', value: sha(statementText) }] });
  if (main) {
    const found = /theorem main/.test(all);
    emit({ role: 'main', name: main[1], found, kind: 'theorem', type: `Not ${main[2]}`, provesClaim: false, refutesClaim: found, axioms: /sorry/.test(all) ? ['sorryAx'] : ['propext'] });
  }
  for (const endpoint of audit.matchAll(/^#harness_audit_endpoint (\S+)/gm)) emit({ role: 'endpoint', name: endpoint[1], found: all.includes(endpoint[1].split('.').at(-1)), kind: 'theorem', type: '?', axioms: [] });
  if (/run_cmd/.test(audit)) emit({ role: 'owned', name: 'EulerPoly.Claim', kind: 'def', axioms: [] });
  process.exit(0);
}
console.error(`fake lake: unsupported ${argv.join(' ')}`);
process.exit(1);
