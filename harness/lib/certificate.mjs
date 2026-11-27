import { spawn } from 'node:child_process';
import { existsSync, readdirSync } from 'node:fs';
import { join } from 'node:path';
import { HarnessError, fileSha256, listFiles, now, readJson, readText, rel, sha256, writeJson, writeText } from './util.mjs';

// A computed claim is accepted only through this path:
// search output -> exact certificate file -> independently written verifier ->
// receipt binding the hashes of every file that was checked.
const RECEIPT = 'receipt.json', OUTPUT = 'output.txt', MANIFEST = 'manifest.json', PROVENANCE = 'provenance.json';
const FLOAT_PATTERNS = [/\bfloat\s*\(/, /^\s*import\s+math\b/m, /^\s*from\s+math\s+import\b/m, /\bnumpy\b/, /\bscipy\b/, /\d\.\d+e[-+]?\d+/];

export function certificateDir(campaign, name) {
  if (!/^[a-z0-9][a-z0-9-]{0,63}$/.test(name)) throw new HarnessError('BAD_NAME', `Certificate names are lowercase slugs: ${name}`);
  return join(campaign.paths.certificates, name);
}

export function hashTree(dir) {
  return Object.fromEntries(listFiles(dir)
    .filter(path => ![RECEIPT, OUTPUT, PROVENANCE].includes(rel(dir, path)))
    .map(path => [rel(dir, path), fileSha256(path)]));
}

// Who generated and who verified a certificate is recorded by the orchestrator
// in provenance.json, outside the model-editable manifest.
export function provenance(campaign, name, update = null) {
  const path = join(certificateDir(campaign, name), PROVENANCE);
  const current = existsSync(path) ? readJson(path) : {};
  if (!update) return current;
  const next = { ...current, ...update };
  writeJson(path, next);
  return next;
}

export function lintVerifier(dir, manifest, origin = {}) {
  const warnings = [];
  const verifier = manifest?.verifier?.path && join(dir, manifest.verifier.path);
  if (verifier && existsSync(verifier) && !manifest.floatJustification) {
    const source = readText(verifier);
    for (const pattern of FLOAT_PATTERNS) {
      if (pattern.test(source)) { warnings.push(`verifier matches ${pattern}: use integers/Fraction or record floatJustification in manifest.json`); break; }
    }
  }
  if (!origin.verifier) warnings.push('no recorded verifier run: write the verifier through the verify-certificate step');
  else if (origin.generator && origin.generator.runId === origin.verifier.runId) warnings.push('generator and verifier were written in the same run; an independent verifier must come from a separate run');
  else if (origin.generator && origin.generator.identity === origin.verifier.identity) warnings.push('generator and verifier share an identity; prefer the other runtime');
  if (!manifest?.claim) warnings.push('manifest.json does not name the PROOF.md claim this certificate supports');
  if (!manifest?.notClaimed) warnings.push('manifest.json should say what the certificate does NOT establish (notClaimed)');
  return warnings;
}

function runCommand(command, cwd, timeoutMs) {
  return new Promise(resolve => {
    const child = spawn(command[0], command.slice(1), { cwd, shell: false, stdio: ['ignore', 'pipe', 'pipe'] });
    const chunks = [];
    let bytes = 0;
    const capture = chunk => { bytes += chunk.length; if (bytes <= 8_388_608) chunks.push(chunk); };
    child.stdout.on('data', capture);
    child.stderr.on('data', capture);
    const timer = setTimeout(() => child.kill('SIGKILL'), timeoutMs);
    child.once('error', error => { clearTimeout(timer); resolve({ code: null, error: error.message, output: '' }); });
    child.once('close', (code, signal) => {
      clearTimeout(timer);
      resolve({ code, signal, output: Buffer.concat(chunks).toString('utf8'), truncated: bytes > 8_388_608 });
    });
  });
}

export async function certify(campaign, name, command, { timeoutMinutes = 30 } = {}) {
  const dir = certificateDir(campaign, name);
  if (!existsSync(dir)) throw new HarnessError('NO_CERTIFICATE', `No directory certificates/${name}`);
  const manifest = existsSync(join(dir, MANIFEST)) ? readJson(join(dir, MANIFEST)) : null;
  const argv = command?.length ? command : manifest?.verifier?.command;
  if (!argv?.length) throw new HarnessError('NO_COMMAND', 'Give the verifier command after -- or set verifier.command in manifest.json');
  const before = hashTree(dir);
  const startedAt = now();
  const result = await runCommand(argv, dir, timeoutMinutes * 60_000);
  const after = hashTree(dir);
  const mutated = Object.keys({ ...before, ...after }).filter(path => before[path] !== after[path]);
  writeText(join(dir, OUTPUT), result.output);
  const lines = result.output.trimEnd().split('\n');
  const receipt = {
    version: 1, name, claim: manifest?.claim || null,
    command: argv, cwd: rel(campaign.config.root, dir),
    exitCode: result.code, signal: result.signal || null, startError: result.error || null,
    status: result.code === 0 && !mutated.length ? 'pass' : 'fail',
    mutatedDuringRun: mutated,
    output: { file: OUTPUT, sha256: sha256(result.output), tail: lines.slice(-20), truncated: Boolean(result.truncated) },
    files: after,
    provenance: provenance(campaign, name),
    warnings: lintVerifier(dir, manifest, provenance(campaign, name)),
    environment: { node: process.version, platform: process.platform },
    startedAt, finishedAt: now()
  };
  writeJson(join(dir, RECEIPT), receipt);
  return receipt;
}

export function certificateStatus(campaign) {
  const root = campaign.paths.certificates;
  if (!existsSync(root)) return {};
  const status = {};
  for (const entry of readdirSync(root, { withFileTypes: true })) {
    if (!entry.isDirectory()) continue;
    const dir = join(root, entry.name);
    const receiptPath = join(dir, RECEIPT);
    if (!existsSync(receiptPath)) { status[entry.name] = { ok: false, reason: 'missing (run discover certify)' }; continue; }
    const receipt = readJson(receiptPath);
    if (receipt.status !== 'pass') { status[entry.name] = { ok: false, reason: `failing (exit ${receipt.exitCode})`, receipt }; continue; }
    const current = hashTree(dir);
    const changed = Object.keys({ ...current, ...receipt.files }).filter(path => current[path] !== receipt.files[path]);
    status[entry.name] = changed.length
      ? { ok: false, reason: `stale (changed since receipt: ${changed.join(', ')})`, receipt }
      : { ok: true, reason: 'pass', receipt, warnings: receipt.warnings || [] };
  }
  return status;
}
