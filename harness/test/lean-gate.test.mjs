import assert from 'node:assert/strict';
import { cpSync, mkdtempSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { test } from 'node:test';
import { leanGate } from '../lib/lean.mjs';
import { HERE, REAL_LAKE } from './helpers.mjs';

// Negative controls for the Lean gate against the real toolchain. Each bad
// fixture compiles; only the gate's checks stand between it and acceptance.
const skip = REAL_LAKE ? false : 'needs a Lean toolchain (lake on PATH); CI runs this';
const project = join(HERE, 'lean-fixtures');
const gate = (name, extra = {}) => leanGate({
  project, modules: [`Fixtures.${name}`], claim: `Fixtures.${name}.Claim`, main: `Fixtures.${name}.main`,
  polarity: 'disprove', namespace: `Fixtures.${name}`, buildTimeoutMs: 20 * 60_000, ...extra
});

test('lean gate: an honest proof passes with only standard axioms', { skip, timeout: 1_800_000 }, async () => {
  const receipt = await gate('Good', { endpoints: ['Fixtures.Good.value_at_forty'], replay: 'auto' });
  assert.equal(receipt.status, 'pass', `${receipt.failures.join('\n')}\nreplay: ${JSON.stringify(receipt.replay)}`);
  assert.equal(receipt.mainTheorem.refutesClaim, true);
  assert.ok(receipt.mainTheorem.axioms.every(axiom => ['propext', 'Classical.choice', 'Quot.sound'].includes(axiom)));
  assert.deepEqual(receipt.statement.closure.sort(), ['Fixtures.Good.Claim', 'Fixtures.Good.IsPrime']);
  assert.match(receipt.statement.fingerprint, /^[0-9a-f]{64}$/);
  assert.equal(receipt.owned.checked, true, JSON.stringify(receipt.owned));
  // Either the toolchain ships leanchecker and the replay passed, or the receipt says it was not run.
  assert.ok(receipt.replay && (receipt.replay.exitCode === 0 || (receipt.replay.command === null && /not run/.test(receipt.replay.tail))), JSON.stringify(receipt.replay));
  assert.ok(receipt.replay.command || receipt.scope.notClaimed.some(item => /kernel replay/i.test(item)));
});

test('lean gate: sorry is rejected', { skip, timeout: 1_800_000 }, async () => {
  const receipt = await gate('Sorry');
  assert.equal(receipt.status, 'fail');
  assert.ok(receipt.failures.some(item => /source rule sorry/.test(item)));
  assert.ok(receipt.failures.some(item => /sorryAx/.test(item)), receipt.failures.join('\n'));
});

test('lean gate: a custom axiom is rejected', { skip, timeout: 1_800_000 }, async () => {
  const receipt = await gate('CustomAxiom');
  assert.equal(receipt.status, 'fail');
  assert.ok(receipt.failures.some(item => /axiom-declaration/.test(item)));
  assert.ok(receipt.failures.some(item => /euler_fails/.test(item)), receipt.failures.join('\n'));
});

test('lean gate: native_decide is rejected', { skip, timeout: 1_800_000 }, async () => {
  const receipt = await gate('Native');
  assert.equal(receipt.status, 'fail');
  assert.ok(receipt.failures.some(item => /native-decide/.test(item)));
  assert.ok(receipt.failures.some(item => /ofReduceBool|trustCompiler/.test(item)), receipt.failures.join('\n'));
});

test('lean gate: a theorem of a weaker statement is rejected', { skip, timeout: 1_800_000 }, async () => {
  const receipt = await gate('Weakened');
  assert.equal(receipt.status, 'fail');
  assert.ok(receipt.failures.some(item => /must be exactly ¬ Fixtures\.Weakened\.Claim/.test(item)), receipt.failures.join('\n'));
});

test('lean gate: a silently weakened statement breaks the pinned fingerprint', { skip, timeout: 1_800_000 }, async () => {
  const pinned = (await gate('Good')).statement.fingerprint;
  const copy = mkdtempSync(join(tmpdir(), 'discover-tamper-'));
  cpSync(project, copy, { recursive: true, filter: path => !path.includes('.lake') });
  cpSync(join(project, 'tampered', 'Good.lean'), join(copy, 'Fixtures', 'Good.lean'));
  const receipt = await leanGate({ project: copy, modules: ['Fixtures.Good'], claim: 'Fixtures.Good.Claim', main: 'Fixtures.Good.main', polarity: 'disprove', expectFingerprint: pinned, buildTimeoutMs: 20 * 60_000 });
  assert.equal(receipt.status, 'fail');
  assert.ok(receipt.failures.some(item => /fingerprint .* differs/.test(item)), receipt.failures.join('\n'));
  assert.equal(receipt.mainTheorem.refutesClaim, true, 'the tampered proof itself is fine; only the pin catches it');
});
