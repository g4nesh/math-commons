import { existsSync, readdirSync } from 'node:fs';
import { join } from 'node:path';
import { HARNESS_ROOT, readText, writeText } from './util.mjs';

// The canonical runtime files live in harness/runtime/. They are copied to the
// locations each runtime reads: Claude Code (.claude/skills, .claude/agents)
// and Codex (.agents/skills, .codex/agents). `--check` fails on drift so the
// copies never silently diverge from the source.
export function runtimeTargets(root) {
  const source = join(HARNESS_ROOT, 'runtime');
  const targets = [];
  for (const skill of readdirSync(join(source, 'skills'))) {
    const file = join(source, 'skills', skill, 'SKILL.md');
    if (!existsSync(file)) continue;
    targets.push({ from: file, to: join(root, '.claude', 'skills', skill, 'SKILL.md') });
    targets.push({ from: file, to: join(root, '.agents', 'skills', skill, 'SKILL.md') });
  }
  for (const name of readdirSync(join(source, 'claude-agents')).filter(name => name.endsWith('.md'))) {
    targets.push({ from: join(source, 'claude-agents', name), to: join(root, '.claude', 'agents', name) });
  }
  for (const name of readdirSync(join(source, 'codex-agents')).filter(name => name.endsWith('.toml'))) {
    targets.push({ from: join(source, 'codex-agents', name), to: join(root, '.codex', 'agents', name) });
  }
  return targets;
}

export function syncRuntime(root, { check = false } = {}) {
  const drift = [];
  for (const { from, to } of runtimeTargets(root)) {
    const text = readText(from);
    if (existsSync(to) && readText(to) === text) continue;
    drift.push(to);
    if (!check) writeText(to, text);
  }
  return drift;
}
