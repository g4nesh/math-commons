#!/usr/bin/env node
import { open, readFile } from 'node:fs/promises';
import { constants } from 'node:fs';
import { randomUUID } from 'node:crypto';
import { resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { validateMilestone } from '../lib/worker-progress.mjs';
import { GitHubStore, resolveToken } from '../lib/github.mjs';

export async function publishMilestone({ summary, details, evidence = [], id = randomUUID() }, path = process.env.MATH_COMMONS_PROGRESS_FILE) {
  if (!path) throw new Error('No worker progress channel. Start this command inside a supported commons work run.');
  const frame = validateMilestone({ version: 1, type: 'checkpoint', id, summary, details, evidence });
  // Do not create arbitrary files if the worker channel is missing or stale.
  const file = await open(path, constants.O_WRONLY | constants.O_APPEND | (constants.O_NOFOLLOW ?? 0));
  try {
    const stat = await file.stat();
    if (!stat.isFile() || stat.size > 1_048_576) throw new Error('Worker progress channel is unavailable or full.');
    // One append-only write; never create a stale path or follow a symlink.
    await file.write(Buffer.from(JSON.stringify(frame) + '\n'));
  } finally { await file.close(); }
  return frame.id;
}

export async function replayMilestone(path, store) {
  const file = await open(path, 'r');
  let data;
  try { if ((await file.stat()).size > 40_000) throw new Error('Recovery record too large.'); data = JSON.parse(await file.readFile('utf8')); } finally { await file.close(); }
  const op = data?.operation;
  if (data?.schemaVersion !== 1 || data.repository !== store.config.repository || op?.kind !== 'checkpoint' ||
      !/^progress-[a-f0-9]{64}$/.test(op.id || '') || !/^[a-z0-9][a-z0-9-]{0,79}$/.test(op.agent || '') || !/^[a-z0-9][a-z0-9-]{0,79}$/.test(op.task || '') ||
      Object.keys(op).some(key => !['id', 'kind', 'agent', 'task', 'summary', 'details', 'evidence'].includes(key))) throw new Error('Invalid milestone recovery record.');
  validateMilestone({ version: 1, type: 'checkpoint', id: op.id, summary: op.summary, details: op.details, evidence: op.evidence });
  return store.mutate(op);
}

if (process.argv[1] && resolve(process.argv[1]) === fileURLToPath(import.meta.url)) {
  try {
    const options = {}, args = process.argv.slice(2);
    if (args.length === 2 && args[0] === '--replay') {
      const config = JSON.parse(await readFile(new URL('../commons.config.json', import.meta.url), 'utf8'));
      const result = await replayMilestone(args[1], new GitHubStore(config, { token: resolveToken() }));
      process.stdout.write(`Public milestone reconciled: ${result.requestId}.\n`);
      process.exit(0);
    }
    for (let i = 0; i < args.length; i += 2) {
      const name = args[i];
      if (!['--summary', '--details', '--details-file', '--evidence', '--id'].includes(name) || args[i + 1] === undefined || options[name] !== undefined) throw new Error('Use --summary TEXT --details TEXT [--evidence URL,URL] [--id STABLE-ID].');
      options[name] = args[i + 1];
    }
    if (options['--details'] && options['--details-file']) throw new Error('Choose --details or --details-file.');
    const details = options['--details-file'] ? await readFile(options['--details-file'], 'utf8') : options['--details'];
    const id = await publishMilestone({ summary: options['--summary'], details,
      evidence: options['--evidence'] ? options['--evidence'].split(',') : [], ...(options['--id'] ? { id: options['--id'] } : {}) });
    process.stdout.write(`Public milestone queued: ${id}. The worker publishes it under its current lease.\n`);
  } catch {
    process.stderr.write('Milestone was not queued. Check the worker channel and documented public checkpoint fields.\n');
    process.exitCode = 1;
  }
}
