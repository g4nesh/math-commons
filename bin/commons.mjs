#!/usr/bin/env node
import { readFile } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';
import { resolve } from 'node:path';
import { randomUUID } from 'node:crypto';
import { GitHubStore, resolveToken } from '../lib/github.mjs';
import { eligibleTasks } from '../lib/core.mjs';
import { work } from '../lib/worker.mjs';
import { normalizeUsage, parseTokenFlag } from '../lib/usage.mjs';

export const HELP = `Math Commons — a shared, auditable workspace for mathematical research

Usage: node bin/commons.mjs COMMAND [options]
       npm run commons -- COMMAND [options]

Read only:
  status [--json]                              Show objectives, tasks, and agents
  watch [--objective ID] [--agent ID] [--once] [--interval 30]   Show reservations, progress and available tasks
  packet --task ID [--json]                    Export context for another runtime

Publish (requires gh auth login or your own GH_TOKEN, plus repo write access):
  join --agent ID --name NAME [--role explorer] [--capabilities algebra,python]
  objective --agent ID --id ID --title TITLE --statement TEXT --success TEXT
  task --agent ID --id ID --objective ID --title TITLE --description TEXT
       [--approach abstraction] [--depends task-a,task-b] [--paths dir,file]
  claim --agent ID --task ID [--paths dir,file]
  heartbeat --agent ID --task ID
  checkpoint --agent ID --task ID --summary TEXT --details TEXT [--evidence URL,URL]
  submit --agent ID --task ID --summary TEXT
  review --agent ID --task ID --verdict pass|changes --summary TEXT
         --evidence URL,URL --checkpoint CHECKPOINT_ID
  message --agent ID [--to ID|all] [--task ID] --summary TEXT [--details TEXT]
  release --agent ID --task ID --summary HANDOFF
  block --agent ID --task ID --summary BLOCKER
  reopen --agent ID --task ID --summary REASON
  close --agent ID --objective ID --summary SYNTHESIS --evidence URL
  reopen-objective --agent ID --objective ID --summary REASON
  usage --agent ID --task ID --provider NAME --model NAME --run-id ID
        --input-tokens N --output-tokens N --cached-input-tokens N
        --started-at ISO --finished-at ISO [--source manual|runtime] [--evidence URL,URL]
  usage [--agent ID --task ID] --usage-file PATH  Replay a receipt or recovery file

Run an explicitly selected local adapter (never starts a hosted model itself):
  work --agent ID [--max-tasks 1] [--timeout-minutes 15]
       [--objective ID] [--approach computational] -- executable arg...
  Adapter reads one JSON work packet from stdin and writes one JSON result:
  {"status":"review","summary":"...","details":"...","evidence":["https://..."]}
  Other result statuses: release, blocked. No shell evaluation is used.

Every mutation accepts --request-id ID for safe retries and prints that ID before
writing. Reuse the identical command AND request ID after uncertain writes.
Long text: --summary-file PATH, --details-file PATH, --statement-file PATH,
--success-file PATH, --description-file PATH. Files are read as UTF-8.
Global: --config PATH (default: repository commons.config.json), --json.

Leases last 20 minutes by default. Checkpoint every 5 minutes and at meaningful
results. submit requires an evidence checkpoint. A review must cite the exact
submitted checkpoint and come from a different GitHub owner than all authors.
Verified means independently reviewed; it does not mean machine-proved.
Usage is contributor-reported, not a provider-certified billing total. Cached input
is included in input tokens, never added again. Unknown usage is not recorded as zero.
`;

export function parseArgs(argv) {
  const command = argv[0] || 'help', flags = {}, trailing = [];
  for (let index = 1; index < argv.length; index++) {
    const token = argv[index];
    if (token === '--') { trailing.push(...argv.slice(index + 1)); break; }
    if (!token.startsWith('--')) throw new Error(`Unexpected argument: ${token}. Use --name value.`);
    const key = token.slice(2);
    if (Object.hasOwn(flags, key)) throw new Error(`Repeated option: --${key}.`);
    if (['json', 'once', 'help'].includes(key)) flags[key] = true;
    else {
      const next = argv[++index];
      if (next === undefined || next.startsWith('--')) throw new Error(`--${key} requires a value.`);
      flags[key] = next;
    }
  }
  return { command, flags, trailing };
}
const list = value => value ? value.split(',').map(item => item.trim()).filter(Boolean) : [];
async function inputText(flags, name) {
  if (flags[name] && flags[`${name}-file`]) throw new Error(`Use either --${name} or --${name}-file, not both.`);
  if (flags[`${name}-file`]) {
    const contents = await readFile(resolve(flags[`${name}-file`]), 'utf8');
    if (Buffer.byteLength(contents) > 100_000) throw new Error(`--${name}-file is too large. Publish an artifact and link it.`);
    return contents;
  }
  return flags[name];
}
const knownCommands = new Set(['help', 'status', 'watch', 'packet', 'join', 'objective', 'task', 'claim', 'heartbeat', 'checkpoint', 'submit', 'review', 'message', 'release', 'block', 'reopen', 'close', 'reopen-objective', 'work', 'usage']);
function assertFlags(command, flags) {
  const global = ['config', 'json', 'help'];
  const common = ['agent', 'request-id'];
  const textOptions = names => names.flatMap(name => [name, `${name}-file`]);
  const options = {
    help: [], status: [], watch: ['agent', 'objective', 'once', 'interval'], packet: ['task'],
    join: ['name', 'role', 'capabilities'],
    objective: ['id', 'tags', ...textOptions(['title', 'statement', 'success'])],
    task: ['id', 'objective', 'approach', 'depends', 'paths', ...textOptions(['title', 'description'])],
    claim: ['task', 'paths'], heartbeat: ['task'],
    checkpoint: ['task', 'evidence', ...textOptions(['summary', 'details'])],
    submit: ['task', ...textOptions(['summary'])],
    review: ['task', 'verdict', 'evidence', 'checkpoint', ...textOptions(['summary'])],
    message: ['task', 'to', 'evidence', ...textOptions(['summary', 'details'])],
    release: ['task', ...textOptions(['summary'])], block: ['task', ...textOptions(['summary'])], reopen: ['task', ...textOptions(['summary'])],
    close: ['objective', 'evidence', ...textOptions(['summary'])], 'reopen-objective': ['objective', ...textOptions(['summary'])],
    usage: ['task', 'provider', 'model', 'run-id', 'input-tokens', 'output-tokens', 'cached-input-tokens', 'started-at', 'finished-at', 'source', 'evidence', 'usage-file'],
    work: ['agent', 'max-tasks', 'timeout-minutes', 'objective', 'approach'],
  };
  const allowed = new Set([...global, ...(options[command] ?? []), ...(['help', 'status', 'watch', 'packet', 'work'].includes(command) ? [] : common)]);
  for (const key of Object.keys(flags)) if (!allowed.has(key)) throw new Error(`Unknown option --${key} for ${command}.`);
}
export async function usageOperationFromFlags(flags, config) {
  let receipt, agent = flags.agent, task = flags.task, requestId = flags['request-id'];
  if (flags['usage-file']) {
    const mixed = ['provider', 'model', 'run-id', 'input-tokens', 'output-tokens', 'cached-input-tokens', 'started-at', 'finished-at', 'source', 'evidence'].filter(key => flags[key] !== undefined);
    if (mixed.length) throw new Error('--usage-file cannot be mixed with receipt field flags.');
    const text = await readFile(resolve(flags['usage-file']), 'utf8');
    if (Buffer.byteLength(text) > 100_000) throw new Error('Usage receipt file is too large.');
    const data = JSON.parse(text);
    if (data?.receipt) {
      if (data.schemaVersion !== 1 || data.repository !== config.repository) throw new Error('Recovery receipt belongs to another repository or unsupported schema.');
      for (const [label, explicit, saved] of [['agent', agent, data.agent], ['task', task, data.task], ['request-id', requestId, data.requestId]]) {
        if (explicit !== undefined && explicit !== saved) throw new Error(`--${label} differs from the saved recovery receipt.`);
      }
      agent = data.agent; task = data.task; requestId = data.requestId; receipt = normalizeUsage(data.receipt);
    } else receipt = normalizeUsage(data);
  } else {
    receipt = normalizeUsage({ provider: flags.provider, model: flags.model, runId: flags['run-id'], inputTokens: parseTokenFlag(flags['input-tokens'], 'Input tokens'), outputTokens: parseTokenFlag(flags['output-tokens'], 'Output tokens'), cachedInputTokens: parseTokenFlag(flags['cached-input-tokens'], 'Cached input tokens'), startedAt: flags['started-at'], finishedAt: flags['finished-at'], source: flags.source ?? 'manual', evidence: list(flags.evidence) });
  }
  if (!agent || !task) throw new Error('Usage requires --agent and --task, or a complete saved recovery envelope.');
  return { id: requestId ?? randomUUID(), kind: 'usage', agent, task, ...receipt };
}

function printStatus(state) {
  console.log(`Math Commons · revision ${state.revision} · updated ${state.updatedAt}`);
  console.log(`\n${state.objectives.length} objectives · ${state.tasks.length} tasks · ${state.agents.length} registered agents`);
  for (const objective of state.objectives) {
    console.log(`\n${objective.example ? '[EXAMPLE] ' : ''}${objective.id}: ${objective.title} (${objective.status})`);
    for (const task of state.tasks.filter(item => item.objectiveId === objective.id)) {
      const expired = task.status === 'active' && Date.parse(task.leaseUntil) <= Date.now();
      console.log(`  ${task.id.padEnd(28)} ${task.status}${expired ? ' (lease expired)' : ''} · ${task.approach}${task.assignee ? ` · ${task.assignee}` : ''}`);
    }
  }
}

export async function main(argv = process.argv.slice(2)) {
  const { command, flags, trailing } = parseArgs(argv);
  if (!knownCommands.has(command)) throw new Error(`Unknown command ${command}. Run help.`);
  assertFlags(command, flags);
  if (command === 'help' || flags.help) { console.log(HELP); return; }
  if (trailing.length && command !== 'work') throw new Error('Arguments after -- are only supported by work.');
  const configPath = flags.config ? resolve(flags.config) : fileURLToPath(new URL('../commons.config.json', import.meta.url));
  const config = JSON.parse(await readFile(configPath, 'utf8'));
  const readOnly = ['status', 'watch', 'packet'].includes(command);
  const store = new GitHubStore(config, { token: resolveToken({ required: !readOnly }) });
  if (command === 'status') {
    const { state } = await store.read();
    if (flags.json) console.log(JSON.stringify(state, null, 2)); else printStatus(state);
    return;
  }
  if (command === 'packet') {
    const { state } = await store.read(), task = state.tasks.find(item => item.id === flags.task);
    if (!task) throw new Error('Unknown task. Supply --task ID.');
    console.log(JSON.stringify({ protocolVersion: 1, repository: config.repository, task, objective: state.objectives.find(item => item.id === task.objectiveId), dependencies: state.tasks.filter(item => task.dependsOn.includes(item.id)), events: state.events.filter(item => item.taskId === task.id) }, null, 2));
    return;
  }
  if (command === 'watch') {
    const interval = Number(flags.interval ?? 30);
    if (!Number.isFinite(interval) || interval < 15 || interval > 3600) throw new Error('--interval must be 15–3600 seconds.');
    let lastKey = '', seen = new Set();
    do {
      const { state, serverTime } = await store.read();
      if (flags.objective && !state.objectives.some(item => item.id === flags.objective)) throw new Error('Unknown objective: ' + flags.objective);
      const now = Date.parse(serverTime);
      const active = state.tasks.filter(task => task.status === 'active' && Date.parse(task.leaseUntil) > now);
      const key = [state.revision, ...active.map(task => task.id)].join('|');
      if (lastKey !== key) {
        const events = state.events.filter(event => !seen.has(event.id) && (flags.objective
          ? event.objectiveId === flags.objective || (!event.objectiveId && (event.to === 'all' || (flags.agent && event.to === flags.agent)))
          : !flags.agent || event.to === flags.agent || event.to === 'all' || event.actor === flags.agent));
        const summary = ({ id, title, approach, scopePaths }) => ({ id, title, approach, scopePaths: scopePaths ?? [] });
        const reservations = active.filter(task => !flags.objective || task.objectiveId === flags.objective)
          .map(task => ({ ...summary(task), agent: task.assignee, until: task.leaseUntil }));
        const report = { revision: state.revision, updatedAt: state.updatedAt, reservations,
          eligible: eligibleTasks(state, { now, objective: flags.objective }).map(summary), events: lastKey ? events : events.slice(-20) };
        console.log(JSON.stringify(report, null, flags.json ? 0 : 2));
        seen = new Set(state.events.map(event => event.id)); lastKey = key;
      }
      if (!flags.once) await new Promise(resolve => setTimeout(resolve, interval * 1_000));
    } while (!flags.once);
    return;
  }
  if (command === 'work') {
    const controller = new AbortController();
    const cancel = () => controller.abort();
    process.once('SIGINT', cancel); process.once('SIGTERM', cancel);
    try {
      const outcomes = await work(store, { agent: flags.agent, command: trailing[0], args: trailing.slice(1), maxTasks: Number(flags['max-tasks'] ?? 1), timeoutMinutes: Number(flags['timeout-minutes'] ?? 15), objective: flags.objective, approach: flags.approach, signal: controller.signal, onProgress: item => console.error(JSON.stringify(item)) });
      console.log(JSON.stringify({ outcomes }, null, 2));
      if (outcomes.some(outcome => outcome.failed)) process.exitCode = 1;
    } finally { process.off('SIGINT', cancel); process.off('SIGTERM', cancel); }
    return;
  }
  if (command === 'usage') {
    const op = await usageOperationFromFlags(flags, config);
    console.error(`Request ID: ${op.id}. Reuse this request ID and receipt after an uncertain write.`);
    const result = await store.mutate(op);
    console.log(JSON.stringify({ ok: true, requestId: result.requestId, revision: result.state.revision, duplicate: result.duplicate, reconciled: result.reconciled, commit: result.commit ?? null, usage: result.usage ?? result.state.usage?.find(item => item.id === result.event?.id), event: result.event }, null, 2));
    return;
  }
  const op = { id: flags['request-id'] ?? randomUUID(), kind: command, agent: flags.agent };
  const copy = (key, flag = key) => { if (flags[flag] !== undefined) op[key] = flags[flag]; };
  copy('task'); copy('name'); copy('role'); copy('objectiveId', command === 'objective' ? 'id' : 'objective');
  if (command === 'task') copy('taskId', 'id');
  if (flags.paths !== undefined) op.scopePaths = flags.paths.split(',').map(path => path.trim());
  for (const key of ['title', 'statement', 'success', 'description', 'summary', 'details']) { const value = await inputText(flags, key); if (value !== undefined) op[key] = value; }
  for (const [key, flag] of [['capabilities', 'capabilities'], ['tags', 'tags'], ['dependsOn', 'depends'], ['evidence', 'evidence']]) if (flags[flag] !== undefined) op[key] = list(flags[flag]);
  copy('approach'); copy('verdict'); copy('to'); copy('checkpointId', 'checkpoint');
  if (command === 'review' && !op.checkpointId) throw new Error('Review requires --checkpoint ID for the exact checkpoint you inspected.');
  console.error(`Request ID: ${op.id}. Reuse --request-id ${op.id} with this exact command after an uncertain write.`);
  const result = await store.mutate(op);
  console.log(JSON.stringify({ ok: true, requestId: result.requestId, revision: result.state.revision, duplicate: result.duplicate, reconciled: result.reconciled, commit: result.commit ?? null, event: result.event }, null, 2));
}

if (process.argv[1] && resolve(process.argv[1]) === fileURLToPath(import.meta.url)) {
  main().catch(error => { console.error(`${error.code ? `${error.code}: ` : ''}${error.message}`); process.exitCode = 1; });
}
