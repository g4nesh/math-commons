/** Compact, same-objective origin context for an explicitly requested task. */
export function requestSourceContext(state, task) {
  const id = task?.request?.fromTaskId;
  if (!id) return null;
  const origin = state.tasks.find(item => item.id === id && item.objectiveId === task.objectiveId && item.id !== task.id);
  if (!origin) return null;
  const checkpoint = origin.checkpoints?.at(-1);
  const evidence = (Array.isArray(checkpoint?.evidence) ? checkpoint.evidence : []).filter(value => {
    try { const url = new URL(value); return typeof value === 'string' && value.length <= 2000 && url.protocol === 'https:' && !url.username && !url.password; } catch { return false; }
  }).slice(0, 3);
  return { id: origin.id, title: String(origin.title || origin.id).slice(0, 240), latestCheckpoint: checkpoint ? {
    id: checkpoint.id, evidenceCount: checkpoint.evidence?.length || 0, summary: String(checkpoint.summary || '').slice(0, 500), evidence
  } : null };
}
