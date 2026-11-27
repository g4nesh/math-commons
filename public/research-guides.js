// Shared by the public site and both local worker adapters.
export const fourierObjective = 'fourier-principal-minors';
const fourierRoot = 'research/' + fourierObjective;
const fourierPages = new Map([
  ['fourier-154-replay', 'order-154/replay.md'],
  ['fourier-154-prior-art', 'order-154/prior-art.md'],
  ['fourier-154-synthesis', 'order-154/synthesis.md'],
  ['fourier-14p-audit', 'family-14p/README.md']
]);

export function taskGuide(objective, task) {
  const page = objective === fourierObjective && fourierPages.get(task);
  return page ? fourierRoot + '/' + page : null;
}

export function objectiveGuides(id, config = {}) {
  const known = id === fourierObjective || id === 'marcus-de-oliveira';
  const root = known ? 'research/' + id : null;
  const active = id && id === config.defaultObjectiveId;
  return {
    overview: (active && config.defaultObjectiveOverview) || (root && root + '/overview.md'),
    brief: (active && config.defaultObjectiveBrief) || (root && root + '/README.md')
  };
}
