import { loadSnapshot } from './state.js';
import { renderContributionsSvg, renderUnavailableSvg } from '../lib/contributions-svg.mjs';

export function createHandler({ load = loadSnapshot, now = Date.now } = {}) {
  return async function handler(req, res) {
    res.setHeader('Content-Type', 'image/svg+xml; charset=utf-8');
    res.setHeader('X-Content-Type-Options', 'nosniff');
    res.setHeader('Content-Security-Policy', "default-src 'none'; style-src 'none'; sandbox");
    if (req.method !== 'GET') {
      res.setHeader('Allow', 'GET');
      res.setHeader('Cache-Control', 'no-store');
      return res.status(405).end(renderUnavailableSvg({ method: true }));
    }
    try {
      const url = new URL(req.url || '/api/contributions', 'https://math-commons-dun.vercel.app');
      const demo = req.query?.demo === '1' || url.searchParams.get('demo') === '1';
      const snapshot = demo ? null : await load();
      const svg = renderContributionsSvg(snapshot, { now: now(), demo });
      // Separate URLs cache separately. README uses the live URL without demo=1.
      res.setHeader('Cache-Control', 'public, max-age=0, s-maxage=900');
      return res.status(200).end(svg);
    } catch {
      // A valid image keeps GitHub's embed legible during upstream failures.
      // It states the outage rather than turning missing data into a zero chart.
      res.setHeader('Cache-Control', 'no-store');
      return res.status(200).end(renderUnavailableSvg());
    }
  };
}
export default createHandler();
