import { resolve } from 'node:path';
import { HarnessError, estimateTokens, readText, isInside } from './util.mjs';

// Locators are the harness's retrieval keys: `path`, `path#heading-anchor` or
// `path:L10-L40`. Prompts carry locators; agents expand only what they need.

export function anchorOf(heading) {
  return heading.trim().toLowerCase().replace(/[`*_~]/g, '')
    .replace(/[^\p{L}\p{N}\s-]/gu, '').trim().replace(/\s+/g, '-');
}

export function outline(text) {
  const lines = text.split('\n');
  const headings = [];
  let fence = false;
  lines.forEach((line, index) => {
    if (/^\s*(```|~~~)/.test(line)) fence = !fence;
    if (fence) return;
    const match = /^(#{1,6})\s+(.+?)\s*#*\s*$/.exec(line);
    if (match) headings.push({ level: match[1].length, title: match[2], anchor: anchorOf(match[2]), line: index + 1 });
  });
  // Section extent: up to the next heading of the same or higher level.
  headings.forEach((heading, index) => {
    const next = headings.slice(index + 1).find(other => other.level <= heading.level);
    heading.endLine = next ? next.line - 1 : lines.length;
    heading.tokens = estimateTokens(lines.slice(heading.line - 1, heading.endLine).join('\n'));
  });
  return headings;
}

export function parseLocator(locator) {
  const lines = /^(.*):L(\d+)(?:-L?(\d+))?$/.exec(locator);
  if (lines) return { path: lines[1], from: Number(lines[2]), to: Number(lines[3] || lines[2]) };
  const hash = locator.indexOf('#');
  if (hash === -1) return { path: locator };
  return { path: locator.slice(0, hash), anchor: locator.slice(hash + 1) };
}

export function extract(text, { anchor, from, to }) {
  const lines = text.split('\n');
  if (from) return lines.slice(from - 1, to).join('\n');
  if (!anchor) return text;
  const heading = outline(text).find(entry => entry.anchor === anchor || entry.anchor === anchorOf(anchor));
  if (!heading) throw new HarnessError('MISSING_SECTION', `No heading #${anchor}`);
  return lines.slice(heading.line - 1, heading.endLine).join('\n').trimEnd();
}

export function resolveLocator(root, locator) {
  const parsed = parseLocator(locator);
  const path = resolve(root, parsed.path);
  if (!isInside(root, path)) throw new HarnessError('OUTSIDE_ROOT', `Locator leaves the campaign: ${locator}`);
  return { ...parsed, absolute: path, text: extract(readText(path), parsed) };
}

// A digest is the producer's own ≤10-line abstract at the top of a long artifact.
// Context packs send it in place of the full file when the budget is tight.
export function digestOf(text) {
  const match = /<!--\s*digest\s*-->([\s\S]*?)<!--\s*\/digest\s*-->/i.exec(text);
  return match ? match[1].trim() : null;
}

export function outlineText(text, path) {
  return outline(text).map(entry => `${'  '.repeat(entry.level - 1)}- ${entry.title} → ${path}#${entry.anchor} (L${entry.line}-L${entry.endLine}, ~${entry.tokens} tok)`).join('\n');
}

// Parse the first Markdown table under a heading into row objects keyed by header.
export function tableRows(text) {
  const rows = [];
  let header = null;
  for (const line of text.split('\n')) {
    if (!/^\s*\|/.test(line)) {
      if (header) break;
      continue;
    }
    const cells = line.trim().replace(/^\||\|$/g, '').split(/(?<!\\)\|/).map(cell => cell.trim());
    if (!header) { header = cells.map(cell => cell.toLowerCase()); continue; }
    if (cells.every(cell => /^:?-{2,}:?$/.test(cell))) continue;
    rows.push(Object.fromEntries(header.map((key, index) => [key, cells[index] ?? ''])));
  }
  return rows;
}
