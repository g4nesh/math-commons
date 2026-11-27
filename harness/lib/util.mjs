import { createHash } from 'node:crypto';
import { existsSync, readFileSync, writeFileSync, mkdirSync, readdirSync, statSync, appendFileSync, renameSync } from 'node:fs';
import { dirname, join, relative, resolve, sep } from 'node:path';
import { fileURLToPath } from 'node:url';

export const HARNESS_ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '..');
export const REPO_ROOT = resolve(HARNESS_ROOT, '..');

export class HarnessError extends Error {
  constructor(code, message) { super(message); this.code = code; }
}

export const sha256 = data => createHash('sha256').update(data).digest('hex');
export const fileSha256 = path => sha256(readFileSync(path));
export const short = hash => (hash || '').slice(0, 12);
export const now = () => new Date().toISOString();

// Rough but stable: math-heavy Markdown averages close to four characters per token
// for both providers. Budgets are planning aids, not billing.
export const estimateTokens = text => Math.ceil(Buffer.byteLength(text || '', 'utf8') / 4);

export function readText(path, fallback = null) {
  if (!existsSync(path)) {
    if (fallback !== null) return fallback;
    throw new HarnessError('MISSING_FILE', `Missing file: ${path}`);
  }
  return readFileSync(path, 'utf8');
}

export function readJson(path, fallback = null) {
  if (!existsSync(path)) {
    if (fallback !== null) return fallback;
    throw new HarnessError('MISSING_FILE', `Missing file: ${path}`);
  }
  try { return JSON.parse(readFileSync(path, 'utf8')); }
  catch (error) { throw new HarnessError('BAD_JSON', `Invalid JSON in ${path}: ${error.message}`); }
}

// Write through a temporary file so an interrupted run never leaves half a state file.
export function writeText(path, text) {
  mkdirSync(dirname(path), { recursive: true });
  const temporary = `${path}.tmp-${process.pid}`;
  writeFileSync(temporary, text);
  renameSync(temporary, path);
}

export const writeJson = (path, value) => writeText(path, JSON.stringify(value, null, 2) + '\n');

export function appendJsonl(path, value) {
  mkdirSync(dirname(path), { recursive: true });
  appendFileSync(path, JSON.stringify(value) + '\n');
}

export function readJsonl(path) {
  if (!existsSync(path)) return [];
  return readFileSync(path, 'utf8').split('\n').filter(line => line.trim()).map((line, index) => {
    try { return JSON.parse(line); }
    catch { throw new HarnessError('BAD_JSONL', `Invalid JSON on line ${index + 1} of ${path}`); }
  });
}

export function listFiles(root, { ignore = ['.git', '.lake', 'node_modules', '.runs', '__pycache__', '.harness'] } = {}) {
  if (!existsSync(root)) return [];
  const out = [];
  const walk = directory => {
    for (const entry of readdirSync(directory, { withFileTypes: true })) {
      if (ignore.includes(entry.name)) continue;
      const path = join(directory, entry.name);
      if (entry.isDirectory()) walk(path);
      else if (entry.isFile()) out.push(path);
    }
  };
  walk(root);
  return out.sort();
}

export const rel = (from, path) => relative(from, path).split(sep).join('/');
export const isInside = (parent, child) => {
  const path = relative(resolve(parent), resolve(child));
  return path === '' || (!path.startsWith('..') && !path.startsWith(sep) && !/^[a-zA-Z]:/.test(path));
};

export const fileSize = path => (existsSync(path) ? statSync(path).size : 0);

export function slugify(text) {
  return String(text).toLowerCase().normalize('NFKD').replace(/[̀-ͯ]/g, '')
    .replace(/[^a-z0-9]+/g, '-').replace(/^-+|-+$/g, '').slice(0, 64);
}

export function assertSlug(slug) {
  if (!/^[a-z0-9][a-z0-9-]{1,63}$/.test(slug || '')) {
    throw new HarnessError('BAD_SLUG', `Campaign slug must be lowercase letters, digits and hyphens: ${slug}`);
  }
  return slug;
}

// Stable JSON for hashing: sorted keys at every level.
export function canonicalJson(value) {
  if (Array.isArray(value)) return `[${value.map(canonicalJson).join(',')}]`;
  if (value && typeof value === 'object') {
    return `{${Object.keys(value).sort().map(key => `${JSON.stringify(key)}:${canonicalJson(value[key])}`).join(',')}}`;
  }
  return JSON.stringify(value);
}

export function truncate(text, max, marker = ' […]') {
  if (!text || text.length <= max) return text || '';
  return text.slice(0, Math.max(0, max - marker.length)).trimEnd() + marker;
}

export function parseArgs(argv) {
  const positional = [], flags = {};
  for (let index = 0; index < argv.length; index++) {
    const token = argv[index];
    if (token === '--') { flags['--'] = argv.slice(index + 1); break; }
    if (token.startsWith('--')) {
      const [key, inline] = token.slice(2).split(/=(.*)/s, 2);
      if (inline !== undefined) flags[key] = inline;
      else if (index + 1 < argv.length && !argv[index + 1].startsWith('--')) flags[key] = argv[++index];
      else flags[key] = true;
    } else positional.push(token);
  }
  return { positional, flags };
}
