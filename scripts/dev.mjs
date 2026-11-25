import { createServer } from 'node:http';
import { readFile, stat } from 'node:fs/promises';
import { resolve, extname, sep } from 'node:path';
import { fileURLToPath } from 'node:url';
import handler from '../api/state.js';
const root = fileURLToPath(new URL('../dist/', import.meta.url));
const port = Number(process.env.PORT || 4317);
const mime = { '.html': 'text/html', '.css': 'text/css', '.js': 'text/javascript', '.svg': 'image/svg+xml', '.json': 'application/json', '.md': 'text/plain' };
const server = createServer(async (req, res) => {
  res.status = code => { res.statusCode = code; return res; };
  res.json = value => { res.setHeader('Content-Type', 'application/json'); res.end(JSON.stringify(value)); };
  try {
    const url = new URL(req.url, 'http://localhost');
    if (url.pathname === '/api/state') return await handler(req, res);
    if (req.method !== 'GET' && req.method !== 'HEAD') return res.status(405).end();
    const path = resolve(root, '.' + decodeURIComponent(url.pathname === '/' ? '/index.html' : url.pathname));
    if (!path.startsWith(root.endsWith(sep) ? root : root + sep)) return res.status(403).end();
    if (!(await stat(path)).isFile()) return res.status(404).end();
    res.setHeader('Content-Type', `${mime[extname(path)] || 'application/octet-stream'}; charset=utf-8`);
    res.setHeader('Cache-Control', 'no-store');
    res.end(req.method === 'HEAD' ? undefined : await readFile(path));
  } catch { res.status(404).end('Not found'); }
});
server.listen(port, '127.0.0.1', () => console.log(`Math Commons: http://127.0.0.1:${port}`));
