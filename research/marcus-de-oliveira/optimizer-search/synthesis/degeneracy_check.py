"""Degeneracy diagnostics for the near-contacts reported by the search lanes.

Reads saved run files only (no optimization). Run from this folder:
    python -I -B degeneracy_check.py            # writes degeneracy_check.json

For each configuration it reports, in canonical units (rms spread of a and b):
  width_over_diam   minimum-width proxy of H = conv{z_sigma} (spread along the
                    least principal axis of the z cloud) divided by diam(H)
  concyclic_res     smallest relative singular value of [|p|^2, Re p, Im p, 1]
                    over p in {a_i} U {-b_j}; ~0 iff a and -b lie on one circle
                    or line (then all z_sigma share one argument mod pi, H is
                    a segment and d lies on it: the Moebius image of the
                    Hermitian case)
  min_abs_aplusb    min |a_i + b_j| (zero factor)
  min_gap_a/b       smallest eigenvalue gap (collision)
  max_imag_b        max |Im b_j| (for real a: distance to the Hermitian case)
For the zero-factor contacts of joint5 it also tests the nested reduction
described in ../README.md (Section 6.2) and checks the identity
d = |U_ik|^2 (a_i + b_k) d' numerically on the 1 x 1 type.
"""
import json, itertools, os
from collections import Counter
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')


def cx(v):
    v = np.asarray(v, dtype=float)
    return v[:, 0] + 1j * v[:, 1]


def perms(n):
    return np.array(list(itertools.permutations(range(n))))


def zall(a, b, P):
    M = a[:, None] + b[None, :]
    return np.prod(M[np.arange(len(a))[None, :], P], axis=1)


def width_ratio(z):
    c = z - z.mean()
    X = np.stack([c.real, c.imag], 1)
    _, _, Vt = np.linalg.svd(X, full_matrices=False)
    p0, p1 = X @ Vt[0], X @ Vt[1]
    return float((p1.max() - p1.min()) / (p0.max() - p0.min()))


def circ_res(p):
    M = np.stack([np.abs(p) ** 2, p.real, p.imag, np.ones(len(p))], 1)
    M = M / np.linalg.norm(M, axis=0)
    s = np.linalg.svd(M, compute_uv=False)
    return float(s[-1] / s[0])


def diag(a, b, P):
    sc = np.sqrt((np.abs(a - a.mean()) ** 2).mean() + (np.abs(b - b.mean()) ** 2).mean())
    z = zall(a, b, P)
    gap = lambda v: min(abs(x - y) for x, y in itertools.combinations(v, 2)) / sc
    return dict(width_over_diam=width_ratio(z), concyclic_res=circ_res(np.concatenate([a, -b])),
                min_abs_aplusb=float(np.abs(a[:, None] + b[None, :]).min() / sc),
                min_gap_a=float(gap(a)), min_gap_b=float(gap(b)),
                max_imag_b=float(np.abs(b.imag).max() / sc))


def label(dg):
    tags = []
    if dg['width_over_diam'] < 1e-4 or dg['concyclic_res'] < 1e-4:
        tags.append('near-segment hull (a, -b near one circle or line)')
    if min(dg['min_gap_a'], dg['min_gap_b']) < 1e-3:
        tags.append('near eigenvalue collision')
    if dg['min_abs_aplusb'] < 1e-2:
        tags.append('near zero factor a_i + b_j = 0')
    return tags or ['none of the three']


out = {}

# 1. hermitian lane: tau-continued constrained optima
rows = []
for f in sorted(os.listdir(os.path.join(ROOT, 'hermitian', 'runs'))):
    if not f.startswith('cont_'):
        continue
    d = json.load(open(os.path.join(ROOT, 'hermitian', 'runs', f)))
    for i, r in enumerate(d['records']):
        a, b = cx(r['a']), cx(r['b'])
        dg = diag(a, b, perms(len(a)))
        rows.append(dict(file='hermitian/runs/' + f, record=i, cls=r.get('cls'), m=r['m'], n_active=r['n_active'],
                         dist_genperm=r['dist_genperm'], dist_blockdiag=r['dist_blockdiag'],
                         continuation_m=r.get('continuation_m'), **dg, tags=label(dg)))
out['hermitian_continuation'] = rows

# 2. joint67 constrained runs: best records
rows = []
for f in ['n6_constr_s602.json', 'n7_constr_s702.json']:
    d = json.load(open(os.path.join(ROOT, 'joint67', 'runs', f)))
    for i, r in enumerate(d['best'][:10]):
        a, b = cx(r['a']), cx(r['b'])
        dg = diag(a, b, perms(len(a)))
        rows.append(dict(file='joint67/runs/' + f, best_index=i, m=r['m'], dist_genperm=r['dist_genperm'],
                         dist_blockdiag=r['dist_blockdiag'], **dg, tags=label(dg)))
out['joint67_constrained'] = rows

# 3. joint5 zero-factor contacts: nested structure and the reduction identity
c = json.load(open(os.path.join(ROOT, 'joint5', 'contacts', 'classification.json')))
P5 = perms(5)
cache, rows, counts = {}, [], Counter()
for rec in c['records']:
    if not rec['cls'].startswith('zero'):
        continue
    f = rec['file']
    if f not in cache:
        cache[f] = json.load(open(os.path.join(ROOT, 'joint5', f)))
    r = cache[f]['best'][rec['best_index']]
    U = np.array(r['U_re']) + 1j * np.array(r['U_im'])
    a, b = cx(r['a']), cx(r['b'])
    d = np.linalg.det(np.diag(a) + U @ np.diag(b) @ U.conj().T)
    z = zall(a, b, P5)
    sc = np.abs(z).max()
    u = d / abs(d)
    on = (np.abs((z * np.conj(u)).imag) < 1e-9 * sc) & (np.abs(z) > 1e-9 * sc)
    S = P5[on]
    found = None
    for size in (1, 2):
        for I in itertools.combinations(range(5), size):
            Ks = set(tuple(sorted(s[list(I)])) for s in S)
            if len(Ks) == 1:
                found = ([int(x) for x in I], [int(x) for x in list(Ks)[0]])
                break
        if found:
            break
    row = dict(file='joint5/' + f, best_index=rec['best_index'], zero_pairs=rec['zero_pairs'],
               n_nonzero_z_on_support_line=int(on.sum()), n_small_U_entries=int((np.abs(U) ** 2 < 1e-12).sum()))
    if found:
        I, K = found
        rr = [x for x in range(5) if x not in I]
        cc = [x for x in range(5) if x not in K]
        row.update(block_rows=I, block_cols=K, complement_concyclic_res=circ_res(np.concatenate([a[rr], -b[cc]])))
    # reduction identity for a single zero pair (i, j) with row i supported on {j, k}
    # or column j supported on {i, m}
    if len(rec['zero_pairs']) == 1:
        i, j = int(rec['zero_pairs'][0][0]), int(rec['zero_pairs'][0][1])
        t = -a[i]
        a2, b2 = a + t, b - t                     # a2[i] = 0, b2[j] ~ 0
        small = np.abs(U) ** 2 < 1e-12
        rowsup = [l for l in range(5) if not small[i, l]]
        colsup = [l for l in range(5) if not small[l, j]]
        if len(rowsup) == 2 and j in rowsup:
            k = [l for l in rowsup if l != j][0]
            X = np.diag(a2) + sum(b2[l] * np.outer(U[:, l], U[:, l].conj()) for l in range(5) if l not in (j, k))
            dprime = np.linalg.det(np.delete(np.delete(X, i, 0), i, 1))
            pred = abs(U[i, k]) ** 2 * b2[k] * dprime
            row.update(reduction='row %d supported on columns {%d,%d}' % (i, j, k),
                       reduction_rel_error=float(abs(pred - d) / abs(d)))
        elif len(colsup) == 2 and i in colsup:
            m = [l for l in colsup if l != i][0]
            V = U.conj().T
            X = np.diag(b2) + sum(a2[l] * np.outer(V[:, l], V[:, l].conj()) for l in range(5) if l not in (i, m))
            dprime = np.linalg.det(np.delete(np.delete(X, j, 0), j, 1))
            pred = abs(U[m, j]) ** 2 * a2[m] * dprime
            row.update(reduction='column %d supported on rows {%d,%d}' % (j, i, m),
                       reduction_rel_error=float(abs(pred - d) / abs(d)))
    key = 'pairs=%d, nonzero z on line=%d, block=%s, complement concyclic<1e-8=%s' % (
        len(rec['zero_pairs']), on.sum(), ('%dx%d' % (len(found[0]), len(found[0]))) if found else 'none',
        (row.get('complement_concyclic_res', 1) < 1e-8))
    counts[key] += 1
    rows.append(row)
out['joint5_zero_factor'] = dict(counts=dict(counts), records=rows)

# 4. joint67 n = 6 near-contacts with a zero factor: same nested-structure test
rows = []
P6 = perms(6)
for r in json.load(open(os.path.join(ROOT, 'joint67', 'runs', 'n6_constr_s602.json')))['near_contacts']:
    a, b = cx(r['a']), cx(r['b'])
    U = np.array(r['U_re']) + 1j * np.array(r['U_im'])
    sc = np.sqrt((np.abs(a - a.mean()) ** 2).mean() + (np.abs(b - b.mean()) ** 2).mean())
    mab = np.abs(a[:, None] + b[None, :])
    if mab.min() / sc > 1e-6:
        continue
    i, j = (int(x) for x in np.unravel_index(mab.argmin(), mab.shape))
    d = np.linalg.det(np.diag(a) + U @ np.diag(b) @ U.conj().T)
    z = zall(a, b, P6)
    sc_z = np.abs(z).max()
    on = (np.abs((z * np.conj(d / abs(d))).imag) < 1e-9 * sc_z) & (np.abs(z) > 1e-9 * sc_z)
    S = P6[on]
    row = dict(restart=r.get('restart'), zero_pair=[i, j], min_abs_aplusb=float(mab.min() / sc),
               n_nonzero_z_on_support_line=int(on.sum()))
    for size in (1, 2, 3):
        hit = None
        for I in itertools.combinations(range(6), size):
            Ks = set(tuple(sorted(s[list(I)])) for s in S)
            if len(Ks) == 1:
                hit = ([int(x) for x in I], [int(x) for x in list(Ks)[0]])
                break
        if hit:
            I, K = hit
            rr = [x for x in range(6) if x not in I]
            cc = [x for x in range(6) if x not in K]
            row.update(block_rows=I, block_cols=K, complement_concyclic_res=circ_res(np.concatenate([a[rr], -b[cc]])))
            break
    rows.append(row)
out['joint67_n6_zero_factor'] = rows
print('joint67 n = 6 zero-factor near-contacts:', rows)

json.dump(out, open(os.path.join(HERE, 'degeneracy_check.json'), 'w'), indent=1)
for k in ('hermitian_continuation', 'joint67_constrained'):
    print(k)
    for r in out[k]:
        print('  %-40s m=%.3g wid=%.1e conc=%.1e |a+b|=%.1e gaps=%.1e,%.1e Imb=%.1e  %s' % (
            r['file'].split('/')[-1] + ' #' + str(r.get('record', r.get('best_index'))), r['m'], r['width_over_diam'],
            r['concyclic_res'], r['min_abs_aplusb'], r['min_gap_a'], r['min_gap_b'], r['max_imag_b'], '; '.join(r['tags'])))
print('joint5 zero-factor:')
for k, v in counts.items():
    print('  %3d  %s' % (v, k))
red = [r for r in out['joint5_zero_factor']['records'] if 'reduction' in r]
print('  reduction identity tested on %d records, max relative error %.2e' % (
    len(red), max(r['reduction_rel_error'] for r in red)))
