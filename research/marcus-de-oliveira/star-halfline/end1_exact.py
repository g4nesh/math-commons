"""Exact certificate on [a, 1] (every t >= a/(1-a)), from the numerical solution of end1.py.

Unknowns: the entries of T_{b,k} for k < d, and alpha, beta (8 numbers) that fix the last Bernstein coefficient
T_inf = sum_v c_v w_v w_v^T,  c_{4i+j} = alpha_i + beta_j  (see end1.py).  Equations: the identity at d+3 rational nodes
and at rational unitaries, and sum alpha + sum beta = 1.  Solved as in interval_exact.py: pivots modulo one prime, free
unknowns rounded from the numerical solution, pivot unknowns lifted from several primes with two more as cross-check.
Then: every c_v >= 0 and every T_{b,k} (k < d) positive definite, in exact arithmetic.
The output is only a candidate; verify_interval.py checks it independently.

Usage: python3 end1_exact.py <deg> <a> <d>      Output: interval<deg>_<a>_1.json"""
import sys, json, time, pickle, random
from fractions import Fraction as Fr
import numpy as np, flint
from exactn import *
from interval_exact import PRIMES, ratrec, modq, matvec_mod
from end1 import square_coordinates
if __name__ == '__main__':
    deg = int(sys.argv[1]); a = Fr(sys.argv[2]); b = Fr(1); d = int(sys.argv[3]); t0 = time.time(); tag = '%d_%s_1' % (deg, str(a).replace('/', '-'))
    fam = json.load(open('family%d.json' % deg)); A0, A1 = fam['A0'], fam['A1']; MONS, blocks = space(deg)
    idx = unknown_index(A0); N = len(idx); NT = d * N + 8; pos = {key: q for q, key in enumerate(idx)}; sq = square_coordinates(deg)
    nU = 300 if deg == 3 else 1100; rng = random.Random(4001); Us = [rational_unitary(rng) for _ in range(nU)]; assert all(is_unitary(U) for U in Us)
    nodes = [a + (b - a) * Fr(j, d + 2) for j in range(d + 3)]; Fv = [[F_exact(U, sj) for U in Us] for sj in nodes]
    # rank-one Gram vectors of the 16 squares, as sparse rational columns over the unknown index
    gv = []
    for bi, w in sq:
        col = {}
        for i in range(len(w)):
            for j in range(i, len(w)):
                if w[i] * w[j] != 0: col[pos[(bi, i, j)]] = w[i] * w[j]
        gv.append(col)
    def system(p):
        R_, I_ = entries_mod(Us, p); ma, mb = monomials_mod(R_, I_, p, MONS); rows = []
        Gm = np.zeros((N, 16), dtype=np.int64)
        for v, col in enumerate(gv):
            for q, val in col.items(): Gm[q, v] = modq(val, p)
        Lm = np.zeros((16, 8), dtype=np.int64)
        for v in range(16): Lm[v, v // 4] = 1; Lm[v, 4 + v % 4] = 1
        for sj, fv in zip(nodes, Fv):
            G = gram_columns_mod(ma, mb, A0, A1, sj, blocks, p); tau = (sj - a) / (b - a)
            GG = np.zeros((nU, 16), dtype=np.int64)
            for v in range(16): GG[:, v] = matvec_mod(G, Gm[:, v], p)
            last = (modq(bernstein(d, d, tau), p) * ((GG @ Lm) % p)) % p
            rows.append(np.hstack([(modq(bernstein(d, k, tau), p) * G) % p for k in range(d)] + [last, np.array([modq(f, p) for f in fv], dtype=np.int64)[:, None]]))
        rows.append(np.array([[0] * (d * N) + [1] * 8 + [1]], dtype=np.int64))
        return np.vstack(rows)
    p = PRIMES[0]; Mnp = system(p); Rm, rank = flint.nmod_mat(Mnp.shape[0], NT + 1, Mnp.reshape(-1).tolist(), p).rref(); del Mnp; pivs = []
    for r in range(rank):
        c = 0 if not pivs else pivs[-1] + 1
        while int(Rm[r, c]) == 0: c += 1
        pivs.append(c)
    del Rm; assert NT not in pivs, "inconsistent system"
    pset = set(pivs); free = np.array([c for c in range(NT) if c not in pset], dtype=np.int64); parr = np.array(pivs, dtype=np.int64)
    print("degree %d, s in [%s, 1], Bernstein degree %d: unknowns %d ; rank %d ; free unknowns %d, of which %d among alpha, beta  [%.0fs]" % (deg, a, d, NT, rank, len(free), int((free >= d * N).sum()), time.time() - t0), flush=True)
    num = pickle.load(open('interval%s_d%d.pkl' % (tag, d), 'rb')); xnum = np.concatenate(list(num['x']) + [num['alpha'], num['beta']]); assert len(xnum) == NT
    ok = False
    for bits in (32, 44, 56):
        DEN = 2**bits; xfree = [round(float(xnum[c]) * DEN) for c in free]; res = []; used = []
        for p in PRIMES[1:]:
            Mnp = system(p); inv = pow(DEN, -1, p); xf = np.array([v % p * inv % p for v in xfree], dtype=np.int64)
            rhs = (Mnp[:, NT] - matvec_mod(Mnp[:, free], xf, p)) % p; sub = np.hstack([Mnp[:, parr], rhs[:, None]]); del Mnp
            Rs, rk = flint.nmod_mat(sub.shape[0], rank + 1, sub.reshape(-1).tolist(), p).rref(); del sub
            assert rk == rank and all(int(Rs[i, i]) == 1 for i in range(rank))
            res.append([int(Rs[i, rank]) for i in range(rank)]); used.append(p); del Rs
            if len(used) < 4: continue
            lift = used[:-2]; mod = 1
            for q in lift: mod *= q
            xp = []; good = True
            for i in range(rank):
                v = 0; m = 1
                for q, rr in zip(lift, res): v += m * ((rr[i] - v) * pow(m, -1, q) % q); m *= q
                fr = ratrec(v, mod)
                if fr is None or any(fr.denominator % q == 0 or modq(fr, q) != rr[i] for q, rr in zip(used[-2:], res[-2:])): good = False; break
                xp.append(fr)
            if good: break
        if not good: sys.exit("rational reconstruction did not stabilise")
        x = [None] * NT
        for c, v in zip(free, xfree): x[c] = Fr(v, DEN)
        for c, v in zip(pivs, xp): x[c] = v
        dev = max(abs(float(x[k]) - float(xnum[k])) for k in range(NT)); al = x[d * N:d * N + 4]; be = x[d * N + 4:]; cv = [al[v // 4] + be[v % 4] for v in range(16)]
        T = [[[[None] * len(a0) for _ in a0] for a0 in A0] for _ in range(d)]
        for k in range(d):
            for q, (bi, i, j) in enumerate(idx): T[k][bi][i][j] = x[k * N + q]; T[k][bi][j][i] = x[k * N + q]
        ok = all(c >= 0 for c in cv) and sum(al) + sum(be) == 1; minpm = None
        for Tk in T:
            for Tb in Tk:
                for m in range(1, len(Tb) + 1):
                    dt = flint.fmpq_mat(m, m, [flint.fmpq(Tb[i][j].numerator, Tb[i][j].denominator) for i in range(m) for j in range(m)]).det(); dt = Fr(int(dt.p), int(dt.q))
                    ok &= dt > 0; minpm = dt if minpm is None else min(minpm, dt)
        print("  bits %d: lifted from %d primes (+2 as check) ; max deviation from the numerical point %.2e ; all c_v >= 0 and T_k (k < d) positive definite: %s (smallest leading minor %.3e)  [%.0fs]" % (bits, len(used) - 2, dev, ok, float(minpm), time.time() - t0), flush=True)
        if ok: break
    if not ok: sys.exit("no valid exact solution found; no certificate written")
    json.dump({'degree': deg, 'a': str(a), 'b': '1', 'd': d,
               'statement': 'for every s in [a,1] and U in U(4): (1-s)*D[12,12] + s*(1-s)*g0 + s^2*D[13,14] = sum_b (A_b(s) m_b)^* T_b(s) (A_b(s) m_b), A_b(s) = A0_b + s*A1_b, T_b(s) = sum_{k<d} Bernstein_k^d((s-a)/(1-a)) T[k][b] + Bernstein_d^d * T_inf, T[k][b] positive definite, T_inf = sum_v c[v] w_v w_v^T with c[v] >= 0 (w_v supported on block end.block[v])',
               'A0': A0, 'A1': A1, 'T': [[[[str(q) for q in row] for row in Tb] for Tb in Tk] for Tk in T],
               'end': {'alpha': [str(q) for q in al], 'beta': [str(q) for q in be], 'c': [str(q) for q in cv], 'block': [bi for bi, w in sq], 'w': [[str(q) for q in w] for bi, w in sq]}},
              open('interval%s.json' % tag, 'w'), separators=(',', ':'))
    print("alpha = %s ; beta = %s ; certificate written: interval%s.json" % ([str(q) for q in al], [str(q) for q in be], tag))
