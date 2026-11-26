"""Exact rational certificate on an interval of s, from the numerical solution of interval_sdp.py.

Unknowns: the entries of the Bernstein coefficients T_{b,k} (k = 0..d).  Equations: the identity
    sum_b (A_b(s) m_b(U))^* T_b(s) (A_b(s) m_b(U)) = F_s(U),      T_b(s) = sum_k B_k(s) T_{b,k},
at d+3 rational nodes s and at rational unitaries U: a linear system over Q.

 1. Modulo one prime, row reduction gives a set of pivot unknowns; the others are free.
 2. Free unknowns := numerical values rounded to a common denominator 2^bits.
 3. The pivot unknowns are then determined.  They are computed modulo several primes and lifted by the Chinese remainder
    theorem and rational reconstruction; two further primes must agree with the lifted values.
 4. Every T_{b,k} is tested for positive definiteness in exact arithmetic; if this fails, bits is increased.
The output is only a candidate.  verify_interval.py checks it independently, in exact arithmetic, with its own nodes.

Usage: python3 interval_exact.py <deg> <a> <b> <d> [number of unitaries]      Output: interval<deg>_<a>_<b>.json"""
import sys, json, time, pickle, random
from fractions import Fraction as Fr
import numpy as np, flint
from exactn import *
PRIMES = [1073741789, 1073741783, 1073741741, 1073741723, 1073741719, 1073741717, 1073741689, 1073741671, 1073741663, 1073741651, 1073741621, 1073741567]
def ratrec(a, m):
    a %= m; bound = int((m // 2) ** 0.5); r0, r1, s0, s1 = m, a, 0, 1
    while r1 > bound: q = r0 // r1; r0, r1, s0, s1 = r1, r0 - q * r1, s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound: return None
    return Fr(r1, s1) if s1 > 0 else Fr(-r1, -s1)
def modq(q, p): return q.numerator % p * pow(q.denominator, -1, p) % p
def matvec_mod(M, x, p):
    """(M @ x) mod p for int64 arrays with entries in [0, p), p < 2^30, without overflow"""
    xl = x & 32767; xh = x >> 15
    return (((M @ xh) % p) * 32768 + (M @ xl)) % p
if __name__ == '__main__':
    deg = int(sys.argv[1]); a, b = Fr(sys.argv[2]), Fr(sys.argv[3]); d = int(sys.argv[4]); t0 = time.time()
    tag = '%d_%s_%s' % (deg, str(a).replace('/', '-'), str(b).replace('/', '-'))
    fam = json.load(open('family%d.json' % deg)); A0, A1 = fam['A0'], fam['A1']; MONS, blocks = space(deg)
    idx = unknown_index(A0); N = len(idx); NT = (d + 1) * N
    nU = int(sys.argv[5]) if len(sys.argv) > 5 else (300 if deg == 3 else 1100)
    rng = random.Random(4001); Us = [rational_unitary(rng) for _ in range(nU)]; assert all(is_unitary(U) for U in Us)
    nodes = [a + (b - a) * Fr(j, d + 2) for j in range(d + 3)]; Fv = [[F_exact(U, sj) for U in Us] for sj in nodes]
    print("degree %d, s in [%s, %s], Bernstein degree %d: unknowns %d ; equations %d (nodes %s)  [%.0fs]" % (deg, a, b, d, NT, nU * len(nodes), ", ".join(map(str, nodes)), time.time() - t0), flush=True)
    def system(p):
        R_, I_ = entries_mod(Us, p); ma, mb = monomials_mod(R_, I_, p, MONS); rows = []
        for sj, fv in zip(nodes, Fv):
            G = gram_columns_mod(ma, mb, A0, A1, sj, blocks, p); tau = (sj - a) / (b - a)
            rows.append(np.hstack([(modq(bernstein(d, k, tau), p) * G) % p for k in range(d + 1)] + [np.array([modq(f, p) for f in fv], dtype=np.int64)[:, None]]))
        return np.vstack(rows)
    # 1. pivot unknowns (cached in pivots<tag>_d<d>.pkl: the row reduction is the slow step at degree 4)
    import os
    cache = 'pivots%s_d%d.pkl' % (tag, d)
    if os.path.exists(cache): pivs = pickle.load(open(cache, 'rb')); rank = len(pivs)
    else:
        p = PRIMES[0]; Mnp = system(p); M = flint.nmod_mat(Mnp.shape[0], NT + 1, Mnp.reshape(-1).tolist(), p); del Mnp
        Rm, rank = M.rref(); del M; pivs = []
        for r in range(rank):
            c = 0 if not pivs else pivs[-1] + 1
            while int(Rm[r, c]) == 0: c += 1
            pivs.append(c)
        del Rm; assert NT not in pivs, "inconsistent system"
        pickle.dump(pivs, open(cache, 'wb'))
    pset = set(pivs); free = np.array([c for c in range(NT) if c not in pset], dtype=np.int64); parr = np.array(pivs, dtype=np.int64)
    print("rank %d ; free unknowns %d  [%.0fs]" % (rank, len(free), time.time() - t0), flush=True)
    if not os.path.exists('interval%s_d%d.pkl' % (tag, d)): sys.exit('pivots cached; no numerical solution yet')
    num = pickle.load(open('interval%s_d%d.pkl' % (tag, d), 'rb')); xnum = np.concatenate(num['x']); assert len(xnum) == NT
    ok = False
    for bits in (32, 44, 56):
        DEN = 2**bits; xfree = [round(float(xnum[c]) * DEN) for c in free]; res = []; used = []
        # 3. pivot unknowns modulo primes
        for p in PRIMES[1:]:
            Mnp = system(p); inv = pow(DEN, -1, p); xf = np.array([v % p * inv % p for v in xfree], dtype=np.int64)
            rhs = (Mnp[:, NT] - matvec_mod(Mnp[:, free], xf, p)) % p; sub = np.hstack([Mnp[:, parr], rhs[:, None]]); del Mnp
            Rs, rk = flint.nmod_mat(sub.shape[0], rank + 1, sub.reshape(-1).tolist(), p).rref(); del sub
            assert rk == rank and all(int(Rs[i, i]) == 1 for i in range(rank)), "pivot set is not independent modulo %d, or the system is inconsistent" % p
            res.append([int(Rs[i, rank]) for i in range(rank)]); used.append(p); del Rs
            if len(used) < 4: continue
            # lift with all but the last two primes, cross-check with the last two
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
        if not good: sys.exit("rational reconstruction did not stabilise with %d primes" % len(used))
        print("  bits %d: pivot unknowns lifted from %d primes, confirmed by 2 more ; largest denominator has %d bits  [%.0fs]" % (bits, len(used) - 2, max(f.denominator for f in xp).bit_length(), time.time() - t0), flush=True)
        x = [None] * NT
        for c, v in zip(free, xfree): x[c] = Fr(v, DEN)
        for c, v in zip(pivs, xp): x[c] = v
        dev = max(abs(float(x[k]) - float(xnum[k])) for k in range(NT))
        T = [[[[None] * len(a0) for _ in a0] for a0 in A0] for _ in range(d + 1)]
        for k in range(d + 1):
            for q, (bi, i, j) in enumerate(idx): T[k][bi][i][j] = x[k * N + q]; T[k][bi][j][i] = x[k * N + q]
        ok = True; minpm = None
        for Tk in T:
            for Tb in Tk:
                for m in range(1, len(Tb) + 1):
                    dt = flint.fmpq_mat(m, m, [flint.fmpq(Tb[i][j].numerator, Tb[i][j].denominator) for i in range(m) for j in range(m)]).det(); dt = Fr(int(dt.p), int(dt.q))
                    ok &= dt > 0; minpm = dt if minpm is None else min(minpm, dt)
        print("  bits %d: max deviation from the numerical point %.2e ; all Bernstein coefficients positive definite: %s (smallest leading minor %.3e)  [%.0fs]" % (bits, dev, ok, float(minpm), time.time() - t0), flush=True)
        if ok: break
    if not ok: sys.exit("no positive definite exact solution found; no certificate written")
    json.dump({'degree': deg, 'a': str(a), 'b': str(b), 'd': d,
               'statement': 'for every s in [a,b] and U in U(4): (1-s)*D[12,12] + s*(1-s)*g0 + s^2*D[13,14] = sum_b (A_b(s) m_b)^* T_b(s) (A_b(s) m_b), A_b(s) = A0_b + s*A1_b, T_b(s) = sum_k Bernstein_k^d((s-a)/(b-a)) T[k][b], every T[k][b] positive definite',
               'A0': A0, 'A1': A1, 'T': [[[[str(q) for q in row] for row in Tb] for Tb in Tk] for Tk in T]}, open('interval%s.json' % tag, 'w'), separators=(',', ':'))
    print("certificate written: interval%s.json" % tag)
