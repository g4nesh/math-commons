"""Exact rational certificate on an interval of s, from the numerical solution of interval_sdp.py.

Unknowns: the entries of the Bernstein coefficients T_{b,k} (k = 0..d).  Equations: the identity
    sum_b (A_b(s) m_b(U))^* T_b(s) (A_b(s) m_b(U)) = F_s(U),      T_b(s) = sum_k B_k(s) T_{b,k},
at d+3 rational nodes s and at rational unitaries U: a linear system over Q.  Its reduced row echelon form is computed
modulo three primes and lifted by rational reconstruction (a fourth prime cross-checks every entry).  Free unknowns are
set to the rounded numerical values, pivot unknowns follow exactly, and every T_{b,k} is tested for positive definiteness
in exact arithmetic.  The result is only a candidate: verify_interval.py checks it independently, with its own nodes.

Usage: python3 interval_exact.py <deg> <a> <b> <d>      Output: interval<deg>_<a>_<b>.json"""
import sys, json, time, pickle, random
from fractions import Fraction as Fr
import numpy as np, flint
from exactn import *
PRIMES = [1073741789, 1073741783, 1073741741, 1073741723]
def ratrec(a, m):
    a %= m; bound = int((m // 2) ** 0.5); r0, r1, s0, s1 = m, a, 0, 1
    while r1 > bound: q = r0 // r1; r0, r1, s0, s1 = r1, r0 - q * r1, s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound: return None
    return Fr(r1, s1)
if __name__ == '__main__':
    deg = int(sys.argv[1]); a, b = Fr(sys.argv[2]), Fr(sys.argv[3]); d = int(sys.argv[4]); t0 = time.time()
    tag = '%d_%s_%s' % (deg, str(a).replace('/', '-'), str(b).replace('/', '-'))
    fam = json.load(open('family%d.json' % deg)); A0, A1 = fam['A0'], fam['A1']; MONS, blocks = space(deg)
    idx = unknown_index(A0); N = len(idx); NT = (d + 1) * N
    nU = int(sys.argv[5]) if len(sys.argv) > 5 else (300 if deg == 3 else 1100)
    rng = random.Random(4001); Us = [rational_unitary(rng) for _ in range(nU)]; assert all(is_unitary(U) for U in Us)
    nodes = [a + (b - a) * Fr(j, d + 2) for j in range(d + 3)]
    print("degree %d, s in [%s, %s], Bernstein degree %d: unknowns %d ; equations %d (nodes %s)  [%.0fs]" % (deg, a, b, d, NT, nU * len(nodes), ", ".join(map(str, nodes)), time.time() - t0), flush=True)
    rrefs = []; pivs = None
    for p in PRIMES:
        R_, I_ = entries_mod(Us, p); ma, mb = monomials_mod(R_, I_, p, MONS); rows = []
        for sj in nodes:
            G = gram_columns_mod(ma, mb, A0, A1, sj, blocks, p); tau = (sj - a) / (b - a)
            w = [bernstein(d, k, tau) for k in range(d + 1)]; wm = [x.numerator % p * pow(x.denominator, -1, p) % p for x in w]
            rhs = np.array([(lambda f: f.numerator % p * pow(f.denominator, -1, p) % p)(F_exact(U, sj)) for U in Us], dtype=np.int64)
            rows.append(np.hstack([(wk * G) % p for wk in wm] + [rhs[:, None]]))
        Mnp = np.vstack(rows); M = flint.nmod_mat(Mnp.shape[0], NT + 1, Mnp.reshape(-1).tolist(), p); del Mnp, rows
        Rm, rank = M.rref(); piv = []
        for r in range(rank):
            c = 0 if not piv else piv[-1] + 1
            while int(Rm[r, c]) == 0: c += 1
            piv.append(c)
        assert NT not in piv, "inconsistent system modulo %d" % p
        if pivs is None: pivs = piv; free = [c for c in range(NT) if c not in set(pivs)]; cols = free + [NT]
        assert piv == pivs, "pivot pattern differs between primes"
        rrefs.append(np.array([[int(Rm[r, c]) for c in cols] for r in range(rank)], dtype=np.int64)); del Rm, M
        print("  prime %d: rank %d  [%.0fs]" % (p, rank, time.time() - t0), flush=True)
    rank = len(pivs); m12 = PRIMES[0] * PRIMES[1]; m123 = m12 * PRIMES[2]; inv1 = pow(PRIMES[0], -1, PRIMES[1]); inv12 = pow(m12, -1, PRIMES[2])
    def crt(r1, r2, r3):
        x = r1 + PRIMES[0] * ((r2 - r1) * inv1 % PRIMES[1]); return x + m12 * ((r3 - x) * inv12 % PRIMES[2])
    lifted = []; maxh = 0; bad = 0
    for r in range(rank):
        ent = {}
        for k in np.nonzero(rrefs[0][r] | rrefs[1][r] | rrefs[2][r])[0]:
            q = ratrec(crt(int(rrefs[0][r, k]), int(rrefs[1][r, k]), int(rrefs[2][r, k])), m123)
            if q is None or (q.numerator % PRIMES[3]) * pow(q.denominator, -1, PRIMES[3]) % PRIMES[3] != int(rrefs[3][r, k]): bad += 1; continue
            ent[cols[k]] = q; maxh = max(maxh, abs(q.numerator), q.denominator)
        lifted.append(ent)
    print("echelon form lifted: rank %d, free unknowns %d ; entries failing reconstruction or the cross-check prime: %d ; largest numerator/denominator %d  [%.0fs]" % (rank, len(free), bad, maxh, time.time() - t0), flush=True)
    assert bad == 0
    num = pickle.load(open('interval%s_d%d.pkl' % (tag, d), 'rb')); xnum = np.concatenate(num['x']); assert len(xnum) == NT
    for bits in (24, 32, 40, 48):
        x = [None] * NT; DEN = 2**bits
        for c in free: x[c] = Fr(round(float(xnum[c]) * DEN), DEN)
        for r, c in enumerate(pivs):
            ent = lifted[r]; val = ent.get(NT, Fr(0))
            for fc, q in ent.items():
                if fc != NT: val -= q * x[fc]
            x[c] = val
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
        print("free unknowns rounded to denominator 2^%d: max deviation from the numerical point %.2e ; all Bernstein coefficients positive definite: %s (smallest leading minor %.3e)  [%.0fs]" % (bits, dev, ok, float(minpm), time.time() - t0), flush=True)
        if ok: break
    if not ok: sys.exit("no positive definite exact solution found; no certificate written")
    json.dump({'degree': deg, 'a': str(a), 'b': str(b), 'd': d,
               'statement': 'for every s in [a,b] and U in U(4): (1-s)*D[12,12] + s*(1-s)*g0 + s^2*D[13,14] = sum_b (A_b(s) m_b)^* T_b(s) (A_b(s) m_b), A_b(s) = A0_b + s*A1_b, T_b(s) = sum_k Bernstein_k^d((s-a)/(b-a)) T[k][b], every T[k][b] positive definite',
               'A0': A0, 'A1': A1, 'T': [[[[str(q) for q in row] for row in Tb] for Tb in Tk] for Tk in T]}, open('interval%s.json' % tag, 'w'), separators=(',', ':'))
    print("certificate written: interval%s.json" % tag)
