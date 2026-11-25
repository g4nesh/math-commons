"""Exact solve for the degree-4 certificate at a rational t.
The identity f_t = sum_b (A_b m_b)^* T_b (A_b m_b), evaluated at rational unitaries, is a linear system over Q in the
entries of the T_b.  Its reduced row echelon form is computed modulo three primes and lifted by rational reconstruction.
Free unknowns are set to rounded values of the numerical strictly feasible point; pivot unknowns follow exactly.
Then every T_b is checked to be positive definite in exact arithmetic.  Output: certificate4_t<t>.json.
(The lifted echelon form is only a means: the certificate is checked independently by verify4.py.)"""
import sys, json, time, pickle, random
from fractions import Fraction as Fr
import numpy as np, flint
from exact4 import *
PRIMES = [1073741789, 1073741783, 1073741741, 1073741723]        # primes below 2^30; the last one is used as a cross-check
def ratrec(a, m):
    """rational reconstruction: n/d = a mod m with |n|, d <= sqrt(m/2); None if impossible"""
    a %= m; bound = int((m // 2) ** 0.5); r0, r1, s0, s1 = m, a, 0, 1
    while r1 > bound: q = r0 // r1; r0, r1, s0, s1 = r1, r0 - q * r1, s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound: return None
    return Fr(r1, s1)
if __name__ == '__main__':
    tq = Fr(sys.argv[1]); t = float(tq); t0 = time.time(); nU = int(sys.argv[2]) if len(sys.argv) > 2 else 1100
    tq2, amp = load_amplitudes('amplitudes4_t%g.json' % t); assert tq2 == tq
    idx = unknown_index(amp); N = len(idx)
    rng = random.Random(4001); Us = [rational_unitary(rng) for _ in range(nU)]; assert all(is_unitary(U) for U in Us)
    fvals = [f_exact(U, tq) for U in Us]
    print("unknowns %d ; rational unitaries %d  [%.0fs]" % (N, nU, time.time() - t0), flush=True)
    rrefs = []; pivs = None
    for p in PRIMES:
        R_, I_ = entries_mod(Us, p); a, b = monomials_mod(R_, I_, p); A = amplitude_columns_mod(a, b, amp, p)
        rhs = np.array([fv.numerator % p * pow(fv.denominator, -1, p) % p for fv in fvals], dtype=np.int64)
        M = flint.nmod_mat(nU, N + 1, np.hstack([A, rhs[:, None]]).reshape(-1).tolist(), p)
        Rm, rank = M.rref()
        tab = np.array([[int(x) for x in row] for row in Rm.table()[:rank]], dtype=np.int64)
        piv = [int(np.nonzero(row)[0][0]) for row in tab]
        assert N not in piv, "inconsistent system modulo %d" % p
        if pivs is None: pivs = piv
        assert piv == pivs, "pivot pattern differs between primes"
        rrefs.append(tab); print("  prime %d: rank %d  [%.0fs]" % (p, rank, time.time() - t0), flush=True)
    rank = len(pivs); free = [c for c in range(N) if c not in set(pivs)]
    # CRT over the first three primes, cross-check with the fourth
    m12 = PRIMES[0] * PRIMES[1]; m123 = m12 * PRIMES[2]; inv1 = pow(PRIMES[0], -1, PRIMES[1]); inv12 = pow(m12, -1, PRIMES[2])
    def crt(r1, r2, r3):
        x = r1 + PRIMES[0] * ((r2 - r1) * inv1 % PRIMES[1]); return x + m12 * ((r3 - x) * inv12 % PRIMES[2])
    rows = []; maxh = 0; bad = 0
    cols = free + [N]
    for r in range(rank):
        ent = {}
        nz = np.nonzero(rrefs[0][r, cols] | rrefs[1][r, cols] | rrefs[2][r, cols])[0]
        for k in nz:
            c = cols[k]; q = ratrec(crt(int(rrefs[0][r, c]), int(rrefs[1][r, c]), int(rrefs[2][r, c])), m123)
            if q is None or (q.numerator % PRIMES[3]) * pow(q.denominator, -1, PRIMES[3]) % PRIMES[3] != int(rrefs[3][r, c]): bad += 1; continue
            ent[c] = q; maxh = max(maxh, abs(q.numerator), q.denominator)
        rows.append(ent)
    print("echelon form lifted: rank %d, free unknowns %d ; entries failing reconstruction or the cross-check prime: %d ; largest numerator/denominator %d  [%.0fs]" % (rank, len(free), bad, maxh, time.time() - t0), flush=True)
    assert bad == 0
    # numerical point -> exact solution
    st = pickle.load(open('strict4_t%g.pkl' % t, 'rb')); xnum = st['x']; assert len(xnum) == N
    x = [None] * N
    DEN = 2**20                                                   # one common denominator keeps the pivot values short
    for c in free: x[c] = Fr(round(float(xnum[c]) * DEN), DEN)
    for r, c in enumerate(pivs):
        ent = rows[r]; val = ent.get(N, Fr(0))
        for fc, q in ent.items():
            if fc != N: val -= q * x[fc]
        x[c] = val
    dev = max(abs(float(x[k]) - float(xnum[k])) for k in range(N))
    print("exact solution: max deviation from the numerical point %.2e  [%.0fs]" % (dev, time.time() - t0), flush=True)
    # blocks and exact positive definiteness
    T = [[[None] * len(a) for _ in a] for a in amp]
    for k, (bi, i, j) in enumerate(idx): T[bi][i][j] = x[k]; T[bi][j][i] = x[k]
    ok = True; minpm = None
    for Tb in T:
        r = len(Tb)
        for k in range(1, r + 1):
            d = flint.fmpq_mat(k, k, [flint.fmpq(Tb[i][j].numerator, Tb[i][j].denominator) for i in range(k) for j in range(k)]).det(); d = Fr(int(d.p), int(d.q))
            ok &= d > 0; minpm = d if minpm is None else min(minpm, d)
    print("all leading principal minors positive: %s ; smallest %.3e  [%.0fs]" % (ok, float(minpm), time.time() - t0))
    json.dump({'t': str(tq), 'degree': 4, 'statement': 't*g0 + (1+t)*D[12,12] + t^2*D[13,14] = sum_b (A_b m_b)^* T_b (A_b m_b) on U(4)',
               'amplitudes': amp, 'T': [[[str(q) for q in row] for row in Tb] for Tb in T]}, open('certificate4_t%g.json' % t, 'w'))
    print("certificate written: certificate4_t%g.json (positive definite: %s)" % (t, ok))
