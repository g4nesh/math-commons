"""Second exact checker for an interval certificate, written separately from verify_interval.py.

Origin: written by a separate reviewing agent (same GitHub owner, so this is not an independent review in the sense of
docs/verification.md) that was asked to find a flaw.  It shares no code with exactn.py, verify_interval.py or
../star-exact/exact_space.py, and differs from verify_interval.py in every choice it could:

 (P') every Bernstein coefficient block: symmetric, positive pivots in exact Gaussian elimination (=> positive definite);
      'end' variant: c_v >= 0 and the last coefficient is used as sum_v c_v |w_v . q|^2.
 (I') for each of nU rational unitaries of a DIFFERENT family (Cayley transforms of Gaussian-rational skew-Hermitian
      matrices) the identity is checked as a polynomial identity in s: per-coefficient quadratic forms are computed
      exactly, then both sides are compared at d+4 rational values of s chosen here (one more than needed).
 (E') rank modulo a different prime of (evaluation matrix on W) x (random matrix) for these unitaries.

Usage: python3 recheck_interval.py <certificate.json> <number of unitaries> [seed]     (320 at degree 3, 1100 at degree 4)
"""
import sys, json, itertools, random, time
from fractions import Fraction as Fr
from math import comb, gcd
import numpy as np, flint

def lcm(a, b): return a * b // gcd(a, b)

def cayley(rng):
    """U = (I - K)(I + K)^{-1}, K skew-Hermitian with Gaussian rational entries; returned as (X, Y) lists of 16 Fractions, U = X + iY"""
    den = lambda: rng.choice([1, 2, 3, 4, 5, 7])
    X = [[Fr(0)] * 4 for _ in range(4)]; Y = [[Fr(0)] * 4 for _ in range(4)]       # K = X + iY, X antisymmetric, Y symmetric
    for i in range(4):
        Y[i][i] = Fr(rng.randint(-6, 6), den())
        for j in range(i + 1, 4):
            x = Fr(rng.randint(-6, 6), den()); y = Fr(rng.randint(-6, 6), den())
            X[i][j] = x; X[j][i] = -x; Y[i][j] = y; Y[j][i] = y
    def real8(X, Y, sign):   # real 8x8 representation of I + sign*K
        M = flint.fmpq_mat(8, 8)
        for i in range(4):
            for j in range(4):
                e = 1 if i == j else 0
                x = e + sign * X[i][j]; y = sign * Y[i][j]
                M[i, j] = flint.fmpq(x.numerator, x.denominator); M[i + 4, j + 4] = flint.fmpq(x.numerator, x.denominator)
                M[i + 4, j] = flint.fmpq(y.numerator, y.denominator); M[i, j + 4] = flint.fmpq(-y.numerator, y.denominator)
        return M
    Um = real8(X, Y, -1) * real8(X, Y, 1).inv()
    f = lambda q: Fr(int(q.p), int(q.q))
    UX = [f(Um[i, j]) for i in range(4) for j in range(4)]; UY = [f(Um[i + 4, j]) for i in range(4) for j in range(4)]
    # exact unitarity
    for i in range(4):
        for j in range(4):
            re = sum(UX[4 * i + k] * UX[4 * j + k] + UY[4 * i + k] * UY[4 * j + k] for k in range(4))
            im = sum(UY[4 * i + k] * UX[4 * j + k] - UX[4 * i + k] * UY[4 * j + k] for k in range(4))
            assert re == (1 if i == j else 0) and im == 0
    return UX, UY

def cmul(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def abs2(a): return a[0] * a[0] + a[1] * a[1]

def F_coeffs(Z, L, deg):
    """integer coefficients (f0, f1, f2) of L^(2 deg) * F_s = f0 + f1 s + f2 s^2, Z[v] = L * u_v as Gaussian integers, v = 4 row + col (0-based)"""
    u = lambda r, c: Z[4 * (r - 1) + (c - 1)]                                                   # 1-based entry
    def D(S, T):
        a = cmul(u(S[0], T[0]), u(S[1], T[1])); b = cmul(u(S[0], T[1]), u(S[1], T[0]))
        return abs2((a[0] - b[0], a[1] - b[1]))                                                 # L^4 |det|^2
    Bv = lambda r, c: abs2(u(r, c)) * L * L                                                     # L^4 |u|^2
    g0 = Bv(1, 1) + Bv(3, 3) + D((1, 4), (3, 4)) - D((1, 2), (1, 4))
    D1 = D((1, 2), (1, 2)); D2 = D((1, 3), (1, 4)); sc = L ** (2 * deg - 4)
    # (1-s) D1 + s(1-s) g0 + s^2 D2 = D1 + s (g0 - D1) + s^2 (D2 - g0)
    return sc * D1, sc * (g0 - D1), sc * (D2 - g0)

if __name__ == '__main__':
    path = sys.argv[1]; nU = int(sys.argv[2]); seed = int(sys.argv[3]) if len(sys.argv) > 3 else 20261125; t0 = time.time()
    C = json.load(open(path)); deg = C['degree']; a = Fr(C['a']); b = Fr(C['b']); d = C['d']; end = C.get('end')
    MONS = list(itertools.combinations_with_replacement(range(16), deg)); blocks = {}
    for i, m in enumerate(MONS):
        blocks.setdefault((tuple(sorted(v // 4 for v in m)), tuple(sorted(v % 4 for v in m))), []).append(i)
    blocks = list(blocks.values()); nb = len(blocks)
    A0, A1 = C['A0'], C['A1']; assert len(A0) == len(A1) == nb
    for a0, a1, bl in zip(A0, A1, blocks):
        assert len(a0) == len(a1)
        for r in a0 + a1: assert len(r) == len(bl) and all(type(c) is int for c in r)
    nT = len(C['T']); assert nT == (d if end else d + 1) and a < b and (not end or b == 1)
    # ---- (P') ----
    okP = True; Tint = []      # Tint[k][bi] = (den, integer upper-triangular list of (i, j, value with multiplicity))
    for Tk in C['T']:
        assert len(Tk) == nb; row_k = []
        for a0, Tb in zip(A0, Tk):
            r = len(a0); assert len(Tb) == r and all(len(row) == r for row in Tb)
            M = [[Fr(q) for q in row] for row in Tb]
            okP &= all(M[i][j] == M[j][i] for i in range(r) for j in range(r))
            W = [row[:] for row in M]                                   # exact elimination, pivots must be > 0
            for i in range(r):
                if not W[i][i] > 0: okP = False; break
                for j in range(i + 1, r):
                    f = W[j][i] / W[i][i]
                    if f:
                        for k in range(i, r): W[j][k] -= f * W[i][k]
            den = 1
            for row in M:
                for q in row: den = lcm(den, q.denominator)
            row_k.append((den, [(i, j, int(M[i][j] * den) * (1 if i == j else 2)) for i in range(r) for j in range(i, r) if M[i][j] != 0]))
        Tint.append(row_k)
    if end:
        cv = [Fr(q) for q in end['c']]; okP &= all(c >= 0 for c in cv); assert len(cv) == len(end['block']) == len(end['w'])
        wv = [[Fr(q) for q in w] for w in end['w']]
        for bi, w in zip(end['block'], wv): assert len(w) == len(A0[bi])
    print("(P') %s: degree %d, [a,b] = [%s, %s], d = %d, %d amplitudes, %d blocks, end variant: %s ; all coefficients symmetric with positive pivots%s: %s  [%.0fs]"
          % (path.split('/')[-1], deg, a, b, d, sum(len(x) for x in A0), nb, bool(end), " and c_v >= 0" if end else "", okP, time.time() - t0), flush=True)
    # sparse integer amplitude rows
    sp = lambda rows, bl: [[(bl[k], c) for k, c in enumerate(r) if c] for r in rows]
    S0 = [sp(a0, bl) for a0, bl in zip(A0, blocks)]; S1 = [sp(a1, bl) for a1, bl in zip(A1, blocks)]
    sblocks = [bi for bi in range(nb) if any(S1[bi])]                   # blocks whose amplitudes depend on s
    # ---- (I') ----
    rng = random.Random(seed); svals = [Fr(-3, 2), Fr(-1, 3), Fr(0), Fr(2, 7), Fr(1, 2), Fr(5, 7), Fr(1), Fr(9, 5), Fr(3), Fr(4)][:d + 4]
    bern = {s: [comb(d, k) * ((s - a) / (b - a)) ** k * (1 - (s - a) / (b - a)) ** (d - k) for k in range(d + 1)] for s in svals}
    bad = 0; Zs = []; maxL = 0
    def qform(den_tri, q):
        den, tri = den_tri; acc = 0
        for i, j, v in tri: acc += v * (q[i][0] * q[j][0] + q[i][1] * q[j][1])
        return Fr(acc, den)
    for n in range(nU):
        UX, UY = cayley(rng); L = 1
        for q in UX + UY: L = lcm(L, q.denominator)
        maxL = max(maxL, L); Z = [(int(x * L), int(y * L)) for x, y in zip(UX, UY)]; Zs.append((UX, UY))
        mon = []
        for m in MONS:
            z = Z[m[0]]
            for v in m[1:]: z = cmul(z, Z[v])
            mon.append(z)
        amp = lambda S: [(sum(c * mon[k][0] for k, c in row), sum(c * mon[k][1] for k, c in row)) for row in S]
        c0 = [Fr(0)] * (d + 1); c1 = [Fr(0)] * (d + 1); c2 = [Fr(0)] * (d + 1)
        for bi in range(nb):
            if not S0[bi]: continue
            q0 = amp(S0[bi]); dep = bi in sblocks; q1 = amp(S1[bi]) if dep else None
            for k in range(nT):
                c0[k] += qform(Tint[k][bi], q0)
                if dep:
                    # (q0 + s q1)^* T (q0 + s q1) = q0^*Tq0 + s * 2Re(q0^* T q1) + s^2 q1^*Tq1
                    den, tri = Tint[k][bi]; acc = 0
                    for i, j, v in tri:                                  # v carries the factor 2 for i<j
                        if i == j: acc += 2 * v * (q0[i][0] * q1[i][0] + q0[i][1] * q1[i][1])
                        else: acc += v * (q0[i][0] * q1[j][0] + q0[i][1] * q1[j][1] + q0[j][0] * q1[i][0] + q0[j][1] * q1[i][1])
                    c1[k] += Fr(acc, den); c2[k] += qform(Tint[k][bi], q1)
        if end:                                                          # last coefficient: sum_v c_v |w_v . q|^2
            for c, bi, w in zip(cv, end['block'], wv):
                q0 = amp(S0[bi]); q1 = amp(S1[bi])
                p0 = (sum(x * q[0] for x, q in zip(w, q0)), sum(x * q[1] for x, q in zip(w, q0)))
                p1 = (sum(x * q[0] for x, q in zip(w, q1)), sum(x * q[1] for x, q in zip(w, q1)))
                c0[d] += c * abs2(p0); c1[d] += 2 * c * (p0[0] * p1[0] + p0[1] * p1[1]); c2[d] += c * abs2(p1)
        f0, f1, f2 = F_coeffs(Z, L, deg)
        for s in svals:
            lhs = f0 + f1 * s + f2 * s * s
            rhs = sum(bern[s][k] * (c0[k] + c1[k] * s + c2[k] * s * s) for k in range(d + 1))
            bad += (lhs != rhs)
        if n % 100 == 99: print("   %d unitaries done, failures so far %d  [%.0fs]" % (n + 1, bad, time.time() - t0), flush=True)
    print("(I') identity as a polynomial in s (checked at %d values of s: %s) fails at %d of %d (unitary, s) pairs ; largest common denominator of a unitary has %d bits  [%.0fs]"
          % (len(svals), ", ".join(map(str, svals)), bad, nU * len(svals), maxL.bit_length(), time.time() - t0), flush=True)
    # ---- (E') ----
    p = 65521; R = np.zeros((nU, 16), dtype=np.int64); I = np.zeros((nU, 16), dtype=np.int64)
    for n, (UX, UY) in enumerate(Zs):
        for v in range(16):
            R[n, v] = UX[v].numerator % p * pow(UX[v].denominator, -1, p) % p; I[n, v] = UY[v].numerator % p * pow(UY[v].denominator, -1, p) % p
    MA = np.array(MONS); ma = R[:, MA[:, 0]].copy(); mb = I[:, MA[:, 0]].copy()
    for k in range(1, deg):
        cR = R[:, MA[:, k]]; cI = I[:, MA[:, k]]; ma, mb = (ma * cR - mb * cI) % p, (ma * cI + mb * cR) % p
    Pi = np.array([bl[i] for bl in blocks for i in range(len(bl)) for j in range(i, len(bl))]); Qi = np.array([bl[j] for bl in blocks for i in range(len(bl)) for j in range(i, len(bl))])
    dimW = len(Pi); ncol = min(dimW, nU); g = np.random.default_rng(seed); Om = g.integers(0, p, size=(dimW, ncol), dtype=np.int64)
    EO = np.zeros((nU, ncol), dtype=np.int64)
    for lo in range(0, dimW, 2000):
        hi = min(lo + 2000, dimW); E = (ma[:, Pi[lo:hi]] * ma[:, Qi[lo:hi]] + mb[:, Pi[lo:hi]] * mb[:, Qi[lo:hi]]) % p
        EO = (EO + E @ Om[lo:hi]) % p                                     # entries < 2000 * p^2 < 2^63
    rE = flint.nmod_mat(nU, ncol, EO.reshape(-1).tolist(), p).rank()
    print("(E') dim W = %d ; rank modulo %d of (evaluations at these %d unitaries) x (random %d x %d matrix): %d  [%.0fs]" % (dimW, p, nU, dimW, ncol, rE, time.time() - t0), flush=True)
    print("RESULT %s: P' %s, I' %s, rank %d" % (path.split('/')[-1], okP, bad == 0, rE))
