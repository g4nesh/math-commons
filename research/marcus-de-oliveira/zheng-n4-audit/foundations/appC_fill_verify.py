"""Independent verification of appendix-c-open-items.md (label appC-verify).

Own code only.  Reads nothing from the supplement and does not import or rerun
appC_fill_checks.py.  From ../lib/moc4lib.py it uses only the transcribed tables:
MATCHINGS (A.1), cubic_C (A.2), n_vec (eq. (4.1)), and small Fraction helpers
(inv, nullspace, matmul, det_leibniz4, QI).

rho~ is built here by a different route from the audited script: as the unique
linear functional with <rho~(V), m(z)> = det((z_i - z_{4+j}) V_ij), solved at 14 points
and certified on the 256-point grid {0,1}^8 (both sides are multi-affine in z), then
compared with the complementary-minor definition on all 70 pairs.

Sections (each prints PASS / FAIL / INFO lines; exit code 1 on any FAIL):
  A  Proposition C.4 and Corollary C.5 (audit Section 1)
  B  Lemma C.11 without (C.2) (audit Section 4)
  C  Lemma C.10 weights without Lam (audit Section 5.1), plus a 50-digit end-to-end test
  D  (C.13), (C.14) and the density step (audit Section 5.2)
  E  (H3) and (C.8) at special stable points; numeric scan for real failures of (H3)
Run:  python -I -B appC_fill_verify.py [A B C D E]     (deterministic seeds)
"""
import itertools
import math
import os
import random
import sys
import time
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L  # noqa: E402
import flint  # noqa: E402
import sympy as sp  # noqa: E402
import numpy as np  # noqa: E402
import mpmath  # noqa: E402

QI = L.QI
FAILS = []


def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + (("  [" + str(info) + "]") if info != "" else ""))
    sys.stdout.flush()
    if not ok:
        FAILS.append(name)


def info(name, msg):
    print("INFO " + name + "  [" + str(msg) + "]")
    sys.stdout.flush()


# ----------------------------------------------------------------------------- helpers
def fq(x):
    x = Fr(x)
    return flint.fmpq(x.numerator, x.denominator)


def frq(x):
    return Fr(int(x.p), int(x.q))


def qrank(rows):
    if not rows:
        return 0
    return flint.fmpq_mat([[fq(x) for x in r] for r in rows]).rank()


def re_(x):
    return x.re if isinstance(x, QI) else Fr(x)


def im_(x):
    return x.im if isinstance(x, QI) else Fr(0)


def crank(rows):
    """rank over C of a matrix with Gaussian-rational entries (real 2x2 block embedding)."""
    R = []
    for r in rows:
        R.append([re_(x) for x in r] + [-im_(x) for x in r])
        R.append([im_(x) for x in r] + [re_(x) for x in r])
    return qrank(R) // 2


def psign(p):
    s = 1
    p = list(p)
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


def gdet(M):
    n = len(M)
    if n == 0:
        return 1
    tot = 0
    for p in itertools.permutations(range(n)):
        t = 1
        for i in range(n):
            t = t * M[i][p[i]]
        tot = tot + t if psign(p) > 0 else tot - t
    return tot


def sub(M, I, J):  # 1-based index tuples
    return [[M[i - 1][j - 1] for j in J] for i in I]


def comp(I):
    return tuple(k for k in range(1, 5) if k not in I)


ALL_PAIRS = [(I, J) for k in range(5) for I in itertools.combinations(range(1, 5), k)
             for J in itertools.combinations(range(1, 5), k)]


def mvec(z):
    """matching vector m(z): B_r(z) = prod over pairs (i,j) of (z_i - z_j)."""
    out = []
    for M in L.MATCHINGS:
        p = 1
        for (i, j) in M:
            p = p * (z[i - 1] - z[j - 1])
        out.append(p)
    return out


def dmvec(z, dz):
    """directional derivative of m at z along dz (product rule)."""
    out = []
    for M in L.MATCHINGS:
        tot = 0
        for k in range(4):
            t = 1
            for idx, (i, j) in enumerate(M):
                t = t * ((dz[i - 1] - dz[j - 1]) if idx == k else (z[i - 1] - z[j - 1]))
            tot = tot + t
        out.append(tot)
    return out


def in_span(v, basis):
    return qrank(basis) == qrank(basis + [v])


# Hessian tensor of C:  T[k][i][j] = d^3 C / du_k du_i du_j  (constants)
MS = sp.symbols("m0:14")
CS = sp.expand(L.cubic_C(MS))
TENS = {}
for k in range(14):
    for i in range(14):
        for j in range(14):
            v = sp.diff(CS, MS[k], MS[i], MS[j])
            if v != 0:
                TENS[(k, i, j)] = int(v)


def Hm(u):
    H = [[Fr(0)] * 14 for _ in range(14)]
    for (k, i, j), v in TENS.items():
        if u[k] != 0:
            H[i][j] += v * u[k]
    return H


def gradC(u):
    g = [Fr(0)] * 14
    for (k, i, j), v in TENS.items():
        g[k] += Fr(v) * u[i] * u[j] / 2
    return g


def qfun(u):
    return u[0] * u[6] - u[1] * u[5]


# ============================================================================ Section A
def section_A():
    print("\n=== A. Proposition C.4 / Corollary C.5 (audit Section 1) ===")
    t0 = time.time()
    ctx = flint.fmpq_mpoly_ctx.get(("v", 16), "lex")
    g = ctx.gens()
    ZERO = 0 * g[0]
    VS = [[g[4 * i + j] for j in range(4)] for i in range(4)]

    def detz(z):
        return ZERO + L.det_leibniz4([[fq(z[i] - z[4 + j]) * VS[i][j] for j in range(4)]
                                      for i in range(4)])
    rng = random.Random(20261125)
    rows, pts = [], []
    while len(rows) < 14:
        z = [Fr(rng.randint(-7, 7)) for _ in range(8)]
        m = mvec(z)
        if qrank(rows + [m]) > len(rows):
            rows.append(m); pts.append(z)
    Minv = L.inv(rows)
    dets = [detz(z) for z in pts]
    rho = [ZERO + sum((fq(Minv[j][k]) * dets[k] for k in range(14)), ZERO) for j in range(14)]
    ok = True
    for bits in itertools.product((0, 1), repeat=8):
        z = [Fr(b) for b in bits]
        m = mvec(z)
        lhs = sum((fq(m[j]) * rho[j] for j in range(14)), ZERO)
        if lhs - detz(z) != ZERO:
            ok = False
            break
    check("A1 rho~ defined by <rho~(V), m(z)> = det((z_i - z_{4+j}) V_ij) at 14 points satisfies it on "
          "all of {0,1}^8, hence identically in (z, V) (both sides multi-affine in z)", ok)
    ok = True
    for (I, J) in ALL_PAIRS:
        nv = L.n_vec(I, J)
        lhs = sum((fq(nv[j]) * rho[j] for j in range(14)), ZERO)
        tgt = ZERO + (-1) ** (sum(I) + sum(J)) * gdet(sub(VS, I, J)) * gdet(sub(VS, comp(I), comp(J)))
        if lhs - tgt != ZERO:
            ok = False
    check("A2 the same rho~ satisfies <rho~(V), n_IJ> = (-1)^(sum I + sum J) det V[I,J] det V[I^c,J^c] "
          "for all 70 (I,J) (polynomial identity; so rho~ = det(V) rho of (C.3))", ok)
    info("A2b", "torus identity rho~(D1 V D2) = det D1 det D2 rho~(V) is immediate from A1: "
         "D_a D1 V D2 + D1 V D2 D_b = D1 (D_a V + V D_b) D2")

    dr = [[rho[k].derivative(a) for a in range(16)] for k in range(14)]

    def ev(p, v):
        return frq(p(*[fq(x) for x in v]))

    def jac(v):
        return [[ev(dr[k][a], v) for a in range(16)] for k in range(14)]

    def omega_conditions(z, w, tag):
        A = [[Fr(1)] * 8, list(z), list(w), [z[i] * w[i] for i in range(8)]]
        A0 = [r[:4] for r in A]; A1 = [r[4:] for r in A]
        if L.det(A0) == 0:
            check(tag + " (i) det A0 != 0", False)
            return None
        A0i = L.inv(A0)
        V = L.matmul(A0i, A1)
        v = [V[i][j] for i in range(4) for j in range(4)]
        J = jac(v)
        rJ = qrank(J)
        Mz = [[(z[i] - z[4 + j]) * V[i][j] for j in range(4)] for i in range(4)]
        Mw = [[(w[i] - w[4 + j]) * V[i][j] for j in range(4)] for i in range(4)]
        rz, rw = qrank(Mz), qrank(Mw)
        N = L.nullspace(L.transpose(J))
        rho0 = [ev(rho[k], v) for k in range(14)]

        def forms(M):
            Rk = L.nullspace(M); Lk = L.nullspace(L.transpose(M))
            out = []
            for l in Lk:
                for r in Rk:
                    row = [Fr(0)] * 8
                    for i in range(4):
                        for j in range(4):
                            c = l[i] * V[i][j] * r[j]
                            row[i] += c; row[4 + j] -= c
                    out.append(row)
            return out
        fz, fw = forms(Mz), forms(Mw)
        Tz, Tw = L.nullspace(fz), L.nullspace(fw)
        img = [mvec(z), mvec(w)] + [dmvec(z, t) for t in Tz] + [dmvec(w, t) for t in Tw]
        rimg = qrank(img)
        lemma12 = True
        for mm in (mvec(z), mvec(w)):
            if any(sum(mm[k] * J[k][a] for k in range(14)) != 0 for a in range(16)):
                lemma12 = False
            if sum(mm[k] * rho0[k] for k in range(14)) != 0:
                lemma12 = False
        allin = all(in_span(x, N) for x in img)
        euler = all(sum(n[k] * rho0[k] for k in range(14)) == 0 for n in N)
        ok = (rJ == 10 and len(N) == 4 and rz == 2 and rw == 2 and qrank(fz) == 4 and qrank(fw) == 4
              and rimg == 4 and lemma12 and allin and euler)
        check(tag + " Omega conditions (i)-(v) hold: det A0 != 0; rank d rho~ = 10; rank M_z = rank M_w = 2; "
              "transversality 4+4; join rank 4; m(z), m(w) annihilate rho~(V) and im d rho~ (Lemma 1.2); "
              "join image inside N(V); Euler", ok,
              "rank J %d, dim N %d, ranks M %d %d, forms %d %d, join %d, all V_ij != 0: %s, det A1 = %s"
              % (rJ, len(N), rz, rw, qrank(fz), qrank(fw), rimg, all(x != 0 for x in v), L.det(A1)))
        return V, v, J, N, rho0, A0i

    z0 = [Fr(x) for x in (0, 1, 3, -2, 5, -4, 7, 2)]
    w0 = [Fr(x) for x in (2, -1, 4, 3, -3, 6, 1, -5)]
    V0, v0, J0, N0, rho00, A0i = omega_conditions(z0, w0, "A3 audit test point:")
    # dominance of (z, w, lambda) -> lambda rho~(V(z, w))
    cols = []
    for p in range(16):
        dA = [[Fr(0)] * 8 for _ in range(4)]
        i = p % 8
        if p < 8:
            dA[1][i] = Fr(1); dA[3][i] = w0[i]
        else:
            dA[2][i] = Fr(1); dA[3][i] = z0[i]
        dA0 = [r[:4] for r in dA]; dA1 = [r[4:] for r in dA]
        rhs = [[dA1[r][c] - sum(dA0[r][k] * V0[k][c] for k in range(4)) for c in range(4)] for r in range(4)]
        dV = L.matmul(A0i, rhs)
        dv = [dV[a // 4][a % 4] for a in range(16)]
        cols.append([sum(J0[k][a] * dv[a] for a in range(16)) for k in range(14)])
    rd = qrank(cols + [rho00])
    check("A4 differential of (z,w,lambda) -> lambda rho~(V(z,w)) has rank 10 at the test point (dense image "
          "in the 10-dim cone M^)", rd == 10, "rank %d" % rd)
    rng2 = random.Random(7781)
    z1 = [Fr(rng2.randint(-9, 9), rng2.randint(1, 3)) for _ in range(8)]
    w1 = [Fr(rng2.randint(-9, 9), rng2.randint(1, 3)) for _ in range(8)]
    omega_conditions(z1, w1, "A5 second random rational point (own seed):")
    # Lemma 1.3 at three random points
    ok = True
    for _ in range(3):
        z = [Fr(rng2.randint(-9, 9)) for _ in range(8)]
        w = [Fr(rng2.randint(-9, 9)) for _ in range(8)]
        A = [[Fr(1)] * 8, z, w, [z[i] * w[i] for i in range(8)]]
        A0 = [r[:4] for r in A]; A1 = [r[4:] for r in A]
        if L.det(A0) == 0:
            continue
        V = L.matmul(L.inv(A0), A1)
        for zz in (z, w):
            if qrank([[(zz[i] - zz[4 + j]) * V[i][j] for j in range(4)] for i in range(4)]) > 2:
                ok = False
    check("A6 Lemma 1.3 (Segre configurations give rank <= 2) at 3 random rational points", ok)
    # rank d rho~ <= 10 also at V with zero entries (closedness argument)
    Vz = [[Fr(1), Fr(0), Fr(2), Fr(-1)], [Fr(0), Fr(3), Fr(1), Fr(1)], [Fr(2), Fr(1), Fr(0), Fr(5)],
          [Fr(1), Fr(-2), Fr(1), Fr(0)]]
    rz = qrank(jac([Vz[i][j] for i in range(4) for j in range(4)]))
    check("A7 rank d rho~ <= 10 at a V with four zero entries (Lemma 1.1 'everywhere')", rz <= 10, "rank %d" % rz)
    # consequences of N(V0) in Sigma2^: C vanishes on N(V0); rank H <= 12 there
    cs = sp.symbols("c0:4")
    u = [sum(sp.Rational(N0[i][k].numerator, N0[i][k].denominator) * cs[i] for i in range(4)) for k in range(14)]
    Cn = sp.expand(L.cubic_C(u))
    rks = []
    rr = random.Random(5)
    for _ in range(3):
        cc = [Fr(rr.randint(-5, 5)) for _ in range(4)]
        uu = [sum(cc[i] * N0[i][k] for i in range(4)) for k in range(14)]
        rks.append(qrank(Hm(uu)))
    check("A8 consequence test of N(V0) in Sigma2^: C restricted to N(V0) is identically 0 and rank H = 12 "
          "at 3 random points of N(V0)", Cn == 0 and all(r <= 12 for r in rks), "ranks %s" % rks)
    # dim Sigma2^ = 12 directly: rank of the join differential (Terracini) at one rational point
    z = [Fr(x) for x in (0, 2, -3, 5, 1, -4, 7, 3)]
    w = [Fr(x) for x in (1, -2, 6, 4, -5, 3, 2, -7)]
    cols = [mvec(z), mvec(w)]
    for k in range(8):
        e = [Fr(int(i == k)) for i in range(8)]
        cols.append(dmvec(z, e)); cols.append(dmvec(w, e))
    rj = qrank(cols)
    check("A9 dim Sigma2^ >= 12 directly: (lam,z,mu,w) -> lam m(z) + mu m(w) has differential rank 12 at a "
          "rational point (so dim Sigma2 = 11, independent of the dual-variety argument)", rj == 12, "rank %d" % rj)
    # dim (M8')^-hat >= 12: tangent of E and projection rank, at V0 and at the second point
    d2 = [[[rho[k].derivative(a).derivative(b) for b in range(16)] for a in range(16)] for k in range(14)]
    for tag, (V, v, J, N) in (("V0", (V0, v0, J0, N0)),):
        rr = random.Random(99)
        cf = [Fr(rr.randint(-6, 6)) for _ in range(4)]
        xi = [sum(cf[i] * N[i][k] for i in range(4)) for k in range(14)]
        Hk = [[[ev(d2[k][a][b], v) for b in range(16)] for a in range(16)] for k in range(14)]
        eqs = []
        for a in range(16):
            eqs.append([sum(xi[k] * Hk[k][a][b] for k in range(14)) for b in range(16)] +
                       [J[k][a] for k in range(14)])
        TE = L.nullspace(eqs)
        rp = qrank([t[16:] for t in TE])
        check("A10 at %s: tangent space of E = {(V, xi): xi annihilates im d rho~_V} has dim 20 and its projection "
              "to xi has rank 12, so dim (M8')^-hat >= 12" % tag, len(TE) == 20 and rp == 12,
              "dim %d, rank %d" % (len(TE), rp))
    print("   section A time %.1f s" % (time.time() - t0))


# ============================================================================ Section B
W4 = [0, 1, 5, 6]
KP = [k for k in range(14) if k not in W4 and k != 13]
KAP = [Fr(0)] * 13 + [Fr(1)]


def detpoly(u):
    """coefficients p0..p14 of det(H(u) + h H(kappa)) by exact interpolation."""
    Hu, HK = Hm(u), Hm(KAP)
    vals = []
    for h in range(15):
        vals.append(flint.fmpq_mat([[fq(Hu[i][j] + h * HK[i][j]) for j in range(14)] for i in range(14)]).det())
    Vd = flint.fmpq_mat([[flint.fmpq(h ** e) for e in range(15)] for h in range(15)])
    sol = Vd.solve(flint.fmpq_mat([[v] for v in vals]))
    return [frq(sol[e, 0]) for e in range(15)]


def smooth_rank(y):
    """rank of {dC_y} together with the differentials at y of all 13x13 minors of H of rank 12."""
    Hy = Hm(y)
    funcs = [gradC(y)]
    for r in range(14):
        for c in range(14):
            I = [i for i in range(14) if i != r]; J = [j for j in range(14) if j != c]
            blk = [[Hy[i][j] for j in J] for i in I]
            if qrank(blk) != 12:
                continue
            rv = L.nullspace(blk)[0]; lv = L.nullspace(L.transpose(blk))[0]
            Ipos = {i: a for a, i in enumerate(I)}; Jpos = {j: b for b, j in enumerate(J)}
            f = [Fr(0)] * 14
            for (k, i, j), val in TENS.items():
                if i in Ipos and j in Jpos:
                    f[k] += val * lv[Ipos[i]] * rv[Jpos[j]]
            funcs.append(f)
    return qrank(funcs), len(funcs) - 1


def yvec(a, b, theta):
    z = [QI(x) for x in a] + [QI(0) - x for x in b]
    return [re_(theta * x) for x in mvec(z)]


def section_B():
    print("\n=== B. Lemma C.11 without (C.2) (audit Section 4) ===")
    t0 = time.time()
    q = MS[0] * MS[6] - MS[1] * MS[5]
    check("B1 dC/du13 = q and d2C/du13^2 = 0", sp.expand(sp.diff(CS, MS[13]) - q) == 0
          and sp.diff(CS, MS[13], 2) == 0)
    ok = all(not (i not in W4 and j == 13) for (k, i, j) in TENS) and \
        all((i in W4 and j in W4) for (k, i, j) in TENS if k == 13)
    check("B2 H(u)_{i,13} = 0 for i not in W={0,1,5,6} (all u) and H(kappa) supported on W x W", ok)
    rng = random.Random(3141)
    ok, degok = True, True
    for _ in range(6):
        u = [Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(14)]
        co = detpoly(u)
        Hu = Hm(u)
        dK = frq(flint.fmpq_mat([[fq(Hu[i][j]) for j in KP] for i in KP]).det())
        if any(c != 0 for c in co[4:]):
            degok = False
        if co[3] != -2 * qfun(u) * dK:
            ok = False
    check("B3 Lemma 4.1 at 6 random rational u (own seed): degree <= 3 in h and p3 = -2 q(u) det H(u)[K',K']",
          ok and degok)
    a = [Fr(0), Fr(1), Fr(3), Fr(7)]
    b = [QI(2, 1), QI(-1, 3), QI(5, -2), QI(4, 4)]
    y = yvec(a, b, QI(1))
    yp = [-36, 108, 189, -569, 344, -126, 378, 288, 168, -784, 496, -512, 512, 16]
    co = detpoly(y)
    rH = qrank(Hm(y))
    check("B4 paper example: y as printed; q(y) = C(y) = 0; rank H(y) = 12; det(H(y)+hH(kappa)) = p2 h^2 with "
          "p2 = -114105877589169227887017984",
          y == [Fr(x) for x in yp] and qfun(y) == 0 and L.cubic_C(y) == 0 and rH == 12
          and co[0] == co[1] == co[3] == 0 and all(c == 0 for c in co[4:])
          and co[2] == -114105877589169227887017984, "rank %d p2 %s" % (rH, co[2]))
    sr, nm = smooth_rank(y)
    check("B5 smoothness criterion at the example: dC_y and the differentials of ALL rank-12 13x13 minors of H span "
          "rank >= 2 (so dim T_y Z <= 12)", sr >= 2, "rank %d from %d minors" % (sr, nm))
    rng = random.Random(2718)
    for t in range(3):
        a = sorted(rng.sample(range(-9, 10), 4))
        a = [Fr(x) for x in a]
        b = [QI(rng.randint(-6, 6), rng.choice([-5, -3, -2, -1, 1, 2, 3, 4])) for _ in range(4)]
        th = QI(rng.randint(-3, 3), rng.choice([-2, -1, 1, 2]))
        y = yvec(a, b, th)
        co = detpoly(y)
        rH = qrank(Hm(y))
        sr, nm = smooth_rank(y)
        ok = (qfun(y) == 0 and L.cubic_C(y) == 0 and rH <= 12 and co[0] == co[1] == co[3] == 0
              and all(c == 0 for c in co[4:]) and co[2] != 0 and sr >= 2)
        check("B6.%d random parameters (a real, b, theta Gaussian rational): q = C = 0, rank H(y) <= 12, "
              "det = p2 h^2 with p2 != 0, smoothness rank >= 2" % t, ok,
              "rank H %d, p2 != 0: %s, smooth rank %d" % (rH, co[2] != 0, sr))
    # special (non-generic) parameters: report only
    specials = [("a1 = a2", [0, 0, 3, 7], [QI(2, 1), QI(-1, 3), QI(5, -2), QI(4, 4)], QI(1)),
                ("b real", [0, 1, 3, 7], [QI(2), QI(-1), QI(5), QI(4)], QI(1, 1)),
                ("a1 + b1 = 0", [0, 1, 3, 7], [QI(0), QI(-1, 3), QI(5, -2), QI(4, 4)], QI(1)),
                ("theta = i", [0, 1, 3, 7], [QI(2, 1), QI(-1, 3), QI(5, -2), QI(4, 4)], QI(0, 1))]
    for name, a, b, th in specials:
        y = yvec([Fr(x) for x in a], b, th)
        co = detpoly(y)
        rH = qrank(Hm(y))
        nz = [e for e in range(15) if co[e] != 0]
        info("B7 special point (%s)" % name, "rank H(y) %d, nonzero h-degrees of det(H(y)+hH(kappa)) %s, q %s, C %s"
             % (rH, nz, qfun(y), L.cubic_C(y)))
    print("   section B time %.1f s" % (time.time() - t0))


# ============================================================================ Section C
def perfect_matchings(pts):
    if not pts:
        yield []
        return
    a = pts[0]
    for k in range(1, len(pts)):
        rest = pts[1:k] + pts[k + 1:]
        for m in perfect_matchings(rest):
            yield [(a, pts[k])] + m


def crosses(p, q):
    a, b = sorted(p); c, d = sorted(q)
    return a < c < b < d or c < a < d < b


def section_C():
    print("\n=== C. Lemma C.10: Vandermonde weights (audit Section 5.1) ===")
    t0 = time.time()
    allm = list(perfect_matchings(list(range(8))))
    NC = [m for m in allm if not any(crosses(p, q) for p, q in itertools.combinations(m, 2))]
    check("C1 105 perfect matchings, 14 noncrossing; each noncrossing one joins odd to even positions",
          len(allm) == 105 and len(NC) == 14 and all((i + j) % 2 == 1 for m in NC for i, j in m))
    subsets = list(itertools.combinations(range(8), 4))
    Mrows = [[1 if all((i in S) != (j in S) for i, j in m) else 0 for m in NC] for S in subsets]
    check("C2 compatibility matrix 70 x 14 has rank 14 (so the weights are unique)", qrank(Mrows) == 14)
    ctx = flint.fmpq_mpoly_ctx.get(("s", 7), "lex")
    s = ctx.gens()
    ZERO = 0 * s[0]
    t = [ZERO]
    for k in range(7):
        t.append(t[-1] + s[k])

    def vand(S):
        p = ZERO + 1
        for i, j in itertools.combinations(sorted(S), 2):
            p = p * (t[j] - t[i])
        return p
    rhs = [vand(S) * vand([k for k in range(8) if k not in S]) for S in subsets]
    chosen, rows = [], []
    for i in reversed(range(70)):
        if qrank(rows + [Mrows[i]]) > len(rows):
            rows.append(Mrows[i]); chosen.append(i)
    Minv = L.inv(rows)
    wts = []
    for j in range(14):
        p = ZERO
        for k in range(14):
            if Minv[j][k] != 0:
                p = p + fq(Minv[j][k]) * rhs[chosen[k]]
        wts.append(p)
    ok = all(sum((wts[j] for j in range(14) if Mrows[i][j]), ZERO) == rhs[i] for i in range(70))
    check("C3 all 70 identities Vand(t_S) Vand(t_S^c) = sum_{pi compatible with S} w_pi hold as polynomials in the "
          "gaps (rows chosen in reverse order, independently of the audit)", ok)
    nonneg = all(w != ZERO and all(c >= 0 for c in w.coeffs()) for w in wts)
    check("C4 every w_pi is a nonzero polynomial with nonnegative coefficients in the gaps", nonneg,
          "terms %s" % [len(w.coeffs()) for w in wts])
    # 50-digit end-to-end test of the use in Lemma C.10
    mpmath.mp.dps = 50
    wd = [w.to_dict() for w in wts]

    def wval(j, gaps):
        tot = mpmath.mpf(0)
        for ex, c in wd[j].items():
            term = mpmath.mpf(int(c.p)) / int(c.q)
            for e, gg in zip(ex, gaps):
                if e:
                    term *= gg ** int(e)
            tot += term
        return tot
    rr = random.Random(424242)
    worst, worst_orth, minent, maxd = mpmath.mpf(0), mpmath.mpf(0), mpmath.mpf(10), mpmath.mpf(0)
    ntests = 0
    for trial in range(6):
        tt = sorted(mpmath.mpf(rr.uniform(-5, 5)) for _ in range(8))
        if trial == 5:
            tt[3] = tt[2] + mpmath.mpf("1e-6")
            tt = sorted(tt)
        for parity in (0, 1):
            S0 = [k for k in range(8) if k % 2 == parity]
            S0c = [k for k in range(8) if k % 2 != parity]
            cc = []
            for k in range(8):
                chi = mpmath.mpf(1)
                for l in range(8):
                    if l != k:
                        chi *= (tt[k] - tt[l])
                cc.append(mpmath.sqrt(1 / abs(chi)))
            M = mpmath.matrix(4, 8)
            for k in range(8):
                for r in range(4):
                    M[r, k] = cc[k] * tt[k] ** r
            MS0 = mpmath.matrix([[M[r, k] for k in S0] for r in range(4)])
            MS0c = mpmath.matrix([[M[r, k] for k in S0c] for r in range(4)])
            O = MS0 ** -1 * MS0c
            E = O * O.T - mpmath.eye(4)
            worst_orth = max(worst_orth, max(abs(E[i, j]) for i in range(4) for j in range(4)))
            minent = min(minent, min(abs(O[i, j]) for i in range(4) for j in range(4)))
            dd = mpmath.det(mpmath.matrix([[O[1, c] * O[2, c] for c in range(3)], [O[1, c] * O[3, c] for c in range(3)],
                                           [O[2, c] * O[3, c] for c in range(3)]]))
            maxd = max(maxd, abs(dd))
            gaps = [tt[k + 1] - tt[k] for k in range(7)]
            D = mpmath.mpf(1)
            for S in (S0, S0c):
                for i, j in itertools.combinations(S, 2):
                    D *= (tt[j] - tt[i])
            p = {}
            for jn, m in enumerate(NC):
                sig = [None] * 4
                for (i, j) in m:
                    ii, jj = (i, j) if i in S0 else (j, i)
                    sig[S0.index(ii)] = S0c.index(jj)
                p[tuple(sig)] = wval(jn, gaps) / D
            for (I, J) in ALL_PAIRS:
                lhs = sum((pv for sg, pv in p.items() if {sg[i - 1] + 1 for i in I} == set(J)), mpmath.mpf(0))
                if I:
                    rhs_ = mpmath.det(mpmath.matrix([[O[i - 1, j - 1] for j in J] for i in I])) ** 2
                else:
                    rhs_ = mpmath.mpf(1)
                worst = max(worst, abs(lhs - rhs_))
            ntests += 1
    check("C5 50-digit end-to-end test (%d Cauchy orthogonal matrices from interlaced nodes, both parities, one with a "
          "node gap 1e-6): O O^T = I, d(O) = 0, and sum_{sigma(I)=J} p_sigma = det O[I,J]^2 for all 70 (I,J) with "
          "p_sigma = w_pi / (Vand(S0) Vand(S0^c))" % ntests,
          worst < mpmath.mpf("1e-30") and worst_orth < mpmath.mpf("1e-30") and maxd < mpmath.mpf("1e-30"),
          "max identity error %s, max |OO^T - I| %s, max |d| %s, min |O_ij| %s"
          % (mpmath.nstr(worst, 3), mpmath.nstr(worst_orth, 3), mpmath.nstr(maxd, 3), mpmath.nstr(minent, 3)))
    print("   section C time %.1f s" % (time.time() - t0))


# ============================================================================ Section D
def dfun(O):
    return gdet([[O[1][c] * O[2][c] for c in range(3)], [O[1][c] * O[3][c] for c in range(3)],
                 [O[2][c] * O[3][c] for c in range(3)]])


def cayley3(p, q, r):
    K = [[Fr(0), -p, q], [p, Fr(0), -r], [-q, r, Fr(0)]]
    I3 = L.eye(3)
    return L.matmul([[I3[i][j] - K[i][j] for j in range(3)] for i in range(3)],
                    L.inv([[I3[i][j] + K[i][j] for j in range(3)] for i in range(3)]))


def stereo(s1, s2):
    D = 1 + s1 * s1 + s2 * s2
    return [2 * s1 / D, 2 * s2 / D, (s1 * s1 + s2 * s2 - 1) / D]


def O_zero00(rng):
    s1, s2 = Fr(rng.randint(-9, 9), rng.randint(1, 7)), Fr(rng.randint(-9, 9), rng.randint(1, 7))
    u = stereo(s1, s2)
    Hu = [[Fr(0)] + u] + [[u[i]] + [Fr(int(i == j)) - u[i] * u[j] for j in range(3)] for i in range(3)]
    Q = cayley3(*[Fr(rng.randint(-9, 9), rng.randint(1, 7)) for _ in range(3)])
    Dm = [[Fr(1), Fr(0), Fr(0), Fr(0)]] + [[Fr(0)] + [-Q[i][j] for j in range(3)] for i in range(3)]
    return L.matmul(Dm, Hu)


def is_SO4(O):
    return L.matmul(O, L.transpose(O)) == L.eye(4) and L.det(O) == 1


def section_D():
    print("\n=== D. (C.13), (C.14), density step (audit Section 5.2) ===")
    t0 = time.time()
    rng = random.Random(1618)
    ok = True
    signs = {}
    for _ in range(10):
        O = O_zero00(rng)
        if not (is_SO4(O) and O[0][0] == 0):
            ok = False; continue
        if dfun(O) != O[0][1] * O[0][2] * O[0][3] * O[1][0] * O[2][0] * O[3][0]:
            ok = False
        for i in range(4):
            for j in range(4):
                P = [[Fr(0)] * 4 for _ in range(4)]
                pr = list(range(4)); pr[0], pr[i] = pr[i], pr[0]
                pc = list(range(4)); pc[0], pc[j] = pc[j], pc[0]
                Op = [[O[pr[a]][pc[b]] for b in range(4)] for a in range(4)]
                if L.det(Op) != 1:
                    k = [x for x in range(4) if x != i][0]
                    Op[k] = [-x for x in Op[k]]
                assert Op[i][j] == 0 and is_SO4(Op)
                prod = Fr(1)
                for c in range(4):
                    if c != j:
                        prod *= Op[i][c]
                for r in range(4):
                    if r != i:
                        prod *= Op[r][j]
                dv = dfun(Op)
                if prod == 0:
                    if dv != 0:
                        ok = False
                elif dv not in (prod, -prod):
                    ok = False
                else:
                    signs.setdefault((i, j), set()).add(1 if dv == prod else -1)
    info("D1a", "samples with an extra zero entry are checked as d = 0 only")
    check("D1 (C.13) exactly at 10 random rational O in SO(4) with O00 = 0; transported to every position (i,j): "
          "d(O) = +-(product of the other six entries of row i and column j) with a sign depending only on (i,j)",
          ok and all(len(v) == 1 for v in signs.values()) and len(signs) == 16,
          "signs %s" % {k: list(v)[0] for k, v in sorted(signs.items())})
    # D2: (C.13) as an identity on all of SO(4) ∩ {O00 = 0}, via the bijective rational parametrization
    #     O = diag(1, -Q) H(u), u in S^2 (stereographic), Q in SO(3) (Cayley).
    ctx = flint.fmpq_mpoly_ctx.get(("s1", "s2", "p", "q", "r"), "lex")
    s1, s2, p, q, r = ctx.gens()
    one = 0 * s1 + 1
    Du = one + s1 ** 2 + s2 ** 2
    un = [2 * s1, 2 * s2, s1 ** 2 + s2 ** 2 - 1]
    HuN = [[0 * s1] + [un[i] * Du for i in range(3)]] + \
          [[un[i] * Du] + [(Du ** 2 if i == j else 0 * s1) - un[i] * un[j] for j in range(3)] for i in range(3)]
    K = [[0 * s1, -p, q], [p, 0 * s1, -r], [-q, r, 0 * s1]]
    IpK = [[(one if i == j else 0 * s1) + K[i][j] for j in range(3)] for i in range(3)]
    ImK = [[(one if i == j else 0 * s1) - K[i][j] for j in range(3)] for i in range(3)]
    adj = [[(-1) ** (i + j) * gdet([[IpK[a][b] for b in range(3) if b != i] for a in range(3) if a != j])
            for j in range(3)] for i in range(3)]
    DQ = one + p ** 2 + q ** 2 + r ** 2
    Nq = [[sum((ImK[i][k] * adj[k][j] for k in range(3)), 0 * s1) for j in range(3)] for i in range(3)]
    Dm = [[DQ, 0 * s1, 0 * s1, 0 * s1]] + [[0 * s1] + [-Nq[i][j] for j in range(3)] for i in range(3)]
    NO = [[sum((Dm[i][k] * HuN[k][j] for k in range(4)), 0 * s1) for j in range(4)] for i in range(4)]
    den = DQ * Du ** 2
    orth = all(sum((NO[i][k] * NO[j][k] for k in range(4)), 0 * s1) == (den ** 2 if i == j else 0 * s1)
               for i in range(4) for j in range(4))
    dN = dfun(NO)
    prodN = NO[0][1] * NO[0][2] * NO[0][3] * NO[1][0] * NO[2][0] * NO[3][0]
    check("D2 (C.13) as a polynomial identity on the parametrization O = diag(1,-Q) H(u) of SO(4) ∩ {O00 = 0} "
          "(bijective onto it; dense rational charts), so it holds on all of SO(4) ∩ {O00 = 0} by continuity",
          orth and NO[0][0] == 0 * s1 and dN == prodN, "orthogonality %s, terms in d %d" % (orth, len(dN.coeffs())))
    # D3: (C.14) with own parametrization: row0 = (a, b, 0, 0), rows (-b x_i, a x_i, y_i, z_i), [x y z] = Q
    ctx2 = flint.fmpq_mpoly_ctx.get(("tau", "p", "q", "r", "e"), "lex")
    tau, p, q, r, e = ctx2.gens()
    z0 = 0 * tau
    one = z0 + 1
    Dt = one + tau ** 2
    K = [[z0, -p, q], [p, z0, -r], [-q, r, z0]]
    IpK = [[(one if i == j else z0) + K[i][j] for j in range(3)] for i in range(3)]
    ImK = [[(one if i == j else z0) - K[i][j] for j in range(3)] for i in range(3)]
    adj = [[(-1) ** (i + j) * gdet([[IpK[a][b] for b in range(3) if b != i] for a in range(3) if a != j])
            for j in range(3)] for i in range(3)]
    DQ = one + p ** 2 + q ** 2 + r ** 2
    Nq = [[sum((ImK[i][k] * adj[k][j] for k in range(3)), z0) for j in range(3)] for i in range(3)]
    Na, Nb = (1 - tau ** 2) * DQ, 2 * tau * DQ
    NO = [[Na, Nb, z0, z0]] + [[-2 * tau * Nq[i][0], (1 - tau ** 2) * Nq[i][0], Dt * Nq[i][1], Dt * Nq[i][2]]
                               for i in range(3)]
    lam = Dt * DQ
    orth = all(sum((NO[i][k] * NO[j][k] for k in range(4)), z0) == (lam ** 2 if i == j else z0)
               for i in range(4) for j in range(4))
    detok = gdet(NO) == lam ** 4
    NX = [Dt * Nq[i][0] for i in range(3)]
    ratios = []
    for i in (1, 2, 3):
        E = [[0] * 4 for _ in range(4)]
        E[0][i] = 1; E[i][0] = -1
        EO = [[sum((E[a][k] * NO[k][b] for k in range(4)), z0) for b in range(4)] for a in range(4)]
        Pe = dfun([[NO[a][b] + e * EO[a][b] for b in range(4)] for a in range(4)])
        lin = ctx2.from_dict({ex[:4] + (0,): c for ex, c in Pe.to_dict().items() if ex[4] == 1})
        jj, kk = [x for x in (1, 2, 3) if x != i]
        tgt = Na * Nb * NX[jj - 1] * NX[kk - 1] * NO[i][2] * NO[i][3]
        ratios.append(1 if lin == tgt else (-1 if lin == -tgt else 0))
    check("D3 (C.14) as a polynomial identity (own parametrization: stereographic (a,b), Cayley [x y z]): "
          "dd(E_i O) = s * a b x_j x_k y_i z_i with one common sign s for i = 1,2,3", orth and detok
          and len(set(ratios)) == 1 and ratios[0] != 0, "orth %s, det 1 %s, s = %s" % (orth, detok, ratios))
    # D4: numeric density near an indecomposable zero-entry point of {d = 0}
    a0, b0 = 3 / 5, 4 / 5
    Q = np.array([[float(x) for x in row] for row in cayley3(Fr(1, 2), Fr(-2, 3), Fr(3, 4))])
    O = np.vstack([[a0, b0, 0, 0]] + [[-b0 * Q[i, 0], a0 * Q[i, 0], Q[i, 1], Q[i, 2]] for i in range(3)])

    def dnp(M):
        return np.linalg.det(np.array([[M[1, c] * M[2, c] for c in range(3)], [M[1, c] * M[3, c] for c in range(3)],
                                       [M[2, c] * M[3, c] for c in range(3)]]))
    basis = []
    for a in range(4):
        for b in range(a + 1, 4):
            X = np.zeros((4, 4)); X[a, b] = 1; X[b, a] = -1; basis.append(X)

    def expm_skew(X):
        w, V = np.linalg.eig(X)
        return np.real(V @ np.diag(np.exp(w)) @ np.linalg.inv(V))

    def grad(M, h=1e-6):
        return np.array([(dnp(expm_skew(h * X) @ M) - dnp(expm_skew(-h * X) @ M)) / (2 * h) for X in basis])
    g0 = grad(O)
    rr = np.random.default_rng(11)
    res = []
    for _ in range(20):
        v = rr.standard_normal(6); v -= (v @ g0) / (g0 @ g0) * g0
        v *= 1e-3 / np.linalg.norm(v)
        M = expm_skew(sum(c * X for c, X in zip(v, basis))) @ O
        for _ in range(30):
            gg = grad(M)
            step = -dnp(M) / (gg @ gg) * gg
            M = expm_skew(sum(c * X for c, X in zip(step, basis))) @ M
        res.append((abs(dnp(M)), np.min(np.abs(M)), np.linalg.norm(M - O), int(np.argmin(np.abs(M)))))
    info("D4a", "per-run min |entry| (position): %s" % ", ".join("%.1e(%d)" % (x[1], x[3]) for x in res))
    check("D4 NUMERIC: from the indecomposable point O (row (3/5,4/5,0,0), d = 0, |grad d| = %.3g) 20 random tangent "
          "steps of size 1e-3 projected back to {d = 0} land on points with all 16 entries nonzero (|entry| > 1e-12, far above rounding; the entry O02 grows only quadratically: H is tangent to {O02 = 0} at O)" % np.linalg.norm(g0),
          all(x[0] < 1e-12 for x in res) and all(x[1] > 1e-12 for x in res) and np.linalg.norm(g0) > 1e-8,
          "max |d| %.1e, runs with min |entry| > 1e-12: %d of 20, max distance to O %.2e"
          % (max(x[0] for x in res), sum(x[1] > 1e-12 for x in res), max(x[2] for x in res)))
    print("   section D time %.1f s" % (time.time() - t0))


# ============================================================================ Section E
CH_PAIRS, CH_NINV = None, None


def setup_rho_eval():
    global CH_PAIRS, CH_NINV
    rows, ch = [], []
    for pr in ALL_PAIRS:
        nv = L.n_vec(*pr)
        if qrank(rows + [nv]) > len(rows):
            rows.append(nv); ch.append(pr)
    CH_PAIRS, CH_NINV = ch, L.inv(rows)


def rho_eval(V):
    """rho~(V) from the complementary-minor definition (certified against A1/A2); any ring."""
    t = []
    for (I, J) in CH_PAIRS:
        val = gdet(sub(V, I, J)) * gdet(sub(V, comp(I), comp(J)))
        t.append(val if (sum(I) + sum(J)) % 2 == 0 else -val)
    return [sum((CH_NINV[j][k] * t[k] for k in range(14)), 0) for j in range(14)]


def mm(A, B):
    return [[sum((A[i][k] * B[k][j] for k in range(len(B))), 0) for j in range(len(B[0]))] for i in range(len(A))]


def madd(*Ms):
    return [[sum((M[i][j] for M in Ms), 0) for j in range(len(Ms[0][0]))] for i in range(len(Ms[0]))]


def mscale(c, M):
    return [[c * x for x in r] for r in M]


def vinv_rows(D):
    return L.inv([[Fr(h) ** e for e in range(D + 1)] for h in range(D + 1)])


VI4, VI8 = None, None


def taylor_coef(fun, D, k):
    VI = VI4 if D == 4 else VI8
    vals = [fun(Fr(h)) for h in range(D + 1)]
    return [sum((VI[k][h] * vals[h][c] for h in range(D + 1)), 0) for c in range(14)]


def EYE4():
    return [[Fr(int(i == j)) for j in range(4)] for i in range(4)]


def h3_data(O):
    """returns (gauge rank, rank[rho0,A], rank[rho0,A,B], rho0, Acols, Bmat, Zs)."""
    I4 = EYE4()
    G = []
    for i in range(4):
        D = [[Fr(0)] * 4 for _ in range(4)]; D[i][i] = Fr(1); G.append(D)
        Ec = [[O[a][i] * O[b][i] for b in range(4)] for a in range(4)]; G.append(Ec)
    flat = lambda M: [M[a][b] for a in range(4) for b in range(a, 4)]
    Grows = [flat(M) for M in G]
    rg = crank(Grows)
    Zs = []
    cur = list(Grows)
    for a in range(4):
        for b in range(a + 1, 4):
            Z = [[Fr(0)] * 4 for _ in range(4)]; Z[a][b] = Z[b][a] = Fr(1)
            if crank(cur + [flat(Z)]) > crank(cur):
                cur.append(flat(Z)); Zs.append(Z)
            if len(Zs) == 3:
                break
        if len(Zs) == 3:
            break
    rho0 = rho_eval(O)
    Acols = []
    for a in range(4):
        for b in range(a + 1, 4):
            X = [[Fr(0)] * 4 for _ in range(4)]; X[a][b] = Fr(1); X[b][a] = Fr(-1)
            XO = mm(X, O)
            Acols.append(taylor_coef(lambda s: rho_eval(madd(O, mscale(s, XO))), 4, 1))

    def Bdiag(Z):
        ZO, Z2O = mm(Z, O), mm(mm(Z, Z), O)
        return taylor_coef(lambda s: rho_eval(madd(O, mscale(s, ZO), mscale(s * s / 2, Z2O))), 8, 2)
    Bd = [Bdiag(Z) for Z in Zs]
    Bm = {}
    for i in range(3):
        Bm[(i, i)] = Bd[i]
    for i in range(3):
        for j in range(i + 1, 3):
            Bs = Bdiag(madd(Zs[i], Zs[j]))
            Bm[(i, j)] = [(Bs[c] - Bd[i][c] - Bd[j][c]) / 2 for c in range(14)]
    Bcols = [Bm[k] for k in sorted(Bm)]
    r7 = crank([rho0] + Acols)
    r13 = crank([rho0] + Acols + Bcols)
    return rg, r7, r13, rho0, Acols, Bm, len(Zs)


def indecomposable(O):
    edges = [(i, j) for i in range(4) for j in range(4) if O[i][j] != 0]
    seen, stack = {("r", 0)}, [("r", 0)]
    while stack:
        t, k = stack.pop()
        for (i, j) in edges:
            nb = ("c", j) if (t == "r" and i == k) else (("r", i) if (t == "c" and j == k) else None)
            if nb and nb not in seen:
                seen.add(nb); stack.append(nb)
    return len(seen) == 8


def c8_alpha(O, rho0, Acols, Bm):
    """C|K_O = alpha det S_O ?  returns (holds, alpha, dimK)."""
    Kb = L.nullspace([rho0] + Acols)
    cs = sp.symbols("k0:%d" % len(Kb))
    R = lambda x: sp.Rational(x.numerator, x.denominator)
    u = [sum(R(Kb[i][k]) * cs[i] for i in range(len(Kb))) for k in range(14)]
    CK = sp.Poly(sp.expand(L.cubic_C(u)), *cs)
    S = sp.Matrix(3, 3, lambda i, j: sum(R(Bm[(min(i, j), max(i, j))][k]) * u[k] for k in range(14)))
    DS = sp.Poly(sp.expand(S.det()), *cs)
    if DS.is_zero:
        return (CK.is_zero, None, len(Kb))
    mon, cD = next((m, c) for m, c in zip(DS.monoms(), DS.coeffs()) if c != 0)
    alpha = CK.coeff_monomial(mon) / cD
    return (sp.expand(CK.as_expr() - alpha * DS.as_expr()) == 0, alpha, len(Kb))


def frame_det(O):
    """det of [7 independent gauge generators; F01, F02, F03] in Sym(4) (10x10), or None."""
    G = []
    for i in range(4):
        D = [[Fr(0)] * 4 for _ in range(4)]; D[i][i] = Fr(1); G.append(D)
        G.append([[O[a][i] * O[b][i] for b in range(4)] for a in range(4)])
    flat = lambda M: [M[a][b] for a in range(4) for b in range(a, 4)]
    sel = []
    for r in [flat(M) for M in G]:
        if qrank(sel + [r]) > len(sel):
            sel.append(r)
    for (a, b) in ((0, 1), (0, 2), (0, 3)):
        Z = [[Fr(0)] * 4 for _ in range(4)]; Z[a][b] = Z[b][a] = Fr(1); sel.append(flat(Z))
    return L.det(sel) if len(sel) == 10 else None


def rational_SO4(rng):
    K = [[Fr(0)] * 4 for _ in range(4)]
    for a in range(4):
        for b in range(a + 1, 4):
            K[a][b] = Fr(rng.randint(-5, 5), rng.randint(1, 4)); K[b][a] = -K[a][b]
    I4 = EYE4()
    return mm(madd(I4, mscale(-1, K)), L.inv(madd(I4, K)))


def section_E():
    global VI4, VI8
    print("\n=== E. (H3) and (C.8) at special stable points (audit Sections 2, 3) ===")
    t0 = time.time()
    setup_rho_eval()
    VI4, VI8 = vinv_rows(4), vinv_rows(8)
    Oc9 = [[Fr(0), Fr(2, 3), Fr(2, 3), Fr(-1, 3)], [Fr(2, 3), Fr(5, 9), Fr(-4, 9), Fr(2, 9)],
           [Fr(2, 3), Fr(-4, 9), Fr(5, 9), Fr(2, 9)], [Fr(1, 3), Fr(-2, 9), Fr(-2, 9), Fr(-8, 9)]]
    K3 = [[0, 4, 3, 2], [-4, 0, 3, -3], [-3, -3, 0, -3], [-2, 3, 3, 0]]
    K = [[Fr(x, 3) for x in r] for r in K3]
    I4 = EYE4()
    Ocay = mm(madd(I4, mscale(-1, K)), L.inv(madd(I4, K)))
    rng = random.Random(5150)
    Orand = rational_SO4(rng)
    Oz = O_zero00(random.Random(77))
    Q = cayley3(Fr(1, 2), Fr(-2, 3), Fr(3, 4))
    Ocd = [[Fr(3, 5), Fr(4, 5), Fr(0), Fr(0)]] + [[Fr(-4, 5) * Q[i][0], Fr(3, 5) * Q[i][0], Q[i][1], Q[i][2]]
                                                  for i in range(3)]
    # complex orthogonal (Gaussian rational) via Cayley with complex skew K
    Kc = [[QI(0)] * 4 for _ in range(4)]
    vals = [QI(1, 2), QI(-1, 1), QI(2, -1), QI(0, 3), QI(1, 1), QI(-2, 1)]
    idx = 0
    for a in range(4):
        for b in range(a + 1, 4):
            Kc[a][b] = vals[idx]; Kc[b][a] = QI(0) - vals[idx]; idx += 1
    IK = [[QI(int(i == j)) + Kc[i][j] for j in range(4)] for i in range(4)]
    # inverse of IK by Gauss-Jordan over QI
    n = 4
    aug = [IK[i][:] + [QI(int(i == j)) for j in range(4)] for i in range(4)]
    for c in range(n):
        k = next(i for i in range(c, n) if not aug[i][c] == QI(0))
        aug[c], aug[k] = aug[k], aug[c]
        piv = aug[c][c]
        aug[c] = [x / piv for x in aug[c]]
        for i in range(n):
            if i != c and not aug[i][c] == QI(0):
                f = aug[i][c]
                aug[i] = [x - f * y for x, y in zip(aug[i], aug[c])]
    IKi = [r[4:] for r in aug]
    Ocx = mm([[QI(int(i == j)) - Kc[i][j] for j in range(4)] for i in range(4)], IKi)
    OOt = mm(Ocx, [list(r) for r in zip(*Ocx)])
    cx_ok = all(OOt[i][j] == QI(int(i == j)) for i in range(4) for j in range(4))
    cases = [("Lemma C.9 example (O00 = 0, on 7 Pluecker divisors)", Oc9, True),
             ("part-1 Cayley matrix (all entries nonzero)", Ocay, True),
             ("random rational SO(4) (own seed)", Orand, True),
             ("random O with O00 = 0 (own generator)", Oz, True),
             ("indecomposable zero-entry point on the Cauchy divisor, row (3/5,4/5,0,0)", Ocd, True),
             ("complex Gaussian-rational O in SO(4,C)", Ocx, False)]
    for name, O, real in cases:
        dval = dfun(O)
        rg, r7, r13, rho0, Acols, Bm, nz = h3_data(O)
        ok = indecomposable(O) and rg == 7 and nz == 3 and r7 == 7 and r13 == 13
        extra = ""
        if real:
            extra = ", det O = %s" % L.det(O)
            holds, alpha, dK = c8_alpha(O, rho0, Acols, Bm)
            extra += ", (C.8) holds: %s, dim K = %d, alpha = %s" % (holds, dK, alpha)
            ok = ok and holds and dK == 7 and ((alpha != 0) == (dval != 0))
            fd = frame_det(O)
            if fd is not None and fd != 0:
                normal = sp.Rational(fd.numerator, fd.denominator) ** 2 * alpha
                extra += ", frame det = %s (= -d: %s), alpha * framedet^2 = %s (= d: %s)" % (
                    fd, fd == -dval, normal, normal == sp.Rational(dval.numerator, dval.denominator))
            else:
                extra += ", frame F01,F02,F03 degenerate here (frame det %s)" % fd
        check("E1 %s: indecomposable, gauge dim 7, rank[rho0,A] = 7, rank[rho0,A,B] = 13 (H3)%s"
              % (name, "; (C.8) holds with alpha != 0 iff d != 0" if real else ""),
              ok and (real or cx_ok), "d = %s%s" % (dval if real else "complex", extra))
    # NUMERIC scan of real SO(4) for near-failures of (H3)
    rnp = np.random.default_rng(2026)
    PAIRS = CH_PAIRS
    NINV = np.array([[float(x) for x in r] for r in CH_NINV])
    sgnp = np.array([(-1) ** (sum(I) + sum(J)) for (I, J) in PAIRS], dtype=float)

    def rho_np(V):
        t = np.empty(14, dtype=complex)
        for k, (I, J) in enumerate(PAIRS):
            d1 = np.linalg.det(V[np.ix_([i - 1 for i in I], [j - 1 for j in J])]) if I else 1.0
            Ic, Jc = comp(I), comp(J)
            d2 = np.linalg.det(V[np.ix_([i - 1 for i in Ic], [j - 1 for j in Jc])]) if Ic else 1.0
            t[k] = sgnp[k] * d1 * d2
        return NINV @ t
    Nfft = 12
    om = np.exp(2j * np.pi * np.arange(Nfft) / Nfft)

    def coef(fun, k):
        vals = np.array([fun(s) for s in om])
        return np.real((vals * om[:, None] ** (-k)).sum(axis=0) / Nfft)
    Xs = []
    for a in range(4):
        for b in range(a + 1, 4):
            X = np.zeros((4, 4)); X[a, b] = 1; X[b, a] = -1; Xs.append(X)
    symb = []
    for a in range(4):
        for b in range(a, 4):
            S = np.zeros((4, 4)); S[a, b] = S[b, a] = 1.0; symb.append(S / np.linalg.norm(S))

    def sigma13(O):
        G = [np.diag(np.eye(4)[i]) for i in range(4)] + [np.outer(O[:, i], O[:, i]) for i in range(4)]
        Gv = np.array([[np.sum(Gm * S) for S in symb] for Gm in G])
        U, sv, Vt = np.linalg.svd(Gv)
        Wb = Vt[7:]                                   # orthonormal complement of the gauge in Sym(4)
        Zs = [sum(c * S for c, S in zip(w, symb)) for w in Wb]
        cols = [np.real(rho_np(O.astype(complex)))]
        for X in Xs:
            cols.append(coef(lambda s: rho_np(O + s * X @ O), 1))
        Bd = [coef(lambda s: rho_np(O + s * Z @ O + s * s / 2 * Z @ Z @ O), 2) for Z in Zs]
        Bc = list(Bd)
        for i in range(3):
            for j in range(i + 1, 3):
                Zij = Zs[i] + Zs[j]
                Bij = coef(lambda s: rho_np(O + s * Zij @ O + s * s / 2 * Zij @ Zij @ O), 2)
                Bc.append((Bij - Bd[i] - Bd[j]) / 2)
        M = np.array(cols + Bc).T
        M = M / np.linalg.norm(M, axis=0)
        return np.linalg.svd(M, compute_uv=False)[-1], sv[6] / sv[0]

    def haar():
        Z = rnp.standard_normal((4, 4))
        Qm, R = np.linalg.qr(Z)
        Qm = Qm @ np.diag(np.sign(np.diag(R)))
        if np.linalg.det(Qm) < 0:
            Qm[:, 0] *= -1
        return Qm
    samples = []
    for _ in range(1000):
        O = haar()
        s13, gsv = sigma13(O)
        samples.append((s13, np.min(np.abs(O)), gsv))
    samples.sort()
    bins = {}
    for s13, me, gsv in samples:
        key = "minentry>=0.1" if me >= 0.1 else ("0.02<=minentry<0.1" if me >= 0.02 else "minentry<0.02")
        bins.setdefault(key, []).append(s13)
    info("E2 NUMERIC Haar scan of 1000 real O in SO(4) (seed 2026): smallest singular value of the column-normalised "
         "14x13 matrix [rho0, A, B]", "; ".join("%s: n=%d min %.3e median %.3e" % (k, len(v), min(v), np.median(v))
                                                for k, v in sorted(bins.items())))
    # local minimisation of sigma13 with a penalty keeping all entries >= 0.1 in modulus
    from scipy.optimize import minimize

    def expm_skew(X):
        w, V = np.linalg.eig(X)
        return np.real(V @ np.diag(np.exp(w)) @ np.linalg.inv(V))
    best = []
    starts = [O for O in [haar() for _ in range(400)] if np.min(np.abs(O)) >= 0.15][:8]
    for O0 in starts:
        def f(th):
            Om = expm_skew(sum(c * X for c, X in zip(th, Xs))) @ O0
            me = np.min(np.abs(Om))
            return sigma13(Om)[0] + (10 * (0.1 - me) if me < 0.1 else 0.0)
        r = minimize(f, np.zeros(6), method="Nelder-Mead", options={"maxiter": 200, "xatol": 1e-6, "fatol": 1e-9})
        best.append(r.fun)
    info("E3 NUMERIC local minimisation of that singular value over real SO(4) with all |O_ij| >= 0.1 (8 starts, "
         "Nelder-Mead, 200 iterations)", "minima %s" % ", ".join("%.3e" % b for b in sorted(best)))
    print("   section E time %.1f s" % (time.time() - t0))


if __name__ == "__main__":
    which = sys.argv[1:] or ["A", "B", "C", "D", "E"]
    T0 = time.time()
    for sec in which:
        {"A": section_A, "B": section_B, "C": section_C, "D": section_D, "E": section_E}[sec]()
    print("\nSUMMARY: %d failure(s)%s; total time %.1f s" % (len(FAILS), (": " + "; ".join(FAILS)) if FAILS else "",
                                                           time.time() - T0))
    sys.exit(1 if FAILS else 0)
