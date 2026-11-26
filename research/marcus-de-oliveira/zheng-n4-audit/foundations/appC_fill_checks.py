"""Exact checks for appendix-c-open-items.md (auditor appCfill).

Own code only. Reads nothing from the supplement. Uses the shared library ../lib/moc4lib.py
(matching basis (A.1), P, kappa, the cubic C of (A.2), minor vectors n_{I,J}).

Sections (each prints PASS/FAIL lines; the script exits 1 on any FAIL):
  A  Proposition C.4: conormal spaces of M8' lie in the secant cone, and dim (M8')* = 11.
  B  Lemma C.11 without (C.2): the h^3 coefficient, uniqueness on the kappa-line, smoothness of y.
  C  Lemma C.10 without Lam: nonnegative double-dimer weights for Vandermonde points.
  D  (C.13) and (C.14) as identities (density step of Lemma C.10).
Run:  python -I -B appC_fill_checks.py      (deterministic; about a minute)
"""
import itertools
import os
import random
import sys
import time
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L  # noqa: E402
import sympy as sp  # noqa: E402
import flint  # noqa: E402

FAILS = []


def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + (("  [" + info + "]") if info else ""))
    if not ok:
        FAILS.append(name)


def qmat(rows):
    return flint.fmpq_mat([[flint.fmpq(int(Fr(x).numerator), int(Fr(x).denominator)) for x in r]
                           for r in rows])


def qrank(rows):
    return qmat(rows).rank()


def R(x):
    x = Fr(x)
    return sp.Rational(x.numerator, x.denominator)


def fr(x):
    """sympy Rational / flint fmpq / int -> Fraction."""
    if isinstance(x, Fr):
        return x
    if isinstance(x, int):
        return Fr(x)
    if isinstance(x, flint.fmpq):
        return Fr(int(x.p), int(x.q))
    x = sp.Rational(x)
    return Fr(int(x.p), int(x.q))


# ---------------------------------------------------------------------------------------------
# Common symbolic objects
# ---------------------------------------------------------------------------------------------
VS = sp.symbols("v0:16")
VM = sp.Matrix(4, 4, VS)
ZS = sp.symbols("z1:9")
MS, CSYM = L.cubic_sympy()


def minor(M, I, J):
    if len(I) == 0:
        return sp.Integer(1)
    return M.extract([i - 1 for i in I], [j - 1 for j in J]).det()


def comp(I):
    return tuple(k for k in range(1, 5) if k not in I)


ALL_PAIRS = [(I, J) for k in range(5) for I in itertools.combinations(range(1, 5), k)
             for J in itertools.combinations(range(1, 5), k)]


def build_rho_tilde():
    """rho~(V) in V_C^* (14 coordinates, dual to B_0..B_13) with
    <rho~(V), n_{I,J}> = (-1)^{sum I + sum J} det V[I,J] det V[I^c,J^c] for all 70 (I,J).
    Solved from 14 independent pairs; consistency on all 70 pairs and the identity
    <rho~(V), m(a,-b)> = det(D_a V + V D_b) are then checked as polynomial identities."""
    nvecs = {p: L.n_vec(*p) for p in ALL_PAIRS}
    chosen, rows = [], []
    for p in ALL_PAIRS:
        if L.rank(rows + [nvecs[p]]) > len(rows):
            rows.append(nvecs[p]); chosen.append(p)
        if len(rows) == 14:
            break
    Ninv = L.inv(rows)

    def target(p):
        I, J = p
        s = (-1) ** (sum(I) + sum(J))
        return sp.expand(s * minor(VM, I, J) * minor(VM, comp(I), comp(J)))

    tg = [target(p) for p in chosen]
    rho = [sp.expand(sum(sp.Rational(Ninv[j][k].numerator, Ninv[j][k].denominator) * tg[k]
                         for k in range(14))) for j in range(14)]
    ok = True
    for p in ALL_PAIRS:
        lhs = sp.expand(sum(sp.Rational(nvecs[p][j].numerator, nvecs[p][j].denominator) * rho[j]
                            for j in range(14)))
        if sp.expand(lhs - target(p)) != 0:
            ok = False
    check("A0 rho~ consistent on all 70 minor pairings (polynomial identity in 16 entries)", ok)
    a = sp.symbols("a1:5"); b = sp.symbols("b1:5")
    mab = [sp.expand(L.eval_poly(L.B_POLYS[r], list(a) + [-x for x in b])) for r in range(14)]
    lhs = sp.expand(sp.Matrix(4, 4, lambda i, j: (a[i] + b[j]) * VM[i, j]).det())
    rhs = sp.expand(sum(rho[j] * mab[j] for j in range(14)))
    check("A0 <rho~(V), m(a,-b)> = det(D_a V + V D_b) (identity in a, b, V: 24 variables)",
          sp.expand(lhs - rhs) == 0)
    return rho


def section_A():
    print("\n=== A. Proposition C.4 ===")
    t0 = time.time()
    rho = build_rho_tilde()
    rho_polys = [sp.Poly(r, *VS, domain="QQ") for r in rho]
    grad = [[p.diff(v) for v in VS] for p in rho_polys]
    m_sym = [sp.expand(L.eval_poly(L.B_POLYS[r], list(ZS))) for r in range(14)]
    m_polys = [sp.Poly(e, *ZS, domain="QQ") for e in m_sym]
    dm_polys = [[p.diff(z) for z in ZS] for p in m_polys]

    def ev(poly, vals):
        return fr(poly.eval(dict(zip(poly.gens, [R(x) for x in vals]))))

    # torus invariance: rho~(D1 V D2) = det D1 det D2 rho~(V)  (gives rank d rho~ <= 10)
    d1 = sp.symbols("p1:5"); d2 = sp.symbols("r1:5")
    Vt = sp.diag(*d1) * VM * sp.diag(*d2)
    sub = dict(zip(VS, list(Vt)))
    fac = sp.prod(d1) * sp.prod(d2)
    ok = all(sp.expand(r.subs(sub, simultaneous=True) - fac * r) == 0 for r in rho)
    check("A1 rho~(D1 V D2) = det(D1) det(D2) rho~(V) (identity)", ok)

    # rational test point (z, w) on the Segre quadric
    z = [Fr(x) for x in (0, 1, 3, -2, 5, -4, 7, 2)]
    w = [Fr(x) for x in (2, -1, 4, 3, -3, 6, 1, -5)]
    A = [[Fr(1)] * 8, z, w, [z[i] * w[i] for i in range(8)]]
    A0 = [r[:4] for r in A]; A1 = [r[4:] for r in A]
    A0i = L.inv(A0)
    V0 = L.matmul(A0i, A1)
    v0 = [V0[i][j] for i in range(4) for j in range(4)]
    check("A2 test point: all 16 entries of V0 nonzero and det V0 != 0",
          all(x != 0 for x in v0) and L.det(V0) != 0)

    # Segre rank lemma at the point (it is proved in general in the report)
    def had(zz, V):
        return [[(zz[i] - zz[4 + j]) * V[i][j] for j in range(4)] for i in range(4)]
    Mz, Mw = had(z, V0), had(w, V0)
    check("A3 rank((z_i - z_{4+j}) V0_ij) = 2 and rank((w_i - w_{4+j}) V0_ij) = 2",
          L.rank(Mz) == 2 and L.rank(Mw) == 2)

    J0 = [[ev(grad[k][a], v0) for a in range(16)] for k in range(14)]
    rJ = qrank(J0)
    check("A4 rank d rho~(V0) = 10", rJ == 10, "rank %d" % rJ)
    N = L.nullspace(L.transpose(J0))           # xi with xi^T J0 = 0
    check("A4 conormal space N(V0) = ann(im d rho~(V0)) has dimension 4", len(N) == 4)
    rho0 = [ev(r, v0) for r in rho_polys]
    check("A4 rho~(V0) in im d rho~ (Euler), so N(V0) annihilates rho~(V0)",
          all(sum(n[k] * rho0[k] for k in range(14)) == 0 for n in N))

    # dominance of (z, w) -> rho~(V(z, w)) onto the 10-dimensional cone
    cols = []
    for p in range(16):
        dA = [[Fr(0)] * 8 for _ in range(4)]
        i = p % 8
        if p < 8:
            dA[1][i] = Fr(1); dA[3][i] = w[i]
        else:
            dA[2][i] = Fr(1); dA[3][i] = z[i]
        dA0 = [r[:4] for r in dA]; dA1 = [r[4:] for r in dA]
        dV = L.matmul(A0i, [[dA1[r][c] - sum(dA0[r][k] * V0[k][c] for k in range(4))
                             for c in range(4)] for r in range(4)])
        dv = [dV[i2][j2] for i2 in range(4) for j2 in range(4)]
        cols.append([sum(J0[k][a] * dv[a] for a in range(16)) for k in range(14)])
    rdom = qrank(cols + [rho0])
    check("A5 rank of d(z,w -> rho~(V(z,w))) together with rho~(V0) = 10 (dominance)",
          rdom == 10, "rank %d" % rdom)

    mz = [ev(p, z) for p in m_polys]; mw = [ev(p, w) for p in m_polys]

    def in_span(vec, basis):
        return qrank(basis) == qrank(basis + [vec])
    check("A6 m(z) and m(w) lie in N(V0)", in_span(mz, N) and in_span(mw, N))

    def inc_tangent(zz, M0):
        # left and right kernels of the rank-2 matrix M0, then the 4 linearised conditions
        Rk = L.nullspace(M0); Lk = L.nullspace(L.transpose(M0))
        conds = []
        for l in Lk:
            for r in Rk:
                row = [Fr(0)] * 8
                for i in range(4):
                    for j in range(4):
                        cf = l[i] * V0[i][j] * r[j]
                        row[i] += cf; row[4 + j] -= cf
                conds.append(row)
        return conds, L.nullspace(conds)
    cz, Tz = inc_tangent(z, Mz); cw, Tw = inc_tangent(w, Mw)
    check("A7 Inc_V0 = {z': rank((z'_i - z'_{4+j}) V0_ij) <= 2} is cut transversally at z and at w"
          " (4 independent linear conditions; tangent dimension 4)",
          qrank(cz) == 4 and qrank(cw) == 4 and len(Tz) == 4 and len(Tw) == 4)
    dmz = [[ev(dm_polys[k][i], z) for i in range(8)] for k in range(14)]
    dmw = [[ev(dm_polys[k][i], w) for i in range(8)] for k in range(14)]
    img = [mz, mw]
    for t in Tz:
        img.append([sum(dmz[k][i] * t[i] for i in range(8)) for k in range(14)])
    for t in Tw:
        img.append([sum(dmw[k][i] * t[i] for i in range(8)) for k in range(14)])
    rj = qrank(img)
    check("A8 join differential at (z, w, 1, 1) has rank 4 = dim N(V0)", rj == 4, "rank %d" % rj)
    check("A8 the join image lies in N(V0) (consistency)", all(in_span(v, N) for v in img))

    # lower bound dim (M8')^-hat >= 12: rank of the projection E -> V_C at (V0, xi0)
    rng = random.Random(7)
    coef = [Fr(rng.randint(-5, 5)) for _ in range(4)]
    xi0 = [sum(coef[i] * N[i][k] for i in range(4)) for k in range(14)]
    hess = [[[ev(rho_polys[k].diff(VS[a]).diff(VS[b]), v0) for b in range(16)]
             for a in range(16)] for k in range(14)]
    eqs = []
    for a in range(16):
        row = [sum(xi0[k] * hess[k][a][b] for k in range(14)) for b in range(16)]
        row += [J0[k][a] for k in range(14)]
        eqs.append(row)
    TE = L.nullspace(eqs)
    proj = [t[16:] for t in TE]
    rp = qrank(proj)
    check("A9 tangent space of E = {(V, xi): xi annihilates im d rho~(V)} at (V0, xi0) has dim 20",
          len(TE) == 20, "dim %d" % len(TE))
    check("A9 projection (V, xi) -> xi has rank 12 there, so dim (M8')^-hat >= 12",
          rp == 12, "rank %d" % rp)
    print("   section A time %.1f s" % (time.time() - t0))


# ---------------------------------------------------------------------------------------------
# B. Lemma C.11 without (C.2)
# ---------------------------------------------------------------------------------------------
def section_B():
    print("\n=== B. Lemma C.11 ===")
    t0 = time.time()
    C = CSYM
    q = MS[0] * MS[6] - MS[1] * MS[5]
    dC13 = sp.expand(sp.diff(C, MS[13]))
    check("B1 dC/du13 = q(u) = u0 u6 - u1 u5 and d2C/du13^2 = 0 (so H(kappa) = Hess q and"
          " H(u) e13 = Hess(q) u)", sp.expand(dC13 - q) == 0 and sp.diff(C, MS[13], 2) == 0)
    Hs = sp.hessian(C, MS)
    Hk = [[sp.diff(Hs[i, j], MS[k]) for j in range(14) for i in range(14)] for k in range(14)]

    def H_at(u):
        return [[sum(Fr(int(Hk[k][j * 14 + i])) * u[k] for k in range(14)) for j in range(14)]
                for i in range(14)]
    HK = H_at([0] * 13 + [1])
    W = [0, 1, 5, 6]
    Kp = [k for k in range(14) if k not in W and k != 13]

    def detpoly(u, deg=14):
        Hu = H_at(u)
        vals = []
        for h in range(deg + 1):
            vals.append(fr(qmat([[Hu[i][j] + h * HK[i][j] for j in range(14)]
                                 for i in range(14)]).det()))
        hh = sp.symbols("h")
        poly = sp.Poly(sp.interpolate(list(zip(range(deg + 1), [R(v) for v in vals])), hh), hh)
        return [fr(c) for c in reversed(poly.all_coeffs())] + [Fr(0)] * 20
    rng = random.Random(11)
    data = []
    for _ in range(6):
        u = [Fr(rng.randint(-9, 9), rng.randint(1, 4)) for _ in range(14)]
        co = detpoly(u)
        Hu = H_at(u)
        dK = fr(qmat([[Hu[i][j] for j in Kp] for i in Kp]).det())
        data.append((u, co, dK))
    check("B2 det(H(u) + h H(kappa)) has degree <= 3 in h at 6 random rational u",
          all(all(c == 0 for c in co[4:]) for _, co, _ in data))
    # p3 = dK * (alpha u0 u6 + beta u1 u5): solve alpha, beta from two points, test on four
    (u1, c1, k1), (u2, c2, k2) = data[0], data[1]
    al, be = sp.symbols("al be")
    sol = sp.solve([R(c1[3]) - R(k1) * (al * R(u1[0] * u1[6]) + be * R(u1[1] * u1[5])),
                    R(c2[3]) - R(k2) * (al * R(u2[0] * u2[6]) + be * R(u2[1] * u2[5]))],
                   [al, be])
    alpha, beta = fr(sol[al]), fr(sol[be])
    ok = all(co[3] == dK * (alpha * u[0] * u[6] + beta * u[1] * u[5]) for u, co, dK in data)
    check("B3 p3(u) = alpha * q(u) * det H(u)[K',K'] with alpha = -beta", ok and alpha == -beta
          and alpha != 0, "alpha=%s beta=%s" % (alpha, beta))

    # the printed example of Lemma C.11
    a = [Fr(0), Fr(1), Fr(3), Fr(7)]
    b = [L.QI(2, 1), L.QI(-1, 3), L.QI(5, -2), L.QI(4, 4)]
    m = L.matching_vector(a, b)
    y = [L.re_part(x) for x in m]
    yp = [-36, 108, 189, -569, 344, -126, 378, 288, 168, -784, 496, -512, 512, 16]
    check("B4 y = Re m(a,-b) equals the printed vector; q(y) = C(y) = 0",
          y == [Fr(x) for x in yp] and y[0] * y[6] - y[1] * y[5] == 0 and L.cubic_C(y) == 0)
    co = detpoly(y)
    Hy = H_at(y)
    rH = qrank(Hy)
    check("B5 rank H(y) = 12 and det(H(y) + h H(kappa)) = p2 h^2 with p2 = "
          "-114105877589169227887017984",
          rH == 12 and co[0] == 0 and co[1] == 0 and co[3] == 0 and all(c == 0 for c in co[4:])
          and co[2] == -114105877589169227887017984, "rank %d, p2 %s" % (rH, co[2]))
    # smoothness: Zariski tangent of Z = V(C, 13x13 minors of H) at y has dimension <= 12
    G = [sp.diff(C, MS[k]).subs(dict(zip(MS, y))) for k in range(14)]
    G = [fr(g) for g in G]
    funcs = [G]
    found = 0
    idx = list(range(14))
    for drop_r in range(14):
        for drop_c in range(14):
            I = [i for i in idx if i != drop_r]; J = [j for j in idx if j != drop_c]
            sub = [[Hy[i][j] for j in J] for i in I]
            if qrank(sub) != 12:
                continue
            # adjugate of the 13x13 block by cofactors
            adj = [[Fr(0)] * 13 for _ in range(13)]
            for r in range(13):
                for c in range(13):
                    mnr = [[sub[i][j] for j in range(13) if j != c] for i in range(13) if i != r]
                    adj[c][r] = (-1) ** (r + c) * fr(qmat(mnr).det())
            fun = []
            for k in range(14):
                Hk_sub = [[Fr(int(Hk[k][j * 14 + i])) for j in J] for i in I]
                fun.append(sum(adj[c][r] * Hk_sub[r][c] for r in range(13) for c in range(13)))
            funcs.append(fun)
            found += 1
            if found == 3:
                break
        if found == 3:
            break
    rf = qrank(funcs)
    check("B6 dC_y and the differentials of three 13x13 minors of H at y span a space of rank >= 2"
          " (so T_y Z has dimension <= 12)", rf >= 2, "rank %d" % rf)
    print("   section B time %.1f s" % (time.time() - t0))


# ---------------------------------------------------------------------------------------------
# C. Lemma C.10, all entries nonzero, without Lam
# ---------------------------------------------------------------------------------------------
def noncrossing_matchings(pts):
    if not pts:
        return [[]]
    out = []
    p0 = pts[0]
    for k in range(1, len(pts), 2):
        inner, outer = pts[1:k], pts[k + 1:]
        for a in noncrossing_matchings(inner):
            for b in noncrossing_matchings(outer):
                out.append([(p0, pts[k])] + a + b)
    return out


def section_C():
    print("\n=== C. Lemma C.10 (Vandermonde weights) ===")
    t0 = time.time()
    NC = noncrossing_matchings(list(range(8)))
    check("C1 there are 14 noncrossing perfect matchings of 8 points on a line", len(NC) == 14)
    subsets = list(itertools.combinations(range(8), 4))
    Mrows = []
    for S in subsets:
        Sset = set(S)
        Mrows.append([1 if all((i in Sset) != (j in Sset) for i, j in pi) else 0 for pi in NC])
    check("C2 compatibility matrix (70 x 14) has rank 14", qrank(Mrows) == 14)
    ctx = flint.fmpq_mpoly_ctx.get(("s", 7), "lex")
    s = ctx.gens()
    t = [ctx.from_dict({}) if False else 0 * s[0]]
    for k in range(7):
        t.append(t[-1] + s[k])

    def vand(S):
        p = 1 + 0 * s[0]
        for i, j in itertools.combinations(sorted(S), 2):
            p = p * (t[j] - t[i])
        return p
    rhs = [vand(S) * vand(tuple(k for k in range(8) if k not in S)) for S in subsets]
    chosen, rows = [], []
    for i, r in enumerate(Mrows):
        if qrank(rows + [r]) > len(rows):
            rows.append(r); chosen.append(i)
        if len(rows) == 14:
            break
    Minv = L.inv(rows)
    wts = []
    for j in range(14):
        p = 0 * s[0]
        for k in range(14):
            c = Minv[j][k]
            if c != 0:
                p = p + flint.fmpq(c.numerator, c.denominator) * rhs[chosen[k]]
        wts.append(p)
    ok = True
    for i, S in enumerate(subsets):
        lhs = 0 * s[0]
        for j in range(14):
            if Mrows[i][j]:
                lhs = lhs + wts[j]
        if lhs != rhs[i]:
            ok = False
    check("C3 the 14 weights solve all 70 equations Vand(t_S) Vand(t_S^c) = sum_{pi compat S} w_pi"
          " (polynomial identity in the gaps s_1..s_7)", ok)
    nonneg = all(all(c >= 0 for c in w.coeffs()) and w != 0 for w in wts)
    nterms = [len(w.coeffs()) for w in wts]
    check("C4 every weight is a nonzero polynomial in s_1..s_7 with nonnegative coefficients",
          nonneg, "terms per weight %s" % nterms)
    print("   section C time %.1f s" % (time.time() - t0))


# ---------------------------------------------------------------------------------------------
# D. (C.13) and (C.14)
# ---------------------------------------------------------------------------------------------
def dfun(O):
    M = sp.Matrix([[O[1, 0] * O[2, 0], O[1, 1] * O[2, 1], O[1, 2] * O[2, 2]],
                   [O[1, 0] * O[3, 0], O[1, 1] * O[3, 1], O[1, 2] * O[3, 2]],
                   [O[2, 0] * O[3, 0], O[2, 1] * O[3, 1], O[2, 2] * O[3, 2]]])
    return M.det()


def section_D():
    print("\n=== D. (C.13), (C.14) ===")
    t0 = time.time()
    a, b, c, d, e, f = sp.symbols("a b c d e f")
    q0, q1, q2, q3 = 1, a, b, c
    p0, p1, p2, p3 = 1, d, e, f
    Lq = sp.Matrix([[q0, -q1, -q2, -q3], [q1, q0, -q3, q2], [q2, q3, q0, -q1], [q3, -q2, q1, q0]])
    Rp = sp.Matrix([[p0, -p1, -p2, -p3], [p1, p0, p3, -p2], [p2, -p3, p0, p1], [p3, p2, -p1, p0]])
    T = Lq * Rp
    n = (1 + a ** 2 + b ** 2 + c ** 2) * (1 + d ** 2 + e ** 2 + f ** 2)
    check("D0 T T^t = n I and det T = n^2 (quaternion parametrisation of SO(4))",
          sp.expand(T * T.T - n * sp.eye(4)) == sp.zeros(4) and sp.expand(T.det() - n ** 2) == 0)
    # (C.13): on T00 = 0, d(T) = T01 T02 T03 T10 T20 T30 (d has degree 6, so O = T/sqrt(n) gives
    # the same identity for O)
    fsol = sp.solve(sp.Eq(T[0, 0], 0), f)
    check("D1 T00 is linear in f with a unique solution", len(fsol) == 1)
    Ts = T.subs(f, fsol[0])
    expr = sp.together(dfun(Ts) - Ts[0, 1] * Ts[0, 2] * Ts[0, 3] * Ts[1, 0] * Ts[2, 0] * Ts[3, 0])
    check("D1 (C.13): d(O) = O01 O02 O03 O10 O20 O30 on {O00 = 0} in SO(4) (rational identity)",
          sp.expand(sp.numer(expr)) == 0)
    # (C.14): O with first row (a1, b1, 0, 0) and rows (-b1 x_i, a1 x_i, y_i, z_i), [x y z] in SO(3)
    tau, u, v, wq = sp.symbols("tau u v wq")
    a1 = (1 - tau ** 2) / (1 + tau ** 2); b1 = 2 * tau / (1 + tau ** 2)
    qq = (1, u, v, wq)
    Nq = sum(x ** 2 for x in qq)
    w0, w1, w2, w3 = qq
    R3 = sp.Matrix([[w0**2 + w1**2 - w2**2 - w3**2, 2*(w1*w2 - w0*w3), 2*(w1*w3 + w0*w2)],
                    [2*(w1*w2 + w0*w3), w0**2 - w1**2 + w2**2 - w3**2, 2*(w2*w3 - w0*w1)],
                    [2*(w1*w3 - w0*w2), 2*(w2*w3 + w0*w1), w0**2 - w1**2 - w2**2 + w3**2]]) / Nq
    x = R3[:, 0]; y = R3[:, 1]; z = R3[:, 2]
    O = sp.Matrix([[a1, b1, 0, 0]] + [[-b1 * x[i], a1 * x[i], y[i], z[i]] for i in range(3)])
    check("D2 the block matrix is in SO(4)", sp.simplify(O * O.T - sp.eye(4)) == sp.zeros(4)
          and sp.simplify(O.det() - 1) == 0)
    Os = sp.Matrix(4, 4, sp.symbols("o0:16"))
    dO = dfun(Os)
    grad = [sp.diff(dO, Os[k]) for k in range(16)]
    sub = dict(zip(list(Os), list(O)))
    vals = []
    for i in (1, 2, 3):
        E = sp.zeros(4)
        E[0, i] = 1; E[i, 0] = -1
        D = E * O
        dd = sum(grad[k].subs(sub) * D[k // 4, k % 4] for k in range(16))
        j, k = [r for r in (1, 2, 3) if r != i]
        target = a1 * b1 * x[j - 1] * x[k - 1] * y[i - 1] * z[i - 1]
        vals.append(sp.cancel(sp.together(dd / target)))
    check("D3 (C.14): dd(E_i O) = s * a b x_j x_k y_i z_i with one common constant s for i = 1,2,3",
          len(set(vals)) == 1 and vals[0].is_number and vals[0] != 0, "s = %s" % vals[0])
    print("   section D time %.1f s" % (time.time() - t0))


if __name__ == "__main__":
    which = sys.argv[1:] or ["A", "B", "C", "D"]
    for sec in which:
        {"A": section_A, "B": section_B, "C": section_C, "D": section_D}[sec]()
    print("\nSUMMARY: %d failure(s)%s" % (len(FAILS), (": " + "; ".join(FAILS)) if FAILS else ""))
    sys.exit(1 if FAILS else 0)
