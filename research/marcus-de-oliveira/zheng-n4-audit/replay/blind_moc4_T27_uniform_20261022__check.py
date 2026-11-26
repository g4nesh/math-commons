"""Replay for blind_moc4_T27_uniform_20261022 (paper case T_28 = {2,3,5,6,9,12}). Auditor code."""
# ---------------------------------------------------------------------------
# Own code (auditor). Rebuilds everything from the paper's definitions via moc4lib:
# generators f_j (Appendix D), the cubic C of (A.2), kappa = e_13.
# Nothing from the supplement is executed or imported.
# ---------------------------------------------------------------------------
import sys, os, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L
import flint

FAIL = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + ((" : " + str(info)) if info != "" else ""))
    if not ok:
        FAIL.append(name)

def fvals_sum(d):
    s = [0] * 24
    for j, c in d.items():
        for t in range(24):
            s[t] += c * L.F_VALUES[j][t]
    return s

def identity(lhs, rhs):
    return fvals_sum({j: 1 for j in lhs}) == fvals_sum({j: 1 for j in rhs})

def f_action(g, js):
    out = {}
    for j in js:
        v = tuple(L.g_act_values(g, L.F_VALUES[j]))
        ks = [k for k in range(82) if tuple(L.F_VALUES[k]) == v]
        out[j] = ks[0] if ks else None
    return out

def stabilizer_has_action(T, claimed):
    """claimed: dict j->k on the generators; True if some g in G fixes T and acts so."""
    G = L.generators(T)
    for g in L.G_ELEMENTS:
        if L.g_act_set(g, T) != frozenset(T):
            continue
        m = f_action(g, G)
        if all(m[j] == k for j, k in claimed.items()):
            return g
    return None

def _split_z(P, ctx):
    out = {}
    for mon, c in zip(P.monoms(), P.coeffs()):
        e = mon[-1]; mm = list(mon); mm[-1] = 0
        out[e] = out.get(e, ctx.from_dict({})) + ctx.from_dict({tuple(mm): c})
    return out

def bareiss_det(M):
    M = [r[:] for r in M]
    n = len(M); sign = 1; prev = None
    for k in range(n - 1):
        piv = None
        for i in range(k, n):
            if not M[i][k].is_zero():
                if piv is None or len(M[i][k].coeffs()) < len(M[piv][k].coeffs()):
                    piv = i
        if piv is None:
            return M[0][0] * 0
        if piv != k:
            M[k], M[piv] = M[piv], M[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                t = M[i][j] * M[k][k] - M[i][k] * M[k][j]
                if prev is not None:
                    q, r = divmod(t, prev)
                    assert r.is_zero()
                    t = q
                M[i][j] = t
            M[i][k] = M[i][k] * 0
        prev = M[k][k]
    return sign * M[n - 1][n - 1]

def radial(weights):
    """weights: list of (f index, variable name). y = sum w_j f_j + z*kappa (matching
    coordinates). Returns ctx, vars, C, Q (with C(y0+z kappa) = C + z Q checked),
    D(z) = det Hess C (y0 + z kappa), N = deg_z D, R = Q^N D(-C/Q)."""
    names = [v for _, v in weights] + ["z"]
    ctx = flint.fmpz_mpoly_ctx.get(tuple(names), 'lex')
    V = dict(zip(names, ctx.gens()))
    y = [ctx.from_dict({}) for _ in range(14)]
    for (j, v) in weights:
        for k in range(14):
            fk = L.F_VECS[j][k]
            assert fk.denominator == 1
            if fk != 0:
                y[k] = y[k] + int(fk) * V[v]
    y[13] = y[13] + V["z"]                       # kappa = e_13 (paper (3.2))
    Cz = L.cubic_C(y)
    parts = _split_z(Cz, ctx)
    assert max(parts) <= 1, "C(y0 + z kappa) is not affine in z"
    Cc = parts.get(0, ctx.from_dict({})); Qc = parts.get(1, ctx.from_dict({}))
    mctx = flint.fmpz_mpoly_ctx.get(tuple("m%d" % i for i in range(14)), 'lex')
    Cm = L.cubic_C(list(mctx.gens()))
    H = []
    for k in range(14):
        dk = Cm.derivative(k); row = []
        for l in range(14):
            dkl = dk.derivative(l); e = ctx.from_dict({})
            for mon, c in zip(dkl.monoms(), dkl.coeffs()):
                idx = [i for i, ex in enumerate(mon) if ex][0]
                e = e + int(c) * y[idx]
            row.append(e)
        H.append(row)
    D = bareiss_det(H)
    Dk = _split_z(D, ctx)
    N = max(Dk)
    R = ctx.from_dict({})
    for k, dk in Dk.items():
        R = R + dk * (-Cc) ** k * Qc ** (N - k)
    return ctx, V, Cc, Qc, D, N, R

def nonneg_coeffs(P):
    return all(int(c) >= 0 for c in P.coeffs())

import sympy
T = {2, 3, 5, 6, 9, 12}
check("contact set equals paper T_28 exactly", frozenset(T) == L.CONTACT_SETS[28])
G = L.generators(T)
check("generators = {0,1,2,5,7,9,39,40,41,67,68}", G == [0, 1, 2, 5, 7, 9, 39, 40, 41, 67, 68], G)
check("non-minor generators are exactly 40, 68", [j for j in G if not L.F_IS_MINOR[j]] == [40, 68])
check("identity f40 + f68 = f0 + f7 + f41 + f67", identity([40, 68], [0, 7, 41, 67]))
check("identity f40 + f1 = f0 + f7 + f39", identity([40, 1], [0, 7, 39]))
check("identity f40 + f9 = f2 + f7 + f41", identity([40, 9], [2, 7, 41]))
g = stabilizer_has_action(T, {40: 68, 68: 40})
check("a group element fixes T_28 and swaps f40, f68", g is not None, g)
def find_ordered(src, dst):
    for gg in L.G_ELEMENTS:
        if all(tuple(L.g_act_values(gg, L.F_VALUES[s])) == tuple(L.F_VALUES[t]) for s, t in zip(src, dst)):
            return gg
    return None
check("group element maps (f59,f53,f55,f67,f54) to (f5,f2,f0,f67,f40) in order", find_ordered([59, 53, 55, 67, 54], [5, 2, 0, 67, 40]) is not None)
check("group element maps (f0,f6,f10,f54) to (f0,f67,f2,f40) in order", find_ordered([0, 6, 10, 54], [0, 67, 2, 40]) is not None)
check("group element maps (f0,f6,f10,f54) to (f67,f41,f39,f40) in order", find_ordered([0, 6, 10, 54], [67, 41, 39, 40]) is not None)

t0 = time.time()
W = [(7, 'q'), (39, 'p'), (5, 'j'), (41, 'c'), (2, 'e'), (40, 'd'), (0, 'b'), (67, 'a')]
ctx, V, Cc, Qc, D, N, R = radial(W)
q, p, j, c, e, d, b, a = [V[k] for k in 'qpjcedba']
print("INFO Hessian determinant: %d terms, degree %d in z, %.0f s" % (len(D.coeffs()), N, time.time() - t0))
check("Q = d^2 + de + dp + ep - ej", Qc == d ** 2 + d * e + d * p + e * p - e * j)
check("C has 25 terms, all coefficients positive, contains d^3", len(Cc.coeffs()) == 25 and all(int(x) > 0 for x in Cc.coeffs()) and int(Cc.to_dict().get((0, 0, 0, 0, 0, 3, 0, 0, 0), 0)) == 1)
LHS = R * Qc ** (4 - N)
quo, rem = divmod(LHS, -4 * Qc ** 2)
check("Q^4 det H(y0 - (C/Q) kappa) divisible by -4Q^2", rem.is_zero())
G9 = quo.sqrt()
check("quotient is a perfect square G9^2; G9 has 2124 terms, homogeneous of degree 9", G9 is not None and len(G9.coeffs()) == 2124 and len(set(sum(m) for m in G9.monoms())) == 1 and G9.total_degree() == 9)
# dehomogenize d = 1
n2 = ('q', 'p', 'j', 'c', 'e', 'b', 'a'); c2 = flint.fmpz_mpoly_ctx.get(n2, 'lex')
def dehom(P):
    out = {}
    for mon, cf in zip(P.monoms(), P.coeffs()):
        assert mon[8] == 0
        m = (mon[0], mon[1], mon[2], mon[3], mon[4], mon[6], mon[7])
        out[m] = out.get(m, 0) + int(cf)
    return c2.from_dict(out)
Gd = dehom(G9); Cr = dehom(Cc)
check("deg_j G9 = 2", Gd.degrees()[2] == 2)
ks = {}
for mon, cf in zip(Gd.monoms(), Gd.coeffs()):
    m = list(mon); ee = m[2]; m[2] = 0
    ks[ee] = ks.get(ee, c2.from_dict({})) + c2.from_dict({tuple(m): cf})
k0, k1, k2 = ks[0], ks[1], ks[2]
G6, r = divmod(k0, Cr)
check("k0 = C * G6 with G6 of 196 terms, bilinear in e, p", r.is_zero() and len(G6.coeffs()) == 196 and G6.degrees()[4] == 1 and G6.degrees()[1] == 1, len(G6.coeffs()))
parts = {}
for mon, cf in zip(G6.monoms(), G6.coeffs()):
    m = list(mon); key = (m[4], m[1]); m[4] = 0; m[1] = 0
    parts[key] = parts.get(key, c2.from_dict({})) + c2.from_dict({tuple(m): cf})
g00, gE, gP, gEP = parts[(0, 0)], parts[(1, 0)], parts[(0, 1)], parts[(1, 1)]
check("g00 has only positive coefficients and constant term 1", all(int(x) > 0 for x in g00.coeffs()) and [cf for mon, cf in zip(g00.monoms(), g00.coeffs()) if sum(mon) == 0] == [1])
k2e = k2.subs({'e': 0}); k1e = k1.subs({'e': 0})
check("k2|e=0 = 0 and k1|e=0 has 375 nonnegative coefficients", k2e.is_zero() and len(k1e.coeffs()) == 375 and nonneg_coeffs(k1e))
# (8) on four chambers
c3 = flint.fmpz_mpoly_ctx.get(('A', 'T', 'U', 'q'), 'lex'); A3, T3, U3, q3 = c3.gens(); one3 = c3.from_dict({(0, 0, 0, 0): 1})
md = [max(P.degrees()[i] for P in (g00, gE, gP, gEP)) for i in range(7)]
def substU(P, a_, b_, c_):
    da, db, dc = md[6], md[5], md[3]; out = c3.from_dict({})
    for mon, cf in zip(P.monoms(), P.coeffs()):
        eq, ep, ej, ec, ee, eb, ea = mon
        out += int(cf) * q3 ** eq * a_[0] ** ea * a_[1] ** (da - ea) * b_[0] ** eb * b_[1] ** (db - eb) * c_[0] ** ec * c_[1] ** (dc - ec)
    return out
for chb in ("ab<1", "ab>1"):
    for chc in ("ac<1", "ac>1"):
        a_ = (A3, one3)
        b_, Le = ((T3, A3 * (1 + T3)), (c3.from_dict({}), one3)) if chb == "ab<1" else ((1 + T3, A3), (T3 * A3, (A3 + 1) * (1 + T3 + A3)))
        c_, Lp = ((U3, A3 * (1 + U3)), (c3.from_dict({}), one3)) if chc == "ac<1" else ((1 + U3, A3), (U3 * A3, (A3 + 1) * (1 + U3 + A3)))
        s00, sE, sP, sEP = [substU(P_, a_, b_, c_) for P_ in (g00, gE, gP, gEP)]
        fs = [s00, sE * Le[1] + Le[0] * s00, sP * Lp[1] + Lp[0] * s00,
              sEP * Le[1] * Lp[1] + Lp[0] * Le[1] * sE + Le[0] * Lp[1] * sP + Le[0] * Lp[0] * s00]
        check("(8) four functions nonnegative-coefficient on chamber %s, %s" % (chb, chc), all(nonneg_coeffs(f) for f in fs), [len(f.coeffs()) for f in fs])
# Bernstein charts (9)
Vn = ['A', 'T', 'W', 'E', 'P', 'Q']; cx = flint.fmpz_mpoly_ctx.get(tuple(Vn), 'lex'); A, Tq, Wv, E, P, Q = cx.gens()
one = cx.from_dict({(0,) * 6: 1})
idx = {'q': 0, 'p': 1, 'j': 2, 'c': 3, 'e': 4, 'b': 5, 'a': 6}; hv = ['p', 'c', 'e', 'b', 'a']
K3 = [k0, k1, k2]
Nv = {v: max(k.degrees()[idx[v]] for k in K3) for v in hv}
ch = flint.fmpz_mpoly_ctx.get(n2 + tuple('h' + v for v in hv), 'lex')
def homog(Pp):
    out = {}
    for mon, cf in zip(Pp.monoms(), Pp.coeffs()):
        out[tuple(list(mon) + [Nv[v] - mon[idx[v]] for v in hv])] = int(cf)
    return ch.from_dict(out)
KH = [homog(k) for k in K3]
def chart(sub, J):
    args = [Q, sub['p'][0], cx.from_dict({}), sub['c'][0], sub['e'][0], sub['b'][0], sub['a'][0]] + [sub[v][1] for v in hv]
    S = [kh.compose(*args) for kh in KH]
    nJ, dJ = J
    return S[0], 2 * S[0] * dJ + nJ * S[1], S[0] * dJ ** 2 + nJ * dJ * S[1] + nJ ** 2 * S[2]
def sturm_groups(dd):
    groups = {}
    for mon, cf in zip(dd.monoms(), dd.coeffs()):
        groups.setdefault(tuple(m for i, m in enumerate(mon) if i != 3), {})[mon[3]] = int(cf)
    x = sympy.Symbol('x'); ng = 0; bad = 0
    for key, cs in groups.items():
        if all(v >= 0 for v in cs.values()):
            continue
        ng += 1; lo = min(cs)
        pol = sympy.Poly(sum(v * x ** (k - lo) for k, v in cs.items()), x)
        if not (pol.eval(0) > 0 and pol.count_roots(0, sympy.oo) == 0):
            bad += 1
    return ng, bad
expected = {"1l": (23160, 30838, 30838), "1g": (62476, 82332, 82332), "2l": (8928, 18266, 25267), "2g": (12719, 25703, 35113),
            "3l": (9857, 20025, 27337), "3g": (33660, 58371, 72870), "7": (4548, 7772, 10983), "8": (30504, 63441, 100982)}
for name in ["1l", "1g", "2l", "2g", "3l", "3g", "7", "8"]:
    kind = name[0]
    if kind == '1':
        a_ = (Tq ** 2 + A * (1 + Tq) ** 2, (1 + A) * (1 + Tq) ** 2); b_ = (Tq ** 2 * (1 + A), Tq ** 2 + A * (1 + Tq) ** 2); e_ = (E, one); J = (E + 1, E)
    elif kind == '2':
        a_ = (1 + A, one); b_ = (Tq ** 2, (1 + A) * (1 + Tq) ** 2); e_ = (E, one); J = ((1 + A + A * Tq) ** 2 + E * (1 + A), E * (1 + A + A * Tq) ** 2)
    elif kind == '3':
        b_ = (1 + A, one); a_ = (Tq ** 2, (1 + A) * (1 + Tq) ** 2); e_ = (E, one); J = ((1 + A + A * Tq) ** 2 + E * (1 + A), E * (1 + A + A * Tq) ** 2)
    if kind in '123':
        if name[1] == 'l':
            c_ = (Wv * a_[1], a_[0] * (1 + Wv)); p_ = (P, one)
        else:
            c_ = ((1 + Wv) * a_[1], a_[0])
            p_ = ((a_[0] + a_[1]) * (c_[0] + c_[1]) * P, a_[1] * c_[1] * Wv * (1 + P))
    elif kind == '7':
        a_ = (A, one); b_ = (1 + Tq, A); Eab = ((A + 1) * (1 + Tq + A), A * Tq)
        e_ = (Eab[0] * E, Eab[1] * (1 + E)); J = (Eab[1], Eab[0] * E); c_ = (Wv, A * (1 + Wv)); p_ = (P, one)
    else:
        U = P; D_ = 1 + U; Dl = (1 + Tq) * D_ - 1; Hh = (1 + Tq) * D_ + A
        a_ = (A, one); b_ = (1 + Tq, A); c_ = (1 + Wv, A)
        p_ = (U * (A + 1) * ((1 + Wv) * D_ + A * U), A * D_ * ((1 + Wv) * D_ - U))
        e_ = ((A + 1) * Hh * E, A * D_ * Dl * (1 + E)); J = (A * Dl, D_ * (A + 1) * Hh * E)
    ds = chart({'a': a_, 'b': b_, 'c': c_, 'e': e_, 'p': p_}, J)
    dens_ok = all(nonneg_coeffs(x[1]) for x in (a_, b_, c_, e_, p_)) and nonneg_coeffs(J[1])
    counts = tuple(len(x.coeffs()) for x in ds); negs = [sum(1 for y in x.coeffs() if int(y) < 0) for x in ds]
    check("chart %s: term counts %s match the document" % (name, counts), counts == expected[name])
    ok = dens_ok
    for i, dd in enumerate(ds):
        if negs[i]:
            ng, bad = sturm_groups(dd)
            print("INFO chart %s datum %d: %d negative coefficients, %d E-groups with a negative coefficient, %d groups failing" % (name, i, negs[i], ng, bad))
            ok = ok and bad == 0
    check("chart %s: Bernstein data nonnegative (coefficients, or exact real-root count in E for mixed groups)" % name, ok, negs)
print("RESULT: %d exact checks failed" % len(FAIL) + ((": " + ", ".join(FAIL)) if FAIL else ""))
sys.exit(1 if FAIL else 0)
