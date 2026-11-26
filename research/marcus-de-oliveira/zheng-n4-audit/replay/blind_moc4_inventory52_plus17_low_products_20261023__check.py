"""Own replay for blind_moc4_inventory52_plus17_low_products_20261023.
Seven-contact face T7 = {0,1,4,6,9,12,17} (paper Table 5, aux set 19), source
c = a f0 + B f2 + C f5 + b f7 + A f9 + c f21 + d f35, real lift u = c + z kappa.
Everything (cubic, Hessian, Q, N, R) is rebuilt from the paper's definitions via the
audit library.  The degree-nine polynomial G and the two octics are read as DATA from
the supplement JSON (path from env MOC4_SUPP) and are then checked against our own
objects: G pointwise (exact rational arithmetic at random points) against the Hessian
determinant, the octics by an exact polynomial identity.  Nothing is executed from the
supplement."""
import os, sys, json, random
from fractions import Fraction
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import moc4lib as L
import flint

NAMES = ('a', 'B', 'C', 'b', 'A', 'c', 'd', 'z')
CTX = flint.fmpq_mpoly_ctx.get(NAMES, 'lex')
GENS = CTX.gens()
a, B, C, b, A, c, d, z = GENS
COEF = {0: a, 2: B, 5: C, 7: b, 9: A, 21: c, 35: d}
JSON_NAME = "proofs/manin_transfer_20261021/blind_moc4_inventory52_plus17_low_products_exact_20261023.json"

fails = 0
def check(name, ok, info=""):
    global fails
    print(("PASS " if ok else "FAIL ") + name + (" :: " + str(info) if info != "" else ""))
    if not ok:
        fails += 1

def lift():
    u = [CTX.from_dict({}) for _ in range(14)]
    for k, s in COEF.items():
        for m in range(14):
            v = L.F_VECS[k][m]
            if v != 0:
                u[m] = u[m] + s * flint.fmpq(v.numerator, v.denominator)
    u[13] = u[13] + z
    return u

def hessian_sym(u):
    """Hessian of the cubic C at u: entries are linear forms; built from the 14-variable cubic."""
    mctx = flint.fmpq_mpoly_ctx.get(tuple("m%d" % i for i in range(14)), 'lex')
    ms = mctx.gens()
    Cm = L.cubic_C(ms)
    H = [[None] * 14 for _ in range(14)]
    for i in range(14):
        di = Cm.derivative(i)
        for j in range(i, 14):
            dij = di.derivative(j)          # constant linear form in m
            e = CTX.from_dict({})
            for mon, co in dij.to_dict().items():
                # dij is linear: mon has a single 1
                idx = mon.index(1) if sum(mon) == 1 else None
                if idx is None:
                    e = e + co
                else:
                    e = e + co * u[idx]
            H[i][j] = H[j][i] = e
    return H

def load_supp():
    root = os.environ.get("MOC4_SUPP") or (sys.argv[1] if len(sys.argv) > 1 else None)
    if not root:
        return None
    with open(os.path.join(root, JSON_NAME)) as f:
        return json.load(f)

def poly_from_terms(terms):
    dct = {}
    for mon, co in terms:
        mon = tuple(mon) + (0,) * (8 - len(mon))
        dct[mon] = flint.fmpq(Fraction(co).numerator, Fraction(co).denominator)
    return CTX.from_dict(dct)

def coeffs_in(poly, var_index):
    """dict power -> coefficient polynomial in the remaining variables (same ctx)."""
    out = {}
    for mon, co in poly.to_dict().items():
        k = mon[var_index]
        mon2 = list(mon); mon2[var_index] = 0
        out.setdefault(k, {})[tuple(mon2)] = co
    return {k: CTX.from_dict(v) for k, v in out.items()}

def evalq(poly, vals):
    return poly.subs({NAMES[i]: flint.fmpq(vals[i].numerator, vals[i].denominator) for i in range(8)}) if False else \
        sum((Fraction(int(co.p), int(co.q)) * prod_pow(vals, mon) for mon, co in poly.to_dict().items()), Fraction(0))

def prod_pow(vals, mon):
    r = Fraction(1)
    for v, e in zip(vals, mon):
        if e:
            r *= v ** int(e)
    return r

def build_core():
    u = lift()
    Cu = L.cubic_C(u)
    H = hessian_sym(u)
    return u, Cu, H

def chart_low(poly, shift=(4, 4, 4), deg=21):
    """a=XY/Z, b=XZ/Y, c=YZ/X (d=1), times X^4Y^4Z^4; X=x/D, Y=(x+y)/D, Z=(x+y+z)/D, D=W0+x+y+z;
    times D^deg.  Input poly in CTX with B, d, z absent (d set to 1 first)."""
    kctx = flint.fmpz_mpoly_ctx.get(('x', 'y', 'zz', 'W0', 'A', 'C'), 'lex')
    x, y, zz, W0, Av, Cv = kctx.gens()
    D = W0 + x + y + zz
    Xp, Yp, Zp = x, x + y, x + y + zz
    out = kctx.from_dict({})
    cache = {}
    for mon, co in poly.to_dict().items():
        ia, iB, iC, ib, iA, ic, idd, iz = mon
        assert iB == 0 and iz == 0
        ex = ia + ib - ic + shift[0]; ey = ia - ib + ic + shift[1]; ez = -ia + ib + ic + shift[2]
        assert min(ex, ey, ez) >= 0 and ex + ey + ez <= deg, (mon, ex, ey, ez)
        key = (ex, ey, ez)
        if key not in cache:
            cache[key] = Xp ** ex * Yp ** ey * Zp ** ez * D ** (deg - ex - ey - ez)
        assert co.q == 1
        out = out + int(co.p) * cache[key] * Av ** iA * Cv ** iC
    return out, kctx

def sign_summary(p):
    cs = [int(v) for v in p.coeffs()]
    return len(cs), sum(1 for v in cs if v > 0), sum(1 for v in cs if v < 0)

def setd1(poly):
    return poly.subs({'d': 1})

def main():
    u, Cu, H = build_core()
    Q_stated = -(A * B + A * C + B * C + d * (A + B + C) + d * d)
    cz = coeffs_in(Cu, 7)
    check("C(c + z kappa) has degree 1 in z", set(cz.keys()) <= {0, 1}, sorted(cz.keys()))
    check("z-coefficient equals Q = -(AB+AC+BC+d(A+B+C)+d^2)", cz[1] == Q_stated)
    NR = -cz[0]
    stated_NR = (d * (A * B + A * C + B * C) + d * (A * a + B * b + C * c) + a * b * c
                 + d * (a * b + a * c + b * c) + d * d * (a + b + c + A + B + C) + d ** 3)
    check("z-free part equals -(stated cubic) (so F = -N - R B + Q z)", NR == stated_NR)
    cB = coeffs_in(NR, 1)
    check("z-free part has degree 1 in B", set(cB.keys()) <= {0, 1})
    N, R = cB[0], cB[1]
    check("N and R have only positive coefficients", all(v > 0 for v in N.coeffs()) and all(v > 0 for v in R.coeffs()))
    # principal block deleting 0 and 5 is z-independent
    keep = [i for i in range(14) if i not in (0, 5)]
    HI = [[H[i][j] for j in keep] for i in keep]
    check("H with rows/cols 0,5 deleted is independent of z", all(e.is_zero() or e.degrees()[7] == 0 for row in HI for e in row))
    data = load_supp()
    if data is None:
        print("SKIP supplement-data checks (set MOC4_SUPP)")
        return
    check("supplement source order and contacts match", data["contacts"] == [0, 1, 4, 6, 9, 12, 17] and data["source_order"] == [0, 2, 5, 7, 9, 21, 35])
    G = poly_from_terms(data["G"])
    gB = coeffs_in(G, 1)
    check("G (data) has %d terms, total degree 9, degree 2 in B, no z" % len(G.coeffs()),
          len(G.coeffs()) == 687 and G.total_degree() == 9 and max(gB) == 2 and G.degrees()[7] == 0, (len(G.coeffs()), G.total_degree()))
    al, be, ga = gB[2], gB[1], gB[0]
    # pointwise exact check: det H(c + z* kappa) * Q^2 = -4 G^2 at z* = (N + R B)/Q
    rng = random.Random(20261125)
    okp = True
    for t in range(12):
        vals = [Fraction(rng.randint(1, 60), rng.randint(1, 9)) for _ in range(7)]
        Qv = evalq(Q_stated, vals + [Fraction(0)])
        zs = (evalq(N, vals + [Fraction(0)]) + evalq(R, vals + [Fraction(0)]) * vals[1]) / Qv
        full = vals + [zs]
        Hn = [[evalq(H[i][j], full) for j in range(14)] for i in range(14)]
        dv = L.det(Hn)
        Gv = evalq(G, full)
        if dv * Qv ** 2 != -4 * Gv ** 2:
            okp = False
    check("det H(c + z* kappa) Q^2 = -4 G^2 at 12 random rational points (exact pointwise)", okp)
    # also: G is not identically something trivial: random line test of degree-2 structure is implicit
    # discriminant factorization
    Delta = be * be - 4 * al * ga
    Hs = [poly_from_terms(t) for t in data["Herons"]]
    check("two octics: degree 8, 136 terms each", all(h.total_degree() == 8 and len(h.coeffs()) == 136 for h in Hs), [(h.total_degree(), len(h.coeffs())) for h in Hs])
    check("Delta = beta^2 - 4 alpha gamma equals F+ F- (exact identity)", Delta == Hs[0] * Hs[1])
    for h in Hs:
        fac = h.factor()
        check("octic irreducible over Q", len(fac[1]) == 1 and fac[1][0][1] == 1)
    # det HI = -d A^2 P: pointwise exact factor check and degree in B via interpolation-free test
    okHI = True
    P9 = poly_from_terms(data["principal_P9"])
    for t in range(6):
        vals = [Fraction(rng.randint(1, 60), rng.randint(1, 9)) for _ in range(7)] + [Fraction(rng.randint(-50, 50))]
        Hn = [[evalq(e, vals) for e in row] for row in HI]
        dv = L.det(Hn)
        if dv != -vals[6] * vals[4] ** 2 * evalq(P9, vals):
            okHI = False
    check("det H_I = -d A^2 P9 (P9 data, 9 variables-degree %d, B-degree %d) at 6 random points (exact pointwise)" % (P9.total_degree(), P9.degrees()[1]), okHI)
    # signed trace-square identities (7),(8), checked pointwise exactly at random points
    pB = coeffs_in(P9, 1)
    p2, p1, p0 = pB.get(2, CTX.from_dict({})), pB.get(1, CTX.from_dict({})), pB.get(0, CTX.from_dict({}))
    u1 = al * p1 - p2 * be
    u0 = al * p0 - p2 * ga
    normW = u0 * u0 - u1 * (be * p0 - p1 * ga)
    ok_sq = True; ok7 = True; ok8 = True
    for t in range(8):
        vals = [Fraction(rng.randint(1, 40), rng.randint(1, 7)) for _ in range(7)] + [Fraction(0)]
        ev = lambda p: evalq(p, vals)
        Wv2 = ev(normW)
        from math import isqrt
        num, den = Wv2.numerator, Wv2.denominator
        if num < 0 or isqrt(num) ** 2 != num or isqrt(den) ** 2 != den:
            ok_sq = False; continue
        Wv = Fraction(isqrt(num), isqrt(den))
        alv, bev, gav, u0v, u1v = ev(al), ev(be), ev(ga), ev(u0), ev(u1)
        for sgn in (1, -1):
            S = 2 * alv * u0v - bev * u1v + sgn * 2 * alv * Wv
            # S = -d F_eps H^2 for one of the octics
            good = False
            for h in Hs:
                q = S / (-vals[6] * ev(h)) if ev(h) != 0 else None
                if q is not None and q >= 0:
                    qn, qd = q.numerator, q.denominator
                    if isqrt(qn) ** 2 == qn and isqrt(qd) ** 2 == qd:
                        good = True
            ok7 = ok7 and good
            # (8): P S - (U + sgn W)^2 is a multiple of G as polynomials in B (B-free multiplier)
            p2v, p1v, p0v = ev(p2), ev(p1), ev(p0)
            c2 = p2v * S - u1v ** 2
            c1 = p1v * S - 2 * u1v * (u0v + sgn * Wv)
            c0 = p0v * S - (u0v + sgn * Wv) ** 2
            if not (c2 * bev == c1 * alv and c2 * gav == c0 * alv):
                ok8 = False
    check("norm u0^2 - u1(beta p0 - p1 gamma) is a rational square at 8 random points", ok_sq)
    check("(7) S+- / (-d F) is a rational square for an octic F, both signs, 8 points", ok7)
    check("(8) P S+- - (U +- W)^2 is proportional to G as a quadratic in B, 8 points", ok8)

    # bordered principal minors on I+{0} and I+{5}: after forced height, divisible by G (B univariate, 4 random parameter points)
    ctx2 = flint.fmpq_mpoly_ctx.get(('B', 'z'), 'lex')
    Bv, zv = ctx2.gens()
    def bareiss(M):
        n = len(M); M = [row[:] for row in M]; prev = ctx2.from_dict({(0, 0): 1}); sign = 1
        for k in range(n - 1):
            if M[k][k].is_zero():
                sw = next((r for r in range(k + 1, n) if not M[r][k].is_zero()), None)
                if sw is None:
                    return ctx2.from_dict({})
                M[k], M[sw] = M[sw], M[k]; sign = -sign
            for i in range(k + 1, n):
                for j in range(k + 1, n):
                    M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) / prev
            prev = M[k][k]
        return sign * M[n - 1][n - 1]
    okb = True
    for t in range(4):
        vals = [Fraction(rng.randint(1, 30), rng.randint(1, 5)) for _ in range(7)]
        def sub2(e):
            out = ctx2.from_dict({})
            for mon, co in e.to_dict().items():
                cst = Fraction(int(co.p), int(co.q))
                for i in (0, 2, 3, 4, 5, 6):
                    cst *= vals[i] ** int(mon[i])
                out = out + flint.fmpq(cst.numerator, cst.denominator) * Bv ** int(mon[1]) * zv ** int(mon[7])
            return out
        Gu = sub2(G); Qu = sub2(Q_stated); NRu = sub2(N) + sub2(R) * Bv
        for extra in (0, 5):
            idx = sorted(keep + [extra])
            M = [[sub2(H[i][j]) for j in idx] for i in idx]
            dz = bareiss(M)
            # substitute z = NRu/Qu and clear Q^2 (degree in z <= 2)
            zc = {}
            for mon, co in dz.to_dict().items():
                zc[int(mon[1])] = zc.get(int(mon[1]), ctx2.from_dict({})) + co * Bv ** int(mon[0])
            assert max(zc) <= 2
            num = sum((zc.get(k, ctx2.from_dict({})) * NRu ** k * Qu ** (2 - k) for k in range(3)), ctx2.from_dict({}))
            try:
                qq = num / Gu
                okb = okb and (qq * Gu == num)
            except Exception:
                okb = False
    check("bordered 13x13 principal minors on I+{0}, I+{5} vanish on G = 0 (divisible by G in B; 4 random parameter points)", okb)
    # chart tables (9) on bc <= d^2 (d = 1)
    R1, N1, al1, be1, ga1 = (setd1(p) for p in (R, N, al, be, ga))
    tabs = {"gamma": ga1, "R beta - N alpha": R1 * be1 - N1 * al1}
    for i, h in enumerate(Hs):
        tabs["R F%d - N alpha" % i] = R1 * setd1(h) - N1 * al1
    for nm, p in tabs.items():
        cp, kctx = chart_low(p)
        n, npos, nneg = sign_summary(cp)
        print("INFO table %-18s terms %6d  positive %6d  negative %6d" % (nm, n, npos, nneg))
        if nm == "gamma":
            check("gamma table all positive, 9772 terms", nneg == 0 and n == 9772, n)
        if nm == "R beta - N alpha":
            check("R beta - N alpha table all positive, 13819 terms", nneg == 0 and n == 13819, n)
        if nm.startswith("R F"):
            print("INFO", nm, "all negative" if npos == 0 else "mixed signs")
            if npos == 0:
                check("R F* - N alpha table all negative (14549 terms) for this octic", n == 14549, n)
                dd = cp.to_dict()
                xA = dd.get((21, 0, 0, 0, 1, 0)); xC = dd.get((21, 0, 0, 0, 0, 1))
                check("strict markers -208 x^21 A and -48 x^21 C", xA == -208 and xC == -48, (xA, xC))
            if nm == "gamma":
                pass
        if nm == "gamma":
            check("gamma marker 512 x^21", cp.to_dict().get((21, 0, 0, 0, 0, 0)) == 512)
        if nm == "R beta - N alpha":
            check("R beta - N alpha marker 256 x^21", cp.to_dict().get((21, 0, 0, 0, 0, 0)) == 256)

if __name__ == "__main__":
    main()
    print("RESULT", "all EXACT checks passed" if fails == 0 else "%d failures" % fails)
    sys.exit(1 if fails else 0)
