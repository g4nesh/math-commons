"""Exact replay and numeric tests for blind_moc4_T48_allocation_I_low_j_uniform_20261023 (audit label C26-T49).
Large-lambda, low-j both-chain cell: three charts.
Run: python -I -B blind_moc4_T48_allocation_I_low_j_uniform_20261023__check.py   [optional env MOC4_T49_G_CACHE=<json>, flags --write-cache / --fresh]
Prints PASS/FAIL per EXACT check and NUMERIC lines; exits nonzero if an EXACT check fails."""
# ---------------------------------------------------------------------------------
# Shared core (identical in every C26-T49 check script).  Own code only; nothing from
# the supplement is executed or imported.  The audit library moc4lib is imported from
# ../lib relative to this file.
#
# G is the radial factor of the T_49 (supplement "T48") eight-generator face, rebuilt
# from the paper's definitions: c = sum coef*f_role, C(c + h kappa) = C0 + h Q0 (C the
# cubic (A.2)), D(h) = det Hess C(c + h kappa) by exact Bareiss elimination, and
# -Q0^3 D(-C0/Q0) = 4 G^2 Q0.  Computing G takes about 4 minutes.  If the environment
# variable MOC4_T49_G_CACHE names a JSON file written by an earlier run (--write-cache),
# it is loaded instead, and accepted only if its canonical SHA256 equals G_SHA256, the
# hash printed by the from-scratch computation.
# ---------------------------------------------------------------------------------
import sys, os, json, time, hashlib, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L
import flint
import sympy as sp

FAILS = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + (("  " + str(info)) if info != "" else ""), flush=True)
    if not ok:
        FAILS.append(name)

GV = ["q", "p", "j", "c", "e", "d", "b", "a"]
ROLES = dict(q=7, p=5, j=29, c=39, e=2, d=40, b=41, a=0)
T49 = frozenset({0, 3, 5, 6, 9, 12})
G_SHA256 = "c89c8cbc5331fcfdbdc22587fbcf8a017fc692eb8fd1f38edbe8a59ca0550a49"


def canon_hash_terms(vars_, terms):
    items = sorted(([int(x) for x in m], int(c)) for m, c in terms)
    return hashlib.sha256(json.dumps({"vars": list(vars_), "terms": items}, separators=(",", ":")).encode()).hexdigest()


def compute_G():
    names = GV + ["h"]
    ctx = flint.fmpz_mpoly_ctx.get(tuple(names), 'degrevlex')
    V = dict(zip(names, ctx.gens())); zero = ctx.from_dict({})
    m = []
    for k in range(14):
        s = zero
        for n, r in ROLES.items():
            x = L.F_VECS[r][k]
            assert x.denominator == 1
            if x != 0:
                s = s + V[n] * int(x)
        m.append(s)
    m[13] = m[13] + V["h"]          # kappa = e_13
    ms, Cexpr = L.cubic_sympy()
    def to_flint(expr):
        s = zero
        for mon, co in sp.Poly(expr, *ms).terms():
            t = zero + int(co)
            for k, e in enumerate(mon):
                if e:
                    t = t * m[k] ** e
            s = s + t
        return s
    H = sp.hessian(Cexpr, ms)
    A = [[to_flint(H[i, k]) for k in range(14)] for i in range(14)]
    Cf = to_flint(Cexpr)
    parts = {}
    for mon, co in Cf.to_dict().items():
        parts.setdefault(int(mon[-1]), {})[tuple(list(mon[:-1]) + [0])] = co
    check("C(c + h kappa) is affine in h", set(parts) <= {0, 1})
    C0 = ctx.from_dict(parts[0]); Q0 = ctx.from_dict(parts[1])
    n = 14; prev = zero + 1; sign = 1
    for k in range(n - 1):
        if A[k][k].is_zero():
            sw = next(i for i in range(k + 1, n) if not A[i][k].is_zero())
            A[k], A[sw] = A[sw], A[k]; sign = -sign
        for i in range(k + 1, n):
            for jx in range(k + 1, n):
                num = A[k][k] * A[i][jx] - A[i][k] * A[k][jx]
                if prev.is_one():
                    A[i][jx] = num
                else:
                    qq, rr = divmod(num, prev)
                    assert rr.is_zero()
                    A[i][jx] = qq
            A[i][k] = zero
        prev = A[k][k]
    D = A[n - 1][n - 1] * sign
    byh = {}
    for mon, co in D.to_dict().items():
        byh.setdefault(int(mon[-1]), {})[tuple(list(mon[:-1]) + [0])] = co
    kmax = max(byh)
    Nn = zero
    for k, dct in byh.items():
        Nn = Nn + ctx.from_dict(dct) * (-C0) ** k * Q0 ** (kmax - k)
    check("det Hess C(c + h kappa) has degree 3 in h", kmax == 3)
    qq, rr = divmod(-Nn, Q0 ** (kmax - 2))
    check("-Q0^3 det Hess C(c - (C0/Q0) kappa) divisible by Q0", rr.is_zero())
    G2, rr = divmod(qq, zero + 4)
    check("... and by 4", rr.is_zero())
    G = G2.sqrt()
    check("quotient is 4 G^2 for a polynomial G (det Hess = -4 G^2 / Q0^2)", G * G == G2)
    terms = [([int(x) for x in mon[:-1]], int(co)) for mon, co in G.to_dict().items()]
    return terms


def load_or_compute_G():
    t0 = time.time()
    cache = os.environ.get("MOC4_T49_G_CACHE")
    if cache and os.path.exists(cache) and "--fresh" not in sys.argv:
        d = json.load(open(cache))
        terms = [(m, c) for m, c in d["terms"]]
        check("cached G has the hash of the from-scratch computation", canon_hash_terms(GV, terms) == G_SHA256)
        src = "cache"
    else:
        terms = compute_G()
        src = "computed from the definitions in %.0f s" % (time.time() - t0)
        if cache and "--write-cache" in sys.argv:
            json.dump({"vars": GV, "terms": terms}, open(cache, "w"))
    h = canon_hash_terms(GV, terms)
    check("G: 1772 terms, homogeneous of degree 9, sha256 " + h + " (" + src + ")",
          len(terms) == 1772 and all(sum(m) == 9 for m, _ in terms) and h == G_SHA256)
    return {tuple(m): c for m, c in terms}


def setup_checks():
    check("paper T_49 = {0,3,5,6,9,12}", L.CONTACT_SETS[49] == T49)
    gens = L.generators(T49)
    check("generators(T_49) = f0,f1,f2,f5,f7,f9,f29,f39,f40,f41 (ten)", gens == [0, 1, 2, 5, 7, 9, 29, 39, 40, 41], gens)
    eight = list(ROLES.values())
    Z8 = frozenset(s for s in range(24) if all(L.F_VALUES[g][s] == 0 for g in eight))
    check("the eight role generators vanish together exactly on T_49", Z8 == T49)
    check("d = f40 is the only non-minor role", [n for n, r in ROLES.items() if not L.F_IS_MINOR[r]] == ["d"])
    # f1, f9 in the span of the eight roles
    for extra, rel in ((1, dict(q=1, c=1, a=1, d=-1)), (9, dict(q=1, e=1, b=1, d=-1))):
        ok = all(sum(rel.get(n, 0) * L.F_VALUES[ROLES[n]][s] for n in ROLES) == L.F_VALUES[extra][s] for s in range(24))
        check("f%d = %s (so a ten-generator source has an eight-role form with d replaced by d - x1 - x9)" %
              (extra, " + ".join(("-" if v < 0 else "") + "f%d" % ROLES[n] for n, v in rel.items())), ok)


def fperms():
    fv = {L.F_VALUES[j]: j for j in range(82)}
    out = []
    for t in L._G_MAPS:
        row = []
        for j in range(82):
            new = [None] * 24
            for s in range(24):
                new[t[s]] = L.F_VALUES[j][s]
            row.append(fv[tuple(new)])
        out.append(row)
    return out


def atom_checks(chains=(), oldH=()):
    FP = fperms()
    chain_orbit = {(p[2], frozenset((p[0], p[12])), p[70]) for p in FP}   # representative (2; 0,12; 70)
    oldH_orbit = {(p[3], p[15], p[26]) for p in FP}                       # representative (3, 15; 26)
    check("chain orbit has 576 members, old-H orbit 288", len(chain_orbit) == 576 and len(oldH_orbit) == 288)
    for centre, leaves, k in chains:
        check("chain (centre f%d; leaves f%d, f%d; non-minor f%d) is in the proved chain orbit" % (centre, leaves[0], leaves[1], k),
              (centre, frozenset(leaves), k) in chain_orbit)
    for aux, main, k in oldH:
        check("old-H (aux f%d, main f%d; f%d) is in the proved old-H orbit (reversed roles: %s)" %
              (aux, main, k, (main, aux, k) in oldH_orbit), (aux, main, k) in oldH_orbit)


def zero_face_checks(names=("q", "p", "j", "e", "b", "a"), claimed=None):
    """claimed: dict role -> paper case number.  Checks that the paper case embeds (up to G)
    in the contact set of the face where that coefficient vanishes."""
    eight = list(ROLES.values())
    for n in names:
        rest = [g for g in eight if g != ROLES[n]]
        Z = frozenset(s for s in range(24) if all(L.F_VALUES[g][s] == 0 for g in rest))
        k = claimed[n]
        hit = None
        for six in itertools.combinations(sorted(Z), 6):
            g = L.find_g(L.CONTACT_SETS[k], six)
            if g is not None:
                hit = (six, g); break
        check("zero face %s = 0: contacts %s contain an image of paper T_%d" % (n, sorted(Z), k), hit is not None,
              "image %s via g=%s" % (hit[0], hit[1]) if hit else "")


def q0_checks():
    q, p, j, c, e, d, b, a = sp.symbols(GV)
    sym = dict(zip(GV, (q, p, j, c, e, d, b, a)))
    m = [sum(sym[n] * L.F_VECS[r][k] for n, r in ROLES.items()) for k in range(14)]
    hh = sp.Symbol("h"); m[13] += hh
    Cu = sp.Poly(sp.expand(L.cubic_C(m)), hh)
    C0 = Cu.coeff_monomial(1); Q0 = Cu.coeff_monomial(hh)
    check("Q0 = (c+d)(d+e+j) - p(e+j)", sp.expand(Q0 - ((c + d) * (d + e + j) - p * (e + j))) == 0)
    X = sp.Poly(sp.expand((e + j) * C0 - e * j * Q0), q, p, j, c, e, d, b, a)
    cs = [cc for _, cc in X.terms()]
    check("(e+j)C - ejQ0 has 42 terms, all positive, including d^3 e and d^3 j",
          len(cs) == 42 and min(cs) > 0 and X.coeff_monomial(d**3 * e) > 0 and X.coeff_monomial(d**3 * j) > 0)
    C00 = sp.expand(C0.subs(c, 0)); Q00 = sp.expand(Q0.subs(c, 0))
    Z = sp.Poly(sp.expand(sp.cancel((e + j) * C00.subs(p, d * (d + e + j) / (e + j)))), q, p, j, c, e, d, b, a)
    cz = [cc for _, cc in Z.terms()]
    check("c = 0, Q0 = 0: (e+j)C0 is a 27-term polynomial with nonnegative coefficients containing d^3(e+j)",
          len(cz) == 27 and min(cz) > 0 and Z.coeff_monomial(d**3 * e) > 0)
    return C0, Q0


# ---- certificate engine ------------------------------------------------------------
def subst(G, fr, ctx):
    """fr: var -> (num, den).  Returns (R, Pi) with R = G(num/den) * Pi, Pi = prod den^deg."""
    degs = {v: max(m[i] for m in G) for i, v in enumerate(GV)}
    one = ctx.from_dict({}) + 1
    pw = {}
    def power(v, w, k):
        if (v, w, k) not in pw:
            pw[(v, w, k)] = one if k == 0 else power(v, w, k - 1) * fr[v][w]
        return pw[(v, w, k)]
    def rec(terms, i):
        if i == len(GV):
            return one * sum(c for _, c in terms)
        v = GV[i]; groups = {}
        for mm, c in terms:
            groups.setdefault(mm[i], []).append((mm, c))
        s = ctx.from_dict({})
        for k, sub in groups.items():
            s = s + rec(sub, i + 1) * power(v, 0, k) * power(v, 1, degs[v] - k)
        return s
    R = rec(list(G.items()), 0)
    Pi = one
    for v in GV:
        Pi = Pi * power(v, 1, degs[v])
    return R, Pi


def strip(R, factors):
    removed = {}
    for name, f in factors:
        while True:
            qq, rr = divmod(R, f)
            if not rr.is_zero():
                break
            R = qq; removed[name] = removed.get(name, 0) + 1
    return R, removed


def signs(P):
    cs = [int(c) for c in P.to_dict().values()]
    return len(cs), sum(1 for c in cs if c > 0), sum(1 for c in cs if c < 0)


def coeff(P, names, mon):
    return int(P.to_dict().get(tuple(mon.get(n, 0) for n in names), 0))


def poly_hash(P, names):
    return canon_hash_terms(names, [(m, c) for m, c in P.to_dict().items()])


def square_groups(P, names, var):
    iv = names.index(var); groups = {}
    for mm, c in P.to_dict().items():
        mm = tuple(int(x) for x in mm)
        groups.setdefault(mm[:iv] + mm[iv + 1:], {})[mm[iv]] = int(c)
    bad = 0; ok = True; nneg = 0
    for key, g in groups.items():
        if min(g.values()) >= 0:
            continue
        bad += 1; nneg += sum(1 for c in g.values() if c < 0)
        co = [g.get(k, 0) for k in range(max(g) + 1)]
        for _ in range(2):
            hi = co[::-1]; qd = [hi[0]]
            for x in hi[1:]:
                qd.append(x + qd[-1])
            if qd.pop() != 0:
                ok = False; break
            co = qd[::-1]
        else:
            if min(co) < 0:
                ok = False
    return bad, ok, nneg, groups


class Q_:
    """Reduced fraction of flint polynomials."""
    def __init__(s, n, d=None):
        d = d if d is not None else n * 0 + 1
        g = n.gcd(d)
        if not g.is_one() and not g.is_zero():
            n = n // g; d = d // g
        if d.leading_coefficient() < 0:
            n = -n; d = -d
        s.n = n; s.d = d
    def _c(s, o): return o if isinstance(o, Q_) else Q_(s.n * 0 + o)
    def __add__(s, o): o = s._c(o); return Q_(s.n * o.d + o.n * s.d, s.d * o.d)
    __radd__ = __add__
    def __sub__(s, o): o = s._c(o); return Q_(s.n * o.d - o.n * s.d, s.d * o.d)
    def __rsub__(s, o): return s._c(o) - s
    def __mul__(s, o): o = s._c(o); return Q_(s.n * o.n, s.d * o.d)
    __rmul__ = __mul__
    def __truediv__(s, o): o = s._c(o); return Q_(s.n * o.d, s.d * o.n)
    def __rtruediv__(s, o): return s._c(o) / s
    def eq(s, o): o = s._c(o); return (s.n * o.d - o.n * s.d).is_zero()


def finish():
    print("SUMMARY: %d EXACT check(s) failed" % len(FAILS) + ((": " + ", ".join(FAILS)) if FAILS else ""))
    sys.exit(1 if FAILS else 0)
# ---- end of shared core ---------------------------------------------------------------


def numeric_capacity(tag, sampler, n=12, seed=0, restarts=8):
    """NUMERIC: for n random coefficient vectors from sampler(rng) (dict of the eight roles,
    possibly not spectrally realizable), minimise the Leibniz sum over SO(4) (multistart BFGS,
    lib.min_phi_over_SO4) and report the worst value relative to the coefficient scale."""
    import numpy as np
    rng = np.random.default_rng(seed)
    worst = None
    for _ in range(n):
        co = sampler(rng)
        cvec = [sum(float(co[k]) * float(L.F_VECS[r][i]) for k, r in ROLES.items()) for i in range(14)]
        res = L.min_phi_over_SO4(cvec, restarts=restarts, seed=int(rng.integers(1 << 30)))
        val = res[0] if isinstance(res, (tuple, list)) else res
        val = float(val) / max(1.0, sum(abs(float(x)) for x in co.values()))
        worst = val if worst is None else min(worst, val)
    print("NUMERIC %s: %d samples, seed %d, restarts %d, worst min Phi / scale = %.3e" % (tag, n, seed, restarts, worst), flush=True)
    return worst
ZF = dict(q=39, p=1, j=3, e=24, b=43, a=44)   # historical T38,T0,T3,T23,T42,T43 = paper T_39,T_1,T_3,T_24,T_43,T_44
setup_checks()
C0s, Q0s = q0_checks()
zero_face_checks(claimed=ZF)
# ===== blind_moc4_T48_allocation_I_low_j_uniform_20261023 =====
atom_checks(chains=[(5, (0, 2), 40), (5, (7, 29), 40), (5, (29, 41), 40)])
a, e, j, q, b = sp.symbols("a e j q b", positive=True)
Q = q + b; Ea = a * e / (a + e); JQ = j * Q / (j + Q); lm = 1 / (Ea + JQ); u = lm * Ea; v = lm * JQ
req = u**2 * (1 / a + 1 / e) + (v * q / Q)**2 * (1 / q + 1 / (j * q / Q)) + (v * b / Q)**2 * (1 / b + 1 / (j * b / Q))
check("both-chain allocation: requirements sum to lambda; currents u + v = 1; j shares jq/Q + jb/Q = j",
      sp.simplify(req - lm) == 0 and sp.simplify(u + v - 1) == 0)
# link B0 = 0 (Q = 2v) to the allocation-II cell: allocation II's u is the unique solution of u(1/a+1/e) - (1-u)/j = 1/2
# the identity u/Ea - v/j = lambda - (lambda - v/Q) = v/Q; at Q = 2v this is 1/2
check("u(1/a+1/e) - v/j = v/Q identically (so = 1/2 exactly on Q = 2v)", sp.simplify(u * (1 / a + 1 / e) - v / j - v / Q) == 0)
G = load_or_compute_G(); t0 = time.time()
names = ["X", "A", "B", "C1", "Z", "P", "c"]
ctx = flint.fmpz_mpoly_ctx.get(tuple(names), 'degrevlex')
X, A, B, C1, Z, P, c = [Q_(g) for g in ctx.gens()]
gX, gA, gB, gC1, gZ, gP, gc = ctx.gens()
expect = {"c1": (245193, {}, 10616832), "c2": (166518, {"X": 3}, 131072), "c3": (457338, {"X": 13}, 1296)}
for ch in ("c1", "c2", "c3"):
    if ch == "c1":
        Lr = 1 / (1 + X); B0 = B; Tj = (2 + 2 * X + B * X) / ((2 + X) * (1 + B)); T = Tj + Z
    elif ch == "c2":
        Lr = 1 + X; B0 = 2 / X + B; Tj = 2 / (1 + Lr) - B0 / (1 + B0); T = Z
    else:
        B0 = X; Lr = 1 + 2 * B / (X * (1 + B)); Tj = X / ((X * (1 + B) + B) * (1 + X)); T = Tj + Z
    u = Lr / (1 + Lr); v = 1 / (1 + Lr)
    TjF = 2 / (1 + Lr) - B0 / (1 + B0)
    aa = 2 * u * (1 + A); QQ = 2 * v * (1 + B0)
    ee = 2 * u * (1 + A) / (A + T + A * T); jj = 2 * v * (1 + B0) / (B0 + T + B0 * T)
    pp = (1 + T) / 2 * P / (1 + P); qq = QQ * C1 / (1 + C1); bb = QQ / (1 + C1)
    Ea_ = aa * ee / (aa + ee); JQ_ = jj * QQ / (jj + QQ); lam_ = 1 / (Ea_ + JQ_)
    check(ch + ": Tj = 2/(1+L) - B0/(1+B0); lambda = (1+T)/2; u = lambda Ea; v = lambda JQ; a = 2u(1+A); Q = 2v(1+B0)",
          Tj.eq(TjF) and lam_.eq((1 + T) / 2) and u.eq(lam_ * Ea_) and v.eq(lam_ * JQ_))
    check(ch + ": (1 - j)(B0 + T + B0 T) = (1+B0)(T - Tj)  [so j <= 1 iff T >= Tj]", ((1 - jj) * (B0 + T + B0 * T)).eq((1 + B0) * (T - TjF)))
    if ch == "c2":
        check("c2: Tj = (2 - B0 X)/((2+X)(1+B0)) <= 0 and B0(L-1) = 2 + BX", Tj.eq((2 - B0 * X) / ((2 + X) * (1 + B0))) and (B0 * (Lr - 1)).eq(2 + B * X))
    if ch == "c3":
        check("c3: B0(L-1) = 2B/(1+B)", (B0 * (Lr - 1)).eq(2 * B / (1 + B)))
    fr = dict(a=aa, e=ee, j=jj, q=qq, b=bb, p=pp, c=c, d=X * 0 + 1)
    R, Pi = subst(G, {k: (w.n, w.d) for k, w in fr.items()}, ctx)
    facs = [("X", gX), ("1+X", 1 + gX), ("2+X", 2 + gX), ("1+B", 1 + gB), ("B", gB), ("A", gA), ("1+A", 1 + gA),
            ("1+C1", 1 + gC1), ("1+P", 1 + gP), ("Z", gZ), ("e_den", (A + T + A * T).n), ("j_den", (B0 + T + B0 * T).n),
            ("T_den", T.d), ("L_den", Lr.d), ("B0_den", B0.d), ("XB+X+B", gX * gB + gX + gB)]
    facs = [(nm, f) for nm, f in facs if not f.is_constant()]
    Rs, rem = strip(R, facs)
    nt, npos, nneg = signs(Rs)
    want, mk, mkc = expect[ch]
    print("INFO %s: cleared by %s; remaining factors removed %s; %d terms; sha256 %s; t=%.0f" %
          (ch, "prod den^deg", rem, nt, poly_hash(Rs, names), time.time() - t0), flush=True)
    check(ch + ": substituted G (positive denominators cleared) has all %d coefficients of one sign (supplement table: %d terms)" % (nt, want), nneg == 0)
    check(ch + ": marker monomial %s present with positive coefficient (supplement: -%d)" % (mk or "1", mkc), coeff(Rs, names, mk) > 0,
          "coef %d" % coeff(Rs, names, mk))
def samp(rng):
    while True:
        a_, e_, j_, Qq = rng.uniform(0.05, 6, 4)
        Ea_ = a_ * e_ / (a_ + e_); JQ_ = j_ * Qq / (j_ + Qq); l_ = 1 / (Ea_ + JQ_)
        if l_ >= 0.5 and j_ <= 1 and 2 * l_ * Ea_ <= a_ and 2 * l_ * JQ_ <= Qq:
            fq = rng.uniform()
            return dict(q=Qq * fq, p=l_, j=j_, c=rng.uniform(0, 2), e=e_, d=1, b=Qq * (1 - fq), a=a_)
numeric_capacity("both-chain capacity at p = lambda (large-lambda low-j cell, not necessarily realizable)", samp, n=12, seed=13)
finish()
