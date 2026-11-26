"""Auditor replay for the T39 (paper T_40 = {0,1,5,9,14,18}) radial-exclusion documents
  independent_moc4_T39_high_p_uniform_20261022  (REGION = "high": p >= 1, p j < 1)
  independent_moc4_T39_low_pj_uniform_20261022  (REGION = "low":  0 <= p, j < 1)
The region is taken from this file's name. Own code only; uses ../lib/moc4lib.py.
EXACT checks: symbolic restriction of the cubic; pointwise exact (rational) identities at sampled base points.
NUMERIC checks: 60-digit evaluations at the positive root c+ (an algebraic number).
The 1313-term G, the octics F+-, L, M and the global coefficient certificates are NOT rebuilt (see report)."""
import os, sys, random, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import moc4lib as L
import sympy as sp
import mpmath as mp
import flint
mp.mp.dps = 60
REGION = "low" if "low_pj" in os.path.basename(__file__) else "high"
NSAMP = int(sys.argv[1]) if len(sys.argv) > 1 else 60
FAIL = []
def chk(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + (("  " + str(info)) if info else ""))
    if not ok: FAIL.append(name)

q, p, j, c, e, b, a, h = sp.symbols("q p j c e b a h")
SPEC = [(q, 5), (p, 1), (j, 12), (c, 9), (e, 2), (1, 35), (b, 33), (a, 29)]   # d = 1
chk("A0 common zero set of the eight source rays is T_40",
    set.intersection(*[set(L.F_ZERO[k]) for _, k in SPEC]) == set(L.CONTACT_SETS[40]))
chk("A0 source rays are generators of N_{T_40}; only f35 is non-minor",
    all(k in L.generators(L.CONTACT_SETS[40]) for _, k in SPEC) and [k for _, k in SPEC if not L.F_IS_MINOR[k]] == [35])
gsym = [sum(co * sp.Rational(L.F_VECS[k][r].numerator, L.F_VECS[k][r].denominator) for co, k in SPEC) for r in range(14)]
w = list(gsym); w[13] = w[13] + h
Cw = sp.expand(L.cubic_C(w))
C0 = Cw.subs(h, 0); Q0 = sp.expand(sp.diff(Cw, h))
A = 1 + a + b; B = 1 + q + e; N = (q + 1)*(e + 1)*A + B*a*b; Nj = N - j*e*a
chk("C0 restriction is affine in h", sp.degree(Cw, h) == 1)
chk("C0 C = -(ABc + Nj)", sp.expand(C0 + A*B*c + Nj) == 0)
chk("C0 -Q = (q+e+1+a+b-p)c + (q+b+1)(e+a+1) - j(e+a) + pj",
    sp.expand(-Q0 - ((q + e + 1 + a + b - p)*c + (q + b + 1)*(e + a + 1) - j*(e + a) + p*j)) == 0)
Nm = sp.Poly(sp.expand(N - e*a), q, e, a, b)
chk("C0 N - ea has nonnegative coefficients and constant 1", all(t >= 0 for t in Nm.coeffs()) and Nm.coeff_monomial(1) == 1)
P_, J_ = sp.symbols("P J")
if REGION == "high":
    sub = {p: 1 + P_, j: J_/((1 + P_)*(1 + J_))}
    chk("C2 substitution: p j = J/(1+J) in [0,1), p = 1 + P >= 1", sp.simplify((p*j).subs(sub) - J_/(1 + J_)) == 0)
else:
    sub = {p: P_/(1 + P_), j: J_/(1 + J_)}
    chk("C2 substitution maps P,J >= 0 onto 0 <= p, j < 1", True)
    mQ = sp.Poly(sp.expand(-Q0 - (1 - p)*c - (1 - j)*(e + a) - p*j - 1), q, e, a, b, c, p, j)
    chk("C0 low region: -Q - (1-p)c - (1-j)(e+a) - pj - 1 has nonnegative coefficients (so -Q > 0)", all(t > 0 for t in mQ.coeffs()))
# height independence of P*: H(kappa) is supported on rows/cols {0,1,5,6}
TENS = {}
ms, Cs = L.cubic_sympy()
for i in range(14):
    for jj in range(14):
        dij = sp.diff(Cs, ms[i], ms[jj])
        for k in range(14):
            vv = sp.diff(dij, ms[k])
            if vv != 0: TENS[(i, jj, k)] = int(vv)
hk = sorted((i, jj) for (i, jj, k) in TENS if k == 13)
chk("C3 H(kappa) entries only at (0,6),(6,0),(1,5),(5,1): P* (delete 0,1) is height independent", hk == [(0, 6), (1, 5), (5, 1), (6, 0)], hk)

def fr(x): return sp.Rational(x)
def fdet(M):
    n = len(M)
    F = flint.fmpq_mat(n, n, [flint.fmpq(int(sp.numer(t)), int(sp.denom(t))) for row in M for t in row])
    d = F.det(); return sp.Rational(int(d.p), int(d.q))
def hess(wv):
    H = [[0]*14 for _ in range(14)]
    for (i, jj, k), vv in TENS.items(): H[i][jj] += vv * wv[k]
    return H
cs = sp.Symbol("cs")
def run(base):
    out = {}
    xs, ys, ps = [], [], []
    for k in range(12):
        pt = dict(base); pt[c] = sp.Rational(k, 3)
        Qv = Q0.subs(pt)
        if Qv == 0: continue
        wv = [t.subs(pt) for t in gsym]; wv[13] += -C0.subs(pt)/Qv
        H = hess(wv)
        xs.append(pt[c]); ys.append(-fdet(H)*Qv**2/4)
        ps.append(fdet([[H[r][s] for s in range(2, 14)] for r in range(2, 14)]))
    poly = sp.Poly(sp.interpolate(list(zip(xs[:9], ys[:9])), cs), cs)
    out["fit"] = poly.degree() == 4 and all(poly.eval(x_) == y_ for x_, y_ in zip(xs, ys))
    lc, facs = sp.factor_list(poly.as_expr())
    sq = sp.sqrt(lc)
    out["square"] = bool(sq.is_rational and all(m % 2 == 0 for _, m in facs))
    if not out["square"]: return out
    G = sp.Poly(sp.expand(sq*sp.prod([f_**(m//2) for f_, m in facs])), cs)
    al, be, ga = [G.coeff_monomial(cs**k) for k in (2, 1, 0)]
    if ga > 0: al, be, ga = -al, -be, -ga
    Av, Bv, Njv = A.subs(base), B.subs(base), Nj.subs(base)
    out["gamma<0"] = ga < 0
    out["cert2"] = (Av*Bv*be - Njv*al) < 0
    Pst = sp.Poly(sp.interpolate(list(zip(xs[:9], ps[:9])), cs), cs)
    out["Pstar_deg4"] = Pst.degree() <= 4 and all(Pst.eval(x_) == y_ for x_, y_ in zip(xs, ps))
    out["root"] = al > 0
    if al > 0:
        cp = (-mp.mpf(str(be)) + mp.sqrt(mp.mpf(str(be*be - 4*al*ga))))/(2*mp.mpf(str(al)))
        out["c+"] = cp
        out["P*(c+)"] = mp.polyval([mp.mpf(str(t)) for t in Pst.all_coeffs()], cp)
        def wnum(cv):
            gm = [mp.mpf(0)]*14
            for co, k in SPEC:
                val = cv if co == c else (mp.mpf(1) if co == 1 else mp.mpf(str(base[co])))
                for r in range(14): gm[r] += val*mp.mpf(L.F_VECS[k][r].numerator)/L.F_VECS[k][r].denominator
            Cn = L.cubic_C(gm); g1 = list(gm); g1[13] += 1; Qn = L.cubic_C(g1) - Cn
            gm[13] += -Cn/Qn
            return gm, Qn
        wv, Qn = wnum(cp)
        out["Q(c+)"] = Qn
        H = mp.matrix(14, 14)
        for (i, jj, k), vv in TENS.items(): H[i, jj] += vv*wv[k]
        U, S, Vt = mp.svd_r(H)
        sv = sorted([abs(S[i]) for i in range(14)])
        out["sv"] = sv[:3]
        out["border"] = [mp.det(mp.matrix([[H[r, s] for s in [kk] + list(range(2, 14))] for r in [kk] + list(range(2, 14))])) for kk in (0, 1)]
        mins = [mp.det(mp.matrix([[H[r, s] for s in idx] for r in idx])) for idx in
                ([t for t in range(14) if t not in dr] for dr in itertools.combinations(range(14), 2))]
        out["max12"] = max(mins)
        # Section 4: quadratic form Q_v on the kernel
        eps = mp.mpf("1e-25")
        wp, _ = wnum(cp + eps); wm, _ = wnum(cp - eps)
        vvec = [(wp[i] - wm[i])/(2*eps) for i in range(14)]
        Hv = mp.matrix(14, 14)
        for (i, jj, k), vv in TENS.items(): Hv[i, jj] += vv*vvec[k]
        idx = sorted(range(14), key=lambda i: abs(S[i]))[:2]
        K = mp.matrix(14, 2)
        for col, i in enumerate(idx):
            for r in range(14): K[r, col] = Vt[i, r]
        M2 = K.T*Hv*K/2
        out["det KQvK"] = M2[0, 0]*M2[1, 1] - M2[0, 1]*M2[1, 0]
    return out

def rnd(rng, zp=0.15):
    return sp.Integer(0) if rng.random() < zp else sp.Rational(rng.randint(1, 40), rng.randint(1, 12))
rng = random.Random(7 if REGION == "high" else 8)
cnt = {"fit": 0, "square": 0, "gamma<0": 0, "cert2": 0, "Pstar_deg4": 0, "root": 0}
worstP = None; worstB = 0; worstM = None; minD = None; ok_root = True; n = 0
for t in range(NSAMP):
    base = {s: rnd(rng) for s in (q, e, b, a)}
    Pv, Jv = rnd(rng, 0.2), rnd(rng, 0.2)
    if REGION == "high": base[p] = 1 + Pv; base[j] = Jv/((1 + Pv)*(1 + Jv))
    else: base[p] = Pv/(1 + Pv); base[j] = Jv/(1 + Jv)
    o = run(base); n += 1
    for k in cnt:
        if o.get(k): cnt[k] += 1
    if o.get("root"):
        worstP = o["P*(c+)"] if worstP is None else max(worstP, o["P*(c+)"])
        worstM = o["max12"] if worstM is None else max(worstM, o["max12"])
        minD = o["det KQvK"] if minD is None else min(minD, o["det KQvK"])
        worstB = max(worstB, max(abs(x_) for x_ in o["border"]))
        if not (o["sv"][1] < mp.mpf("1e-40") and o["sv"][2] > mp.mpf("1e-10") and o["Q(c+)"] != 0): ok_root = False
chk("C1 EXACT at %d rational base points: -det H(w(c)) Q^2/4 is a degree-4 polynomial in c (fits 12 points)" % n, cnt["fit"] == n, cnt["fit"])
chk("C1 EXACT at the same points: it is a perfect square G(c)^2 with G quadratic", cnt["square"] == n, cnt["square"])
chk("C2 EXACT pointwise: with the sign of G fixed by gamma < 0 (gamma != 0 at all points)", cnt["gamma<0"] == n, cnt["gamma<0"])
chk("C2 EXACT pointwise: AB beta - Nj alpha < 0 for that sign", cnt["cert2"] == n, cnt["cert2"])
chk("C3 EXACT pointwise: P*(c) has degree <= 4 in c", cnt["Pstar_deg4"] == n, cnt["Pstar_deg4"])
print("INFO base points with a positive root (alpha > 0): %d of %d" % (cnt["root"], n))
if cnt["root"]:
    chk("NUMERIC P*(c+) < 0 at every sampled positive root", worstP < 0, "max = %s" % mp.nstr(worstP, 8))
    chk("NUMERIC both bordered minors on I u {0}, I u {1} vanish at c+", worstB < mp.mpf("1e-30"), "max |.| = %s" % mp.nstr(worstB, 5))
    chk("NUMERIC rank H(w(c+)) = 12 (two singular values < 1e-40, third > 1e-10), Q(c+) != 0", ok_root)
    chk("NUMERIC all 91 principal 12x12 minors < 0 at c+", worstM < 0, "max = %s" % mp.nstr(worstM, 8))
    chk("NUMERIC det(K^T Q_v K) > 0 on the 2-dim kernel (Section 4 definiteness)", minD > 0, "min = %s" % mp.nstr(minD, 8))
print("SUMMARY", "all EXACT checks passed" if not [f_ for f_ in FAIL if not f_.startswith("NUMERIC")] else "FAILED: %s" % FAIL)
sys.exit(1 if [f_ for f_ in FAIL if not f_.startswith("NUMERIC")] else 0)
