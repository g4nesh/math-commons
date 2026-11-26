# Certificate replay (auditor cert-T38) for blind_moc4_T38_uniform_20261022 (paper case T_39).
# Own code only: nothing from the supplement is executed or imported.  Every polynomial is
# rebuilt from the paper's definitions (cubic C of (A.2), generator vectors f_j, kappa) with
# python-flint.  Run:  python -I -B blind_moc4_T38_uniform_20261022__certificates.py
# Optional: MOC4_SUPP=<supplement root> only to list the supplement's T38 data files (none are used).
import os, sys, time, glob
from fractions import Fraction as Fr
_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, "..", "lib"))
from moc4lib import F_VECS, cubic_C, generators, CONTACT_SETS
import flint
import sympy as sp

T0 = time.time()
FAIL = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + ((" : " + str(info)) if info != "" else ""), flush=True)
    if not ok: FAIL.append(name)
def info(name, msg): print("INFO " + name + " : " + str(msg), flush=True)
def nt(P): return len(P.coeffs())
def nonneg(P): return all(int(c) >= 0 for c in P.coeffs())
def nneg(P): return sum(1 for c in P.coeffs() if int(c) < 0)

supp = os.environ.get("MOC4_SUPP")
if supp:
    hits = sorted(os.path.basename(x) for x in glob.glob(os.path.join(supp, "proofs", "manin_transfer_20261021", "*T38*")))
    info("supplement files with T38 in the name", hits)

check("contact set is paper T_39 and its generators are {0,1,2,5,7,12,39,40,67,70}",
      CONTACT_SETS[39] == frozenset({2, 3, 5, 7, 9, 12}) and generators(CONTACT_SETS[39]) == [0, 1, 2, 5, 7, 12, 39, 40, 67, 70])

# ---------------------------------------------------------------- 1. G40 and G70 from the Hessian
mctx = flint.fmpz_mpoly_ctx.get(("m", 14), "lex")
MS = mctx.gens()
CM = cubic_C(list(MS))
HM = [[CM.derivative(i).derivative(j) for j in range(14)] for i in range(14)]

def bareiss(M):
    M = [row[:] for row in M]; n = len(M); sign = 1; prev = None
    for k in range(n):
        best = None
        for i in range(k, n):
            for j in range(k, n):
                if not M[i][j].is_zero():
                    L = nt(M[i][j])
                    if best is None or L < best[0]: best = (L, i, j)
        if best is None: return M[0][0] * 0
        _, i, j = best
        if i != k: M[k], M[i] = M[i], M[k]; sign = -sign
        if j != k:
            for row in M: row[k], row[j] = row[j], row[k]
            sign = -sign
        p = M[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                v = p * M[i][j] - M[i][k] * M[k][j]
                if prev is not None:
                    qq, rr = divmod(v, prev); assert rr.is_zero(); v = qq
                M[i][j] = v
        prev = p
    return sign * M[n - 1][n - 1]

def radial(names, pairs):
    """w = sum coef*f_k (d = 1 for the non-minor generator) + h*kappa; returns ctx (names + h), C0, Q, G
    with C(w) = C0 + h Q and Q^4 det Hess C(w)|_{h=-C0/Q} = -4 Q^2 G^2."""
    ctx = flint.fmpz_mpoly_ctx.get(tuple(names) + ("h",), "lex")
    g = ctx.gens(); hv = g[-1]; vm = dict(zip(names, g))
    w = [ctx.constant(0)] * 14
    for coef, k in pairs:
        cv = vm[coef] if isinstance(coef, str) else ctx.constant(coef)
        for r in range(14):
            assert Fr(F_VECS[k][r]).denominator == 1
            w[r] = w[r] + cv * int(F_VECS[k][r])
    w[13] = w[13] + hv
    C = cubic_C(w)
    assert C.degrees()[-1] == 1
    Q = C.derivative(len(names)); C0 = C - Q * hv
    H = [[HM[i][j].compose(*w, ctx=ctx) if not HM[i][j].is_zero() else ctx.constant(0) for j in range(14)] for i in range(14)]
    D = bareiss(H)
    dh = D.degrees()[-1]
    dd = [dict() for _ in range(dh + 1)]
    for mon, cf in zip(D.monoms(), D.coeffs()):
        e = list(mon); k = e[-1]; e[-1] = 0; dd[k][tuple(e)] = cf
    Dk = [ctx.from_dict(x) for x in dd]
    N = sum((Dk[k] * (-C0) ** k * Q ** (4 - k) for k in range(dh + 1)), ctx.constant(0))
    G2, rem = divmod(N, -4 * Q * Q)
    G = G2.sqrt() if rem.is_zero() else None
    ok = rem.is_zero() and G is not None and G * G == G2
    # drop h from the context
    c2 = flint.fmpz_mpoly_ctx.get(tuple(names), "lex")
    def drop(P): return c2.from_dict({m[:-1]: c for m, c in zip(P.monoms(), P.coeffs())})
    if ok and int(drop(G).coefficient(0) if False else 0) == 0: pass
    Gd = drop(G); Gd = Gd if int(Gd.coeffs()[-1]) > 0 else -Gd   # sign: positive constant term (lex: last)
    return c2, drop(C0), drop(Q), Gd, ok, dh

N40 = ("q", "p", "j", "c", "e", "b", "a")
P40 = [("q", 39), ("p", 12), ("j", 5), ("c", 7), ("e", 2), (1, 40), ("b", 67), ("a", 0)]
N70 = ("p", "r", "a", "b", "c", "q", "s")
P70 = [("p", 0), ("r", 1), ("a", 2), ("b", 12), ("c", 39), ("q", 67), (1, 70), ("s", 5)]
c40, C40, Q40, G40, ok40, dh40 = radial(N40, P40)
c70, C70, Q70, G70, ok70, dh70 = radial(N70, P70)
check("g40: Q^4 detHess(C)|_{h=-C/Q} = -4 Q^2 G40^2 exactly (polynomial identity)", ok40, f"deg_h detHess = {dh40}")
check("g40: G40 has 1455 terms, total degree 9", nt(G40) == 1455 and G40.total_degree() == 9, (nt(G40), G40.total_degree()))
check("g70: Q^4 detHess(C)|_{h=-C/Q} = -4 Q^2 G70^2 exactly (polynomial identity)", ok70, f"deg_h detHess = {dh70}")
check("g70: G70 has 1113 terms, total degree 9", nt(G70) == 1113 and G70.total_degree() == 9, (nt(G70), G70.total_degree()))
info("time after G40, G70", round(time.time() - T0, 1))

# ---------------------------------------------------------------- helpers: chart substitution
def to_flint(expr, ctx, syms):
    P = sp.Poly(sp.expand(expr), *syms)
    d = {}
    for mon, cf in P.terms():
        assert cf.is_integer, expr
        d[tuple(mon)] = int(cf)
    return ctx.from_dict(d) if d else ctx.constant(0)

def numden(expr, ctx, syms):
    n, d = sp.fraction(sp.cancel(sp.together(expr)))
    return to_flint(n, ctx, syms), to_flint(d, ctx, syms)

def subst(P, srcnames, mapping, ctx2, degs=None):
    """P(src vars) with x -> num_x/den_x; returns P(...)*prod den_x^{D_x} with D_x = degs[x] (default deg_x P)."""
    D = degs or P.degrees()
    hn = []
    for nm in srcnames: hn += [nm + "_n", nm + "_d"]
    hctx = flint.fmpz_mpoly_ctx.get(tuple(hn), "lex")
    hd = {}
    for mon, cf in zip(P.monoms(), P.coeffs()):
        e = []
        for k, ek in enumerate(mon): e += [ek, D[k] - ek]
        hd[tuple(e)] = cf
    PH = hctx.from_dict(hd)
    args = []
    for nm in srcnames:
        nu, de = mapping[nm]; args += [nu, de]
    return PH.compose(*args, ctx=ctx2)

def strip(P):
    """divide by integer content and monomial content (positive factors), for size comparison."""
    if P.is_zero(): return P
    c = P.content(); P = P / c if int(c) != 1 else P
    tc = P.term_content()
    q, r = divmod(P, tc); assert r.is_zero()
    return q

def posden(d):
    return (not d.is_zero()) and all(int(c) > 0 for c in d.coeffs())

# exact Sturm count of distinct real roots of an integer univariate polynomial in (0, +inf)
def sturm_pos_roots(coeffs):
    f = flint.fmpq_poly(coeffs)
    seq = [f, f.derivative()]
    while not seq[-1].is_zero():
        r = seq[-2] % seq[-1]
        if r.is_zero(): break
        seq.append(-r)
    def var(vals):
        s = [v for v in vals if v != 0]
        return sum(1 for x, y in zip(s, s[1:]) if (x > 0) != (y > 0))
    at0 = var([p(0) for p in seq])
    atinf = var([p.leading_coefficient() if not p.is_zero() else 0 for p in seq])
    return at0 - atinf, at0, atinf

def group_check(P, varidx):
    """Group P by the monomials of all variables except varidx; each group is a univariate
    polynomial in that variable.  Return (#groups, #groups needing a sign check, #failures)."""
    groups = {}
    for mon, cf in zip(P.monoms(), P.coeffs()):
        key = mon[:varidx] + mon[varidx + 1:]
        groups.setdefault(key, {})[mon[varidx]] = int(cf)
    need = 0; bad = []
    for key, g in groups.items():
        if all(v >= 0 for v in g.values()): continue
        need += 1
        lo = min(g); hi = max(g)
        co = [g.get(k, 0) for k in range(lo, hi + 1)]
        nroots, v0, vinf = sturm_pos_roots(co)
        if not (co[0] > 0 and co[-1] > 0 and nroots == 0):
            bad.append((key, co))
    return len(groups), need, bad

check("Sturm self-test: E^2-3E+2 has 2 roots in (0,inf), E^2-2E+2 none, (E-1)^2 one, E^2+1 none",
      [sturm_pos_roots(x)[0] for x in ([2, -3, 1], [2, -2, 1], [1, -2, 1], [1, 0, 1])] == [2, 0, 1, 0])
_tc = flint.fmpz_mpoly_ctx.get(("x", "E"), "lex"); _x, _E = _tc.gens()
check("group-check self-test: x*(E^2-3E+3) passes, x*(E^2-3E+2) is caught", not group_check(_x * (_E**2 - 3 * _E + 3), 1)[2] and len(group_check(_x * (_E**2 - 3 * _E + 2), 1)[2]) == 1)

# ---------------------------------------------------------------- 2. G40 structure
q, p, j, c, e, b, a = c40.gens()
dg = G40.degrees()
check("G40 is quadratic in j", dg[2] == 2)
def coeff_in(P, idx, k, ctx):
    return ctx.from_dict({m: cf for m, cf in zip(P.monoms(), P.coeffs()) if m[idx] == k}) / (ctx.gens()[idx] ** k) if True else None
def part(P, idx, k, ctx):
    d = {}
    for m, cf in zip(P.monoms(), P.coeffs()):
        if m[idx] == k:
            mm = list(m); mm[idx] = 0; d[tuple(mm)] = cf
    return ctx.from_dict(d) if d else ctx.constant(0)
k0, k1, k2 = (part(G40, 2, i, c40) for i in range(3))
check("C (g40) does not involve j; C has 21 positive coefficients incl. d^3 (constant 1 at d=1)",
      C40.degrees()[2] == 0 and nt(C40) == 21 and nonneg(C40) and int(C40.coefficient(nt(C40) - 1)) == 1 and C40.monoms()[-1] == (0,) * 7)
G6, r6 = divmod(k0, C40)
check("k0 = C * G6 exactly", r6.is_zero())
check("G6 has 128 terms and is bilinear in e, p (degree <= 1 in each)", nt(G6) == 128 and G6.degrees()[1] <= 1 and G6.degrees()[4] <= 1, (nt(G6), G6.degrees()))
g00 = part(part(G6, 1, 0, c40), 4, 0, c40)
gE = part(part(G6, 1, 0, c40), 4, 1, c40)
gP = part(part(G6, 1, 1, c40), 4, 0, c40)
gEP = part(part(G6, 1, 1, c40), 4, 1, c40)
check("G6 = g00 + gE e + gP p + gEP e p", G6 == g00 + gE * e + gP * p + gEP * e * p)
check("g00 has positive coefficients and constant 1", nonneg(g00) and g00.monoms()[-1] == (0,) * 7 and int(g00.coeffs()[-1]) == 1, nt(g00))
k1e0 = k1.subs({"e": 0}) if False else part(k1, 4, 0, c40)
k2e0 = part(k2, 4, 0, c40)
check("k2|_{e=0} = 0 and k1|_{e=0} has nonnegative coefficients (doc: 234 terms)", k2e0.is_zero() and nonneg(k1e0), nt(k1e0))
info("term counts k0, k1, k2", (nt(k0), nt(k1), nt(k2)))

# ---------------------------------------------------------------- 3. the 16 strict-base tables
A, Tt, W, E, P_, Qv = sp.symbols("A T W E P Q", positive=True)
S40 = (A, Tt, W, E, P_, Qv)
x40 = flint.fmpz_mpoly_ctx.get(("A", "T", "W", "E", "P", "Q"), "lex")
X = dict(zip(["A", "T", "W", "E", "P", "Q"], x40.gens()))
# polynomial ring for (q,c,b,a) plus symbols for Le, Lp numerators/denominators
c4 = flint.fmpz_mpoly_ctx.get(("q", "c", "b", "a", "NE", "DE", "NP", "DP"), "lex")
def to4(P):  # from c40 (q,p,j,c,e,b,a) with p=j=e absent
    d = {}
    for m, cf in zip(P.monoms(), P.coeffs()):
        assert m[1] == 0 and m[2] == 0 and m[4] == 0
        d[(m[0], m[3], m[5], m[6], 0, 0, 0, 0)] = cf
    return c4.from_dict(d)
h00, hE, hP, hEP = map(to4, (g00, gE, gP, gEP))
NE_, DE_, NP_, DP_ = c4.gens()[4:]
F1 = h00
F2 = DE_ * hE + NE_ * h00
F3 = DP_ * hP + NP_ * h00
F4 = DE_ * DP_ * hEP + NP_ * DE_ * hE + NE_ * DP_ * hP + NE_ * NP_ * h00
bch = {"ab<1": Tt / (A * (1 + Tt)), "ab>1": (1 + Tt) / A}
cch = {"ac<1": W / (A * (1 + W)), "ac>1": (1 + W) / A}
sizes16 = []
allok = True
for bn, bx in bch.items():
    for cn, cx in cch.items():
        Le = sp.Integer(0) if bn == "ab<1" else (A * bx - 1) / ((A + 1) * (bx + 1))
        Lp = sp.Integer(0) if cn == "ac<1" else (A * cx - 1) / ((A + 1) * (cx + 1))
        mp = {"q": numden(Qv, x40, S40), "c": numden(cx, x40, S40), "b": numden(bx, x40, S40), "a": numden(A, x40, S40)}
        ne, de = numden(Le, x40, S40); npp, dpp = numden(Lp, x40, S40)
        if Le == 0: ne, de = x40.constant(0), x40.constant(1)
        if Lp == 0: npp, dpp = x40.constant(0), x40.constant(1)
        dens_ok = all(posden(mp[k][1]) for k in mp) and posden(de) and posden(dpp) and nonneg(ne) and nonneg(npp)
        row = []
        for Fk in (F1, F2, F3, F4):
            D = list(Fk.degrees())
            mapping = dict(mp); mapping.update({"NE": (ne, x40.constant(1)), "DE": (de, x40.constant(1)), "NP": (npp, x40.constant(1)), "DP": (dpp, x40.constant(1))})
            R = subst(Fk, ["q", "c", "b", "a", "NE", "DE", "NP", "DP"], mapping, x40, D)
            okk = (not R.is_zero()) and nonneg(R)
            allok &= okk and dens_ok
            row.append(nt(strip(R)))
        sizes16.append(((bn, cn), row, dens_ok))
for s in sizes16: info("strict-base table sizes (after removing positive content) " + str(s[0]), s[1])
check("16 strict-base tables: g00, gE+Le g00, gP+Lp g00, gEP+Lp gE+Le gP+Le Lp g00 have only nonnegative coefficients in all four product chambers (denominators positive)", allok)
# the algebra behind k0 >= C g00 (1-Le e)(1-Lp p): symbolic identity
ge, gp, gx, gy, gz, le, lp, ee, pp = sp.symbols("g00 gE gP gEP z le lp e p")
G6s = ge + gp * ee + gx * pp + gy * ee * pp
xx = gp + le * ge; yy = gx + lp * ge; zz = gy + lp * gp + le * gx + le * lp * ge
check("identity G6 - g00(1-Le e)(1-Lp p) = e x (1-Lp p) + p y (1-Le e) + e p z  (x,y,z the shifted functions)",
      sp.expand(G6s - ge * (1 - le * ee) * (1 - lp * pp) - (ee * xx * (1 - lp * pp) + pp * yy * (1 - le * ee) + ee * pp * zz)) == 0)
info("time after 16 tables", round(time.time() - T0, 1))

# ---------------------------------------------------------------- 4. the 24 interval tables (8 domains)
z = Tt / (1 + Tt)
def Pac(ax, cx): return (ax + 1) * (cx + 1) / (ax * cx - 1)
def Eab(ax, bx): return (ax + 1) * (bx + 1) / (ax * bx - 1)
rows = []
for rn, ax, bx, Jx in [
        ("(ab<1, a,b<=1)", (z**2 + A) / (1 + A), z**2 / ((z**2 + A) / (1 + A)), 1 + 1 / E),
        ("(ab<1, a>=1)", 1 + A, z**2 / (1 + A), 1 / E + (1 + A) * (1 - z)**2 / ((1 + A) - z)**2),
        ("(ab<1, b>=1)", z**2 / (1 + A), 1 + A, 1 / E + (1 + A) * (1 - z)**2 / ((1 + A) - z)**2)]:
    for cn in ("ac<1", "ac>1"):
        if cn == "ac<1": cx = W / (ax * (1 + W)); px = P_
        else: cx = (1 + W) / ax; px = Pac(ax, cx) * P_ / (1 + P_)
        rows.append((rn + ", " + cn, dict(a=ax, b=bx, c=cx, e=E, p=px, q=Qv), Jx))
ax = A; bx = (1 + Tt) / A
rows.append(("ab>1, ac<1", dict(a=ax, b=bx, c=W / (A * (1 + W)), e=Eab(ax, bx) * E / (1 + E), p=P_, q=Qv), 1 / (Eab(ax, bx) * E)))
# allocation domain
U_ = sp.symbols("U", positive=True)
S8 = (A, Tt, W, U_, E, Qv)
x8 = flint.fmpz_mpoly_ctx.get(("A", "T", "W", "U", "E", "Q"), "lex")
Dd = 1 + U_; Del = (1 + Tt) * Dd - 1; Hh = (1 + Tt) * Dd + A
row8 = dict(a=A, b=(1 + Tt) / A, c=(1 + W) / A, p=U_ * (A + 1) * ((1 + W) * Dd + A * U_) / (A * Dd * ((1 + W) * Dd - U_)),
            e=(A + 1) * Hh * E / (A * Dd * Del * (1 + E)), q=Qv)
J8 = A * Del / (Dd * (A + 1) * Hh * E)
# symbolic checks of the domain formulas
dl = U_ / (1 + U_); r0 = 1 - dl
aa, bb, cc = row8["a"], row8["b"], row8["c"]
Er0 = r0 * (aa + 1) * (bb + r0) / (aa * bb - r0)
check("domain 8: p = delta(a+1)(c+delta)/(ac-delta) with delta = U/(1+U)", sp.simplify(row8["p"] - dl * (aa + 1) * (cc + dl) / (aa * cc - dl)) == 0)
check("domain 8: e = E_r0 E/(1+E) and J = r0^2/e - r0^2/E_r0", sp.simplify(row8["e"] - Er0 * E / (1 + E)) == 0 and sp.simplify(J8 - (r0**2 / row8["e"] - r0**2 / Er0)) == 0)
check("domain 8: delta -> p is increasing on [0,1) from 0 to P_ac (ac>=1), so U>=0 covers p in [0,P_ac)",
      sp.simplify(sp.limit(dl * (aa + 1) * (cc + dl) / (aa * cc - dl), U_, sp.oo) - Pac(aa, cc)) == 0)
ax7 = A; bx7 = (1 + Tt) / A; e7 = Eab(ax7, bx7) * E / (1 + E)
check("domain 7: 1/(Eab E) = 1/e - 1/Eab", sp.simplify(1 / (Eab(ax7, bx7) * E) - (1 / e7 - 1 / Eab(ax7, bx7))) == 0)
rows.append(("ab>=1, ac>=1 (allocation)", row8, J8))

degs = [max(k0.degrees()[i], k1.degrees()[i], k2.degrees()[i]) for i in range(7)]
srcn = ["q", "p", "j", "c", "e", "b", "a"]
interval_sizes = []
sturm_total = 0; sturm_bad = []; ok_all = True
for rn, mp_, Jx in rows:
    ctx, syms = (x8, S8) if "allocation" in rn else (x40, S40)
    mp = {k: numden(v, ctx, syms) for k, v in mp_.items()}
    mp["j"] = (ctx.constant(0), ctx.constant(1))  # j does not occur in k0,k1,k2
    nJ, dJ = numden(Jx, ctx, syms)
    dens_ok = all(posden(mp[k][1]) for k in mp) and posden(dJ) and nonneg(nJ) and not nJ.is_zero()
    K0, K1, K2 = (subst(K, srcn, mp, ctx, degs) for K in (k0, k1, k2))
    base = K0; mid = 2 * dJ * K0 + nJ * K1; end = dJ * dJ * K0 + nJ * dJ * K1 + nJ * nJ * K2
    res = []
    for nm, Pp in (("base", base), ("middle", mid), ("endpoint", end)):
        Pp = strip(Pp)
        if nonneg(Pp) and not Pp.is_zero():
            res.append((nm, nt(Pp), "nonneg"))
        else:
            # sign check by grouping on E
            ei = ctx.names().index("E") if hasattr(ctx, "names") else 3
            ng, need, bad = group_check(Pp, ei)
            sturm_total += need; sturm_bad += [(rn, nm, x) for x in bad]
            res.append((nm, nt(Pp), f"{nneg(Pp)} negative coefficients; {need} E-groups need Sturm; {len(bad)} fail"))
            if bad: ok_all = False
    ok_all &= dens_ok
    interval_sizes.append((rn, res, dens_ok))
    info("interval tables " + rn, f"denominators positive: {dens_ok}; " + "; ".join(f"{x[0]} {x[1]} terms ({x[2]})" for x in res))
info("time after 24 tables", round(time.time() - T0, 1))
check(f"24 interval tables: nonnegative coefficients or exact Sturm positivity of every E-group on E>0 ({sturm_total} groups Sturm-checked)", ok_all and not sturm_bad,
      sturm_bad[:3])
# Bernstein identity
jj, JJ, K0s, K1s, K2s = sp.symbols("j J k0 k1 k2")
tt = jj / JJ
check("Bernstein identity G40 = (1-j/J)^2 k0 + (j/J)(1-j/J)(2k0+J k1) + (j/J)^2 G40(J)",
      sp.simplify((1 - tt)**2 * K0s + tt * (1 - tt) * (2 * K0s + JJ * K1s) + tt**2 * (K0s + JJ * K1s + JJ**2 * K2s) - (K0s + K1s * jj + K2s * jj**2)) == 0)

# ---------------------------------------------------------------- 5. G70: SOC peeling in the five KKT charts
pv, rv, av, bv, cv, qv, sv = c70.gens()
d70 = G70.degrees()
check("G70 is quadratic in s", d70[6] == 2)
L0, L1, L2 = (part(G70, 6, i, c70) for i in range(3))
check("L0, L1, L2 are (at most) quadratic in a", all(Li.degrees()[2] <= 2 for Li in (L0, L1, L2)), [Li.degrees()[2] for Li in (L0, L1, L2)])
B1 = L1 + 2 * av * L0
B2 = L2 + av * L1 + av * av * L0
L00 = part(L0, 2, 0, c70)
check("L0|_{a=0} has nonnegative coefficients and constant 1 (doc: 231 terms)", nonneg(L00) and L00.monoms()[-1] == (0,) * 7 and int(L00.coeffs()[-1]) == 1, nt(L00))
Ss, aS, Ls = sp.symbols("S a_ s_")
l0, l1, l2 = sp.symbols("l0 l1 l2")
Gs = lambda ss: l0 + ss * l1 + ss**2 * l2
check("identity (1+aS)^2 G70(a, S/(1+aS)) = L0 + S B1 + S^2 B2",
      sp.simplify((1 + aS * Ss)**2 * Gs(Ss / (1 + aS * Ss)) - (l0 + Ss * (l1 + 2 * aS * l0) + Ss**2 * (l2 + aS * l1 + aS**2 * l0))) == 0)
# charts
Tk, Lk, Uk, Wk, Yk, Rk = sp.symbols("T L U W Y R", positive=True)
S70 = (Tk, Lk, Uk, Wk, Yk, Rk)
x70 = flint.fmpz_mpoly_ctx.get(("T", "L", "U", "W", "Y", "R"), "lex")
th = Tk / (1 + Tk)
charts = {}
ph = 2 * th * Uk / (1 + Uk); qh = 2 * (1 - th) * Wk / (1 + Wk)
charts["HH"] = dict(b=th / Lk, c=(1 - th) / Lk, p=ph, q=qh, J=1 + Lk - (ph + qh) / 4, beta=th, d1=ph / 2, d2=qh / 2)
qc = 2 * (1 - th) * (1 + Wk)
charts["HC"] = dict(b=th / Lk, p=ph, q=qc, c=qc / (2 * Lk * (1 + Wk) + Wk), J=Lk + sp.Rational(1, 2) + th / (2 * (1 + Uk)), beta=th, d1=ph / 2, d2=1 - th)
qh2 = 2 * (1 - th) * Uk / (1 + Uk); pc = 2 * th * (1 + Wk)
charts["CH"] = dict(c=(1 - th) / Lk, q=qh2, p=pc, b=pc / (2 * Lk * (1 + Wk) + Wk), J=Lk + sp.Rational(1, 2) + (1 - th) / (2 * (1 + Uk)), beta=th, d1=th, d2=qh2 / 2)
pC = 2 * th * (1 + Uk); qC = 2 * (1 - th) * (1 + Uk + Wk)
charts["CC"] = dict(p=pC, q=qC, b=th / Lk, c=2 * (1 - th) * (1 + Uk) * (1 + Uk + Wk) / (2 * Lk * (1 + Uk) * (1 + Uk + Wk) + Wk), J=Lk + 1 / (2 * (1 + Uk)), beta=th, d1=th, d2=1 - th)
qC2 = 2 * (1 - th) * (1 + Uk); pC2 = 2 * th * (1 + Uk + Wk)
charts["CC'"] = dict(q=qC2, p=pC2, c=(1 - th) / Lk, b=2 * th * (1 + Uk) * (1 + Uk + Wk) / (2 * Lk * (1 + Uk) * (1 + Uk + Wk) + Wk), J=Lk + 1 / (2 * (1 + Uk)), beta=th, d1=th, d2=1 - th)
# energy, feasibility and stationarity of each chart point
bt = sp.symbols("beta")
for nm, ch in charts.items():
    b_, c_, p_, q_, J_ = ch["b"], ch["c"], ch["p"], ch["q"], ch["J"]
    be, d1, d2 = ch["beta"], ch["d1"], ch["d2"]
    energy = 1 - d1 - d2 + be**2 / b_ + (1 - be)**2 / c_ + d1**2 / p_ + d2**2 / q_
    ok_e = sp.simplify(energy - J_) == 0
    # feasibility: 0 <= d1 <= beta, 0 <= d2 <= 1-beta (each difference a positive-coefficient rational function)
    feas = []
    for xx_ in (d1, be - d1, d2, (1 - be) - d2):
        n_, dn_ = sp.fraction(sp.cancel(sp.together(xx_)))
        feas.append(n_ == 0 or (all(cf > 0 for cf in sp.Poly(n_, *S70).coeffs()) and all(cf > 0 for cf in sp.Poly(dn_, *S70).coeffs())) or
                    (all(cf < 0 for cf in sp.Poly(n_, *S70).coeffs()) and all(cf < 0 for cf in sp.Poly(dn_, *S70).coeffs())))
    # stationarity in beta of f(beta) with delta_i = min(p_i/2, beta_i) in the chart's pattern
    D1 = (p_ / 2) if nm in ("HH", "HC") else bt
    D2 = (q_ / 2) if nm in ("HH", "CH") else (1 - bt)
    f = 1 - D1 - D2 + bt**2 / b_ + (1 - bt)**2 / c_ + D1**2 / p_ + D2**2 / q_
    stat = sp.simplify(sp.diff(f, bt).subs(bt, be)) == 0
    check(f"KKT chart {nm}: energy at the chart point equals J; the point is feasible; d/dbeta = 0 there", ok_e and all(feas) and stat, (ok_e, feas, stat))

soc_ok = True
for nm, ch in charts.items():
    mp = {k: numden(ch[k], x70, S70) for k in ("b", "c", "p", "q")}
    nJ, dJ = numden(ch["J"], x70, S70)
    mp["a"] = (nJ * x70.gens()[4], dJ * (1 + x70.gens()[4]))
    mp["r"] = (x70.gens()[5], x70.constant(1))
    mp["s"] = (x70.constant(0), x70.constant(1))
    dens_ok = all(posden(mp[k][1]) for k in mp if k != "s") and nonneg(nJ) and posden(dJ)
    sizes = []
    for nmP, Pp in (("L0", L0), ("B1", B1), ("B2", B2)):
        Rp = subst(Pp, list(N70), mp, x70)
        okp = nonneg(Rp)
        soc_ok &= okp and dens_ok
        # split by powers of Y
        by = {}
        for m_, cf in zip(Rp.monoms(), Rp.coeffs()): by[m_[4]] = by.get(m_[4], 0) + 1
        sizes.append(f"{nmP}: {'nonneg' if okp else str(nneg(Rp)) + ' NEGATIVE'}; terms per Y-power {[by.get(k, 0) for k in range(max(by) + 1 if by else 0)]}")
    info(f"SOC tables chart {nm} (denominators positive: {dens_ok})", " | ".join(sizes))
check("60 SOC tables: cleared Y-coefficients of L0, B1, B2 nonnegative in all five KKT charts (a = J Y/(1+Y), r free)", soc_ok)
info("total time", round(time.time() - T0, 1))
print(f"SUMMARY: {len(FAIL)} failed checks", flush=True)
sys.exit(1 if FAIL else 0)
