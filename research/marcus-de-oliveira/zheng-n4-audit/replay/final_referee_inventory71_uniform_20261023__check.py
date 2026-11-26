"""Own-code replay for the Zheng n=4 audit, document final_referee_inventory71_uniform_20261023
(primary proof of paper case T_72 = {1,3,4,7,8,12}).
Everything is rebuilt from the paper's definitions through the audit library (lib/moc4lib.py).
Nothing from the supplement is read, imported or executed.
Exit code 1 if an EXACT check fails.  Runtime about 5 minutes (the 14x14 symbolic Hessian
determinant dominates).  Run through the job limiter."""
import sys, os, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lib'))
import moc4lib as L
import flint

FAIL = []
def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name, flush=True)
    if not cond:
        FAIL.append(name)

t0 = time.time()
FV = L.F_VALUES
IDX = {tuple(v): j for j, v in enumerate(FV)}
def gf(gi, j):
    return IDX[tuple(L.g_act_values(L.G_ELEMENTS[gi], list(FV[j])))]
def orbit(tup):
    return {tuple(gf(gi, j) for j in tup) for gi in range(1152)}

# ---------- A/B: contact set and transport ----------
T = set(L.CONTACT_SETS[72]); Tp = {0, 1, 2, 10, 14, 17}
g = (L.perm_from_oneline("1324"), L.perm_from_oneline("1342"), 0)
check("stated action rows (0,2,1,3), columns (0,2,3,1), no transpose maps T'={0,1,2,10,14,17} onto paper T_72",
      L.g_act_set(L.G_ELEMENTS.index(g), Tp) == T)
gen = L.generators(Tp)
check("generators(T') = [0,2,10,12,13,15,21,23,26] (doc's active nine-source list)", gen == [0, 2, 10, 12, 13, 15, 21, 23, 26])
check("all generators except f26 are minor; f26 = h4132,- non-minor",
      [L.F_IS_MINOR[j] for j in gen] == [True] * 8 + [False])
prof = [0, 2, 12, 13, 15, 21, 23, 26]
check("the eight profiles of (1) are linearly independent (rank 8 = dim of contact kernel 14-6)",
      L.rank([list(L.F_VECS[j]) for j in prof]) == 8 and L.rank([list(L.P[s]) for s in Tp]) == 6)
rows = [23, 15, 8, 7, 4, 16, 9]
check("identity rows: coefficients q,a,b,v,c,w,t are the gaps at [23,15,8,7,4,16,9]",
      all(FV[prof[i]][rows[k]] == (1 if i == k else 0) for k in range(7) for i in range(8)))
check("g11 = t + d, g20 = q + d, g12 = v + d (so d = g11 - g9)",
      [FV[j][11] for j in prof] == [0, 0, 0, 0, 0, 0, 1, 1] and [FV[j][20] for j in prof] == [1, 0, 0, 0, 0, 0, 0, 1]
      and [FV[j][12] for j in prof] == [0, 0, 0, 1, 0, 0, 0, 1])
check("(2) circuit f10 + f26 = f0 + f13 + f23 (values on all 24 rows)",
      all(FV[10][s] + FV[26][s] == FV[0][s] + FV[13][s] + FV[23][s] for s in range(24)))
check("f10 is a minor vector (so d<=0 side is a nonnegative minor sum)", L.F_IS_MINOR[10])
# atoms (foundations/so4-atoms.md): orbit membership with roles
soc = orbit((1, 9, 22)); H = orbit((3, 15, 26)); ch = orbit((2, 0, 12, 70))
check("SOC pairs (2,12;26),(2,15;26),(12,15;26) lie in the proved SOC orbit of (1,9;22)",
      all(((i, j, k) in soc or (j, i, k) in soc) for (i, j, k) in [(2, 12, 26), (2, 15, 26), (12, 15, 26)]))
check("old-H (aux f23 = t, main f15 = c; f26) lies in the proved old-H orbit, with these roles", (23, 15, 26) in H)
check("old-H reversed roles (15,23;26) NOT in orbit (roles as printed matter)", (15, 23, 26) not in H)
check("chain (centre f15; leaves f21, f23; f26) lies in the proved chain orbit", (15, 21, 23, 26) in ch or (15, 23, 21, 26) in ch)

# ---------- cubic, height linearity, Hessian determinant, G9 ----------
names = ('q', 'a', 'b', 'v', 'c', 'w', 't', 'z')
ctx = flint.fmpz_mpoly_ctx.get(names, 'degrevlex')
q, a, b, v, c, w, t, z = ctx.gens(); d = ctx.constant(1)
coef = [q, a, b, v, c, w, t, d]
m = [ctx.constant(0)] * 14
for j, cf in zip(prof, coef):
    for k in range(14):
        assert L.F_VECS[j][k].denominator == 1
        m[k] = m[k] + int(L.F_VECS[j][k]) * cf
m[13] = m[13] + z     # kappa = e_13
mctx = flint.fmpz_mpoly_ctx.get(tuple('m%d' % i for i in range(14)), 'degrevlex')
Cm = L.cubic_C(list(mctx.gens()))
Cst = (q + d) * (v + d) * (w + t + d) + a * b * c
Qst = a * b + a * c + b * c + d * d - t * (c - d)
check("(8) C(g + z kappa) = C + zQ with C=(q+d)(v+d)(w+t+d)+abc, Q=ab+ac+bc+d^2-t(c-d) (d=1; homogeneous)",
      Cm.compose(*m, ctx=ctx) == Cst + z * Qst)
Hm = [[Cm.derivative(i).derivative(j).compose(*m, ctx=ctx) for j in range(14)] for i in range(14)]
A = [r[:] for r in Hm]; n = 14; sign = 1; prev = ctx.constant(1)
for k in range(n - 1):
    if A[k][k].is_zero():
        cand = [i for i in range(k + 1, n) if not A[i][k].is_zero()]
        i = min(cand, key=lambda i: len(A[i][k].coeffs()))
        A[k], A[i] = A[i], A[k]; sign = -sign
    for i in range(k + 1, n):
        for j in range(k + 1, n):
            qq, rr = divmod(A[k][k] * A[i][j] - A[i][k] * A[k][j], prev)
            assert rr.is_zero()
            A[i][j] = qq
        A[i][k] = ctx.constant(0)
    prev = A[k][k]
detH = sign * A[n - 1][n - 1]
print("info det Hess C(g+z kappa): terms", len(detH.coeffs()), "z-degree", detH.degrees()[7], "time", round(time.time() - t0, 1), flush=True)
dd = [{} for _ in range(detH.degrees()[7] + 1)]
for mon, cf in zip(detH.monoms(), detH.coeffs()):
    e = list(mon); kz = e[7]; e[7] = 0; dd[kz][tuple(e)] = cf
Pk = [ctx.from_dict(x) for x in dd]
LHS = sum((Pk[k] * (-Cst) ** k * Qst ** (4 - k) for k in range(len(Pk))), ctx.constant(0))
G9q, rem = divmod(LHS, -4 * Qst * Qst)
check("Q^4 det Hess C(g - C kappa/Q) is divisible by -4Q^2", rem.is_zero())
G9 = G9q.sqrt()
if G9.to_dict().get((0,) * 8, 0) < 0:
    G9 = -G9
check("(9) Q^4 det Hess C(g - C kappa/Q) = -4 Q^2 G9^2 exactly, G9 polynomial", -4 * Qst * Qst * G9 * G9 == LHS)
check("G9: 1615 terms, degree 9, [d^9]G9 = 1, degrees in (q,a,b,v,c,w,t) = (3,2,2,3,2,3,2)",
      len(G9.coeffs()) == 1615 and G9.total_degree() == 9 and G9.to_dict().get((0,) * 8) == 1
      and G9.degrees()[:7] == (3, 2, 2, 3, 2, 3, 2))
check("(10) G9(q,a,b,v,c,w,t) = G9(v,b,a,q,c,w,t)", G9.compose(v, b, a, q, c, w, t, z, ctx=ctx) == G9)
zero = ctx.constant(0)
GK = G9.compose(q, zero, zero, v, c, w, zero, z, ctx=ctx)
check("face a=b=t=0: G9 = ((q+1)(v+1)+c(q+v+1))(q+1)^2(v+1)^2(w+1)^3 (used by T_74 document, eq. (17))",
      GK == ((q + 1) * (v + 1) + c * (q + v + 1)) * (q + 1) ** 2 * (v + 1) ** 2 * (w + 1) ** 3)

# ---------- the 17 positivity tables ----------
S = flint.fmpz_mpoly_ctx.get(('Q0', 'V', 'X', 'Y', 'Z', 'T', 'W'), 'degrevlex')
SN = ('Q0', 'V', 'X', 'Y', 'Z', 'T', 'W')
Q0, V, X, Y, Z, TT, W = S.gens(); one = S.constant(1); szero = S.constant(0)
groups = {}
for mon, cf in zip(G9.monoms(), G9.coeffs()):
    groups.setdefault((mon[1], mon[2], mon[4], mon[5], mon[6]), {})[(mon[0], mon[3])] = int(cf)
MM = {'a': 2, 'b': 2, 'c': 2, 'w': 3, 't': 2}
def clear(sub):
    pw = {}
    for k in MM:
        N, D = sub[k]; pw[k] = ([one], [one])
        for i in range(MM[k]):
            pw[k][0].append(pw[k][0][-1] * N); pw[k][1].append(pw[k][1][-1] * D)
    tot = szero
    for key, cd in groups.items():
        term = sum((cf * Q0 ** e0 * V ** e3 for (e0, e3), cf in cd.items()), szero)
        for k, e in zip(['a', 'b', 'c', 'w', 't'], key):
            term = term * pw[k][0][e] * pw[k][1][MM[k] - e]
        tot = tot + term
    return tot
def coef_of(P, mon):
    e = [0] * 7
    for k, val in mon.items():
        e[SN.index(k)] = val
    return P.to_dict().get(tuple(e), 0)
Lx = 1 + X; D11 = 1 + X + Y + Z
charts = [
    ('(11)', {'a': (Y + Z, D11), 'b': (X + Z, D11), 'c': (X + Y, D11)}, (28693, 10012, 93570),
     ([({'X': 2}, 1), ({'Y': 2}, 1)], [({}, 1)], [({'T': 3, 'X': 3}, 1), ({'T': 3, 'Y': 3}, 1)])),
    ('(12)', {'a': (1 + Y * Lx * Lx, Lx * (1 + Y)), 'b': (Z, Lx * (1 + Z)), 'c': (one, Lx * (1 + Z))}, (27186, 10665, 70122),
     ([({'X': 1}, 8)], [({}, 8)], [({'T': 3}, 8)])),
    ('(13)', {'a': (Z, Lx * (1 + Y) * (1 + Z)), 'b': (one, Lx * (1 + Y) * (1 + Z)), 'c': (one, Lx)}, (36773, 13175, 70650),
     ([({'X': 1}, 8), ({'Y': 1}, 4)], [({}, 8)], [({'T': 3}, 8)])),
    ('(14)', {'a': (Y * Z, Lx * (1 + Y) * (1 + Z)), 'b': (Y, Lx * (1 + Y) * (1 + Z)), 'c': (Lx, one)}, (38081, 12700, 79550),
     ([({}, 2)], [({}, 10)], [({'T': 5}, 3)])),
    ('solo', {'a': (szero, one), 'b': (szero, one), 'c': (one, Lx)}, (703, 271, 1674),
     ([({}, 2)], [({}, 10)], [({'T': 3}, 10)])),
]
total = 0
for nm, sub, counts, strict in charts:
    Nc, Dc = sub['c']
    rowsub = [((TT * Dc, Nc * (1 + TT)), (W, one)), ((Dc, Nc), (W, one))]
    if nm != '(14)':
        rowsub.append((((1 + TT) * Dc, Nc), ((1 + TT) * Dc * W, Nc * TT * (1 + W))))
    else:   # (16)
        rowsub.append(((X * (1 + TT) + TT, X * (1 + X) * (1 + TT)), ((X * (1 + TT) + TT) * W, (1 + X) * TT * (1 + W))))
    for rn, (ts, ws), cnt, st in zip(('ct<1', 'ct=1', 'ct>1'), rowsub, counts, strict):
        s = dict(sub); s['t'] = ts; s['w'] = ws
        P = clear(s); cs = P.coeffs(); total += len(cs)
        check(f"table {nm} {rn}: {len(cs)} terms (doc {cnt}), all coefficients > 0, strict marker {st}",
              len(cs) == cnt and all(x > 0 for x in cs) and all(coef_of(P, mo) == val for mo, val in st))
for nm, sub, cnt, cst in [('(19a)', {'a': (X, 1 + X), 'b': (Y, 1 + Y)}, 1293, 1),
                          ('(19b)', {'a': (1 + X, one), 'b': (Y, (1 + X) * (1 + Y))}, 1680, 2)]:
    s = dict(sub); s['c'] = (szero, one); s['w'] = (W, one); s['t'] = (TT, one)
    P = clear(s); cs = P.coeffs(); total += len(cs)
    check(f"table c=0 {nm}: {len(cs)} terms (doc {cnt}), all coefficients > 0, constant term {cst}",
          len(cs) == cnt and all(x > 0 for x in cs) and coef_of(P, {}) == cst)
check("total number of positive coefficients = 496798", total == 496798)

# ---------- NUMERIC (brief) ----------
import numpy as np
print("NUMERIC section", flush=True)
Gd = {tuple(int(x) for x in e): int(co) for e, co in G9.to_dict().items()}
def G9num(cf):
    qv, av, bv, vv, cv, wv, tv = cf
    val = 0.0; scale = 0.0
    for e, co in Gd.items():
        tm = float(co) * qv**e[0] * av**e[1] * bv**e[2] * vv**e[3] * cv**e[4] * wv**e[5] * tv**e[6]
        val += tm; scale += abs(tm)
    return val, scale
def H3(a_, b_, c_):
    mx = max(a_, b_, c_); s_ = a_ + b_ + c_
    return s_ / 2 if 2 * mx <= s_ else (mx * (s_ - mx)) ** 0.5
worst_rel = 0.0; inneg = 0; found = 0; minphi = np.inf
for seed in range(8):
    pt = L.random_point_on_contact_stratum(sorted(Tp), seed=seed, tries=30, margin=1e-2)
    if pt is None:
        continue
    found += 1
    gaps = pt['gaps']; cf = [gaps[r_] for r_ in rows]; dv = gaps[11] - gaps[9]
    if dv <= 1e-9:
        continue
    nc = [x / dv for x in cf]
    val, scale = G9num(nc); worst_rel = max(worst_rel, abs(val) / scale)
    qv, av, bv, vv, cv, wv, tv = nc
    neg = H3(av, bv, cv) < 1 and tv * (cv - 1) < 1 and (cv * tv <= 1 or wv < tv / (cv * tv - 1))
    inneg += neg
    c_vec = np.linalg.lstsq(np.array(L.P_numpy(), float), np.array(pt['gaps'], float), rcond=None)[0]
    mv, _ = L.min_phi_over_SO4(c_vec, restarts=8, seed=seed)
    minphi = min(minphi, mv)
print(f"NUMERIC T' stratum points found {found}/8 (seeds 0-7, margin 1e-2): max |G9|/sum|terms| at d>0 points = {worst_rel:.2e} "
      f"(G9 must vanish at every source: radial necessity); points in the negative region: {inneg}; min Phi over SO(4) (8 restarts) = {minphi:.4g}")
print("time", round(time.time() - t0, 1))
print("SUMMARY", "all EXACT checks passed" if not FAIL else f"{len(FAIL)} FAILED: {FAIL}")
sys.exit(1 if FAIL else 0)
