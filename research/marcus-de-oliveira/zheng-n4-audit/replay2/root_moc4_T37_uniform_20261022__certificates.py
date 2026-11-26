# Certificate replay for root_moc4_T37_uniform_20261022 (paper case T_38), auditor label cert-T37.
# Own code only; nothing from the supplement is executed or imported. The supplement's radial
# JSON is loaded as DATA (only to compare the rebuilt G9 with the stored factor).
# Usage: MOC4_SUPP=<supplement root> python -I -B root_moc4_T37_uniform_20261022__certificates.py
#        optional env G9_CACHE=<file> : cache the rebuilt G9 (it is always rebuilt if the file is absent)
# Prints one line per check; exits nonzero if any EXACT check fails.
import os, sys, json, time, random
from fractions import Fraction as Fr
_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, "..", "lib"))
from moc4lib import F_VECS, F_VALUES, CONTACT_SETS, cubic_C, matching_vector, values
from flint import fmpz_mpoly_ctx, fmpz, fmpq

FAIL = []
T0 = time.time()
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + ((" : " + str(info)) if info != "" else ""), flush=True)
    if not ok: FAIL.append(name)
def info(name, msg): print("INFO " + name + " : " + str(msg) + f"  [{time.time()-T0:.1f}s]", flush=True)

# ---------------------------------------------------------------- 0. sanity: X in Sing C (paper App. A)
ms14 = fmpz_mpoly_ctx.get([f"m{i}" for i in range(14)])
M = ms14.gens()
Cm = cubic_C(list(M))
grad = [Cm.derivative(i) for i in range(14)]
rng = random.Random(3807)
ok = True
for _ in range(5):
    a = [Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(4)]
    b = [Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(4)]
    m = matching_vector(a, b)
    mq = [fmpq(x.numerator, x.denominator) for x in m]
    for gp in grad:
        val = fmpq(0)
        for ex, co in gp.to_dict().items():
            t = fmpq(int(co))
            for i, k in enumerate(ex):
                if k: t *= mq[i] ** k
            val += t
        ok &= (val == 0)
check("grad C vanishes at 5 random rational matching vectors (X in Sing C)", ok)

# ---------------------------------------------------------------- 1. rebuild G9 from (A.2)
GEN = [1, 2, 15, 9, 10, 56, 0, 53]          # q p j c e d b a   (document eq. (1))
NAMES = ["q", "p", "j", "c", "e", "d", "b", "a"]
check("contact set of the document == paper T_38", CONTACT_SETS[38] == frozenset({1, 2, 3, 8, 11, 12}))
check("generators vanish on T_38", all(F_VALUES[g][s] == 0 for g in GEN for s in CONTACT_SETS[38]))
# integer scaling of the source coordinates
den = 1
for g in GEN:
    for x in F_VECS[g]:
        den = den * Fr(x).denominator // __import__("math").gcd(den, Fr(x).denominator)
info("common denominator of generator coordinates", den)
R9 = fmpz_mpoly_ctx.get(NAMES + ["h"])
V9 = R9.gens(); hh = V9[8]
w = [R9.constant(0)] * 14
for s_, g in zip(V9[:8], GEN):
    for r in range(14):
        x = Fr(F_VECS[g][r]) * den
        assert x.denominator == 1
        if x: w[r] = w[r] + int(x) * s_
w[13] = w[13] + hh                      # w = den*g + h*kappa  (kappa = e_13)
Cw = cubic_C(w)
Cd = Cw.to_dict()
check("C(den*g + h kappa) is affine in h", all(ex[8] <= 1 for ex in Cd))
C0 = R9.from_dict({ex: co for ex, co in Cd.items() if ex[8] == 0})
Q1 = R9.from_dict({ex[:8] + (0,): co for ex, co in Cd.items() if ex[8] == 1})
info("C0 terms / Q terms", (len(C0), len(Q1)))
# Hessian: constant third-derivative tensor applied to w
Hm = [[Cm.derivative(i).derivative(j) for j in range(14)] for i in range(14)]
def lin(poly):  # linear form in m -> polynomial in R9 via w
    out = R9.constant(0)
    for ex, co in poly.to_dict().items():
        k = ex.index(1)
        out = out + int(co) * w[k]
    return out
H = [[lin(Hm[i][j]) for j in range(14)] for i in range(14)]

def bareiss_det(A):
    A = [row[:] for row in A]; n = len(A); sg = 1; prev = None
    for k in range(n - 1):
        cand = [(len(A[i][j]), i, j) for i in range(k, n) for j in range(k, n) if not A[i][j].is_zero()]
        if not cand: return R9.constant(0)
        _, ri, cj = min(cand)
        if ri != k: A[k], A[ri] = A[ri], A[k]; sg = -sg
        if cj != k:
            for row in A: row[k], row[cj] = row[cj], row[k]
            sg = -sg
        pv = A[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                x = pv * A[i][j] - A[i][k] * A[k][j]
                if prev is not None:
                    qq, rr = divmod(x, prev)
                    assert rr.is_zero()
                    x = qq
                A[i][j] = x
            A[i][k] = R9.constant(0)
        prev = pv
    return sg * A[n - 1][n - 1]

R8 = fmpz_mpoly_ctx.get(NAMES)
def to8(p):  # drop the h variable (must not occur)
    return R8.from_dict({ex[:8]: co for ex, co in p.to_dict().items()})
cache = os.environ.get("G9_CACHE")
if cache and os.path.exists(cache):
    G9 = R8.from_dict({tuple(e): int(c) for e, c in json.load(open(cache))})
    info("G9 loaded from own cache (built by this script)", len(G9))
else:
    det = bareiss_det(H)
    info("det Hess C(den*g + h kappa): terms, h-degree", (len(det), det.degrees()[8]))
    check("det Hess is of degree <= 4 in h", det.degrees()[8] <= 4)
    # E = Q^4 det(H at h = -C0/Q)
    E = R9.constant(0)
    dd = det.to_dict()
    for k in range(det.degrees()[8] + 1):
        pk = R9.from_dict({ex[:8] + (0,): co for ex, co in dd.items() if ex[8] == k})
        E = E + pk * (-C0) ** k * Q1 ** (4 - k)
    qq, rr = divmod(E, Q1 ** 2)
    check("Q^2 divides Q^4 det Hess(den g - (C/Q) kappa)", rr.is_zero())
    cont = qq.content()
    lead = qq.to_dict().get((0, 0, 0, 0, 0, 18, 0, 0, 0))
    sgn = -1 if lead is not None and lead < 0 else 1
    base = R9.from_dict({ex: co // (sgn * cont) for ex, co in qq.to_dict().items()})
    G = base.sqrt()
    check("E / Q^2 = (negative constant) * (integer square)", G is not None and G * G == base and sgn == -1,
          f"constant = {sgn * cont}")
    G9 = to8(G)
    if G9.to_dict().get((0, 0, 0, 0, 0, 9, 0, 0), 0) < 0: G9 = -G9
    if cache:
        json.dump([[list(map(int, e)), int(c)] for e, c in G9.to_dict().items()], open(cache, "w"))
d9 = G9.to_dict().get((0, 0, 0, 0, 0, 9, 0, 0), 0)
check("G9 rebuilt: 1154 terms, d^9 coefficient +1, homogeneous of degree 9",
      len(G9) == 1154 and d9 == 1 and all(sum(e) == 9 for e in G9.to_dict()), f"{len(G9)} terms, d^9 coeff {d9}")
# compare with the supplement's stored factor (DATA only)
supp = os.environ.get("MOC4_SUPP") or (sys.argv[1] if len(sys.argv) > 1 else None)
if supp:
    rad = json.load(open(os.path.join(supp, "proofs", "manin_transfer_20261021", "root_moc4_T37_full_radial_20261022.json")))
    assert rad["variables"][:8] == NAMES
    st = next(f["terms"] for f in rad["radial_factors"] if f["degree"] == 9)
    Gs = R8.from_dict({tuple(ex[:8]): -int(co) for ex, co in st})
    check("rebuilt G9 equals the negative of the supplement's stored degree-9 factor", Gs == G9)
    st2 = next(f["terms"] for f in rad["radial_factors"] if f["degree"] == 2)
    info("supplement degree-2 factor terms (expected to be +-Q up to coordinates)", len(st2))

# ---------------------------------------------------------------- 2. dehomogenize d = 1
R7 = fmpz_mpoly_ctx.get(["q", "p", "j", "c", "e", "b", "a"])
gd = {}
for ex, co in G9.to_dict().items():
    k = (ex[0], ex[1], ex[2], ex[3], ex[4], ex[6], ex[7])
    gd[k] = gd.get(k, 0) + int(co)
G = R7.from_dict(gd)
DEG = G.degrees()
info("G (d=1): terms, degrees in q p j c e b a", (len(G), DEG))

# ---------------------------------------------------------------- 3. chart substitutions
CV = ["T", "J", "U", "A", "C", "B", "Q"]
O = fmpz_mpoly_ctx.get(CV)
T, J, U, A, Cc, B, Q = O.gens()
one = O.constant(1)
def clear(vals):
    """vals: list of (num, den) for q p j c e b a; returns sum co prod num^k den^(deg-k)."""
    pw = []
    for (n, dn), m in zip(vals, DEG):
        pn = [one]; pd = [one]
        for _ in range(m): pn.append(pn[-1] * n); pd.append(pd[-1] * dn)
        pw.append((pn, pd, m))
    out = O.constant(0)
    for ex, co in G.to_dict().items():
        t = O.constant(int(co))
        for i, k in enumerate(ex):
            pn, pd, m = pw[i]
            t = t * pn[k] * pd[m - k]
        out = out + t
    return out
def pos_on_domain(poly, strict):
    """nonnegative coefficients and some positive monomial using only strict variables -> > 0 on domain"""
    d = poly.to_dict()
    return all(c > 0 for c in d.values()) and any(all(ex[i] == 0 for i in range(7) if i not in strict) for ex in d)

def absorb(poly, direction):
    """Absorb each negative coefficient m at x^ex by a binomial square
       L x^(ex-dir) + m x^ex + R x^(ex+dir) >= 0  (valid for nonnegative variables iff m^2 <= 4 L R).
       Own allocation: take L = R scaled from available mass. Returns (remainder dict, list) or None."""
    rem = {ex: Fr(int(c)) for ex, c in poly.to_dict().items()}
    used = []
    negs = sorted([ex for ex, c in rem.items() if c < 0], key=lambda e: rem[e])
    for ex in negs:
        m = -rem[ex]
        lo = tuple(x - y for x, y in zip(ex, direction)); hi = tuple(x + y for x, y in zip(ex, direction))
        if min(lo) < 0: return None, f"lower neighbour of {ex} not a monomial"
        L0 = rem.get(lo, Fr(0)); R0 = rem.get(hi, Fr(0))
        if L0 <= 0 or R0 <= 0 or 4 * L0 * R0 < m * m: return None, f"insufficient mass at {ex}: m={m}, L={L0}, R={R0}"
        # choose L = m^2/(4 R0) * s with s = 1 if possible ... take the balanced split l = m/2 * sqrt(L0/R0) approx:
        # use the smallest L with L >= m^2/(4 R0) (uses all of R0) unless that exceeds L0; then rescale.
        L = m * m / (4 * R0); Rr = R0
        if L > L0: return None, "allocation failed"
        rem[lo] -= L; rem[hi] -= Rr; rem[ex] = Fr(0)
        used.append((ex, str(L), str(-m), str(Rr)))
    return rem, used

CHARTS = []
def chart(name, vals, strict, expect_terms, direction=None, expect_neg=0):
    CHARTS.append((name, vals, strict, expect_terms, direction, expect_neg))
def fr(n, d=one): return (n, d)
# --- both small: p = T/(1+T), j = J/(1+J)
pS = fr(T, 1 + T); jS = fr(J, 1 + J)
chart("both_small_u<1", [fr(Q), pS, jS, fr(Cc), fr(U * (1 + J), (1 + U) * J), fr(B), fr(A)], {1, 2, 3, 6}, 3186, (0, 0, 0, 0, 0, 1, 0), 2)
chart("both_small_u=1", [fr(Q), pS, jS, fr(Cc), fr(1 + J, J), fr(B), fr(A)], {1, 2, 3, 6}, 1309, (0, 0, 0, 0, 0, 1, 0), 2)
chart("both_small_u>1_b>=j", [fr(Q), pS, jS, fr(Cc), fr((1 + U) * (1 + J), J),
      fr(J * J * U + (1 + U) * (1 + J) ** 2 * B, J * U * (1 + J) * (1 + B)), fr(A)], {1, 2, 3, 6}, 12522)
# --- j > 1: j = 1+J, p = T/((1+T)(1+J))
pH = fr(T, (1 + T) * (1 + J)); jH = fr(1 + J)
chart("j>1_u<1", [fr(Q), pH, jH, fr(Cc), fr(U, (1 + U) * (1 + J)), fr(B), fr((1 + U + J) * A, J * (1 + U) * (1 + J) * (1 + A))], {1, 2, 3, 6}, 12416)
chart("j>1_u=1", [fr(Q), pH, jH, fr(Cc), fr(one, 1 + J), fr(B), fr(A, J * (1 + J) * (1 + A))], {1, 2, 3, 6}, 2564)
D_ = 1 + U + A
chart("j>1_u>1_simplex", [fr(Q), pH, jH, fr(Cc), fr(J * D_ + U, J * D_ * (1 + J)), fr((J * D_ + U) * B, (1 + J) * U * (1 + B)),
      fr(A, J * D_ * (1 + J))], {1, 2, 3, 6}, 31706, (1, 0, 0, 0, 0, 1, 0), 20)
# --- p > 1: p = 1+J, j = T/((1+T)(1+J)), c = C/(J(1+C))
pP = fr(1 + J); jP = fr(T, (1 + T) * (1 + J)); cP = fr(Cc, J * (1 + Cc))
chart("p>1_u<1", [fr(Q), pP, jP, cP, fr(U * (1 + T) * (1 + J), (1 + U) * T), fr(B), fr(A)], {0, 1, 2, 3, 6}, 6354)
chart("p>1_u=1", [fr(Q), pP, jP, cP, fr((1 + T) * (1 + J), T), fr(B), fr(A)], {0, 1, 2, 3, 6}, 2943)
chart("p>1_u>1", [fr(Q), pP, jP, cP, fr((1 + U) * (1 + T) * (1 + J), T), fr((1 + U) * (1 + T) * (1 + J) * B, T * U * (1 + B)), fr(A)], {0, 1, 2, 3, 6}, 23198)
# --- p = 1, j = J/(1+J)
chart("p=1_u<1", [fr(Q), fr(one), jS, fr(Cc), fr(U * (1 + J), (1 + U) * J), fr(B), fr(A)], {1, 2, 3, 6}, 1077)
chart("p=1_u=1", [fr(Q), fr(one), jS, fr(Cc), fr(1 + J, J), fr(B), fr(A)], {1, 2, 3, 6}, 450)
chart("p=1_u>1", [fr(Q), fr(one), jS, fr(Cc), fr((1 + U) * (1 + J), J), fr((1 + U) * (1 + J) * B, J * U * (1 + B)), fr(A)], {1, 2, 3, 6}, 2843)
# --- j = 1, p = T/(1+T)
chart("j=1_u<1", [fr(Q), pS, fr(one), fr(Cc), fr(U, 1 + U), fr(B), fr(A)], {2, 3, 6}, 833)
chart("j=1_u=1", [fr(Q), pS, fr(one), fr(Cc), fr(one), fr(B), fr(A)], {2, 3, 6}, 296)
chart("j=1_u>1_b>=1", [fr(Q), pS, fr(one), fr(Cc), fr(1 + U), fr(U * (1 + B) + B, U * (1 + B)), fr(A)], {2, 3, 6}, 1359)

results = {}
for name, vals, strict, exp_terms, direction, exp_neg in CHARTS:
    # denominators must be positive on the chart domain
    dens_ok = all(pos_on_domain(dn, strict) for _, dn in vals)
    P_ = clear(vals)
    nneg = sum(1 for c in P_.to_dict().values() if c < 0)
    line = f"{len(P_)} terms (document: {exp_terms}), {nneg} negative coefficients (document: {exp_neg})"
    if nneg == 0:
        rem = {ex: Fr(int(c)) for ex, c in P_.to_dict().items()}; used = []
    else:
        rem, used = absorb(P_, direction)
    if rem is None:
        check(f"chart {name}: certificate", False, line + " ; absorption failed: " + str(used)); continue
    nonneg = all(c >= 0 for c in rem.values())
    wit = [ex for ex, c in rem.items() if c > 0 and all(ex[i] == 0 for i in range(7) if i not in strict)]
    w0 = min(wit, key=lambda e: (sum(e), e)) if wit else None
    results[name] = (len(P_), nneg)
    check(f"chart {name}: denominators > 0 on domain; cleared G9 has nonneg. coefficients after {len(used)} binomial squares; strict witness {w0} (coeff {rem.get(w0) if w0 else None})",
          dens_ok and nonneg and w0 is not None and len(P_) == exp_terms and nneg == exp_neg, line)
    if used and name.startswith("both_small"):
        info(f"chart {name} squares (center, L, middle, R)", used)
    if name == "j>1_u>1_simplex":
        info(f"chart {name} first squares", used[:3])
        t = (0, 1, 2, 1, 0, 0, 0)  # 3 J U^4 A^2 claimed:  exponent (T,J,U,A,C,B,Q) = (0,1,4,2,0,0,0)
        info("coefficient of J U^4 A^2 in the remainder", rem.get((0, 1, 4, 2, 0, 0, 0)))

# ---------------------------------------------------------------- 4. lower-b high-u0 chart and its j=1 layer
Zv = [fr(Q), pS, jS, fr(Cc), fr((1 + U) * (1 + J), J), fr(J * B, (1 + J) * (1 + B)), fr(A)]
Z = clear(Zv)
check("lower-b chart: denominators > 0 on domain", all(pos_on_domain(dn, {1, 2, 6}) for _, dn in Zv))
Zd = Z.to_dict()
Z0 = O.from_dict({ex: c for ex, c in Zd.items() if ex[3] == 0})
ZA = Z - Z0
check("lower-b chart Z: 7134 terms; A-part ZA has nonneg. coefficients; Z0 has 2330 terms",
      len(Z) == 7134 and all(c >= 0 for c in ZA.to_dict().values()) and len(Z0) == 2330, (len(Z), len(Z0), sum(c < 0 for c in Zd.values())))
unit, facs = Z0.factor()
info("factorization of Z0 (total degree, multiplicity, terms)", (unit, [(f.total_degree(), int(m), len(f)) for f, m in facs]))
P5 = next(f for f, m in facs if f.total_degree() == 5)
F11 = next(f for f, m in facs if f.total_degree() == 11)
if P5.to_dict().get((0, 2, 0, 0, 0, 0, 0), 0) < 0: P5 = -P5; F11 = -F11
sq = (J + 1) ** 2 * (T * J * B + T * B + T - J) ** 2
R11 = F11 - U * Cc * Q * sq
check("Z = ZA + (J+1) P5 (R11 + U C Q (J+1)^2 (TJB+TB+T-J)^2) exactly", Z == ZA + (J + 1) * P5 * (R11 + U * Cc * Q * sq))
check("P5 (40 terms) and R11 (309 terms) have positive coefficients; P5 has 2J^2, R11 has 5J^4",
      len(P5) == 40 and len(R11) == 309 and all(c > 0 for c in P5.to_dict().values()) and all(c > 0 for c in R11.to_dict().values())
      and P5.to_dict().get((0, 2, 0, 0, 0, 0, 0)) == 2 and R11.to_dict().get((0, 4, 0, 0, 0, 0, 0)) == 5)
# j = 1 lower chart, computed directly
Z1 = clear([fr(Q), pS, fr(one), fr(Cc), fr(1 + U), fr(B, 1 + B), fr(A)])
def coeffJ(p, k): return O.from_dict({(e[0], 0) + e[2:]: c for e, c in p.to_dict().items() if e[1] == k})
check("deg_J Z = 7 and Z1 (direct substitution j=1, e=1+U, b=B/(1+B)) equals [J^7] Z", Z.degrees()[1] == 7 and coeffJ(Z, 7) == Z1, f"{len(Z1)} terms")
Al, Pl, Rl = coeffJ(ZA, 7), coeffJ(P5, 2), coeffJ(R11, 4)
check("Z1 = Alead + Plead (Rlead + U C Q (TB-1)^2), factors nonneg., constants 2 and 5",
      Z1 == Al + Pl * (Rl + U * Cc * Q * (T * B - 1) ** 2) and all(c >= 0 for p in (Al, Pl, Rl) for c in p.to_dict().values())
      and Pl.to_dict().get((0,) * 7) == 2 and Rl.to_dict().get((0,) * 7) == 5)
check("P5 and R11 have J-degree exactly 2 and 4 (so the leading-coefficient extraction is the J^7 coefficient)",
      P5.degrees()[1] == 2 and R11.degrees()[1] == 4 and ZA.degrees()[1] <= 7)

# ---------------------------------------------------------------- 5. covering: exact grid test
def hard(q, p, j, c, e, b, a):
    if p * j >= 1: return False
    if c > 0 and p >= 1 + 1 / c: return False
    if j >= 1 + 1 / (a + e): return False
    if b > 0 and j >= 1 / e + 1 / b: return False
    return True
def ev(nd, x):
    n, dn = nd
    def E(p):
        s = Fr(0)
        for ex, co in p.to_dict().items():
            t = Fr(int(co))
            for i, k in enumerate(ex):
                if k: t *= x[i] ** int(k)
            s += t
        return s
    return E(n) / E(dn)
def locate(q, p, j, c, e, b, a):
    """invert the chart maps; return (chart name, chart point (T,J,U,A,C,B,Q))"""
    u = j * e
    if p < 1 and j < 1:
        Tt = p / (1 - p); Jj = j / (1 - j)
        if u < 1: return "both_small_u<1", (Tt, Jj, u / (1 - u), a, c, b, q)
        if u == 1: return "both_small_u=1", (Tt, Jj, Fr(1), a, c, b, q)
        Uu = u - 1; cap = e / Uu
        if b >= j: s = (b - j) / (cap - j); return "both_small_u>1_b>=j", (Tt, Jj, Uu, a, c, s / (1 - s), q)
        s = b / j; return "lower_b", (Tt, Jj, Uu, a, c, s / (1 - s), q)
    if j > 1:
        Jj = j - 1; s = p * j; Tt = s / (1 - s)
        if u <= 1:
            s2 = a / (1 / Jj - e); Aa = s2 / (1 - s2)
            return ("j>1_u<1", (Tt, Jj, u / (1 - u), Aa, c, b, q)) if u < 1 else ("j>1_u=1", (Tt, Jj, Fr(1), Aa, c, b, q))
        x, y = u - 1, j * a; sx = Jj * (x + y); Dd = 1 / (1 - sx)
        Uu, Aa = Jj * Dd * x, Jj * Dd * y
        cap = (Jj * Dd + Uu) / (j * Uu); s3 = b / cap
        return "j>1_u>1_simplex", (Tt, Jj, Uu, Aa, c, s3 / (1 - s3), q)
    if p > 1:
        Jj = p - 1; s = j * p; Tt = s / (1 - s); s4 = c * Jj; Cc_ = s4 / (1 - s4)
        if u < 1: return "p>1_u<1", (Tt, Jj, u / (1 - u), a, Cc_, b, q)
        if u == 1: return "p>1_u=1", (Tt, Jj, Fr(1), a, Cc_, b, q)
        Uu = u - 1; s5 = b * Uu / e; return "p>1_u>1", (Tt, Jj, Uu, a, Cc_, s5 / (1 - s5), q)
    if p == 1:
        Jj = j / (1 - j)
        if u < 1: return "p=1_u<1", (Fr(1), Jj, u / (1 - u), a, c, b, q)
        if u == 1: return "p=1_u=1", (Fr(1), Jj, Fr(1), a, c, b, q)
        Uu = u - 1; s5 = b * Uu / e; return "p=1_u>1", (Fr(1), Jj, Uu, a, c, s5 / (1 - s5), q)
    assert j == 1
    Tt = p / (1 - p)
    if u < 1: return "j=1_u<1", (Tt, Fr(1), u / (1 - u), a, c, b, q)
    if u == 1: return "j=1_u=1", (Tt, Fr(1), Fr(1), a, c, b, q)
    Uu = u - 1
    if b >= 1: s = (b - 1) * Uu; return "j=1_u>1_b>=1", (Tt, Fr(1), Uu, a, c, s / (1 - s), q)
    return "j=1_lower_b", (Tt, Fr(1), Uu, a, c, b / (1 - b), q)
CH = {n: (v, st) for n, v, st, *_ in CHARTS}
CH["lower_b"] = (Zv, {1, 2, 6})
CH["j=1_lower_b"] = ([fr(Q), pS, fr(one), fr(Cc), fr(1 + U), fr(B, 1 + B), fr(A)], {2, 6})
def Gval(x):
    s = Fr(0)
    for ex, co in G.to_dict().items():
        t = Fr(int(co))
        for i, k in enumerate(ex):
            if k: t *= x[i] ** int(k)
        s += t
    return s
grid_vals = [Fr(0), Fr(1, 7), Fr(1, 3), Fr(1, 2), Fr(2, 3), Fr(1), Fr(3, 2), Fr(2), Fr(3), Fr(7)]
pos_vals = grid_vals[1:]
rng = random.Random(38)
npts, bad, counts, gneg = 0, [], {}, 0
pts = []
for _ in range(20000):
    q = rng.choice(pos_vals); p = rng.choice(grid_vals); j = rng.choice(pos_vals); c = rng.choice(grid_vals)
    e = rng.choice(pos_vals + [Fr(1) / j, Fr(3, 2) / j]); b = rng.choice(grid_vals + [j]); a = rng.choice(pos_vals)
    if rng.random() < 0.3: p = Fr(1) / j if rng.random() < 0.5 else Fr(1)   # boundary layers
    if rng.random() < 0.2: j = Fr(1)
    if hard(q, p, j, c, e, b, a): pts.append((q, p, j, c, e, b, a))
for x in pts[:3000]:
    npts += 1
    name, y = locate(*x)
    counts[name] = counts.get(name, 0) + 1
    vals, st = CH[name]
    okdom = all(v >= 0 for v in y) and all(y[i] > 0 for i in st)
    okmap = tuple(ev(nd, y) for nd in vals) == x
    if not (okdom and okmap): bad.append((x, name, y))
    if Gval(x) <= 0: gneg += 1
check(f"cover sanity: {npts} exact hard-region grid points (incl. layers p=1, j=1, u0=1, b=j) each lie in the image of the chart chosen by the covering argument, with chart parameters in its domain",
      not bad and npts > 0, f"per chart {sorted(counts.items())}; failures {bad[:2]}")
check(f"G9 > 0 at all {npts} grid points (exact evaluation)", gneg == 0)
info("elapsed", f"{time.time()-T0:.1f}s")
print(f"SUMMARY: {len(FAIL)} failed checks", flush=True)
sys.exit(1 if FAIL else 0)
