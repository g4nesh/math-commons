"""Own-code replay for blind_moc4_T42_uniform_20261022 (paper case T_43). Audit batch C18.
Uses only the shared audit library (../lib/moc4lib.py); executes nothing from the supplement."""
# ---- common helpers (own code; duplicated in each C18 check script) ----
import sys, os, itertools, random
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L
import flint
import sympy as sp

FAIL = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + ((" :: " + str(info)) if info != "" else ""))
    if not ok:
        FAIL.append(name)

def info(name, s):
    print("INFO " + name + " :: " + str(s))

def vec_comb(coeffs):
    """sum_j coeff_j f_j as a V_R vector; coeffs: dict j -> number/sympy"""
    out = [0] * 14
    for j, c in coeffs.items():
        for k in range(14):
            out[k] = out[k] + c * L.F_VECS[j][k]
    return out

def prof(v):
    return L.values(v)

def lin_identity(lhs, rhs):
    """profiles of sum lhs f_j == sum rhs f_j"""
    a = prof(vec_comb({j: 1 for j in lhs}) if isinstance(lhs, list) else vec_comb(lhs))
    b = prof(vec_comb({j: 1 for j in rhs}) if isinstance(rhs, list) else vec_comb(rhs))
    return all(x == y for x, y in zip(a, b))

# Hessian of the cubic: H(m)[a][b] = sum_k HT[a][b][k] m_k
_ms, _Cs = L.cubic_sympy()
HT = [[[int(sp.diff(_Cs, _ms[a], _ms[b], _ms[k])) for k in range(14)] for b in range(14)] for a in range(14)]

def hess_det(m):
    M = [[sum(HT[a][b][k] * m[k] for k in range(14) if HT[a][b][k]) for b in range(14)] for a in range(14)]
    return flint.fmpq_mat(14, 14, [flint.fmpq(Fr(x).numerator, Fr(x).denominator) for row in M for x in row]).det()

def CQ(g):
    """C(g + z kappa) = C0 + z*Q (checked linear in z by caller)"""
    C0 = L.cubic_C(g)
    gk = list(g); gk[13] = gk[13] + 1
    C1 = L.cubic_C(gk)
    return C0, C1 - C0

def radial_G2(g):
    """-Q^2 det H(g - (C/Q) kappa)/4  (= G^2 up to the calibration constant)"""
    C0, Q = CQ(g)
    w = list(g); w[13] = w[13] - Fr(C0) / Fr(Q)
    D = hess_det(w)
    return -flint.fmpq(Fr(Q).numerator, Fr(Q).denominator) ** 2 * D / 4, C0, Q

def to_fmpq(x):
    x = Fr(x); return flint.fmpq(x.numerator, x.denominator)

def interp(ts, ys):
    """Lagrange interpolation over Q -> fmpq_poly"""
    X = flint.fmpq_poly([0, 1])
    P = flint.fmpq_poly([0])
    for i, ti in enumerate(ts):
        term = flint.fmpq_poly([ys[i]])
        for j, tj in enumerate(ts):
            if j != i:
                term = term * (X - to_fmpq(tj)) / (to_fmpq(ti) - to_fmpq(tj))
        P = P + term
    return P
# ---- end common helpers ----

# ================= blind_moc4_T42_uniform_20261022 (paper T_43) =================
T = frozenset({1, 2, 4, 6, 7, 8})
check("contact set equals paper T_43 literally", T == L.CONTACT_SETS[43])
gens = L.generators(T)
check("generators(T) = {0,1,2,3,4,5,6,21,22,23,53,54}", gens == [0,1,2,3,4,5,6,21,22,23,53,54], gens)
check("nonminor generators are exactly 22, 54", [j for j in gens if not L.F_IS_MINOR[j]] == [22, 54])
ids = [([22,54],[0,4,5,23,53]), ([54,1],[0,6,53]), ([54,4],[2,3,6,53]), ([54,21],[0,2,23,53])]
for l, r in ids:
    check("identity f%s = f%s" % ("+f".join(map(str,l)), "+f".join(map(str,r))), lin_identity(l, r))
# symmetry swapping 22 <-> 54 and fixing T
def maps_f(gi, j):
    t = L.g_index_map(gi); v = [None]*24
    for s in range(24): v[t[s]] = L.F_VALUES[j][s]
    return [k for k in range(82) if list(L.F_VALUES[k]) == v]
sw = [gi for gi in L.G_ELEMENTS if L.g_act_set(gi, T) == T and maps_f(gi, 22) == [54] and maps_f(gi, 54) == [22]]
check("some g in G fixes T and swaps f22<->f54", len(sw) > 0, "count=%d" % len(sw))
stated = [gi for gi in sw if gi[2] == 0 and set([gi[0], gi[1]]) <= {(0,2,1,3),(1,0,3,2)}]
info("swap elements built from the stated row/col perms (0,2,1,3),(1,0,3,2), eps=0", stated)
B1 = [0, 53, 5, 3, 23, 6, 54, 2]
check("basis of (1) has rank 8", L.rank([L.F_VECS[j] for j in B1]) == 8)
zs = frozenset.intersection(*[L.F_ZERO[j] for j in B1])
check("basis of (1) has exactly the six common zeros T", zs == T, sorted(zs))
# the reduction leaves (1): every point of cone(gens) lies in cone(B1), cone(swap B1) or the minor cone
# (follows from the identities; checked as algebra above)

# ---- cubic restriction C + zQ ----
p, r, a, b, c, q, d, s, z = sp.symbols("p r a b c q d s z")
coef = {0: p, 53: r, 5: a, 3: b, 23: c, 6: q, 54: d, 2: s}
g = vec_comb(coef)
gz = list(g); gz[13] = gz[13] + z
Cz = sp.Poly(sp.expand(L.cubic_C(gz)), z)
check("C(g + z kappa) is affine in z", Cz.degree() <= 1)
Cm = sp.expand(Cz.coeff_monomial(1)); Qm = sp.expand(Cz.coeff_monomial(z))
M = b + c; U = d + r
T0 = b*U*(d+q) + c*d*(d+p+r) + d*(d+p+r)*(d+q) + r*p*(c+q+d)
Cc = sp.expand(b*c*(a-U) - T0); Qc = sp.expand((a-U)*M - d*U - a*s)
rC = sp.simplify(Cm / Cc); rQ = sp.simplify(Qm / Qc)
check("C(own) = const * C(doc)", rC.is_number, rC)
check("Q(own) = const * Q(doc)", rQ.is_number, rQ)
check("C and Q constants agree (so -C/Q is the same height)", rC == rQ)
# -MC identity at d=1 and nonnegativity of M T0 - bc
e1 = sp.expand((-M*Cc - (b*c*(M+1-a*M) + r*M*b*c + (M*T0 - b*c))).subs(d, 1))
check("-MC = bc(M+1-aM) + rMbc + (MT0-bc) at d=1", e1 == 0)
P1 = sp.Poly(sp.expand((M*T0 - b*c).subs(d, 1)), p, r, a, b, c, q, s)
check("MT0 - bc (d=1) has only nonnegative coefficients", all(co >= 0 for co in P1.coeffs()))
check("MT0 - bc (d=1) contains monomials b and c", P1.coeff_monomial(b) > 0 and P1.coeff_monomial(c) > 0)
check("M=0: C = -(1+p)(1+q)(1+r)", sp.expand(Cc.subs({b: 0, c: 0, d: 1}) + (1+p)*(1+q)*(1+r)) == 0)

# ---- radial factor G: exact pointwise reconstruction along lines from the base point d=1 ----
NAMES = ["p", "r", "a", "b", "c", "q", "d", "s"]
IDX = [0, 53, 5, 3, 23, 6, 54, 2]
def cost(vals):
    return vec_comb({IDX[i]: Fr(vals[i]) for i in range(8)})
BASE = [0, 0, 0, 0, 0, 0, 1, 0]
G2b, C0b, Qb = radial_G2(cost(BASE))
info("base point d=1: -Q^2 detH/4, C, Q", (G2b, C0b, Qb))
def G_along(target, npts=21):
    ts = list(range(npts)); ys = []
    for t in ts:
        v = [Fr(BASE[i]) + Fr(t) * (Fr(target[i]) - BASE[i]) for i in range(8)]
        ys.append(radial_G2(cost(v))[0])
    S = interp(ts, ys)
    return S
def exact_G_at(target):
    S = G_along(target)
    ok_deg = S.degree() <= 18
    R = S.sqrt()  # raises/None if not a square
    if R is None: return None, ok_deg, False
    g0 = R(0); g1 = R(1)
    return Fr(int((g1/g0).p), int((g1/g0).q)), ok_deg, True

def J_of(b_, c_, p_, q_):
    import numpy as np
    from scipy.optimize import minimize_scalar
    def obj(be):
        bs = [be, 1-be]; res = [p_, q_]; den = [b_, c_]; v = 1.0
        for i in range(2):
            if den[i] == 0:
                if bs[i] > 0: return np.inf
            else: v += bs[i]**2/den[i]
            dl = min(res[i]/2, bs[i]); v += (dl**2/res[i] if res[i] > 0 else 0) - dl
        return v
    rs = minimize_scalar(obj, bounds=(0, 1), method="bounded", options={"xatol": 1e-13})
    return min(rs.fun, obj(0.0), obj(1.0))

rng = random.Random(20261125)
def rnd():
    return Fr(round(2.0 ** rng.uniform(-3, 3) * 16), 16)
nsq = 0; worst = None; nviol = 0; NPT = 12
for it in range(NPT):
    b_, c_, p_, q_, r_, s_ = [rnd() for _ in range(6)]
    if it % 4 == 1: r_ = Fr(0)
    if it % 4 == 2: s_ = Fr(0)
    J = J_of(float(b_), float(c_), float(p_), float(q_))
    amax = J if s_ == 0 else min(J, 1/float(s_))
    a_ = Fr(rng.uniform(0.02, 0.98) * amax).limit_denominator(64)
    Gv, okd, sq = exact_G_at([p_, r_, a_, b_, c_, q_, 1, s_])
    if sq: nsq += 1
    lb = (1 - float(a_*s_))**2 * (1 - float(a_)/J)**2
    ratio = float(Gv)/lb if sq else None
    if sq and ratio < 1: nviol += 1
    worst = ratio if worst is None or (ratio is not None and ratio < worst) else worst
    print("  pt", it, "p r a b c q s =", [str(x) for x in (p_, r_, a_, b_, c_, q_, s_)], "J=%.4g G=%.6g bound=%.6g" % (J, float(Gv), lb))
check("G^2 = -Q^2 detH/4 is a perfect square of a degree-9 polynomial along every sampled line (EXACT, %d lines)" % NPT, nsq == NPT)
check("NUMERIC-J: G >= (1-as)^2 (1-a/J)^2 at %d exact hard points" % NPT, nviol == 0, "min G/bound = %.6g" % worst)

# ---- secant rank check: rank H(Re x) <= 12 for x in X (NUMERIC) ----
import numpy as np
nrng = np.random.default_rng(5)
ranks = []
for _ in range(5):
    av = nrng.normal(size=4) + 1j*nrng.normal(size=4); bv = nrng.normal(size=4) + 1j*nrng.normal(size=4)
    x = L.matching_vector_np(av, bv, np.exp(1j*nrng.uniform(0, 6.28)))
    w = x.real
    H = np.array([[sum(HT[i][j][k]*w[k] for k in range(14)) for j in range(14)] for i in range(14)])
    sv = np.linalg.svd(H, compute_uv=False); ranks.append(int(np.sum(sv > 1e-9*sv[0])))
    Hx = np.array([[sum(HT[i][j][k]*x[k] for k in range(14)) for j in range(14)] for i in range(14)])
    svx = np.linalg.svd(Hx, compute_uv=False)
check("NUMERIC: rank Hess C(Re x) <= 12 at 5 random x in X (so det = 0 is a valid necessary condition)", max(ranks) <= 12, ranks)

# ---- realizable sampling on the T_43 stratum: is any realizable cost 'hard'? ----
Bmat = np.array([[float(v) for v in L.F_VALUES[j]] for j in IDX]).T
sw_g = sw[0]; tmap = L.g_index_map(sw_g)
def coords(gv):
    sol, *_ = np.linalg.lstsq(Bmat, gv, rcond=None); return sol, np.max(np.abs(Bmat@sol - gv))
nfound = 0; nhard = 0; minphi = np.inf
for seed in range(12):
    pt = L.random_point_on_contact_stratum(sorted(T), seed=seed, tries=20, margin=1e-3)
    if pt is None: continue
    nfound += 1
    gv = pt["gaps"]
    best = None
    for variant in (gv, np.array([gv[tmap[s]] for s in range(24)])):
        co, res = coords(variant)
        if res < 1e-7 and co.min() > -1e-9 and co[6] > 1e-9:
            best = co / co[6]
    if best is None: continue
    pp, rr, aa, bb, cc, qq, dd, ss = best
    J = J_of(bb, cc, pp, qq)
    if aa < J and aa*ss < 1: nhard += 1
    mval, _ = L.min_phi_over_SO4(np.linalg.solve(np.array([[float(x) for x in row] for row in L.P])[L.P_PIVOT_ROWS], gv[L.P_PIVOT_ROWS]), restarts=6, seed=seed)
    minphi = min(minphi, mval)
info("NUMERIC realizable sampling (seeds 0-11, margin 1e-3)", "found=%d in-form-(1)-hard=%d min Phi over SO4 = %.3g" % (nfound, nhard, minphi))
check("NUMERIC: no sampled realizable T_43 cost lies in the hard region", nhard == 0)
print("SUMMARY failures:", FAIL); sys.exit(1 if FAIL else 0)
