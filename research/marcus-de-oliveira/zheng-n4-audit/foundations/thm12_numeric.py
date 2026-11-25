"""Independent numerical and exact checks of Theorem 1.2 (Zheng, MOC n=4, Section 2).

Own code only; uses no file of the supplement.  Usage:
    python -I -B thm12_numeric.py [quick|full]
Deterministic seeds.  Sections:
  E0  exact (sympy) checks of eq. (2.2) signs and of Gale invariance of q_A
  E1  numeric checks on a complex elliptic normal quartic: the q-vectors of x(L), L in Pic^4(E),
      span a 2-dimensional space (deg W = 2, f maps to a line); Gale duality g(x(L)) = x(iota L)
  E2  the paper's construction on a real elliptic curve with 8 real marked points: the circle B,
      the sign pattern of the Gale self-conjugacy relation, unitarity along B, real/orthogonal
      endpoints at the ramification points, and r(U_L) = lam r(O-) + (1-lam) r(O+) along B
  N1  multistart least squares for (O-, O+, lam) on Haar and structured unitaries
  N2  min over U(4) vs min over SO(4) of linear functionals of r (random, and the 48 h_{pi,eps})
"""
import sys, time, itertools
import numpy as np
from scipy.linalg import expm
from scipy.optimize import least_squares, minimize

MODE = sys.argv[1] if len(sys.argv) > 1 else "full"
FULL = MODE == "full"
T0 = time.time()
def log(*a):
    print(*a, flush=True)
def stamp(s):
    log(f"\n=== {s}   [t={time.time()-T0:.1f}s]")

# ---------------------------------------------------------------- minors and r(V)
SUBS = {k: list(itertools.combinations(range(4), k)) for k in range(5)}
PAIRS = [(I, J) for k in range(5) for I in SUBS[k] for J in SUBS[k]]          # 70 pairs, (I,J)
assert len(PAIRS) == 70
IDX = {}
for k in (1, 2, 3):
    Is = np.array([I for I in SUBS[k] for J in SUBS[k]]); Js = np.array([J for I in SUBS[k] for J in SUBS[k]])
    IDX[k] = (Is, Js)

def r_vec(V):
    """vector of |det V[I,J]|^2 in PAIRS order"""
    out = [np.ones(1)]
    for k in (1, 2, 3):
        Is, Js = IDX[k]
        sub = V[Is[:, :, None], Js[:, None, :]]
        out.append(np.abs(np.linalg.det(sub)) ** 2)
    out.append(np.array([abs(np.linalg.det(V)) ** 2]))
    return np.concatenate(out)

A_OF = []   # A = I^c u (4+J) as 0-based subset of range(8)
SIGN = []   # (-1)^{k + sum I + sum J} with 1-based indices
for (I, J) in PAIRS:
    Ic = [i for i in range(4) if i not in I]
    A_OF.append(tuple(sorted(Ic + [4 + j for j in J])))
    SIGN.append((-1) ** (len(I) + sum(i + 1 for i in I) + sum(j + 1 for j in J)))
SIGN = np.array(SIGN, dtype=float)
ALLA = list(itertools.combinations(range(8), 4))
assert sorted(A_OF) == sorted(ALLA)

def q_vec(V, subsets=A_OF):
    out = []
    for A in subsets:
        Ac = [i for i in range(8) if i not in A]
        out.append(np.linalg.det(V[:, list(A)]) * np.linalg.det(V[:, Ac]))
    return np.array(out)

def haar(rng, n=4):
    Z = (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    return Q * (np.diag(R) / np.abs(np.diag(R)))

def so4(th):
    K = np.zeros((4, 4)); K[np.triu_indices(4, 1)] = th; return expm(K - K.T)

def u4(th):
    H = np.zeros((4, 4), complex)
    H[np.triu_indices(4, 1)] = th[:6] + 1j * th[6:12]
    H = H + H.conj().T + np.diag(th[12:16])
    return expm(1j * H)

# ---------------------------------------------------------------- E0 exact checks
def section_E0():
    stamp("E0  exact checks (sympy)")
    import sympy as sp
    rng = np.random.default_rng(1)
    # Gaussian-rational unitary via Cayley transform of a skew-Hermitian K
    for trial in range(2 if FULL else 1):
        a = [sp.Rational(int(x), 7) for x in rng.integers(-9, 10, 16)]
        K = sp.zeros(4, 4)
        it = iter(a)
        for i in range(4):
            K[i, i] = sp.I * next(it)
        for i in range(4):
            for j in range(i + 1, 4):
                z = next(it) + sp.I * sp.Rational(int(rng.integers(-9, 10)), 5)
                K[i, j] = z; K[j, i] = -sp.conjugate(z)
        Id = sp.eye(4)
        U = sp.simplify((Id - K) * (Id + K).inv())
        assert sp.simplify(U * U.H - Id) == sp.zeros(4, 4)
        V = Id.row_join(U)
        q0 = V[:, [0, 1, 2, 3]].det() * V[:, [4, 5, 6, 7]].det()
        bad = 0
        for n, (I, J) in enumerate(PAIRS):
            A = list(A_OF[n]); Ac = [i for i in range(8) if i not in A]
            qA = sp.expand(V[:, A].det() * V[:, Ac].det())
            m = U.extract(list(I), list(J)).det() if I else sp.Integer(1)
            rhs = int(SIGN[n]) * sp.expand(m * sp.conjugate(m))
            if sp.nsimplify(sp.expand(qA / q0 - rhs)) != 0:
                bad += 1
        log(f"  (2.2) exact, Gaussian-rational unitary #{trial}: mismatches among 70 = {bad}")
    # Gale invariance of q_A for a random integer 4x8 matrix: q_A(ker basis) = c * q_A(V)
    Vi = sp.Matrix(rng.integers(-5, 6, (4, 8)).tolist())
    Wk = sp.Matrix.hstack(*Vi.nullspace()).T          # rows span ker V
    assert (Vi * Wk.T) == sp.zeros(4, 4)
    ratios = set()
    for A in ALLA:
        Ac = [i for i in range(8) if i not in A]
        qv = Vi[:, list(A)].det() * Vi[:, Ac].det(); qw = Wk[:, list(A)].det() * Wk[:, Ac].det()
        if qv != 0:
            ratios.add(sp.nsimplify(qw / qv))
        elif qw != 0:
            ratios.add("inf")
    log(f"  Gale invariance (exact, integer V): distinct ratios q_A(g V)/q_A(V) over 70 A: {len(ratios)} -> {list(ratios)[:2]}")
    # dimension of span of q-vectors (HMSV: 14)
    M = np.array([q_vec(rng.standard_normal((4, 8)) + 1j * rng.standard_normal((4, 8)), ALLA) for _ in range(60)])
    s = np.linalg.svd(M, compute_uv=False)
    log(f"  span of q-vectors of 60 random configurations: singular values s14={s[13]:.3e}, s15={s[14]:.3e} (numeric rank 14 expected)")
    # (2.2) numerically for Haar unitaries
    err = 0
    for _ in range(200):
        U = haar(rng); V = np.hstack([np.eye(4), U])
        q = q_vec(V); err = max(err, np.max(np.abs(q / np.linalg.det(U) - SIGN * r_vec(U))))
    log(f"  (2.2) numeric, 200 Haar U: max |q_A/q_0 - sign*|minor|^2| = {err:.2e}")

# ---------------------------------------------------------------- elliptic curve tools
class Curve:
    """y^2 = x^3 + a x + b; points are (x,y) complex tuples or None (= O)"""
    def __init__(self, a, b): self.a, self.b = a, b
    def f(self, x): return x ** 3 + self.a * x + self.b
    def neg(self, P): return None if P is None else (P[0], -P[1])
    def add(self, P, Q):
        if P is None: return Q
        if Q is None: return P
        x1, y1 = P; x2, y2 = Q
        if abs(x1 - x2) < 1e-13 * (1 + abs(x1)):
            if abs(y1 + y2) < 1e-11 * (1 + abs(y1)):
                return None
            m = (3 * x1 ** 2 + self.a) / (2 * y1)
        else:
            m = (y2 - y1) / (x2 - x1)
        x3 = m * m - x1 - x2
        return (x3, -(y1 + m * (x3 - x1)))
    def sum(self, Ps):
        S = None
        for P in Ps: S = self.add(S, P)
        return S

def config(C, pts, Q):
    """4x8 matrix of x(L_Q), L_Q = O(3O + Q): basis 1, x, y, (y + y_Q)/(x - x_Q)  (Q != O)"""
    xQ, yQ = Q
    return np.array([[1, x, y, (y + yQ) / (x - xQ)] for (x, y) in pts], dtype=complex).T

def gale_relation(V, W):
    """least singular values of the system  V diag(d) W^T = 0  (unknown d in C^8)"""
    M = np.array([np.outer(V[:, i], W[:, i]).ravel() for i in range(8)]).T
    _, s, vh = np.linalg.svd(M)
    return s, vh[-1].conj()

def section_E1():
    stamp("E1  complex elliptic curve: deg W = 2 (q-vectors on a line) and Gale duality (2.4)")
    rng = np.random.default_rng(2)
    C = Curve(complex(rng.standard_normal(), rng.standard_normal()), complex(rng.standard_normal(), rng.standard_normal()))
    pts = []
    for _ in range(8):
        x = complex(rng.standard_normal(), rng.standard_normal()); pts.append((x, np.sqrt(C.f(x))))
    S = C.sum(pts)
    Qs = []
    for _ in range(40):
        x = complex(rng.standard_normal(), rng.standard_normal()); Qs.append((x, np.sqrt(C.f(x))))
    M = []
    for Q in Qs:
        q = q_vec(config(C, pts, Q), ALLA); M.append(q / np.linalg.norm(q))
    s = np.linalg.svd(np.array(M), compute_uv=False)
    log(f"  40 random L in Pic^4(E): singular values of stacked q-vectors s1..s4 = {s[:4]}")
    log(f"  -> numeric rank 2 (s3/s1 = {s[2]/s[0]:.2e}): all f(L) lie on one projective line")
    worst = 0; gap = np.inf; minabs = np.inf
    for Q in Qs[:20]:
        V = config(C, pts, Q); Wd = config(C, pts, C.add(S, C.neg(Q)))
        sv, d = gale_relation(V, Wd)
        worst = max(worst, sv[-1] / sv[0]); gap = min(gap, sv[-2] / sv[0]); minabs = min(minabs, np.min(np.abs(d)) / np.max(np.abs(d)))
        # and q-vectors of x(L), x(iota L) proportional
    log(f"  Gale (2.4): max s8/s1 = {worst:.2e}, min s7/s1 = {gap:.2e} (1-dim solution), min |d_i|/max|d_i| = {minabs:.2e} (all d_i nonzero)")
    # ramification: f(Q) = f(S-Q)
    Q = Qs[0]
    q1 = q_vec(config(C, pts, Q), ALLA); q2 = q_vec(config(C, pts, C.add(S, C.neg(Q))), ALLA)
    log(f"  f(L) = f(iota L): |q1/|q1| - phase*q2/|q2|| = {np.linalg.norm(q1/np.linalg.norm(q1) - (np.vdot(q2,q1)/abs(np.vdot(q2,q1)))*q2/np.linalg.norm(q2)):.2e}")

def real_halve(C, S, e1):
    """real H on the identity component with 2H = S"""
    a, b = C.a, C.b; xs = S[0].real
    # x(2H) = ((3x^2+a)/(2y))^2 - 2x = xs  <=>  (3x^2+a)^2 - 4 f(x)(2x + xs) = 0
    poly = np.polysub(np.polymul([3, 0, a], [3, 0, a]), np.polymul([4, 0, 4 * a, 4 * b], [2, xs]))
    best = None
    for x in np.roots(poly):
        if abs(x.imag) < 1e-9 and x.real >= e1:
            for sg in (1, -1):
                H = (complex(x.real), sg * np.sqrt(complex(C.f(x.real))))
                D = C.add(H, H)
                err = abs(D[0] - S[0]) + abs(D[1] - S[1])
                if best is None or err < best[0]: best = (err, H)
    return best

def unitary_from_config(V, d):
    """V D V^* = 0 with D = diag(d) real, signs ++++----: return U unitary representing V"""
    Z0 = np.linalg.solve(V[:, :4], V[:, 4:])
    DI, DJ = d[:4], -d[4:]
    return np.diag(DI ** -0.5) @ Z0 @ np.diag(DJ ** 0.5)

def section_E2():
    stamp("E2  the paper's construction on a real curve with 8 real marked points")
    rng = np.random.default_rng(3)
    allres = []
    ncurves = 6 if FULL else 2
    for cnum in range(ncurves):
        e = np.sort(rng.uniform(-2, 2, 3))[::-1]; e = e - e.mean()
        e1, e2, e3 = e
        a = e1 * e2 + e1 * e3 + e2 * e3; b = -e1 * e2 * e3
        C = Curve(a, b)
        pts = []
        while len(pts) < 8:
            x = rng.uniform(e3, e1 + 3)
            if C.f(x) > 0.05:
                pts.append((complex(x), complex(np.sqrt(C.f(x)) * rng.choice([-1, 1]))))
        S = C.sum(pts)
        if S is not None and S[0].real < e1:          # S on the egg: move it to the identity component
            pts[7] = C.add(pts[7], (complex(e3), 0j)); S = C.sum(pts)
        err, H = real_halve(C, S, e1)
        # imaginary real form {W : conj W = -W}: W = (x, i t), t^2 = -f(x); components x<=e3 (with O) and e2<=x<=e1
        comps = {"C0 (x<=e3, through O)": (None, (complex(e3), 0j), lambda s: e3 - s),
                 "C1 (oval e2<=x<=e1)": ((complex(e2), 0j), (complex(e1), 0j), lambda s: e2 + s * (e1 - e2))}
        for name, (Wa, Wb, xmap) in comps.items():
            if name.startswith("C0"):
                svals = np.concatenate([np.geomspace(1e-3, 1e3, 120)])
            else:
                svals = np.linspace(0.002, 0.998, 120)
            Ws = [(complex(xmap(s)), 1j * np.sqrt(-C.f(xmap(s)))) for s in svals]
            Qend = [C.add(H, Wa), C.add(H, Wb)]
            Qs = [C.add(H, W) for W in Ws]
            # sign pattern of the Gale self-conjugacy relation V D V^* = 0 along the component
            pats = set(); dims = []; imag = 0; ds = []
            for Q in Qend + Qs:
                V = config(C, pts, Q)
                sv, d = gale_relation(V, V.conj())
                dims.append((sv[-1] / sv[0], sv[-2] / sv[0]))
                d = d / d[0]; imag = max(imag, np.max(np.abs(d.imag)) / np.max(np.abs(d.real))); d = d.real
                ds.append(d); pats.add(tuple(np.sign(d).astype(int)))
            log(f"  curve {cnum} e={np.round(e,3)} {name}: |H-halving err|={err:.1e}; sign patterns along component: {pats}; "
                f"max s8/s1={max(x[0] for x in dims):.1e}, min s7/s1={min(x[1] for x in dims):.1e}, max imag/real of d={imag:.1e}, "
                f"min |d_i|={min(np.min(np.abs(d)) for d in ds):.2e}")
            if len(pats) != 1:
                log("   !! sign pattern NOT constant on the component"); continue
            pat = np.array(next(iter(pats)))
            if (pat > 0).sum() != 4:
                log("   !! pattern is not 4+/4-"); continue
            perm = list(np.where(pat > 0)[0]) + list(np.where(pat < 0)[0])
            Us = []; uerr = 0
            for Q, d in zip(Qend + Qs, ds):
                V = config(C, pts, Q)[:, perm]; dd = d[perm]
                U = unitary_from_config(V, dd); uerr = max(uerr, np.max(np.abs(U @ U.conj().T - np.eye(4)))); Us.append(U)
            Om, Op = Us[0], Us[1]
            realerr = max(np.max(np.abs(Om.imag)), np.max(np.abs(Op.imag)))
            rm, rp = r_vec(Om), r_vec(Op)
            lam_all = []; res = 0
            for U in Us[2:]:
                r = r_vec(U)
                lam = np.dot(r - rp, rm - rp) / np.dot(rm - rp, rm - rp)
                res = max(res, np.max(np.abs(r - lam * rm - (1 - lam) * rp))); lam_all.append(lam)
            lam_all = np.array(lam_all)
            allres.append(res)
            log(f"     relabel perm={perm}: max unitarity err={uerr:.1e}; endpoints imag part={realerr:.1e} "
                f"(orthogonality err {max(np.max(np.abs(Om.real@Om.real.T-np.eye(4))),np.max(np.abs(Op.real@Op.real.T-np.eye(4)))):.1e}); "
                f"|r(O-)-r(O+)|={np.linalg.norm(rm-rp):.3f}")
            log(f"     along B: max |r(U_L) - lam r(O-) - (1-lam) r(O+)|_inf = {res:.2e}; lam in [{lam_all.min():.6f}, {lam_all.max():.6f}]; "
                f"lam monotone: {bool(np.all(np.diff(lam_all) > 0) or np.all(np.diff(lam_all) < 0))}")
    log(f"  E2 summary: worst collinearity residual over all circles = {max(allres) if allres else float('nan'):.2e}")

# ---------------------------------------------------------------- N1 optimizer
def fit_two_orthogonal(rU, rng, starts):
    best = (np.inf, None)
    def F(p):
        lam = np.sin(p[12]) ** 2
        return rU - lam * r_vec(so4(p[:6])) - (1 - lam) * r_vec(so4(p[6:12]))
    for s in range(starts):
        p0 = np.concatenate([rng.uniform(-np.pi, np.pi, 12), [rng.uniform(0, np.pi / 2)]])
        sol = least_squares(F, p0, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
        v = np.max(np.abs(sol.fun))
        if v < best[0]: best = (v, sol.x)
        if v < 1e-12: break
    return best

def structured(rng):
    out = []
    w = np.exp(2j * np.pi / 4)
    out.append(("Fourier F4", np.array([[w ** (i * j) for j in range(4)] for i in range(4)]) / 2))
    out.append(("F4 x diag phases", np.diag(np.exp(1j * rng.uniform(0, 6.3, 4))) @ out[0][1] @ np.diag(np.exp(1j * rng.uniform(0, 6.3, 4)))))
    h2 = haar(rng, 2); h2b = haar(rng, 2)
    B = np.zeros((4, 4), complex); B[:2, :2] = h2; B[2:, 2:] = h2b
    out.append(("block diag U(2)xU(2)", B))
    B3 = np.zeros((4, 4), complex); B3[:3, :3] = haar(rng, 3); B3[3, 3] = np.exp(0.7j)
    out.append(("block diag U(3)xU(1)", B3))
    out.append(("Kronecker U(2)(x)U(2)", np.kron(haar(rng, 2), haar(rng, 2))))
    P = np.eye(4)[[1, 3, 0, 2]]
    out.append(("permutation", P.astype(complex)))
    for eps in (1e-3, 1e-1):
        out.append((f"near permutation eps={eps}", P @ u4(eps * rng.standard_normal(16))))
    Ost = np.array([[1, 4, 4, 4], [-4, 1, 4, -4], [-4, -4, 1, 4], [-4, 4, -4, 1]]) / 7.0
    out.append(("O* of Example 1.3", Ost.astype(complex)))
    for eps in (1e-2, 1e-1, 0.3):
        out.append((f"O* exp(i eps H) eps={eps}", Ost @ u4(eps * rng.standard_normal(16))))
    out.append(("diag phases", np.diag(np.exp(1j * rng.uniform(0, 6.3, 4)))))
    return out

def section_N1():
    stamp("N1  multistart least squares for (O-, O+, lam), residual = max over the 70 equations")
    rng = np.random.default_rng(4)
    nh = 120 if FULL else 10
    res = []
    hard = []
    for k in range(nh):
        U = haar(rng); v, _ = fit_two_orthogonal(r_vec(U), rng, 30)
        res.append(v)
        if v > 1e-10: hard.append((k, U, v))
    res = np.array(res)
    log(f"  Haar U(4), n={nh} (seed 4): residual quantiles min/median/90%/max = {res.min():.1e}/{np.median(res):.1e}/{np.quantile(res,0.9):.1e}/{res.max():.1e}; #>1e-10: {(res>1e-10).sum()}")
    for k, U, v in hard:
        v2, _ = fit_two_orthogonal(r_vec(U), rng, 200)
        log(f"   retry Haar #{k}: first {v:.2e} -> with 200 starts {v2:.2e}")
    for name, U in structured(rng):
        v, _ = fit_two_orthogonal(r_vec(U), rng, 60)
        if v > 1e-10:
            v = min(v, fit_two_orthogonal(r_vec(U), rng, 300)[0])
        log(f"  {name:28s}: unitarity err {np.max(np.abs(U@U.conj().T-np.eye(4))):.1e}; best residual {v:.2e}")
    # does the found segment consist of unitary points? (prediction of the construction, not of the statement)
    log("  segment test: fit a unitary to interior points of the found segment [r(O-), r(O+)]")
    for k in range(2 if FULL else 1):
        U = haar(rng); v, p = fit_two_orthogonal(r_vec(U), rng, 30)
        rm, rp = r_vec(so4(p[:6])), r_vec(so4(p[6:12])); lam = np.sin(p[12]) ** 2
        out = []
        for mu in (0.25, 0.5, 0.75):
            tgt = mu * rm + (1 - mu) * rp
            b = np.inf
            for s in range(15):
                sol = least_squares(lambda t: r_vec(u4(t)) - tgt, rng.uniform(-np.pi, np.pi, 16), method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
                b = min(b, np.max(np.abs(sol.fun)))
                if b < 1e-11: break
            out.append(b)
        log(f"   Haar #{k}: fit residual {v:.1e}, lam={lam:.4f}, |r(O-)-r(O+)|={np.linalg.norm(rm-rp):.3f}; unitary-fit residual at mu=.25/.5/.75: {out[0]:.1e}/{out[1]:.1e}/{out[2]:.1e}")

# ---------------------------------------------------------------- N2 linear functionals
PIDX = {p: n for n, p in enumerate(PAIRS)}
def coeff(I, J):
    return PIDX[(tuple(sorted(I)), tuple(sorted(J)))]

def h_functional(pi, eps):
    c = np.zeros(70)
    for i in range(4): c[coeff((i,), (pi[i],))] += 0.5
    c[coeff((0, 1), (pi[0], pi[2]))] += eps * 0.5
    c[coeff((0, 1), (pi[1], pi[3]))] += eps * 0.5
    c[coeff((0, 2), (pi[0], pi[1]))] -= eps * 0.5
    c[coeff((0, 2), (pi[2], pi[3]))] -= eps * 0.5
    return c

def minimize_over(c, group, rng, starts, want_arg=False):
    if group == "U":
        f = lambda t: c @ r_vec(u4(t)); n = 16
    else:
        f = lambda t: c @ r_vec(so4(t)); n = 6
    best = (np.inf, None)
    for s in range(starts):
        sol = minimize(f, rng.uniform(-np.pi, np.pi, n), method="BFGS", options={"gtol": 1e-10})
        if sol.fun < best[0]: best = (sol.fun, sol.x)
    return best if want_arg else best[0]

def section_N2():
    stamp("N2  conv r(U(4)) = conv r(O(4)): minima of linear functionals")
    rng = np.random.default_rng(5)
    Ost = np.array([[1, 4, 4, 4], [-4, 1, 4, -4], [-4, -4, 1, 4], [-4, 4, -4, 1]]) / 7.0
    psi = h_functional((0, 1, 2, 3), -1) @ r_vec(Ost)
    log(f"  sanity: h_(id,-1)(O*) = {psi:.10f}  vs  -158/2401 = {-158/2401:.10f}")
    perms = list(itertools.permutations(range(4)))
    vals = [h_functional(p, e) @ r_vec(np.eye(4)[list(s)]) for p in perms for e in (1, -1) for s in perms]
    log(f"  sanity: values of the 48 h at the 24 permutation matrices: min {min(vals):.3g}, max {max(vals):.3g}")
    worst = -np.inf; mins_U = []; mins_O = []
    sel = [(p, e) for p in perms for e in (1, -1)]
    if not FULL: sel = sel[:4]
    for p, e in sel:
        c = h_functional(p, e)
        mU = minimize_over(c, "U", rng, 12 if FULL else 4); mO = minimize_over(c, "O", rng, 12 if FULL else 4)
        mins_U.append(mU); mins_O.append(mO); worst = max(worst, mO - mU)
    log(f"  48 h_(pi,eps): min over U(4) in [{min(mins_U):.10f}, {max(mins_U):.10f}], min over SO(4) in [{min(mins_O):.10f}, {max(mins_O):.10f}];"
        f" max (minO - minU) = {worst:.2e}   (paper: -0.0674422489 for both)")
    worst = -np.inf; worst_cert = -np.inf; nr = 40 if FULL else 4; nfail = 0
    for k in range(nr):
        c = rng.standard_normal(70)
        mU, thU = minimize_over(c, "U", rng, 10 if FULL else 3, want_arg=True); mO = minimize_over(c, "O", rng, 10 if FULL else 3)
        # certificate: write r(U_min) = lam r(O-) + (1-lam) r(O+); then min(c.r(O-), c.r(O+)) <= c.r(U_min)
        v, p = fit_two_orthogonal(r_vec(u4(thU)), rng, 30)
        mcert = min(c @ r_vec(so4(p[:6])), c @ r_vec(so4(p[6:12])))
        if mO - mU > 1e-8:
            nfail += 1
            log(f"   functional #{k}: SO(4)-BFGS min {mO:.8f} > U(4) min {mU:.8f}; two-orthogonal fit of the U(4) minimizer (residual {v:.1e}) gives an O with c.r(O) = {mcert:.8f}")
        worst = max(worst, mO - mU); worst_cert = max(worst_cert, min(mO, mcert) - mU)
    log(f"  {nr} Gaussian random functionals: raw max (BFGS min over SO(4) - min over U(4)) = {worst:.2e} ({nfail} cases where SO(4) multistart missed);"
        f" after adding the orthogonal points from the two-orthogonal fit: max = {worst_cert:.2e} (>0 beyond tolerance would contradict Thm 1.2)")

if __name__ == "__main__":
    log(f"thm12_numeric.py mode={MODE}; numpy {np.__version__}")
    which = sys.argv[2].split(",") if len(sys.argv) > 2 else ["E0", "E1", "E2", "N1", "N2"]
    for s in which:
        globals()["section_" + s]()
    stamp("done")
