"""Second, independent numerical test of Theorem 1.2 (label thm12b).

Theorem 1.2 (Zheng, n=4): for every U in U(4) there are O-, O+ in O(4) and lam in [0,1]
with |det U[I,J]|^2 = lam |det O-[I,J]|^2 + (1-lam) |det O+[I,J]|^2 for all |I|=|J|.

Own code only; nothing from the supplement is read or executed.
Usage:  python thm12b_numeric.py PART [NAME]   PART in {sign, pencil, ls, trace, trace1 NAME, mp1 NAME, theta, all}
Every random choice is seeded; the seeds are printed.
"""
import sys, itertools, time
import numpy as np
from scipy.linalg import expm

np.set_printoptions(precision=6, suppress=True, linewidth=150)

# ---------------------------------------------------------------- index sets
SUBS = {k: list(itertools.combinations(range(4), k)) for k in range(5)}
IJ = [(I, J) for k in range(5) for I in SUBS[k] for J in SUBS[k]]          # 70 pairs
assert len(IJ) == 70
_BATCH = {}
for k in (1, 2, 3):
    pairs = [(I, J) for (I, J) in IJ if len(I) == k]
    R = np.array([I for I, J in pairs]); C = np.array([J for I, J in pairs])
    _BATCH[k] = (R, C)
POS = {k: [n for n, (I, J) in enumerate(IJ) if len(I) == k] for k in range(5)}


def r_vec(V):
    """70-vector of squared moduli of minors (empty minor = 1)."""
    out = np.empty(70)
    out[POS[0]] = 1.0
    for k in (1, 2, 3):
        R, C = _BATCH[k]
        sub = V[R[:, :, None], C[:, None, :]]
        out[POS[k]] = np.abs(np.linalg.det(sub)) ** 2
    out[POS[4]] = abs(np.linalg.det(V)) ** 2
    return out


def haar_unitary(rng, n=4):
    Z = (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    return Q * (np.diag(R) / np.abs(np.diag(R)))


def haar_orth(rng, n=4):
    Z = rng.standard_normal((n, n))
    Q, R = np.linalg.qr(Z)
    return Q * np.sign(np.diag(R))


def skew(a):
    S = np.zeros((4, 4)); iu = np.triu_indices(4, 1)
    S[iu] = a; return S - S.T


def herm(h):
    H = np.zeros((4, 4), complex)
    iu = np.triu_indices(4, 1)
    H[np.diag_indices(4)] = h[:4]
    H[iu] = h[4:10] + 1j * h[10:16]
    return H + np.triu(H, 1).conj().T


def realness_defect(V, tol=1e-9):
    """max |Im(V_ij V_kl conj(V_il V_kj))| over all i<k, j<l: zero iff V is real up to
    left/right diagonal phases (when the bipartite support graph is connected)."""
    m = 0.0
    for i, k in itertools.combinations(range(4), 2):
        for j, l in itertools.combinations(range(4), 2):
            m = max(m, abs((V[i, j] * V[k, l] * np.conj(V[i, l] * V[k, j])).imag))
    return m


# ---------------------------------------------------------------- test set
def test_matrices():
    L = []
    for k in range(30):
        seed = 1000 + k
        L.append((f"haar[{seed}]", haar_unitary(np.random.default_rng(seed))))
    w = np.exp(2j * np.pi / 4)
    F = np.array([[w ** (j * k) for k in range(4)] for j in range(4)]) / 2
    L.append(("fourier4", F))
    L.append(("orth[2000]", haar_orth(np.random.default_rng(2000)).astype(complex)))
    Ostar = np.array([[1, 4, 4, 4], [-4, 1, 4, -4], [-4, -4, 1, 4], [-4, 4, -4, 1]]) / 7
    L.append(("O*_ex1.3", Ostar.astype(complex)))
    rng = np.random.default_rng(2001)
    W2a, W2b = haar_unitary(rng, 2), haar_unitary(rng, 2)
    B22 = np.zeros((4, 4), complex); B22[:2, :2] = W2a; B22[2:, 2:] = W2b
    L.append(("blockdiag2+2[2001]", B22))
    W3 = haar_unitary(np.random.default_rng(2002), 3)
    B31 = np.eye(4, dtype=complex); B31[:3, :3] = W3
    L.append(("blockdiag3+1[2002]", B31))
    rng = np.random.default_rng(2003)
    X = np.eye(4, dtype=complex); X[:2, :2] = haar_unitary(rng, 2)
    Y = np.eye(4, dtype=complex); Y[1:, 1:] = haar_unitary(rng, 3)
    L.append(("twozeros[2003]", X @ Y))
    rng = np.random.default_rng(2004)
    x = rng.standard_normal(4) + 1j * rng.standard_normal(4); x[0] = 0; x /= np.linalg.norm(x)
    M = np.column_stack([x, haar_unitary(rng)[:, :3]])
    Q, R = np.linalg.qr(M); Q[:, 0] *= R[0, 0] / abs(R[0, 0])
    L.append(("onezero[2004]", Q))
    rng = np.random.default_rng(2005)
    P = np.eye(4)[[1, 3, 0, 2]].astype(complex)
    H = herm(rng.standard_normal(16))
    for eps in (0.3, 0.1, 0.03, 0.01):
        L.append((f"nearperm eps={eps}[2005]", P @ expm(1j * eps * H)))
    rng = np.random.default_rng(2006)
    O = haar_orth(rng); H = herm(rng.standard_normal(16))
    for eps in (0.3, 0.1, 0.03, 0.01):
        L.append((f"nearorth eps={eps}[2006]", O @ expm(1j * eps * H)))
    return L


# ---------------------------------------------------------------- part: signs (2.2), Gale
def part_sign():
    print("== PART sign: eq. (2.2) sign rule and Gale invariance of q_A (own derivation) ==")
    worst22, worstg = 0.0, 0.0
    for seed in range(3000, 3020):
        U = haar_unitary(np.random.default_rng(seed))
        V = np.hstack([np.eye(4), U])
        q0 = np.linalg.det(U)
        r = r_vec(U)
        for n, (I, J) in enumerate(IJ):
            k = len(I)
            Ic = [i for i in range(4) if i not in I]
            A = sorted(Ic + [4 + j for j in J]); Ac = sorted(set(range(8)) - set(A))
            qA = np.linalg.det(V[:, A]) * np.linalg.det(V[:, Ac]) if A else 1
            sgn = (-1) ** (k + sum(i + 1 for i in I) + sum(j + 1 for j in J))
            worst22 = max(worst22, abs(qA / q0 - sgn * r[n]))
        # Gale dual: rows spanning ker V, q-vector proportional to q(V)
        W = np.linalg.svd(V)[2][4:].conj()          # V W^T = 0
        assert np.abs(V @ W.T).max() < 1e-12
        A4 = list(itertools.combinations(range(8), 4))
        qV = np.array([np.linalg.det(V[:, list(A)]) * np.linalg.det(V[:, sorted(set(range(8)) - set(A))]) for A in A4])
        qW = np.array([np.linalg.det(W[:, list(A)]) * np.linalg.det(W[:, sorted(set(range(8)) - set(A))]) for A in A4])
        c = qW @ qV.conj() / (qV @ qV.conj())
        worstg = max(worstg, np.abs(qW - c * qV).max() / np.abs(qW).max())
        # Gale dual of [I,U] equals conj([I,U]) up to projective equivalence: W ~ M conj(V) D
    print(f"  20 Haar U (seeds 3000-3019): max |q_A/q_0 - (-1)^(k+sum I+sum J) |det U[I,J]|^2| = {worst22:.2e}")
    print(f"  max relative deviation of q(Gale dual) from a multiple of q(V) = {worstg:.2e}")


# ---------------------------------------------------------------- part: pencil
def pencil_info(U):
    V = np.hstack([np.eye(4), U])
    mons = [(a, b) for a in range(4) for b in range(a, 4)]
    M = np.array([[V[a, i] * V[b, i] for (a, b) in mons] for i in range(8)])
    s = np.linalg.svd(M, compute_uv=False)
    _, _, Vh = np.linalg.svd(M)
    ker = Vh[8:].conj()             # 2 x 10, M @ ker.T = 0
    def sym(c):
        Q = np.zeros((4, 4), complex)
        for (a, b), x in zip(mons, c):
            if a == b: Q[a, a] = x
            else: Q[a, b] = Q[b, a] = x / 2
        return Q
    Q1, Q2 = sym(ker[0]), sym(ker[1])
    # binary quartic det(s Q1 + t Q2): sample and fit, roots in P^1
    ts = np.exp(2j * np.pi * np.arange(5) / 5)
    vals = [np.linalg.det(Q1 + t * Q2) for t in ts]
    coef = np.linalg.solve(np.vander(ts, 5), vals)   # highest degree first
    roots = np.roots(coef)
    # chordal distances between roots on P^1 (include infinity if leading coef small)
    pts = list(roots) + ([np.inf] * (4 - len(roots)))
    def chord(z, w):
        if z == np.inf and w == np.inf: return 0.0
        if z == np.inf: return 1 / np.sqrt(1 + abs(w) ** 2)
        if w == np.inf: return 1 / np.sqrt(1 + abs(z) ** 2)
        return abs(z - w) / np.sqrt((1 + abs(z) ** 2) * (1 + abs(w) ** 2))
    sep = min(chord(pts[i], pts[j]) for i, j in itertools.combinations(range(4), 2))
    # columns distinct: min chordal distance between the 8 points of P^3
    Vn = V / np.linalg.norm(V, axis=0)
    dmin = min(np.sqrt(max(0.0, 1 - abs(np.vdot(Vn[:, i], Vn[:, j])) ** 2)) for i, j in itertools.combinations(range(8), 2))
    return s[7] / s[0], s[8] / s[0] if len(s) > 8 else 0.0, sep, dmin


def part_pencil(mats):
    print("== PART pencil: quadrics through the 8 columns of [I,U] ==")
    print("  cols: sigma8/sigma1 of the 8x10 evaluation matrix (8 independent conditions iff >0),")
    print("        min chordal root separation of det(sQ1+tQ2) on P^1 (smooth base quartic iff >0),")
    print("        min chordal distance between the 8 column points (distinct iff >0)")
    for name, U in mats:
        a, _, sep, dmin = pencil_info(U)
        flag = "general position" if (a > 1e-8 and sep > 1e-6 and dmin > 1e-8) else "DEGENERATE"
        print(f"  {name:28s} s8/s1={a:.2e}  rootsep={sep:.2e}  ptdist={dmin:.2e}  {flag}")


# ---------------------------------------------------------------- part: least squares
def lm_retract(res, state, retract, n, iters=300, h=1e-7):
    """Levenberg-Marquardt on a manifold: the step delta in R^n is applied by retract(state, delta)
    (exponential map), Jacobian by central differences at delta = 0."""
    f = res(state); cost = f @ f; mu = 1e-3
    for _ in range(iters):
        J = np.empty((f.size, n))
        for m in range(n):
            e = np.zeros(n); e[m] = h
            J[:, m] = (res(retract(state, e)) - res(retract(state, -e))) / (2 * h)
        g = J.T @ f; Hs = J.T @ J
        while True:
            d = -np.linalg.solve(Hs + mu * np.diag(np.diag(Hs) + 1e-12), g)
            new = retract(state, d); fn = res(new); cn = fn @ fn
            if cn < cost:
                state, f, mu = new, fn, max(mu / 5, 1e-12)
                done = cn < 1e-30
                cost = cn; break
            mu *= 4
            if mu > 1e12: return state, f
        if done: break
    return state, f


def decompose(U, nstart, seed, tol=1e-26):
    rU = r_vec(U)
    rng = np.random.default_rng(seed)
    sols = []
    def res(st):
        A, B, m = st; lam = np.sin(m) ** 2
        return lam * r_vec(A) + (1 - lam) * r_vec(B) - rU
    def retract(st, d):
        A, B, m = st
        return (A @ expm(skew(d[:6])), B @ expm(skew(d[6:12])), m + d[12])
    for st in range(nstart):
        O1, O2 = haar_orth(rng), haar_orth(rng)
        state, f = lm_retract(res, (O1, O2, rng.uniform(0.2, 1.4)), retract, 13)
        A, B, m = state; lam = np.sin(m) ** 2
        if lam < 0.5:  # convention: lam >= 1/2 is the weight of O-
            A, B, lam = B, A, 1 - lam
        sols.append((np.abs(f).max(), lam, A, B))
    sols.sort(key=lambda t: t[0])
    return rU, sols


def cluster(sols, thr=1e-4, tol=1e-11):
    reps = []
    for res, lam, A, B in sols:
        if res > tol: continue
        ra, rb = r_vec(A), r_vec(B)
        for rep in reps:
            if np.abs(rep[0] - ra).max() < thr and np.abs(rep[1] - rb).max() < thr:
                rep[2] += 1; break
            if np.abs(rep[0] - rb).max() < thr and np.abs(rep[1] - ra).max() < thr:
                rep[2] += 1; break
        else:
            reps.append([ra, rb, 1, lam, A, B])
    return reps


def decomp_jac_sv(A, B, lam, h=1e-6):
    """Singular values of d/d(params) [lam r(A) + (1-lam) r(B)] at a solution (13 params:
    6 + 6 Lie-algebra directions and lam). Rank < 13 means a positive-dimensional solution set."""
    def G(d):
        return (lam + d[12]) * r_vec(A @ expm(skew(d[:6]))) + (1 - lam - d[12]) * r_vec(B @ expm(skew(d[6:12])))
    J = np.empty((70, 13))
    for m in range(13):
        e = np.zeros(13); e[m] = h
        J[:, m] = (G(e) - G(-e)) / (2 * h)
    return np.linalg.svd(J, compute_uv=False)


def tangent_check(U, d):
    """Is the segment direction d tangent to r(U(4)) at U?  J = dr(e^{iH}U)/dH."""
    h = 1e-6; J = np.empty((70, 16))
    for m in range(16):
        e = np.zeros(16); e[m] = h
        J[:, m] = (r_vec(expm(1j * herm(e)) @ U) - r_vec(expm(-1j * herm(e)) @ U)) / (2 * h)
    u, s, vt = np.linalg.svd(J)
    rank = int((s > 1e-7 * s[0]).sum())
    Ub = u[:, :rank]
    resid = np.linalg.norm(d - Ub @ (Ub.T @ d)) / np.linalg.norm(d)
    return rank, resid


def part_ls(mats, nstart_rand=6, nstart_struct=12):
    print("== PART ls: multistart least squares for (O-, O+, lam), lam = sin^2(mu) in [0,1] ==")
    print("  residual = max_{I,J} |lam r(O-) + (1-lam) r(O+) - r(U)|; convergent = residual < 1e-12")
    print("  nclass = number of distinct unordered pairs {r(O-), r(O+)} (up to 1e-4) among convergent (residual < 1e-11) starts")
    rows = []
    for idx, (name, U) in enumerate(mats):
        ns = nstart_rand if name.startswith("haar") else nstart_struct
        t0 = time.time()
        rU, sols = decompose(U, ns, 5000 + idx)
        reps = cluster(sols)
        best = sols[0]
        nconv = sum(1 for s in sols if s[0] < 1e-12)
        line = f"  {name:28s} best={best[0]:.1e} conv={nconv}/{ns} nclass={len(reps)}"
        sub = []
        for rp in reps:
            d = rp[0] - rp[1]
            rank, tres = tangent_check(U, d) if np.abs(d).max() > 1e-9 else (0, float("nan"))
            rp.append(tres)
            sub.append(f"      class: starts={rp[2]:2d} lam={rp[3]:.6f} |r(O-)-r(O+)|={np.abs(d).max():.3e}"
                       f" rank(dr)={rank} tangent_resid={tres:.1e}")
        order = sorted(range(len(reps)), key=lambda i: (np.nan_to_num(reps[i][6], nan=9.0)))
        reps = [reps[i] for i in order]; sub = [sub[i] for i in order]
        if reps:
            sv = decomp_jac_sv(reps[0][4], reps[0][5], reps[0][3])
            line += f" jac_sv(13)={' '.join(f'{x:.0e}' for x in sv)}"
        line += f" ({time.time() - t0:.1f}s)"
        print(line, flush=True)
        for s_ in sub: print(s_, flush=True)
        rows.append((name, U, best, reps))
    return rows


# ---------------------------------------------------------------- part: trace along the segment
def solve_unitary(target, V0, rng=None, nstart=0):
    res = lambda V: r_vec(V) - target
    retract = lambda V, d: expm(1j * herm(d)) @ V
    best = None
    starts = [V0] + [haar_unitary(rng) for _ in range(nstart)]
    for S in starts:
        V, f = lm_retract(res, S, retract, 16)
        val = np.abs(f).max()
        if best is None or val < best[0]:
            best = (val, V)
    return best


def part_trace(rows, names):
    print("== PART trace: continuation of unitaries along the segment, and the extension test ==")
    print("  P(t) = t r(O-) + (1-t) r(O+); t=lam is U. We solve r(V)=P(t) for V in U(4) by continuation from U.")
    print("  realdef = max|Im V_ij V_kl conj(V_il V_kj)| (0 iff V is real orthogonal up to diagonal phases)")
    for name, U, best, reps in rows:
        if name not in names or not reps: continue
        ra, rb, cnt, lam, A, B = reps[0][:6]
        print(f"  -- {name}: lam={lam:.6f}, realdef(U)={realness_defect(U):.2e}")
        for t_end in (1.0, 0.0):
            V = U.copy(); worst = 0.0; mid = []
            ts = np.linspace(lam, t_end, 12)[1:]
            for t in ts:
                val, V = solve_unitary(t * ra + (1 - t) * rb, V)
                worst = max(worst, val); mid.append(realness_defect(V))
            print(f"     path lam->{t_end:.0f}: worst residual {worst:.1e}; realdef along path {mid[0]:.2e} .. {mid[-2]:.2e};"
                  f" at endpoint realdef={mid[-1]:.2e},"
                  f" |r(V)-r(O_end)|={np.abs(r_vec(V) - (ra if t_end == 1 else rb)).max():.1e}")
            rng = np.random.default_rng(7000)
            for delta in (0.02, 0.1):
                t = t_end + (delta if t_end == 1 else -delta)
                val, _ = solve_unitary(t * ra + (1 - t) * rb, V, rng, nstart=7)
                print(f"       beyond endpoint, t={t:+.3f}: best residual over 8 starts = {val:.2e}"
                      f"  (segment length scale {np.abs(ra - rb).max():.2e}, offset {delta * np.abs(ra - rb).max():.2e})")


# ---------------------------------------------------------------- part: mpmath refinement
def part_mp(rows, k=3, dps=30):
    import mpmath as mp
    mp.mp.dps = dps + 10
    print(f"== PART mp: Gauss-Newton refinement at {dps} digits for the {k} worst Haar cases ==")
    cand = sorted([r for r in rows if r[0].startswith("haar") and r[3]], key=lambda r: -r[2][0])[:k]

    def polar(X):
        X = mp.matrix(X)
        for _ in range(60):
            Xn = (X + mp.inverse(X.H)) / 2
            if mp.mnorm(Xn - X, 1) < mp.mpf(10) ** (-(dps + 8)): X = Xn; break
            X = Xn
        return X

    def rmp(V):
        out = []
        for (I, J) in IJ:
            if not I: out.append(mp.mpf(1)); continue
            S = mp.matrix([[V[i, j] for j in J] for i in I])
            out.append(abs(mp.det(S)) ** 2)
        return mp.matrix(out)

    def skewmp(a):
        S = mp.zeros(4, 4); n = 0
        for i in range(4):
            for j in range(i + 1, 4):
                S[i, j] = a[n]; S[j, i] = -a[n]; n += 1
        return S

    for name, U, best, reps in cand:
        ra, rb, cnt, lam, A, B = reps[0][:6]
        Ump = polar(U.tolist()); A0 = polar(A.tolist()); B0 = polar(B.tolist())
        rU = rmp(Ump)
        def F(x):
            Am = A0 * mp.expm(skewmp(x[:6])); Bm = B0 * mp.expm(skewmp(x[6:12])); l = x[12]
            return rmp(Am) * l + rmp(Bm) * (1 - l) - rU
        x = [mp.mpf(0)] * 12 + [mp.mpf(lam)]
        nf0 = mp.norm(F(x), mp.inf)
        h = mp.mpf(10) ** (-(dps // 2 + 8))
        eps = mp.mpf(10) ** (-(dps + 4))
        for it in range(10):
            f = F(x)
            J = mp.zeros(70, 13)
            for m in range(13):
                xp = list(x); xp[m] += h; xm = list(x); xm[m] -= h
                col = (F(xp) - F(xm)) / (2 * h)
                for i in range(70): J[i, m] = col[i]
            # damped (minimum-norm) Gauss-Newton step: J has a one-dimensional kernel
            N = J.T * J
            for m in range(13): N[m, m] += eps
            dx = mp.lu_solve(N, -(J.T * f))
            x = [x[i] + dx[i] for i in range(13)]
            nf = mp.norm(F(x), mp.inf)
            if nf < mp.mpf(10) ** (-dps): break
        sv = mp.svd_r(J, compute_uv=False)
        sv = sorted([abs(s) for s in sv])
        print(f"  {name}: double-precision start residual {mp.nstr(nf0, 3)} (U re-unitarised at {dps + 10} digits)"
              f" -> residual {mp.nstr(nf, 3)} after {it + 1} steps; lam={mp.nstr(x[12], 20)};"
              f" Jacobian sigma_min, sigma_2/sigma_max = {mp.nstr(sv[0] / sv[-1], 3)}, {mp.nstr(sv[1] / sv[-1], 3)}")
        Am = A0 * mp.expm(skewmp(x[:6])); Bm = B0 * mp.expm(skewmp(x[6:12]))
        print(f"     orthogonality defects |O^tO - I|: {mp.nstr(mp.mnorm(Am.T * Am - mp.eye(4), 1), 3)},"
              f" {mp.nstr(mp.mnorm(Bm.T * Bm - mp.eye(4), 1), 3)}; lam in [0,1]: {0 <= x[12] <= 1}")


# ---------------------------------------------------------------- part: theta functions
def part_theta():
    """Numerical check of the explicit formulas of the own proof (theorem-1-2-second.md, Sec. A)."""
    import mpmath as mp
    mp.mp.dps = 30
    print("== PART theta: Frobenius-Stickelberger form of Delta_A, the line, iota, and a real circle ==")
    rng = np.random.default_rng(9000)
    tau = mp.mpc(0.13, 1.1)
    q = mp.exp(1j * mp.pi * tau)
    th = lambda z: mp.jtheta(1, mp.pi * z, q)
    u = [mp.mpc(*rng.uniform(0, 1, 2)) for _ in range(8)]
    Sig = sum(u)
    A4 = list(itertools.combinations(range(8), 4))

    def config(s, seed):
        g = np.random.default_rng(seed)
        rows = []
        for k in range(4):
            a = [mp.mpc(*g.uniform(0, 1, 2)) for _ in range(3)]
            a.append(s - sum(a))
            rows.append([mp.fprod(th(z - aa) for aa in a) for z in u])
        return mp.matrix(rows)

    def Delta(V, A): return mp.det(mp.matrix([[V[k, i] for i in A] for k in range(4)]))

    worst = 0
    for t, s in enumerate([mp.mpc(0.3, 0.2), mp.mpc(-0.1, 0.7), mp.mpc(0.45, -0.33)]):
        V = config(s, 9100 + t)
        rat = []
        for A in A4:
            pred = th(sum(u[i] for i in A) - s) * mp.fprod(th(u[i] - u[j]) for i, j in itertools.combinations(A, 2))
            rat.append(Delta(V, A) / pred)
        worst = max(worst, max(abs(x / rat[0] - 1) for x in rat))
    print(f"  Delta_A(s) = C(s) theta1(u_A - s) prod_(i<j in A) theta1(u_i-u_j): max rel. spread of C over 70 A, 3 s: {mp.nstr(worst, 3)}")

    S = [mp.mpc(*rng.uniform(-1, 1, 2)) for _ in range(6)]
    G = []
    for t, s in enumerate(S):
        V = config(s, 9200 + t)
        row = [Delta(V, A) * Delta(V, tuple(sorted(set(range(8)) - set(A)))) for A in A4]
        nrm = mp.norm(mp.matrix(row)); G.append([x / nrm for x in row])
    sv = mp.svd_c(mp.matrix(G), compute_uv=False)
    print("  singular values of the 6x70 matrix of normalised q-vectors q(s_1..s_6):",
          ", ".join(mp.nstr(x, 3) for x in sv), " (rank 2 = a projective line)")
    s = mp.mpc(0.21, 0.37)
    V1, V2 = config(s, 9300), config(Sig - s, 9301)
    q1 = mp.matrix([Delta(V1, A) * Delta(V1, tuple(sorted(set(range(8)) - set(A)))) for A in A4])
    q2 = mp.matrix([Delta(V2, A) * Delta(V2, tuple(sorted(set(range(8)) - set(A)))) for A in A4])
    c = q2[0] / q1[0]
    print(f"  iota: max |q(Sigma-s) - c q(s)| / |q| = {mp.nstr(mp.norm(q2 - c * q1, mp.inf) / mp.norm(q2, mp.inf), 3)}")

    # ---- real circle demo: rectangular lattice, eight real points, s = Sigma/2 + i y
    print("  real-structure demo: tau = 1.2 i, eight real points (seed 9400), s on the circle Re s = Re(Sigma)/2")
    tau = mp.mpc(0, 1.2); q = mp.exp(1j * mp.pi * tau)
    th = lambda z: mp.jtheta(1, mp.pi * z, q)
    g = np.random.default_rng(9400)
    u = [mp.mpf(x) for x in np.sort(g.uniform(0, 1, 8))]
    Sig = sum(u)
    for half in (0, mp.mpf(1) / 2):
        pats = []; reald = []; vals = []
        for y in np.linspace(0, 1.2, 25)[:-1]:
            s = Sig / 2 + half + 1j * mp.mpf(y)
            a = [mp.mpf(0.1), mp.mpf(0.37), mp.mpf(0.71)]
            rowsV = []
            for k in range(4):
                b = [x + mp.mpf(0.05) * k for x in a] + [None]
                b[3] = s - sum(b[:3])
                rowsV.append([mp.fprod(th(z - bb) for bb in b) for z in u])
            V = np.array([[complex(rowsV[k][i]) for i in range(8)] for k in range(4)])
            # real-structure relation V D V^* = 0 with D diagonal: 16 complex equations, 8 unknowns
            Mx = np.array([[V[a_, i] * np.conj(V[b_, i]) for i in range(8)] for a_ in range(4) for b_ in range(4)])
            _, sv_, vh = np.linalg.svd(Mx)
            d = vh[-1].conj(); d = d / d[0]
            pats.append("".join("+" if x.real > 0 else "-" for x in d))
            reald.append((np.abs(d.imag).max(), sv_[-1] / sv_[0], sv_[-2] / sv_[0]))
            Z = np.linalg.solve(V[:, :4], V[:, 4:])
            cr = max(abs((Z[i, j] * Z[k, l] / (Z[i, l] * Z[k, j])).imag)
                     for i, k in itertools.combinations(range(4), 2) for j, l in itertools.combinations(range(4), 2))
            vals.append((y, cr))
        print(f"   circle Re s = Re Sigma/2 + {float(half)}: sign patterns of D along the circle: {sorted(set(pats))};"
              f" max Im(d)={max(x[0] for x in reald):.1e}; kernel dims ok (sigma8/sigma1 max {max(x[1] for x in reald):.1e},"
              f" sigma7/sigma1 min {min(x[2] for x in reald):.1e})")
        print("    realness (max |Im| of cross-ratios of Z) at y=0 (ramification), y=0.6 (ramification), y=0.25:",
              f"{vals[0][1]:.1e}, {vals[12][1]:.1e}, {vals[5][1]:.1e}")


# ---------------------------------------------------------------- main
def main():
    part = sys.argv[1] if len(sys.argv) > 1 else "all"
    mats = test_matrices()
    if part in ("trace1", "mp1"):   # one matrix per run: argv[2] = test-matrix name
        sel = [m for m in mats if m[0] == sys.argv[2]]
        rows = part_ls(sel)
        if part == "trace1": part_trace(rows, {sys.argv[2]})
        else: part_mp(rows, k=1)
        return
    if part == "quick":
        part_ls(mats[:3], nstart_rand=6); return
    if part in ("sign", "all"): part_sign()
    if part in ("pencil", "all"): part_pencil(mats)
    if part in ("ls", "trace", "mp", "all"):
        sel = mats
        if part == "trace":
            keep = {"haar[1000]", "haar[1001]", "haar[1002]", "fourier4", "onezero[2004]", "twozeros[2003]",
                    "blockdiag3+1[2002]", "nearperm eps=0.1[2005]", "nearorth eps=0.1[2006]"}
            sel = [m for m in mats if m[0] in keep]
        if part == "mp":
            sel = [m for m in mats if m[0].startswith("haar")]
        rows = part_ls(sel)
        if part in ("trace", "all"):
            part_trace(rows, {m[0] for m in sel if m[0] in {"haar[1000]", "haar[1001]", "haar[1002]", "fourier4",
                              "onezero[2004]", "twozeros[2003]", "blockdiag3+1[2002]", "nearperm eps=0.1[2005]",
                              "nearorth eps=0.1[2006]"}})
        if part in ("mp", "all"): part_mp(rows)
    if part in ("theta", "all"): part_theta()


if __name__ == "__main__":
    main()
