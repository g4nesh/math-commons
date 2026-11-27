"""Fixed-U support-gap search over structured unitary families, n = 5..8.

Usage: python -I -B structured_search.py SEED SECONDS_PER_FAMILY OUT.json FAMILY [FAMILY ...]
       python -I -B structured_search.py --selftest
       python -I -B structured_search.py --list

For each family U(p) (p = family parameters, possibly none) we maximize the scale- and shift-invariant gap
    h(a, b, theta, p) = [Re(e^{-i theta} d) - max_sigma Re(e^{-i theta} z_sigma)] / S,
d = det(D_a + U D_b U*), z_sigma = prod_i (a_i + b_sigma(i)), S = rho^n with rho^2 = (sum|a - mean a|^2 +
sum|b - mean b|^2) / 2n (so h is invariant under (a, b) -> (c(a + t), c(b - t)), c > 0).
h > 0 for some theta iff d lies outside conv{z_sigma}. The max is replaced by a log-sum-exp at relative
temperature tau (annealed 5e-2 -> 1e-4; LSE >= max, so the smoothed value is a lower bound), with an
analytic gradient in (a, b, theta) and central differences in p. L-BFGS-B, restart r of family f uses
default_rng([SEED, f, r]); even restarts are free, odd restarts carry a quadratic penalty keeping all
|a_i - a_j|, |b_i - b_j|, |a_i + b_j| >= DELTA * rms (non-degenerate spectra). The final value is the
robust edge margin m (lift_lib.robust_margin, sign-correct), and every record with m > -1e-9 (up to 3 per
family and mode) is recomputed at 50 digits with moc_hp.hp_margin.
"""
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
OPT = os.path.join(HERE, "..", "optimizer-search")
sys.path.insert(0, os.path.join(OPT, "lib"))
sys.path.insert(0, os.path.join(OPT, "lifted"))
import numpy as np
from scipy.optimize import minimize
from moc_core import perms
from moc_param import random_spectra, random_unitary
from moc_hp import hp_margin
from lift_lib import robust_margin

DELTA = 0.1
TAUS = (5e-2, 1e-2, 2e-3, 4e-4, 1e-4)

# ---------------------------------------------------------------- families
OSTAR = np.array([[1, 4, 4, 4], [-4, 1, 4, -4], [-4, -4, 1, 4], [-4, 4, -4, 1]], float) / 7


def _vmat():  # orthogonal-hull-five/exact_ustar.py V, entries (a + c i sqrt15)/20 and 4 sqrt6/20
    E = lambda a, c=0: (a + c * 1j * np.sqrt(15)) / 20
    s = 4 * np.sqrt(6) / 20
    return np.array([[E(4), s, s, s, s],
                     [s, E(-1, -1), E(-6, -2), E(-6, 2), E(9, 1)],
                     [s, E(-6, -2), E(-1, 1), E(9, -1), E(-6, 2)],
                     [s, E(-6, 2), E(9, -1), E(-1, 1), E(-6, -2)],
                     [s, E(9, 1), E(-6, 2), E(-6, -2), E(-1, -1)]])


V5 = _vmat()


def F(n):
    j = np.arange(n)
    return np.exp(2j * np.pi * np.outer(j, j) / n) / np.sqrt(n)


def dsum(*Bs):
    n = sum(len(B) for B in Bs)
    U = np.zeros((n, n), complex)
    k = 0
    for B in Bs:
        U[k:k + len(B), k:k + len(B)] = B
        k += len(B)
    return U


def R2(p):  # general SU(2)-type rotation with a phase
    c, s, e = np.cos(p[0]), np.sin(p[0]), np.exp(1j * p[1])
    return np.array([[c, -s * e], [s / e, c]])


def dita(A, B):  # (F_A (x) I_B) blockdiag(diag(e^{i p_k}) F_B): Dita construction, phases p in R^(A*B)
    K = np.kron(F(A), np.eye(B))
    return lambda p: K @ dsum(*[np.diag(np.exp(1j * p[k * B:(k + 1) * B])) @ F(B) for k in range(A)])


def const(U):
    return lambda p: U


FAMILIES = {  # name: (n, U_of(p), number of parameters)
    "F5": (5, const(F(5)), 0),
    "Ostar+1": (5, const(dsum(OSTAR, np.eye(1))), 0),
    "V": (5, const(V5), 0),
    "F6": (6, const(F(6)), 0),
    "Dita2x3": (6, dita(2, 3), 6),
    "Dita3x2": (6, dita(3, 2), 6),
    "Ostar+I2": (6, const(dsum(OSTAR, np.eye(2))), 0),
    "Ostar+F2": (6, const(dsum(OSTAR, F(2))), 0),
    "Ostar+R2": (6, lambda p: dsum(OSTAR, R2(p)), 2),
    "V+1": (6, const(dsum(V5, np.eye(1))), 0),
    "F7": (7, const(F(7)), 0),
    "Ostar+F3": (7, const(dsum(OSTAR, F(3))), 0),
    "V+F2": (7, const(dsum(V5, F(2))), 0),
    "V+I2": (7, const(dsum(V5, np.eye(2))), 0),
    "F8": (8, const(F(8)), 0),
    "Dita2x4": (8, dita(2, 4), 8),
    "Dita4x2": (8, dita(4, 2), 8),
    "Ostar+Ostar": (8, const(dsum(OSTAR, OSTAR)), 0),
    "Ostar+F4": (8, const(dsum(OSTAR, F(4))), 0),
    "Ostar+I4": (8, const(dsum(OSTAR, np.eye(4))), 0),
    "OstarxR2": (8, lambda p: np.kron(OSTAR, R2(p)), 2),
    "R2xOstar": (8, lambda p: np.kron(R2(p), OSTAR), 2),
    "V+F3": (8, const(dsum(V5, F(3))), 0),
    "V+I3": (8, const(dsum(V5, np.eye(3))), 0),
}


# ---------------------------------------------------------------- objective
def unpack(x, n):
    a = x[:n] + 1j * x[n:2 * n]
    b = x[2 * n:3 * n] + 1j * x[3 * n:4 * n]
    return a, b, x[4 * n], x[4 * n + 1:]


def h_and_grad(x, n, U_of, tau, P):
    a, b, th, p = unpack(x, n)
    U = U_of(p)
    c = np.exp(-1j * th)
    M = np.diag(a) + U @ np.diag(b) @ U.conj().T
    d = np.linalg.det(M)
    Mi = np.linalg.inv(M)
    da = d * np.diag(Mi)
    db = d * np.einsum("ji,jk,ki->i", U.conj(), Mi, U)  # d (U* M^-1 U)_jj
    z = (a[None, :] + b[P]).prod(axis=1)
    rho2 = (np.sum(np.abs(a - a.mean()) ** 2) + np.sum(np.abs(b - b.mean()) ** 2)) / (2 * n)
    S = rho2 ** (n / 2)
    r = (c * z).real
    rd = (c * d).real
    t = S * tau
    e = np.exp((r - r.max()) / t)
    w = e / e.sum()
    h = rd / S - (r.max() / S + tau * np.log(e.sum()))
    sw = (w * r).sum()
    idx = np.flatnonzero(w > 1e-18)  # the rest contribute below rounding
    Ps = P[idx]
    wQ = (w[idx] * z[idx])[:, None] / (a[None, :] + b[Ps])
    sb = (np.bincount(Ps.ravel(), wQ.real.ravel(), minlength=n)
          + 1j * np.bincount(Ps.ravel(), wQ.imag.ravel(), minlength=n))  # sum over sigma(i) = j
    gSa = S * (a - a.mean()) / (2 * rho2)
    gSb = S * (b - b.mean()) / (2 * rho2)
    ga = np.conj(c * (da - wQ.sum(axis=0))) / S - (rd - sw) * gSa / S ** 2
    gb = np.conj(c * (db - sb)) / S - (rd - sw) * gSb / S ** 2
    gth = ((c * d).imag - (w * (c * z).imag).sum()) / S
    gp = []
    for k in range(len(p)):
        ep = np.zeros(len(p)); ep[k] = 1e-6
        xp, xm = x.copy(), x.copy()
        xp[4 * n + 1:] += ep; xm[4 * n + 1:] -= ep
        dp = np.linalg.det(np.diag(a) + U_of(p + ep) @ np.diag(b) @ U_of(p + ep).conj().T)
        dm = np.linalg.det(np.diag(a) + U_of(p - ep) @ np.diag(b) @ U_of(p - ep).conj().T)
        gp.append(((c * dp).real - (c * dm).real) / 2e-6 / S)
    g = np.concatenate([ga.real, ga.imag, gb.real, gb.imag, [gth], gp])
    return h, g


def seps(x, n):
    a, b, _, _ = unpack(x, n)
    mb = b.mean()
    rho = np.sqrt((np.sum(np.abs(a + mb - (a + mb).mean()) ** 2) + np.sum(np.abs(b - mb) ** 2)) / (2 * n))
    ia, ja = np.triu_indices(n, 1)
    return (np.abs(a[ia] - a[ja]).min() / rho, np.abs(b[ia] - b[ja]).min() / rho,
            np.abs(a[:, None] + b[None, :]).min() / rho)


def penalty(x, n):
    a, b, _, _ = unpack(x, n)
    mb = b.mean()
    rho = np.sqrt((np.sum(np.abs(a + mb - (a + mb).mean()) ** 2) + np.sum(np.abs(b - mb) ** 2)) / (2 * n))
    ia, ja = np.triu_indices(n, 1)
    g = np.concatenate([np.abs(a[ia] - a[ja]), np.abs(b[ia] - b[ja]), np.abs(a[:, None] + b[None, :]).ravel()]) / rho
    return float(np.sum(np.maximum(0.0, DELTA - g) ** 2))


def penalty_grad(x, n):
    g = np.zeros_like(x)
    for k in range(4 * n):
        e = np.zeros_like(x); e[k] = 1e-7
        g[k] = (penalty(x + e, n) - penalty(x - e, n)) / 2e-7
    return g


def normalize(x, n):  # shift mean b to 0 and joint rms to 1 (no rotation, so theta stays valid); h is invariant
    a, b, th, p = unpack(x, n)
    t = b.mean(); a, b = a + t, b - t
    rho = np.sqrt((np.sum(np.abs(a - a.mean()) ** 2) + np.sum(np.abs(b) ** 2)) / (2 * n))
    a, b = a / rho, b / rho  # positive real scaling: d and z scale by rho^-n > 0, theta unchanged
    return np.concatenate([a.real, a.imag, b.real, b.imag, [th], p])


def exact_eval(x, n, U_of, P):
    a, b, th, p = unpack(x, n)
    U = U_of(p)
    M = np.diag(a) + U @ np.diag(b) @ U.conj().T
    d = np.linalg.det(M)
    z = (a[None, :] + b[P]).prod(axis=1)
    m = float(robust_margin(d, z))
    r = (np.exp(-1j * th) * z).real
    S = ((np.sum(np.abs(a - a.mean()) ** 2) + np.sum(np.abs(b - b.mean()) ** 2)) / (2 * n)) ** (n / 2)
    hexact = float(((np.exp(-1j * th) * d).real - r.max()) / S)
    diam = np.ptp(z.real) + np.ptp(z.imag)
    active = np.flatnonzero(r >= r.max() - 1e-7 * diam)
    return m, hexact, active, U


def run_one(n, U_of, k, rng, pen, P, tlim):
    a, b = random_spectra(n, rng, ["gaussian", "circle", "mixed"][int(rng.integers(3))])
    p0 = rng.uniform(-np.pi, np.pi, k)
    U = U_of(p0)
    d = np.linalg.det(np.diag(a) + U @ np.diag(b) @ U.conj().T)
    z = (a[None, :] + b[P]).prod(axis=1)
    th0 = np.angle(d - z.mean())
    x = normalize(np.concatenate([a.real, a.imag, b.real, b.imag, [th0], p0]), n)
    for tau in TAUS:
        if time.time() > tlim:
            break

        def f(y):
            h, g = h_and_grad(y, n, U_of, tau, P)
            if pen:
                return -h + 100 * penalty(y, n), -g + 100 * penalty_grad(y, n)
            return -h, -g
        with np.errstate(all="ignore"):
            res = minimize(f, x, jac=True, method="L-BFGS-B", options=dict(maxiter=250))
        if np.all(np.isfinite(res.x)):
            x = normalize(res.x, n)
    return x


def record(name, n, x, U_of, P, pen, restart):
    m, hx, active, U = exact_eval(x, n, U_of, P)
    a, b, th, p = unpack(x, n)
    sa, sb, sab = seps(x, n)
    return dict(family=name, n=n, restart=restart, penalized=bool(pen), m=m, h_exact=hx, theta=float(th),
                n_active=int(len(active)), active=[list(map(int, perms(n)[i])) for i in active[:8]],
                sep_a=float(sa), sep_b=float(sb), sep_ab=float(sab), params=[float(v) for v in p],
                a=[[float(v.real), float(v.imag)] for v in a], b=[[float(v.real), float(v.imag)] for v in b],
                unitarity_defect=float(np.abs(U @ U.conj().T - np.eye(n)).max()))


def search(seed, secs, names):
    out = {}
    for fi, name in enumerate(names):
        n, U_of, k = FAMILIES[name]
        P = np.asarray(perms(n), np.intp)
        t0 = time.time(); tlim = t0 + secs
        recs = []
        r = 0
        while time.time() < tlim:
            rng = np.random.default_rng([seed, list(FAMILIES).index(name), r])
            pen = r % 2 == 1
            x = run_one(n, U_of, k, rng, pen, P, tlim)
            recs.append(record(name, n, x, U_of, P, pen, r))
            r += 1
        for pen in (False, True):
            sub = sorted([t for t in recs if t["penalized"] == pen], key=lambda t: -t["m"])
            for t in [t for t in sub if t["m"] > -1e-9][:3]:
                a = np.array([complex(*v) for v in t["a"]]); b = np.array([complex(*v) for v in t["b"]])
                hp = hp_margin(a, b, U_of(np.array(t["params"])), dps=50)
                t["hp50_m"] = hp["m_str"]
        free = sorted([t for t in recs if not t["penalized"]], key=lambda t: -t["m"])
        pend = sorted([t for t in recs if t["penalized"]], key=lambda t: -t["m"])
        feas = [t for t in pend if min(t["sep_a"], t["sep_b"], t["sep_ab"]) >= 0.9 * DELTA]
        out[name] = dict(n=n, params=k, restarts=len(recs), seconds=round(time.time() - t0, 1),
                         n_free_above_m1e9=sum(t["m"] > -1e-9 for t in free),
                         best_free=free[:2], best_penalized_feasible=feas[:2],
                         m_quantiles_free=np.quantile([t["m"] for t in free], [0, .5, .9, 1]).tolist() if free else [])
        bf = free[0]["m"] if free else None
        bp = feas[0]["m"] if feas else None
        print("%-12s n=%d R=%3d best_free=%s best_pen=%s" % (name, n, len(recs), bf, bp), flush=True)
    return out


def selftest():
    rng = np.random.default_rng(7)
    for name in ("Dita2x3", "Ostar+R2", "V"):
        n, U_of, k = FAMILIES[name]
        Up = U_of(rng.normal(size=k))
        assert np.abs(Up @ Up.conj().T - np.eye(n)).max() < 1e-12, name
        P = np.asarray(perms(n), np.intp)
        x = rng.normal(size=4 * n + 1 + k)
        h, g = h_and_grad(x, n, U_of, 0.05, P)
        for j in range(len(x)):
            e = np.zeros_like(x); e[j] = 1e-6
            fd = (h_and_grad(x + e, n, U_of, 0.05, P)[0] - h_and_grad(x - e, n, U_of, 0.05, P)[0]) / 2e-6
            assert abs(fd - g[j]) < 1e-5 * max(1, abs(fd)), (name, j, fd, g[j])
        assert abs(h_and_grad(normalize(x, n), n, U_of, 0.05, P)[0] - h) < 1e-9
    for name, (n, U_of, k) in FAMILIES.items():
        U = U_of(np.zeros(k))
        assert np.abs(U @ U.conj().T - np.eye(n)).max() < 1e-12, name
    # U = I gives d = z_identity, which is in the hull
    a, b = random_spectra(5, rng, "gaussian")
    x = np.concatenate([a.real, a.imag, b.real, b.imag, [0.0]])
    m, _, _, _ = exact_eval(x, 5, const(np.eye(5)), np.asarray(perms(5), np.intp))
    assert m < 1e-12
    print("selftest ok")


if __name__ == "__main__":
    if sys.argv[1] == "--selftest":
        selftest()
    elif sys.argv[1] == "--list":
        print(" ".join(FAMILIES))
    else:
        seed, secs, path = int(sys.argv[1]), float(sys.argv[2]), sys.argv[3]
        res = search(seed, secs, sys.argv[4:])
        json.dump(dict(seed=seed, seconds_per_family=secs, delta=DELTA, taus=TAUS, families=res),
                  open(path, "w"), indent=1)
