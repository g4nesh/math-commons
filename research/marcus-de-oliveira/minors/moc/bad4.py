"""n=4 targeted search.  Counterexamples can only occur for U whose compound vector violates one of
the 48 exotic facets.  Phase 1: collect deep violators.  Phase 2: for each, maximise J over (a,b,theta)
with U fixed (batched, analytic gradients).  U is never varied after phase 1: no joint refinement is implemented.
Half of the unitaries are 'deepened' in phase 1; these all converge to the one extremal orbit (see reanalyse4.py)."""
import sys, json, time, numpy as np
sys.path.insert(0, 'moc')
from bsearch import Problem, haar, normalise, expm_skew
from reduced4 import PAIRS, PERMS, VERTS
F = np.load('moc/facets4.npy')
onv = np.array([(np.abs(VERTS @ f[:-1] + f[-1]) < 1e-9) for f in F])
EX = F[onv.sum(1) == 13]                       # 48 exotic facets, inside: EX[:, :-1] x + EX[:, -1] <= 0
EXn = EX / np.linalg.norm(EX[:, :-1], axis=1, keepdims=True)

def red_batch(U):
    M1 = np.abs(U)**2
    N = np.zeros((U.shape[0], 3, 3))
    for p, (S, Sc) in enumerate(PAIRS):
        for q, (T, Tc) in enumerate(PAIRS):
            A = U[:, list(S)][:, :, list(T)]; B = U[:, list(S)][:, :, list(Tc)]
            N[:, p, q] = np.abs(np.linalg.det(A))**2 + np.abs(np.linalg.det(B))**2
    return np.concatenate([M1[:, :3, :3].reshape(-1, 9), N[:, :2, :2].reshape(-1, 4)], 1)

def viol(U):
    x = red_batch(U)
    return (x @ EXn[:, :-1].T + EXn[:, -1]).max(1)   # >0 : outside polytope (Euclidean, reduced coords)

def deepen(U, steps=300, lr=0.05, rng=None):
    """gradient ascent on violation via finite differences in the 16-dim Lie algebra"""
    n = 4; basis = []
    for i in range(n):
        H = np.zeros((n, n), complex); H[i, i] = 1j; basis.append(H)
    for i in range(n):
        for j in range(i + 1, n):
            H = np.zeros((n, n), complex); H[i, j] = 1; H[j, i] = -1; basis.append(H)
            H = np.zeros((n, n), complex); H[i, j] = 1j; H[j, i] = 1j; basis.append(H)
    basis = np.array(basis); eps = 1e-6
    for t in range(steps):
        v0 = viol(U)
        g = np.zeros((U.shape[0], 16))
        for k in range(16):
            Ep = expm_skew(np.broadcast_to(eps * basis[k], U.shape).copy())
            g[:, k] = (viol(U @ Ep) - v0) / eps
        Hs = np.einsum('nk,kij->nij', g, basis)
        nrm = np.linalg.norm(g, axis=1)[:, None, None] + 1e-30
        U = U @ expm_skew(lr * (1 - t / steps) * Hs / nrm)
    return U

if __name__ == '__main__':
    seed = int(sys.argv[1]); nU = int(sys.argv[2]); N = int(sys.argv[3]); steps = int(sys.argv[4])
    rng = np.random.default_rng(seed)
    pool = []
    while sum(len(p) for p in pool) < nU:
        U = haar(20000, 4, rng); v = viol(U); pool.append(U[v > 0])
    U0 = np.concatenate(pool)[:nU]
    v_before = viol(U0)
    Ud = deepen(U0.copy())
    v_after = viol(Ud)
    print("phase1: %d bad U; violation before: max %.4f; after deepening: max %.4f median %.4f" %
          (nU, v_before.max(), v_after.max(), np.median(v_after)), flush=True)
    # mix: half deepened, half raw
    Us = np.concatenate([Ud[: nU // 2], U0[nU // 2:]])
    pr = Problem(4)
    results = []
    t0 = time.time()
    for u_idx in range(nU):
        U = np.broadcast_to(Us[u_idx], (N, 4, 4)).copy()
        a = rng.standard_normal((N, 4)) + 1j * rng.standard_normal((N, 4))
        b = rng.standard_normal((N, 4)) + 1j * rng.standard_normal((N, 4))
        a, b = normalise(a, b); th = rng.uniform(0, 2 * np.pi, N)
        ma = np.zeros_like(a); mb = np.zeros_like(b); mt = np.zeros(N)
        for t in range(1, steps + 1):
            frac = t / steps; tau = 0.05 * (1 - frac) + 0.0005 * frac
            J, (gU, gA, gB, gth) = pr.evaluate(U, a, b, th, tau)
            gn = np.sqrt((np.abs(gA)**2).sum(1) + (np.abs(gB)**2).sum(1) + gth**2) + 1e-30
            ma = 0.9 * ma + 0.1 * gA / gn[:, None]; mb = 0.9 * mb + 0.1 * gB / gn[:, None]; mt = 0.9 * mt + 0.1 * gth / gn
            st = 0.03 * (1 - 0.9 * frac)
            a, b, th = a + st * ma, b + st * mb, th + st * mt
            a, b = normalise(a, b)
        Je, _ = pr.evaluate(U, a, b, th, None, grad=False)
        k = int(np.nanargmax(Je))
        results.append({'viol': float(viol(Us[u_idx:u_idx+1])[0]), 'J': float(Je[k]),
                         'U_re': Us[u_idx].real.tolist(), 'U_im': Us[u_idx].imag.tolist(),
                         'a': [a[k].real.tolist(), a[k].imag.tolist()], 'b': [b[k].real.tolist(), b[k].imag.tolist()], 'th': float(th[k])})
        if (u_idx + 1) % 10 == 0 or u_idx == nU - 1:
            Js = np.array([r['J'] for r in results])
            print("  %d/%d U done, max J so far %.4e, median of per-U maxima %.4e  (%.0fs)" %
                  (u_idx + 1, nU, Js.max(), np.median(Js), time.time() - t0), flush=True)
            json.dump(results, open('moc/runs/bad4_s%d.json' % seed, 'w'))
