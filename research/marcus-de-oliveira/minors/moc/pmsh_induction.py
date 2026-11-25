"""Test candidate induction steps for PMSH: does  sum_j p_j m(M_j) = m(H')  and  sum_j p_j c_j m(M_j) = m_K
hold for canonical choices of M_j (spectrum c minus c_j), p_j = |U_nj|^2 ?"""
import sys, itertools, numpy as np
sys.path.insert(0, 'moc')
from moclib import haar_unitary
def minors(M):
    k = M.shape[0]; return np.array([np.linalg.det(M[np.ix_(T, T)]).real for r in range(1, k + 1) for T in itertools.combinations(range(k), r)])
def householder_to_en(v):
    n = len(v); th = np.angle(v[-1]) if abs(v[-1]) > 0 else 0.0
    e = np.zeros(n, complex); e[-1] = np.exp(1j * th)
    u = v - e
    if np.linalg.norm(u) < 1e-14: return np.eye(n, dtype=complex) * 1
    R = np.eye(n) - 2 * np.outer(u, u.conj()) / (u.conj() @ u)
    return R          # R v = e  (unitary, Hermitian)
def min_rotation_to_en(v):
    n = len(v); th = np.angle(v[-1]); e = np.zeros(n, complex); e[-1] = np.exp(1j * th)
    # rotation in span{v, e} taking v -> e, identity on the orthogonal complement
    w = e - (v.conj() @ e) * v; nw = np.linalg.norm(w)
    if nw < 1e-14: return np.eye(n, dtype=complex)
    w /= nw; cosA = (v.conj() @ e).real; sinA = np.sqrt(max(0, 1 - cosA**2))
    P = np.outer(v, v.conj()) + np.outer(w, w.conj())
    Rplane = cosA * (np.outer(v, v.conj()) + np.outer(w, w.conj())) + sinA * (np.outer(w, v.conj()) - np.outer(v, w.conj()))
    return np.eye(n) - P + Rplane
rng = np.random.default_rng(3)
for n in [3, 4]:
    for name, f in [('householder', householder_to_en), ('min-rotation', min_rotation_to_en)]:
        errs = []
        for _ in range(20):
            U = haar_unitary(n, rng); c = rng.uniform(0.5, 3, n); H = (U * c) @ U.conj().T
            p = np.abs(U[-1, :])**2
            G = minors(H[:-1, :-1])
            K = np.array([np.linalg.det(H[np.ix_(list(T) + [n - 1], list(T) + [n - 1])]).real for r in range(1, n) for T in itertools.combinations(range(n - 1), r)])
            sG = 0; sK = 0
            for j in range(n):
                R = f(U[:, j]); assert np.allclose(R @ U[:, j], np.exp(1j * np.angle(U[-1, j])) * np.eye(n)[-1])
                Mj = (R @ H @ R.conj().T)[:-1, :-1]
                assert np.allclose(np.sort(np.linalg.eigvalsh(Mj)), np.sort(np.delete(c, j)))
                mj = minors(Mj); sG = sG + p[j] * mj; sK = sK + p[j] * c[j] * mj
            errs.append(max(np.abs(sG - G).max(), np.abs(sK - K).max()))
        print("n=%d %-13s max error in identities (G),(K): %.2e" % (n, name, max(errs)))
