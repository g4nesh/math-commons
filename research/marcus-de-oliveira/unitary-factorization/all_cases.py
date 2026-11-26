"""Theorem F on structured unitaries: every branch of Lemmas 3 and 5 of README.md, in floating point.

construct.py covers the generic branch (all off-diagonal entries of S non-zero, proper triangle).  This script implements
the whole case analysis: the explicit phases of the remark after Lemma 3 (which need no case distinction for degenerate
triangles) and the four graph cases of Lemma 5 (six edges, triangle, one edge, no edge), and runs it on permutation, monomial, diagonal and real orthogonal matrices, direct sums, the Fourier
matrix, real and complex Hadamard matrices, tensor products and matrices with prescribed zeros.  Entries below TOL in
modulus are treated as zero.  For each matrix it checks D1 U D2 = O W with O real orthogonal and the squared-minor formula
of Corollary F1.  A matrix for which the case list of Lemma 5 did not apply would raise an error."""
import itertools, sys, numpy as np
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0); TOL = 1e-9
P2 = list(itertools.combinations(range(4), 2)); IDX = {S: k for k, S in enumerate(P2)}
MATCH = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
def compound(U): return np.array([[np.linalg.det(U[np.ix_(S, T)]) for T in P2] for S in P2])
def csqrt(z): return np.sqrt(z + 0j)
def triangle_kind(c):
    """shape of the triangle with side lengths |c_M| (for the report only)"""
    a, b, g = np.sort(np.abs(c))[::-1]
    return 'all sides zero' if a < TOL else 'one side zero' if g < TOL else 'flat triangle' if a - b - g > -1e-12 else 'proper triangle'
def lemma4(S, idx, x):
    """S[idx, idx] = diag(Lam) + x x^T with all x_j non-zero: returns Omega (on idx), mu, y (on idx)"""
    Lam = np.array([S[j, j] for j in idx]) - x * x; h = Lam * np.conj(x) / x
    assert np.abs(h - h[0]).max() < 1e-6 and abs(abs(h[0]) - 1) < 1e-6, "Lemma 4 fails"
    h = h.mean(); h = h / abs(h); n2 = (np.abs(x)**2).sum()
    return csqrt(h) * x / np.abs(x), n2 * np.conj(h), np.abs(x) / np.sqrt(n2)
def lemma5(S):
    """symmetric unitary S with S12 S34 = S13 S24 = S14 S23: returns Omega, mu, y, branch with S = Omega (I + mu y y^T) Omega"""
    edges = [(j, k) for j, k in P2 if abs(S[j, k]) >= TOL]; Om = np.ones(4, complex); y = np.zeros(4); mu = 0.0
    if len(edges) == 6:
        x1 = csqrt(S[0, 1] * S[0, 2] / S[1, 2]); x = np.array([x1, S[0, 1] / x1, S[0, 2] / x1, S[0, 3] / x1]); Om, mu, y = lemma4(S, range(4), x); branch = 'six edges'
    elif len(edges) == 3 and len({v for e in edges for v in e}) == 3:
        i, j, k = sorted({v for e in edges for v in e}); m = ({0, 1, 2, 3} - {i, j, k}).pop()
        xi_ = csqrt(S[i, j] * S[i, k] / S[j, k]); x = np.array([xi_, S[i, j] / xi_, S[i, k] / xi_]); o3, mu, y3 = lemma4(S, (i, j, k), x)
        Om[[i, j, k]] = o3; y[[i, j, k]] = y3; Om[m] = csqrt(S[m, m]); branch = 'triangle'
    elif len(edges) == 1:
        j, k = edges[0]; a = S[j, j]; dl = S[j, j] * S[k, k] - S[j, k]**2
        if abs(a) >= TOL: d1 = csqrt(np.conj(a) / abs(a)); d2 = csqrt(-np.conj(dl * np.conj(a) / abs(a))); branch = 'one edge'
        else: d1 = 1.0; d2 = csqrt(-np.conj(dl)); branch = 'one edge, zero diagonal entry'
        G = np.array([[d1 * d1 * S[j, j], d1 * d2 * S[j, k]], [d1 * d2 * S[j, k], d2 * d2 * S[k, k]]])        # symmetric unitary with eigenvalue 1
        w, Q = np.linalg.eigh(G.real + 0.618 * G.imag); lam = np.array([Q[:, t] @ G @ Q[:, t] for t in range(2)]); t = int(np.argmax(np.abs(lam - 1)))
        assert abs(lam[1 - t] - 1) < 1e-6, "no eigenvalue 1"
        mu = lam[t] - 1; y[[j, k]] = Q[:, t]; Om[j] = 1 / d1; Om[k] = 1 / d2
        for m in set(range(4)) - {j, k}: Om[m] = csqrt(S[m, m])
    elif len(edges) == 0:
        Om = csqrt(np.diag(S)); y[0] = 1.0; branch = 'no edge'
    else: raise AssertionError("graph outside the case list of Lemma 5: %s" % edges)
    return Om, mu, y, branch
def factorise(U):
    C = compound(U); d = np.linalg.det(U)
    a_ = np.array([C[IDX[S], IDX[(0, 1)]] * C[IDX[S], IDX[(2, 3)]] for S, _ in MATCH]); b_ = np.array([C[IDX[S], IDX[(0, 2)]] * C[IDX[S], IDX[(1, 3)]] for S, _ in MATCH])
    r = np.abs(a_); assert max(r[M] - (r.sum() - r[M]) for M in range(3)) < 1e-9, "Lemma 1 fails"
    # explicit phases (remark after Lemma 3): a_M = x_M yv_M with sum x_M conj(yv_M) = 0 by orthogonality of the columns of U
    x = np.array([C[IDX[S], IDX[(0, 1)]] for S, _ in MATCH]); yv = np.array([C[IDX[S], IDX[(2, 3)]] for S, _ in MATCH])
    assert abs((x * np.conj(yv)).sum()) < 1e-12, "column orthogonality fails"
    xi = np.array([np.conj(v) / v if abs(v) >= TOL else 1.0 for v in yv]); kind = triangle_kind(a_); q = (xi * b_).sum()
    if abs(q) >= TOL: xi = xi * 1j * np.conj(q) / abs(q)
    e1 = csqrt(d * d * xi.prod()); e = np.array([e1, xi[0] / e1, xi[1] / e1, xi[2] / e1]); S = U.T @ np.diag(e) @ U           # Lemma 2
    assert max(abs(np.linalg.det(S[np.ix_(T, Tc)])) for T, Tc in MATCH) < 1e-9, "Lemma 3 fails"
    Om, mu, y, branch = lemma5(S); assert np.abs(Om[:, None] * (np.eye(4) + mu * np.outer(y, y)) * Om[None, :] - S).max() < 1e-7 and abs(abs(1 + mu) - 1) < 1e-7
    psi = np.angle(1 + mu) / 2; W = np.eye(4) + (np.exp(1j * psi) - 1) * np.outer(y, y); D1 = csqrt(e); D2 = 1 / Om
    V = D1[:, None] * U * D2[None, :]; O = V @ np.linalg.inv(W)
    return D1, D2, O, y, psi, kind, branch
def sq_minors(Mx):
    return np.array([abs(np.linalg.det(Mx[np.ix_(S, T)]))**2 for k in (1, 2, 3) for S in itertools.combinations(range(4), k) for T in itertools.combinations(range(4), k)])
def haar(n):
    z = (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))) / np.sqrt(2); q, r = np.linalg.qr(z); dd = np.diag(r); return q * (dd / np.abs(dd))
def haar_o(n):
    q, r = np.linalg.qr(rng.standard_normal((n, n))); return q * np.sign(np.diag(r))
def phases(n=4): return np.diag(np.exp(1j * rng.uniform(0, 2 * np.pi, n)))
def perm(): return np.eye(4)[rng.permutation(4)]
def dsum(*blocks):
    M = np.zeros((4, 4), complex); p = 0
    for B in blocks: k = B.shape[0]; M[p:p + k, p:p + k] = B; p += k
    return M
H2 = np.array([[1, 1], [1, -1]]) / np.sqrt(2); F4 = np.array([[1j**(j * k) for k in range(4)] for j in range(4)]) / 2
def zero_corner():
    """Haar-like unitary with u_11 = 0: first column orthogonal to e_1"""
    U = haar(4); v = U[:, 0].copy(); v[0] = 0; v /= np.linalg.norm(v); M = np.column_stack([v, haar(4)[:, 1:]]); q, _ = np.linalg.qr(M); return q
FAMILIES = {
    'Haar unitary': lambda: haar(4),
    'permutation': perm,
    'monomial (permutation times phases)': lambda: phases() @ perm(),
    'diagonal': phases,
    'real orthogonal': lambda: haar_o(4).astype(complex),
    'phases x real orthogonal x phases': lambda: phases() @ haar_o(4) @ phases(),
    'U(2) + U(2), permuted': lambda: perm() @ dsum(haar(2), haar(2)) @ perm(),
    'U(3) + U(1), permuted': lambda: perm() @ dsum(haar(3), haar(1)) @ perm(),
    'U(1) + U(1) + U(2), permuted': lambda: perm() @ dsum(haar(1), haar(1), haar(2)) @ perm(),
    'O(2) + U(2)': lambda: dsum(haar_o(2), haar(2)),
    'Fourier matrix, permuted, with phases': lambda: phases() @ perm() @ F4 @ perm(),
    'real Hadamard, permuted, with phases': lambda: phases() @ perm() @ np.kron(H2, H2) @ perm() @ phases(),
    'complex Hadamard family': lambda: (lambda t: np.array([[1, 1, 1, 1], [1, np.exp(1j * t), -1, -np.exp(1j * t)], [1, -1, 1, -1], [1, -np.exp(1j * t), -1, np.exp(1j * t)]]) / 2)(rng.uniform(0, 2 * np.pi)),
    'tensor product U(2) x U(2)': lambda: np.kron(haar(2), haar(2)),
    'tensor product O(2) x U(2)': lambda: np.kron(haar_o(2), haar(2)),
    'u_11 = 0': zero_corner,
    'one Givens rotation with phases': lambda: phases() @ dsum(haar(2), np.eye(2)) @ phases(),
}
if __name__ == '__main__':
    NT = int(sys.argv[2]) if len(sys.argv) > 2 else 200; worst = 0; total = 0
    print("%-40s %5s  %-9s %-9s %-9s  branches of Lemma 5 / triangles of Lemma 3" % ("family", "n", "|V - OW|", "|Im O|", "F1"))
    for name, gen in FAMILIES.items():
        e1 = e2 = e3 = 0; br = {}; tr = {}
        for _ in range(NT):
            U = gen(); D1, D2, O, y, psi, kind, branch = factorise(U); total += 1; br[branch] = br.get(branch, 0) + 1; tr[kind] = tr.get(kind, 0) + 1
            W = np.eye(4) + (np.exp(1j * psi) - 1) * np.outer(y, y); Or = O.real; c2 = np.cos(psi / 2)**2
            e1 = max(e1, np.abs(D1[:, None] * U * D2[None, :] - Or @ W).max()); e2 = max(e2, np.abs(O.imag).max(), np.abs(Or @ Or.T - np.eye(4)).max())
            e3 = max(e3, np.abs(sq_minors(U) - c2 * sq_minors(Or) - (1 - c2) * sq_minors(Or @ (np.eye(4) - 2 * np.outer(y, y)))).max())
        worst = max(worst, e1, e2, e3); print("%-40s %5d  %.1e   %.1e   %.1e   %s / %s" % (name, NT, e1, e2, e3, br, tr))
    # Lemma 5 directly, on symmetric unitaries S = Om (I + mu y y^T) Om and on the 2x2 block with zero diagonal
    e5 = 0; br = {}
    for t in range(4 * NT):
        if t % 4 == 0: al = rng.uniform(0, 2 * np.pi); S = dsum(np.exp(1j * al) * np.array([[0, 1], [1, 0]]), phases(2)); P = perm(); S = P @ S @ P.T
        else:
            y0 = np.zeros(4); k = (2, 3, 4)[t % 4 - 1]; y0[:k] = rng.standard_normal(k); y0 /= np.linalg.norm(y0); y0 = perm() @ y0; Om0 = np.exp(1j * rng.uniform(0, 2 * np.pi, 4))
            S = Om0[:, None] * (np.eye(4) + (np.exp(1j * rng.uniform(0, 2 * np.pi)) - 1) * np.outer(y0, y0)) * Om0[None, :]
        Om, mu, y, branch = lemma5(S); br[branch] = br.get(branch, 0) + 1
        e5 = max(e5, np.abs(Om[:, None] * (np.eye(4) + mu * np.outer(y, y)) * Om[None, :] - S).max(), abs(abs(1 + mu) - 1))
    worst = max(worst, e5); print("Lemma 5 on %d symmetric unitaries: largest error %.1e ; branches %s" % (4 * NT, e5, br))
    print("%d factorisations; largest error %.1e ; %s" % (total, worst, "ALL CASES COVERED" if worst < 1e-6 else "A CHECK FAILED"))
