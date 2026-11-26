"""Constructive factorisation of a 4x4 unitary, following the proof of Theorem F in README.md step by step.

Input: a generic U in U(4) for which the off-diagonal rank-one reconstruction is numerically resolved.
See Lemma 5 in README.md for the other cases. Output: diagonal unitaries D1, D2, a real orthogonal O, a real unit vector y and psi with
        D1 U D2 = O (I + (e^{i psi} - 1) y y^T).
No optimisation is used: the phases come from a triangle with side lengths |C[S_M,12] C[S_M,34]| (which satisfy the
triangle inequality by the Pluecker relation for columns 1, 2 and Jacobi's complementary-minor identity) and one rotation.
Floating point; the script reports the errors of every intermediate claim of the proof."""
import itertools, sys, numpy as np
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0); NU = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
P2 = list(itertools.combinations(range(4), 2)); IDX = {S: k for k, S in enumerate(P2)}
MATCH = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
def haar():
    z = (rng.standard_normal((4, 4)) + 1j * rng.standard_normal((4, 4))) / np.sqrt(2); q, r = np.linalg.qr(z); d = np.diag(r); return q * (d / np.abs(d))
def compound(U): return np.array([[np.linalg.det(U[np.ix_(S, T)]) for T in P2] for S in P2])
def triangle_phases(c):
    """unit complex numbers xi with sum xi_M c_M = 0, for complex c whose moduli satisfy the triangle inequality"""
    r = np.abs(c); o = np.argsort(-r); a, b, g = r[o]          # a >= b >= g, a <= b + g
    if a < 1e-300: return np.ones(3, complex)
    # place a along the positive real axis; b and g at angles pi - B and pi + G with the law of cosines
    cosB = np.clip((a * a + b * b - g * g) / (2 * a * b), -1, 1) if b > 0 else 1.0
    vb = -b * np.exp(-1j * np.arccos(cosB)); vg = -(a + vb)                     # a + vb + vg = 0, |vg| = g
    target = np.zeros(3, complex); target[o[0]] = a; target[o[1]] = vb; target[o[2]] = vg
    xi = np.ones(3, complex)
    for M in range(3):
        if r[M] > 0: xi[M] = target[M] / c[M]
    return xi / np.where(np.abs(xi) > 0, np.abs(xi), 1)
def factorise(U):
    U = np.asarray(U, dtype=complex)
    if U.shape != (4, 4) or not np.isfinite(U).all() or not np.allclose(U.conj().T @ U, np.eye(4), atol=1e-12, rtol=0):
        raise ValueError("Expected a finite 4x4 unitary matrix")
    info = {}; C = compound(U); d = np.linalg.det(U)
    c1 = np.array([C[IDX[S], IDX[(0, 1)]] * C[IDX[S], IDX[(2, 3)]] for S, _ in MATCH]); c2 = np.array([C[IDX[S], IDX[(0, 2)]] * C[IDX[S], IDX[(1, 3)]] for S, _ in MATCH])
    r = np.abs(c1); info['triangle'] = max(r[M] - (r.sum() - r[M]) for M in range(3))          # <= 0 is the triangle inequality
    xi = triangle_phases(c1); info['closure'] = abs((xi * c1).sum())
    q = (xi * c2).sum(); rot = 1j * np.conj(q) / abs(q) if abs(q) > 1e-300 else 1.0; xi = rot * xi           # Re(sum xi c2) = 0, sum xi c1 = 0
    th = np.angle(xi) + np.angle(d)                                                # (a_i + a_j - a_k - a_l)/2 = theta_M + arg d
    a = np.array([th.sum(), th[0] - th[1] - th[2], -th[0] + th[1] - th[2], -th[0] - th[1] + th[2]]) / 2
    E = np.exp(1j * a); S = U.T @ np.diag(E) @ U
    info['minors'] = max(abs(np.linalg.det(S[np.ix_(T, Tc)])) for T, Tc in MATCH)
    # Below roundoff, the generic chart cannot resolve the nonzero entries used as divisors.
    if min(abs(S[i, j]) for i, j in P2) <= 32 * np.finfo(float).eps:
        raise ValueError("The generic numerical chart is singular or unresolved; see Lemma 5 for the other cases")
    x1 = np.sqrt(S[0, 1] * S[0, 2] / S[1, 2] + 0j); x = np.array([x1, S[0, 1] / x1, S[0, 2] / x1, S[0, 3] / x1])
    Lam = np.diag(S - np.outer(x, x)); info['offdiag'] = np.abs(S - np.outer(x, x) - np.diag(Lam)).max(); info['unimodular'] = np.abs(np.abs(Lam) - 1).max()
    h = Lam * np.conj(x) / x; info['h_equal'] = np.abs(h - h[0]).max(); h = h.mean(); n2 = (np.abs(x)**2).sum(); info['Re_h'] = abs(2 * h.real + n2)
    y = np.abs(x) / np.sqrt(n2); om = np.sqrt(h) * x / np.abs(x); mu = n2 * np.conj(h); psi = np.angle(1 + mu) / 2
    D1 = np.exp(1j * a / 2); D2 = np.conj(om); V = D1[:, None] * U * D2[None, :]
    W = np.eye(4) + (np.exp(1j * psi) - 1) * np.outer(y, y)
    if np.abs(V.T @ V - W @ W).max() > 1e-6: psi += np.pi; W = np.eye(4) + (np.exp(1j * psi) - 1) * np.outer(y, y)
    info['VtV'] = np.abs(V.T @ V - W @ W).max(); O = V @ np.linalg.inv(W); info['imO'] = np.abs(O.imag).max(); info['orth'] = np.abs(O.real @ O.real.T - np.eye(4)).max()
    if not all(np.isfinite(value) for value in info.values()):
        raise ValueError("The generic factorization produced a nonfinite diagnostic")
    return D1, D2, O.real, y, psi, info
def sq_minors(Mx):
    return np.array([abs(np.linalg.det(Mx[np.ix_(S, T)]))**2 for k in (1, 2, 3) for S in itertools.combinations(range(4), k) for T in itertools.combinations(range(4), k)])
if NU < 1:
    raise SystemExit("The number of Haar samples must be positive")
worst = {}; wf = 0
for n in range(NU):
    U = haar(); D1, D2, O, y, psi, info = factorise(U)
    for k, v in info.items(): worst[k] = max(worst.get(k, -np.inf), v)
    c2 = np.cos(psi / 2)**2
    error = np.abs(sq_minors(U) - c2 * sq_minors(O) - (1 - c2) * sq_minors(O @ (np.eye(4) - 2 * np.outer(y, y)))).max()
    if not np.isfinite(error):
        raise ValueError("The squared-minor reconstruction produced a nonfinite error")
    wf = max(wf, error)
print("%d Haar unitaries, construction by the proof (no optimisation). Largest values of:" % NU)
print("   triangle defect max_M (|c_M| - sum of the other two): %.2e (<= 0 expected)" % worst['triangle'])
print("   |sum xi c| after closing the triangle: %.1e ; complementary minors of S: %.1e ; off-diagonal part minus rank one: %.1e" % (worst['closure'], worst['minors'], worst['offdiag']))
print("   ||Lambda_j| - 1|: %.1e ; spread of h_j: %.1e ; |2 Re h + |x|^2|: %.1e ; |V^T V - W^2|: %.1e" % (worst['unimodular'], worst['h_equal'], worst['Re_h'], worst['VtV']))
print("   |Im O|: %.1e ; |O O^T - I|: %.1e ; error of the squared-minor formula over all 68 minors of sizes 1 to 3: %.1e" % (worst['imO'], worst['orth'], wf))
if max(max(worst.values()), wf) > 1e-11:
    raise SystemExit("A numerical residual exceeds the reported tolerance of 1e-11")
