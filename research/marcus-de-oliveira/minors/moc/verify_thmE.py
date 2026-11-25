"""Direct check of the construction in Theorem E (three distinct eigenvalues, one simple).
B = U diag(b) U*,  b = beta on block J0, gamma on block J (|J|=k), delta at one index m.
Coupling (X subset Y): X ~ DPP(range U[:,J]),  Y ~ DPP(range U[:, J+{m}]),  found by an LP
(exists by Lyons Thm 6.2 + Strassen).  Then  det(A+B) = sum P(X, Y) z(X, Y\\X)  must hold exactly."""
import sys, itertools, numpy as np
from scipy.optimize import linprog
sys.path.insert(0, 'moc')
from moclib import haar_unitary, det_sum

def check(n, k, rng):
    U = haar_unitary(n, rng)
    J = list(range(k)); m = k; J2 = J + [m]
    Xs = list(itertools.combinations(range(n), k)); Ys = list(itertools.combinations(range(n), k + 1))
    pX = np.array([abs(np.linalg.det(U[np.ix_(X, J)]))**2 for X in Xs])
    pY = np.array([abs(np.linalg.det(U[np.ix_(Y, J2)]))**2 for Y in Ys])
    pairs = [(i, j) for i, X in enumerate(Xs) for j, Y in enumerate(Ys) if set(X) <= set(Y)]
    A = np.zeros((len(Xs) + len(Ys), len(pairs)))
    for c, (i, j) in enumerate(pairs): A[i, c] = 1; A[len(Xs) + j, c] = 1
    r = linprog(np.zeros(len(pairs)), A_eq=A, b_eq=np.concatenate([pX, pY]), bounds=[(0, None)] * len(pairs), method='highs')
    assert r.status == 0, "coupling LP infeasible"
    lam = r.x
    a = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    beta, gamma, delta = rng.standard_normal(3) + 1j * rng.standard_normal(3)
    b = np.full(n, beta); b[J] = gamma; b[m] = delta
    d = det_sum(a, b, U)
    conv = 0
    for c, (i, j) in enumerate(pairs):
        X = set(Xs[i]); jj = (set(Ys[j]) - X).pop()
        z = np.prod([a[t] + (gamma if t in X else delta if t == jj else beta) for t in range(n)])
        conv += lam[c] * z
    return abs(d - conv) / abs(d), lam.min(), abs(pX.sum() - 1) + abs(pY.sum() - 1)

rng = np.random.default_rng(42)
for (n, k) in [(4, 1), (4, 2), (5, 2), (5, 3), (6, 2), (6, 3), (6, 4), (7, 3)]:
    errs = [check(n, k, rng) for _ in range(20)]
    e = np.array(errs)
    print("n=%d, multiplicities (beta,gamma,delta)=(%d,%d,1): max rel |det - sum lam z| = %.1e over 20 trials, min weight %.1e" % (n, n - k - 1, k, e[:, 0].max(), e[:, 1].min()))
