"""Numerical factorisation of 4x4 unitaries (floating point; evidence, not proof).

Claim tested:  every U in U(4) can be written   D1 U D2 = O (I + (e^{i psi} - 1) y y^T)
with D1, D2 diagonal unitary, O real orthogonal, y a real unit vector.  Equivalently: Im(D1 U D2) has rank <= 1.
For such U, Proposition S of README.md gives, for all minors of all sizes,
        |det U[S,T]|^2 = cos^2(psi/2) det O[S,T]^2 + sin^2(psi/2) det (O R_y)[S,T]^2 ,     R_y = I - 2 y y^T.

For each test matrix: solve for the phases by least squares on the 2x2 minors of Im(D1 U D2), recover (O, y, psi), and
report  (i) how far O is from real orthogonal,  (ii) the largest error in the displayed formula over all 1x1, 2x2 and
3x3 minors,  (iii) the margin of (star) at U and at the two orthogonal endpoints."""
import itertools, sys, numpy as np
from scipy.optimize import least_squares
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0); NU = int(sys.argv[2]) if len(sys.argv) > 2 else 200
P2 = list(itertools.combinations(range(4), 2))
def haar():
    z = (rng.standard_normal((4, 4)) + 1j * rng.standard_normal((4, 4))) / np.sqrt(2); q, r = np.linalg.qr(z); d = np.diag(r); return q * (d / np.abs(d))
def dressed(x, U): return np.exp(1j * x[:4])[:, None] * U * np.exp(1j * np.r_[0, x[4:]])[None, :]
def resid(x, U):
    Y = dressed(x, U).imag; return np.array([Y[S[0], T[0]] * Y[S[1], T[1]] - Y[S[0], T[1]] * Y[S[1], T[0]] for S in P2 for T in P2])
def factorise(U, starts=200):
    best = None
    for _ in range(starts):
        r = least_squares(resid, rng.uniform(-np.pi, np.pi, 7), args=(U,), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=400)
        if best is None or np.abs(r.fun).max() < np.abs(best.fun).max(): best = r
        if np.abs(best.fun).max() < 1e-13: break
    V = dressed(best.x, U); u_, s_, vt = np.linalg.svd(V.imag); y = vt[0]
    Vy = V @ y; k = np.argmax(np.abs(Vy)); psi = np.angle(Vy[k]); xr = (Vy * np.exp(-1j * psi)).real
    if xr[k] < 0: psi += np.pi
    O = V @ (np.eye(4) + (np.exp(-1j * psi) - 1) * np.outer(y, y))
    return V, O, y, psi, np.abs(best.fun).max()
def minors(Mx):
    out = []
    for k in (1, 2, 3):
        for S in itertools.combinations(range(4), k):
            for T in itertools.combinations(range(4), k): out.append(abs(np.linalg.det(Mx[np.ix_(S, T)]))**2)
    return np.array(out)
def margin(Mx):
    D = lambda S, T: abs(np.linalg.det(Mx[np.ix_(S, T)]))**2
    return D((0, 2), (0, 2)) + D((0, 3), (0, 3)) + D((0, 2), (1, 2)) + (np.sqrt(D((0, 1), (0, 1))) + np.sqrt(D((0, 2), (0, 3))))**2 - D((0, 1), (2, 3))
tests = [("Haar", haar()) for _ in range(NU)]
S4 = np.array([[0, 1, 1, -1], [-1, 0, -1, -1], [-1, 1, 0, 1], [1, 1, -1, 0]], float)
F = np.exp(2j * np.pi * np.outer(range(4), range(4)) / 4) / 2
tests += [("Fourier matrix", F), ("near a permutation", np.linalg.qr(np.eye(4)[[1, 2, 3, 0]] + 0.05 * (rng.standard_normal((4, 4)) + 1j * rng.standard_normal((4, 4))))[0]),
          ("Drury point times phases", np.diag(np.exp(1j * rng.uniform(0, 6, 4))) @ (0.168231 * np.eye(4) + np.sqrt((1 - 0.168231**2) / 3) * S4) @ np.diag(np.exp(1j * rng.uniform(0, 6, 4))))]
worst = dict(res=0, imO=0, orth=0, formula=0); minmarg = np.inf; lam = []; fails = 0
for name, U in tests:
    V, O, y, psi, res = factorise(U)
    Or = O.real; Op = Or @ (np.eye(4) - 2 * np.outer(y, y)); c2 = np.cos(psi / 2)**2
    err = np.abs(minors(U) - (c2 * minors(Or) + (1 - c2) * minors(Op))).max(); imO = np.abs(O.imag).max(); orth = np.abs(Or @ Or.T - np.eye(4)).max()
    ok = max(orth, err) < 1e-9 and imO < 1e-6; fails += not ok      # |Im O| is ill-conditioned when U is already orthogonal up to phases
    worst['res'] = max(worst['res'], res); worst['imO'] = max(worst['imO'], imO); worst['orth'] = max(worst['orth'], orth); worst['formula'] = max(worst['formula'], err)
    minmarg = min(minmarg, margin(U), margin(Or), margin(Op)); lam.append(c2)
    if name != "Haar": print("   %-26s: psi = %+.4f, cos^2(psi/2) = %.4f, |Im O| = %.1e, formula error = %.1e, margins of (star): U %.4f, O %.4f, O R_y %.4f" % (name, psi, c2, imO, err, margin(U), margin(Or), margin(Op)))
print("%d unitaries (%d Haar + 3 special): factorisation failed (formula or orthogonality error above 1e-9, or |Im O| above 1e-6) for %d" % (len(tests), NU, fails))
print("   largest: rank-one residual %.1e ; |Im O| %.1e ; |O O^T - I| %.1e ; error in the squared-minor formula over all 164 minors %.1e" % (worst['res'], worst['imO'], worst['orth'], worst['formula']))
print("   cos^2(psi/2) ranged over [%.4f, %.4f] ; smallest margin of (star) seen at U or at an endpoint: %.3e" % (min(lam), max(lam), minmarg))
