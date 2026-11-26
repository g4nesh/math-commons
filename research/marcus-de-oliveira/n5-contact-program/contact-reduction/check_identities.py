"""Exact and numerical checks of the identities used in the contact reduction.

Run from the repository root:
    python -I -B research/marcus-de-oliveira/n5-contact-program/contact-reduction/check_identities.py

Every claim prints one line "PASS"/"FAIL" with a label EXACT (exact rational
arithmetic, python-flint fmpq/fmpz or sympy) or NUMERIC (floating point).

(a) Lemma 1: for U unitary, |det U[K,J]|^2 = sum_{sigma(K)=J} lambda_sigma(U),
    M(mu(U)) = (|U_ij|^2), sum mu = 1, |mu_sigma| <= 1 and
    det(D_a + U D_b U^*) = sum_sigma mu_sigma(U) z_sigma(a,b).
(b) rank of the Jacobian of Phi at a random integer configuration = 2n-2.
(c) slice formulas (sympy).
(d) kappa-line identity and rank of A = [vec P_sigma].
(e) NUMERIC Kempf-Ness descent.
"""

import itertools
import random
import sys
import time

import numpy as np
import sympy as sp
from flint import fmpq, fmpq_mat, fmpz_mat

RESULTS = []


def report(ok, label, text):
    RESULTS.append(ok)
    print(("PASS" if ok else "FAIL") + f" [{label}] {text}", flush=True)


# ---------------------------------------------------------------- Gaussian rationals
class G:
    """Gaussian rational re + i*im with fmpq parts."""

    __slots__ = ("re", "im")

    def __init__(self, re=0, im=0):
        self.re = re if isinstance(re, fmpq) else fmpq(re)
        self.im = im if isinstance(im, fmpq) else fmpq(im)

    def __add__(self, o):
        return G(self.re + o.re, self.im + o.im)

    def __sub__(self, o):
        return G(self.re - o.re, self.im - o.im)

    def __neg__(self):
        return G(-self.re, -self.im)

    def __mul__(self, o):
        return G(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)

    def conj(self):
        return G(self.re, -self.im)

    def abs2(self):
        return self.re * self.re + self.im * self.im

    def inv(self):
        d = self.abs2()
        return G(self.re / d, -self.im / d)

    def iszero(self):
        return self.re == 0 and self.im == 0

    def __eq__(self, o):
        return self.re == o.re and self.im == o.im

    def scale(self, q):
        return G(self.re * q, self.im * q)


ZERO, ONE = G(0, 0), G(1, 0)


def gdet(M):
    """Determinant of a square list-of-lists of G by fraction-full elimination."""
    n = len(M)
    if n == 0:
        return ONE
    A = [row[:] for row in M]
    det = ONE
    for c in range(n):
        p = next((r for r in range(c, n) if not A[r][c].iszero()), None)
        if p is None:
            return ZERO
        if p != c:
            A[c], A[p] = A[p], A[c]
            det = -det
        det = det * A[c][c]
        ic = A[c][c].inv()
        for r in range(c + 1, n):
            if A[r][c].iszero():
                continue
            f = A[r][c] * ic
            A[r] = [A[r][k] - f * A[c][k] for k in range(n)]
    return det


def sign(p):
    s = 1
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


def rand_q(rng, num=3, den=3):
    return fmpq(rng.randint(-num, num), rng.randint(1, den))


def cayley_unitary(n, rng):
    """U = (I-K)(I+K)^{-1}, K skew-Hermitian with Gaussian-rational entries.
    Computed in the real 2n x 2n block representation Z=X+iY -> [[X,-Y],[Y,X]]."""
    X = [[fmpq(0)] * n for _ in range(n)]
    Y = [[fmpq(0)] * n for _ in range(n)]
    for i in range(n):
        Y[i][i] = rand_q(rng)
        for j in range(i + 1, n):
            x, y = rand_q(rng), rand_q(rng)
            X[i][j], X[j][i] = x, -x  # real part antisymmetric
            Y[i][j], Y[j][i] = y, y  # imaginary part symmetric
    R = fmpq_mat(2 * n, 2 * n)
    for i in range(n):
        for j in range(n):
            R[i, j] = X[i][j]
            R[i, j + n] = -Y[i][j]
            R[i + n, j] = Y[i][j]
            R[i + n, j + n] = X[i][j]
    I2 = fmpq_mat(2 * n, 2 * n)
    for i in range(2 * n):
        I2[i, i] = 1
    Ub = (I2 - R) * (I2 + R).inv()
    return [[G(Ub[i, j], Ub[i + n, j]) for j in range(n)] for i in range(n)]


def check_unitary(U):
    n = len(U)
    for i in range(n):
        for j in range(n):
            s = ZERO
            for k in range(n):
                s = s + U[i][k] * U[j][k].conj()
            if not s == (ONE if i == j else ZERO):
                return False
    return True


def lemma1(n, seed):
    rng = random.Random(seed)
    U = cayley_unitary(n, rng)
    ok_unit = check_unitary(U)
    perms = list(itertools.permutations(range(n)))
    cdet = gdet(U).conj()
    lam = []
    for p in perms:
        t = G(sign(p), 0)
        for i in range(n):
            t = t * U[i][p[i]]
        lam.append(t * cdet)
    mu = [l.re for l in lam]
    # accumulate sum_{sigma(K)=J} lambda_sigma over bitmasks
    acc = {}
    for p, l in zip(perms, lam):
        for K in range(1 << n):
            J = 0
            for i in range(n):
                if K >> i & 1:
                    J |= 1 << p[i]
            key = (K, J)
            acc[key] = acc[key] + l if key in acc else l
    npairs, ok_minor = 0, True
    for K in range(1 << n):
        rows = [i for i in range(n) if K >> i & 1]
        for J in range(1 << n):
            cols = [j for j in range(n) if J >> j & 1]
            if len(cols) != len(rows):
                continue
            npairs += 1
            d = gdet([[U[i][j] for j in cols] for i in rows])
            lhs = G(d.abs2(), 0)
            rhs = acc.get((K, J), ZERO)
            if not lhs == rhs:
                ok_minor = False
    ok_M = all(
        sum((mu[k] for k, p in enumerate(perms) if p[i] == j), fmpq(0)) == U[i][j].abs2()
        for i in range(n) for j in range(n))
    ok_sum = sum(mu, fmpq(0)) == 1
    ok_abs = all(abs(m) <= 1 for m in mu)
    # determinant identity
    a = [G(rand_q(rng, 5, 4), rand_q(rng, 5, 4)) for _ in range(n)]
    b = [G(rand_q(rng, 5, 4), rand_q(rng, 5, 4)) for _ in range(n)]
    Mx = []
    for i in range(n):
        row = []
        for j in range(n):
            s = a[i] if i == j else ZERO
            for k in range(n):
                s = s + U[i][k] * b[k] * U[j][k].conj()
            row.append(s)
        Mx.append(row)
    lhs = gdet(Mx)
    rhs_mu, rhs_lam = ZERO, ZERO
    for p, l, m in zip(perms, lam, mu):
        z = ONE
        for i in range(n):
            z = z * (a[i] + b[p[i]])
        rhs_mu = rhs_mu + z.scale(m)
        rhs_lam = rhs_lam + z * l
    ok_det = lhs == rhs_mu and lhs == rhs_lam
    return ok_unit, npairs, ok_minor, ok_M, ok_sum, ok_abs, ok_det, sum(1 for m in mu if m < 0)


def part_a():
    print("== (a) Lemma 1, rational Cayley unitaries ==")
    for n in (3, 4, 5, 6):
        for seed in (1, 2, 3):
            t0 = time.time()
            ok_u, npairs, okm, okM, oks, oka, okd, nneg = lemma1(n, 1000 * n + seed)
            dt = time.time() - t0
            report(ok_u, "EXACT", f"n={n} seed={seed}: U U^* = I")
            report(okm, "EXACT", f"n={n} seed={seed}: |det U[K,J]|^2 = sum_(sigma(K)=J) lambda_sigma for all {npairs} pairs (K,J) incl. empty/full (imag parts cancel)")
            report(okM, "EXACT", f"n={n} seed={seed}: M(mu(U)) = (|U_ij|^2)")
            report(oks, "EXACT", f"n={n} seed={seed}: sum mu = 1")
            report(oka, "EXACT", f"n={n} seed={seed}: |mu_sigma| <= 1 ({nneg} of {len(list(itertools.permutations(range(n))))} mu_sigma negative)")
            report(okd, "EXACT", f"n={n} seed={seed}: det(D_a + U D_b U^*) = sum mu_sigma z_sigma = sum lambda_sigma z_sigma (Gaussian-rational a,b) [{dt:.1f}s]")


# ---------------------------------------------------------------- (b) Jacobian rank
def part_b():
    print("== (b) rank of the Jacobian of Phi at a random integer point ==")
    rng = random.Random(7)
    for n in range(3, 8):
        V = [[rng.randint(-20, 20), rng.randint(-20, 20)] for _ in range(n)]
        W = [[rng.randint(-20, 20), rng.randint(-20, 20)] for _ in range(n)]
        perms = list(itertools.permutations(range(n)))
        rows = []
        for p in perms:
            br = [V[i][0] * W[p[i]][1] - V[i][1] * W[p[i]][0] for i in range(n)]
            row = [0] * (4 * n)
            for i in range(n):
                rest = 1
                for k in range(n):
                    if k != i:
                        rest *= br[k]
                j = p[i]
                # d/dv_i1, d/dv_i2, d/dw_j1, d/dw_j2 of v_i1 w_j2 - v_i2 w_j1
                row[2 * i] += W[j][1] * rest
                row[2 * i + 1] += -W[j][0] * rest
                row[2 * n + 2 * j] += -V[i][1] * rest
                row[2 * n + 2 * j + 1] += V[i][0] * rest
            rows.append(row)
        r = fmpz_mat(rows).rank()
        report(r == 2 * n - 2, "EXACT", f"n={n}: rank of d Phi ({len(perms)} x {4*n}) at a random integer configuration = {r} (expected 2n-2 = {2*n-2})")


# ---------------------------------------------------------------- (c) slice formulas
def part_c():
    print("== (c) slice formulas (sympy) ==")
    for n in (4, 5):
        s = sp.Symbol("s")
        p = {i: sp.Symbol(f"p{i}") for i in range(3, n + 1)}
        q = {j: sp.Symbol(f"q{j}") for j in range(2, n + 1)}
        v = {1: (0, s), 2: (1, 0)}
        for i in range(3, n + 1):
            v[i] = (p[i], 1)
        w = {1: (1, 1)}
        for j in range(2, n + 1):
            w[j] = (q[j], 1)

        def br(u, up):
            return u[0] * up[1] - u[1] * up[0]

        def F(sig):  # sig: dict 1..n -> 1..n
            return sp.Mul(*[br(v[i], w[sig[i]]) for i in range(1, n + 1)])

        perms = [dict(zip(range(1, n + 1), (x + 1 for x in pp)))
                 for pp in itertools.permutations(range(n))]
        ok1, c1 = True, 0
        for j in range(2, n + 1):
            for sig in perms:
                if sig[1] == 1 and sig[2] == j:
                    tau = dict(sig)
                    tau[1], tau[2] = sig[2], sig[1]  # tau = sigma o (1 2)
                    c1 += 1
                    if sp.simplify(F(sig) / F(tau) - 1 / q[j]) != 0:
                        ok1 = False
        report(ok1, "EXACT", f"n={n}: F_sigma/F_tau = 1/q_j for tau = sigma o (1 2), sigma(1)=1, sigma(2)=j ({c1} cases)")
        ok2, c2 = True, 0
        for i in range(3, n + 1):
            for sig in perms:
                if sig[1] == 1 and sig[i] == 2:
                    tau = dict(sig)
                    tau[1], tau[i] = sig[i], sig[1]  # tau = sigma o (1 i)
                    c2 += 1
                    target = (p[i] - q[2]) / (q[2] * (p[i] - 1))
                    if sp.simplify(F(sig) / F(tau) - target) != 0:
                        ok2 = False
        report(ok2, "EXACT", f"n={n}: F_sigma/F_tau = (p_i - q_2)/(q_2 (p_i - 1)) for tau = sigma o (1 i), sigma(1)=1, sigma(i)=2 ({c2} cases)")
        ident = {i: i for i in range(1, n + 1)}
        target = -s * sp.Mul(*[p[i] - q[i] for i in range(3, n + 1)])
        report(sp.expand(F(ident) - target) == 0, "EXACT", f"n={n}: F_id = -s prod_(i>=3) (p_i - q_i)")


# ---------------------------------------------------------------- (d) kappa line, rank A
def part_d():
    print("== (d) kappa line and rank of A ==")
    rng = random.Random(11)
    for n in (4, 5):
        for trial in range(3):
            def rg():
                return G(rand_q(rng, 9, 5), rand_q(rng, 9, 5))
            c, an = rg(), rg()
            b = [rg() for _ in range(n)]
            t = rg()
            Q = ONE
            for j in range(n):
                Q = Q * (c + b[j])
            assert not (Q + t).iszero()
            a = [c] * (n - 1) + [an]
            anp = c + Q * (an - c) * (Q + t).inv()
            ap = [c] * (n - 1) + [anp]
            alpha = (Q + t) * Q.inv()
            ok = True
            for pp in itertools.permutations(range(n)):
                z, zp = ONE, ONE
                for i in range(n):
                    z = z * (a[i] + b[pp[i]])
                    zp = zp * (ap[i] + b[pp[i]])
                if not alpha * zp == z + t:
                    ok = False
            report(ok, "EXACT", f"n={n} trial={trial}: alpha' z(a',b) = z(a,b) + t(1,...,1) for all sigma (Gaussian-rational c, a_n, b, t)")
    for n in (3, 4, 5):
        perms = list(itertools.permutations(range(n)))
        A = [[1 if pp[i] == j else 0 for pp in perms] for i in range(n) for j in range(n)]
        r = fmpz_mat(A).rank()
        report(r == (n - 1) ** 2 + 1, "EXACT", f"n={n}: rank [vec P_sigma] ({n*n} x {len(perms)}) = {r} (expected (n-1)^2+1 = {(n-1)**2+1})")


# ---------------------------------------------------------------- (e) Kempf-Ness
def phi(V, n):
    """V: 2n x 2 complex array (rows v_1..v_n, w_1..w_n)."""
    out = []
    for pp in itertools.permutations(range(n)):
        t = 1.0 + 0j
        for i in range(n):
            u, up = V[i], V[n + pp[i]]
            t *= u[0] * up[1] - u[1] * up[0]
        out.append(t)
    return np.array(out)


def expm_herm(H):
    w, Q = np.linalg.eigh(H)
    return (Q * np.exp(w)) @ Q.conj().T


def kempf_ness(V, iters=200000, eta=0.5, tol=1e-13):
    """Gradient descent of ||g.v||^2 on G = SL2(C) x {prod t = 1}.
    Each step applies an exact group element (det = 1, prod t = 1)."""
    V = V.copy()
    m = V.shape[0]
    hist = []
    for it in range(iters):
        N = np.sum(np.abs(V) ** 2)
        C = V.T @ V.conj()  # sum_k u_k u_k^*
        H = C - np.trace(C).real / 2 * np.eye(2)
        r = np.sum(np.abs(V) ** 2, axis=1) - N / m
        gnorm = np.sqrt(np.sum(np.abs(H) ** 2) + np.sum(r ** 2)) / N
        if it % 1000 == 0 or gnorm < tol:
            hist.append((it, N, gnorm))
        if gnorm < tol:
            break
        g = expm_herm(-eta * H / N)
        V = (V @ g.T) * np.exp(-eta * r / N)[:, None]
    return V, hist, it


def part_e():
    print("== (e) NUMERIC Kempf-Ness descent ==")
    rng = np.random.default_rng(5)
    for n in (4, 5):
        for trial in range(2):
            V = rng.normal(size=(2 * n, 2)) + 1j * rng.normal(size=(2 * n, 2))
            V *= np.exp(rng.normal(size=(2 * n, 1)))  # unbalanced norms
            P0 = phi(V, n)
            U, hist, it = kempf_ness(V)
            P1 = phi(U, n)
            norms = np.sum(np.abs(U) ** 2, axis=1)
            C = U.T @ U.conj()
            rel_phi = np.max(np.abs(P1 - P0)) / np.max(np.abs(P0))
            spread = (norms.max() - norms.min()) / norms.mean()
            Cdev = np.max(np.abs(C - np.trace(C).real / 2 * np.eye(2))) / np.trace(C).real
            ok = rel_phi < 1e-10 and spread < 1e-10 and Cdev < 1e-10
            report(ok, "NUMERIC", f"n={n} trial={trial}: converged in {it} steps, ||v||^2 {np.sum(np.abs(V)**2):.4g} -> {hist[-1][1]:.6g}; |u_k|^2 spread {spread:.1e}; sum u u^* - (tr/2)I rel {Cdev:.1e}; Phi(u) vs Phi(v) rel {rel_phi:.1e}")
    # strictly semistable: v1, v2, w1, w2 on one line (n = 4)
    n = 4
    V = rng.normal(size=(2 * n, 2)) + 1j * rng.normal(size=(2 * n, 2))
    line = np.array([1.0, 0.3 - 0.7j])
    for k in (0, 1, 4, 5):  # v1, v2, w1, w2
        V[k] = (rng.normal() + 1j * rng.normal()) * line
    P0 = phi(V, n)
    nz = int(np.sum(np.abs(P0) > 1e-12 * np.max(np.abs(P0))))
    U, hist, it = kempf_ness(V, iters=400000)
    P1 = phi(U, n)
    rel_phi = np.max(np.abs(P1 - P0)) / np.max(np.abs(P0))
    norms = np.sum(np.abs(U) ** 2, axis=1)
    # distance of the other 4 vectors from being parallel (polystable limit)
    others = U[[2, 3, 6, 7]]
    sv = np.linalg.svd(others, compute_uv=False)
    print(f"  semistable n=4 (v1,v2,w1,w2 parallel): Phi(v) has {nz} nonzero entries of 24, max |Phi| = {np.max(np.abs(P0)):.3g}")
    for (k, N, gn) in hist[:: max(1, len(hist) // 8)] + [hist[-1]]:
        print(f"    step {k:7d}: ||g.v||^2 = {N:.10f}  gradient norm = {gn:.3e}")
    print(f"    after {it} steps: |u_k|^2 = {np.array2string(norms, precision=5)}; singular values of (v3,v4,w3,w4) = {np.array2string(sv, precision=3)}")
    rate = np.polyfit(np.log([h[0] for h in hist[2:]]), np.log([h[2] for h in hist[2:]]), 1)[0]
    report(rel_phi < 1e-10 and nz > 0 and hist[-1][2] > 1e-13 and rate > -1.5,
           "NUMERIC",
           f"n=4 strictly semistable: Phi(v) != 0 preserved to rel {rel_phi:.1e}; descent does NOT reach the tolerance (gradient decays only polynomially, fitted exponent {rate:.2f}); the remaining 4 vectors collapse toward one line (smallest/largest singular value {sv[1]/sv[0]:.2e}), i.e. the infimum is approached only in the limit (polystable limit: two lines of 4)")


if __name__ == "__main__":
    t0 = time.time()
    part_a()
    part_b()
    part_c()
    part_d()
    part_e()
    print(f"== summary: {sum(RESULTS)} PASS, {len(RESULTS) - sum(RESULTS)} FAIL, {time.time() - t0:.1f}s ==")
    sys.exit(0 if all(RESULTS) else 1)
