#!/usr/bin/env python3
"""Seeded floating-point spot checks of identities used in ../reductions/README.md
and in Kovacec p2611 v2 section 2. These are regression checks for sign and index
errors, NOT proofs: every identity here is also argued symbolically in adversarial/README.md
or in the cited source. Bounds: n in 2..6, 200 trials per identity, seed 20261125.

Usage: python adversarial/identity_spotchecks.py
"""
import numpy as np

SEED, TRIALS, TOL = 20261125, 200, 1e-9
rng = np.random.default_rng(SEED)


def crandn(*shape):
    return rng.standard_normal(shape) + 1j * rng.standard_normal(shape)


def haar(n):
    q, r = np.linalg.qr(crandn(n, n))
    return q * (np.diag(r) / abs(np.diag(r)))


def d(a, b, U):
    return np.linalg.det(np.diag(a) + U @ np.diag(b) @ U.conj().T)


def rel(x, y):
    return abs(x - y) / max(1.0, abs(x), abs(y))


worst = {}
def record(name, err):
    worst[name] = max(worst.get(name, 0.0), err)


for _ in range(TRIALS):
    n = int(rng.integers(2, 7))
    a, b, U = crandn(n), crandn(n), haar(n)
    base = d(a, b, U)

    # LUR phase invariance (reductions §3, last table row)
    L = np.diag(np.exp(1j * rng.uniform(0, 2 * np.pi, n)))
    R = np.diag(np.exp(1j * rng.uniform(0, 2 * np.pi, n)))
    record("LUR phase invariance", rel(d(a, b, L @ U @ R), base))

    # swap a,b and U -> U*  (reductions §3)
    record("swap spectra, U->U*", rel(d(b, a, U.conj().T), base))

    # opposite shifts (reductions §3)
    t = crandn(1)[0]
    record("opposite shift", rel(d(a + t, b - t, U), base))

    # complex scaling multiplies by c^n (used for theta=0 sufficiency)
    c = crandn(1)[0]
    record("scaling by c^n", rel(d(c * a, c * b, U), c ** n * base))

    # reciprocal transform (reductions §3; Shitov Thm 12)
    cc = crandn(1)[0]
    A = np.diag(a); B = U @ np.diag(b) @ U.conj().T
    X, Y = A - cc * np.eye(n), B + cc * np.eye(n)
    lhs = np.linalg.det(np.linalg.inv(X) + np.linalg.inv(Y))
    rhs = base / (np.linalg.det(X) * np.linalg.det(Y))
    record("reciprocal transform", rel(lhs, rhs))
    sigma = rng.permutation(n)
    zp = np.prod(1 / (a - cc) + 1 / (b[sigma] + cc))
    z = np.prod(a + b[sigma]) / (np.prod(a - cc) * np.prod(b + cc))
    record("reciprocal vertex scaling", rel(zp, z))

    # Cayley transform (reductions §4, Lemmas 1-2)
    K = crandn(n, n); K = (K - K.conj().T) / 2
    I = np.eye(n)
    C = (I - K) @ np.linalg.inv(I + K)
    record("Cayley unitary", np.abs(C.conj().T @ C - I).max())
    Kinv = (I - U) @ np.linalg.inv(I + U)
    record("Cayley inverse skew", np.abs(Kinv + Kinv.conj().T).max())
    record("Cayley round trip", np.abs((I - Kinv) @ np.linalg.inv(I + Kinv) - U).max())

    # Kovacec p2611 v2 §2 Corollary: generalized Givens image formula
    i, j = sorted(rng.choice(n, 2, replace=False))
    v = rng.standard_normal(4); cq, bq, sq, dq = v / np.linalg.norm(v)
    G = np.eye(n, dtype=complex)
    G[i, i], G[i, j] = cq + 1j * bq, sq + 1j * dq
    G[j, i], G[j, j] = -sq + 1j * dq, cq - 1j * bq
    record("gGivens unitary", np.abs(G.conj().T @ G - I).max())
    M = np.prod([a[k] + b[k] for k in range(n) if k not in (i, j)])
    pred = M * ((a[i] + b[i]) * (a[j] + b[j]) + (a[i] - a[j]) * (b[i] - b[j]) * (dq ** 2 + sq ** 2))
    record("gGivens determinant formula", rel(d(a, b, G), pred))

ok = all(v < TOL for v in worst.values())
for k, v in worst.items():
    print(f"{k:32s} worst error {v:.2e}")
print(f"seed={SEED} trials={TRIALS} tol={TOL}: {'all within tolerance' if ok else 'FAILURES PRESENT'}")
raise SystemExit(0 if ok else 1)
