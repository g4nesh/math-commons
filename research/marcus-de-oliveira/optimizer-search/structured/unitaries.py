"""Structured unitaries for the fixed-U search."""
import numpy as np

O_STAR = np.array([[1, 4, 4, 4], [-4, 1, 4, -4], [-4, -4, 1, 4], [-4, 4, -4, 1]]) / 7.0


def fourier(n):
    j = np.arange(n)
    return np.exp(2j * np.pi * np.outer(j, j) / n) / np.sqrt(n)


def circulant(n, rng):
    F = fourier(n)
    return F @ np.diag(np.exp(2j * np.pi * rng.random(n))) @ F.conj().T


def twisted_hadamard(p, q, rng):
    """(F_p (x) I_q) D (I_p (x) F_q): a complex Hadamard matrix for any phases D."""
    D = np.diag(np.exp(2j * np.pi * rng.random(p * q)))
    return np.kron(fourier(p), np.eye(q)) @ D @ np.kron(np.eye(p), fourier(q))


def sylvester8():
    H = np.array([[1, 1], [1, -1]])
    return np.kron(np.kron(H, H), H) / np.sqrt(8)


def givens(n, i, j, t, phi=0.0):
    G = np.eye(n, dtype=complex)
    c, s = np.cos(t), np.sin(t)
    G[i, i] = c
    G[j, j] = c
    G[i, j] = -s * np.exp(1j * phi)
    G[j, i] = s * np.exp(-1j * phi)
    return G


def block(*Bs):
    n = sum(len(B) for B in Bs)
    W = np.zeros((n, n), complex)
    k = 0
    for B in Bs:
        m = len(B)
        W[k:k + m, k:k + m] = B
        k += m
    return W


def ostar_sum(n, rng, mix=0.5, complex_phase=False, second=None):
    """O* (+) V rotated out of block form by plane rotations that couple the two
    blocks on both sides. V is I, F_{n-4}, or O* (n = 8)."""
    k = n - 4
    if second == "ostar" and k == 4:
        V = O_STAR
    elif second == "fourier" and k >= 2:
        V = fourier(k)
    else:
        V = np.eye(k)
    W = block(O_STAR, V).astype(complex)
    L = np.eye(n, dtype=complex)
    R = np.eye(n, dtype=complex)
    for _ in range(2):
        i, j = rng.integers(0, 4), 4 + rng.integers(0, k)
        ph = 2 * np.pi * rng.random() if complex_phase else 0.0
        L = givens(n, i, j, mix * rng.standard_normal(), ph) @ L
        i, j = rng.integers(0, 4), 4 + rng.integers(0, k)
        ph = 2 * np.pi * rng.random() if complex_phase else 0.0
        R = R @ givens(n, i, j, mix * rng.standard_normal(), ph)
    return L @ W @ R


def ostar_kron(rng, t=None, complex_phase=False):
    """n = 8: O* (x) R(t), rotated out of the tensor form by a random plane rotation pair."""
    t = np.pi / 4 if t is None else t
    R2 = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    W = np.kron(O_STAR, R2).astype(complex)
    if complex_phase:
        W = W @ np.diag(np.exp(2j * np.pi * rng.random(8)))
        W = W @ np.kron(np.eye(4), fourier(2)) @ W.T.conj() @ W
    return W


def cycle_rotations(n, rng, length=None, reps=1, complex_phase=True, big=True):
    """Product of plane rotations on consecutive pairs of a random cycle
    (c0 c1)(c1 c2)...(c_{L-1} c0), repeated reps times, random angles."""
    L = n if length is None else length
    cyc = rng.permutation(n)[:L]
    U = np.eye(n, dtype=complex)
    for _ in range(reps):
        for k in range(L):
            i, j = cyc[k], cyc[(k + 1) % L]
            t = rng.uniform(0.2, np.pi / 2 - 0.2) if big else 0.4 * rng.standard_normal()
            ph = 2 * np.pi * rng.random() if complex_phase else 0.0
            U = givens(n, i, j, t, ph) @ U
    return U


def sq_minors(U, ks=None):
    """Vector of |det U[S,T]|^2 for |S| = |T| = k, k in ks."""
    import itertools
    n = len(U)
    ks = range(1, n) if ks is None else ks
    out = []
    for k in ks:
        subs = list(itertools.combinations(range(n), k))
        for S in subs:
            for T in subs:
                out.append(abs(np.linalg.det(U[np.ix_(S, T)])) ** 2)
    return np.array(out)
