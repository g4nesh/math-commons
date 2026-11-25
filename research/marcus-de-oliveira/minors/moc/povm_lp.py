"""POVM generalisation: for E_0..E_{m-1} >= 0 with sum I (n x n), is there a random colouring
kappa:[n]->[m] with, for every S, the law of the colour-count vector of kappa on S equal to the
coefficients of det(sum_j c_j E_j[S,S]) (a probability vector: mixed discriminants)?
Unitary/pattern case = E_j orthogonal projections."""
import sys, itertools, numpy as np
from scipy.optimize import linprog

def rand_povm(n, m, rng, ranks=None, proj=False):
    # Naimark: isometry V: C^n -> C^D, blocks of C^D
    if ranks is None: ranks = [n] * m
    ranks = list(ranks)
    while sum(ranks) < n: ranks[rng.integers(len(ranks))] += 1
    D = sum(ranks)
    Z = rng.standard_normal((D, n)) + 1j * rng.standard_normal((D, n))
    Q, _ = np.linalg.qr(Z)            # D x n isometry
    E = []; s = 0
    for r in ranks:
        B = Q[s:s+r, :]; s += r
        E.append(B.conj().T @ B)
    return E

def count_stats(E, n, m):
    """dict (S, counts) -> coefficient"""
    out = {}
    for k in range(1, n + 1):
        for S in itertools.combinations(range(n), k):
            S = list(S)
            for tau in itertools.product(range(m), repeat=k):
                Mx = np.array([E[tau[r]][S[r], S] for r in range(k)])
                val = np.linalg.det(Mx)
                cnt = tuple(np.bincount(tau, minlength=m))
                out[(tuple(S), cnt)] = out.get((tuple(S), cnt), 0) + val
    return out

def colouring_matrix(n, m, keys):
    cols = list(itertools.product(range(m), repeat=n))
    kidx = {k: i for i, k in enumerate(keys)}
    V = np.zeros((len(cols), len(keys)))
    for c, kap in enumerate(cols):
        for k in range(1, n + 1):
            for S in itertools.combinations(range(n), k):
                cnt = tuple(np.bincount([kap[i] for i in S], minlength=m))
                V[c, kidx[(S, cnt)]] = 1
    return cols, V

def gap(E, n, m, V, keys):
    st = count_stats(E, n, m)
    p = np.array([st[k_] for k_ in keys])
    assert np.abs(p.imag).max() < 1e-9
    p = p.real
    M, d = V.shape
    c = np.zeros(M + 1); c[-1] = 1
    A_ub = np.vstack([np.hstack([V.T, -np.ones((d, 1))]), np.hstack([-V.T, -np.ones((d, 1))])])
    r = linprog(c, A_ub=A_ub, b_ub=np.concatenate([p, -p]), A_eq=np.hstack([np.ones((1, M)), [[0]]]),
                b_eq=[1.0], bounds=[(0, None)] * (M + 1), method='highs')
    return r.fun

if __name__ == '__main__':
    rng = np.random.default_rng(0)
    for (n, m, trials) in [(2, 3, 300), (3, 3, 300), (4, 3, 150), (5, 3, 40), (3, 4, 200), (4, 4, 60)]:
        st = count_stats(rand_povm(n, m, rng), n, m); keys = sorted(st.keys())
        cols, V = colouring_matrix(n, m, keys)
        g = [gap(rand_povm(n, m, rng, ranks=[rng.integers(1, n + 1) for _ in range(m)]), n, m, V, keys) for _ in range(trials)]
        g = np.array(g)
        print("n=%d colours=%d trials=%d: infeasible=%d, worst gap=%.2e" % (n, m, trials, (g > 1e-9).sum(), g.max()), flush=True)
