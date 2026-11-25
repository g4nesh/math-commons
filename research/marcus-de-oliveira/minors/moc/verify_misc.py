import sys, itertools, numpy as np
sys.path.insert(0, 'moc')
from moclib import haar_unitary, det_sum
from pattern_lp import perm_stats_matrix, stats, make_blocks
from scipy.optimize import linprog
rng = np.random.default_rng(11)
# (1) explicit two-value weights: det(D_a + U D_b U*) = sum_X |det U[X,J]|^2 z_X
worst = 0
for n in range(2, 8):
    for k in range(1, n):
        U = haar_unitary(n, rng); a = rng.standard_normal(n) + 1j * rng.standard_normal(n)
        be, ga = rng.standard_normal(2) + 1j * rng.standard_normal(2)
        J = list(range(k)); b = np.full(n, be); b[J] = ga
        s = 0; tot = 0
        for X in itertools.combinations(range(n), k):
            w = abs(np.linalg.det(U[np.ix_(X, J)]))**2; tot += w
            s += w * np.prod([a[i] + (ga if i in X else be) for i in range(n)])
        worst = max(worst, abs(s - det_sum(a, b, U)) / abs(s), abs(tot - 1))
print("(1) two-value explicit Pluecker weights: max rel error over n=2..7, all k: %.1e" % worst)
# (2) embedding of Drury's 4x4 violator kills every pattern with >= 4 parts
U4 = np.load('moc/runs/deepest4.npy')
def lp_gap(pattern, U):
    reps, keys, V, blocks = perm_stats_matrix(pattern); n = sum(pattern)
    st = stats(U, blocks, n); p = np.array([st.get(k_, 0.0) for k_ in keys]); m, d = V.shape
    c = np.zeros(m + 1); c[-1] = 1
    A_ub = np.vstack([np.hstack([V.T, -np.ones((d, 1))]), np.hstack([-V.T, -np.ones((d, 1))])])
    r = linprog(c, A_ub=A_ub, b_ub=np.concatenate([p, -p]), A_eq=np.hstack([np.ones((1, m)), [[0]]]), b_eq=[1.0], bounds=[(0, None)] * (m + 1), method='highs')
    return r.fun
for pattern in [(1, 1, 1, 1), (2, 1, 1, 1), (2, 2, 1, 1), (3, 1, 1, 1), (1, 1, 1, 1, 1)]:
    n = sum(pattern); blocks = make_blocks(pattern)
    cols = [Bk[0] for Bk in blocks[:4]]              # one column from each of 4 blocks
    rest = [j for j in range(n) if j not in cols]
    U = np.zeros((n, n), complex)
    U[np.ix_(range(4), cols)] = U4
    for i, j in zip(range(4, n), rest): U[i, j] = 1
    assert np.allclose(U @ U.conj().T, np.eye(n))
    print("(2) pattern %-15s with Drury's U4 embedded: LP infeasibility gap %.3e (>0 => no certificate)" % (pattern, lp_gap(pattern, U)))
