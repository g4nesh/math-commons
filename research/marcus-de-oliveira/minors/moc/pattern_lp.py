"""Pattern compound-Birkhoff test.
B's eigenvalues take value c_j on column block B_j (sizes = pattern).  d is a polynomial
in (a, c) whose coefficients are grouped statistics
   G_U(S, kappa) = sum_{T : |T cap B_j| = kappa_j} |det U[S,T]|^2 .
A certificate for MOC is a probability lambda on S_n with  sum lambda_s G_s = G_U."""
import sys, itertools, numpy as np
from scipy.optimize import linprog
sys.path.insert(0, 'moc')
from moclib import haar_unitary

def stats(U, blocks, n):
    lab = np.empty(n, int)
    for j, Bk in enumerate(blocks):
        lab[list(Bk)] = j
    m = len(blocks)
    out = {}
    for k in range(1, n):
        for S in itertools.combinations(range(n), k):
            for T in itertools.combinations(range(n), k):
                kap = tuple(np.bincount(lab[list(T)], minlength=m))
                v = abs(np.linalg.det(U[np.ix_(S, T)]))**2
                out[(S, kap)] = out.get((S, kap), 0.0) + v
    return out

def make_blocks(pattern):
    blocks, s = [], 0
    for p in pattern:
        blocks.append(tuple(range(s, s+p))); s += p
    return blocks

def perm_stats_matrix(pattern):
    n = sum(pattern); blocks = make_blocks(pattern)
    # classes of permutations: sigma maps row i to column sigma(i); only block of sigma(i) matters
    lab = np.empty(n, int)
    for j, Bk in enumerate(blocks): lab[list(Bk)] = j
    classes = {}
    for s in itertools.permutations(range(n)):
        key = tuple(lab[list(s)])
        if key not in classes: classes[key] = s
    reps = list(classes.values())
    keys = None; rows = []
    for s in reps:
        P = np.zeros((n, n)); P[range(n), s] = 1
        st = stats(P, blocks, n)
        if keys is None: keys = sorted(st.keys())
        rows.append([st.get(k_, 0.0) for k_ in keys])
    return reps, keys, np.array(rows), blocks

def test(pattern, trials, seed=0):
    n = sum(pattern)
    reps, keys, V, blocks = perm_stats_matrix(pattern)
    rng = np.random.default_rng(seed)
    m, d = V.shape
    worst = 0; bad = 0
    for _ in range(trials):
        U = haar_unitary(n, rng)
        st = stats(U, blocks, n)
        p = np.array([st.get(k_, 0.0) for k_ in keys])
        c = np.zeros(m+1); c[-1] = 1
        A_ub = np.vstack([np.hstack([V.T, -np.ones((d,1))]), np.hstack([-V.T, -np.ones((d,1))])])
        b_ub = np.concatenate([p, -p])
        A_eq = np.hstack([np.ones((1,m)), np.zeros((1,1))])
        r = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=[1.0], bounds=[(0,None)]*(m+1), method='highs')
        worst = max(worst, r.fun); bad += r.fun > 1e-9
    return worst, bad, m

if __name__ == '__main__':
    pats = [tuple(map(int, a.split(','))) for a in sys.argv[2:]]
    trials = int(sys.argv[1])
    for pat in pats:
        w, bad, m = test(pat, trials)
        print("pattern %-14s classes=%4d  trials=%d  infeasible=%d  worst Linf gap=%.2e" % (pat, m, trials, bad, w), flush=True)
