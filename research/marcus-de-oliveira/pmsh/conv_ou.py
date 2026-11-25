"""Is conv{p_U : U in U(4)} = conv{p_O : O in O(4)} ?  Compare support functions h(l) = max l.p over U(4) and over O(4)
for random linear functionals l on the compound coordinates (B, D)."""
import sys, numpy as np
from scipy.optimize import minimize
from scipy.linalg import expm
from qtriv import xvec_U, NX
iu = np.triu_indices(4, 1)
def U_of(x):
    H = np.zeros((4, 4), complex); H[iu] = x[:6] + 1j * x[6:12]; H = H - H.conj().T; H[np.diag_indices(4)] = 1j * x[12:16]; return expm(H)
def O_of(x):
    K = np.zeros((4, 4)); K[iu] = x; K = K - K.T; return expm(K).astype(complex)
def support(l, f, dim, R, rng):
    best = -np.inf
    for _ in range(R):
        r = minimize(lambda x: -(l @ xvec_U(f(x))), rng.uniform(-np.pi, np.pi, dim), method='BFGS', options={'gtol': 1e-10})
        best = max(best, -r.fun)
    # also all permutation-matrix values are attained by both groups (signed permutations are orthogonal)
    return best
seed = int(sys.argv[1]); K = int(sys.argv[2]); rng = np.random.default_rng(seed); gaps = []
for k in range(K):
    l = rng.standard_normal(NX)
    hU = support(l, U_of, 16, 30, rng); hO = support(l, O_of, 6, 60, rng)
    gaps.append(hU - hO)
    print("functional %d: max_U %.10f  max_O %.10f  gap %.2e" % (k, hU, hO, hU - hO), flush=True)
print("RESULT seed %d: max gap (U minus O) over %d functionals = %.3e" % (seed, K, max(gaps)))
