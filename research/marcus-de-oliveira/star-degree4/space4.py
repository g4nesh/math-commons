"""Degree-4 bookkeeping: monomials in the 16 entries of U, grouped by torus weight (row multiset, column multiset)."""
import itertools, numpy as np
DEG = 4
MONS = list(itertools.combinations_with_replacement(range(16), DEG))
MARR = np.array(MONS)
def weight(m):
    r = [0] * 4; c = [0] * 4
    for v in m: r[v // 4] += 1; c[v % 4] += 1
    return (tuple(r), tuple(c))
BLOCKS = {}
for i, m in enumerate(MONS): BLOCKS.setdefault(weight(m), []).append(i)
BLOCK_LIST = list(BLOCKS.values())
def mono_batch(Us):
    """Us: array (K,4,4) complex -> (K, 3876) values of all degree-4 monomials"""
    flat = np.asarray(Us).reshape(len(Us), 16)
    out = flat[:, MARR[:, 0]].copy()
    for k in range(1, DEG): out *= flat[:, MARR[:, k]]
    return out
def haar(rng, K):
    z = (rng.standard_normal((K, 4, 4)) + 1j * rng.standard_normal((K, 4, 4))) / np.sqrt(2)
    q, r = np.linalg.qr(z); d = np.diagonal(r, axis1=1, axis2=2); return q * (d / np.abs(d))[:, None, :]
def parts_batch(Us):
    B = np.abs(Us)**2
    def D(S, T): return np.abs(Us[:, S[0], T[0]] * Us[:, S[1], T[1]] - Us[:, S[0], T[1]] * Us[:, S[1], T[0]])**2
    g0 = B[:, 0, 0] + B[:, 2, 2] + D((0, 3), (2, 3)) - D((0, 1), (0, 3))
    return g0, D((0, 1), (0, 1)), D((0, 2), (0, 3))
if __name__ == '__main__':
    sizes = np.array([len(b) for b in BLOCK_LIST])
    print("degree-4 monomials %d ; weight blocks %d ; sizes %s ; symmetric pairs %d" % (len(MONS), len(BLOCK_LIST), dict(zip(*np.unique(sizes, return_counts=True))), int((sizes * (sizes + 1) // 2).sum())))
