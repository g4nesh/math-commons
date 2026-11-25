"""Cross-check with Zheng (Oct 2026), Example 1.3, as reproduced in ../literature/README.md:
the functional psi coincides with one of our 48 exotic facets of the n=4 compound polytope,
and our facet code reproduces psi(O*) = -158/2401 for Zheng's rational orthogonal O*."""
import sys, numpy as np
sys.path.insert(0, 'moc')
from reduced4 import VERTS, PERMS
from bad4 import red_batch
F = np.load('moc/facets4.npy'); onv = np.array([(np.abs(VERTS @ f[:-1] + f[-1]) < 1e-9) for f in F])
O = np.array([[1, 4, 4, 4], [-4, 1, 4, -4], [-4, -4, 1, 4], [-4, 4, -4, 1]]) / 7.0
print("O^T O = I:", np.allclose(O.T @ O, np.eye(4)), " det =", round(np.linalg.det(O), 12))
m = lambda V, I, J: abs(np.linalg.det(V[np.ix_([i - 1 for i in I], [j - 1 for j in J])]))**2
psi = lambda V: 0.5 * (sum(abs(V[i, i])**2 for i in range(4)) - m(V, [1, 2], [1, 3]) - m(V, [1, 2], [2, 4])
                       + m(V, [1, 3], [1, 2]) + m(V, [1, 3], [3, 4]))
pv = [psi(np.eye(4)[list(s)]) for s in PERMS]
print("psi on permutation matrices: values/multiplicities", np.unique(np.round(pv, 9), return_counts=True))
print("psi(O*) = %.10f ; -158/2401 = %.10f" % (psi(O), -158 / 2401))
x = red_batch(O[None].astype(complex))[0]; gs = []
for i in np.where(onv.sum(1) == 13)[0]:
    f = F[i]; vals = VERTS @ f[:-1] + f[-1]; unit = np.abs(vals[np.abs(vals) > 1e-9]).min()
    gs.append(-(x @ f[:-1] + f[-1]) / unit)
    if np.allclose(-vals / unit, pv):
        print("psi coincides with exotic facet #%d on all 24 vertices (hence on the whole affine span)" % i)
print("min over the 48 exotic facets at O*: %.10f (facets violated: %d)" % (min(gs), sum(g < -1e-12 for g in gs)))
