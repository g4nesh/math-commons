import sys; sys.path.insert(0,'moc')
import numpy as np
from scipy.spatial import ConvexHull
from reduced4 import *
from moclib import haar_unitary
h = ConvexHull(VERTS)
eq = h.equations
# dedupe facets
eqr = np.round(eq/np.abs(eq[:,:-1]).max(1,keepdims=True), 8)
_, idx = np.unique(eqr, axis=0, return_index=True)
F = eq[idx]
print("vertices", len(VERTS), "affine dim", np.linalg.matrix_rank(VERTS-VERTS[0]), "facets", len(F))
# count vertices on each facet
onv = [(np.abs(VERTS@f[:-1]+f[-1])<1e-9).sum() for f in F]
print("vertices-per-facet histogram", np.unique(onv, return_counts=True))
np.save('moc/facets4.npy', F)
rng = np.random.default_rng(1)
viol = []
for t in range(20000):
    U = haar_unitary(4, rng)
    x = vec(*red_coords(U))
    v = F[:,:-1]@x + F[:,-1]
    viol.append(v.max())
viol=np.array(viol)
print("fraction of Haar U(4) outside compound polytope: %.4f, max violation %.4f"%((viol>1e-12).mean(), viol.max()))
