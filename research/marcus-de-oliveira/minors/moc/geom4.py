import sys, numpy as np, itertools
sys.path.insert(0, 'moc')
from bad4 import red_batch, EX, EXn, haar, deepen
from reduced4 import VERTS, PERMS
F = np.load('moc/facets4.npy')
onv = np.array([(np.abs(VERTS @ f[:-1] + f[-1]) < 1e-9) for f in F])
exo = np.where(onv.sum(1) == 13)[0]
# 1. are the 48 exotic facets one orbit under left/right permutations (and transpose)?
def permute_vec_action(Pl, Pr):
    # action on p_U of U -> Pl U Pr; compute via its effect on vertices: sigma -> composite
    pass
rng = np.random.default_rng(7)
U = haar(400000, 4, rng)
X = red_batch(U)
vals = X @ EXn[:, :-1].T + EXn[:, -1]
v = vals.max(1); bad = v > 0
print("bad fraction %.4f; violating >=2 exotic facets simultaneously: %d of %d" % (bad.mean(), ((vals > 0).sum(1) >= 2).sum(), bad.sum()))
Xb = X[bad]; vb = v[bad]; fb = vals[bad].argmax(1)
# distance to nearest vertex vs global scale
dv = np.min(np.linalg.norm(Xb[:, None, :] - VERTS[None, :, :], axis=2), axis=1)
diam = np.max(np.linalg.norm(VERTS[:, None] - VERTS[None], axis=2))
print("min distance from a bad p_U to a vertex: %.4f (polytope diameter %.3f)" % (dv.min(), diam))
# 2. barycentric coordinates of projection onto violated facet simplex
mins = []; maxs = []
for k in range(len(Xb)):
    fi = exo[fb[k]]; f = F[fi]; idx = np.where(onv[fi])[0]
    nrm = f[:-1] / np.linalg.norm(f[:-1]); off = f[-1] / np.linalg.norm(f[:-1])
    y = Xb[k] - (Xb[k] @ nrm + off) * nrm          # projection onto facet hyperplane
    A = np.vstack([VERTS[idx].T, np.ones(13)]); lam = np.linalg.lstsq(A, np.append(y, 1), rcond=None)[0]
    mins.append(lam.min()); maxs.append(lam.max())
mins = np.array(mins); maxs = np.array(maxs)
print("projection onto violated facet: inside the facet simplex for %.1f%% of bad points; min barycentric coord over all: %.4f; max coord over all: %.3f" % (100 * (mins >= -1e-12).mean(), mins.min(), maxs.max()))
print("max violation in this sample %.4f" % vb.max())
# deepest point structure
Ud = deepen(U[bad][np.argsort(-vb)[:50]].copy(), steps=600)
vd = (red_batch(Ud) @ EXn[:, :-1].T + EXn[:, -1]).max(1)
k = vd.argmax(); print("deepest after refinement: %.6f" % vd[k]); np.set_printoptions(precision=4, suppress=True)
print("|U|^2 of deepest violator:\n", np.abs(Ud[k])**2)
np.save('moc/runs/deepest4.npy', Ud[k])
