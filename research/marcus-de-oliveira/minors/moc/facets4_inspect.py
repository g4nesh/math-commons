import sys; sys.path.insert(0,'moc')
import numpy as np
from fractions import Fraction
from reduced4 import *
from moclib import haar_unitary
F = np.load('moc/facets4.npy')
def full_from_vec_coeffs(f):
    return f
onv = np.array([(np.abs(VERTS@f[:-1]+f[-1])<1e-9) for f in F])
for cnt in [20,13]:
    rows = [i for i in range(len(F)) if onv[i].sum()==cnt]
    print("=== facets with",cnt,"vertices:",len(rows))
    for i in rows[:3]:
        off = [PERMS[j] for j in range(24) if not onv[i,j]]
        vals = VERTS@F[i,:-1]+F[i,-1]
        print(" excluded perms:", off, " values on excluded:", np.round(vals[~onv[i]]/np.abs(vals[~onv[i]]).min(),4))
# which facets violated
rng = np.random.default_rng(1)
cnts = np.zeros(len(F))
for t in range(20000):
    U = haar_unitary(4, rng)
    x = vec(*red_coords(U))
    v = F[:,:-1]@x + F[:,-1]
    cnts += v>1e-12
print("violation counts by facet type (vertices on facet -> total violations):")
for cnt in [13,18,20]:
    rows = onv.sum(1)==cnt
    print(cnt, cnts[rows].sum())
