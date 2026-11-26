"""Are the six spectrum-free inequalities attached to g0 equivalent under the stabiliser of g0 in the symmetry group
(row permutations x column permutations x transpose)?  Exact integer computation on functionals."""
import numpy as np
from base import *
g0 = np.zeros(NX, dtype=int); g0[iB(0,0)] += 1; g0[iB(2,2)] += 1; g0[iD((0,3),(2,3))] += 1; g0[iD((0,1),(0,3))] -= 1
PAIRS = [(((0,1),(0,1)), ((0,2),(0,3))), (((0,1),(0,1)), ((0,2),(1,2))), (((0,3),(0,3)), ((0,1),(0,2))),
         (((0,3),(0,3)), ((0,1),(1,3))), (((0,2),(0,2)), ((0,3),(0,1))), (((0,2),(0,2)), ((0,3),(2,3)))]
key0 = tuple(VX @ g0)
def jac(S, T):   # canonical representative under Jacobi D[S,T] = D[S^c,T^c]
    return min((tuple(sorted(S)), tuple(sorted(T))), (comp(S), comp(T)))
stab = []
for al in PERMS:
    for be in PERMS:
        for tr in (False, True):
            if tuple(VX @ act(g0, al, be, tr)) == key0: stab.append((al, be, tr))
print("stabiliser of g0 (as a functional on the polytope):", len(stab), "elements")
def image(pair, al, be, tr):
    out = []
    for (S, T) in pair:
        S2, T2 = [al[s] for s in S], [be[t] for t in T]
        out.append(jac(T2, S2) if tr else jac(S2, T2))
    return tuple(out)
canon = [tuple(jac(S, T) for (S, T) in p) for p in PAIRS]
orbits = []
for i, p in enumerate(canon):
    orb = {image(PAIRS[i], *g) for g in stab}
    orbits.append(orb)
classes = []
for i in range(6):
    for c in classes:
        if canon[i] in orbits[c[0]]: c.append(i); break
    else: classes.append([i])
print("orbits of the six (a,b) pairs under the stabiliser:", classes)
