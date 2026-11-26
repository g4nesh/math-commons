"""Describe the sampled zero set: which entries and 2x2 minors vanish at each zero, grouped into patterns."""
import pickle, itertools, sys, numpy as np
Z = pickle.load(open(sys.argv[1], 'rb'))['Z']
P2 = list(itertools.combinations(range(4), 2))
def pattern(U, tol=1e-5):
    ent = tuple(int(abs(U[i, j]) < tol) for i in range(4) for j in range(4))
    mn = tuple(int(abs(np.linalg.det(U[np.ix_(S, T)])) < tol) for S in P2 for T in P2)
    return ent, mn
pats = {}
for U in Z: pats.setdefault(pattern(U), []).append(U)
print("%d zeros, %d distinct vanishing patterns" % (len(Z), len(pats)))
nm = lambda S: "".join(str(i + 1) for i in S)
for (ent, mn), Us in sorted(pats.items(), key=lambda kv: -len(kv[1])):
    E = [(i + 1, j + 1) for i in range(4) for j in range(4) if ent[4 * i + j]]
    Mz = ["%s|%s" % (nm(S), nm(T)) for a, S in enumerate(P2) for b, T in enumerate(P2) if mn[6 * a + b]]
    print("  x%-4d zero entries %s ; zero minors: %s" % (len(Us), E, " ".join(Mz)))
