"""Exact check (sympy, integers) that the proven unitary identities (i)-(iv) cut out exactly the affine span of the
24 permutation compound vectors: rank(E) = 39, so the solution set has affine dimension 52 - 39 = 13 = dim aff(Q)."""
import sympy as sp, numpy as np
from qtriv import E, e, VX
Ei = sp.Matrix(E.astype(int)); assert np.allclose(E, np.array(Ei.tolist(), float))
rE = Ei.rank()
Vi = sp.Matrix(VX.astype(int)); D = Vi[1:, :] - sp.ones(23, 1) * Vi[0, :]; rQ = D.rank()
ok = all(list(Ei * Vi.row(k).T) == [int(v) for v in e] for k in range(24))
print("exact rank of identities:", rE, "-> affine dim", 52 - rE, "; exact affine dim of Q:", rQ, "; all 24 vertices satisfy identities:", ok)
