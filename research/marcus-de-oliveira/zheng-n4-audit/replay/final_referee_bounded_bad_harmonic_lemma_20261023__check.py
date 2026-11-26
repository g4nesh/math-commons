"""Replay for "Bounded coefficients and removable negative sublevel sets"
(final_referee_bounded_bad_harmonic_lemma_20261023).  Own code only; nothing from the
supplement is read or executed.  One printed line per check; exit status nonzero iff an
EXACT check fails.

Theorem 1 is an analytic lemma; it is audited by reading (see the report).  This script
replays the finite algebraic facts the document uses in its "height exception" paragraph,
and the one-line algebra of the application paragraph:
  1. grad C(v0 + f kappa) = grad C(v0) + f H_kappa v0 for the paper's cubic (A.2),
     so v0 + f kappa in Sing C for two values of f forces grad C(v0) = 0 = H_kappa v0
  2. ker H_kappa = span of the 16 single-entry vectors ("singleton-assignment span")
  3. a cost sum_ij a_ij n_ij with all permutation values >= 0 has Phi >= 0 (Birkhoff), exact
     on rational O via the identity Phi = sum a_ij O_ij^2 and doubly stochastic (O_ij^2)
  4. application: c = sum_k lam_k f_k (minor, lam_k >= 0) + d f_bad  =>  Phi_c(O) >= d Phi_fbad(O)
     (checked exactly on random rational O for every non-minor f_j as f_bad)
  5. NUMERIC illustration of the pole step: Im of z^-k (c + O(z)) takes both signs on small circles
"""
import os
import random
import sys
from fractions import Fraction as Fr

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L  # noqa: E402

FAIL = []


def check(name, ok, level="EXACT", info=""):
    print(f"[{level}] {'PASS' if ok else 'FAIL'}  {name}" + (f"  ({info})" if info else ""))
    if not ok and level == "EXACT":
        FAIL.append(name)


ms, Cexpr = L.cubic_sympy()
grad = [sp.diff(Cexpr, m) for m in ms]
kap = [sp.Integer(int(k)) for k in L.KAPPA]
Hk = sp.hessian(Cexpr, ms).subs(dict(zip(ms, kap)))
f = sp.symbols("f")
Hkv = Hk * sp.Matrix(ms)
ok = all(sp.expand(g.xreplace({m: m + f * k for m, k in zip(ms, kap)}) - g - f * Hkv[i]) == 0
         for i, g in enumerate(grad))
check("grad C(v0 + f kappa) = grad C(v0) + f H_kappa v0 identically (f^2 term is grad C(kappa) = 0)", ok)

Nmat = sp.Matrix([[sp.Rational(x) for x in L.n_vec((i,), (j,))] for i in range(1, 5) for j in range(1, 5)])
check("ker H_kappa (dim 10) = span of the 16 single-entry vectors",
      Hk.rank() == 4 and Nmat.rank() == 10 and Hk * Nmat.T == sp.zeros(14, 16))

rng = random.Random(23)
ok3 = True
ok4 = True
nonminor = [j for j in range(82) if not L.F_IS_MINOR[j]]
minor = [j for j in range(82) if L.F_IS_MINOR[j]]
for trial in range(10):
    O = L.random_rational_SO4(rng)
    sq = [[x * x for x in r] for r in O]
    ok3 &= all(sum(r) == 1 for r in sq) and all(sum(sq[i][j] for i in range(4)) == 1 for j in range(4))
    for jb in nonminor:
        lam = {j: Fr(rng.randint(0, 5), rng.randint(1, 3)) for j in rng.sample(minor, 5)}
        dd = Fr(rng.randint(0, 5), rng.randint(1, 3))
        c = [dd * x for x in L.F_VECS[jb]]
        for j, lv in lam.items():
            c = [x + lv * y for x, y in zip(c, L.F_VECS[j])]
        if not (L.Phi(c, O) >= dd * L.Phi(L.F_VECS[jb], O)):
            ok4 = False
    # minor vectors have Phi = det(O[I,J])^2 >= 0
    ok4 &= all(L.Phi(L.F_VECS[j], O) >= 0 for j in minor)
check("(O_ij^2) doubly stochastic on 10 exact rational O in SO(4)", ok3)
check("Phi_c(O) >= d Phi_fbad(O) for c = (nonneg. comb. of minor f_j) + d f_bad, d >= 0; all 48 non-minor f_bad, 10 exact O",
      ok4)

import numpy as np  # noqa: E402

rs = np.random.default_rng(1)
okpole = True
for _ in range(200):
    k = int(rs.integers(1, 5))
    cc = complex(rs.normal(), rs.normal())
    coeffs = rs.normal(size=6) + 1j * rs.normal(size=6)
    r = 1e-3
    th = np.linspace(0, 2 * np.pi, 2001)
    z = r * np.exp(1j * th)
    Hval = z ** (-k) * (cc + sum(coeffs[i] * z ** (i + 1) for i in range(6)))
    okpole &= (Hval.imag.max() > 0) and (Hval.imag.min() < 0)
check("NUMERIC pole step: Im H takes both signs on |z| = 1e-3 for 200 random Laurent heads (seed 1)", okpole, "NUMERIC")

print("SUMMARY", "ALL EXACT CHECKS PASSED" if not FAIL else f"FAILED: {FAIL}")
sys.exit(1 if FAIL else 0)
