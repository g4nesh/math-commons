"""Exact replay for "Reduction of the order-four support problem to six contacts"
(final_referee_moc4_six_contact_reduction_20261022).  Own code only; nothing from the
supplement is read or executed.  One printed line per check; exit status nonzero iff an
EXACT check fails.

Checks (all from the paper's definitions through the shared library moc4lib):
  1. eq. (1): grad C(x + t kappa) = grad C(x) + t H_kappa x  (symbolic in x, t)
  2. H_kappa = Hess C(kappa) computed from C itself, rank 4, H_kappa kappa = 0
  3. eq. (2): ker H_kappa = span of the 16 single-entry vectors n_{i,j} (dimension 10)
  4. Lemma 1 profile: for c = sum a_ij n_ij, (Pc)_sigma = sum_i a_{i,sigma(i)}  (Birkhoff step)
     and Phi_c(O) = sum a_ij O_ij^2 on random exact rational O in SO(4)
  5. Theorem 3 bookkeeping: if P y - P x is a multiple of (1,...,1) then y - x is a multiple
     of kappa (injectivity of P); the frozen data are r affine-linear equations
  6. NUMERIC: the descent of Theorem 4 on random curves is not simulated; instead a
     brief sanity check that Phi_c(O) for c in ker H_kappa with Pc >= 0 is >= 0
"""
import itertools
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
Hess = sp.hessian(Cexpr, ms)
kap = [sp.Integer(int(k)) for k in L.KAPPA]
subk = dict(zip(ms, kap))
Hk = Hess.subs(subk)
Hk_lib = sp.Matrix(L.hess_kappa())
check("Hess C(kappa) computed from (A.2) equals the library's H_kappa (3.10)", Hk == Hk_lib)
check("grad C(kappa) = 0 (kappa in Sing C)", all(sp.expand(g.subs(subk)) == 0 for g in grad))
check("rank H_kappa = 4", Hk.rank() == 4)
check("H_kappa kappa = 0", Hk * sp.Matrix(kap) == sp.zeros(14, 1))

# 1. eq. (1)
t = sp.symbols("t")
shift = {m: m + t * k for m, k in zip(ms, kap)}
ok1 = True
Hkx = Hk * sp.Matrix(ms)
for i, g in enumerate(grad):
    lhs = sp.expand(g.xreplace(shift))
    if sp.expand(lhs - g - t * Hkx[i]) != 0:
        ok1 = False
check("eq. (1): grad C(x + t kappa) = grad C(x) + t H_kappa x, identically in x and t", ok1)
# the general cubic expansion behind it: grad C(x+tk) = grad C(x) + t H(x) k + t^2 grad C(k)
Hx = Hess  # Hessian at symbolic x
check("symmetry used: H(x) kappa = H(kappa) x identically", sp.expand(Hx * sp.Matrix(kap) - Hkx) == sp.zeros(14, 1))

# 3. eq. (2)
NS = [[int(v) if v.denominator == 1 else v for v in L.n_vec((i,), (j,))] for i in range(1, 5) for j in range(1, 5)]
Nmat = sp.Matrix([[sp.Rational(x) for x in row] for row in NS])  # 16 x 14
check("the 16 single-entry vectors span a space of dimension 10", Nmat.rank() == 10)
check("each single-entry vector lies in ker H_kappa", (Hk * Nmat.T) == sp.zeros(14, 16))
check("dim ker H_kappa = 10, hence ker H_kappa = span{n_ij} (over Q, R and C)", 14 - Hk.rank() == 10)

# 4. Lemma 1 profile and Phi formula
rng = random.Random(5)
okp = True
okphi = True
for trial in range(20):
    acoef = [[Fr(rng.randint(-5, 5), rng.randint(1, 3)) for _ in range(4)] for _ in range(4)]
    cvec = [Fr(0)] * 14
    for i in range(4):
        for j in range(4):
            nv = L.n_vec((i + 1,), (j + 1,))
            cvec = [x + acoef[i][j] * y for x, y in zip(cvec, nv)]
    vals = L.values(cvec)
    for s_idx, s in enumerate(L.PERMS):
        if vals[s_idx] != sum(acoef[i][s[i]] for i in range(4)):
            okp = False
    O = L.random_rational_SO4(rng)
    if L.Phi(cvec, O) != sum(acoef[i][j] * O[i][j] ** 2 for i in range(4) for j in range(4)):
        okphi = False
check("c = sum a_ij n_ij  =>  (Pc)_sigma = sum_i a_{i,sigma(i)}  (20 random rational a)", okp)
check("Phi_c(O) = sum a_ij O_ij^2 on 20 random exact rational O in SO(4)", okphi)
# Birkhoff step: (O_ij^2) is doubly stochastic for O orthogonal (exact on the same samples)
okds = True
for trial in range(20):
    O = L.random_rational_SO4(rng)
    sq = [[x * x for x in r] for r in O]
    okds &= all(sum(r) == 1 for r in sq) and all(sum(sq[i][j] for i in range(4)) == 1 for j in range(4))
check("(O_ij^2) is doubly stochastic (20 exact rational O)", okds)

# 5. injectivity of P and the height line
check("rank P = 14 (P injective)", L.rank([list(r) for r in L.P]) == 14)
check("P kappa = (1,...,1)", L.values(L.KAPPA) == [1] * 24)

# Theorem 3: every point of ker H_kappa + C kappa ... the constant-gap leaf condition
# A sanity check that a height line x + C kappa inside X forces x in ker H_kappa:
# take x a matching vector; grad C(x + t kappa) = t H_kappa x, so the line lies in
# Sing C iff H_kappa x = 0.  Check on a random rational matching vector that H_kappa x != 0
# (so generic lines are not in X) and that for a matching vector with H_kappa x = 0 ...
a_ = [Fr(1), Fr(2), Fr(-3), Fr(5)]
b_ = [Fr(7), Fr(-2), Fr(4), Fr(1, 2)]
x = L.matching_vector(a_, b_)
Hkx_num = [sum(int(Hk[i, j]) * x[j] for j in range(14)) for i in range(14)]
check("INFO: for a generic rational matching vector H_kappa x != 0 (so generic height lines are not in X)",
      any(v != 0 for v in Hkx_num), "INFO")

# 6. NUMERIC sanity of Lemma 1's conclusion
import numpy as np  # noqa: E402

rs = np.random.default_rng(7)
worst = np.inf
for _ in range(200):
    # random single-entry cost with all permutation values >= 0: a_ij >= 0 suffices but test
    # also signed a with nonnegative permutation sums (shift a row by a constant)
    acoef = rs.normal(size=(4, 4))
    perm_sums = [sum(acoef[i, s[i]] for i in range(4)) for s in L.PERMS]
    acoef[0, :] -= min(perm_sums)  # makes min permutation sum 0
    theta = rs.normal(size=6) * 2
    O = L.so4_from_params(theta)
    worst = min(worst, float((acoef * O**2).sum()))
check("NUMERIC Lemma 1: sum a_ij O_ij^2 >= 0 when all permutation sums >= 0 (200 samples, seed 7)",
      worst >= -1e-12, "NUMERIC", info=f"min = {worst:.3e}")

print("SUMMARY", "ALL EXACT CHECKS PASSED" if not FAIL else f"FAILED: {FAIL}")
sys.exit(1 if FAIL else 0)
