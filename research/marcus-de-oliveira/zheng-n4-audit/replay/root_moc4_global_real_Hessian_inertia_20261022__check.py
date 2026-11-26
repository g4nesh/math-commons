"""Exact replay for "Hessian inertia for spectrally realizable real matching vectors"
(root_moc4_global_real_Hessian_inertia_20261022).  Own code only; nothing from the
supplement is read or executed.  One printed line per check; exit status nonzero iff an
EXACT check fails.

The cubic is the paper's C of (A.2) (moc4lib.cubic_C).  Inertia and ranks are invariant
under a real linear change of coordinates and under C -> -C, so any normalization of the
"matching cubic" gives the same counts.

Checks (EXACT = exact rational arithmetic on specific instances; these are instances of the
theorem, not a proof of it):
  1. smooth point of X: rank H(m) = 8 at real rational matching vectors
  2. Lemma 1 instances: rank H(m1 + m2) <= 12 for pairs of rational matching vectors
     (exactly 12 for generic pairs)
  3. Theorem 2 instances: A = H(Re m) for Gaussian-rational spectra and scales has inertia
     (6, 6, 2), via the exact characteristic polynomial (all roots real) and Descartes' rule
  4. degenerate instances (repeated eigenvalues, a_i + b_j = 0, zero eigenvalues): p, n <= 6
  5. all 91 principal 12-minors of A are >= 0, at least one > 0, on the generic instances
  6. NUMERIC: 2000 random complex spectra and scales, eigenvalue sign counts
"""
import itertools
import os
import random
import sys
from fractions import Fraction as Fr

import flint
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
Hs = sp.hessian(Cexpr, ms)
# third-derivative tensor: H(m)_{ij} = sum_k T[i][j][k] m_k
T = [[[int(sp.diff(Hs[i, j], ms[k])) for k in range(14)] for j in range(14)] for i in range(14)]


def H_at(m):
    return [[sum(T[i][j][k] * m[k] for k in range(14) if T[i][j][k]) for j in range(14)] for i in range(14)]


def to_fmpq_mat(M):
    n_ = len(M)
    return flint.fmpq_mat(n_, n_, [flint.fmpq(Fr(x).numerator, Fr(x).denominator) for r in M for x in r])


def inertia(M):
    """(p, n, z) of a real symmetric rational matrix from its exact characteristic polynomial."""
    cp = to_fmpq_mat(M).charpoly()
    co = [Fr(int(c.p), int(c.q)) for c in cp.coeffs()]  # low -> high degree
    z = 0
    while co[z] == 0:
        z += 1
    co = co[z:]
    deg = len(co) - 1

    def sign_changes(seq):
        s = [x for x in seq if x != 0]
        return sum(1 for u, v in zip(s, s[1:]) if (u > 0) != (v > 0))

    p_ = sign_changes(co)
    n_neg = sign_changes([c * (-1) ** k for k, c in enumerate(co)])
    assert p_ + n_neg == deg  # all roots real and nonzero: Descartes is exact
    return p_, n_neg, z


def rank_q(M):
    return to_fmpq_mat(M).rank()


rng = random.Random(2022)


def rq(lo=-6, hi=6, den=4):
    return Fr(rng.randint(lo, hi), rng.randint(1, den))


# 1. smooth real points
r8 = []
for _ in range(5):
    a_ = [rq() for _ in range(4)]
    b_ = [rq() for _ in range(4)]
    m = L.matching_vector(a_, b_)
    r8.append(rank_q(H_at(m)))
check("rank H(m) = 8 at 5 random real rational matching vectors", all(r == 8 for r in r8), info=f"{r8}")

# 2. Lemma 1 instances
r12 = []
for _ in range(6):
    m1 = L.matching_vector([rq() for _ in range(4)], [rq() for _ in range(4)], rq(1, 5))
    m2 = L.matching_vector([rq() for _ in range(4)], [rq() for _ in range(4)], rq(-5, -1))
    r12.append(rank_q(H_at([x + y for x, y in zip(m1, m2)])))
check("Lemma 1 instances: rank H(m1 + m2) <= 12 for 6 random pairs of rational matching vectors",
      all(r <= 12 for r in r12), info=f"{r12}")
check("  and = 12 for all of them (generic pairs)", all(r == 12 for r in r12), "INFO")
# contrast: a random point of the ambient space has full rank (the bound is not vacuous)
rr = rank_q(H_at([rq() for _ in range(14)]))
check("contrast: H at a random rational point of V has rank 14", rr == 14, "INFO", info=f"{rr}")

# 3. Theorem 2 generic instances
generic = []
for _ in range(6):
    a_ = [L.QI(rq(), rq()) for _ in range(4)]
    b_ = [L.QI(rq(), rq()) for _ in range(4)]
    al = L.QI(rq(), rq(1, 6))
    m = L.matching_vector(a_, b_, al)
    A = H_at([L.re_part(x) for x in m])
    generic.append(A)
ins = [inertia(A) for A in generic]
# density input without Terracini: at one Gaussian-rational instance, rank_C H(m) = 8 and
# rank [A; B] = 14 (A = H(Re m), B = H(Im m)), i.e. ker H(m) cap ker H(conj m) = 0.  Each is a
# nonvanishing minor, i.e. a real polynomial in the 18 real parameters that is not identically
# zero, so both hold on a dense open set of parameters.
a_ = [L.QI(1, 2), L.QI(-3, 1), L.QI(Fr(1, 2), -1), L.QI(2, 3)]
b_ = [L.QI(0, 1), L.QI(4, -2), L.QI(-1, Fr(1, 3)), L.QI(Fr(5, 2), 0)]
m = L.matching_vector(a_, b_, L.QI(1, 1))
Am = H_at([x.re for x in m])
Bm = H_at([x.im for x in m])
real_emb = [ra + [-x for x in rb] for ra, rb in zip(Am, Bm)] + [rb + ra for ra, rb in zip(Am, Bm)]
check("instance: rank_C H(m) = 8 (real embedding [[A,-B],[B,A]] has rank 16)",
      flint.fmpq_mat(28, 28, [flint.fmpq(Fr(x).numerator, Fr(x).denominator) for r in real_emb for x in r]).rank() == 16)
stack = Am + Bm
check("instance: rank [A; B] = 14, so ker H(m) cap ker H(conj m) = 0",
      flint.fmpq_mat(28, 14, [flint.fmpq(Fr(x).numerator, Fr(x).denominator) for r in stack for x in r]).rank() == 14)
check("Theorem 2: inertia of H(Re m) = (6,6,2) at 6 Gaussian-rational spectra and scales",
      all(x == (6, 6, 2) for x in ins), info=f"{ins}")

# 4. degenerate instances
deg_cases = {
    "a1 = a2": ([L.QI(1, 1), L.QI(1, 1), L.QI(-2, 3), L.QI(5, -1)], [L.QI(0, 2), L.QI(3, 1), L.QI(-1, -1), L.QI(2, 0)]),
    "a1 = a2, b3 = b4": ([L.QI(1, 1), L.QI(1, 1), L.QI(-2, 3), L.QI(5, -1)], [L.QI(0, 2), L.QI(3, 1), L.QI(-1, -1), L.QI(-1, -1)]),
    "a1 + b1 = 0": ([L.QI(1, 1), L.QI(2, -1), L.QI(-2, 3), L.QI(5, -1)], [L.QI(-1, -1), L.QI(3, 1), L.QI(-1, -2), L.QI(2, 0)]),
    "a1+b1 = 0 and a2+b2 = 0": ([L.QI(1, 1), L.QI(2, -1), L.QI(-2, 3), L.QI(5, -1)], [L.QI(-1, -1), L.QI(-2, 1), L.QI(-1, -2), L.QI(2, 0)]),
    "a = (x,x,x,y)": ([L.QI(1, 2), L.QI(1, 2), L.QI(1, 2), L.QI(-3, 0)], [L.QI(0, 2), L.QI(3, 1), L.QI(-1, -1), L.QI(2, 1)]),
    "a scalar": ([L.QI(1, 2)] * 4, [L.QI(0, 2), L.QI(3, 1), L.QI(-1, -1), L.QI(2, 1)]),
    "real a, b (m real)": ([L.QI(1), L.QI(2), L.QI(-3), L.QI(5)], [L.QI(0), L.QI(7), L.QI(-1), L.QI(2)]),
    "one circle a": ([L.QI(1), L.QI(0, 1), L.QI(-1), L.QI(Fr(3, 5), Fr(4, 5))], [L.QI(0, 2), L.QI(3, 1), L.QI(-1, -1), L.QI(2, 1)]),
}
okdeg = True
rows = []
for name, (a_, b_) in deg_cases.items():
    m = L.matching_vector(a_, b_, L.QI(2, 1))
    A = H_at([L.re_part(x) for x in m])
    pn = inertia(A)
    rows.append(f"{name}: {pn}")
    okdeg &= pn[0] <= 6 and pn[1] <= 6
check("degenerate instances: positive and negative index <= 6", okdeg, info="; ".join(rows))

# 5. principal 12-minors
okpm = True
info = []
for A in generic[:4]:
    vals = []
    for drop in itertools.combinations(range(14), 2):
        keep = [i for i in range(14) if i not in drop]
        sub = [[A[i][j] for j in keep] for i in keep]
        vals.append(to_fmpq_mat(sub).det())
    npos = sum(1 for v in vals if v > 0)
    nneg = sum(1 for v in vals if v < 0)
    info.append(f"{npos} pos/{nneg} neg/{91 - npos - nneg} zero")
    okpm &= nneg == 0 and npos >= 1
check("all 91 principal 12-minors >= 0 and at least one > 0 (4 generic instances)", okpm, info="; ".join(info))

# 6. NUMERIC
import numpy as np  # noqa: E402

rs = np.random.default_rng(31)
Tn = np.array(T, dtype=float)
P_np = L.P_numpy()
worst_p = worst_n = 0
counts = {}
NS = 2000
for _ in range(NS):
    a_ = rs.normal(size=4) + 1j * rs.normal(size=4)
    b_ = rs.normal(size=4) + 1j * rs.normal(size=4)
    if rs.random() < 0.3:  # force a coincidence sometimes
        a_[1] = a_[0]
    if rs.random() < 0.3:
        b_[2] = -a_[3]
    al = rs.normal() + 1j * rs.normal()
    m = np.array([complex(x) for x in L.matching_vector(list(a_), list(b_), al)])
    A = np.tensordot(Tn, m.real, axes=([2], [0]))
    ev = np.linalg.eigvalsh(A)
    tol = 1e-9 * max(1.0, np.abs(ev).max())
    pcount = int((ev > tol).sum())
    ncount = int((ev < -tol).sum())
    worst_p = max(worst_p, pcount)
    worst_n = max(worst_n, ncount)
    counts[(pcount, ncount)] = counts.get((pcount, ncount), 0) + 1
check(f"NUMERIC: {NS} random complex spectra/scales (seed 31, 30% forced coincidences): max p, max n <= 6",
      worst_p <= 6 and worst_n <= 6, "NUMERIC", info=f"inertia counts {counts}")

print("SUMMARY", "ALL EXACT CHECKS PASSED" if not FAIL else f"FAILED: {FAIL}")
sys.exit(1 if FAIL else 0)
