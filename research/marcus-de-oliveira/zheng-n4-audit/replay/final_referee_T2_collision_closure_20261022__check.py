"""Audit replay for final_referee_T2_collision_closure_20261022 (branch of paper case
T_6 = {0,1,2,4,6,8}: cross-spectrum collisions).

Own code only; uses the shared library ../lib/moc4lib.py.  Nothing from the supplement is
read or executed.  Prints one line per check; exits nonzero if an EXACT check fails.
Run:  python -I -B final_referee_T2_collision_closure_20261022__check.py
"""
import os
import sys
import random
import itertools
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L  # noqa: E402

FAIL = []


def check(name, ok, info="", level="EXACT"):
    print(f"{'PASS' if ok else 'FAIL'} {level} {name} {info}".rstrip())
    if not ok and level == "EXACT":
        FAIL.append(name)


T = frozenset({0, 1, 2, 4, 6, 8})
check("T equals paper T_6", T == L.CONTACT_SETS[6])

# --- the document's m_j table: m_j(sigma) = 1[sigma(I) = J], 0-based labels ---
MTAB = {0: ("0", "3"), 3: ("03", "23"), 4: ("0", "2"), 5: ("01", "13"), 6: ("3", "0"),
        9: ("03", "12"), 10: ("02", "13"), 21: ("3", "1"), 23: ("03", "01")}
for j, (I, J) in MTAB.items():
    I0 = [int(ch) for ch in I]
    J0 = set(int(ch) for ch in J)
    prof = tuple(1 if set(s[i] for i in I0) == J0 else 0 for s in L.PERMS)
    check(f"m_{j} equals paper f_{j} (minor {L.f_name(j)})",
          prof == L.F_VALUES[j] and L.F_IS_MINOR[j])

# --- entries (i,j) used / unused by T ---
used = set((i, L.PERMS[s][i]) for s in T for i in range(4))
unused = sorted(set((i, j) for i in range(4) for j in range(4)) - used)
check("exactly four unused entries (0,2),(0,3),(3,0),(3,1)",
      unused == [(0, 2), (0, 3), (3, 0), (3, 1)], str(unused))

# --- used entry: the six-block sigma(i)=j is a contact set whose generators are all minor ---
allminor = True
for (i, j) in sorted(used):
    block = frozenset(s for s in range(24) if L.PERMS[s][i] == j)
    gens = L.generators(block | T)
    allminor &= all(L.F_IS_MINOR[k] for k in gens)
    gens_b = L.generators(block)
    allminor &= all(L.F_IS_MINOR[k] for k in gens_b)
check("used entry: generators of (six-block) and of (six-block u T) are all minor", allminor)

# --- the linear subspace W = {Pc : (Pc)_T = 0} and the coordinate claim ---
basis = []   # basis of W as 24-vectors, via nullspace of the T rows of P
NS = L.nullspace([L.P[s] for s in sorted(T)])
Wvals = [L.values(c) for c in NS]
check("dim {c : (Pc)_T = 0} = 8", len(NS) == 8)
COORD = [3, 7, 9, 10, 12, 13, 16, 18]
check("(g3,g7,g9,g10,g12,g13,g16,g18) are coordinates on W",
      L.rank([[w[k] for k in COORD] for w in Wvals]) == 8)


def W_from_q(q):
    """the unique g in W with g[COORD] = q"""
    M = [[Wvals[a][k] for a in range(8)] for k in COORD]
    lam = L.matvec(L.inv(M), q)
    return [sum((lam[a] * Wvals[a][s] for a in range(8)), Fraction(0)) for s in range(24)]


IDS = {
    (0, 2): {0: 19, 4: 12, 5: 10, 10: 7, 23: 3},
    (0, 3): {0: 18, 4: 12, 9: 7, 23: 3},
    (3, 0): {3: 12, 6: 9, 9: 7, 21: 3},
    (3, 1): {3: 12, 5: 10, 6: 15, 10: 7, 21: 3},
}
DIMS = {(0, 2): 5, (0, 3): 4, (3, 0): 4, (3, 1): 5}
for (i, j), dec in IDS.items():
    block = sorted(s for s in range(24) if L.PERMS[s][i] == j)
    # subspace of W where the six block gaps are equal
    rows = []
    for s in block[1:]:
        rows.append([Wvals[a][s] - Wvals[a][block[0]] for a in range(8)])
    K = L.nullspace(rows)
    check(f"A{i}=C{j}: dimension of equal-block subspace = {DIMS[(i, j)]}", len(K) == DIMS[(i, j)],
          f"got {len(K)}")
    ok = True
    for kv in K:
        g = [sum((kv[a] * Wvals[a][s] for a in range(8)), Fraction(0)) for s in range(24)]
        rhs = [sum((g[cs] * L.F_VALUES[f][s] for f, cs in dec.items()), Fraction(0))
               for s in range(24)]
        ok &= (g == rhs)
    check(f"A{i}=C{j}: decomposition holds on the whole subspace, all 24 coordinates", ok)
    check(f"A{i}=C{j}: all generators used are minor", all(L.F_IS_MINOR[f] for f in dec))

# --- realizability: an actual spectral collision lands in that subspace ---
rng = random.Random(7)


def rq():
    return L.QI(Fraction(rng.randint(-9, 9), rng.randint(1, 5)), Fraction(rng.randint(-9, 9), rng.randint(1, 5)))


# with the supplement's convention F = prod(A_i - C_sigma(i)), i.e. b = -C in the library
ok_real = True
for (i, j) in IDS:
    for trial in range(3):
        A = [rq() for _ in range(4)]
        C = [rq() for _ in range(4)]
        A[i] = C[j]
        Fv = L.sigma_points(A, [-x for x in C])
        blockvals = [Fv[s] for s in range(24) if L.PERMS[s][i] == j]
        ok_real &= all(x == 0 for x in blockvals)
check("a collision A_i = C_j makes the six block products zero (Gaussian rational spectra)", ok_real)

print(f"SUMMARY failures={len(FAIL)}")
sys.exit(1 if FAIL else 0)
