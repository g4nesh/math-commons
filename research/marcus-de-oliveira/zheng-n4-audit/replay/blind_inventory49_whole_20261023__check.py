"""Exact replay for blind_inventory49_whole_20261023 (paper case T_50 = {1,3,5,6,10,11}).

The document assembles: the conditional reduction (blind_inventory49_bounded_native_conditional,
replayed in its own script), the seven-element cases F, G, A, B, the six-element cases cited
for new ties, and the common lemmas.  This script checks the combinatorial content stated in
the document and the stage order of Appendix G (Table 6) for every cited contact set.
Run: python -I -B <this file>.  Exits 1 if an EXACT check fails.
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
from moc4lib import *            # noqa

FAIL = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + ("  " + str(info) if info else ""))
    if not ok:
        FAIL.append(name)

# Table 6 of the paper (stages), transcribed from the paper text (p. 63)
STAGES6 = {1: "1 2 3 4 5 6 7 15 18 19 20 21 22 23 25 27 28 29 31 34 35 36 37 39 40 42 43 44 72 88",
           2: "10 26 46 90", 3: "9 11 12 47", 4: "8 13 17 49", 5: "14 16 24 33 38 70",
           6: "52 57 58 65 68 71 76 78", 7: "45 51 66 69 75 81", 8: "48 50 53 59 62 67 73 79 86",
           9: "30 41 54 55 60 61 63 74 77 80 83 84", 10: "56 64 82 85 87", 11: "32", 12: "89"}
STAGE = {int(k): st for st, v in STAGES6.items() for k in v.split()}
STAGE7 = {j: 1 for j in (1, 2, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 21, 22, 24)}
STAGE7.update({3: 2, 4: 2, 23: 2, 26: 2, 6: 4, 18: 4, 20: 5, 27: 5, 19: 7, 25: 7})

def act(row, col, tr):
    return (tuple(col), tuple(row), 1 if tr else 0)

T = frozenset({1, 3, 5, 6, 10, 11})
check("literal contact tuple is paper T_50", T == CONTACT_SETS[50])
check("Table 6 lists every T_k exactly once", sorted(STAGE) == list(range(1, 91)))
check("representative {1,3,4,6,7,8} -> T (row (0,2,3,1), column (1,0,2,3), no transpose)",
      g_act_set(act((0, 2, 3, 1), (1, 0, 2, 3), 0), {1, 3, 4, 6, 7, 8}) == T)
check("index.tsv style canonical form {0,1,2,7,9,10} is in the orbit", find_g({0, 1, 2, 7, 9, 10}, T) is not None)
check("statement 'contact rank 6, kernel dimension 8'", rank([P[s] for s in sorted(T)]) == 6)
check("kernel basis [f0,f1,f2,f3,f8,f10,f27,f59] independent and vanishing on T",
      rank([F_VECS[j] for j in (0, 1, 2, 3, 8, 10, 27, 59)]) == 8 and
      all(F_VALUES[j][s] == 0 for j in (0, 1, 2, 3, 8, 10, 27, 59) for s in T))
check("ten-source list equals generators(T)", generators(T) == [0, 1, 2, 3, 4, 8, 10, 27, 30, 59])
check("every source except 30 is a minor vector", [j for j in generators(T) if not F_IS_MINOR[j]] == [30])
check("d <= g4 and d <= g15 (f30 = 1 at 4 and 15, all sources >= 0)",
      F_VALUES[30][4] == 1 and F_VALUES[30][15] == 1 and all(F_VALUES[j][k] >= 0 for j in generators(T) for k in range(24)))
check("possible extra ties with f30 = 0 are exactly 0,7,12,13,14,16,20", sorted(F_ZERO[30] - T) == [0, 7, 12, 13, 14, 16, 20])

F7, G7 = frozenset({0, 1, 2, 7, 9, 10, 14}), frozenset({0, 1, 2, 7, 9, 10, 23})
A7, B7 = frozenset({0, 1, 2, 7, 9, 10, 17}), frozenset({0, 1, 2, 7, 9, 10, 16})
aux_index = {"F": 7, "G": 10, "A": 9, "B": 8}
for nm, S in (("F", F7), ("G", G7), ("A", A7), ("B", B7)):
    check("%s is Table 5 set T%d^(7), stage %d < stage(T_50)=%d" % (nm, aux_index[nm], STAGE7[aux_index[nm]], STAGE[50]),
          AUX7[aux_index[nm]] == S and STAGE7[aux_index[nm]] < STAGE[50])
table = [(0, 10, {0, 1, 3, 5, 7, 8}, (0, 1, 2, 3), (0, 1, 3, 2)),
         (7, 43, {1, 2, 4, 6, 7, 8}, (0, 1, 3, 2), (0, 1, 3, 2)),
         (12, "B", B7, (0, 2, 3, 1), (0, 1, 2, 3)),
         (13, "F", F7, (0, 2, 3, 1), (0, 1, 2, 3)),
         (14, "A", A7, (0, 2, 3, 1), (0, 1, 2, 3)),
         (16, 66, {2, 4, 5, 6, 9, 12}, (0, 3, 1, 2), (0, 1, 2, 3)),
         (20, "G", G7, (0, 2, 3, 1), (0, 1, 2, 3))]
for j, nm, S, row, col in table:
    S = frozenset(S)
    if isinstance(nm, int):
        check("tie %d: document's 'T%d' set is paper T_%d, stage %d < %d" % (j, nm - 1, nm, STAGE[nm], STAGE[50]),
              S == CONTACT_SETS[nm] and STAGE[nm] < STAGE[50])
    img = g_act_set(act(row, col, 0), S)
    check("tie %d: action image of %s is inside T+{%d}" % (j, nm, j), img <= T | {j}, sorted(img))
# Bound (1) on exact random instances
import random
rng = random.Random(5)
ok = True
for _ in range(30):
    co = {j: Fr(rng.randint(0, 9)) for j in generators(T)}
    c = [sum(co[j] * F_VECS[j][k] for j in co) for k in range(14)]
    v = values(c)
    O = random_rational_SO4(rng)
    ok &= -Phi(c, O) <= max(Fr(0), -Phi(F_VECS[30], O)) * min(v[4], v[15])
check("(1) -L_O <= M min(g4,g15) on 30 exact random instances", ok)
print("EXACT checks failed:", len(FAIL))
sys.exit(1 if FAIL else 0)
