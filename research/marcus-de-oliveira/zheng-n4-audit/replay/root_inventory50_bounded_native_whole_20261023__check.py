"""Exact replay for root_inventory50_bounded_native_whole_20261023 (paper case T_51 = {1,3,5,10,11,12}).

Checks the combinatorial claims of the one-page assembly: the contact set, the seven-element
inputs B, C, D (Table 5) and their stage order, that every possible new tie is covered by a
contact set of an earlier stage, and the 'eight exclusive coefficients' claim.
The analytic content is replayed in root_inventory50_bounded_native_conditional_20261023__check.py.
Run: python -I -B <this file>.  Exits 1 if an EXACT check fails.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
from moc4lib import *            # noqa

FAIL = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + ("  " + str(info) if info else ""))
    if not ok:
        FAIL.append(name)

STAGES6 = {1: "1 2 3 4 5 6 7 15 18 19 20 21 22 23 25 27 28 29 31 34 35 36 37 39 40 42 43 44 72 88",
           2: "10 26 46 90", 3: "9 11 12 47", 4: "8 13 17 49", 5: "14 16 24 33 38 70",
           6: "52 57 58 65 68 71 76 78", 7: "45 51 66 69 75 81", 8: "48 50 53 59 62 67 73 79 86",
           9: "30 41 54 55 60 61 63 74 77 80 83 84", 10: "56 64 82 85 87", 11: "32", 12: "89"}
STAGE = {int(k): st for st, v in STAGES6.items() for k in v.split()}
STAGE7 = {j: 1 for j in (1, 2, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 21, 22, 24)}
STAGE7.update({3: 2, 4: 2, 23: 2, 26: 2, 6: 4, 18: 4, 20: 5, 27: 5, 19: 7, 25: 7})

T = frozenset({1, 3, 5, 10, 11, 12})
NU = 51
check("literal contacts {1,3,5,10,11,12} = paper T_51", T == CONTACT_SETS[NU])
B7, C7, D7 = map(frozenset, ({0, 1, 2, 7, 9, 10, 16}, {0, 1, 2, 7, 9, 16, 17}, {0, 1, 2, 7, 9, 16, 23}))
for nm, S, k in (("B", B7, 8), ("C", C7, 11), ("D", D7, 12)):
    check("%s = Table 5 T%d^(7), stage %d < %d" % (nm, k, STAGE7[k], STAGE[NU]), AUX7[k] == S and STAGE7[k] < STAGE[NU])
gens = generators(T)
check("generators(T) = f0,f1,f2,f8,f10,f12,f30,f59 (eight, independent)",
      gens == [0, 1, 2, 8, 10, 12, 30, 59] and rank([F_VECS[j] for j in gens]) == 8, gens)
check("only non-minor generator is f30 and its coefficient is g4 (f30 is the only generator nonzero at 4)",
      [j for j in gens if not F_IS_MINOR[j]] == [30] and [j for j in gens if F_VALUES[j][4]] == [30] and F_VALUES[30][4] == 1)
check("ties not killing d: {0,6,7,13,14,16,20}", sorted(F_ZERO[30] - T) == [0, 6, 7, 13, 14, 16, 20])
for j in sorted(F_ZERO[30] - T):
    T2 = T | {j}
    prior6 = [nu for nu, S in CONTACT_SETS.items() if STAGE[nu] < STAGE[NU] and any(g_act_set(g, S) <= T2 for g in range(1152))]
    prior7 = [k for k, S in AUX7.items() if STAGE7[k] < STAGE[NU] and any(g_act_set(g, S) == T2 for g in range(1152))]
    check("tie %d covered by an earlier-stage contact set" % j, bool(prior6 or prior7), "six: %s seven: %s" % (prior6, prior7))
check("tie 6 and 13 are images of B, 14 of C, 20 of D",
      find_g(B7, T | {6}) is not None and find_g(B7, T | {13}) is not None and
      find_g(C7, T | {14}) is not None and find_g(D7, T | {20}) is not None)
print("EXACT checks failed:", len(FAIL))
sys.exit(1 if FAIL else 0)
