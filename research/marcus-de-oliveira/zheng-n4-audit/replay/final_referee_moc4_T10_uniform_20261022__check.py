# Replay for final_referee_moc4_T10_uniform_20261022 (primary proof of paper case T_11).
# Own code only; nothing from the supplement is read or run.  All checks EXACT (integer).
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
from moc4lib import (F_VALUES, F_IS_MINOR, P, rank, generators, CONTACT_SETS, find_g,
                     G_ELEMENTS, g_act_values, g_act_set)

FAILS = []
def chk(name, ok, level="EXACT"):
    print(("PASS" if ok else "FAIL"), level, name, flush=True)
    if level == "EXACT" and not ok:
        FAILS.append(name)
def lin(terms):
    out = [0]*24
    for c, j in terms:
        for s in range(24):
            out[s] += c*F_VALUES[j][s]
    return out
FV = lambda j: F_VALUES[j]

T = {0, 1, 2, 4, 7, 9}
g = find_g(T, CONTACT_SETS[11])
chk("document set {0,1,2,4,7,9} is mapped onto paper T_11 = {0,1,3,5,6,8} by g = %s" % (g,), g is not None)
GT = generators(T)
chk("generators(T) = {0,1,2,3,4,5,7,11,12,21,22,26}", GT == [0, 1, 2, 3, 4, 5, 7, 11, 12, 21, 22, 26])
chk("non-minor generators are exactly f22, f26", [j for j in GT if not F_IS_MINOR[j]] == [22, 26])
Av = {j: FV(j)[16] - FV(j)[3] for j in GT}
Bv = {j: -FV(j)[3] + FV(j)[6] - FV(j)[12] + FV(j)[13] + FV(j)[18] for j in GT}
chk("A = g16 - g3: +1 on f1,f4,f7,f22; -1 on f26; 0 otherwise",
    {j for j in GT if Av[j] == 1} == {1, 4, 7, 22} and {j for j in GT if Av[j] == -1} == {26} and all(Av[j] in (0, 1, -1) for j in GT))
chk("B: +1 on f0,f2,f12,f26; -1 on f22; 0 otherwise",
    {j for j in GT if Bv[j] == 1} == {0, 2, 12, 26} and {j for j in GT if Bv[j] == -1} == {22} and all(Bv[j] in (0, 1, -1) for j in GT))
chk("A + B >= 0 on every generator (sum is c0+c1+c2+c4+c7+c12)", all(Av[j] + Bv[j] >= 0 for j in GT))
ids = [("f22+f26 = f3+f5+f11+2f21", [(1, 22), (1, 26)], [(1, 3), (1, 5), (1, 11), (2, 21)]),
       ("f26+f1 = f0+f11+f21", [(1, 26), (1, 1)], [(1, 0), (1, 11), (1, 21)]),
       ("f26+f4 = f2+f3+f11+f21", [(1, 26), (1, 4)], [(1, 2), (1, 3), (1, 11), (1, 21)]),
       ("f26+f7 = f3+f5+f12+f21", [(1, 26), (1, 7)], [(1, 3), (1, 5), (1, 12), (1, 21)]),
       ("f22+f0 = f1+f3+f5+f21", [(1, 22), (1, 0)], [(1, 1), (1, 3), (1, 5), (1, 21)]),
       ("f22+f2 = f4+f5+f21", [(1, 22), (1, 2)], [(1, 4), (1, 5), (1, 21)]),
       ("f22+f12 = f7+f11+f21", [(1, 22), (1, 12)], [(1, 7), (1, 11), (1, 21)])]
for name, L, R in ids:
    chk("minor replacement " + name + " (right side all minor)", lin(L) == lin(R) and all(F_IS_MINOR[j] for _, j in R))
cA = [(lambda g: g[23], 0), (lambda g: g[15], 2), (lambda g: g[14], 3), (lambda g: g[10], 5), (lambda g: g[17], 11),
      (lambda g: g[8], 12), (lambda g: g[16], 21), (lambda g: g[3] - g[16], 26)]
cB = [(lambda g: g[23], 1), (lambda g: g[18], 3), (lambda g: g[15], 4), (lambda g: g[21], 5), (lambda g: g[8], 7),
      (lambda g: g[6], 11), (lambda g: g[19], 21), (lambda g: g[3] - g[6] + g[12] - g[13] - g[18], 22)]
dec = lambda co: all(lin([(cf(list(FV(j))), k) for cf, k in co]) == list(FV(j)) for j in GT)
chk("A-branch identity g = g23 f0 + g15 f2 + g14 f3 + g10 f5 + g17 f11 + g8 f12 + g16 f21 + e f26 on all generators", dec(cA))
chk("B-branch identity g = g23 f1 + g18 f3 + g15 f4 + g21 f5 + g8 f7 + g6 f11 + g19 f21 + d f22 on all generators", dec(cB))
chk("generators span the 8-dim space {Pc : (Pc)_T = 0}",
    14 - rank([P[s] for s in sorted(T)]) == 8 and rank([list(FV(j)) for j in GT]) == 8)
chk("real-v stratum: generators of {0..5} u T all minor", all(F_IS_MINOR[j] for j in generators(set(range(6)) | T)))
chk("b = g8 = 0 boundary: T u {8} contains paper T_7 = {0,1,2,4,7,8}", CONTACT_SETS[7] <= frozenset(T | {8}))
fidx = {F_VALUES[j]: j for j in range(82)}
def gmap_f(gg):
    return {j: fidx[tuple(g_act_values(gg, list(F_VALUES[j])))] for j in range(82)}
def find_tuple(src, dst):
    for gg in G_ELEMENTS:
        mp = gmap_f(gg)
        if [mp[j] for j in src] == list(dst):
            return gg
    return None
g1 = find_tuple([1, 11, 28], [2, 12, 26])
chk("SOC orbit: some g maps (f1,f11,f28) to (f2,f12,f26): %s" % (g1,), g1 is not None)
g2 = find_tuple([4, 6, 2, 28], [4, 7, 3, 22])
chk("four-coefficient orbit: some g maps (f4,f6,f2,f28) to (f4,f7,f3,f22) with roles (a,B,c,d): %s" % (g2,), g2 is not None)
chk("image of paper T_10 under that g: %s" % sorted(g_act_set(g2, CONTACT_SETS[10])), True, "INFO")
print("RESULT", "FAIL" if FAILS else "ALL EXACT CHECKS PASS", FAILS)
sys.exit(1 if FAILS else 0)
