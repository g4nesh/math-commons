"""Own-code replay for the Zheng n=4 audit, document final_referee_inventory75_whole_assembly_20261023
(assembly of paper case T_76). Checks the combinatorial interfaces between the two component proofs.
The component proofs are replayed in their own scripts. Nothing from the supplement is used."""
import sys, os, itertools, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L
NUMERIC = "--numeric" in sys.argv
FAIL = []
def chk(name, ok, extra=""):
    print(("PASS " if ok else "FAIL ") + name + ((" : " + str(extra)) if extra else ""), flush=True)
    if not ok: FAIL.append(name)
# Table 6 of the paper (order of proof), transcribed from the paper.
_ST6 = {1: [1,2,3,4,5,6,7,15,18,19,20,21,22,23,25,27,28,29,31,34,35,36,37,39,40,42,43,44,72,88],
        2: [10,26,46,90], 3: [9,11,12,47], 4: [8,13,17,49], 5: [14,16,24,33,38,70],
        6: [52,57,58,65,68,71,76,78], 7: [45,51,66,69,75,81], 8: [48,50,53,59,62,67,73,79,86],
        9: [30,41,54,55,60,61,63,74,77,80,83,84], 10: [56,64,82,85,87], 11: [32], 12: [89]}
_ST7 = {1: [1,2,5,7,8,9,10,11,12,13,14,15,16,17,21,22,24], 2: [3,4,23,26], 4: [6,18], 5: [20,27], 7: [19,25]}
STAGE6 = {k: s for s, ks in _ST6.items() for k in ks}
STAGE7 = {k: s for s, ks in _ST7.items() for k in ks}
REP6 = {L.orbit_rep(T): k for k, T in L.CONTACT_SETS.items()}
REP7 = {L.orbit_rep(T): k for k, T in L.AUX7.items()}
def closers(S, max_stage):
    """Paper cases (stage <= max_stage) or trivial (all-minor) orbits among the 6-subsets of S,
    and Table 5 sets equal to S up to G."""
    out = set()
    for sub in itertools.combinations(sorted(S), 6):
        if all(L.F_IS_MINOR[j] for j in L.generators(sub)):
            out.add(("trivial", 0)); continue
        k = REP6[L.orbit_rep(sub)]
        if STAGE6[k] <= max_stage: out.add(("T_%d" % k, STAGE6[k]))
    if len(S) == 7 and L.orbit_rep(S) in REP7:
        k = REP7[L.orbit_rep(S)]
        if STAGE7[k] <= max_stage: out.add(("T7_%d" % k, STAGE7[k]))
    return sorted(out)
def exclusive(gens, j):
    return [s for s in range(24) if L.F_VALUES[j][s] != 0 and all(L.F_VALUES[i][s] == 0 for i in gens if i != j)]
T = frozenset({0, 1, 3, 10, 11, 12}); S16 = T | {16}; S14 = T | {14}
chk("contact set equals paper T_76 literally", T == L.CONTACT_SETS[76])
g = L.generators(T)
chk("complete source [0,1,2,8,10,12,15,30], non-minor f30 with coefficient g2", g == [0, 1, 2, 8, 10, 12, 15, 30] and exclusive(g, 30) == [2] and [j for j in g if not L.F_IS_MINOR[j]] == [30])
chk("f1 has exclusive position 16 in the parent source (coefficient g16)", exclusive(g, 1) == [16])
chk("S16 source = parent source without f1", L.generators(S16) == [j for j in g if j != 1])
chk("S16 = Table 5 j=18 (stage 4 < 6)", S16 == L.AUX7[18] and STAGE7[18] < STAGE6[76])
gg = L.find_g(L.AUX7[27], S14)
chk("S14 = g.(Table 5 j=27) (stage 5 < 6)", gg is not None and STAGE7[27] < STAGE6[76], gg)
others = sorted(L.F_ZERO[30] - T - {14, 16})
chk("all other seventh zeros of f30 closed by earlier six-cases", all(any(c[0].startswith("T_") for c in closers(T | {x}, 5)) for x in others), others)
pos = [x for x in range(24) if x not in T and L.F_VALUES[30][x] > 0]
chk("every other extra tie has f30 > 0 there, forcing d = 0", all(L.F_VALUES[30][x] > 0 for x in pos) and sorted(pos + others + [14, 16]) == sorted(set(range(24)) - T))
print("SUMMARY failures:", FAIL)
sys.exit(1 if FAIL else 0)
