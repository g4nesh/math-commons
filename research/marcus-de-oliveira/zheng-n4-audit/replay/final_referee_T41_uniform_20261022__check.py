"""Auditor replay for final_referee_T41_uniform_20261022 (paper case T_42 = {1,2,3,6,7,8}).
Own code only; uses ../lib/moc4lib.py. Prints one line per check; exits 1 if an EXACT check fails."""
import os, sys, random
from fractions import Fraction as Fr
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import moc4lib as L

FAIL = []
def chk(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + (("  " + str(info)) if info else ""))
    if not ok: FAIL.append(name)

T = frozenset({1, 2, 3, 6, 7, 8})
chk("A1 contact set equals paper T_42", T == L.CONTACT_SETS[42])
gens = L.generators(T)
chk("B1 generators = 0,1,2,3,4,5,6,15,16,17,53,54", gens == [0, 1, 2, 3, 4, 5, 6, 15, 16, 17, 53, 54], gens)
chk("B1 non-minor generators exactly 17, 54", [k for k in gens if not L.F_IS_MINOR[k]] == [17, 54])
V = L.F_VALUES
def comb(lst):
    out = [0] * 24
    for k in lst:
        for s in range(24): out[s] += V[k][s]
    return out
IDS = [([17, 54], [3, 6, 6, 15, 53]), ([17, 0], [1, 3, 6, 15]), ([17, 2], [4, 6, 15]), ([17, 5], [3, 6, 16]),
       ([54, 1], [0, 6, 53]), ([54, 4], [2, 3, 6, 53]), ([54, 16], [5, 15, 6, 53])]
for lhs, rhs in IDS:
    chk("C2 identity %s = %s on all 24 values" % ("+".join("f%d" % k for k in lhs), "+".join("f%d" % k for k in rhs)),
        comb(lhs) == comb(rhs))
S17 = {1, 3, 4, 6, 15, 16, 17, 53}; S54 = {0, 2, 3, 5, 6, 15, 53, 54}
chk("C3 every right side of the f17 absorptions lies in the f17 canonical support",
    all(set(r) <= S17 for l, r in IDS[:4]))
chk("C3 every right side of the f54 absorptions lies in the f54 canonical support",
    all(set(r) <= S54 for l, r in [IDS[0]] + IDS[4:]))
chk("C3 f17 canonical support = gens minus {0,2,5,54}", S17 == set(gens) - {0, 2, 5, 54})
chk("C3 f54 canonical support = gens minus {1,4,16,17}", S54 == set(gens) - {1, 4, 16, 17})

# Run the reduction algorithm on random nonnegative coefficient vectors and check the final form
def reduce(lam):
    lam = dict(lam)
    def absorb(lhs, rhs):
        m = min(lam[k] for k in lhs)
        for k in lhs: lam[k] -= m
        for k in rhs: lam[k] += m
    absorb(*IDS[0])
    if lam[17] > 0:
        for l, r in IDS[1:4]: absorb(l, r)
    elif lam[54] > 0:
        for l, r in IDS[4:]: absorb(l, r)
    return lam
rng = random.Random(5)
ok = True
for trial in range(3000):
    lam = {k: (Fr(rng.randint(0, 6)) if rng.random() < 0.7 else Fr(0)) for k in gens}
    new = reduce(lam)
    if any(x < 0 for x in new.values()): ok = False; break
    before = [sum(lam[k] * V[k][s] for k in gens) for s in range(24)]
    after = [sum(new[k] * V[k][s] for k in gens) for s in range(24)]
    if before != after: ok = False; break
    supp = {k for k in gens if new[k] != 0}
    if not (new[17] == 0 and new[54] == 0) and not (supp <= S17 or supp <= S54): ok = False; break
chk("C3 reduction on 3000 random integer sources: coefficients stay >= 0, values unchanged, end in a minor sum or a canonical source", ok)

# Is there a symmetry of T_42 exchanging the two canonical sources (17 <-> 54)?
def g_facet(g, k):
    vals = L.g_act_values(g, list(V[k]))
    for kk in range(82):
        if list(V[kk]) == vals: return kk
    return None
stab = [g for g in L.G_ELEMENTS if L.g_act_set(g, T) == T]
swap = [g for g in stab if g_facet(g, 17) == 54]
info = ""
if swap:
    g = swap[0]
    info = "g=%s maps support S17 -> %s" % (g, sorted(g_facet(g, k) for k in S17))
chk("INFO stabilizer of T_42 has order %d; elements sending f17 to f54: %d" % (len(stab), len(swap)), True, info)
print("SUMMARY", "all EXACT checks passed" if not FAIL else "FAILED: %s" % FAIL)
sys.exit(1 if FAIL else 0)
