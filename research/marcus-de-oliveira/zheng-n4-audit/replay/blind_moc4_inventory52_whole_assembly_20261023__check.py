"""Own replay for blind_moc4_inventory52_whole_assembly_20261023 (paper case T_53).
Combinatorial and cone checks from the definitions via the shared audit library.
Nothing from the supplement is executed or read."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import moc4lib as L

fails = 0
def check(name, ok, info=""):
    global fails
    print(("PASS " if ok else "FAIL ") + name + (" :: " + str(info) if info != "" else ""))
    if not ok:
        fails += 1

T6 = frozenset({0, 1, 4, 6, 9, 12})
T7 = T6 | {17}
FV = L.F_VALUES
check("T6 equals paper T_53 literally", T6 == L.CONTACT_SETS[53])
check("T7 equals Table 5 aux set 19", T7 == L.AUX7[19])
gen6 = L.generators(T6)
check("generators(T6) == [0,1,2,5,7,9,21,35]", gen6 == [0, 1, 2, 5, 7, 9, 21, 35], gen6)
check("only f35 non-minor among them", [j for j in gen6 if not L.F_IS_MINOR[j]] == [35])
excl = {0: 18, 1: 17, 2: 15, 5: 11, 7: 14, 9: 7, 21: 5, 35: 2}
M = [[FV[j][excl[k]] for k in gen6] for j in gen6]
check("evaluation at exclusive perms is identity", M == [[int(i == j) for j in range(8)] for i in range(8)])
# rank of contact rows of P and the cone dimension
check("rank of P restricted to T6 rows is 6", L.rank([L.P[s] for s in sorted(T6)]) == 6)
nz17 = [j for j in gen6 if FV[j][17] != 0]
check("at perm 17 only f1 is nonzero (tie 17 kills exactly f1)", nz17 == [1], nz17)
check("generators(T7) == [0,2,5,7,9,21,35]", L.generators(T7) == [0, 2, 5, 7, 9, 21, 35])
# simplicial cone: 8 generators independent
check("8 generators linearly independent", L.rank([L.F_VECS[j] for j in gen6]) == 8)
# noncontact perms where f35 > 0 vs zero
pos35 = sorted(s for s in range(24) if s not in T6 and FV[35][s] > 0)
zero35 = sorted(s for s in range(24) if s not in T6 and FV[35][s] == 0)
check("non-contact perms with f35 = 0 are exactly {5,7,11,14,15,17,18}", zero35 == [5, 7, 11, 14, 15, 17, 18], (zero35, pos35))

# all-minor orbits among 6-sets
orbit_case = {L.orbit_rep(T): k for k, T in L.CONTACT_SETS.items()}
def six_subsets_info(S):
    import itertools
    out = []
    for sub in itertools.combinations(sorted(S), 6):
        rep = L.orbit_rep(sub)
        k = orbit_case.get(rep)
        allminor = all(L.F_IS_MINOR[j] for j in L.generators(sub))
        out.append((sub, k, allminor))
    return out

print("--- extra ties on the six-contact face T6 (native_harmonic (2)); paper cases contained in T6+{j}")
for j in [5, 7, 11, 14, 15, 17, 18]:
    info = six_subsets_info(T6 | {j})
    cases = sorted({k for (_, k, am) in info if k is not None and k != 53})
    am = [sub for (sub, k, a) in info if a]
    print("INFO tie", j, ": paper cases (other than 53) contained:", cases, "; all-minor 6-subsets:", len(am))
    if j != 17:
        check("tie %d: T6+{j} contains an earlier paper case or an all-minor 6-set" % j, bool(cases) or bool(am))

print("--- zero coefficients on the seven-contact face T7 (whole_assembly table)")
zero_tab = {"a(f0)": 18, "B(f2)": 15, "C(f5)": 11, "b(f7)": 14, "A(f9)": 7, "c(f21)": 5}
for nm, j in zero_tab.items():
    gk = {"a(f0)": 0, "B(f2)": 2, "C(f5)": 5, "b(f7)": 7, "A(f9)": 9, "c(f21)": 21}[nm]
    check("coefficient %s is the exclusive gap at %d on T7 cone" % (nm, j),
          FV[gk][j] == 1 and all(FV[k][j] == 0 for k in L.generators(T7) if k != gk))
    info = six_subsets_info(T7 | {j})
    cases = sorted({k for (_, k, am) in info if k is not None and k != 53})
    print("INFO zero", nm, "tie", j, ": paper cases contained in T7+{j}:", cases)
    check("zero %s: T7+{%d} contains an earlier paper case" % (nm, j), bool(cases))

# four-coefficient atom transports (rep (f0,f6;f10;f54), roles (a,b unordered; c; d))
def g_act_f(g, j):
    t = L.g_index_map(g)
    w = [None] * 24
    for s in range(24):
        w[t[s]] = FV[j][s]
    w = tuple(w)
    return FV.index(w) if w in FV else None

rep = (0, 6, 10, 54)
targets = [(7, 21, 9, 35), (0, 21, 2, 35), (0, 7, 5, 35)]
for tgt in targets:
    found = None
    for g in L.G_ELEMENTS:
        im = [g_act_f(g, j) for j in rep]
        if set(im[:2]) == set(tgt[:2]) and im[2] == tgt[2] and im[3] == tgt[3]:
            found = g
            break
    check("four-coefficient atom (0,6;10;54) maps to %s with roles (a,b;c;d)" % (tgt,), found is not None, found)

# S3 stabilizer of T7 permuting pairs (f0,f9),(f7,f2),(f21,f5) and fixing f35
pairs = [(0, 9), (7, 2), (21, 5)]
induced = set()
ok_pairs = True
for g in L.G_ELEMENTS:
    if L.g_act_set(g, T7) != T7:
        continue
    im = {j: g_act_f(g, j) for j in L.generators(T7)}
    if im[35] != 35:
        ok_pairs = False
        continue
    perm = []
    for (m, o) in pairs:
        tgt = [i for i, (m2, o2) in enumerate(pairs) if im[m] == m2 and im[o] == o2]
        if len(tgt) != 1:
            ok_pairs = False
        perm.append(tgt[0] if tgt else None)
    induced.add(tuple(perm))
check("stabilizer of T7 fixes f35 and permutes the pairs (main, opposite) consistently", ok_pairs)
check("induced action is the full S3", len(induced) == 6, sorted(induced))

# NUMERIC: realizable points on the seven-contact face T7; region claims (2)-(3) and min Phi
import numpy as np
excl7 = {0: 18, 2: 15, 5: 11, 7: 14, 9: 7, 21: 5, 35: 2}
nfound = 0; worst_reg = None; worst_phi = None; viol = []
for seed in range(40):
    pt = L.random_point_on_contact_stratum(T7, seed=seed, tries=8, margin=1e-3)
    if pt is None:
        continue
    g = np.asarray(pt["gaps"]).real
    co = {k: g[e] for k, e in excl7.items()}
    a_, B_, C_, b_, A_, c_, d_ = co[0], co[2], co[5], co[7], co[9], co[21], co[35]
    if min(co.values()) <= 1e-6:
        continue
    nfound += 1
    trip = sorted([(a_, A_), (b_, B_), (c_, C_)])   # sort main coefficients with their opposite partners
    (a1, A1), (b1, B1), (c1, C1) = trip
    dd = d_ * d_
    ok_low = max(a1 * b1, b1 * c1, a1 * c1) > dd
    capA = d_ * (b1 + d_) * (c1 + d_) / (b1 * c1 - dd) if b1 * c1 > dd else np.inf
    capB = d_ * (a1 + d_) * (c1 + d_) / (a1 * c1 - dd) if a1 * c1 > dd else np.inf
    capC = d_ * (a1 + d_) * (b1 + d_) / (a1 * b1 - dd) if a1 * b1 > dd else np.inf
    if a1 * c1 <= dd < b1 * c1:
        ok = A1 >= capA * (1 - 1e-7)
    elif a1 * b1 <= dd < a1 * c1:
        ok = A1 >= capA * (1 - 1e-7) or B1 >= capB * (1 - 1e-7)
    elif dd < a1 * b1:
        ok = A1 >= capA * (1 - 1e-7) or B1 >= capB * (1 - 1e-7) or C1 >= capC * (1 - 1e-7)
    else:
        ok = False
    if not (ok_low and ok):
        viol.append((seed, co))
    m14 = L.matching_vector_np(pt["a"], pt["b"], pt["alpha"])
    cvec = np.real(m14) - np.real(L.F_values_np(pt["a"], pt["b"], pt["alpha"])[pt["pi_star"]]) * np.eye(14)[13]
    val, _ = L.min_phi_over_SO4(cvec / max(1.0, np.abs(cvec).max()), restarts=6, seed=seed)
    worst_phi = val if worst_phi is None else min(worst_phi, val)
print("NUMERIC T7 stratum: %d strictly positive realizable samples (seeds 0..39, margin 1e-3); region-claim violations: %d; min Phi over SO(4) (6 restarts each, scaled) = %s" % (nfound, len(viol), worst_phi))
if viol:
    print("NUMERIC VIOLATION DATA", viol[:3])
print("RESULT", "all EXACT checks passed" if fails == 0 else "%d failures" % fails)
sys.exit(1 if fails else 0)
