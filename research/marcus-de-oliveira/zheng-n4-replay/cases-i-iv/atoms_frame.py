"""The audit's universal inequalities on SO(4) in the coordinates of ../../star-orthogonal (exact).

A linear functional of the squared minors of a 4x4 orthogonal matrix is determined by its values at the 24 permutation
matrices (the audit's matrix P is injective and the squared-minor vectors of SO(4) and of the permutation matrices have the
same linear span).  So two such functionals agree on SO(4) as soon as they agree at the permutation matrices, and that is a
finite integer check.

g0(O) := D[13,13] + D[14,14] + D[13,14] + D[13,23] - D[12,34],  D[S,T] = det O[S,T]^2  (1-based index sets),
is the functional of inequality (star) in ../../pmsh/equivalence.md.  For each family of zheng-n4-audit/PROOF.md 7.6 the script
finds every relabelling g of rows and columns (with or without transposition) for which the family's non-minor functional
Phi_k equals g0(g.O), and prints which minors of g.O the family's minor terms become.  Uses the audit library for the 82
vectors and nothing else.  Exit code 0 iff every family is matched."""
import itertools, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "zheng-n4-audit", "lib"))
import moc4lib as L
P2 = list(itertools.combinations(range(4), 2))
def perm_matrix(s): return [[1 if s[i] == j else 0 for j in range(4)] for i in range(4)]       # 1 in row i, column sigma(i)
def relabel(M, tr, p, q):
    if tr: M = [[M[j][i] for j in range(4)] for i in range(4)]
    return [[M[p[i]][q[j]] for j in range(4)] for i in range(4)]
def D(M, S, T): return (M[S[0]][T[0]] * M[S[1]][T[1]] - M[S[0]][T[1]] * M[S[1]][T[0]])**2
def g0(M): return D(M, (0, 2), (0, 2)) + D(M, (0, 3), (0, 3)) + D(M, (0, 2), (0, 3)) + D(M, (0, 2), (1, 2)) - D(M, (0, 1), (2, 3))
PM = [perm_matrix(s) for s in L.PERMS]
FRAMES = [(tr, p, q) for tr in (0, 1) for p in itertools.permutations(range(4)) for q in itertools.permutations(range(4))]
G0VAL = {fr: tuple(g0(relabel(M, *fr)) for M in PM) for fr in FRAMES}
def minor_name(j, fr):
    target = L.F_VALUES[j]
    for i in range(4):
        for k in range(4):
            if tuple(relabel(M, *fr)[i][k]**2 for M in PM) == target: return "O[%d,%d]" % (i + 1, k + 1)
    for S in P2:
        for T in P2:
            if tuple(D(relabel(M, *fr), S, T) for M in PM) == target: return "D[%d%d,%d%d]" % (S[0] + 1, S[1] + 1, T[0] + 1, T[1] + 1)
    return None
FAMILIES = [("SOC", "Phi_22 + 2|x_1 x_9| >= 0", 22, (1, 9)),
            ("old-H", "Phi_70 + Phi_2 + 2|x_12 x_2| >= 0", 70, (2, 12)),
            ("chain", "Phi_70 + 2|x_2|(|x_0| + |x_12|) >= 0", 70, (2, 0, 12)),
            ("four-coefficient", "Phi_54 + Phi_10 + 2(|x_0 x_6| + |x_0 x_10| + |x_6 x_10|) >= 0", 54, (10, 0, 6)),
            ("rank-two", "Phi_35 + Phi_12/4 + |x_12| sqrt(Phi_2 + Phi_29) >= 0", 35, (12, 2, 29))]
ok = True
nonminor = [j for j in range(82) if not L.F_IS_MINOR[j]]
covered = {j for j in nonminor if any(G0VAL[fr] == L.F_VALUES[j] for fr in FRAMES)}
print("non-minor vectors whose functional is g0 after a relabelling: %d of %d" % (len(covered), len(nonminor))); ok &= len(covered) == len(nonminor) == 48
for name, text, k, minors in FAMILIES:
    frames = [fr for fr in FRAMES if G0VAL[fr] == L.F_VALUES[k]]; forms = sorted({tuple(minor_name(j, fr) for j in minors) for fr in frames})
    ok &= len(frames) > 0 and all(None not in f for f in forms)
    print("%s: %s" % (name, text)); print("   relabellings with Phi_%d = g0: %d ; the minors x_%s become:" % (k, len(frames), ", x_".join(map(str, minors))))
    for f in forms: print("      " + ", ".join(f))
star = ("D[12,12]", "D[13,14]"); hit = star in {tuple(minor_name(j, fr) for j in (2, 12)) for fr in FRAMES if G0VAL[fr] == L.F_VALUES[70]}
print("old-H in the frame (x_2^2, x_12^2) = (D[12,12], D[13,14]) reads g0 + D[12,12] + 2 sqrt(D[12,12] D[13,14]) >= 0, which is inequality (star):", hit); ok &= hit
print("ALL MATCHED" if ok else "MISMATCH"); sys.exit(0 if ok else 1)
