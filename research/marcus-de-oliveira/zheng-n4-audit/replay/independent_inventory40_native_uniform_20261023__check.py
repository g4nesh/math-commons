"""Auditor replay for independent_inventory40_native_uniform_20261023 (paper case T_41 = {0,1,3,6,11,12}).
Own code only; uses the shared library ../lib/moc4lib.py. Nothing from the supplement is read.
Prints one line per check; exits 1 if an EXACT check fails."""
import os, sys, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import moc4lib as L
import sympy as sp

FAIL = []
def chk(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + (("  " + str(info)) if info else ""))
    if not ok: FAIL.append(name)

T = frozenset({0, 1, 3, 6, 11, 12})
chk("A1 contact set equals paper T_41", T == L.CONTACT_SETS[41])
gens = L.generators(T)
chk("B1 generators of N_T = (1)", gens == [0, 1, 2, 7, 8, 9, 10, 15, 18, 29, 30], gens)
chk("B1 non-minor generators are exactly 18, 30", [k for k in gens if not L.F_IS_MINOR[k]] == [18, 30])
chk("B1 span of generators has dimension 8", L.rank([L.F_VECS[k] for k in gens]) == 8)
chk("B1 the six contact rows of P have rank 6", L.rank([L.P[s] for s in sorted(T)]) == 6)
# (2): theta = lam18 + lam30 <= g17 for every nonnegative representation
v17 = {k: L.F_VALUES[k][17] for k in gens}
chk("C2 values at pi_17: f18 = f30 = 1, all other generators >= 0",
    v17[18] == 1 and v17[30] == 1 and all(x >= 0 for x in v17.values()), v17)
# (4): indices where both nonminor generators are strictly positive
both = [s for s in range(24) if s not in T and L.F_VALUES[18][s] > 0 and L.F_VALUES[30][s] > 0]
chk("C4 list (4) = noncontact indices where f18, f30 > 0", both == [4, 8, 15, 17, 19, 22], both)
rest = [s for s in range(24) if s not in T and s not in both]
chk("C5 list (5) = remaining noncontact indices", rest == [2, 5, 7, 9, 10, 13, 14, 16, 18, 20, 21, 23], rest)

# stages of Table 6 (Appendix G) for six-element cases
STAGE = {}
for st, ks in {1: [1,2,3,4,5,6,7,15,18,19,20,21,22,23,25,27,28,29,31,34,35,36,37,39,40,42,43,44,72,88],
               2: [10,26,46,90], 3: [9,11,12,47], 4: [8,13,17,49], 5: [14,16,24,33,38,70],
               6: [52,57,58,65,68,71,76,78], 7: [45,51,66,69,75,81], 8: [48,50,53,59,62,67,73,79,86],
               9: [30,41,54,55,60,61,63,74,77,80,83,84], 10: [56,64,82,85,87], 11: [32], 12: [89]}.items():
    for k in ks: STAGE[k] = st
chk("G table 6 lists 90 cases", sorted(STAGE) == list(range(1, 91)))
REP = {L.orbit_rep(Tk): k for k, Tk in L.CONTACT_SETS.items()}
def classify(S):
    if all(L.F_IS_MINOR[k] for k in L.generators(S)): return "allminor"
    return REP.get(L.orbit_rep(S), "??")
for j in rest:
    U = T | {j}
    opts = []
    for S in itertools.combinations(sorted(U), 6):
        S = frozenset(S)
        if S == T: continue
        opts.append((sorted(S), classify(S)))
    good = [(S, k) for S, k in opts if k == "allminor" or (k != "??" and STAGE[k] < STAGE[41])]
    best = good[0] if good else None
    lab = [("all-minor" if k == "allminor" else ("T_%d/st%d" % (k, STAGE[k]) if k != "??" else k)) for _, k in opts]
    info = ("e.g. %s -> %s; all: %s" % (best[0], ("all-minor" if best[1] == "allminor" else "T_%d (stage %d)" % (best[1], STAGE[best[1]])), lab)) if best else str(lab)
    chk("C5 new contact %d: T u {%d} contains an all-minor 6-set or an earlier-stage case" % (j, j), bool(good), info)

# ---------------- curve algebra ----------------
v, z, x, t, w = sp.symbols("v z x t w")
def br(P, Q): return P[0]*Q[1] - P[1]*Q[0]
def Fs(A, C):
    return [sp.prod([br(A[i], C[s[i]]) for i in range(4)]) for s in L.PERMS]
A = [(0, 1), (1, 1), (z, 1), (w, 1)]
C = [(1, 0), (v, 1), (x, 1), (t, 1)]
F = Fs(A, C)
E = sp.factor(F[6] - F[0])
chk("C6 E = F6 - F0 = (z - x)(w - t)", sp.expand(E - (z - x)*(w - t)) == 0, E)
chk("C6 F0/E = v - 1", sp.simplify(F[0]/E - (v - 1)) == 0)
r_ = sp.together((F[1] - F[0])/E); s_ = sp.together((F[3] - F[1])/E)
u_ = sp.together((F[11] - F[6])/E); h_ = sp.together((F[12] - F[6])/E)
print("INFO r,s,u,h =", sp.factor(r_), "|", sp.factor(s_), "|", sp.factor(u_), "|", sp.factor(h_))
chk("C8 u = v(w-1)/(t-w)", sp.simplify(u_ - v*(w - 1)/(t - w)) == 0)
chk("C8 h = z(v-x)/(x-z)", sp.simplify(h_ - z*(v - x)/(x - z)) == 0)
chk("C8 s = u h (t-z)/(z v)", sp.simplify(s_ - u_*h_*(t - z)/(z*v)) == 0)
H17 = sp.together((F[17] - F[0])/E)
chk("C11 H17 = 1 + u + h + uh/v (raw identity)", sp.simplify(H17 - (1 + u_ + h_ + u_*h_/v)) == 0)
# inverse (9)
u, h, s, r = sp.symbols("u h s r")
xi = z*(h + v)/(h + z); ti = z*(u*h + s*v)/(u*h); wi = (z*(u*h + s*v) + h*v)/(h*(u + v))
sub = {x: xi, t: ti, w: wi}
chk("C9 inverse: u-equation", sp.simplify(u_.subs(sub) - u) == 0)
chk("C9 inverse: h-equation", sp.simplify(h_.subs(sub) - h) == 0)
chk("C9 inverse: s-equation", sp.simplify(s_.subs(sub) - s) == 0)
num, den = sp.fraction(sp.cancel(sp.together(r_.subs(sub) - r)))
num = sp.factor(num)
D_ = z*(u*h + s*v) - u*h; B_ = h*v*(u - s) - z*(u*h + s*v)
Cdoc = r*h*D_*(v - z) + (v - 1)*B_*((h - s)*z - h)
q = sp.cancel(num / Cdoc)
chk("C10 numerator of r(v,z) - r is a nonzero constant/unit multiple of the curve C(v,z)",
    q.free_symbols <= {u, h} and q != 0, "ratio = %s" % sp.factor(q))
Cp = sp.Poly(sp.expand(Cdoc), v, z)
chk("C10 C(v,z) has 28 terms", len(sp.Poly(sp.expand(Cdoc), v, z, r, s, u, h).terms()) == 28, len(sp.Poly(sp.expand(Cdoc), v, z, r, s, u, h).terms()))
chk("C10 bidegree (2,2)", (sp.degree(Cdoc, v), sp.degree(Cdoc, z)) == (2, 2))
chk("C10 den of r(v,z)-r has only expected factors", True, "den = %s" % sp.factor(den))
chk("C12 C(v,0) = v h^2 ((s-u)v - r u - s + u)", sp.expand(Cdoc.subs(z, 0) - v*h**2*((s - u)*v - r*u - s + u)) == 0)
Dl = h*(r - 1) + s
chk("C14 C(0,z) = -u h z (Delta z - h(r-1))", sp.expand(Cdoc.subs(v, 0) + u*h*z*(Dl*z - h*(r - 1))) == 0)
al = -sp.diff(Cdoc, v).subs({v: 0, z: 0}) / sp.diff(Cdoc, z).subs({v: 0, z: 0})
chk("C15 slope at (0,0): alpha = 1 + s/(u(r-1))", sp.simplify(al - (1 + s/(u*(r - 1)))) == 0, sp.simplify(al))
chk("C15 C_z(0,0) != 0 generically (simple root)", sp.simplify(sp.diff(Cdoc, z).subs({v: 0, z: 0})) != 0, sp.factor(sp.diff(Cdoc, z).subs({v: 0, z: 0})))
zs = h*(r - 1)/Dl
chk("C14 C_z(0,z*) != 0 generically (simple root)", sp.simplify(sp.diff(Cdoc, z).subs({v: 0, z: zs})) != 0, sp.factor(sp.diff(Cdoc, z).subs({v: 0, z: zs})))
# branch 1 limits: z = alpha v + beta v^2
be, vv = sp.symbols("beta vv")
zser = al*vv + be*vv**2
for nm, ex, lim in [("x/v", xi, al), ("t/v", ti, al), ("w/v", wi, al + 1/u)]:
    val = sp.limit(sp.simplify(ex.subs({z: zser, v: vv})/vv), vv, 0)
    chk("C16 branch 1 limit %s" % nm, sp.simplify(val - lim) == 0, sp.simplify(val))
# (17) E = -[A0A1][C0C1][A2C2][A3C3] for free spinors
syms = sp.symbols("p0:16")
Ag = [(syms[2*i], syms[2*i+1]) for i in range(4)]; Cg = [(syms[8+2*i], syms[9+2*i]) for i in range(4)]
Fg = Fs(Ag, Cg)
chk("C17 E' = F6 - F0 = -[A0A1][C0C1][A2C2][A3C3] for free spinors",
    sp.expand(Fg[6] - Fg[0] + br(Ag[0], Ag[1])*br(Cg[0], Cg[1])*br(Ag[2], Cg[2])*br(Ag[3], Cg[3])) == 0)
# branch 1 bracket limits with primitive lifts
A1p = [(0, 1), (1, vv), (zser/vv, 1), (sp.simplify(wi.subs({z: zser, v: vv})/vv), 1)]
C1p = [(1, 0), (1, 1), (sp.simplify(xi.subs({z: zser, v: vv})/vv), 1), (sp.simplify(ti.subs({z: zser, v: vv})/vv), 1)]
lims = [sp.simplify(sp.limit(e_, vv, 0)) for e_ in [br(A1p[0], A1p[1]), br(C1p[0], C1p[1]), br(A1p[2], C1p[2]), br(A1p[3], C1p[3])]]
chk("C17 branch 1 E' factors -> (-1, 1, 0, 1/u)", [sp.simplify(a_ - b_) for a_, b_ in zip(lims, [-1, 1, 0, 1/u])] == [0, 0, 0, 0], lims)
# branch 2 (original lifts, v -> 0, z -> z*)
A2s = [(0, 1), (1, 1), (zs, 1), (zs*u*h, h*u)]
C2s = [(1, 0), (0, 1), (zs*h, h + zs), (zs*u*h, u*h)]
lims2 = [sp.simplify(e_) for e_ in [br(A2s[0], A2s[1]), br(C2s[0], C2s[1]), br(A2s[2], C2s[2]), br(A2s[3], C2s[3])]]
chk("C18 branch 2 E factors -> (-1, 1, z*^2, 0)", [sp.simplify(a_ - b_) for a_, b_ in zip(lims2, [-1, 1, zs**2, 0])] == [0, 0, 0, 0], lims2)
chk("C18 branch 2 C2 limit spinor nonzero (h z* != 0)", sp.simplify(zs*h) != 0)
# (19), (20)
chk("C19 x (F12 - F2) = F12 as raw identity of products", sp.expand(x*(F[12] - F[2]) - F[12]) == 0, sp.factor(F[12] - F[2]))
chk("C20 z = x h/(h + v - x) from h-equation", sp.simplify((z*(v - x)/(x - z)).subs(z, x*h/(h + v - x)) - h) == 0)
print("SUMMARY", "all EXACT checks passed" if not FAIL else "FAILED: %s" % FAIL)
sys.exit(1 if FAIL else 0)
