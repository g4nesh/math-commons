"""Second-pass replay for blind_inventory86_matching_section_uniform_20261023 (paper case T_87).

Own code only. Uses the shared audit library moc4lib (paper coordinates B_0..B_13,
kappa = e_13, cubic C of (A.2), facets f_j of Appendix D). Nothing from the supplement
is executed, imported or read. The document's identities (5), (7), (8) are stated in the
supplement's own matching coordinates; they are replayed here in a coordinate-free form
(see the report): an element of the linear span of the 14 partials of C whose restriction
to the affine leaf is a polynomial in the non-minor coefficient d alone.

Run: python -I -B this_file.py [--numeric]
Prints one line per check; exits 1 if any EXACT check fails.
"""
import os
import sys
import itertools
import random
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L  # noqa: E402
import sympy  # noqa: E402

FAIL = []


def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + ((" :: " + str(info)) if info else ""))
    if not ok:
        FAIL.append(name)


T = frozenset({4, 5, 7, 9, 12, 14})
GOOD = [0, 1, 2, 5, 12, 45, 59]          # the document's seven minor generators
BAD, BAD2 = 70, 35                       # its two non-minor generators
BASIS8 = GOOD + [BAD]                    # order (a, b, c, e, f, g, k, d)
NAMES = ["a", "b", "c", "e", "f", "g", "k", "d"]

# ---------------------------------------------------------------- A. statement / case identity
check("T_is_paper_T87_literally", L.CONTACT_SETS[87] == T)
gens = L.generators(T)
check("generators(T)=={0,1,2,5,12,35,45,59,70}", sorted(gens) == sorted(GOOD + [BAD, BAD2]), gens)
check("minor_flags", all(L.F_IS_MINOR[j] for j in GOOD) and not L.F_IS_MINOR[BAD]
      and not L.F_IS_MINOR[BAD2])
PT = [L.P[s] for s in sorted(T)]
check("rank_P[T,:]==6 (kernel dim 8)", L.rank(PT) == 6)
check("8 profiles independent", L.rank([L.F_VECS[j] for j in BASIS8]) == 8)
check("profiles in kernel of P[T,:]",
      all(L.F_VALUES[j][s] == 0 for j in BASIS8 + [BAD2] for s in T))

# circuit (3)
lhs = [L.F_VECS[BAD2][i] + L.F_VECS[BAD][i] for i in range(14)]
rhs = [sum(L.F_VECS[j][i] for j in [1, 2, 5, 12, 45, 59]) for i in range(14)]
check("circuit (3) f35+f70=f1+f2+f5+f12+f45+f59", lhs == rhs)

# coefficients = exclusive gaps (18,16,13,10,8,3,2): at perm p, among the 8 profiles only
# the matching minor is nonzero, with value 1
excl = [18, 16, 13, 10, 8, 3, 2]
ok = True
for j, p in zip(GOOD, excl):
    row = [L.F_VALUES[q][p] for q in BASIS8]
    want = [1 if q == j else 0 for q in BASIS8]
    ok &= (row == want)
check("good coefficients are the gaps at (18,16,13,10,8,3,2)", ok)

# d<0: gaps (17,15,11,6,1,0) give b,c,e,f,g,k >= delta: row = e_x + e_d exactly
neg = [17, 15, 11, 6, 1, 0]
ok = True
rows = {}
for j, p in zip([1, 2, 5, 12, 45, 59], neg):
    row = [L.F_VALUES[q][p] for q in BASIS8]
    rows[p] = row
    want = [1 if q in (j, BAD) else 0 for q in BASIS8]
    ok &= (row == want)
check("d<0 rows: gap_p = coeff + d at p in (17,15,11,6,1,0)", ok)

# ---------------------------------------------------------------- stabilizer
stab = [g for g in L.G_ELEMENTS if L.g_act_set(g, T) == T]
acts = {}
for g in stab:
    acts.setdefault(L.g_index_map(g), g)
check("stabilizer has 12 distinct actions on the 24 indices", len(acts) == 12, len(acts))


def fimage(g, j):
    img = L.g_act_c(g, L.F_VECS[j])
    for q in range(82):
        if L.F_VECS[q] == img:
            return q
    return None


swap = [g for g in acts.values() if fimage(g, BAD2) == BAD and fimage(g, BAD) == BAD2
        and sorted(fimage(g, j) for j in GOOD) == sorted(GOOD)]
check("some stabilizer element swaps f35<->f70 and permutes the 7 minors", len(swap) > 0,
      "count=%d" % len(swap))
pres = [g for g in acts.values() if fimage(g, BAD) == BAD]
check("6 stabilizer actions preserve f70", len(pres) == 6)
check("f70-preserving actions permute the 7 minors",
      all(sorted(fimage(g, j) for j in GOOD) == sorted(GOOD) for g in pres))
edges = {frozenset(L.g_act_set(g, {4, 9})) for g in pres}
listed = {frozenset(e) for e in [(4, 9), (4, 14), (5, 7), (5, 12), (7, 14), (9, 12)]}
check("orbit of edge {4,9} = the six listed edges", edges == listed)
# connectivity of the edge graph on T
adj = {v: set() for v in T}
for e in listed:
    x, y = tuple(e)
    adj[x].add(y); adj[y].add(x)
seen, st = {4}, [4]
while st:
    v = st.pop()
    for w in adj[v]:
        if w not in seen:
            seen.add(w); st.append(w)
check("edge graph connected on the six contacts", seen == set(T))
# ordered version: for every listed ordered edge (p,q) some f70-preserving g maps (4,9)->(p,q)
# or (q,p); since H is a group, g^{-1} maps the edge back to (4,9).
check("f70-preserving actions form a group",
      all(tuple(L.g_index_map(g)[t] for t in L.g_index_map(h)) in
          {L.g_index_map(k) for k in pres} for g in pres for h in pres))

# ---------------------------------------------------------------- extra-tie table
table = {2: 68, 10: 68, 13: 68, 3: 51, 8: 51, 16: 51, 18: 61}
lit = {68: {1, 4, 5, 7, 9, 12}, 51: {1, 3, 5, 10, 11, 12}, 61: {1, 2, 3, 9, 11, 12}}
check("listed tuples are literally paper T_68, T_51, T_61",
      all(L.CONTACT_SETS[k] == frozenset(v) for k, v in lit.items()))
posrow = [0, 1, 6, 11, 15, 17, 19, 20, 21, 22, 23]
check("table covers all 18 noncontacts exactly once",
      sorted(list(table) + posrow) == sorted(set(range(24)) - T))
ok = True
wit = {}
for j, k in table.items():
    tgt = T | {j}
    found = None
    for g in L.G_ELEMENTS:
        if L.g_act_set(g, L.CONTACT_SETS[k]) <= tgt:
            found = g; break
    wit[j] = found
    ok &= found is not None
check("each extra tie j: some g maps T_k into T u {j}", ok,
      {j: (L.perm_to_oneline(g[0]), L.perm_to_oneline(g[1]), g[2]) for j, g in wit.items()})
check("positive source row: f70(j) > 0 for the 11 remaining j",
      all(L.F_VALUES[BAD][j] > 0 for j in posrow))
# alternative resolution (info): which paper cases are contained in T u {j} at all
info = {}
for j in sorted(set(range(24)) - T):
    tgt = T | {j}
    hits = []
    for k, S in L.CONTACT_SETS.items():
        if k == 87:
            continue
        if any(L.g_act_set(t, S) <= tgt for t in range(1152)):
            hits.append(k)
    info[j] = hits
print("INFO six-sets contained in T u {j} up to G:", info)
print("INFO aux7 equal to T u {j} up to G:",
      {j: [k for k, S in L.AUX7.items() if L.orbit_rep(S) == L.orbit_rep(T | {j})]
       for j in sorted(set(range(24)) - T)})

# ---------------------------------------------------------------- the affine leaf (6)
# variables: w = (a,b,c,e,f,g,k,d), z (kappa), and the complex differences
# (E, r, s, u, h) = (H9, H5, H7, H12, H14) * E with H4 = 0. E kept as a homogenizing parameter.
a, b, c, e, f, g_, k, d, z = sympy.symbols("a b c e f g k d z")
E, r, s, u, h = sympy.symbols("E r s u h")
W = [a, b, c, e, f, g_, k, d]
# particular solution v0: the unique vector with values (4:0, 5:r, 7:s, 9:E, 12:u, 14:h) on T
# that is orthogonal (standard pairing on V = Q^14) to the 8 profiles; it is rational and
# linear in (E,r,s,u,h), hence real for real parameters. Any other real choice differs by a
# real element of span(profiles), i.e. shifts w by a real vector, which leaves Im w unchanged.
target = {4: 0, 5: r, 7: s, 9: E, 12: u, 14: h}
Ms = sympy.Matrix([[sympy.Rational(x) for x in L.P[sidx]] for sidx in sorted(T)] +
                  [[sympy.Rational(x) for x in L.F_VECS[j]] for j in BASIS8])
rhsv = sympy.Matrix([target[sidx] for sidx in sorted(T)] + [0] * 8)
check("v0 system (6 T-rows + 8 orthogonality rows) is invertible", Ms.rank() == 14)
v0 = [sympy.expand(x) for x in Ms.LUsolve(rhsv)]
y = [sympy.expand(sum(W[t] * L.F_VECS[BASIS8[t]][i] for t in range(8)) + z * L.KAPPA[i] + v0[i])
     for i in range(14)]
# leaf check: values on T
vals = [sympy.expand(sum(L.P[sidx][i] * y[i] for i in range(14))) for sidx in range(24)]
check("leaf: H4=0? no -- F4 = z, F5-F4=r, F7-F4=s, F9-F4=E, F12-F4=u, F14-F4=h",
      sympy.expand(vals[5] - vals[4] - r) == 0 and sympy.expand(vals[7] - vals[4] - s) == 0
      and sympy.expand(vals[9] - vals[4] - E) == 0 and sympy.expand(vals[12] - vals[4] - u) == 0
      and sympy.expand(vals[14] - vals[4] - h) == 0 and sympy.expand(vals[4] - z) == 0)
# real coefficients: c(x) = Re(i eta (y - F4 kappa)); its decomposition in the 8 profiles is
# -eta*Im(w) + Re(i eta)*... ; with real (E,r,s,u,h) v0 is real, so c = sum_t (-eta Im w_t) f_t.
check("v0 is a rational-linear function of (E,r,s,u,h) (real for real parameters)",
      all(sympy.Poly(x, E, r, s, u, h).total_degree() <= 1 and
          all(co.is_rational for co in sympy.Poly(x, E, r, s, u, h).coeffs()) for x in v0 if x != 0))

ms = sympy.symbols("m0:14")
Cexpr = sympy.expand(L.cubic_C(ms))
partials = [sympy.diff(Cexpr, mi) for mi in ms]
sub = dict(zip(ms, y))
G = [sympy.Poly(sympy.expand(p.subs(sub)), *W, z, E, r, s, u, h) for p in partials]
others = [a, b, c, e, f, g_, k, z]
allmons = sorted({mon for p in G for mon in p.monoms()})
idx_other = [W.index(x) if x in W else 8 for x in others]
bad_mons = [mon for mon in allmons if any(mon[i] > 0 for i in idx_other)]
A = sympy.Matrix([[p.as_dict().get(mon, 0) for p in G] for mon in bad_mons])
ns = A.nullspace()
check("span of partials restricted to leaf contains a polynomial in (d;E,r,s,u,h) only",
      len(ns) >= 1, "dim=%d" % len(ns))
quads = []
for vec in ns:
    vec = vec * sympy.lcm([sympy.fraction(x)[1] for x in vec])
    pol = sympy.expand(sum(vec[i] * G[i].as_expr() for i in range(14)))
    quads.append((list(vec), pol))
for vec, pol in quads:
    print("INFO combination", vec, "->", sympy.factor(pol))
ok_monic = any(sympy.Poly(pol, d).degree() == 2 and sympy.Poly(pol, d).LC().is_number
               and sympy.Poly(pol, d).LC() != 0 for _, pol in quads)
check("(8) coordinate-free: some combination is quadratic in d with constant nonzero d^2 coefficient",
      ok_monic)
# (5): full-complex-equal stratum is E=r=s=u=h=0; the polynomial reduces to const*d^2
ok5 = any(sympy.expand(pol.subs({E: 0, r: 0, s: 0, u: 0, h: 0})) != 0 and
          sympy.Poly(sympy.expand(pol.subs({E: 0, r: 0, s: 0, u: 0, h: 0})), d).as_expr()
          == sympy.Poly(pol, d).LC() * d ** 2 for _, pol in quads)
check("(5) coordinate-free: on E=r=s=u=h=0 the same combination is const*d^2", ok5)

# (7): z appears at most linearly in each partial; find partials whose z-coefficient has a
# w-part with all coefficients of one sign (so Im is nonzero at strict points).
zcoef = [sympy.Poly(sympy.expand(sympy.diff(p.as_expr(), z)), *W, E, r, s, u, h) for p in G]
check("each partial has z-degree <= 1",
      all(sympy.Poly(p.as_expr(), z).degree() <= 1 for p in G))
good7 = []
for i, zc in enumerate(zcoef):
    lin = [zc.as_expr().coeff(x) for x in W]
    if zc.total_degree() <= 1 and any(cf != 0 for cf in lin) and (
            all(cf >= 0 for cf in lin) or all(cf <= 0 for cf in lin)):
        good7.append((i, sympy.expand(zc.as_expr())))
check("(7) coordinate-free: some partial is linear in z with a one-signed w-part of its z-coefficient",
      len(good7) > 0, good7)


# ---------------------------------------------------------------- exact test on actual matching vectors
def decompose(yv):
    """solve yv = sum w_t f_t + z kappa + v0(E,r,s,u,h) for exact complex data (QI)."""
    vv = [sum((L.P[sidx][i] * yv[i] for i in range(14)), L.QI(0)) for sidx in range(24)]
    F4 = vv[4]
    par = {E: vv[9] - F4, r: vv[5] - F4, s: vv[7] - F4, u: vv[12] - F4, h: vv[14] - F4}
    return vv, par


rng = random.Random(87)
ok = True
ntest = 0
qpoly = [pol for _, pol in quads if sympy.Poly(pol, d).degree() == 2][0]
for trial in range(6):
    A_ = [L.QI(Fraction(rng.randint(-9, 9), rng.randint(1, 5)), Fraction(rng.randint(-9, 9), rng.randint(1, 5))) for _ in range(4)]
    B_ = [L.QI(Fraction(rng.randint(-9, 9), rng.randint(1, 5)), Fraction(rng.randint(-9, 9), rng.randint(1, 5))) for _ in range(4)]
    m = L.matching_vector(A_, B_)
    vv, par = decompose(m)
    # w and z by solving sympy linear system with Gaussian rationals
    to_s = lambda q: sympy.Rational(q.re) + sympy.I * sympy.Rational(q.im) if isinstance(q, L.QI) else sympy.Rational(q)
    parS = {kk: to_s(vv_) for kk, vv_ in par.items()}
    eqs = [sympy.expand(y[i].subs(parS) - to_s(m[i])) for i in range(14)]
    sl = sympy.solve(eqs, W + [z], dict=True)
    if len(sl) != 1:
        ok = False; continue
    sl = sl[0]
    val = sympy.expand(qpoly.subs(parS).subs(sl))
    ok &= (sympy.simplify(val) == 0)
    ntest += 1
check("quadratic relation vanishes on %d exact Gaussian-rational matching vectors" % ntest,
      ok and ntest == 6)

# discriminant of the d-quadratic is invariant under the real shift of d caused by a
# different choice of v0: compare with the document's displayed quadratic (8) at E = 1
qd = sympy.Poly(qpoly.subs(E, 1), d)
A2, B1, C0 = qd.all_coeffs()
disc_mine = sympy.expand((B1 ** 2 - 4 * A2 * C0) / A2 ** 2)
doc_q = -d ** 2 + (h + r + s - u - 1) * d + r * (u - h)
dq = sympy.Poly(doc_q, d)
A2d, B1d, C0d = dq.all_coeffs()
disc_doc = sympy.expand((B1d ** 2 - 4 * A2d * C0d) / A2d ** 2)
check("(8): discriminant of own d-quadratic equals that of the displayed quadratic",
      sympy.expand(disc_mine - disc_doc) == 0)
shift = sympy.expand(B1 / (2 * A2) - B1d / (2 * A2d))
check("(8): own quadratic = displayed quadratic after a real affine shift of d",
      sympy.expand(qd.as_expr() / A2 - doc_q.subs(d, d + shift) / A2d) == 0, "shift=%s" % shift)

# optional literal replay of (8) and (7) in the document's own leaf data (JSON read as data)
supp = os.environ.get("MOC4_SUPP") or next((x.split("=", 1)[1] for x in sys.argv if x.startswith("--supp=")), None)
if supp:
    import json
    jp = os.path.join(supp, "proofs", "manin_transfer_20261021",
                      "blind_inventory86_matching_section_exact_20261023.json")
    J = json.load(open(jp))
    lift = J["lift14x8"]
    check("data: document's lift14x8 columns equal the paper-coordinate profiles f0..f70",
          all([Fraction(lift[i][t]) for i in range(14)] == L.F_VECS[BASIS8[t]] for t in range(8)))
    check("data: document's height direction is kappa", [Fraction(x) for x in J["height_direction"]] == L.KAPPA)
    loc = {str(x): x for x in (r, s, u, h)}
    v0d = [sympy.sympify(x, locals=loc) for x in J["real_offset14"]]
    yd = [sympy.expand(sum(W[t] * lift[i][t] for t in range(8)) + z * L.KAPPA[i] + v0d[i]) for i in range(14)]
    vd = [sympy.expand(sum(L.P[sidx][i] * yd[i] for i in range(14))) for sidx in range(24)]
    check("data: document's v0 satisfies the leaf equations (F9-F4=1, r, s, u, h)",
          all(sympy.expand(vd[p] - vd[4] - q) == 0 for p, q in [(9, 1), (5, r), (7, s), (12, u), (14, h)]))
    subd = dict(zip(ms, yd))
    G13m2 = sympy.expand(partials[13].subs(subd) - partials[2].subs(subd))
    check("(8) literal: G13 - G2 on the document's leaf = -d^2 + (h+r+s-u-1)d + r(u-h)",
          sympy.expand(G13m2 - doc_q) == 0)
    G1 = sympy.expand(partials[1].subs(subd))
    check("(7) literal: z-coefficient of G1 = -d-e-g+h-u",
          sympy.expand(sympy.diff(G1, z) - (-d - e - g_ + h - u)) == 0)
    # (5) literal: on the homogeneous full-equal stratum y = L8 w + z kappa
    y0 = [sympy.expand(sum(W[t] * lift[i][t] for t in range(8)) + z * L.KAPPA[i]) for i in range(14)]
    sub0 = dict(zip(ms, y0))
    check("(5) literal: (d13 C - d2 C)(L8 w + z kappa) = -d^2",
          sympy.expand(partials[13].subs(sub0) - partials[2].subs(sub0) + d ** 2) == 0)
else:
    print("SKIP literal data checks (set MOC4_SUPP to the supplement root)")

# ---------------------------------------------------------------- numeric falsification
if "--numeric" in sys.argv:
    import numpy as np
    worst = np.inf
    found = 0
    tested = 0
    for seed in range(int(os.environ.get("NSEEDS", "12"))):
        pt = L.random_point_on_contact_stratum(T, seed=seed, tries=int(os.environ.get("TRIES", "40")), margin=float(os.environ.get("MARGIN", "0")))
        if pt is None:
            continue
        found += 1
        mvec = L.matching_vector_np(pt["a"], pt["b"], pt["alpha"])
        cvec = mvec.real.copy()
        Fps = (L.P_numpy() @ mvec)[pt["pi_star"]]
        cvec[13] -= Fps.real
        vals_ = L.P_numpy() @ cvec
        # decompose into the 8 profiles (least squares; exact up to fp)
        Bm = np.array([[float(x) for x in L.F_VECS[j]] for j in BASIS8]).T
        coef, *_ = np.linalg.lstsq(Bm, cvec, rcond=None)
        val, O = L.min_phi_over_SO4([float(x) for x in cvec], restarts=12, seed=seed)
        tested += 1
        worst = min(worst, val / max(1e-12, np.max(np.abs(vals_))))
        print("NUMERIC seed=%d contact=%s d=%.4g minPhi=%.6g (normalized %.3g)"
              % (seed, pt["contact"], coef[7], val, val / np.max(np.abs(vals_))))
    print("NUMERIC margin=%s" % os.environ.get("MARGIN", "0"))
    print("NUMERIC summary: stratum points found=%d tested=%d worst normalized min Phi=%.3g"
          % (found, tested, worst))

print("SUMMARY: %d EXACT check(s) failed" % len(FAIL) if FAIL else "SUMMARY: all EXACT checks passed")
sys.exit(1 if FAIL else 0)
