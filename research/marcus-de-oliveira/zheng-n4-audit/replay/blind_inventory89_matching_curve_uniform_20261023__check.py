"""Own-code replay for the Zheng n=4 audit, document
blind_inventory89_matching_curve_uniform_20261023 (primary proof of paper case T_90).

Everything is rebuilt from the paper's definitions through the shared library moc4lib
(matching basis, P, cubic C of (A.2), generators f_j of Appendix D, group G).
Nothing from the supplement is imported, executed or read.
Exact checks always run; floating-point falsification runs with --numeric.
Prints one line per check; exit code 1 if an EXACT check fails."""
import sys, os, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L
import sympy as sp

NUMERIC = "--numeric" in sys.argv[1:]
FAIL = []


def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + "EXACT " + name + ((" : " + str(info)) if info != "" else ""))
    if not ok:
        FAIL.append(name)


T = [1, 5, 6, 9, 14, 18]
Tset = frozenset(T)
OFF = [j for j in range(24) if j not in Tset]
B8 = [1, 2, 5, 9, 29, 33, 59, 35]          # (f1,f2,f5,f9,f29,f33,f59,f35) <-> (a,b,c,e,f,g,k,d)
FV = L.F_VALUES

# ---------------- Section 2: contact set, generators, source basis ----------------
check("T equals paper T_90 (Appendix E)", L.CONTACT_SETS[90] == Tset, sorted(L.CONTACT_SETS[90]))
gens = L.generators(T)
check("generators(T) = {1,2,5,9,28,29,33,35,59}", gens == [1, 2, 5, 9, 28, 29, 33, 35, 59], gens)
check("non-minor generators are exactly f28, f35",
      [j for j in gens if not L.F_IS_MINOR[j]] == [28, 35], [L.f_name(j) for j in gens])
kdim = 14 - L.rank([L.P[i] for i in T])
check("dim ker P[T,:] = 8", kdim == 8, kdim)
check("the 8 basis vectors are independent (rank 8) and vanish on T",
      L.rank([L.F_VECS[j] for j in B8]) == 8 and all(FV[j][i] == 0 for j in B8 for i in T))
ex_ok = True
for idx, j in zip((17, 15, 11, 7, 4, 12, 0), B8[:7]):
    row = [FV[q][idx] for q in B8]
    want = [1 if q == j else 0 for q in B8]
    ex_ok &= (row == want)
check("exclusive gaps: coefficient of f1,f2,f5,f9,f29,f33,f59 = gap at 17,15,11,7,4,12,0", ex_ok)
lhs = [FV[28][i] + FV[35][i] for i in range(24)]
rhs = [2 * FV[1][i] + FV[2][i] + FV[5][i] + FV[29][i] + FV[33][i] for i in range(24)]
check("circuit (3): f28 + f35 = 2f1 + f2 + f5 + f29 + f33", lhs == rhs)
# d < 0 bounds: gap rows 22,19,20,8,3 in the basis
want = {22: {1: 1, 35: 2}, 19: {2: 1, 35: 1}, 20: {5: 1, 35: 1}, 8: {29: 1, 35: 1}, 3: {33: 1, 35: 1}}
ok = all({q: FV[q][i] for q in B8 if FV[q][i]} == w for i, w in want.items())
check("d<0: gaps 22,19,20,8,3 are a+2d, b+d, c+d, f+d, g+d", ok)
# with gaps>=0 the rewritten source is a nonnegative combination: implied by the two checks above.

# cone N_T equals cone of the 9 generators (Lemma 4.2 in the paper); here a direct check that
# every c in ker P[T] with P c >= 0 has nonnegative minor coefficients (exclusive gaps) -- done.

# ---------------- transposition and stabilizer ----------------
tmap = L.g_index_map((L.IDENTITY, L.IDENTITY, 1))
check("transposition sigma -> sigma^{-1} maps T to T", frozenset(tmap[j] for j in T) == Tset)


def act_on_f(t):
    """permutation of generator indices induced by index map t (values moved)."""
    out = {}
    for j in gens:
        img = [None] * 24
        for s in range(24):
            img[t[s]] = FV[j][s]
        out[j] = [q for q in range(82) if list(FV[q]) == img][0]
    return out


ft = act_on_f(tmap)
check("transposition swaps f28 <-> f35 and permutes the 7 minor generators",
      ft[28] == 35 and ft[35] == 28 and sorted(ft[q] for q in B8[:7]) == sorted(B8[:7]), ft)

stab = [g for g, t in zip(L.G_ELEMENTS, L._G_MAPS) if frozenset(t[j] for j in T) == Tset]
check("stabilizer of T in G has 16 elements", len(stab) == 16, len(stab))
keep35 = [g for g in stab if act_on_f(L.g_index_map(g))[35] == 35]
check("exactly 8 stabilizer elements fix f35 (hence f28, and permute the minors)", len(keep35) == 8
      and all(act_on_f(L.g_index_map(g))[28] == 28 for g in keep35), len(keep35))
edges = sorted(tuple(sorted((L.g_index_map(g)[1], L.g_index_map(g)[5]))) for g in keep35)
want_edges = sorted([(1, 5), (1, 9), (1, 14), (1, 18), (5, 6), (6, 9), (6, 14), (6, 18)])
check("images of the pair (1,5) under those 8 are the listed K_{2,4} edges", edges == want_edges, edges)
# connectivity of the listed graph on T
adj = {v: set() for v in T}
for x, y in want_edges:
    adj[x].add(y); adj[y].add(x)
seen, st = {1}, [1]
while st:
    v = st.pop()
    for w2 in adj[v]:
        if w2 not in seen:
            seen.add(w2); st.append(w2)
check("the edge graph is connected on all six contacts", seen == Tset)

# ---------------- boundary table ----------------
z35 = sorted(j for j in OFF if FV[35][j] == 0)
check("off-T indices where f35 vanishes are {0,4,7,11,12,15,17}", z35 == [0, 4, 7, 11, 12, 15, 17], z35)
pos = sorted(j for j in OFF if FV[35][j] > 0)
check("f35 > 0 exactly at {2,3,8,10,13,16,19,20,21,22,23} (d=0 rows)",
      pos == [2, 3, 8, 10, 13, 16, 19, 20, 21, 22, 23], pos)
T40 = L.CONTACT_SETS[40]
check("cited six-set {0,1,5,9,14,18} is paper T_40", T40 == frozenset({0, 1, 5, 9, 14, 18}))
ok = True
wit = {}
for j in (0, 4, 7, 11, 12, 15):
    U = Tset | {j}
    found = None
    for g, t in zip(L.G_ELEMENTS, L._G_MAPS):
        if frozenset(t[q] for q in T40) <= U:
            found = g; break
    wit[j] = (L.perm_to_oneline(found[0]), L.perm_to_oneline(found[1]), found[2]) if found else None
    ok &= found is not None
check("extra contact j in {0,4,7,11,12,15}: some g maps T_40 into T u {j}", ok, wit)
A24 = L.AUX7[24]
check("cited seven-set {0,3,7,11,12,16,20} is paper T^(7)_24", A24 == frozenset({0, 3, 7, 11, 12, 16, 20}))
g17 = L.find_g(A24, Tset | {17})
check("extra contact 17: T u {17} is in the G-orbit of T^(7)_24", g17 is not None,
      (L.perm_to_oneline(g17[0]), L.perm_to_oneline(g17[1]), g17[2]) if g17 else None)
# stage order (Table 6): T_40 and T^(7)_24 are Stage 1, T_90 is Stage 2.
print("INFO (read from the paper, not computed): Table 6 puts T_40 and T^(7)_24 in Stage 1 and T_90 in Stage 2")

# ---------------- Section 3: identities on the cubic ----------------
a, b, c, e, f, g, k, d, z, r, s, u, h = sp.symbols('a b c e f g k d z r s u h')
W = [a, b, c, e, f, g, k, d]
ms = sp.symbols('m0:14')
C = sp.expand(L.cubic_C(ms))
Gs = [sp.diff(C, ms[i]) for i in range(14)]
L8m = [sum(sp.Rational(L.F_VECS[jj][i]) * W[t] for t, jj in enumerate(B8)) + (z if i == 13 else 0)
       for i in range(14)]
e5 = sp.expand((Gs[12] - Gs[2]).subs(dict(zip(ms, L8m)), simultaneous=True))
check("(5): (d12 C - d2 C)(L8 w + z kappa) = -d(d+f+g)", sp.expand(e5 + d * (d + f + g)) == 0, sp.factor(e5))

v0 = [r - s + 1, -1, -h - r + s + u - 1, r - u + 1, -r] + [0] * 9
Pv0 = [sp.expand(sum(int(L.P[i][j]) * v0[j] for j in range(14))) for i in range(24)]
diffs = [sp.expand(Pv0[i] - Pv0[1]) for i in T]
check("(6): P v0 has differences (F_i - F_1)_{i in T} = (0,1,r,s,u,h)", diffs == [0, 1, r, s, u, h], diffs)
mm = [L8m[i] + v0[i] for i in range(14)]
sub = dict(zip(ms, mm))
G = [sp.expand(Gi.subs(sub, simultaneous=True)) for Gi in Gs]
lin = [G[4] - G[8], G[11] - G[3], G[10] - G[7], G[12] - G[2]]
vars4 = [b, c, f, g]
deg_ok = all(sp.Poly(q, *vars4).total_degree() <= 1 for q in lin) and \
    all(sp.Poly(q, e, k, z).total_degree() == 0 for q in lin)
check("four differences are linear in (b,c,f,g) and free of (e,k,z)", deg_ok)
J4 = sp.Matrix([[sp.diff(q, v) for v in vars4] for q in lin])
det7 = sp.expand(J4.det())
want7 = sp.expand(2 * (d + r) * (-d - h + s) * (-d + u - 1) * (d + h + r - s - u + 1))
check("(7): determinant of the 4x4 system", sp.expand(det7 - want7) == 0 or sp.expand(det7 + want7) == 0,
      "det = " + str(sp.factor(det7)) + ("  (equals (7))" if sp.expand(det7 - want7) == 0 else "  (equals -(7))"))
sol = sp.solve(lin, vars4, dict=True)
B0 = -h + r**2 - r * s - r * u + r + s * u
want8 = {b: (r * a - d**2 + (-h + 2 * r - s - u + 1) * d + B0) / (2 * (d + r)),
         c: -(-r * a - d**2 + (-h - s + u - 1) * d + B0) / (2 * (-d + u - 1)),
         f: -(-r * a + d**2 + (h - s - u + 1) * d + B0) / (2 * (d + h + r - s - u + 1)),
         g: -(-r * a - d**2 + (-h + s - u - 1) * d + B0) / (2 * (-d - h + s))}
ok8 = len(sol) == 1 and all(sp.simplify(sol[0][v] - want8[v]) == 0 for v in vars4)
check("(8): unique solution (b,c,f,g) as printed, with B0 = -h+r^2-rs-ru+r+su", ok8)
# (9)
q1, q2, q3 = G[13] - G[12], G[2] - G[3], G[6]
c1 = sp.expand(sp.diff(q1, e)); c2 = sp.expand(sp.diff(q2, k)); c3 = sp.expand(sp.diff(q3, z))
ok9 = (sp.expand(c1 - (-d - f - g)) == 0 and sp.expand(c2 - (-c - d - f)) == 0
       and sp.expand(c3 - (e + k + r - s + 1)) == 0)
tri = (sp.diff(q1, k) == 0 and sp.diff(q1, z) == 0 and sp.diff(q2, z) == 0
       and sp.Poly(q1, e).degree() == 1 and sp.Poly(q2, k).degree() == 1 and sp.Poly(q3, z).degree() == 1)
check("(9): coefficients of e in G13-G12, k in G2-G3, z in G6 are -d-f-g, -c-d-f, e+k+r-s+1", ok9,
      [c1, c2, c3])
check("(9): the three equations are triangular in (e,k,z) and each linear in its variable", tri,
      "dq1/dk=%s dq1/dz=%s dq2/dz=%s dq2/de=%s" % (sp.diff(q1, k), sp.diff(q1, z), sp.diff(q2, z),
                                                     sp.factor(sp.diff(q2, e))))
seven = lin + [q1, q2, q3]
vars7 = [b, c, f, g, e, k, z]
J7 = sp.Matrix([[sp.diff(q, v) for v in vars7] for q in seven]).det()
prod9 = (-d - f - g) * (-c - d - f) * (e + k + r - s + 1)
okJ = sp.expand(J7 - det7 * prod9) == 0
check("selected 7x7 Jacobian = det(4x4) * product of (9)", okJ)
# (10)
Af = h + r - s - u + 1; Ag = h - s
P2 = ((r - u - 1) * d**2 + (h * r - h * u - h - r * s - r * u - r + s * u + s + u**2 - 1) * d
      - h * r + h * u - h + r**3 - r**2 * s - 2 * r**2 * u + 2 * r**2 + 2 * r * s * u - r * s
      + r * u**2 - 2 * r * u + r - s * u**2 + s * u)
lhs10 = 2 * (d + Af) * (d + Ag) * (want8[f] + want8[g] + d)
check("(10): 2(d+Af)(d+Ag)(f+g+d) = -(r-u+1) r a + P2(d)",
      sp.simplify(lhs10 - (-(r - u + 1) * r * a + P2)) == 0)
L0 = -h + 2 * r - s - u + 1
check("b-equation of (8) is 2(d+r)b = ra - d^2 + L0 d + B0 with L0 as printed",
      sp.simplify(2 * (d + r) * want8[b] - (r * a - d**2 + L0 * d + B0)) == 0)
check("denominators of f and g are 2(d+Af) and -2(d+Ag)",
      sp.simplify(want8[f] * 2 * (d + Af) - (-(-r * a + d**2 + (h - s - u + 1) * d + B0))) == 0
      and sp.simplify(want8[g] * (-2) * (d + Ag) - (-(-r * a - d**2 + (-h + s - u - 1) * d + B0))) == 0)
# Lemma 3 imaginary-part identity
x, y, A_, al, X_, be = sp.symbols('x y A alpha X beta', real=True)
rr, L0s, B0s = sp.symbols('r L0 B0', real=True)
dd = x + sp.I * y; aa = A_ + sp.I * al; bb = X_ + sp.I * be
expr = sp.expand(2 * (dd + rr) * bb - (rr * aa - dd**2 + L0s * dd + B0s))
imI = sp.expand(sp.im(expr))
check("Lemma 3: Im of 2(d+r)b - (ra - d^2 + L0 d + B0) = 2[(x+r)beta + yX] - (r alpha - 2xy + L0 y)",
      sp.expand(imI - (2 * ((x + rr) * be + y * X_) - (rr * al - 2 * x * y + L0s * y))) == 0)
Xsol = sp.solve(sp.Eq(2 * ((x + rr) * be + y * X_), rr * al - 2 * x * y + L0s * y), X_)[0]
check("Lemma 3: X = -x - (x+r) beta/y + (r/2) alpha/y + L0/2",
      sp.simplify(Xsol - (-x - (x + rr) * be / y + rr * al / (2 * y) + L0s / 2)) == 0)
# ideal check: X subset {grad C = 0} on a random exact matching vector in the section family
import random
rng = random.Random(7)
av = [L.QI(rng.randint(-5, 5), rng.randint(-5, 5)) for _ in range(4)]
bv = [L.QI(rng.randint(-5, 5), rng.randint(-5, 5)) for _ in range(4)]
mv = L.matching_vector(av, bv)
Gv = [sp.nsimplify(0)]
okg = True
for i in range(14):
    val = sp.diff(C, ms[i]).subs({ms[j]: sp.Rational(mv[j].re) + sp.I * sp.Rational(mv[j].im) for j in range(14)})
    okg &= sp.expand(val) == 0
check("grad C = 0 at an exact Gaussian-rational matching vector (sanity of coordinates)", okg)

# ---------------- NUMERIC ----------------
# (a) Coordinate-alignment test on actual matching vectors of X that satisfy the five real
#     contact equations Re(F_i - F_1) = 0 (i in T) but not necessarily positivity: the
#     identities (6), (8), (10) hold on all of X, so they must hold at these points.
# (b) Brief attempt at the main statement on points with contact set exactly T.
if NUMERIC:
    import numpy as np
    from scipy.optimize import least_squares
    Pn = L.P_numpy()
    L8n = np.array([[float(L.F_VECS[jj][i]) for jj in B8] for i in range(14)])
    kap = np.zeros(14); kap[13] = 1.0
    M = np.hstack([L8n, kap[:, None]]).astype(complex)
    rng = np.random.default_rng(2026)
    found = 0; worst_res = 0.0; worst8 = 0.0; worst10 = 0.0; imE = 0.0; worstH = 0.0; worstgap = 0.0
    def unpack(x):
        return x[0:4] + 1j * x[4:8], x[8:12] + 1j * x[12:16], np.exp(1j * x[16])
    def res(x):
        A, Bv, al = unpack(x)
        Fv = L.F_values_np(A, Bv, al).real
        return np.array([Fv[i] - Fv[1] for i in T[1:]])
    for trial in range(30):
        x0 = rng.normal(size=17)
        sol = least_squares(res, x0, method="trf", xtol=1e-15, ftol=1e-15, gtol=1e-15)
        if np.max(np.abs(res(sol.x))) > 1e-11:
            continue
        found += 1
        A, Bv, al = unpack(sol.x)
        m = L.matching_vector_np(A, Bv, al)
        Fv = Pn @ m
        E = Fv[5] - Fv[1]
        imE = max(imE, abs(E.real) / abs(E))
        mp = m / E
        Fp = Pn @ mp
        H = [Fp[i] - Fp[1] for i in T]
        worstH = max(worstH, max(abs(Hh.imag) for Hh in H))
        rv, sv, uv, hv = H[2].real, H[3].real, H[4].real, H[5].real
        v0n = np.array([rv - sv + 1, -1, -hv - rv + sv + uv - 1, rv - uv + 1, -rv] + [0] * 9, complex)
        w_, _, _, _ = np.linalg.lstsq(M, mp - v0n, rcond=None)
        worst_res = max(worst_res, np.linalg.norm(M @ w_ - (mp - v0n)) / (1 + np.linalg.norm(mp)))
        aw, bw, cw, ew, fw, gw, kw, dw = w_[:8]
        subsn = {a: aw, d: dw, r: rv, s: sv, u: uv, h: hv}
        for vv, val in zip([b, c, f, g], [bw, cw, fw, gw]):
            worst8 = max(worst8, abs(complex(want8[vv].subs(subsn).evalf()) - val) / (1 + abs(val)))
        l10 = complex(lhs10.subs(subsn).evalf()); r10 = complex((-(r - u + 1) * r * a + P2).subs(subsn).evalf())
        worst10 = max(worst10, abs(l10 - r10) / (1 + abs(r10)))
        # real coefficient = -eta Im w with eta = Im E, compared with the real gaps at 17,15,11,7,4,12,0
        realcoef = -E.imag * np.imag(w_[:8])
        gaps = (Fv - Fv[1]).real
        worstgap = max(worstgap, max(abs(realcoef[q] - gaps[idx]) / (1 + abs(gaps[idx]))
                                     for q, idx in enumerate((17, 15, 11, 7, 4, 12, 0))))
    print("NUMERIC (a) points of X with the five real contact equations (seed 2026, 30 starts): %d" % found)
    print("NUMERIC (a) max |Re E|/|E|  (E must be purely imaginary): %.2e" % imE)
    print("NUMERIC (a) max |Im H_i|  (H_i must be real): %.2e" % worstH)
    print("NUMERIC (a) max rel. residual of the affine form (6): %.2e" % worst_res)
    print("NUMERIC (a) max rel. error of (8): %.2e" % worst8)
    print("NUMERIC (a) max rel. error of (10): %.2e" % worst10)
    print("NUMERIC (a) max rel. error of 'real coefficient = -eta Im w = gap at 17,15,11,7,4,12,0': %.2e" % worstgap)
    got = 0; worstphi = np.inf
    for seed in range(4):
        pt = L.random_point_on_contact_stratum(T, seed=seed, tries=10, margin=0.0, pi_star=1)
        if pt is None:
            continue
        got += 1
        cvec = np.linalg.lstsq(Pn, np.array(pt["gaps"], float), rcond=None)[0]
        val, _ = L.min_phi_over_SO4(cvec, restarts=6, seed=seed)
        worstphi = min(worstphi, val)
    print("NUMERIC (b) lib stratum points with contacts containing T (seeds 0..3, 10 tries, margin 0): %d; "
          "worst min_O Phi (6 BFGS restarts): %s" % (got, ("%.6e" % worstphi) if got else "n/a"))

print("SUMMARY: %d EXACT failures" % len(FAIL), FAIL)
sys.exit(1 if FAIL else 0)
