"""Independent verification of the three-plus-three generalized-circle theorem
(foundations/three-plus-three-circle.md).  Own code, written from the definitions; it
does not reuse six_circle_replay.py.  It imports only the audit library lib/moc4lib.py
(P, the minor vectors, the Appendix D table) and reads three supplement JSON files as
data (json.load), never executing anything from the supplement.

Usage:  python -I -B six_circle_verify.py STAGE [SUPP]
  STAGE in {atoms, dd, chambers, random, numeric, all}; SUPP defaults to $MOC4_SUPP.
Seeds are fixed; output is deterministic up to floating-point timing lines.
"""
import itertools
import json
import os
import random
import sys
import time
from fractions import Fraction as Fr
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L  # noqa: E402

PERMS, SGN = L.PERMS, L.SGN
FAILS = []


def report(ok, msg):
    print(("[PASS] " if ok else "[FAIL] ") + msg, flush=True)
    if not ok:
        FAILS.append(msg)


def supp_dir():
    if len(sys.argv) > 2:
        return sys.argv[2]
    p = os.environ.get("MOC4_SUPP")
    if not p:
        sys.exit("give the supplement path as argv[2] or MOC4_SUPP")
    return p


def load(name):
    with open(os.path.join(supp_dir(), "data", name)) as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------------
# Generators of K and the value-space form of a functional on V_R
# ---------------------------------------------------------------------------------
MINOR_VECS = [L.n_vec(I, J) for (I, J) in L.MINOR_KEYS]
PIV = L.P_PIVOT_ROWS
PPIV_INV = L.inv([L.P[r] for r in PIV])          # c = PPIV_INV * vals[PIV]
MASTER_VALS = [a + 2 * b + c for a, b, c in zip(L.F_VALUES[70], L.F_VALUES[2], L.F_VALUES[12])]


def atoms_vals():
    d = load("final_referee_six_circle_12atoms_partition_data_20261022.json")
    return [[Fr(x) for x in a] for a in d["fixed_atoms"]]


def to_value_functional(w):
    """w_hat in Q^24 (supported on pivot rows) with w.c = w_hat . (P c) for c in V_R."""
    coeff = [sum((Fr(w[k]) * PPIV_INV[k][j] for k in range(14)), Fr(0)) for j in range(14)]
    out = [Fr(0)] * 24
    for j, r in enumerate(PIV):
        out[r] = coeff[j]
    return out


def in_imP(vals):
    c = [sum((PPIV_INV[k][j] * vals[PIV[j]] for j in range(14)), Fr(0)) for k in range(14)]
    return all(sum((L.P[s][k] * c[k] for k in range(14)), Fr(0)) == vals[s] for s in range(24))


# ---------------------------------------------------------------------------------
# Exact double description (own implementation): extreme rays of {w : A w >= 0}
# ---------------------------------------------------------------------------------
def _prim(v):
    den = 1
    for x in v:
        den = den * Fr(x).denominator // gcd(den, Fr(x).denominator)
    iv = [int(Fr(x) * den) for x in v]
    g = 0
    for x in iv:
        g = gcd(g, abs(x))
    return tuple(x // g for x in iv) if g else tuple(iv)


def dd_extreme_rays(A):
    """A: list of rational rows (length d), rank d (pointed cone).  Returns primitive
    integer extreme rays of {w : a.w >= 0 for all rows a}.  Standard double description
    with the combinatorial adjacency test (Fukuda-Prodon)."""
    m, d = len(A), len(A[0])
    A = [_prim(a) for a in A]
    # pick d independent rows greedily
    basis, rows = [], []
    for i, a in enumerate(A):
        if L.rank(rows + [list(a)]) > len(rows):
            rows.append(list(a)); basis.append(i)
            if len(basis) == d:
                break
    assert len(basis) == d, "not full rank"
    Binv = L.inv([[Fr(x) for x in A[i]] for i in basis])
    rays = []   # (vector, zero-mask)
    for j in range(d):
        v = _prim([Binv[k][j] for k in range(d)])
        mask = 0
        for t, i in enumerate(basis):
            if t != j:
                mask |= 1 << i
        rays.append((v, mask))
    done = set(basis)
    for i in range(m):
        if i in done:
            continue
        a = A[i]
        vals = [sum(x * y for x, y in zip(a, r[0])) for r in rays]
        pos = [k for k, s in enumerate(vals) if s > 0]
        neg = [k for k, s in enumerate(vals) if s < 0]
        zer = [k for k, s in enumerate(vals) if s == 0]
        new = []
        masks = [r[1] for r in rays]
        for kp in pos:
            for kn in neg:
                Z = masks[kp] & masks[kn]
                if bin(Z).count("1") < d - 2:
                    continue
                adj = True
                for k2 in range(len(rays)):
                    if k2 != kp and k2 != kn and masks[k2] & Z == Z:
                        adj = False; break
                if not adj:
                    continue
                vp, vn = vals[kp], vals[kn]
                nv = _prim([vp * x - vn * y for x, y in zip(rays[kn][0], rays[kp][0])])
                new.append((nv, Z | (1 << i)))
        rays = ([(rays[k][0], rays[k][1] | (1 << i)) for k in zer] +
                [rays[k] for k in pos] + new)
        done.add(i)
    return [r[0] for r in rays]


def dual_vertices(gens):
    """Vertices of {w : w(kappa) = 1, w(g) >= 0 for g in gens} as Fraction 14-vectors."""
    rays = dd_extreme_rays(gens)
    out = []
    for r in rays:
        assert r[13] > 0, "kappa not interior?"
        out.append(tuple(Fr(x, r[13]) for x in r))
    return out


# ---------------------------------------------------------------------------------
# Stage atoms: quaternion identities, the inequality, the 12 atoms as images of m
# ---------------------------------------------------------------------------------
def stage_atoms():
    import sympy as sp
    p = sp.symbols("p0:4", real=True); q = sp.symbols("q0:4", real=True)

    def qmul(a, b):
        a0, a1, a2, a3 = a; b0, b1, b2, b3 = b
        return (a0*b0 - a1*b1 - a2*b2 - a3*b3, a0*b1 + a1*b0 + a2*b3 - a3*b2,
                a0*b2 - a1*b3 + a2*b0 + a3*b1, a0*b3 + a1*b2 - a2*b1 + a3*b0)
    qc = (q[0], -q[1], -q[2], -q[3])
    E = [tuple(1 if k == j else 0 for k in range(4)) for j in range(4)]
    cols = [qmul(qmul(p, e), qc) for e in E]
    O = sp.Matrix(4, 4, lambda i, j: sp.expand(cols[j][i]))
    N = sum(t**2 for t in p); M = sum(t**2 for t in q)
    report(sp.expand(O.T * O - N * M * sp.eye(4)) == sp.zeros(4, 4) and
           sp.expand(O.det() - N**2 * M**2) == 0,
           "O(p,q): x -> p x conj(q) satisfies O^T O = |p|^2|q|^2 I and det O = |p|^4|q|^4")

    def Phi_poly(vals, Om):
        tot = 0
        for s, sg, v in zip(PERMS, SGN, vals):
            if v:
                tot += sg * v * Om[0, s[0]] * Om[1, s[1]] * Om[2, s[2]] * Om[3, s[3]]
        return sp.expand(tot)
    names = (L.F_TABLE[2], L.F_TABLE[12], L.F_TABLE[70])
    report(names == (("n", "13", "34"), ("n", "14", "24"), ("h", "4231", 1)),
           "lib Appendix D: f2 = n_{13,34}, f12 = n_{14,24}, f70 = h_{4231,+}")
    # Hadamard coordinates P = H p, Q = H q
    Pv = sp.symbols("P0:4", real=True); Qv = sp.symbols("Q0:4", real=True)
    H = sp.Matrix([[1, 1, 1, 1], [1, -1, -1, 1], [1, 1, -1, -1], [1, -1, 1, -1]]) / 2
    report(sp.simplify(H * H.T - sp.eye(4)) == sp.zeros(4, 4), "H is orthogonal")
    pP = H.T * sp.Matrix(Pv); qQ = H.T * sp.Matrix(Qv)
    sub = {**{p[i]: pP[i] for i in range(4)}, **{q[i]: qQ[i] for i in range(4)}}
    Os = O.subs(sub, simultaneous=True).applyfunc(sp.expand)
    NP = sum(t**2 for t in Pv); NQ = sum(t**2 for t in Qv)
    x2 = sp.expand(Os.extract([0, 2], [2, 3]).det())
    x3 = sp.expand(Os.extract([0, 3], [1, 3]).det())
    v = [sp.expand(Pv[i]**2 * NQ - Qv[i]**2 * NP) for i in range(4)]
    dl = Pv[0]*Qv[2] + Pv[1]*Qv[3]
    DL = Pv[2]*Qv[0] + Pv[3]*Qv[1]
    Phi_f2 = Phi_poly(L.F_VALUES[2], Os); Phi_f12 = Phi_poly(L.F_VALUES[12], Os)
    Phi_f2f70 = Phi_poly([a + b for a, b in zip(L.F_VALUES[2], L.F_VALUES[70])], Os)
    Phi_m = Phi_poly(MASTER_VALS, Os)
    report(sp.expand(Phi_f2 - x2**2) == 0 and sp.expand(Phi_f12 - x3**2) == 0,
           "Phi_{f2} = x2^2, Phi_{f12} = x3^2 (x2 = det O[13,34], x3 = det O[14,24]), polynomial identities")
    report(sp.expand(x2 - (v[0] + v[1])) == 0 and sp.expand(x3 - (v[1] + v[3])) == 0
           and sp.expand(sum(v)) == 0,
           "x2 = v0+v1, x3 = v1+v3, v0+v1+v2+v3 = 0, with v_i = P_i^2|Q|^2 - Q_i^2|P|^2")
    ident = sp.expand(Phi_f2f70 - (2*(dl**2 + DL**2)*NP*NQ - 2*(v[0]*v[3] + v[1]*v[2])))
    report(ident == 0, "KEY: Phi(f2+f70) = 2(delta^2+Delta^2)|P|^2|Q|^2 - 2(v0 v3 + v1 v2),"
           " delta = P0Q2+P1Q3, Delta = P2Q0+P3Q1 (bihomogeneous polynomial identity)")
    report(sp.expand(x2*x3 - (v[0]*v[3] - v[1]*v[2])) == 0, "x2 x3 = v0 v3 - v1 v2")
    report(sp.expand(Phi_m - (Phi_f2f70 + x2**2 + x3**2)) == 0 and
           sp.expand(Phi_m - (2*(dl**2 + DL**2)*NP*NQ + (v[0]-v[3])**2 - 4*v[1]*v[2])) == 0,
           "Phi_m = Phi(f2+f70) + x2^2 + x3^2 = 2(delta^2+Delta^2) + (v0-v3)^2 - 4 v1 v2")
    # numeric sanity of (E1) and of the trigonometric lemma
    import numpy as np
    rng = np.random.default_rng(7)
    f = sp.lambdify(list(Pv) + list(Qv), [Phi_f2f70, x2, x3, Phi_m], "numpy")
    X = rng.normal(size=(200000, 8))
    X[:, :4] /= np.linalg.norm(X[:, :4], axis=1)[:, None]
    X[:, 4:] /= np.linalg.norm(X[:, 4:], axis=1)[:, None]
    a, b, c, dm = f(*X.T)
    e1 = a + 2 * np.abs(b * c)
    report(e1.min() > -1e-12 and dm.min() > -1e-12,
           f"NUMERIC: 200000 random (p,q): min of Phi(f2+f70)+2|x2 x3| = {e1.min():.3e}, min Phi_m = {dm.min():.3e}")
    ph = rng.uniform(-np.pi/2, np.pi/2, size=(2, 10**6))
    lhs = np.cos((ph[0] + ph[1]) / 2)**2 - np.minimum(np.cos(ph[0])**2, np.cos(ph[1])**2)
    report(lhs.min() > -1e-14, f"NUMERIC: lemma cos^2((a+b)/2) >= min(cos^2 a, cos^2 b) on [-pi/2,pi/2]^2: min {lhs.min():.2e}")

    # the 12 atoms: in im P, and each equals Phi_m(S1 O^(T) S2) with det S1 det S2 = +1
    Ov = sp.Matrix(4, 4, sp.symbols("o0:16"))
    Phim_sym = Phi_poly(MASTER_VALS, Ov)
    AT = atoms_vals()
    report(len(AT) == 12 and len({tuple(a) for a in AT}) == 12 and all(in_imP(a) for a in AT),
           "12 distinct stored atoms, each value vector in im P")
    allok = True
    for k, a in enumerate(AT):
        found = None
        for tr in (0, 1):
            for r in PERMS:
                for cc in PERMS:
                    # Phi_m(R O' C)(..) with (R O' C)_{ij} = O'_{r(i), cc(j)}
                    ok = True
                    for t_idx, tau in enumerate(PERMS):
                        # Phi_m(R O C) = sgn(r)sgn(cc) Phi_{m'}(O), m'_tau = m_{cc^{-1} tau r}
                        s = L.compose(L.inverse(cc), L.compose(tau, r))
                        if tr:
                            s = L.compose(L.inverse(cc), L.compose(L.inverse(tau), r))
                        if MASTER_VALS[L.perm_index(s)] != a[t_idx]:
                            ok = False; break
                    if ok:
                        found = (tr, r, cc); break
                if found:
                    break
            if found:
                break
        if not found:
            allok = False; print("   atom", k, "not found in orbit"); continue
        tr, r, cc = found
        Ob = Ov.T if tr else Ov
        Rm = sp.Matrix(4, 4, lambda i, j: Ob[r[i], cc[j]])
        sg = L.sgn(r) * L.sgn(cc)
        if sg < 0:                                   # fix det: flip first column sign
            Rm = Rm * sp.diag(-1, 1, 1, 1)
        lhs = Phi_poly(a, Ov)
        rhs = Phi_poly(MASTER_VALS, Rm)
        ok = sp.expand(lhs - rhs) == 0
        allok &= ok
    report(allok, "each atom c_k satisfies Phi_{c_k}(O) = Phi_m(S1 O^(T) S2) identically in 16 entries,"
           " with signed permutation matrices of det product +1 (so Phi_{c_k} >= 0 on SO(4))")


# ---------------------------------------------------------------------------------
# Stage dd
# ---------------------------------------------------------------------------------
_VCACHE = {}


def vertex_lists():
    if "R" in _VCACHE:
        return _VCACHE["R"], _VCACHE["Rp"]
    t0 = time.time()
    R = dual_vertices([list(v) for v in MINOR_VECS])
    AT = [L.solve_P(a) for a in atoms_vals()]
    Rp = dual_vertices([list(v) for v in MINOR_VECS] + [list(a) for a in AT])
    _VCACHE["R"], _VCACHE["Rp"] = R, Rp
    _VCACHE["t"] = time.time() - t0
    return R, Rp


def stage_dd():
    s = [sum((v[k] for v in MINOR_VECS), Fr(0)) for k in range(14)]
    report(s == [7 * x for x in L.KAPPA] and L.rank([list(v) for v in MINOR_VECS]) == 14,
           "sum of the 34 minor vectors = 7 kappa and they have rank 14 (kappa interior to K_0 and K)")
    R, Rp = vertex_lists()
    print(f"   own DD: |R| = {len(R)}, |R'| = {len(Rp)}  ({_VCACHE['t']:.1f}s)")
    rows = {tuple(Fr(x) for x in L.P[s]) for s in range(24)}
    report(rows <= set(R) and rows <= set(Rp), "the 24 evaluation functionals are vertices of R and R'")
    rel = load("final_referee_minor_relaxation_boundary_20261022.json")
    pos = load("final_referee_six_circle_positive_partition_data_20261022.json")
    alt = load("final_referee_six_circle_12atoms_partition_data_20261022.json")
    J72 = [tuple(Fr(x) for x in w) for w in pos["target_coordinates"]["minor_relaxation"]]
    J72b = [tuple(Fr(x) for x in w) for w in rel["facet_normals"]]
    J86 = [tuple(Fr(x) for x in w) for w in alt["vertices"]]
    report(set(J72) == set(R) and len(J72) == len(R) == 72, "own DD vertex set of R equals the stored 72-list")
    report(set(J72b) == set(R), "the 72 facet normals of the relaxation file are the same set")
    report(set(J86) == set(Rp) and len(J86) == len(Rp) == 86, "own DD vertex set of R' equals the stored 86-list")
    cut = [i for i, w in enumerate(J72) if w not in set(Rp)]
    print(f"   vertices of R removed by the atoms (indices in the 72-list): {cut}")
    # independent re-check of each vertex: feasibility and active rank 13
    AT = [L.solve_P(a) for a in atoms_vals()]
    gens = [list(v) for v in MINOR_VECS] + [list(a) for a in AT]
    ok = True
    for w in Rp:
        vals = [sum(x * y for x, y in zip(w, g)) for g in gens]
        if min(vals) < 0:
            ok = False
        act = [g for g, x in zip(gens, vals) if x == 0]
        if L.rank(act) != 13:
            ok = False
    report(ok, "every vertex of R' is feasible with an active set of rank 13 (re-check)")
    # NUMERIC cross-check of the DD implementation: random LPs over R and R'
    import numpy as np
    from scipy.optimize import linprog
    rng = np.random.default_rng(11)
    for name, G, V in [("R", [list(v) for v in MINOR_VECS], R), ("R'", gens, Rp)]:
        Gm = np.array([[float(x) for x in g] for g in G]); Vm = np.array([[float(x) for x in v] for v in V])
        worst = 0.0
        for _ in range(500):
            obj = rng.normal(size=14)
            res = linprog(obj, A_ub=-Gm, b_ub=np.zeros(len(G)), A_eq=np.eye(14)[13:14], b_eq=[1.0],
                          bounds=[(None, None)] * 14, method="highs")
            worst = max(worst, abs(res.fun - (Vm @ obj).min()))
        report(worst < 1e-7, f"NUMERIC: 500 random LPs over {name}: LP optimum = min over own vertex list (max gap {worst:.1e})")


# ---------------------------------------------------------------------------------
# Configurations: homogeneous Gaussian-rational points, F_sigma, profile space
# ---------------------------------------------------------------------------------
def cm(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def cs(a, b):
    return (a[0] - b[0], a[1] - b[1])


def br(u, v):   # [u, v] = u0 v1 - u1 v0 (lib convention, v_i = (z_i, 1))
    return cs(cm(u[0], v[1]), cm(u[1], v[0]))


def Fvec(A, C):
    B = [[br(A[i], C[j]) for j in range(4)] for i in range(4)]
    out = []
    for s in PERMS:
        t = B[0][s[0]]
        for i in range(1, 4):
            t = cm(t, B[i][s[i]])
        out.append(t)
    return out


def gap_vals(F, alpha):
    re = [cm(alpha, f)[0] for f in F]
    mn = min(re)
    return [x - mn for x in re]


ONE, ZERO = (Fr(1), Fr(0)), (Fr(0), Fr(0))
E_INF, E_ZERO = (ONE, ZERO), (ZERO, ONE)


def profile_basis(A, C, ia, jc):
    """Value vectors spanning {P c : c in L(Gamma)}: Re, Im of the four coefficient vectors
    (free points A[ia], C[jc] replaced by the homogeneous basis) and the all-ones vector."""
    vecs = []
    for za in (E_INF, E_ZERO):
        for wc in (E_INF, E_ZERO):
            A2 = list(A); C2 = list(C); A2[ia] = za; C2[jc] = wc
            F = Fvec(A2, C2)
            vecs.append([f[0] for f in F]); vecs.append([f[1] for f in F])
    vecs.append([Fr(1)] * 24)
    # exact basis by row reduction
    basis = []
    for v in vecs:
        if L.rank(basis + [v]) > len(basis):
            basis.append(v)
    return basis


def profile_rays(basis):
    """Extreme rays (as 24-value vectors) of {sum u_k basis_k : values >= 0}."""
    d = len(basis)
    rows = [[basis[k][s] for k in range(d)] for s in range(24)]
    rays = dd_extreme_rays(rows)
    return [[sum((Fr(u[k]) * basis[k][s] for k in range(d)), Fr(0)) for s in range(24)] for u in rays]


class Tester:
    def __init__(self):
        R, Rp = vertex_lists()
        self.WR = [to_value_functional(w) for w in R]
        self.WRp = [to_value_functional(w) for w in Rp]

    @staticmethod
    def minval(W, vals):
        return min(sum((a * b for a, b in zip(w, vals) if a), Fr(0)) for w in W)


def chamber_of(order_colors):
    """order_colors: cyclic list of 'A'/'C' (6 entries, 3 each).  Returns I, II or III."""
    n = len(order_colors)
    ia = [k for k in range(n) if order_colors[k] == "A"]
    gaps = []
    for t in range(3):
        a, b = ia[t], ia[(t + 1) % 3]
        gaps.append(((b - a) % n) - 1)
    g = sorted(gaps)
    return {(0, 0, 3): "I", (0, 1, 2): "II", (1, 1, 1): "III"}[tuple(g)]


def relabelings_III(order):
    """order: cyclic list of ('A', idx) / ('C', idx) in alternating pattern.  Returns the
    six relabelings (rowsA[0..2], colsC[0..2]) with C_k on the arc between A_k and A_{k+1}
    not containing the third A."""
    out = []
    pos = {lab: k for k, lab in enumerate(order)}
    As = [lab for lab in order if lab[0] == "A"]
    n = len(order)
    for perm in itertools.permutations(range(3)):
        a = [As[perm[0]], As[perm[1]], As[perm[2]]]
        cs_ = []
        for k in range(3):
            x, y, z = pos[a[k]], pos[a[(k + 1) % 3]], pos[a[(k + 2) % 3]]
            # the arc from x to y (either direction) not containing z
            fwd = [(x + t) % n for t in range(1, (y - x) % n)]
            bwd = [(y + t) % n for t in range(1, (x - y) % n)]
            arc = fwd if z not in fwd else bwd
            cc = [order[t] for t in arc if order[t][0] == "C"]
            assert len(cc) == 1
            cs_.append(cc[0])
        out.append(([l[1] for l in a], [l[1] for l in cs_]))
    return out


# ---------------------------------------------------------------------------------
# Stage chambers: replay of the stored certificates with own F_sigma; own samples
# ---------------------------------------------------------------------------------
def stage_chambers():
    import sympy as sp
    x, y, z, zeta, eta, s_, t_, u_ = sp.symbols("x y z zeta eta s t u")
    # representative frame: A = (inf, 0, 1, zeta), C = (s, t, u, eta), [A,C] = A0 C1 - A1 C0
    Ah = [(1, 0), (0, 1), (1, 1), (zeta, 1)]
    Ch = [(s_, 1), (t_, 1), (u_, 1), (eta, 1)]
    Fs = []
    for sig in PERMS:
        f = 1
        for i in range(4):
            a, c = Ah[i], Ch[sig[i]]
            f *= a[0] * c[1] - a[1] * c[0]
        Fs.append(sp.expand(f))
    cols = []
    for mon in [(0, 0), (1, 0), (0, 1), (1, 1)]:
        cols.append([sp.Poly(f, zeta, eta).coeff_monomial(zeta**mon[0] * eta**mon[1]) for f in Fs])
    report(all(sp.expand(Fs[k] - (-Ch[PERMS[k][1]][0]) * (1 - Ch[PERMS[k][2]][0]) * (zeta - Ch[PERMS[k][3]][0])) == 0
               for k in range(24)), "own F_sigma equals the document's -C_s(2)(1-C_s(3))(zeta-C_s(4))")
    # columns lie in im P identically: c_k = PPIV_INV vals[PIV], check all 24 rows
    ccols = []
    ok = True
    for col in cols:
        c = [sp.expand(sum(sp.Rational(PPIV_INV[k][j].numerator, PPIV_INV[k][j].denominator) * col[PIV[j]]
                           for j in range(14))) for k in range(14)]
        for srow in range(24):
            if sp.expand(sum(L.P[srow][k] * c[k] for k in range(14)) - col[srow]) != 0:
                ok = False
        ccols.append(c)
    report(ok, "the four coefficient columns are P-images of real vectors, identically in s,t,u")
    chart = {"I": {s_: x / (1 + x + y + z), t_: (x + y) / (1 + x + y + z), u_: (x + y + z) / (1 + x + y + z)},
             "II": {s_: x / (1 + x + y), t_: (x + y) / (1 + x + y), u_: 1 + z},
             "III": {s_: -x, t_: y / (1 + y), u_: 1 + z}}
    pos = load("final_referee_six_circle_positive_partition_data_20261022.json")
    alt = load("final_referee_six_circle_12atoms_partition_data_20261022.json")
    J72 = [tuple(Fr(v) for v in w) for w in pos["target_coordinates"]["minor_relaxation"]]
    J86 = [tuple(Fr(v) for v in w) for w in alt["vertices"]]
    R, Rp = vertex_lists()
    rows = {tuple(Fr(v) for v in L.P[s]) for s in range(24)}

    def poly_of(terms):
        tot = 0
        nonneg = True
        for mon, co in terms:
            q = sp.Rational(co)
            if q < 0:
                nonneg = False
            tot += q * x**mon[0] * y**mon[1] * z**mon[2]
        return sp.expand(tot), nonneg

    def check_record(rec, Vlist, cham, override_target=None):
        tgt = rec["target"] if override_target is None else override_target
        w = Vlist[tgt]
        Rp_, nn = poly_of(rec["R"])
        if Rp_ == 0 or not nn:
            return False, tgt
        Ws = {}
        for key, terms in rec["W"].items():
            pw, nn2 = poly_of(terms)
            if not nn2:
                return False, tgt
            Ws[int(key)] = pw
        if sp.expand(sum(Ws.values()) - Rp_) != 0:          # kappa column
            return False, tgt
        for k in range(4):
            lhs = sum(Ws[sg] * cols[k][sg] for sg in Ws)
            wk = sum(sp.Rational(w[j].numerator, w[j].denominator) * ccols[k][j] for j in range(14))
            diff = sp.together((lhs - Rp_ * wk).subs(chart[cham], simultaneous=True))
            if sp.expand(sp.numer(diff)) != 0:
                return False, tgt
        return True, tgt

    for cham, recs, Vlist, own in [("I", pos["cells"]["I"]["records"], J72, R),
                                   ("II", pos["cells"]["II"]["records"], J72, R),
                                   ("III", alt["records"], J86, Rp)]:
        t0 = time.time()
        good = True; tg = []
        for rec in recs:
            ok, tgt = check_record(rec, Vlist, cham)
            good &= ok; tg.append(tgt)
        nonord = {i for i, w in enumerate(Vlist) if w not in rows}
        cover = set(tg) == nonord and len(tg) == len(nonord)
        report(good and cover and set(Vlist) == set(own),
               f"chamber {cham}: {len(recs)} certificates replayed with own F_sigma (nonnegative coefficients,"
               f" R != 0, kappa + 4 column identities exact); targets = all {len(nonord)} non-ordinary vertices"
               f" of own DD list  ({time.time()-t0:.1f}s)")
        bad_t = next(i for i in sorted(nonord) if i != recs[0]["target"])
        ok, _ = check_record(recs[0], Vlist, cham, override_target=bad_t)
        report(not ok, f"chamber {cham}: negative control (record 0 against vertex {bad_t}) rejected")
    # failed list in the 12atoms file must be empty
    report(alt.get("failed", []) == [], "12-atom file lists no failed vertex")

    # own sample test: exact extreme rays of the profile cone at rational chamber points
    T = Tester()
    rng = random.Random(20261125)
    worst = {"I": None, "II": None, "III": None}
    nrays = 0; npts = 0; need_atom = 0
    def rq():
        return Fr(rng.choice([1, 1, 2, 3, 5, 7, 11, 13, 50, 200]), rng.choice([1, 2, 3, 7, 19, 100, 1000]))
    pts = []
    for cham in ("I", "II", "III"):
        for _ in range(25):
            pts.append((cham, rq(), rq(), rq()))
    for cham, X, Y, Z in pts:
        if cham == "I":
            D = 1 + X + Y + Z; S, Tt, U = X / D, (X + Y) / D, (X + Y + Z) / D
        elif cham == "II":
            D = 1 + X + Y; S, Tt, U = X / D, (X + Y) / D, 1 + Z
        else:
            S, Tt, U = -X, Y / (1 + Y), 1 + Z
        A = [E_INF, E_ZERO, (ONE, ONE), E_INF]
        C = [((S, Fr(0)), ONE), ((Tt, Fr(0)), ONE), ((U, Fr(0)), ONE), E_INF]
        basis = profile_basis(A, C, 3, 3)
        rays = profile_rays(basis)
        npts += 1; nrays += len(rays)
        W = T.WR if cham in ("I", "II") else T.WRp
        mv = min(T.minval(W, r) for r in rays)
        if cham == "III":
            need_atom += sum(1 for r in rays if T.minval(T.WR, r) < 0)
        worst[cham] = mv if worst[cham] is None else min(worst[cham], mv)
    report(all(v >= 0 for v in worst.values()),
           f"own exact sample: {npts} rational chamber points, {nrays} profile extreme rays, every ray in K_0"
           f" (I, II) resp. K (III); min w(ray) per chamber = {[str(worst[c]) for c in ('I','II','III')]};"
           f" chamber-III rays outside the minor cone: {need_atom}")


# ---------------------------------------------------------------------------------
# Stage random: exact test on random rational configurations (incl. degenerate)
# ---------------------------------------------------------------------------------
def stage_random(n_conf=3900, n_strong=3900):
    T = Tester()
    rng = random.Random(4242)

    def rq(a=9, den=(1, 2, 3, 4, 5, 7)):
        return Fr(rng.randint(-a, a), rng.choice(den))

    def gq():
        return (rq(), rq())

    def moebius():
        while True:
            m = [[(Fr(rng.randint(-3, 3)), Fr(rng.randint(-3, 3))) for _ in range(2)] for _ in range(2)]
            det = cs(cm(m[0][0], m[1][1]), cm(m[0][1], m[1][0]))
            if det != ZERO:
                return m

    def cadd(a, b):
        return (a[0] + b[0], a[1] + b[1])

    def apply(m, h):   # homogeneous Moebius action h -> m h
        return (cadd(cm(m[0][0], h[0]), cm(m[0][1], h[1])),
                cadd(cm(m[1][0], h[0]), cm(m[1][1], h[1])))

    kinds = ["generic", "AA", "CC", "AC", "inf", "free_eq", "free_real", "AC+inf", "AA+CC", "free_circle",
             "alt", "alt+AC", "alt+inf+free_eq"]
    stats = {k: [0, 0] for k in kinds}
    cham_count = {"I": 0, "II": 0, "III": 0}
    rank_hist = {}
    fails = []
    strong_rays = 0; strong_fail = 0; need_atom_c = 0; need_atom_rays = 0
    t0 = time.time()
    for it in range(n_conf):
        kind = kinds[it % len(kinds)]
        # six real parameters (theta) with labels; None = infinity
        th = []
        used = set()
        while len(th) < 6:
            v = rq(12, (1, 2, 3, 5))
            if v not in used:
                used.add(v); th.append(v)
        lab = [("A", 0), ("A", 1), ("A", 2), ("C", 0), ("C", 1), ("C", 2)]
        if kind.startswith("alt"):            # force the alternating pattern A C A C A C
            srt = sorted(th)
            th = [srt[0], srt[2], srt[4], srt[1], srt[3], srt[5]]
            if kind == "alt+AC":
                th[3] = th[0]                     # C0 = A0 (adjacent on the circle)
            if kind == "alt+inf+free_eq":
                th[5] = None
        if kind in ("AA", "AA+CC"):
            th[1] = th[0]
        if kind in ("CC", "AA+CC"):
            th[4] = th[3]
        if kind in ("AC", "AC+inf"):
            th[3] = th[rng.randrange(3)]
        if kind in ("inf", "AC+inf"):
            th[rng.choice([2, 5])] = None
        frame = moebius() if (kind not in ("inf", "AC+inf", "free_real", "alt+inf+free_eq") and it % 3 != 0) else None
        def h_of(t):
            hp = E_INF if t is None else ((t, Fr(0)), ONE)
            return apply(frame, hp) if frame else hp
        six = [h_of(t) for t in th]
        # free points
        if kind in ("free_eq", "alt+inf+free_eq"):
            fa = six[rng.randrange(6)]; fc = six[rng.randrange(6)]
        elif kind == "free_real":
            fa = ((rq(), Fr(0)), ONE); fc = ((rq(), Fr(0)), ONE)
        elif kind == "free_circle":
            fa = h_of(rq()); fc = h_of(rq())
        else:
            fa = (gq(), ONE) if rng.random() < 0.9 else E_INF
            fc = (gq(), ONE) if rng.random() < 0.9 else E_INF
        # random placement of the three A rows / C columns among 0..3
        ra = rng.sample(range(4), 4); rc = rng.sample(range(4), 4)
        A = [None] * 4; C = [None] * 4
        for k in range(3):
            A[ra[k]] = six[k]; C[rc[k]] = six[3 + k]
        A[ra[3]] = fa; C[rc[3]] = fc
        # cyclic order (ties broken by a random tiny perturbation)
        key = []
        for k in range(6):
            t = th[k]
            key.append((float("inf") if t is None else float(t), rng.random(), k))
        order = [lab[k] for _, _, k in sorted(key)]
        cham = chamber_of([o[0] for o in order])
        cham_count[cham] += 1
        F = Fvec(A, C)
        if all(f == ZERO for f in F):
            continue
        # prescribed relabelings
        if cham == "III":
            rels = relabelings_III(order)
        else:
            rels = [([0, 1, 2], [0, 1, 2])]
        W = T.WR if cham != "III" else T.WRp
        ok = True
        alphas = [gq() for _ in range(3)]
        alphas = [a for a in alphas if a != ZERO] or [ONE]
        for (ia, jc) in rels:
            A2 = [A[ra[i]] for i in ia] + [A[ra[3]]]
            C2 = [C[rc[j]] for j in jc] + [C[rc[3]]]
            F2 = Fvec(A2, C2)
            for al in alphas:
                vals = gap_vals(F2, al)
                if T.minval(W, vals) < 0:
                    ok = False
                    fails.append((kind, cham, th, A2, C2, al))
                if cham == "III" and T.minval(T.WR, vals) < 0:
                    need_atom_c += 1
        stats[kind][0] += 1; stats[kind][1] += (0 if ok else 1)
        if it < n_strong:
            ia, jc = rels[0]
            A2 = [A[ra[i]] for i in ia] + [A[ra[3]]]
            C2 = [C[rc[j]] for j in jc] + [C[rc[3]]]
            basis = profile_basis(A2, C2, 3, 3)
            rank_hist[len(basis)] = rank_hist.get(len(basis), 0) + 1
            if len(basis) > 5:
                fails.append(("rank>5", kind, th))
                strong_fail += 1
                continue
            rays = profile_rays(basis)
            for r in rays:
                strong_rays += 1
                if T.minval(W, r) < 0:
                    strong_fail += 1
                    fails.append(("ray", kind, cham, th))
                if cham == "III" and T.minval(T.WR, r) < 0:
                    need_atom_rays += 1
    tot = sum(v[0] for v in stats.values()); nf = sum(v[1] for v in stats.values())
    print("   per kind (configs, failures):", {k: tuple(v) for k, v in stats.items()})
    print("   chambers by perturbed cyclic order:", cham_count)
    report(nf == 0, f"EXACT: {tot} random rational concyclic 3+3 configurations (Gaussian-rational circle frames,"
           f" {tot - stats['generic'][0]} degenerate/special), up to 3 directions each and all prescribed"
           f" relabelings: gap vector in g.K (K_0 in chambers I/II, K in III); failures {nf}  ({time.time()-t0:.0f}s)")
    print(f"   chamber-III gap vectors outside the minor cone (atoms needed): {need_atom_c}")
    print("   profile-space rank histogram (strong test):", rank_hist)
    report(strong_fail == 0, f"EXACT strong form: profile cones of the first {n_strong} configurations,"
           f" {strong_rays} extreme rays, all in g.K; failures {strong_fail}; chamber-III rays needing atoms: {need_atom_rays}")
    for f in fails[:5]:
        print("   FAILURE:", f)


# ---------------------------------------------------------------------------------
# Stage numeric: floating-point attack on the statement
# ---------------------------------------------------------------------------------
def stage_numeric(n_conf=20000):
    import numpy as np
    from scipy.optimize import minimize
    import mpmath as mp
    rng = np.random.default_rng(20261126 + int(os.environ.get("SIXCIRCLE_BATCH", "0")))
    perms = np.array(PERMS); sg = np.array(SGN, dtype=float)

    def O_from(x):
        p = x[:4] / np.linalg.norm(x[:4]); q = x[4:] / np.linalg.norm(x[4:])
        p0, p1, p2, p3 = p; q0, q1, q2, q3 = q
        Lp = np.array([[p0, -p1, -p2, -p3], [p1, p0, -p3, p2], [p2, p3, p0, -p1], [p3, -p2, p1, p0]])
        Rq = np.array([[q0, q1, q2, q3], [-q1, q0, -q3, q2], [-q2, q3, q0, -q1], [-q3, -q2, q1, q0]])
        return Lp @ Rq

    def weights_batch(X):
        p = X[:, :4] / np.linalg.norm(X[:, :4], axis=1)[:, None]
        q = X[:, 4:] / np.linalg.norm(X[:, 4:], axis=1)[:, None]
        Os = np.stack([O_from(np.concatenate([a, b])) for a, b in zip(p, q)])
        Wt = np.ones((len(X), 24))
        for i in range(4):
            Wt *= Os[:, i, perms[:, i]]
        return Wt * sg

    # sanity: O_from is in SO(4)
    Ot = O_from(rng.normal(size=8))
    assert np.allclose(Ot @ Ot.T, np.eye(4)) and abs(np.linalg.det(Ot) - 1) < 1e-12
    X0 = rng.normal(size=(3000, 8))
    W0 = weights_batch(X0)

    def phi(x, vals):
        O = O_from(x)
        return float(np.dot(sg * vals, np.prod(O[np.arange(4)[None, :], perms], axis=1)))

    def Fnp(A, C):
        B = A[:, None] - C[None, :]
        return np.prod(B[np.arange(4)[None, :], perms], axis=1)

    worst = []
    t0 = time.time()
    n_local = 0
    for it in range(n_conf):
        th = rng.standard_cauchy(6)
        if it % 5 == 0:                         # near-coincidence
            i, j = rng.choice(6, 2, replace=False); th[j] = th[i] + rng.normal() * 1e-3
        m = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        pts = (m[0, 0] * th + m[0, 1]) / (m[1, 0] * th + m[1, 1])
        fa = complex(rng.normal(), rng.normal()) * rng.choice([1, 3]); fc = complex(rng.normal(), rng.normal()) * rng.choice([1, 3])
        if it % 7 == 0:
            fa = pts[rng.integers(6)]
        if it % 11 == 0:
            fa, fc = rng.normal(), rng.normal()
        ra = rng.permutation(4); rc = rng.permutation(4)
        A = np.zeros(4, complex); C = np.zeros(4, complex)
        A[ra[:3]] = pts[:3]; C[rc[:3]] = pts[3:]; A[ra[3]] = fa; C[rc[3]] = fc
        al = np.exp(2j * np.pi * rng.random())
        re = (al * Fnp(A, C)).real
        vals = re - re.min()
        if vals.max() <= 0:
            continue
        vals = vals / vals.max()
        ph = W0 @ vals
        best = np.argsort(ph)[:2]
        cands = [ph.min()]
        starts = [X0[b] for b in best]
        if ph.min() < 0.02:
            starts += [rng.normal(size=8) for _ in range(2)]
        for st in starts:
            r = minimize(phi, st, args=(vals,), method="BFGS", options={"gtol": 1e-10, "maxiter": 300})
            cands.append(r.fun); n_local += 1
            worst.append((r.fun, it, r.x.copy(), (A.copy(), C.copy(), al)))
        worst.sort(key=lambda t: t[0]); worst = worst[:10]
    print(f"   {n_conf} configurations, {n_local} local minimizations, {time.time()-t0:.0f}s")
    wmin = worst[0][0]
    # mpmath re-evaluation of the 10 worst
    mp.mp.dps = 50
    mpw = []
    for val, it, x, (A, C, al) in worst:
        xs = [mp.mpf(float(t)) for t in x]
        pn = mp.sqrt(sum(t**2 for t in xs[:4])); qn = mp.sqrt(sum(t**2 for t in xs[4:]))
        p = [t / pn for t in xs[:4]]; q = [t / qn for t in xs[4:]]
        p0, p1, p2, p3 = p; q0, q1, q2, q3 = q
        Lp = mp.matrix([[p0, -p1, -p2, -p3], [p1, p0, -p3, p2], [p2, p3, p0, -p1], [p3, -p2, p1, p0]])
        Rq = mp.matrix([[q0, q1, q2, q3], [-q1, q0, -q3, q2], [-q2, q3, q0, -q1], [-q3, -q2, q1, q0]])
        O = Lp * Rq
        Am = [mp.mpc(complex(a)) for a in A]; Cm = [mp.mpc(complex(c)) for c in C]
        alm = mp.mpc(complex(al))
        Fm = []
        for s in PERMS:
            t = mp.mpc(1)
            for i in range(4):
                t *= Am[i] - Cm[s[i]]
            Fm.append((alm * t).real)
        mn = min(Fm); vm = [f - mn for f in Fm]; mx = max(vm); vm = [v / mx for v in vm]
        ph = mp.mpf(0)
        for s, sgn_, v in zip(PERMS, SGN, vm):
            ph += sgn_ * v * O[0, s[0]] * O[1, s[1]] * O[2, s[2]] * O[3, s[3]]
        mpw.append(ph)
    print("   10 worst float minima:", ["%.2e" % w[0] for w in worst])
    print("   same points in mpmath (50 digits):", [mp.nstr(v, 3) for v in mpw])
    report(wmin > -1e-9 and min(mpw) > -1e-9,
           f"NUMERIC: {n_conf} random six-on-a-circle configurations (normalized max gap 1): min over SO(4)"
           f" worst value {wmin:.3e}; mpmath re-evaluation worst {mp.nstr(min(mpw), 3)}")


def main():
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    t0 = time.time()
    if stage in ("atoms", "all"):
        print("== stage atoms"); stage_atoms()
    if stage in ("dd", "all"):
        print("== stage dd"); stage_dd()
    if stage in ("chambers", "all"):
        print("== stage chambers"); stage_chambers()
    if stage in ("random", "all"):
        print("== stage random"); stage_random()
    if stage in ("numeric", "all"):
        print("== stage numeric"); stage_numeric(int(os.environ.get("SIXCIRCLE_N", "20000")))
    print(f"stage {stage}: {'ALL PASS' if not FAILS else str(len(FAILS)) + ' FAIL'}  ({time.time()-t0:.1f}s)")
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
