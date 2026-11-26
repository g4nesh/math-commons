"""Own-code replay of the "three-plus-three generalized-circle" certificate (sixcircle referee).

Usage:  python -I -B six_circle_replay.py STAGE [SUPP]
  STAGE in {atoms, dd, chambers, numeric, all}; SUPP defaults to $MOC4_SUPP.

Nothing from the supplement is imported or executed.  Its three JSON data files are
loaded with json.load as data only:
  data/final_referee_minor_relaxation_boundary_20261022.json   (72 vertices, list only)
  data/final_referee_six_circle_positive_partition_data_20261022.json  (chambers I, II)
  data/final_referee_six_circle_12atoms_partition_data_20261022.json   (atoms, 86 vertices, chamber III)
Exact objects (P, minor vectors, f_j, group G) come from the audit library ../lib/moc4lib.py.
EXACT stages: atoms, dd, chambers.  NUMERIC stage: numeric (seeds printed).
"""
import os
import sys
import json
import time
import itertools
from fractions import Fraction as Fr
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L  # noqa: E402

T0 = time.monotonic()
FAILS = []


def check(name, ok, info=""):
    tag = "PASS" if ok else "FAIL"
    if not ok:
        FAILS.append(name)
    print(f"[{tag}] {name}" + (f"  ({info})" if info else ""), flush=True)


def supp_path():
    if len(sys.argv) > 2:
        return sys.argv[2]
    p = os.environ.get("MOC4_SUPP")
    if not p:
        sys.exit("give the supplement path as argv[2] or MOC4_SUPP")
    return p


def load(name):
    with open(os.path.join(supp_path(), "data", name), encoding="utf-8") as fh:
        return json.load(fh)


def frvec(v):
    return tuple(Fr(str(a)) for a in v)


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), Fr(0))


# ---------------------------------------------------------------------------------
# Common objects
# ---------------------------------------------------------------------------------
MINOR_C = [L.n_vec(I, J) for (I, J) in L.MINOR_KEYS]            # 34 vectors in V_R (lib coords)
MINOR_VALS = [tuple(int(x) for x in L.values(c)) for c in MINOR_C]
ROWS_P = [tuple(Fr(x) for x in r) for r in L.P]                  # the 24 evaluation functionals
KAPPA = tuple(L.KAPPA)


def master_values():
    """f70 + 2 f2 + f12 (H-threshold atom a = d + d^2/b with d = 1, b = 1, c = 0)."""
    v70, v2, v12 = L.F_VALUES[70], L.F_VALUES[2], L.F_VALUES[12]
    return tuple(v70[s] + 2 * v2[s] + v12[s] for s in range(24))


def g_orbit(vals):
    out = {}
    for g in L.G_ELEMENTS:
        im = tuple(L.g_act_values(g, list(vals)))
        out.setdefault(im, g)
    return out


def atom_data():
    alt = load("final_referee_six_circle_12atoms_partition_data_20261022.json")
    atoms_vals = [tuple(Fr(str(a)) for a in av) for av in alt["fixed_atoms"]]
    return alt, atoms_vals


# ---------------------------------------------------------------------------------
# Stage atoms
# ---------------------------------------------------------------------------------
def stage_atoms():
    print("== stage atoms (EXACT)")
    check("f2 = n_{13,34}, f12 = n_{14,24}, f39 = n_{14,13}, f70 = h_{4231,+} (lib Appendix D)",
          L.F_TABLE[2] == ("n", "13", "34") and L.F_TABLE[12] == ("n", "14", "24")
          and L.F_TABLE[39] == ("n", "14", "13") and L.F_TABLE[70] == ("h", "4231", +1))
    m = master_values()
    orb = g_orbit(m)
    print("   master atom f70 + 2 f2 + f12 values:", m, " G-orbit size", len(orb))
    # the supplement script's own base: 2 n_{12,34} + n_{14,23} + h_{3421,-}
    hv = L.values(L.h_vec(L.perm_from_oneline("3421"), -1))
    sb = tuple(int(2 * a + b + c) for a, b, c in zip(L.n_values((1, 2), (3, 4)), L.n_values((1, 4), (2, 3)), hv))
    check("supplement base 2n_{12,34}+n_{14,23}+h_{3421,-} lies in the G-orbit of f70+2f2+f12", sb in orb)
    alt, atoms = atom_data()
    check("12 distinct stored atoms, each a G-image of the master atom f70+2f2+f12",
          len(set(atoms)) == 12 and all(tuple(int(x) for x in a) in orb for a in atoms))
    for a in atoms:
        L.solve_P(list(a))  # raises if not in im P
    check("each atom value vector lies in im P (exact solve on all 24 rows)", True)
    # minors are G-permuted (validity of the cone under relabelling)
    mset = set(MINOR_VALS)
    check("G permutes the 34 minor value vectors",
          all(tuple(L.g_act_values(g, list(v))) in mset for g in L.G_ELEMENTS[::7] for v in MINOR_VALS))
    # kappa = (1/7) sum of minors, minors span 14  => kappa interior to K
    s = [sum(c[k] for c in MINOR_C) for k in range(14)]
    check("sum of the 34 minor vectors = 7 kappa; minors have rank 14",
          all(s[k] == 7 * KAPPA[k] for k in range(14)) and L.rank([list(c) for c in MINOR_C]) == 14)
    stage_quaternion_identity()


def stage_quaternion_identity():
    """Replay of the KEY identity behind (E1):  Phi(f2+f70) = 2[delta^2+Delta^2-(v0 v2+v1 v3)],
    x2 = det O[13,34] = v0+v1, x3 = det O[14,24] = v1+v2, with O = (x -> p x conj q)."""
    import sympy as S
    p = S.symbols("p0:4")
    q = S.symbols("q0:4")

    def qmul(a, b):
        a0, a1, a2, a3 = a
        b0, b1, b2, b3 = b
        return (a0*b0 - a1*b1 - a2*b2 - a3*b3, a0*b1 + a1*b0 + a2*b3 - a3*b2,
                a0*b2 - a1*b3 + a2*b0 + a3*b1, a0*b3 + a1*b2 - a2*b1 + a3*b0)
    qc = (q[0], -q[1], -q[2], -q[3])
    cols = []
    for k in range(4):
        e = [0, 0, 0, 0]
        e[k] = 1
        cols.append([S.expand(t) for t in qmul(qmul(p, e), qc)])
    O = [[cols[j][i] for j in range(4)] for i in range(4)]
    W = S.Rational(1, 2) * S.Matrix([[1, 1, 1, 1], [1, -1, -1, 1], [1, -1, 1, -1], [1, 1, -1, -1]])
    Pv = list(W * S.Matrix(p))
    Qv = list(W * S.Matrix(q))
    NP = sum(t**2 for t in p)
    NQ = sum(t**2 for t in q)
    v = [S.expand(Pv[i]**2 * NQ - Qv[i]**2 * NP) for i in range(4)]       # homogenized a_i - A_i
    delta = Pv[0]*Qv[3] + Pv[1]*Qv[2]
    Delta = Pv[2]*Qv[1] + Pv[3]*Qv[0]
    phi = S.expand(sum(a * b for a, b in zip([L.F_VALUES[2][s] + L.F_VALUES[70][s] for s in range(24)],
                                               L.leibniz_weights(O))))
    rhs = S.expand(2 * ((delta**2 + Delta**2) * NP * NQ - (v[0]*v[2] + v[1]*v[3])))
    x2 = S.expand(L.minor_det(O, (1, 3), (3, 4)))
    x3 = S.expand(L.minor_det(O, (1, 4), (2, 4)))
    ok_o = all(S.expand(sum(O[k][i]*O[k][j] for k in range(4)) - (NP*NQ if i == j else 0)) == 0
               for i in range(4) for j in range(4))
    print("   O^T O = |p|^2|q|^2 I:", ok_o)
    best = None
    for s1, s2 in itertools.product([1, -1], repeat=2):
        if S.expand(phi - rhs) == 0 and S.expand(x2 - s1*(v[0]+v[1])) == 0 and S.expand(x3 - s2*(v[1]+v[2])) == 0:
            best = (s1, s2)
    check("KEY quaternion identity Phi(f2+f70) = 2[(d^2+D^2)|P|^2|Q|^2 - (v0v2+v1v3)] and "
          "x2 = +-(v0+v1), x3 = +-(v1+v2), as polynomial identities", ok_o and best is not None,
          f"signs {best}")


# ---------------------------------------------------------------------------------
# Stage dd: own double description of the dual polytopes
# ---------------------------------------------------------------------------------
def to_int(v):
    den = 1
    for x in v:
        den = den * Fr(x).denominator // gcd(den, Fr(x).denominator)
    w = [int(Fr(x) * den) for x in v]
    g = 0
    for x in w:
        g = gcd(g, abs(x))
    return tuple(x // g for x in w) if g else tuple(w)


def dd_extreme_rays(cons, d=14):
    """Extreme rays of {w in R^d : <a, w> >= 0 for a in cons} (pointed). Exact integers.
    Incremental double description with the combinatorial adjacency test."""
    cons = [to_int(a) for a in cons]
    # initial basis of d independent constraints
    basis = []
    for i, a in enumerate(cons):
        if L.rank([list(cons[j]) for j in basis + [i]]) == len(basis) + 1:
            basis.append(i)
        if len(basis) == d:
            break
    A0 = [[Fr(x) for x in cons[i]] for i in basis]
    Ainv = L.inv(A0)
    rays = []
    for j in range(d):
        r = to_int([Ainv[i][j] for i in range(d)])
        z = 0
        for k, i in enumerate(basis):
            if k != j:
                z |= 1 << i
        rays.append((r, z))
    order = basis + [i for i in range(len(cons)) if i not in basis]
    done = set(basis)
    maxr = len(rays)
    for i in order[d:]:
        a = cons[i]
        vals = [sum(x * y for x, y in zip(a, r)) for r, _ in rays]
        pos = [k for k, s in enumerate(vals) if s > 0]
        neg = [k for k, s in enumerate(vals) if s < 0]
        zer = [k for k, s in enumerate(vals) if s == 0]
        new = []
        zsets = [z for _, z in rays]
        for kp in pos:
            for kn in neg:
                Z = zsets[kp] & zsets[kn]
                if bin(Z).count("1") < d - 2:
                    continue
                adj = True
                for k in range(len(rays)):
                    if k != kp and k != kn and (zsets[k] & Z) == Z:
                        adj = False
                        break
                if not adj:
                    continue
                sp, sn = vals[kp], vals[kn]
                r = to_int([sp * y - sn * x for x, y in zip(rays[kp][0], rays[kn][0])])
                new.append((r, Z | (1 << i)))
        rays = [(rays[k][0], rays[k][1] | (1 << i)) for k in zer] + [rays[k] for k in pos] + new
        done.add(i)
        maxr = max(maxr, len(rays))
    return [r for r, _ in rays], maxr


def normalize_vertex(r):
    k = r[13]  # <r, kappa> with kappa = e_13
    assert k > 0
    return tuple(Fr(x, k) for x in r)


def check_vertex_set(name, cons, listed):
    t = time.monotonic()
    rays, maxr = dd_extreme_rays(cons)
    mine = {normalize_vertex(r) for r in rays}
    print(f"   {name}: DD gives {len(mine)} vertices (max intermediate rays {maxr}, "
          f"{time.monotonic() - t:.1f}s)")
    # each listed point is feasible and a vertex (active constraints of rank 13)
    consF = [tuple(Fr(x) for x in a) for a in cons]
    okv = all(all(dot(a, v) >= 0 for a in consF) and dot(KAPPA, v) == 1
              and L.rank([list(a) for a in consF if dot(a, v) == 0]) == 13 for v in listed)
    check(f"{name}: every stored vertex is feasible with an active set of rank 13", okv)
    # independent NUMERIC cross-check of completeness: random LPs over the polytope
    import numpy as np
    from scipy.optimize import linprog
    rng = np.random.default_rng(7)
    Af = np.array([[float(x) for x in a] for a in consF])
    Vf = np.array([[float(x) for x in v] for v in listed])
    kap = np.array([float(x) for x in KAPPA])
    worst = 0.0
    for _ in range(400):
        obj = rng.standard_normal(14)
        res = linprog(obj, A_ub=-Af, b_ub=np.zeros(len(Af)), A_eq=kap[None, :], b_eq=[1.0],
                      bounds=(None, None), method="highs")
        worst = max(worst, abs(res.fun - (Vf @ obj).min()))
    check(f"{name}: 400 random LPs attain their optimum at a listed vertex (NUMERIC)", worst < 1e-7,
          f"max discrepancy {worst:.1e}")
    return mine, mine == set(listed)


def stage_dd():
    print("== stage dd (EXACT double description, own code)")
    rel = load("final_referee_minor_relaxation_boundary_20261022.json")
    pos = load("final_referee_six_circle_positive_partition_data_20261022.json")
    alt, atoms_vals = atom_data()
    v72 = [frvec(v) for v in rel["facet_normals"]]
    v72b = [frvec(v) for v in pos["target_coordinates"]["minor_relaxation"]]
    v86 = [frvec(v) for v in alt["vertices"]]
    check("stored 72-list (relaxation file) equals the chamber I/II target list, same order", v72 == v72b)
    check("the 24 rows of P (evaluation functionals) are among the 72 and the 86",
          set(ROWS_P) <= set(v72) and set(ROWS_P) <= set(v86))
    mine72, ok = check_vertex_set("R = {w : w(kappa)=1, w(n) >= 0 for 34 minors}", MINOR_C, v72)
    check("R has exactly the 72 stored vertices (own DD, completeness included)", ok and len(v72) == 72)
    atoms_c = [L.solve_P(list(a)) for a in atoms_vals]
    mine86, ok = check_vertex_set("R' = R cut by the 12 atoms", MINOR_C + atoms_c, v86)
    check("R' has exactly the 86 stored vertices (own DD, completeness included)", ok and len(v86) == 86)
    # which atoms are needed? and the two harmful vertices of R
    harmful = sorted({i for i, v in enumerate(v72) for a in atoms_c if dot(a, v) < 0})
    print("   vertices of R cut off by some atom (indices in the 72-list):", harmful)


# ---------------------------------------------------------------------------------
# Stage chambers: exact replay of the 158 positive partitions
# ---------------------------------------------------------------------------------
def build_V():
    """F_sigma for A = (inf, 0, 1, zeta), C = (s, t, u, eta), homogeneous brackets in the
    lib convention [i,4+j] = x_i y_j - x_j y_i with C_j represented by (C_j, 1).
    Returns the 24 x 4 matrix of coefficients of (1, zeta, eta, zeta*eta)."""
    import sympy as S
    s, t, u, ze, et = S.symbols("s t u zeta eta")
    Apts = [(1, 0), (0, 1), (1, 1), (ze, 1)]
    Cpts = [(s, 1), (t, 1), (u, 1), (et, 1)]
    V = []
    for p in L.PERMS:
        F = 1
        for i in range(4):
            (x1, y1), (x2, y2) = Apts[i], Cpts[p[i]]
            F *= (x1 * y2 - x2 * y1)
        F = S.Poly(S.expand(F), ze, et)
        V.append([F.coeff_monomial(m) for m in (1, ze, et, ze * et)])
    return V, (s, t, u)


def stage_chambers():
    import sympy as S
    import flint
    print("== stage chambers (EXACT polynomial identities)")
    V, (s, t, u) = build_V()
    # columns of V lie in im P identically: V_col = P (Ppiv^{-1} V_col[piv])
    piv = L.P_PIVOT_ROWS
    Pinv = L.inv([[Fr(x) for x in L.P[r]] for r in piv])
    okim = True
    for col in range(4):
        vc = [S.Rational(x) if not isinstance(x, S.Basic) else x for x in [V[r][col] for r in range(24)]]
        c = [sum(S.Rational(Pinv[a][b].numerator, Pinv[a][b].denominator) * vc[piv[b]] for b in range(14))
             for a in range(14)]
        for r in range(24):
            if S.expand(sum(L.P[r][a] * c[a] for a in range(14)) - vc[r]) != 0:
                okim = False
    check("each coefficient column of (F_sigma) is P of a real vector, identically in s,t,u", okim)
    x, y, z = S.symbols("x y z")
    chambers = {"I": ({s: x/(1+x+y+z), t: (x+y)/(1+x+y+z), u: (x+y+z)/(1+x+y+z)}, (1+x+y+z)**3),
                "II": ({s: x/(1+x+y), t: (x+y)/(1+x+y), u: 1+z}, (1+x+y)**2),
                "III": ({s: -x, t: y/(1+y), u: 1+z}, 1+y)}
    ctx = flint.fmpq_mpoly_ctx.get(("x", "y", "z"), "lex")
    X, Y, Z = ctx.gens()

    def to_flint(expr):
        e = S.cancel(expr)
        num, den = S.fraction(e)
        assert den.free_symbols == set(), ("not polynomial", e)
        pol = S.Poly(num, x, y, z)
        out = ctx.from_dict({}) if False else 0 * X
        for (a, b, c), co in pol.terms():
            co = S.Rational(co) / den
            out += flint.fmpq(int(co.p), int(co.q)) * X**a * Y**b * Z**c
        return out

    def stored(terms):
        out = 0 * X
        seen = set()
        for power, value in terms:
            power = tuple(power)
            val = Fr(value)
            if not (len(power) == 3 and all(k in (0, 1) for k in power) and val >= 0 and power not in seen):
                return None
            seen.add(power)
            out += flint.fmpq(val.numerator, val.denominator) * X**power[0] * Y**power[1] * Z**power[2]
        return out

    pos = load("final_referee_six_circle_positive_partition_data_20261022.json")
    alt, _ = atom_data()
    targets_initial = [frvec(v) for v in pos["target_coordinates"]["minor_relaxation"]]
    targets_final = [frvec(v) for v in alt["vertices"]]
    rowsP = set(ROWS_P)
    total = 0
    for name, (sub, den) in chambers.items():
        tt = time.monotonic()
        polys = [[to_flint(den * S.sympify(V[r][c]).subs(sub)) for c in range(4)] for r in range(24)]
        targets = targets_final if name == "III" else targets_initial
        recs = alt["records"] if name == "III" else pos["cells"][name]["records"]
        need = {k for k, w in enumerate(targets) if w not in rowsP}
        seen = set()
        bad = 0
        nterms = 0
        for rec in recs:
            k = rec["target"]
            if k not in need or k in seen:
                bad += 1
                continue
            seen.add(k)
            R = stored(rec["R"])
            Wd = {int(j): stored(v) for j, v in rec["W"].items()}
            if R is None or R == 0 * X or any(w is None for w in Wd.values()) or not all(0 <= j < 24 for j in Wd):
                bad += 1
                continue
            ssum = 0 * X
            for w in Wd.values():
                ssum += w
            if ssum != R:              # kappa column: sum of weights = R  (w(kappa) = v_sigma(kappa) = 1)
                bad += 1
                continue
            w = targets[k]
            alpha = [sum(w[a] * Pinv[a][b] for a in range(14)) for b in range(14)]  # w = alpha . P[piv]
            for c in range(4):
                tgt = 0 * X
                for b in range(14):
                    if alpha[b]:
                        tgt += flint.fmpq(alpha[b].numerator, alpha[b].denominator) * polys[piv[b]][c]
                lhs = 0 * X
                for j, wj in Wd.items():
                    lhs += wj * polys[j][c]
                if lhs != R * tgt:
                    bad += 1
                    break
            nterms += len(rec["R"]) + sum(len(v) for v in rec["W"].values())
        # negative control: the first record checked against the wrong target must fail
        rec = recs[0]
        R = stored(rec["R"]); Wd = {int(j): stored(v) for j, v in rec["W"].items()}
        wrong = targets[(rec["target"] + 1) % len(targets)]
        alpha = [sum(wrong[a] * Pinv[a][b] for a in range(14)) for b in range(14)]
        detect = False
        for c in range(4):
            tgt = 0 * X
            for b in range(14):
                if alpha[b]:
                    tgt += flint.fmpq(alpha[b].numerator, alpha[b].denominator) * polys[piv[b]][c]
            lhs = 0 * X
            for j, wj in Wd.items():
                lhs += wj * polys[j][c]
            detect = detect or (lhs != R * tgt)
        check(f"chamber {name}: negative control (record 0 against a wrong target vertex) is rejected", detect)
        total += len(recs)
        check(f"chamber {name}: {len(recs)} partitions cover all {len(need)} non-ordinary vertices, "
              f"nonnegative multiaffine weights, sum W = R != 0, 4 column identities exact",
              bad == 0 and seen == need, f"{nterms} stored terms, {time.monotonic() - tt:.1f}s")
    check("158 partitions in total (48 + 48 + 62)", total == 158)


# ---------------------------------------------------------------------------------
# Stage numeric
# ---------------------------------------------------------------------------------
def stage_numeric():
    import numpy as np
    from scipy.optimize import linprog
    print("== stage numeric (NUMERIC; floating point)")
    Pn = np.array(L.P, float)
    piv = L.P_PIVOT_ROWS
    Ppinv = np.linalg.inv(Pn[piv])
    alt, atoms_vals = atom_data()
    m = master_values()
    orbit = np.array(sorted(g_orbit(m).keys()), float)                 # all valid atoms (value space)
    MV = np.array(MINOR_VALS, float)
    v86 = np.array([[float(x) for x in frvec(v)] for v in alt["vertices"]])
    v72 = np.array([[float(x) for x in frvec(v)] for v in load(
        "final_referee_six_circle_positive_partition_data_20261022.json")["target_coordinates"]["minor_relaxation"]])
    S4 = np.array(L.PERMS)
    SG = np.array(L.SGN, float)

    def gaps_from_points(A, C, alpha):
        F = alpha * np.prod(A[None, :] - C[S4], axis=1)
        re = F.real
        return re - re.min()

    def in_cone(g, gens):
        res = linprog(np.zeros(len(gens)), A_eq=gens.T, b_eq=g, bounds=(0, None), method="highs")
        return res.status == 0

    def qmat(th):
        # O = (x -> p x conj q) for unit quaternions p, q (onto SO(4))
        p = th[:4] / np.linalg.norm(th[:4])
        q = th[4:] / np.linalg.norm(th[4:])
        a, b, c, d = p
        Lp = np.array([[a, -b, -c, -d], [b, a, -d, c], [c, d, a, -b], [d, -c, b, a]])
        a, b, c, d = q
        Rq = np.array([[a, b, c, d], [-b, a, -d, c], [-c, d, a, -b], [-d, -c, b, a]])  # x -> x conj q
        return Lp @ Rq

    def phimin(g, rng, restarts=12):
        from scipy.optimize import minimize
        f = lambda th: float(g @ L.leibniz_weights_np(qmat(th)))
        best = np.inf
        for _ in range(restarts):
            r = minimize(f, rng.standard_normal(8), method="BFGS")
            best = min(best, r.fun)
        return best

    Ot = qmat(rng_tmp := np.random.default_rng(1).standard_normal(8))
    assert np.allclose(Ot.T @ Ot, np.eye(4)) and abs(np.linalg.det(Ot) - 1) < 1e-12

    def mobius(pts, M):
        a, b, c, d = M
        return (a * pts + b) / (c * pts + d)

    gens_all = np.vstack([MV, orbit])
    rng = np.random.default_rng(20261125)
    # (N1) realizable configurations: six concyclic points (random circle), random
    #      positions of the two free points, random direction; reference = argmin.
    N1, fails, worst_phi, nphi, need_atoms = 4000, 0, np.inf, 0, 0
    worst_v86 = np.inf
    for it in range(N1):
        real6 = rng.standard_normal(6) * rng.choice([0.3, 1, 3])
        if it % 5 == 0:                      # near-coincidences
            real6[rng.integers(6)] = real6[rng.integers(6)] + 1e-3 * rng.standard_normal()
        M = rng.standard_normal(4) + 1j * rng.standard_normal(4)
        six = mobius(real6.astype(complex), M)                  # on a random generalized circle
        zf, wf = rng.standard_normal(2) + 1j * rng.standard_normal(2)
        ia, ic = rng.integers(4), rng.integers(4)                # which A and C are free
        A = np.insert(six[:3], ia, zf)
        C = np.insert(six[3:], ic, wf)
        alpha = np.exp(2j * np.pi * rng.random())
        g = gaps_from_points(A, C, alpha)
        g = g / max(1e-300, np.abs(g).max())
        if not in_cone(g, gens_all):
            fails += 1
        elif not in_cone(g, MV):
            need_atoms += 1
        if it < 200:
            worst_phi = min(worst_phi, phimin(g, rng, restarts=8))
            nphi += 1
        # chamber-frame test with the 86 vertices when the configuration is already in
        # normal form A = (inf,0,1,.) is not available here; done in N2.
    check(f"N1: {N1} random concyclic 3+3 configurations (random circle, random free slots, "
          f"random direction): gap vector in cone(minors + G-orbit of master atom) by LP",
          fails == 0, f"failures {fails}; {need_atoms} of them are not in the minor cone alone")
    check(f"N1: min Phi over SO(4) (8 BFGS starts, quaternion chart) on the first {nphi}: worst {worst_phi:.3e}",
          worst_phi > -1e-9)
    # (N2) the full 5-dim profile family in the three chambers and on degenerate parameters:
    #      enumerate extreme rays of {c in L(s,t,u) : Pc >= 0} and test the 86 (resp. 72) vertices.
    def Vnum(s, t, u):
        out = np.zeros((24, 4))
        for k, p in enumerate(L.PERMS):
            c1, c2, c3 = [(s, t, u, None)[p[i]] for i in (1, 2, 3)]
            # F = (-C_{p1})(1 - C_{p2})(zeta - C_{p3}) with C_4 = eta (index 3)
            # expand in (1, zeta, eta, zeta*eta)
            def lin(ci):  # returns coefficients (const, eta) of the factor value C
                return (0.0, 1.0) if ci is None else (ci, 0.0)
            a1, a2, a3 = lin(c1), lin(c2), lin(c3)
            # f1 = -C1 ; f2 = 1 - C2 ; f3 = zeta - C3 ; each C = const + e*eta
            f1 = np.array([-a1[0], 0, -a1[1], 0])          # basis (1, zeta, eta, zeta*eta)
            f2 = np.array([1 - a2[0], 0, -a2[1], 0])
            f3 = np.array([-a3[0], 1, -a3[1], 0])

            def mul(f, g):
                r = np.zeros(4)
                for i in range(4):
                    for j in range(4):
                        if f[i] == 0 or g[j] == 0:
                            continue
                        zi, ei = i in (1, 3), i in (2, 3)
                        zj, ej = j in (1, 3), j in (2, 3)
                        if (zi and zj) or (ei and ej):
                            raise ValueError("degree > 1")
                        r[(1 if (zi or zj) else 0) + (2 if (ei or ej) else 0)] += f[i] * g[j]
                return r
            out[k] = mul(mul(f1, f2), f3)
        return out

    def profile_rays(s, t, u):
        Vm = np.hstack([Vnum(s, t, u), np.ones((24, 1))])      # values of c = B r + k kappa
        Vm = Vm / np.linalg.norm(Vm, axis=0, keepdims=True).clip(1e-300)
        rays = []
        for S4sub in itertools.combinations(range(24), 4):
            M = Vm[list(S4sub)]
            _, sv, vt = np.linalg.svd(M)
            if sv[-1] < 1e-9 * sv[0]:
                continue
            r = vt[-1]
            gv = Vm @ r
            if gv.min() < -1e-9:
                r, gv = -r, -gv
            if gv.min() < -1e-9 or np.abs(gv).max() < 1e-12:
                continue
            rays.append(gv / np.abs(gv).max())
        return rays

    def test_rays(rays, verts):
        worst = np.inf
        for gv in rays:
            c = Ppinv @ gv[piv]
            worst = min(worst, float((verts @ c).min()))
        return worst

    samples = []
    for _ in range(12):
        a, b, c = np.sort(rng.random(3))
        samples.append(("I", a, b, c))
        a, b = np.sort(rng.random(2))
        samples.append(("II", a, b, 1 + rng.exponential()))
        samples.append(("III", -rng.exponential(), rng.random(), 1 + rng.exponential()))
    samples += [("I", 0.2, 0.2, 0.7), ("I", 0.0, 0.5, 1.0), ("II", 0.3, 1.0, 1.0), ("III", 0.0, 0.5, 2.0),
                ("III", -1.0, 0.0, 1.0), ("III", -2.0, 0.5, 0.5), ("II", 0.4, 0.4, 3.0), ("I", 0.5, 0.5, 0.5)]
    worst = {"I": np.inf, "II": np.inf, "III": np.inf}
    worst_phi2, nrays = np.inf, 0
    for ch, a, b, c in samples:
        rays = profile_rays(a, b, c)
        nrays += len(rays)
        verts = v86 if ch == "III" else v72
        worst[ch] = min(worst[ch], test_rays(rays, verts))
        for gv in rays[:: max(1, len(rays) // 3)][:3]:
            worst_phi2 = min(worst_phi2, phimin(gv, rng, restarts=6))
    check(f"N2: extreme rays of the profile cone at {len(samples)} chamber parameters "
          f"({nrays} rays, incl. closure/degenerate points): min over chamber vertices of w(c)",
          min(worst.values()) > -1e-7, ", ".join(f"{k}: {v:.2e}" for k, v in worst.items()))
    check(f"N2: min Phi over SO(4) on a subset of those rays (6 starts each): {worst_phi2:.3e}",
          worst_phi2 > -1e-9)
    # (N3) control: chamber III needs the atoms: vertices of R (72) fail on some profile ray
    w3 = np.inf
    hard = []
    for ch, a, b, c in samples:
        if ch == "III":
            for gv in profile_rays(a, b, c):
                cc = Ppinv @ gv[piv]
                m72 = float((v72 @ cc).min())
                w3 = min(w3, m72)
                if m72 < -1e-6:
                    hard.append(gv)
    print(f"   control N3: chamber III rays tested against the 72 minor-only vertices: min {w3:.3e}"
          f"  (negative means minors alone do not suffice there); {len(hard)} such rays")
    wh = np.inf
    for gv in hard[:60]:
        wh = min(wh, phimin(gv, rng, restarts=8))
    check(f"N3: min Phi over SO(4) on {min(60, len(hard))} profile rays outside the minor cone "
          f"(8 starts each)", wh > -1e-9, f"worst {wh:.3e}")
    # (N4) control: non-concyclic configurations are not all in the cone
    out = 0
    for it in range(500):
        A = rng.standard_normal(4) + 1j * rng.standard_normal(4)
        C = rng.standard_normal(4) + 1j * rng.standard_normal(4)
        g = gaps_from_points(A, C, np.exp(2j * np.pi * rng.random()))
        g = g / np.abs(g).max()
        out += not in_cone(g, gens_all)
    print(f"   control N4: 500 generic (non-concyclic) configurations: {out} outside cone(minors + atom orbit)")


def main():
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    if stage in ("atoms", "all"):
        stage_atoms()
    if stage in ("dd", "all"):
        stage_dd()
    if stage in ("chambers", "all"):
        stage_chambers()
    if stage in ("numeric", "all"):
        stage_numeric()
    print(f"stage {stage}: {'ALL PASS' if not FAILS else 'FAILURES: ' + '; '.join(FAILS)}"
          f"  ({time.monotonic() - T0:.1f}s)")
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
