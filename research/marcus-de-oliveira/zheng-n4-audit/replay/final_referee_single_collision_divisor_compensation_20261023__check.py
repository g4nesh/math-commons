"""Exact replay for "Uniform collision identities along a single denominator divisor"
(final_referee_single_collision_divisor_compensation_20261023).  Own code only.  Uses the
sibling own-code script for the shared-zero law (cross collisions).  One line per check;
exit status nonzero iff an EXACT check fails.
"""
import itertools
import os
import random
import runpy
import sys
from fractions import Fraction as Fr

import sympy as sp

here = os.path.dirname(os.path.abspath(__file__))
SZ = runpy.run_path(os.path.join(here, "final_referee_shared_zero_birkhoff_20261022__check.py"))
gaussian_unitary = SZ["gaussian_unitary"]
real_rational_SO4 = SZ["real_rational_SO4"]
abs2 = SZ["abs2"]
law_for_cross_pair = SZ["law_for_cross_pair"]

FAIL = []
PERMS = list(itertools.permutations(range(4)))


def check(name, ok, level="EXACT"):
    print(f"[{level}] {'PASS' if ok else 'FAIL'}  {name}")
    if not ok and level == "EXACT":
        FAIL.append(name)


def sgn(p):
    s = 1
    for i in range(4):
        for j in range(i + 1, 4):
            if p[i] > p[j]:
                s = -s
    return s


# ---------------------------------------------------------------- exact max flow
def max_flow(cap, src, snk):
    nodes = list(cap)
    flow = {u: {v: Fr(0) for v in cap[u]} for u in cap}
    for u in list(cap):
        for v in cap[u]:
            flow.setdefault(v, {}).setdefault(u, Fr(0))
            cap.setdefault(v, {}).setdefault(u, Fr(0))
    total = Fr(0)
    while True:
        par = {src: None}
        q = [src]
        while q and snk not in par:
            u = q.pop(0)
            for v in cap[u]:
                if v not in par and cap[u][v] - flow[u][v] > 0:
                    par[v] = u
                    q.append(v)
        if snk not in par:
            return total, flow
        path = []
        v = snk
        while par[v] is not None:
            path.append((par[v], v))
            v = par[v]
        aug = min(cap[u][v] - flow[u][v] for u, v in path)
        for u, v in path:
            flow[u][v] += aug
            flow[v][u] -= aug
        total += aug


def internal_law(U, cols=(0, 1)):
    """Theorem 2: law for e = (alpha, beta, gamma, gamma) where alpha, beta sit at the columns
    cols of U and gamma at the other two.  Returns (law on S4 in U's own column labels, data)."""
    p = [U[i, cols[0]] for i in range(4)]
    q = [U[i, cols[1]] for i in range(4)]
    a = [abs2(x) for x in p]
    b = [abs2(x) for x in q]
    w = {(i, j): abs2(p[i]*q[j] - p[j]*q[i]) for i in range(4) for j in range(i + 1, 4)}
    big = Fr(10)
    cap = {"s": {}}
    for e, we in w.items():
        cap["s"][("e",) + e] = we
        cap[("e",) + e] = {("v", e[0]): big, ("v", e[1]): big}
    for i in range(4):
        cap[("v", i)] = {"t": a[i]}
    cap["t"] = {}
    val, flow = max_flow(cap, "s", "t")
    theta = {}
    for (i, j) in w:
        theta[(i, j)] = flow[("e", i, j)][("v", i)]   # alpha at i, beta at j
        theta[(j, i)] = flow[("e", i, j)][("v", j)]
    others = [c for c in range(4) if c not in cols]
    law = {}
    for (i, j), th in theta.items():
        if th == 0:
            continue
        rest = [k for k in range(4) if k not in (i, j)]
        for g in (others, others[::-1]):
            sg = [None]*4
            sg[i], sg[j] = cols[0], cols[1]
            sg[rest[0]], sg[rest[1]] = g[0], g[1]
            law[tuple(sg)] = law.get(tuple(sg), 0) + th / 2
    return law, dict(a=a, b=b, w=w, theta=theta, val=val)


def det_poly(U, a, c):
    return sp.expand((sp.diag(*a) - U * sp.diag(*c) * U.H).det(method="berkowitz"))


def law_poly(law, a, c):
    return sp.expand(sum(sp.Rational(w.numerator, w.denominator) *
                         sp.prod([a[i] - c[sg[i]] for i in range(4)]) for sg, w in law.items()))


def main():
    rng = random.Random(1123)
    al, be, ga = sp.symbols("alpha beta gamma")
    x = list(sp.symbols("x0:4"))

    # ---- Lemma 1: L - sum p g = -eta Im(R/e) given sum p = 1 (formal algebra) ----
    Fs = sp.symbols("F0:3")
    Ft, Du, Ee, eta = sp.symbols("Ftau D E eta")
    ps = sp.symbols("p0:3")
    # complex-linear part: (D - Ft)/E - sum p (F - Ft)/E = (D - sum p F)/E when sum p = 1
    expr = (Du - Ft) - sum(ps[k]*(Fs[k] - Ft) for k in range(3)) - (Du - sum(ps[k]*Fs[k] for k in range(3)))
    check("Lemma 1: (D-F_tau) - sum p (F_sigma - F_tau) = D - sum p F_sigma when sum p = 1",
          sp.expand(expr.subs(ps[2], 1 - ps[0] - ps[1])) == 0)

    # ---- Theorem 2 on exact unitaries ----
    Us = [gaussian_unitary(rng) for _ in range(4)] + [real_rational_SO4(rng) for _ in range(2)]
    for idx, U in enumerate(Us):
        law, D = internal_law(U)
        a, b, w, theta = D["a"], D["b"], D["w"], D["theta"]
        ok = sum(w.values()) == 1
        ok &= all(sum(w[tuple(sorted((i, j)))] for j in range(4) if j != i) == a[i] + b[i] for i in range(4))
        check(f"U#{idx}: total edge mass 1 and incident mass a_i + b_i", ok)
        okS = True
        for r in range(1, 5):
            for S in itertools.combinations(range(4), r):
                ins = sum(w[(i, j)] for i, j in itertools.combinations(S, 2))
                ps_ = [U[i, 0] for i in S]
                qs_ = [U[i, 1] for i in S]
                lag = Fr(sum(abs2(z) for z in ps_)) * Fr(sum(abs2(z) for z in qs_)) - \
                    abs2(sum(ps_[k]*sp.conjugate(qs_[k]) for k in range(len(S))))
                okS &= ins == lag and ins <= sum(a[i] for i in S)
        check(f"U#{idx}: for all S, sum_(i<j in S) w_ij = |p_S|^2|q_S|^2 - |<p_S,q_S>|^2 <= sum_S a_i", okS)
        okf = D["val"] == 1 and all(t >= 0 for t in theta.values())
        okf &= all(theta[(i, j)] + theta[(j, i)] == w[(min(i, j), max(i, j))] for i in range(4) for j in range(4) if i != j)
        okf &= all(sum(theta[(i, j)] for j in range(4) if j != i) == a[i] for i in range(4))
        okf &= all(sum(theta[(j, i)] for j in range(4) if j != i) == b[i] for i in range(4))
        check(f"U#{idx}: exact max flow has value 1; theta >= 0, theta_ij+theta_ji = w_ij, out = a_i, in = b_i", okf)
        e = [al, be, ga, ga]
        M = U * sp.diag(*e) * U.H
        okm = True
        for r in range(0, 5):
            for I in itertools.combinations(range(4), r):
                lhs = sp.expand(M.extract(list(I), list(I)).det()) if r else 1
                rhs = sp.expand(sum(sp.Rational(wt.numerator, wt.denominator) * sp.prod([e[sg[i]] for i in I])
                                    for sg, wt in law.items()))
                okm &= sp.expand(lhs - rhs) == 0
        check(f"U#{idx}: law represents all 16 principal minors of U diag(alpha,beta,gamma,gamma) U^* "
              f"({len(law)} weights, sum {sum(law.values())})", okm and sum(law.values()) == 1
              and all(v >= 0 for v in law.values()))
        check(f"U#{idx}: det(diag(x) - U diag(alpha,beta,gamma,gamma) U^*) = sum p prod(x_i - e_sigma(i)) (polynomial identity)",
              sp.expand(det_poly(U, x, e) - law_poly(law, x, e)) == 0)

    # ---- homogeneous divisibility (Corollary 3 / remark on primitive spinors) ----
    O = real_rational_SO4(rng)
    Asym = [sp.symbols(f"A{k}_1 A{k}_2") for k in range(4)]
    Csym = [sp.symbols(f"C{k}_1 C{k}_2") for k in range(4)]
    br = lambda v, w_: v[0]*w_[1] - v[1]*w_[0]
    Pi = {sg: sp.expand(sp.prod([br(Asym[k], Csym[sg[k]]) for k in range(4)])) for sg in PERMS}
    Delta = sp.expand(sum(sgn(sg) * Pi[sg] * sp.prod([O[k, sg[k]] for k in range(4)]) for sg in PERMS))
    # affine chart identification: Delta(A=(a,1), C=(e,1)) = det(diag(a) - O diag(e) O^t)
    ee = list(sp.symbols("e0:4"))
    aff = {}
    for k in range(4):
        aff[Asym[k][0]] = x[k]; aff[Asym[k][1]] = 1
        aff[Csym[k][0]] = ee[k]; aff[Csym[k][1]] = 1
    check("Leibniz form: Delta_O(A,C) on the affine chart equals det(diag(a) - O diag(e) O^t)",
          sp.expand(Delta.xreplace(aff) - det_poly(O, x, ee)) == 0)
    tt = sp.Symbol("t")

    def vanish_on(law, kind, i, j):
        f = Delta - sum(sp.Rational(wt.numerator, wt.denominator) * Pi[sg] for sg, wt in law.items())
        if kind == "cross":      # [A_i, C_j] = 0  <=> A_i = t C_j  (C_j != 0)
            rep = {Asym[i][0]: tt*Csym[j][0], Asym[i][1]: tt*Csym[j][1]}
        elif kind == "C":        # [C_i, C_j] = 0
            rep = {Csym[j][0]: tt*Csym[i][0], Csym[j][1]: tt*Csym[i][1]}
        else:                    # [A_i, A_j] = 0
            rep = {Asym[j][0]: tt*Asym[i][0], Asym[j][1]: tt*Asym[i][1]}
        return sp.expand(f.xreplace(rep)) == 0

    ok = True
    for i0 in range(4):
        for j0 in range(4):
            ok &= vanish_on(law_for_cross_pair(O, i0, j0), "cross", i0, j0)
    check("homogeneous: Delta_O - sum mu Pi vanishes on [A_i,C_j]=0 for all 16 cross pairs (shared-zero law)", ok)
    ok = True
    for c1, c2 in itertools.combinations(range(4), 2):
        others = [c for c in range(4) if c not in (c1, c2)]
        law, _ = internal_law(O, cols=tuple(others))   # alpha, beta at the other two columns
        ok &= vanish_on(law, "C", c1, c2)
    check("homogeneous: vanishes on [C_i,C_j]=0 for all 6 internal pairs (Theorem 2 law)", ok)
    ok = True
    for r1, r2 in itertools.combinations(range(4), 2):
        others = [c for c in range(4) if c not in (r1, r2)]
        lawT, _ = internal_law(O.T, cols=tuple(others))
        law = {}
        for sg, wt in lawT.items():
            inv = [None]*4
            for k in range(4):
                inv[sg[k]] = k
            law[tuple(inv)] = law.get(tuple(inv), 0) + wt
        ok &= vanish_on(law, "A", r1, r2)
    check("homogeneous: vanishes on [A_i,A_j]=0 for all 6 internal pairs (interchanged spectra, inverse perms)", ok)
    print(f"EXACT failures: {len(FAIL)}")
    return O


def numeric_curve(O):
    """Along a_0(t) = e_0 + t^k the quotient (D - sum p F)/t^k stays bounded, while (D - F_tau)/t^k
    does not; mpmath 50 digits."""
    import mpmath as mp
    mp.mp.dps = 50
    rng = random.Random(5)
    law = law_for_cross_pair(O, 0, 0)
    Om = [[mp.mpf(int(O[i, j].p)) / int(O[i, j].q) for j in range(4)] for i in range(4)]
    worst_ratio = 0
    for k in (1, 2, 3):
        base = [mp.mpc(rng.uniform(-2, 2), rng.uniform(-2, 2)) for _ in range(8)]
        e = base[4:]
        vals = []
        naive = []
        for m in (4, 8, 12):
            t = mp.mpf(10)**(-m)
            a = [e[0] + t**k] + base[1:4]
            M = mp.matrix(4, 4)
            for i in range(4):
                for j in range(4):
                    M[i, j] = (a[i] if i == j else 0) - sum(Om[i][l]*e[l]*Om[j][l] for l in range(4))
            Dv = mp.det(M)
            Fv = {sg: mp.fprod([a[i] - e[sg[i]] for i in range(4)]) for sg in PERMS}
            vals.append(abs((Dv - sum(mp.mpf(w.numerator)/w.denominator * Fv[sg] for sg, w in law.items())) / t**k))
            naive.append(abs((Dv - Fv[(1, 0, 2, 3)]) / t**k))
        print(f"[NUMERIC] k={k}: |(D - sum p F)/t^k| at t=1e-4,1e-8,1e-12: {[mp.nstr(v, 6) for v in vals]}; "
              f"without the law |(D - F_tau)/t^k|: {[mp.nstr(v, 3) for v in naive]}")
        worst_ratio = max(worst_ratio, float(vals[-1] / max(vals[0], mp.mpf(10)**-40)))
    check(f"bounded quotient along curves with vanishing order 1,2,3 (ratio last/first {worst_ratio:.3f})",
          worst_ratio < 2, level="NUMERIC")


if __name__ == "__main__":
    O = main()
    numeric_curve(O)
    sys.exit(1 if FAIL else 0)
