"""Audit replay for independent_moc4_T8_double_boundary_20261022
(contact set T_9 u {18,20} = {0,1,2,3,11,12,18,20}; supplement label T8 = paper T_9).

Own code only; uses the shared library ../lib/moc4lib.py.  Nothing from the supplement is
read or executed.  Prints one line per check; exits nonzero if an EXACT check fails.
Run:  python -I -B independent_moc4_T8_double_boundary_20261022__check.py
The chart helpers here are also imported by the T8_uniform replay (Section 3).
"""
import os
import sys
import random

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L  # noqa: E402
import sympy as sp  # noqa: E402

FAIL = []


def check(name, ok, info=""):
    print(f"{'PASS' if ok else 'FAIL'} EXACT {name} {info}".rstrip())
    if not ok:
        FAIL.append(name)


z, w, v, vb, x, t, r, s = sp.symbols("z w v vb x t r s")
Pp = 1 - t
M = s + x - 1
Hh = r + s - Pp
Ll = Hh + x
Bch = r * (x - t) + Pp * (s - Pp)
Dch = r * (x - t) - t * (s - Pp)
Qch = (r + t - 1) * (x - 1) - s
Q0 = (r - Pp) * (x - 1) - s
ZINV = (t * Ll * v - x * (r + t)) / (s + t * v - 1)
WINV = ZINV * v * (s + t * v - 1) / (v * (s + t + x - 1) - x)
Avec = [z, 0, 1, x]
Cvec = [w, v, t * v, None]


def Fsig(sidx):
    p = sp.Integer(1)
    for i in range(4):
        j = L.PERMS[sidx][i]
        if Cvec[j] is None:
            continue
        p *= (Avec[i] - Cvec[j])
    return sp.expand(p)


HRAW = [sp.cancel((Fsig(k) - Fsig(0)) / (v * (z - w))) for k in range(24)]
HCH = [sp.cancel(h.subs({z: ZINV, w: WINV})) for h in HRAW]     # functions of v on the chart


def DD(j):
    """g_j / rho where rho = mu Im v:  -(H(v) - H(vb)) / (v - vb)."""
    h = HCH[j]
    return sp.cancel(-(h - h.subs(v, vb)) / (v - vb))


N = v * vb
O = (v - 1) * (vb - 1)
Dt = (t * v - x) * (t * vb - x)
Dp = (t * Ll * v - x * (r + t)) * (t * Ll * vb - x * (r + t))


def zero(e):
    return sp.cancel(sp.together(e)) == 0


def run_chart_checks(prefix=""):
    check(prefix + "chart: H12 = 1-t-r and H18 = s identically under the inverse formulas",
          zero(HCH[12] - (1 - t - r)) and zero(HCH[18] - s))
    check(prefix + "chart: r = z(v-1)(tv-w)/(v(z-w)) (same r as the polar chart)",
          zero((Pp - HRAW[12]) - z * (v - 1) * (t * v - w) / (v * (z - w))))
    check(prefix + "g4 = rho x/N", zero(DD(4) - x / N))
    check(prefix + "g7 = rho B/O", zero(DD(7) - Bch / O))
    check(prefix + "g13 = rho r(x-1)/O", zero(DD(13) - r * (x - 1) / O))
    check(prefix + "g23 = rho r t x M H/Dp", zero(DD(23) - r * t * x * M * Hh / Dp))
    check(prefix + "M g4 + s g7 - Q g20 - M g11 = 0",
          zero(M * DD(4) + s * DD(7) - Qch * DD(20) - M * DD(11)))
    check(prefix + "(s-P) g7 - M g23 - D g20 = 0", zero((s - Pp) * DD(7) - M * DD(23) - Dch * DD(20)))


if __name__ == "__main__":
    T = frozenset({0, 1, 2, 3, 11, 12})
    run_chart_checks()
    check("Q (double-boundary) = Q0 (uniform Thm 3)", sp.expand(Qch - Q0) == 0)
    # SOC profile identity with -A f18, A = g18 - g4 (equivalently +a f18, a = g4 - g18)
    BASIS = [1, 2, 7, 8, 9, 12, 15, 18]
    coef = {1: 23, 2: 21, 7: 20, 8: 9, 9: 7, 12: 6, 15: 18}
    gs = sp.symbols("g0:24")
    gen = [0] * 24
    for k, gi in coef.items():
        for q in range(24):
            gen[q] += L.F_VALUES[k][q] * gs[gi]
    for q in range(24):
        gen[q] += L.F_VALUES[18][q] * (gs[4] - gs[18])
    # gen must reproduce every coordinate of a vector in N_T: check on the 8 pivots
    piv = [4, 6, 7, 8, 9, 13, 16, 18]
    # express pivots via coordinates and verify consistency on the 16 other gap forms
    Kc, Lc, Jc, Qc, Cc, Bc, hc, Ac = sp.symbols("K L J Q C B h A")
    SY = [Kc, Lc, Jc, Qc, Cc, Bc, hc, Ac]
    FORMS = [sum(sp.Integer(L.F_VALUES[k][q]) * SY[i] for i, k in enumerate(BASIS)) for q in range(24)]
    sub = {gs[q]: FORMS[q] for q in range(24)}
    check("profile identity g = g23 f1 + g21 f2 + g20 f7 + g9 f8 + g7 f9 + g6 f12 + g18 f15 + (g4-g18) f18 on N_T",
          all(sp.expand(gen[q].subs(sub) - FORMS[q]) == 0 for q in range(24)))
    # boundary relations (g18 = g20 = 0): g14 = g9 + a, b = g9 - x g14 (chart), g21 = g13 - a - b
    check("forms: g14 = g9 + g4 when g18 = g20 = 0; g21 = g13 - g4 - g7 when g18 = 0",
          sp.expand((FORMS[14] - FORMS[9] - FORMS[4]).subs({hc: 0, Jc: 0})) == 0
          and sp.expand((FORMS[21] - FORMS[13] + FORMS[4] + FORMS[7]).subs({hc: 0})) == 0)
    # chart: g7 - g9 + x g14 is a multiple of g20 (so b = g9 - x g14 on the boundary)
    lin = sp.cancel((DD(7) - DD(9) + x * DD(14)) / DD(20))
    check("chart: g7 - g9 + x g14 = lambda * g20 with lambda free of v, vb",
          not lin.has(v) and not lin.has(vb), f"lambda={sp.factor(lin)}")
    # on the boundary g20 = g11 = 0: from the two linear relations a/b = -s/M, g23/b = (s-P)/M,
    # g14/b = -1/M, g9/b = (s-1)/M.  Check with g20 = 0 eliminated: express g14, g9 via a, b.
    a_, b_ = sp.symbols("a b")
    sols = {"a": -s * b_ / M}
    g14v = sp.cancel((b_ / x) * 0 + (sp.Symbol("g9") + sols["a"]))
    # b = g9 - x g14 = g9 - x (g9 + a)  =>  g9 = (b + x a)/(1 - x)
    g9v = sp.cancel((b_ + x * sols["a"]) / (1 - x))
    check("boundary: g9/b = (s-1)/M and g14/b = -1/M",
          zero(g9v / b_ - (s - 1) / M) and zero((g9v + sols["a"]) / b_ + 1 / M))
    # g21/b = (x-1)(r-P)(s-P)/(M B) on the boundary: g21 = g13 - a - b with g13 = rho r(x-1)/O,
    # g7 = rho B/O => g13/b = r(x-1)/B
    check("boundary: g21/b = r(x-1)/B + s/M - 1 = (x-1)(r-P)(s-P)/(M B)",
          zero(r * (x - 1) / Bch + s / M - 1 - (x - 1) * (r - Pp) * (s - Pp) / (M * Bch)))
    # polynomial compensation identities (4), (5)
    Cc2 = Pp * (1 - x) - s * (2 - t - x)
    Jp = (r * s ** 2 * x - r * s ** 2 + r * s * t * x - 2 * r * s * x + r * s + s ** 2 * t * x - s ** 2 * t
          + s * t ** 2 * x - s * t ** 2 - 2 * s * t * x + 2 * s * t + s * x ** 2 - s * x - t ** 2 * x + t ** 2
          + t * x ** 2 - t - x ** 2 + x)
    B0 = r * ((s - Pp - 1) * x + 1 - s) + Pp * (Pp - s) * x
    A1 = M * x * O + s * Bch * N
    E20 = sp.cancel(DD(20) * O * Dp)
    E6 = sp.cancel(DD(6) * O * Dt)
    check("E20 = g20 O Dp / rho is a polynomial", sp.denom(sp.together(E20)).is_number)
    check("E6 = g6 O Dt / rho is a polynomial", sp.denom(sp.together(E6)).is_number)
    check("A1 = (M g4 + s g7) N O / rho", zero(A1 - (M * DD(4) + s * DD(7)) * N * O))
    check("(4): B{(r+t)(s-P)Dt - t J N} = (x-t)E20 - t B0 A1",
          sp.expand(Bch * ((r + t) * (s - Pp) * Dt - t * Jp * N) - (x - t) * E20 + t * B0 * A1) == 0)
    check("(5): J E6 - P{C + s(P-r)} B Dt = tP(x-1)(s-P) M A1 + sP(x-1) E20",
          sp.expand(Jp * E6 - Pp * (Cc2 + s * (Pp - r)) * Bch * Dt
                    - t * Pp * (x - 1) * (s - Pp) * M * A1 - s * Pp * (x - 1) * E20) == 0)
    check("capacity: (s-P)M - s^2 = C = P(1-x) - s(2-t-x)", sp.expand((s - Pp) * M - s ** 2 - Cc2) == 0)
    # SOC lemma for (f1,f9,f18) numerically: Phi_f18 >= -2 sqrt(Phi_f1 Phi_f9) on SO(4)
    import numpy as np
    rng = np.random.default_rng(7)
    worst = 1e9
    f1, f9, f18 = [np.array([float(c) for c in L.F_VECS[j]]) for j in (1, 9, 18)]
    for _ in range(20000):
        Oq = L.so4_from_params(rng.normal(size=6) * 2)
        p1, p9, p18 = (L.Phi_np(f, Oq) for f in (f1, f9, f18))
        worst = min(worst, p18 + 2 * np.sqrt(max(p1, 0) * max(p9, 0)))
    print(f"NUMERIC SOC(f1,f9,f18): min over 20000 random SO(4) (seed 7) of Phi18 + 2 sqrt(Phi1 Phi9) = {worst:.3e}")
    print("SUMMARY", "FAILED: " + ", ".join(FAIL) if FAIL else "all EXACT checks passed")
    sys.exit(1 if FAIL else 0)
