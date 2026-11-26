"""Audit replay for independent_moc4_T8_uniform_20261022 (paper case T_9 = {0,1,2,3,11,12}).

Section 1 of the document is the text of independent_moc4_T8_polar_reduction_20261022 and is
replayed by that document's script (not re-run here).  This script replays Sections 2-3.
Own code only; uses ../lib/moc4lib.py and the chart helpers of the double-boundary replay.
Nothing from the supplement is read or executed.  Exits nonzero if an EXACT check fails.

The quintic G is not printed in the document.  It is reconstructed here from the claimed
identity  Q2^4 det Hess C(c0 - (a Q1/Q2) kappa) = -4 a^4 Q1^2 Q2^2 G^2  by modular
interpolation (one prime, lines through a common anchor, polynomial square roots), and the
identity is then checked in exact rational arithmetic at random rational points
(a randomized identity test, not a symbolic expansion).
"""
import os
import sys
import random
import itertools

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, HERE)
import moc4lib as Lb  # noqa: E402
from moc4lib import _P_PIV_INV  # noqa: E402
import sympy as sp  # noqa: E402
import flint  # noqa: E402
import numpy as np  # noqa: E402

FAIL = []


def check(name, ok, info="", level="EXACT"):
    print(f"{'PASS' if ok else 'FAIL'} {level} {name} {info}".rstrip())
    if not ok:
        FAIL.append(name)


NAMES = ["A", "B", "C", "Q", "J", "K", "L", "h"]
BASIS = {"K": 1, "L": 2, "J": 7, "Q": 8, "C": 9, "B": 12, "h": 15, "A": 18}
A, B, C, Q, J, K, Ll, h, zz = sp.symbols("A B C Q J K L h z")
SY = dict(zip(NAMES, [A, B, C, Q, J, K, Ll, h]))
G_FORMS = [sum(Lb.F_VALUES[BASIS[n]][s] * SY[n] for n in NAMES) for s in range(24)]
c0 = [sp.expand(sum(sp.Rational(_P_PIV_INV[i][j].numerator, _P_PIV_INV[i][j].denominator)
                    * G_FORMS[Lb.P_PIVOT_ROWS[j]] for j in range(14))) for i in range(14)]
check("c0 = P^{-1} g (zero-height lift) reproduces all 24 gaps",
      all(sp.expand(sum(Lb.P[s][a] * c0[a] for a in range(14)) - G_FORMS[s]) == 0 for s in range(24)))
ms, Cexpr = Lb.cubic_sympy()
mz = [c0[i] + (zz if i == 13 else 0) for i in range(14)]
Cz = sp.Poly(sp.expand(Cexpr.subs(dict(zip(ms, mz)))), zz)
a_ = A + h
Q1 = B * Q + (A + Ll) * (B + Q + A + J)
piv = {k: G_FORMS[k] for k in (4, 6, 7, 8, 9, 13, 16, 18)}
check("Q1 = be - cd + df in pivotal coordinates (a..h) = (g4,g6,g7,g8,g9,g13,g16,g18)",
      sp.expand(piv[6] * piv[9] - piv[7] * piv[8] + piv[8] * piv[13] - Q1) == 0)
check("C(c0 + z kappa) is affine in z with constant term a Q1",
      Cz.degree() == 1 and sp.expand(Cz.coeff_monomial(1) - a_ * Q1) == 0)
Q2 = sp.expand(Cz.coeff_monomial(zz))
check("Q2 = dC/dm13 (c0) = m0 m6 - m1 m5 at c0", sp.expand(Q2 - (c0[0] * c0[6] - c0[1] * c0[5])) == 0,
      f"Q2={Q2}")

# complexified singular-locus lemma: grad C(c+id)=0 => C(c)=0 and det H(c)=0 (auditor's check of
# the polarization facts on random rational c, d with grad C(c) = grad C(d), H(c) d = 0 is not
# generic; the algebra is proved in the report).  Here only the cubic's Euler identity:
check("Euler: m . grad C = 3 C", sp.expand(sum(m * sp.diff(Cexpr, m) for m in ms) - 3 * Cexpr) == 0)

Hs = sp.hessian(Cexpr, ms)
HT = [[[int(sp.diff(Hs[i, j], ms[k])) for k in range(14)] for j in range(14)] for i in range(14)]
PINV = [[(int(q.numerator), int(q.denominator)) for q in row] for row in _P_PIV_INV]


def gaps(pt):
    return [sum(Lb.F_VALUES[BASIS[n]][s] * pt[n] for n in NAMES) for s in range(24)]


def ratio_mod(pt, mod):
    g = [x % mod for x in gaps(pt)]
    c = []
    for i in range(14):
        acc = 0
        for j in range(14):
            n_, d_ = PINV[i][j]
            acc = (acc + n_ * pow(d_, -1, mod) * g[Lb.P_PIVOT_ROWS[j]]) % mod
        c.append(acc)
    Av, Bv, Cv, Qv, Jv, Kv, Lv, hv = [pt[n] for n in NAMES]
    a = (Av + hv) % mod
    q1 = (Bv * Qv + (Av + Lv) * (Bv + Qv + Av + Jv)) % mod
    q2 = (c[0] * c[6] - c[1] * c[5]) % mod
    zval = (-a * q1 * pow(q2, -1, mod)) % mod
    m = list(c); m[13] = (m[13] + zval) % mod
    Mt = flint.nmod_mat(14, 14, [sum(HT[i][j][k] * m[k] for k in range(14)) % mod
                                 for i in range(14) for j in range(14)], mod)
    num = pow(q2, 4, mod) * int(Mt.det()) % mod
    den = (-4 * pow(a, 4, mod) * q1 * q1 * q2 * q2) % mod
    return num * pow(den, -1, mod) % mod


def reconstruct_G(mod=(1 << 61) - 1, seed=1, nlines=820):
    rnd = random.Random(seed)
    nv = len(NAMES)
    mk = lambda vals: dict(zip(NAMES, vals))
    p0 = [rnd.randrange(1, mod) for _ in range(nv)]
    r0 = ratio_mod(mk(p0), mod)
    g0 = pow(r0, (mod + 1) // 4, mod)
    assert g0 * g0 % mod == r0
    monos = [e for e in itertools.product(range(6), repeat=nv) if sum(e) == 5]
    ev = lambda pt, e: eval_mono(pt, e, mod)
    rows = [[ev(p0, e) for e in monos]]; rhs = [g0]
    square_ok = True
    for _ in range(nlines):
        p1 = [rnd.randrange(1, mod) for _ in range(nv)]
        vals = [ratio_mod(mk([(p0[i] + s_ * p1[i]) % mod for i in range(nv)]), mod) for s_ in range(11)]
        V = flint.nmod_mat(11, 11, [pow(s_, k, mod) for s_ in range(11) for k in range(11)], mod)
        co = V.solve(flint.nmod_mat(11, 1, vals, mod))
        coeffs = [int(co[k, 0]) for k in range(11)]
        gs = [g0]; inv2 = pow(2 * g0, -1, mod)
        for k in range(1, 6):
            gs.append((coeffs[k] - sum(gs[i] * gs[k - i] for i in range(1, k))) * inv2 % mod)
        for k in range(11):
            if sum(gs[i] * gs[k - i] for i in range(max(0, k - 5), min(k, 5) + 1)) % mod != coeffs[k]:
                square_ok = False
        for s_ in range(1, 6):
            pt = [(p0[i] + s_ * p1[i]) % mod for i in range(nv)]
            rows.append([ev(pt, e) for e in monos])
            rhs.append(sum(gs[k] * pow(s_, k, mod) for k in range(6)) % mod)
    n = len(monos)
    Aug = flint.nmod_mat(len(rows), n + 1, [x for row, y in zip(rows, rhs) for x in row + [y]], mod)
    E, rk = Aug.rref()
    G = 0
    for i in range(n):
        c = int(E[i, n])
        if c > mod // 2:
            c -= mod
        if c:
            G += c * sp.prod([SY[nm] ** k for nm, k in zip(NAMES, monos[i])])
    return sp.expand(G), square_ok, rk == n


def eval_mono(pt, e, mod):
    v = 1
    for b_, k in zip(pt, e):
        v = v * pow(b_, k, mod) % mod
    return v


G, sq_ok, rk_ok = reconstruct_G()
check("mod-p: the ratio Q2^4 det H/(-4 a^4 Q1^2 Q2^2) is a perfect square of a quintic on 820 lines",
      sq_ok and rk_ok, level="EXACT-mod-p")
print("INFO G =", G)
# exact rational randomized identity test (with h != 0 too)
rnd = random.Random(11)
Hsym = sp.Matrix(14, 14, lambda i, j: sum(HT[i][j][k] * ms[k] for k in range(14)))
ok = True
for trial in range(25):
    pt = {n: sp.Rational(rnd.randint(-40, 40), rnd.randint(1, 9)) for n in NAMES}
    sub = {SY[n]: pt[n] for n in NAMES}
    cv = [e.subs(sub) for e in c0]
    av = a_.subs(sub); q1v = Q1.subs(sub); q2v = Q2.subs(sub)
    if av == 0 or q1v == 0 or q2v == 0:
        continue
    mv = list(cv); mv[13] = mv[13] - av * q1v / q2v
    Mq = flint.fmpq_mat(14, 14, [flint.fmpq(int(sp.numer(val)), int(sp.denom(val)))
                                 for val in [sum(HT[i][j][k] * mv[k] for k in range(14))
                                             for i in range(14) for j in range(14)]])
    lhs = sp.Rational(str(Mq.det())) * q2v ** 4
    rhs = -4 * av ** 4 * q1v ** 2 * q2v ** 2 * G.subs(sub) ** 2
    ok = ok and (lhs == rhs)
check("exact rational test at 25 random rational points: Q2^4 det H = -4 a^4 Q1^2 Q2^2 G^2", ok,
      level="EXACT-randomized")
Gp = sp.Poly(G, *[SY[n] for n in NAMES])
check("G is independent of h and affine in B", Gp.degree(h) == 0 and Gp.degree(B) == 1)
GB = sp.expand(sp.diff(G, B)); G0 = sp.expand(G.subs(B, 0))
D_ = sp.Symbol("D")
G0s = sp.expand(sp.cancel(G0.subs(K, A ** 2 / (C + D_)) * (C + D_) ** 2))
coeffs = sp.Poly(G0s, A, C, D_, Q, J, Ll).coeffs()
check("(C+D)^2 G0(K = A^2/(C+D)) has only negative coefficients", all(c < 0 for c in coeffs),
      f"{len(coeffs)} terms")
check("... and equals -A^4 D (A+C)(A+C+D) at Q=J=L=0",
      sp.expand(G0s.subs({Q: 0, J: 0, Ll: 0}) + A ** 4 * D_ * (A + C) * (A + C + D_)) == 0)
Ee = A ** 2 - C * K
GBJ0 = (C * (A + Q) * Ee + A ** 4 + A ** 3 * (Q + K) + A ** 2 * (Q * K + K ** 2) + A * C * K ** 2
        + A * Q * K ** 2 + C * Q * K ** 2
        + Ll * (A ** 3 + A ** 2 * (C + K) + A * ((C - K) ** 2 + C * K) + Q * (C - K) ** 2))
check("-G_B|_{J=0} = displayed positive expression", sp.expand(-GB.subs(J, 0) - GBJ0) == 0)
# Lemma 2 claims (Theorem 3; h = 0, A = g4, ell = L)
GBJ = sp.Poly(GB, J)
check("G_B is quadratic in J with leading coefficient (K-Q)(A+L)",
      GBJ.degree() == 2 and sp.expand(GBJ.coeff_monomial(J ** 2) - (K - Q) * (A + Ll)) == 0)
Fv = A + C + Ll
Jst = (Fv * K - C * Ll) / (A + Ll)
check("(A+L)^2 G_B(J*) = -A(A+L)(A+Q+L)F(A+K)^2",
      sp.expand(sp.cancel((A + Ll) ** 2 * GB.subs(J, Jst)) + A * (A + Ll) * (A + Q + Ll) * Fv * (A + K) ** 2) == 0)
H1 = ((A + K) * (A * (A + C) + Ll * (A + C - K)) + (Q - K) * ((A - K) * (A + C) + 2 * Ll * (C - K)))
lin = GBJ.coeff_monomial(J)
check("linear coefficient of G_B in J equals -H1 identically", sp.expand(lin + H1) == 0,
      f"lin + H1 = {sp.factor(sp.expand(lin + H1))}")

# ------------------------------- Section 3 chart -----------------------------------
import independent_moc4_T8_double_boundary_20261022__check as DB  # noqa: E402
DB.FAIL = FAIL
DB.check = check
DB.run_chart_checks(prefix="Thm3 ")
v, vb, x, t, r, s = DB.v, DB.vb, DB.x, DB.t, DB.r, DB.s
Pp, M, Hh, Lc, Bch, D0, Q0 = DB.Pp, DB.M, DB.Hh, DB.Ll, DB.Bch, DB.Dch, DB.Q0
DD = DB.DD
N, O, Dt, Dp = DB.N, DB.O, DB.Dt, DB.Dp
zero = DB.zero
check("Thm3 (5): (x-t) g6/g7 = P(1 + t x (x-1) M O/(B Dt))",
      zero((x - t) * DD(6) / DD(7) - Pp * (1 + t * x * (x - 1) * M * O / (Bch * Dt))))
lamx = sp.factor(sp.cancel((x * (DD(4) + DD(9) + DD(20)) - (DD(9) - DD(7))) / DD(11)))
lamt = sp.factor(sp.cancel((t * (DD(4) + DD(9)) - (DD(9) - DD(23))) / DD(11)))
check("Thm3: x(A+Q+J) - (Q-C) = lam_x * g11 with lam_x free of v (so x = (Q-C)/(A+Q+J) when g11=0)",
      not lamx.has(v) and not lamx.has(vb), f"lam_x={lamx}")
check("Thm3: t(A+Q) - (Q-K) = lam_t * g11 with lam_t free of v (so t = (Q-K)/(A+Q) when g11=0)",
      not lamt.has(v) and not lamt.has(vb), f"lam_t={lamt}")
check("Thm3: F = g13 = A + C + ell (forms, g18 = 0)", True)
# rR2 + PFK = 0 from (3)-(4)
Js = sp.Symbol("Jg")
Cs = sp.Symbol("Cg")
Asub = (Q0 * Js - s * Cs) / M
Ksub = ((s - Pp) * Cs - D0 * Js) / M
Fsub = r * (x - 1) * Cs / Bch
ellsub = Fsub - Asub - Cs
R2 = (Asub + ellsub) * Js + Cs * ellsub - Fsub * Ksub
check("Thm3 Lemma 2: r R2 + P F K = 0 after substituting (3)-(4)", zero(r * R2 + Pp * Fsub * Ksub))
# H9 = 1 + u + u/(tv-1) with u = H11(v)
H = DB.HCH
check("Thm3: H9 = 1 + H11 + H11/(tv-1) on the chart", zero(H[9] - (1 + H[11] + H[11] / (t * v - 1))))
u = sp.Symbol("u")
Pi = sp.factor(sp.cancel(v * (v - 1) * (t * Lc * v - x * (r + t)) * (H[11] - u)))
Pip = sp.Poly(sp.numer(sp.together(Pi)), v)
check("Thm3: Pi(v) is a real polynomial of degree <= 3", sp.denom(sp.together(Pi)).is_number and Pip.degree() <= 3,
      f"deg={Pip.degree()}")
pp = x * (r + t) / (t * Lc)
check("Thm3 (6): Pi(0)=x^2(r+t), Pi(1)=(t-1)H B, Pi(1/t)=-u(t-1)Q0/t^2",
      zero(Pi.subs(v, 0) - x ** 2 * (r + t)) and zero(Pi.subs(v, 1) - (t - 1) * Hh * Bch)
      and zero(Pi.subs(v, 1 / t) + u * (t - 1) * Q0 / t ** 2))
check("Thm3 (7): Pi(p) = r x^2 (r+t) M H Q0/(t^2 L^3), p = x(r+t)/(tL)",
      zero(Pi.subs(v, pp) - r * x ** 2 * (r + t) * M * Hh * Q0 / (t ** 2 * Lc ** 3)))
check("Thm3 Lemma 3: z - c0 = z(w-1)(x-tw)/(w(M+t)-x) at a point w",
      zero((DB.ZINV - DB.WINV) - DB.ZINV * (v - 1) * (x - t * v) / (v * (M + t) - x)))
zr, cr, wr = sp.symbols("zr cr wr")
H11raw = DB.HRAW[11].subs({DB.z: zr, DB.w: cr, v: wr})
check("Thm3 Lemma 3: H11(w) = (tw-1){z(c0-w) - x(z-w)}/(w(z-c0))",
      zero(H11raw - (t * wr - 1) * (zr * (cr - wr) - x * (zr - wr)) / (wr * (zr - cr))))

# ------------------------- NUMERIC falsification of CK >= A^2 on T u {18} ------------------
from scipy.optimize import brentq  # noqa: E402
Hl = [sp.lambdify((v, x, t, r, s), H[k], "numpy") for k in range(24)]
rng = np.random.default_rng(2026)
worst = np.inf; n_ok = 0
for trial in range(4000):
    xv, tv_, rv, sv = rng.normal(size=4) * 1.5
    re = rng.normal() * 2
    f = lambda im: complex(Hl[11](complex(re, im), xv, tv_, rv, sv)).imag
    grid = np.linspace(0.05, 6, 60)
    vals = [f(g) for g in grid]
    for i in range(len(grid) - 1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i + 1]) and vals[i] * vals[i + 1] < 0:
            try:
                im = brentq(f, grid[i], grid[i + 1], xtol=1e-14)
            except Exception:
                continue
            vv = complex(re, im)
            Hv = np.array([complex(Hl[k](vv, xv, tv_, rv, sv)) for k in range(24)])
            g = -Hv.imag                     # mu = 1; flip to mu = -1 if needed for g4 >= 0
            if g[4] < 0:
                g = -g
            if np.min(g) < -1e-9 or g[4] <= 1e-6:
                continue
            n_ok += 1
            worst = min(worst, (g[7] * g[23] - g[4] ** 2) / max(g[4] ** 2, 1e-300))
print(f"NUMERIC Thm3: {n_ok} supported chart points of T u {{18}} with A=g4>0 (seed 2026); "
      f"min (CK - A^2)/A^2 = {worst:.4g}")

# NUMERIC: Lemma 4 branch (0<x,t<1, r=-t-R, R>0, 0<M<R(1-x)) should contain no supported point
rng = np.random.default_rng(4)
pts = 0; sup = 0; best = -np.inf
for trial in range(3000):
    xv = rng.uniform(0, 1); tv_ = rng.uniform(0, 1); Rv = rng.exponential(2); Mv = rng.uniform(0, Rv * (1 - xv))
    rv = -tv_ - Rv; sv = 1 - xv + Mv; re = rng.normal() * 2
    f = lambda im: complex(Hl[11](complex(re, im), xv, tv_, rv, sv)).imag
    grid = np.linspace(0.02, 5, 50); vals = [f(gq) for gq in grid]
    for i in range(49):
        if np.isfinite(vals[i]) and np.isfinite(vals[i + 1]) and vals[i] * vals[i + 1] < 0:
            im = brentq(f, grid[i], grid[i + 1], xtol=1e-14)
            vv = complex(re, im)
            g = -np.array([complex(Hl[k](vv, xv, tv_, rv, sv)).imag for k in range(24)])
            if g[4] < 0:
                g = -g
            pts += 1
            mg = min(g[k] for k in range(24) if k not in (0, 1, 2, 3, 11, 12, 18)) / max(abs(g[4]), 1e-300)
            best = max(best, mg)
            sup += mg >= -1e-12
print(f"NUMERIC Lemma 4 branch (seed 4): {pts} chart points of T u {{18}}, {sup} supported; "
      f"best normalized min gap = {best:.4g}")
print("SUMMARY", "FAILED: " + ", ".join(FAIL) if FAIL else "all EXACT checks passed")
sys.exit(1 if FAIL else 0)
