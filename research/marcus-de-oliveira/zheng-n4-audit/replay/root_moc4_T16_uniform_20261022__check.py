"""Exact replay and brief numeric tests for the audit of root_moc4_T16_uniform_20261022
(supplement label T16 = paper case T_17 = {0,1,4,5,9,12}).
Own code only; nothing from the supplement is executed, imported or read.
Run from this directory:  python -I -B root_moc4_T16_uniform_20261022__check.py [--no-numeric]
Prints one line per check; exits nonzero if an EXACT check fails."""
import os, sys, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lib'))
import moc4lib as L
import sympy as sp
from fractions import Fraction as Fr

FAIL = []
def chk(name, ok, info=''):
    print(('PASS ' if ok else 'FAIL ') + name + ((' : ' + str(info)) if info != '' else ''))
    if not ok:
        FAIL.append(name)

T = {0, 1, 4, 5, 9, 12}
NU = 17
# historical supplement labels cited in the document's inputs -> paper cases (via the
# supplement's case index, which lists the paper contact set for each historical label;
# here we only use the paper contact sets themselves)
CITED_PAPER_CASES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 16]
STAGE = {1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1, 15: 1, 10: 2, 9: 3, 11: 3, 12: 3, 8: 4, 13: 4, 17: 4, 16: 5}

# ---------------------------------------------------------------- combinatorics
def stage_combinatorics():
    gens = L.generators(T)
    chk('contact set equals paper T_17 (Appendix E)', L.CONTACT_SETS[NU] == frozenset(T))
    chk('generators F_T = [0,1,2,5,7,9,12,35]', gens == [0, 1, 2, 5, 7, 9, 12, 35], gens)
    chk('only f35 is nonminor', [j for j in gens if not L.F_IS_MINOR[j]] == [35], L.f_name(35))
    chk('8 generators of rank 8 (simplicial cone); dim(im P cap {g_T=0}) = 8',
        L.rank([list(L.F_VALUES[j]) for j in gens]) == 8 and 14 - L.rank([L.P[s] for s in sorted(T)]) == 8)
    # the document's expansion: coefficient functionals (dict index -> coeff) for each generator
    coef = {0: {18: 1}, 1: {6: 1, 7: -1, 8: -1, 16: 1}, 2: {2: -1, 7: -1, 13: 1},
            5: {6: 1, 7: -1, 8: -1, 10: 1}, 7: {2: -1, 6: -1, 8: 1}, 9: {7: 1}, 12: {6: 1}, 35: {2: 1}}
    ok = True
    for m in gens:
        for j in gens:
            val = sum(c * L.F_VALUES[m][s] for s, c in coef[j].items())
            if val != (1 if j == m else 0):
                ok = False
    chk('expansion g = sum coef_j(g) f_j holds on the 8-dim span (coef forms are the dual basis)', ok)
    zer = sorted(L.F_ZERO[35] - T)
    chk('nontied zeros of f35 are {6,7,11,14,15,17,18}', zer == [6, 7, 11, 14, 15, 17, 18], zer)
    # reductions of each extra contact j to an earlier case
    allok = True
    for j in zer:
        S = T | {j}
        found = []
        for nu in CITED_PAPER_CASES:
            for sub in itertools.combinations(sorted(S), 6):
                g = L.find_g(L.CONTACT_SETS[nu], set(sub))
                if g is not None:
                    found.append((nu, sub))
                    break
        early = [f for f in found if STAGE.get(f[0], 99) <= 3]
        print('INFO extra contact %d: cited cases embedding in T u {j}: %s' % (j, found))
        if not early:
            allok = False
    chk('every nontied zero j of f35: T u {j} contains a G-image of a cited case of stage <= 3', allok)
    E = {s for s in range(24) if L.PERMS[s][0] == 0}
    gE = L.generators(E | T)
    chk('real-v layer: entry face {sigma(1)=1} u T has only minor generators', all(L.F_IS_MINOR[j] for j in gE), gE)

# ---------------------------------------------------------------- algebra
x, h, r, u, v, vb, z, w, t, s = sp.symbols('x h r u v vb z w t s')
INF = None

def Fprod(A, C, sig_idx):
    sig = L.PERMS[sig_idx]
    out = sp.Integer(1)
    for i in range(4):
        if C[sig[i]] is INF:
            out *= -1            # homogeneous bracket [(a,1),(1,0)] = -1
            continue
        out *= (A[i] - C[sig[i]])
    return out

def is_zero(e):
    return sp.cancel(sp.together(e)) == 0

def stage_algebra():
    A = [z, sp.Integer(1), x, sp.Integer(0)]
    C = [w, v, INF, t * v]
    F = [sp.expand(Fprod(A, C, j)) for j in range(24)]
    R = [sp.cancel((F[j] - F[0]) / (v * (z - w))) for j in range(24)]
    chk('R1 = x(v-1)/v', is_zero(R[1] - x * (v - 1) / v))
    chk('R4 = R1 + (t-1)(x-1)', is_zero(R[4] - R[1] - (t - 1) * (x - 1)))
    chk('R5 = 1 - t', is_zero(R[5] - (1 - t)))
    chk('(R4-R1)/R5 = 1 - x', is_zero((R[4] - R[1]) / R[5] - (1 - x)))
    chk('R1/R5 = x/hh with hh = v(1-t)/(v-1)', is_zero(R[1] / R[5] - x / (v * (1 - t) / (v - 1))))
    tt = 1 - h * (v - 1) / v
    chk('t = 1 - h(v-1)/v  <=>  tv = (1-h)v + h', is_zero(tt * v - ((1 - h) * v + h)))
    H = [sp.cancel(R[j].subs(t, tt) / (1 - tt)) for j in range(24)]
    chk('H1 = x/h, H4 = x/h - x + 1, H5 = 1',
        is_zero(H[1] - x / h) and is_zero(H[4] - (x / h - x + 1)) and is_zero(H[5] - 1))
    chk('(2) H2 = (x-1)(v/(h(v-1)) - 1)', is_zero(H[2] - (x - 1) * (v / (h * (v - 1)) - 1)))
    Imf = lambda f: (f - f.subs(v, vb)) / (2 * sp.I)
    chk('(2) d = g2 = -mu Im H2 = mu (x-1) Im v / (h |v-1|^2)',
        is_zero(-Imf(H[2]) - (x - 1) * ((v - vb) / (2 * sp.I)) / (h * (v - 1) * (vb - 1))))
    a = tt * (v - 1) - u * (1 - tt)
    D = tt * (x - 1) - u * (1 - tt)
    b = r * (1 - tt)
    E9 = w * (v * (b + tt - x) + z * (x - tt * v)) + z * v * (-b + tt * v - tt)
    E12 = a * z - D * w - tt * (v - x)
    n9 = sp.numer(sp.together(H[9].subs(t, tt) - r)) if False else sp.numer(sp.together(H[9] - r))
    n12 = sp.numer(sp.together(H[12] - u))
    q9 = sp.factor(sp.cancel(n9 / sp.numer(sp.together(E9))))
    q12 = sp.factor(sp.cancel(n12 / sp.numer(sp.together(E12))))
    ok9 = not sp.together(q9).has(z) and not sp.together(q9).has(w)
    ok12 = not sp.together(q12).has(z) and not sp.together(q12).has(w)
    chk('H9 = r  <=>  E9 = 0 (numerators proportional, factor free of z,w)', ok9, q9)
    chk('H12 = u  <=>  E12 = 0 (numerators proportional, factor free of z,w)', ok12, q12)
    chk('a = -(v-1)(hu + hv - h - v)/v', is_zero(a + (v - 1) * (h * u + h * v - h - v) / v))
    chk('D = (1-h)(x-1) - hu + h(x-1+u)/v', is_zero(D - ((1 - h) * (x - 1) - h * u + h * (x - 1 + u) / v)))
    chk('x - tv = x - h - (1-h)v', is_zero(x - tt * v - (x - h - (1 - h) * v)))
    zsol = (D * w + tt * (v - x)) / a
    chk('E12 solved: z = (Dw + t(v-x))/a', is_zero(E12.subs(z, zsol)))
    curve = D * (x - tt * v) * w**2 + (a * v * (b + tt - x) + tt * (v - x) * (x - tt * v) + D * v * (-b + tt * v - tt)) * w \
        + tt * (v - x) * v * (-b + tt * v - tt)
    chk('(3) a * E9(z = zsol) = curve (3)', is_zero(a * E9.subs(z, zsol) - curve))
    chk('z = w divisor: E12 = t(v-x)(z-1), E9 = z(v-z)(tv-x)',
        is_zero(E12.subs(w, z) - tt * (v - x) * (z - 1)) and is_zero(E9.subs(w, z) - z * (v - z) * (tt * v - x)))
    chk('F0/(v(z-w)(1-t)) = ((h-1)v - h)/h', is_zero(F[0].subs(t, tt) / (v * (z - w) * (1 - tt)) - ((h - 1) * v - h) / h))
    # D's reciprocal coefficient: D = d0 + d1/v with d1 = h(x-1+u); if d1 = 0 then d0 = x-1
    chk('D with u = 1-x is the constant x-1', is_zero(D.subs(u, 1 - x) - (x - 1)))
    # ---------------- fibre v = 1
    vs = 1 + s
    sub = {v: vs}
    curve_s = sp.factor(sp.together(curve.subs(sub)))
    cnum = sp.numer(curve_s); cden = sp.denom(curve_s)
    P0 = sp.factor(sp.together(curve.subs(v, 1)))
    chk('curve at v = 1 is (x-1)^2 w (w-1)', is_zero(P0 - (x - 1)**2 * w * (w - 1)), P0)
    a1 = 1 - h * u
    chk('a = a1 (v-1) + O((v-1)^2), a1 = 1 - hu', is_zero(sp.diff(a, v).subs(v, 1) - a1) and is_zero(a.subs(v, 1)))
    stage_branches(F)


# ---------------------------------------------------------------- truncated Laurent series
# Exact series in s = v - 1 with coefficients in the rational function field QQ(x,h,r,u).
from sympy import QQ
from sympy.polys.fields import field as _field
KF, KX, KH, KR, KU = _field("x,h,r,u", QQ)

class Ser:
    """Truncated Laurent series: coefficient dict {exponent: K element}, valid for exponents < prec."""
    def __init__(self, c, prec):
        self.prec = prec
        self.c = {e: a for e, a in c.items() if e < prec and a != 0}
    @staticmethod
    def const(a, prec):
        return Ser({0: KF(a) if not hasattr(a, 'numer') else a}, prec)
    def __add__(self, o):
        o = o if isinstance(o, Ser) else Ser.const(o, 10**6)
        p = min(self.prec, o.prec); c = dict(self.c)
        for e, a in o.c.items():
            c[e] = c.get(e, 0) + a
        return Ser(c, p)
    __radd__ = __add__
    def __neg__(self):
        return Ser({e: -a for e, a in self.c.items()}, self.prec)
    def __sub__(self, o):
        return self + (-(o if isinstance(o, Ser) else Ser.const(o, 10**6)))
    def __rsub__(self, o):
        return (-self) + o
    def val(self):
        return min(self.c) if self.c else self.prec
    def __mul__(self, o):
        o = o if isinstance(o, Ser) else Ser.const(o, 10**6)
        p = min(self.prec + o.val(), o.prec + self.val())
        c = {}
        for e1, a1 in self.c.items():
            for e2, a2 in o.c.items():
                if e1 + e2 < p:
                    c[e1 + e2] = c.get(e1 + e2, 0) + a1 * a2
        return Ser(c, p)
    __rmul__ = __mul__
    def inv(self):
        v0 = self.val(); c0 = self.c[v0]
        rel = self.prec - v0                       # relative precision
        e = Ser({k - v0: a / c0 for k, a in self.c.items() if k != v0}, rel)
        out = Ser({0: KF(1)}, rel); term = Ser({0: KF(1)}, rel)
        for _ in range(rel):
            term = term * (-e)
            out = out + term
        return Ser({k - v0: a / c0 for k, a in out.c.items()}, rel - v0)
    def __truediv__(self, o):
        o = o if isinstance(o, Ser) else Ser.const(o, 10**6)
        return self * o.inv()
    def coeff(self, e):
        assert e < self.prec, ('precision', e, self.prec)
        return self.c.get(e, KF(0))

def poly_to_ser(expr, gens, sers, prec):
    P = sp.Poly(sp.expand(expr), *gens, x, h, r, u)
    out = Ser({}, prec)
    for mon, cf in P.terms():
        term = Ser({0: KF(int(cf)) * KX**mon[-4] * KH**mon[-3] * KR**mon[-2] * KU**mon[-1]}, prec)
        for g_, e_ in zip(sers, mon[:-4]):
            for _ in range(e_):
                term = term * g_
        out = out + term
    return out

def stage_branches(F, M=6):
    S = Ser({1: KF(1)}, M)                    # s
    one = Ser({0: KF(1)}, M)
    V = one + S
    Vinv = V.inv()
    Tt = one - KH * S * Vinv
    A = Tt * (V - 1) - KU * (one - Tt)
    Dd = Tt * (KX - 1) - KU * (one - Tt)
    B = KR * (one - Tt)
    c2 = Dd * (KX - Tt * V)
    c1 = A * V * (B + Tt - KX) + Tt * (V - KX) * (KX - Tt * V) + Dd * V * (-B + Tt * V - Tt)
    c0 = Tt * (V - KX) * V * (-B + Tt * V - Tt)
    chk('series: a = (1-hu) s + O(s^2)', A.val() == 1 and A.coeff(1) == 1 - KH * KU)
    ok_all = True
    for w0, partner in ((0, 7), (1, 14)):
        W = Ser({0: KF(w0)}, M)
        for _ in range(3):                    # Newton iteration (precision 1 -> 2 -> 4 -> M), exact truncated
            Pw = c2 * W * W + c1 * W + c0
            dP = 2 * c2 * W + c1
            W = W - Pw / dP
            W = Ser(W.c, M)
        resid = c2 * W * W + c1 * W + c0
        chk('branch w=%d: Newton series solves curve (3) to O(s^%d)' % (w0, resid.prec), resid.val() >= resid.prec, resid.prec)
        Nn = Dd * W + Tt * (V - KX)
        Z = Nn / A
        res = {}
        for j in (2, partner):
            Fj = sp.expand(F[j] - F[0])
            num = poly_to_ser(Fj, (z, w, v, t), (Z, W, V, Tt), M)
            Hj = num / (V * (Z - W) * (one - Tt))
            if Hj.val() < -1:
                ok_all = False
            res[j] = Hj
            print('INFO branch w=%d: H%d has valuation %d, residue %s' % (w0, j, Hj.val(), Hj.coeff(-1)))
        chk('branch w=%d: res H2 = (x-1)/h and res H%d = -(x-1)/h' % (w0, partner),
            res[2].coeff(-1) == (KX - 1) / KH and res[partner].coeff(-1) == -(KX - 1) / KH)
        tot = res[2] + res[partner]
        chk('branch w=%d: H2 + H%d is regular (real-coefficient series), value %s' % (w0, partner, tot.coeff(0)),
            tot.val() >= 0)
        if w0 == 1:
            chk('branch w=1: w = 1 + (1-hu)/(x-1) s^2 + O(s^3)', W.coeff(1) == 0 and W.coeff(2) == (1 - KH * KU) / (KX - 1))
            chk('branch w=1: z = 1 + s + O(s^2) (uses a1 != 0)', Z.coeff(0) == 1 and Z.coeff(1) == 1)
        else:
            chk('branch w=0: w = O(s), inverse numerator -> 1 - x, z has a simple pole',
                W.coeff(0) == 0 and Nn.coeff(0) == 1 - KX and Z.val() == -1)
    chk('no pole of order >= 2 for H2 and its compensator on either branch', ok_all)

# ---------------------------------------------------------------- numeric
def stage_numeric(seed=20261125, N=3000, nphi=12):
    import numpy as np
    rng = np.random.default_rng(seed)
    Pn = L.P_numpy()
    def F_num(A, C):
        out = []
        for j in range(24):
            sig = L.PERMS[j]; p = 1 + 0j
            for i in range(4):
                if C[sig[i]] is None: p *= -1; continue
                p *= (A[i] - C[sig[i]])
            out.append(p)
        return np.array(out)
    nsupp = 0; worst_coef = np.inf; phis = []
    for _ in range(N):
        X, Hh, Rr, U = rng.normal(scale=2, size=4)
        V = complex(rng.normal(scale=2), rng.normal(scale=2))
        tt = 1 - Hh * (V - 1) / V
        a = tt * (V - 1) - U * (1 - tt); D = tt * (X - 1) - U * (1 - tt); b = Rr * (1 - tt)
        c2 = D * (X - tt * V); c1 = a * V * (b + tt - X) + tt * (V - X) * (X - tt * V) + D * V * (-b + tt * V - tt)
        c0 = tt * (V - X) * V * (-b + tt * V - tt)
        for W in np.roots([c2, c1, c0]):
            Z = (D * W + tt * (V - X)) / a
            Fv = F_num([Z, 1, X, 0], [W, V, None, tt * V])
            for mus in (1, -1):
                kap = 1j * mus / (V * (Z - W) * (1 - tt))
                gg = np.real(kap * (Fv - Fv[0]))
                sc = np.max(np.abs(gg))
                if np.max(np.abs(gg[[0, 1, 4, 5, 9, 12]])) > 1e-8 * sc or np.min(gg) < -1e-9 * sc:
                    continue
                nsupp += 1
                coefs = [gg[18], gg[6] - gg[7] - gg[8] + gg[16], -gg[2] - gg[7] + gg[13],
                         gg[6] - gg[7] - gg[8] + gg[10], -gg[2] - gg[6] + gg[8], gg[7], gg[6], gg[2]]
                worst_coef = min(worst_coef, min(coefs) / sc)
                if len(phis) < nphi:
                    cvec = np.linalg.lstsq(Pn, gg / sc, rcond=None)[0]
                    m = L.min_phi_over_SO4(cvec, restarts=6, seed=len(phis))
                    phis.append(m if np.isscalar(m) else m[0])
    print('NUMERIC seed %d, %d chart samples (both roots): %d supported with exact six contacts; min coefficient/max|g| = %.3e'
          % (seed, N, nsupp, worst_coef))
    if phis:
        print('NUMERIC min over SO(4) of Phi (6 restarts each) at %d supported points: worst = %.3e' % (len(phis), min(phis)))

if __name__ == '__main__' and '--branches-only' in sys.argv:
    _A = [z, sp.Integer(1), x, sp.Integer(0)]; _C = [w, v, INF, t * v]
    stage_branches([sp.expand(Fprod(_A, _C, j)) for j in range(24)])
    print('EXACT FAILURES:', FAIL if FAIL else 'none'); sys.exit(1 if FAIL else 0)
if __name__ == '__main__':
    stage_combinatorics()
    stage_algebra()
    if '--no-numeric' not in sys.argv:
        stage_numeric()
    print('EXACT FAILURES:', FAIL if FAIL else 'none')
    sys.exit(1 if FAIL else 0)
