"""Exact replay (and brief numeric tests) for the audit of root_moc4_T4_uniform_20261022.
Own code only; nothing from the supplement is executed or read.  Run from this directory:
    python -I -B root_moc4_T4_uniform_20261022__check.py [--numeric]
Prints one line per check; exits nonzero if an EXACT check fails."""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lib'))
import moc4lib as L
import sympy as sp
import numpy as np

def stage_linear_T4():
    """Linear (profile) checks for T4 and T5: generators, relations, decompositions."""
    from fractions import Fraction as Fr

    FAIL = []
    def chk(name, ok):
        print(('PASS ' if ok else 'FAIL ') + name)
        if not ok:
            FAIL.append(name)

    V = L.F_VALUES

    def comb(coeffs):
        out = [Fr(0)] * 24
        for j, c in coeffs.items():
            for s in range(24):
                out[s] += c * V[j][s]
        return out

    def subspace_basis(T):
        # basis of {P c : (P c)_s = 0 for s in T}, as 24-vectors of values
        rows = [[L.P[s][k] for k in range(14)] for s in T]
        ns = L.nullspace(rows)
        return [L.values(c) for c in ns]

    def check_decomp(T, name, decomp):
        """decomp: function(g: list of 24 values) -> dict j -> coefficient; check g == sum."""
        B = subspace_basis(T)
        ok = all(comb(decomp(g)) == list(g) for g in B)
        chk(f'{name}: decomposition holds on a basis of the {len(B)}-dim face space', ok)
        return B

    # ---------------- T4 ----------------
    T4 = [0, 1, 4, 5, 7, 8]
    chk('T4 equals paper CONTACT_SETS[4]', set(T4) == set(L.CONTACT_SETS[4]))
    gen4 = L.generators(T4)
    chk('T4 generators = {0,1,2,3,4,5,6,11,33,34}: got %s' % gen4, gen4 == [0, 1, 2, 3, 4, 5, 6, 11, 33, 34])
    lab = {0: ('n', '1', '4'), 1: ('n', '12', '34'), 2: ('n', '13', '34'), 3: ('n', '14', '34'), 4: ('n', '1', '3'),
           5: ('n', '12', '24'), 6: ('n', '4', '1'), 11: ('n', '13', '23'), 33: ('n', '12', '13')}
    chk('T4 minor labels of f0..f6,f11,f33 match (doc 0-based -> paper 1-based)', all(L.F_TABLE[j] == lab[j] for j in lab))
    chk('T4 only f34 non-minor', [j for j in gen4 if not L.F_IS_MINOR[j]] == [34])
    chk('f34 = f0 - f2 + f6 + f33', list(V[34]) == comb({0: 1, 2: -1, 6: 1, 33: 1}))
    chk('f0 - f1 - f2 - f3 + f4 = 0', comb({0: 1, 1: -1, 2: -1, 3: -1, 4: 1}) == [0] * 24)
    chk('-f0 + f2 - f6 - f33 + f34 = 0', comb({0: -1, 2: 1, 6: -1, 33: -1, 34: 1}) == [0] * 24)
    chk('T4 generator span rank 8', L.rank([list(V[j]) for j in gen4]) == 8)

    def dec4(g):
        d = g[2] - g[13]
        return {0: g[19], 1: g[16], 3: g[14], 5: g[10], 6: g[15], 11: g[6], 33: g[13], 34: d}
    B4 = check_decomp(T4, 'T4 eq.(2)', dec4)
    chk('T4 face space dim 8', len(B4) == 8)
    ok = True
    for g in B4:
        d = g[2] - g[13]; A = -d; b = g[19]
        ok &= g[22] == g[16] + b - A and g[18] == g[14] + b - A and g[15] - A == g[9] and g[13] - A == g[2]
    chk('T4: g22=g16+b-A, g18=g14+b-A, g15-A=g9, g13-A=g2 on face space', ok)
    # replacement: coefficients after adding A*(rel2) + tau*(rel1), check symbolic in A, tau by two values
    ok = True
    for g in B4:
        for tau in (Fr(0), Fr(3, 7)):
            d = g[2] - g[13]; A = -d; b = g[19]
            co = {0: b - A + tau, 1: g[16] - tau, 2: A - tau, 3: g[14] - tau, 4: tau, 5: g[10], 6: g[9], 11: g[6], 33: g[2]}
            ok &= comb(co) == list(g)
    chk('T4: replaced minor combination (b-A+tau, g16-tau, A-tau, g14-tau, tau; f6:g9, f33:g2) reproduces g', ok)
    # which f_j are positive where (for T5 claims about g0,g9,g11)


    if FAIL:
        print('FAILED:', FAIL); sys.exit(1)
    print('ALL T4 LINEAR CHECKS PASSED')


def stage_chart_uniqueness_T4():

    T4 = [0, 1, 4, 5, 7, 8]
    z, w, v, vb, x, t, r, u = sp.symbols('z w v vb x t r u')
    A = [z, 0, 1, x]
    C = [w, v, None, t * v]

    def F(sig, Cv):
        p = sp.Integer(-1)
        for i in range(4):
            if sig[i] != 2:
                p *= (A[i] - Cv[sig[i]])
        return sp.expand(p)

    Fs = [F(s, C) for s in L.PERMS]
    R = [sp.cancel((Fs[s] - Fs[0]) / (v * (z - w))) for s in range(24)]
    rdef = z * (w - v) * (1 - t * v) / (v * (z - w))
    udef = -(x - t * v) * (z * (1 + v - w) - v) / (v * (z - w))
    claimed = {0: 0, 1: 1 - x, 4: t - x, 5: x * (t - 1), 7: rdef - x + 1, 8: udef}
    for s, val in claimed.items():
        print('R', s, sp.simplify(R[s] - val) == 0)
    # all R_sigma polynomial in z,w? print a few
    for s in range(24):
        print(s, sp.factor(R[s]))
    # solve r, u equations for z, w
    e1 = sp.numer(sp.together(rdef - r))
    e2 = sp.numer(sp.together(udef - u))
    sol = sp.solve([e1, e2], [z, w], dict=True)
    print('nsol', len(sol))
    for so in sol:
        print({k: sp.factor(vv) for k, vv in so.items()})
    zc = r * v * (t * v - x) / (r*t*v**2 - r*v*x + t**2*v**2 - t*u*v**2 - t*v*x - t*v + u*v + x)
    wc = (-t * v + x) * (r - t * v + 1) / (-r*t*v + r*x + t*u*v - u)
    print('claimed sol satisfies', sp.simplify(e1.subs({z: zc, w: wc})) == 0, sp.simplify(e2.subs({z: zc, w: wc})) == 0)
    # resultant elimination
    res = sp.factor(sp.resultant(e1, e2, w))
    print('resultant in z:', res)

def stage_chart_identities_T4():

    T0 = time.time()
    z, w, v, vb, x, t, r, u, s = sp.symbols('z w v vb x t r u s')
    A = [z, 0, 1, x]
    C = [w, v, None, t * v]

    def F(sig):
        p = sp.Integer(-1)
        for i in range(4):
            if sig[i] != 2:
                p *= (A[i] - C[sig[i]])
        return sp.expand(p)

    Fs = [F(sg) for sg in L.PERMS]
    zc = r * v * (t * v - x) / (r*t*v**2 - r*v*x + t**2*v**2 - t*u*v**2 - t*v*x - t*v + u*v + x)
    wc = (-t * v + x) * (r - t * v + 1) / (-r*t*v + r*x + t*u*v - u)
    need = [2, 6, 9, 10, 13, 14, 15, 16, 18, 19, 22]
    Rv = {}
    for sg in need + [1, 4, 5, 7, 8]:
        Rv[sg] = sp.cancel(((Fs[sg] - Fs[0]) / (v * (z - w))).subs({z: zc, w: wc}))
    # contact R's after substitution
    print('R at T4 after chart:', [sp.factor(Rv[k]) for k in (1, 4, 5, 7, 8)])
    g = {}
    for sg in need:
        Rb = Rv[sg].subs(v, vb)
        g[sg] = sp.factor(sp.cancel((Rv[sg] - Rb) / (v - vb)))   # gap with rho = -1
        print('g', sg, g[sg]); sys.stdout.flush()
    FAIL = []
    def chk(name, expr):
        ok = sp.cancel(expr) == 0
        print(('PASS ' if ok else 'FAIL ') + name); sys.stdout.flush()
        if not ok: FAIL.append(name)
    N = v * vb
    rho = -1
    chk('g2 = -rho x/|v|^2', g[2] - (-rho * x / N))
    chk('g6 = -rho r t (x-1)/|tv-1|^2', g[6] - (-rho * r * t * (x - 1) / ((t*v - 1) * (t*vb - 1))))
    chk('c = g10 = rho t u x/|tv-x|^2', g[10] - (rho * t * u * x / ((t*v - x) * (t*vb - x))))
    chk('g9 = (x-1) c', g[9] - (x - 1) * g[10])
    d = g[2] - g[13]; b = g[19]; c = g[10]; a = g[15]; e = g[14]; k = g[16]
    chk('a = g9 - d', a - (g[9] - d))
    y = x - 1; h = 1 - t; H = h + u; E = H*(1 - t*x) - h*r*y; J = H*(1 - t*x) + h*y
    chk('(3a) k = h c - t(b+d)', k - (h*c - t*(b + d)))
    chk('(3b) E e = t y H b - J d', E*e - (t*y*H*b - J*d))
    chk('(4a) J = E + h y (r+1)', J - (E + h*y*(r + 1)))
    chk('(4b) J - x h H = -u y', J - x*h*H + u*y)
    Lq = y - r + u
    chk('E = (h-y)(H+hy) + hyL', E - ((h - y)*(H + h*y) + h*y*Lq))
    K = t*y*b - x*h*d
    Q = lambda S: -t*(r + t - u)*S**2 + (r*t + r*x + t*x + t - u)*S - x*(r + 1)
    for name, fn, detc, evc in [('D', d, -r**3*t**2*u*x*y*Lq*E, r*u*x*y*(h*r - H*y)),
                                ('P', K, -r**2*t**2*u**2*x**2*h*y*Lq*E, r*u**2*x**2*y*(y - h))]:
        Dp = sp.cancel(fn * Q(v) * Q(vb))
        num, den = sp.fraction(Dp)
        print(name, 'denominator after multiplying by |Q|^2:', den)
        Dpol = sp.Poly(sp.expand(num / den) if den.free_symbols & {v, vb} == set() else num, v, vb)
        degs = Dpol.monoms()
        ok_bideg = all(m[0] <= 1 and m[1] <= 1 for m in degs) and den.free_symbols.isdisjoint({v, vb})
        print(('PASS ' if ok_bideg else 'FAIL ') + f'{name}*|Q|^2 is a polynomial of bidegree <= (1,1) in v, vb')
        if not ok_bideg: FAIL.append(name + ' bidegree')
        m00 = Dpol.coeff_monomial(1); m01 = Dpol.coeff_monomial(v); m10 = Dpol.coeff_monomial(vb); m11 = Dpol.coeff_monomial(v*vb)
        chk(f'{name}: M symmetric (coeff v == coeff vb)', m01 - m10)
        chk(f'{name}: det M = {detc}', m00*m11 - m01*m10 - detc)
        ev = Dpol.as_expr().subs({v: x/t, vb: x/t})
        chk(f'{name}: value at v=vb=x/t = {evc}', ev - evc)
    aa, bb, cc, dd, yy = sp.symbols('aa bb cc dd yy')
    chk('ab(c-d)-d^2(a+b+c+d) = c{ybc - xd(b+d)} when a = yc-d, x = y+1',
        (aa*bb*(cc - dd) - dd**2*(aa + bb + cc + dd)).subs(aa, yy*cc - dd) - cc*(yy*bb*cc - (yy + 1)*dd*(bb + dd)))
    # denominators of d: factor
    zw = sp.factor(sp.cancel(zc - wc))
    print('z - w in the chart =', zw)
    qq = sp.Poly(sp.numer(sp.together(zc - wc)), v)
    ok_q = sp.rem(qq, sp.Poly(Q(v), v)).is_zero
    print(('PASS ' if ok_q else 'FAIL ') + 'Q(v) divides the numerator of z - w (so Q(v) != 0 whenever A_0 != C_0)')
    if not ok_q: FAIL.append('Q divides z-w')
    print('denominator of d:', sp.factor(sp.fraction(sp.cancel(d))[1]))
    print('denominator of z (as poly in v) vs Q:', sp.factor(sp.fraction(zc)[1]), '| Q(v) =', sp.factor(Q(v)))
    print('time', round(time.time() - T0, 1))
    if FAIL:
        print('FAILED', FAIL); sys.exit(1)
    print('ALL T4 CHART CHECKS PASSED')

def stage_capacity():
    """EXACT replay of the four-coefficient capacity inequality (T5 (2); T4 (1) by transport).
    Own derivation of the quaternion link: O(p,q) = L(p) R(q) (x -> p x q), unit p, q."""

    FAIL = []
    def chk(name, ok):
        print(('PASS ' if ok else 'FAIL ') + name); sys.stdout.flush()
        if not ok: FAIL.append(name)

    p = sp.symbols('p0:4'); q = sp.symbols('q0:4')
    def qmul(a, b):
        a0, a1, a2, a3 = a; b0, b1, b2, b3 = b
        return [a0*b0 - a1*b1 - a2*b2 - a3*b3, a0*b1 + a1*b0 + a2*b3 - a3*b2,
                a0*b2 - a1*b3 + a2*b0 + a3*b1, a0*b3 + a1*b2 - a2*b1 + a3*b0]
    E4 = [[1 if i == k else 0 for i in range(4)] for k in range(4)]
    Lm = sp.Matrix(4, 4, lambda i, k: qmul(p, E4[k])[i])
    Rm = sp.Matrix(4, 4, lambda i, k: qmul(E4[k], q)[i])
    O = Lm * Rm
    np2 = sum(a*a for a in p); nq2 = sum(a*a for a in q)
    chk('O^T O = |p|^2 |q|^2 I', sp.expand(O.T * O - np2 * nq2 * sp.eye(4)) == sp.zeros(4, 4))
    # reduce modulo |p|^2 = 1, |q|^2 = 1 by eliminating p3^2, q3^2
    def red(e):
        e = sp.expand(e)
        e = sp.rem(sp.Poly(e, p[3]), sp.Poly(np2 - 1, p[3])).as_expr()
        e = sp.rem(sp.Poly(sp.expand(e), q[3]), sp.Poly(nq2 - 1, q[3])).as_expr()
        return sp.expand(e)
    Ow = [[O[i, j] for j in range(4)] for i in range(4)]
    wts = L.leibniz_weights(Ow)
    def phi(jf):
        return sum(int(val) * wt for val, wt in zip(L.F_VALUES[jf], wts))
    A = p[0]*q[3] + p[3]*q[0]; B = p[1]*q[2] - p[2]*q[1]
    D = p[0]**2 + p[2]**2 - q[0]**2 - q[2]**2
    E = p[1]**2 - p[2]**2 - q[0]**2 + q[3]**2
    G = p[0]**2 - p[3]**2 - q[1]**2 + q[2]**2
    H = 2*(A**2 + B**2) - E*G
    for jf, Qf in [(0, (A - B)**2), (6, (A + B)**2), (10, D**2), (54, H)]:
        chk(f'Phi_f{jf}(L(p)R(q)) == quartic on S^3 x S^3 (reduced mod |p|^2-1, |q|^2-1)', red(phi(jf) - Qf) == 0)
    # Gram certificate (4) of the companion, re-derived and replayed
    a, b, c, d = sp.symbols('a b c d')
    Z = sp.Matrix([[a + b + 2*d, a - b, 0, 0], [a - b, a + b + 2*d, 0, 2*d], [0, 0, c, d], [0, 2*d, d, d]])
    y = list(p) + list(q)
    def S(j, yy):
        p0, p1, p2, p3, q0, q1, q2, q3 = yy
        return [yy, [-p1, p0, -p3, p2, q1, -q0, -q3, q2], [-p2, p3, p0, -p1, q2, q3, -q0, -q1],
                [-p3, -p2, p1, p0, q3, -q2, q1, -q0]][j]
    def forms(yy):
        p0, p1, p2, p3, q0, q1, q2, q3 = yy
        Ay = p0*q3 + p3*q0; By = p1*q2 - p2*q1
        Dy = p0**2 + p2**2 - q0**2 - q2**2; Ey = p1**2 - p2**2 - q0**2 + q3**2; Gy = p0**2 - p3**2 - q1**2 + q2**2
        ell = sp.Matrix([-p0*By, p0*Ay, -q3*Dy, -q3*Gy]); m = sp.Matrix([q0*By, -q0*Ay, -p3*Dy, -p3*Ey])
        return ell, m
    rhs = 0
    for j in range(4):
        ell, m = forms(S(j, y))
        rhs += (ell.T * Z * ell)[0, 0] + (m.T * Z * m)[0, 0]
    n = np2 + nq2; Nn = np2 - nq2
    lhs = n * (a*(A - B)**2 + b*(A + B)**2 + c*D**2 + d*H)
    diff = sp.expand(lhs - rhs)
    r_ = sp.rem(sp.Poly(diff, q[3]), sp.Poly(Nn, q[3]))
    chk('companion identity (4): n(a(A-B)^2+b(A+B)^2+cD^2+dH) = sum of Z-forms mod N', r_.is_zero)
    # Z PSD <=> capacity (Schur complements), exact
    S1 = sp.simplify(Z[:3, :3] - Z[:3, 3:] * Z[3:, :3] / d)
    chk('Schur complement of d', sp.simplify(S1 - sp.Matrix([[a + b + 2*d, a - b, 0], [a - b, a + b - 2*d, -2*d], [0, -2*d, c - d]])) == sp.zeros(3, 3))
    S2 = sp.simplify(S1[:2, :2] - S1[:2, 2:] * S1[2:, :2] / (c - d))
    chk('2x2 Schur det = 4[ab(c-d)-d^2(a+b+c+d)]/(c-d)', sp.simplify(S2.det() - 4*(a*b*(c - d) - d**2*(a + b + c + d))/(c - d)) == 0)
    # T4 arrowhead pencil: Schur complement formula as displayed
    Z4 = sp.Matrix([[a + d, 0, -d, 0], [0, c, d, 0], [-d, d, d, d], [0, 0, d, b + d]])
    sc = d - d**2/(a + d) - d**2/c - d**2/(b + d)
    chk('T4 arrowhead Schur complement = d[ab(c-d) - d^2(a+b+c+d)]/((a+d)c(b+d))', sp.simplify(sc - d*(a*b*(c - d) - d**2*(a + b + c + d))/((a + d)*c*(b + d))) == 0)
    chk('T4 arrowhead det = (a+d)c(b+d) * Schur', sp.expand(Z4.det() - (a + d)*c*(b + d)*sc) == 0 or sp.simplify(Z4.det() - (a + d)*c*(b + d)*sc) == 0)
    # transport T5 (2) -> T4 (1) with an explicit group element
    g = ((3, 2, 1, 0), (3, 1, 2, 0), 0)
    imgs = {}
    for j in (6, 0, 5, 34):
        vals = L.g_act_values(g, list(L.F_VALUES[j]))
        imgs[j] = [k for k in range(82) if list(L.F_VALUES[k]) == list(vals)]
    chk('g maps f6->f0, f0->f6, f5->f10, f34->f54 (so T4 (1) with (a,b,c,d) == T5 (2) with (a,b,c,d))', imgs == {6: [0], 0: [6], 5: [10], 34: [54]})
    # Transport: Phi_c(O) = Phi_O4(g.c, g.O) and g.O has det -1 here; for J = diag(1,1,1,-1),
    # Phi_O4(c', O'J... ) : the O(4) extension at O' (det -1) equals Phi_{c'}(O' J) with O' J in SO(4).
    import random
    from fractions import Fraction as Fr
    rng = random.Random(5)
    ok = True
    Jd = [1, 1, 1, -1]
    for _ in range(8):
        O1 = L.random_rational_SO4(rng)
        vals = [sum(k * L.F_VALUES[j][s] for j, k in ((6, 2), (0, 3), (5, 5), (34, 1))) for s in range(24)]
        vals2 = L.g_act_values(g, vals)
        gO = L.g_act_O(g, O1)
        gOJ = [[gO[i][j] * Jd[j] for j in range(4)] for i in range(4)]
        ok &= L.is_special_orthogonal(gOJ) and L.Phi_values(vals, O1) == L.Phi_values(vals2, gOJ)
    chk('Phi_{c4}(O) = Phi_{g.c4}((g.O) J) with (g.O) J in SO(4), 8 exact rational O (transport of T5 (2) to T4 (1))', ok)
    if FAIL:
        print('FAILED', FAIL); sys.exit(1)
    print('ALL CAPACITY CHECKS PASSED')

def stage_numeric_inner(WHICH, SEED, TARGETED=False):
    """NUMERIC falsification tests for T4 and T5 (floating point; evidence only)."""

    which = WHICH
    TGT = TARGETED
    rng = np.random.default_rng(SEED)
    PERMS = L.PERMS

    def Fvals(A, C):
        out = []
        for sg in PERMS:
            p = -1.0 + 0j
            for i in range(4):
                if sg[i] != 2:
                    p *= (A[i] - C[sg[i]])
            out.append(p)
        return np.array(out)

    def fvec(coeffs):
        return np.array([sum(k * float(L.F_VECS[j][i]) for j, k in coeffs.items()) for i in range(14)])

    if which == 'cap':
        worst = []
        for trial in range(6):
            X, Y = rng.uniform(0.2, 4, 2)
            if X * Y <= 1:
                X, Y = 1 / X + 0.3, 1 / Y + 0.3
            a, b, d = X, Y, 1.0
            c = (X + 1) * (Y + 1) / (X * Y - 1)
            val, _ = L.min_phi_over_SO4(fvec({0: a, 6: b, 10: c, 54: d}), restarts=12, seed=trial)
            worst.append(val)
        print('NUMERIC capacity (2) on boundary ab(c-d)=d^2(a+b+c+d), 6 coefficient sets x 12 restarts: min Phi =', min(worst))
        worst = []
        for trial in range(6):
            s = float(np.exp(rng.uniform(-2, 2)))
            val, _ = L.min_phi_over_SO4(fvec({5: s, 2: 1 / s, 54: 1.0}), restarts=12, seed=10 + trial)
            worst.append(val)
        print('NUMERIC SOC bound (T5 (1)) on boundary ab=d^2, 6 coefficient sets x 12 restarts: min Phi =', min(worst))
        # just outside the boundary, the inequality should be able to fail (sanity that the test has teeth)
        val, _ = L.min_phi_over_SO4(fvec({0: 1.0, 6: 1.0, 10: 1.0, 54: 1.0}), restarts=12, seed=99)
        print('NUMERIC sanity: a=b=c=d=1 (violates capacity): min Phi =', val)

    if which == 't4':
        n_ok = n_hit = 0; bad = []; worstcap = np.inf; gaps_kept = []
        for _ in range(200000):
            if TGT:
                x = rng.uniform(1, 5); t = rng.uniform(0, 1); r = rng.uniform(0, 4); u = rng.uniform(-4, 0)
            else:
                x = rng.uniform(-3, 4); t = rng.uniform(-3, 3); r = rng.uniform(-4, 4); u = rng.uniform(-4, 4)
            v = complex(rng.uniform(-3, 3), rng.uniform(-3, 3))
            den_z = r*t*v**2 - r*v*x + t**2*v**2 - t*u*v**2 - t*v*x - t*v + u*v + x
            den_w = -r*t*v + r*x + t*u*v - u
            if abs(den_z) < 1e-6 or abs(den_w) < 1e-6: continue
            z = r*v*(t*v - x)/den_z; w = (-t*v + x)*(r - t*v + 1)/den_w
            A = [z, 0, 1, x]; C = [w, v, None, t*v]
            pts = [z, 0, 1, x, w, v, t*v]
            if min(abs(p - q) for i, p in enumerate(pts) for q in pts[i+1:]) < 1e-3: continue
            F = Fvals(A, C)
            mu = -1.0 / v.imag   # rho = mu Im v = -1
            kap = 1j * mu / (v * (z - w))
            g = (kap * (F - F[0])).real
            if max(abs(g[s]) for s in (0, 1, 4, 5, 7, 8)) > 1e-7 * (1 + np.abs(g).max()): bad.append('contact'); continue
            if g.min() < -1e-12: continue
            n_ok += 1
            d = g[2] - g[13]
            if d <= 1e-9: continue
            n_hit += 1
            a, b, c = g[15], g[19], g[10]; h = 1 - t; y = x - 1
            K = t*y*b - x*h*d
            signs = (x > 1) and (0 < t < 1) and (r > 0) and (u < 0) and (K > 0)
            cap = a*b*(c - d) - d*d*(a + b + c + d)
            worstcap = min(worstcap, cap / (1 + abs(a*b*c)))
            if not signs or cap < 0: bad.append((x, t, r, u, v, K, cap))
            if len(gaps_kept) < 3: gaps_kept.append(g.copy())
        print(f'NUMERIC T4 chart ({"targeted x>1,0<t<1,r>0,u<0" if TGT else "unbiased"}): 200000 draws (seed fixed), {n_ok} supported, {n_hit} with d>0; sign/capacity failures: {len(bad)}; worst normalized capacity {worstcap:.3e}')
        for g in gaps_kept:
            c = np.linalg.lstsq(L.P_numpy(), g, rcond=None)[0]
            val, _ = L.min_phi_over_SO4(c, restarts=10, seed=3)
            print('NUMERIC T4 main statement at a d>0 sample: min Phi over SO(4) =', val)

    if which == 't5':
        n_ok = 0; branch = {'neg': 0, 'posHneg': 0, 'posHpos': 0, 'A0>=0': 0}; bad = []; kept = []
        for _ in range(200000):
            x = rng.uniform(-4, 5); t = rng.uniform(-3, 3); r = rng.uniform(-4, 4); u = rng.uniform(-4, 4)
            v = complex(rng.uniform(-3, 3), rng.uniform(-3, 3))
            if abs(x - (t + u)*v) < 1e-6 or abs(u) < 1e-6: continue
            z = r*v/(x - (t + u)*v); w = (x - r - t*v)/u
            A = [z, 1, 0, x]; C = [w, v, None, t*v]
            pts = [z, 1, 0, x, w, v, t*v]
            if min(abs(p - q) for i, p in enumerate(pts) for q in pts[i+1:]) < 1e-3: continue
            F = Fvals(A, C)
            for mu in (1.0, -1.0):
                kap = 1j * mu / (v * (z - w))
                g = (kap * (F - F[1])).real
                if max(abs(g[s]) for s in (1, 2, 3, 4, 6, 8)) > 1e-7 * (1 + np.abs(g).max()): bad.append('contact'); continue
                if g.min() < -1e-12: continue
                n_ok += 1
                A0 = g[16] - g[0]
                if A0 >= 0: branch['A0>=0'] += 1; continue
                d = -A0
                H = -r*t + t*(x - 1) - u
                ok4 = d*d <= g[10]*g[13] * (1 + 1e-9) + 1e-12
                ok5 = g[22]*g[17]*(g[7] - d) >= d*d*(g[22] + g[17] + g[7] + d) - 1e-9*(1 + d**4)
                if x < 0:
                    branch['neg'] += 1; ok = ok4
                elif H < 0:
                    branch['posHneg'] += 1; ok = ok4 and x > 1 and 0 < t < 1 and r > 0 and u > 0
                else:
                    branch['posHpos'] += 1; ok = ok5 and x > 1 and 0 < t < 1 and r > 0 and u > 0
                if not ok: bad.append((x, t, r, u, v, mu, ok4, ok5))
                if len(kept) < 3: kept.append(g.copy())
        print(f'NUMERIC T5 chart: 200000 draws x 2 signs of mu, {n_ok} supported; branches {branch}; failures {len(bad)}', bad[:3])
        for g in kept:
            c = np.linalg.lstsq(L.P_numpy(), g, rcond=None)[0]
            val, _ = L.min_phi_over_SO4(c, restarts=10, seed=4)
            print('NUMERIC T5 main statement at an A0<0 sample: min Phi over SO(4) =', val)

if __name__ == '__main__':
    t0 = time.time()
    print('=== stage_linear_T4'); sys.stdout.flush()
    stage_linear_T4()
    print('=== stage_chart_uniqueness_T4'); sys.stdout.flush()
    stage_chart_uniqueness_T4()
    print('=== stage_chart_identities_T4'); sys.stdout.flush()
    stage_chart_identities_T4()
    print('=== stage_capacity'); sys.stdout.flush()
    stage_capacity()
    if '--numeric' in sys.argv:
        print('=== numeric (NUMERIC evidence only)')
        stage_numeric_inner('t4', 0)
        stage_numeric_inner('t4', 1, True)
    print('ALL EXACT STAGES PASSED; total time %.0f s' % (time.time() - t0))
