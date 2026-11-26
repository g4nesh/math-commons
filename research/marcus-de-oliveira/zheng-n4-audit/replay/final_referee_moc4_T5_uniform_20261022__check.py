"""Exact replay (and brief numeric tests) for the audit of final_referee_moc4_T5_uniform_20261022.
Own code only; nothing from the supplement is executed or read.  Run from this directory:
    python -I -B final_referee_moc4_T5_uniform_20261022__check.py [--numeric]
Prints one line per check; exits nonzero if an EXACT check fails."""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lib'))
import moc4lib as L
import sympy as sp
import numpy as np

def stage_linear_T5():
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

    # ---------------- T5 ----------------
    T5 = [1, 2, 3, 4, 6, 8]
    chk('T5 equals paper CONTACT_SETS[5]', set(T5) == set(L.CONTACT_SETS[5]))
    gen5 = L.generators(T5)
    chk('T5 generators = {0,1,2,3,4,5,6,9,10,53,54}: got %s' % gen5, gen5 == [0, 1, 2, 3, 4, 5, 6, 9, 10, 53, 54])
    chk('T5 only f54 non-minor', [j for j in gen5 if not L.F_IS_MINOR[j]] == [54])
    print('INFO f54 =', L.f_name(54), ' f53 =', L.f_name(53), ' f9 =', L.f_name(9), ' f10 =', L.f_name(10))
    chk('f54 + f1 = f0 + f6 + f53', comb({54: 1, 1: 1}) == comb({0: 1, 6: 1, 53: 1}))
    chk('f54 + f4 = f2 + f3 + f6 + f53', comb({54: 1, 4: 1}) == comb({2: 1, 3: 1, 6: 1, 53: 1}))
    chk('f54 + f9 = f2 + f5 + f10 + f53', comb({54: 1, 9: 1}) == comb({2: 1, 5: 1, 10: 1, 53: 1}))
    A0 = lambda vals: vals[16] - vals[0]
    chk('A0 = g16 - g0 on generators: -1 on f54, +1 on f1,f4,f9, 0 else',
        all(A0(V[j]) == {54: -1, 1: 1, 4: 1, 9: 1}.get(j, 0) for j in gen5))
    chk('f54 positive at 0, 9, 11', all(V[54][s] > 0 for s in (0, 9, 11)))
    print('INFO f54 values', V[54])

    def dec5(g):
        d = -(g[16] - g[0])
        return {0: g[22], 2: g[13], 3: g[12], 5: g[10], 6: g[17], 10: g[7], 53: g[16], 54: d}
    B5 = check_decomp(T5, 'T5 eq.(3)', dec5)
    chk('T5 face space dim 8', len(B5) == 8)
    chk('T5 generator span rank 8', L.rank([list(V[j]) for j in gen5]) == 8)
    # T5 group check: T4 / T5 are distinct orbits
    chk('T4, T5 in different G-orbits', L.find_g([0, 1, 4, 5, 7, 8], T5) is None)
    # (1) of T5 is "a relabeling of the old SOC master"; (1) of T4 vs (2) of T5: same shape under G?
    found = None
    for g in L.G_ELEMENTS:
        t = L.g_index_map(g)
        imgs = {}
        for j in (6, 0, 5, 34):
            vals = L.g_act_values(g, V[j])
            jj = [k for k in range(82) if list(V[k]) == list(vals)]
            imgs[j] = jj[0] if jj else None
        if imgs == {6: 0, 0: 6, 5: 10, 34: 54}:
            found = g; break
    chk('some g in G maps (f6,f0,f5,f34) of T4-(1) to (f0,f6,f10,f54) of T5-(2): %s' % (found,), found is not None)
    if FAIL:
        print('FAILED:', FAIL); sys.exit(1)
    print('ALL LINEAR CHECKS PASSED')

def stage_chart_T5():

    T0 = time.time()
    FAIL = []
    def chk(name, expr):
        ok = sp.cancel(sp.together(expr)) == 0
        print(('PASS ' if ok else 'FAIL ') + name); sys.stdout.flush()
        if not ok: FAIL.append(name)

    z, w, v, vb, x, t, r, u = sp.symbols('z w v vb x t r u')
    A = [z, 1, 0, x]
    C = [w, v, None, t * v]
    SIGN = -1   # sign of the dropped factor at infinity
    def F(sig):
        p = sp.Integer(SIGN)
        for i in range(4):
            if sig[i] != 2:
                p *= (A[i] - C[sig[i]])
        return sp.expand(p)
    Fs = [F(sg) for sg in L.PERMS]
    T5 = [1, 2, 3, 4, 6, 8]
    Rraw = {s: sp.factor(sp.cancel((Fs[s] - Fs[1]) / (v * (z - w)))) for s in T5}
    print('raw R at T5:', Rraw)
    rdef = z * (v - w) * (x - t * v) / (v * (z - w))
    udef = (z - v) * (x - t * v) / (v * (z - w))
    # transposition differences named in the text (proportional to): print F_s - F_1 factored
    for s in T5:
        print('F_%d - F_1 =' % s, sp.factor(Fs[s] - Fs[1]))
    e1 = sp.numer(sp.together(rdef - r)); e2 = sp.numer(sp.together(udef - u))
    sol = sp.solve([e1, e2], [z, w], dict=True)
    print('solutions:', [{k: sp.factor(q) for k, q in so.items()} for so in sol])
    zc = r * v / (x - (t + u) * v); wc = (x - r - t * v) / u
    chk('(6) satisfies r-equation', e1.subs({z: zc, w: wc}))
    chk('(6) satisfies u-equation', e2.subs({z: zc, w: wc}))
    need = [0, 7, 9, 10, 11, 12, 13, 16, 17, 22] + T5
    Rv = {s: sp.cancel(((Fs[s] - Fs[1]) / (v * (z - w))).subs({z: zc, w: wc})) for s in need}
    print('R at T5 in chart:', {s: sp.factor(Rv[s]) for s in T5})
    g = {}
    for s in need:
        g[s] = sp.cancel((Rv[s] - Rv[s].subs(v, vb)) / (v - vb))   # g_s = -rho * g[s]
    rho, sig = sp.symbols('rho sig')
    N = v * vb; E = (x - t*v) * (x - t*vb); Fd = (x - r - (t + u)*v) * (x - r - (t + u)*vb)
    G = {s: -rho * g[s] for s in need}
    Ls = x - r - t - u; D0 = r*t + u*x; J = r + (t + u - 1)*x; H = -r*t + t*(x - 1) - u; W = t*(x - 1) - r - u
    sg = rho / E
    A0 = G[16] - G[0]
    chk('g0 = -rho x/N', G[0] + rho*x/N)
    chk('g7 = -sigma t x (r+u)', G[7] + sg*t*x*(r + u))
    chk('g9 = -sigma t x (r+ux)', G[9] + sg*t*x*(r + u*x))
    chk('g10 = -sigma r t (x-1)', G[10] + sg*r*t*(x - 1))
    chk('x g11 = (x-1) g9', x*G[11] - (x - 1)*G[9])
    chk('g12 = -rho r u L/F', G[12] + rho*r*u*Ls/Fd)
    chk('g17 = A0 - sigma t u x (x-1)', G[17] - (A0 - sg*t*u*x*(x - 1)))
    chk('H g13 = t L A0 + sigma t u x (t-1) W', H*G[13] - (t*Ls*A0 + sg*t*u*x*(t - 1)*W))
    chk('H g22 = (t-1)(x-r)A0 + sigma r t x (t-1)(x-r-u-1)', H*G[22] - ((t - 1)*(x - r)*A0 + sg*r*t*x*(t - 1)*(x - r - u - 1)))
    chk('D0 A0 = sigma r u [x t (t-1)(x-1) - J H E/F]', D0*A0 - sg*r*u*(x*t*(t - 1)*(x - 1) - J*H*E/Fd))
    for s in T5:
        chk(f'contact {s}: gap vanishes identically in chart', G[s])
    # auxiliary algebra used in the sign arguments
    chk('H = tL + (t-1)(t+u)', H - (t*Ls + (t - 1)*(t + u)))
    chk('W = H + r(t-1)', W - (H + r*(t - 1)))
    # negative branch, (8)
    X, T, R, U, d = sp.symbols('X T R U d')
    negsub = {x: -X, t: -T, r: -R, u: U}
    Hn = T*(X + 1 - R) - U; Wn = R - U + T*(X + 1); Bn = T*U*X*(X + 1); Ln = R + T - X - U
    chk('neg: H', H.subs(negsub) - Hn); chk('neg: W', W.subs(negsub) - Wn); chk('neg: L', Ls.subs(negsub) - Ln)
    # with sigma = 1, A0 = -d, the linear-in-(A0, sigma) identities
    g10n = (-1*r*t*(x - 1)).subs(negsub); chk('neg: g10 = RT(X+1)', g10n - R*T*(X + 1))
    g17n = (-d - t*u*x*(x - 1)).subs(negsub); chk('neg: g17 = B - d', g17n - (Bn - d))
    Hg13n = (t*Ls*(-d) + t*u*x*(t - 1)*W).subs(negsub); chk('neg: Hg13 = T[Ld - UX(T+1)W]', Hg13n - T*(Ln*d - U*X*(T + 1)*Wn))
    Hg22n = ((t - 1)*(x - r)*(-d) + r*t*x*(t - 1)*(x - r - u - 1)).subs(negsub)
    chk('neg: Hg22 = (T+1)[(R-X)d + RTX(R-X-U-1)]', Hg22n - (T + 1)*((R - X)*d + R*T*X*(R - X - U - 1)))
    g7n = (-t*x*(r + u)).subs(negsub); chk('neg: g7 = TX(R-U)', g7n - T*X*(R - U))
    da = R*T*X*(X + U + 1 - R)/(R - X)
    chk('neg: da - B = TX(R+UX)(X+1-R)/(R-X)', da - Bn - T*X*(R + U*X)*(X + 1 - R)/(R - X))
    Psi = lambda dd: R*T*(X + 1)*(T/Hn)*(Ln*dd - U*X*(T + 1)*Wn) - dd**2
    Qn = -(X + 1)*Hn - R*(T + 1)
    chk('neg: Psi(B) = T^2 U X (R+UX)(X+1) Q/H', Psi(Bn) - T**2*U*X*(R + U*X)*(X + 1)*Qn/Hn)
    Pn = R*(R - U) + (X + 1)*(R*T - (T + 1)*X)
    chk('neg: Psi(da) = R T^2 X (R+UX) P/(R-X)^2', Psi(da) - R*T**2*X*(R + U*X)*Pn/(R - X)**2)
    Umax = R + T - X
    chk('neg: Q at U=R+T-X equals -X(T+1)(X+1-R)', Qn.subs(U, Umax) + X*(T + 1)*(X + 1 - R))
    chk('neg: dQ/dU = X+1 (>0)', sp.diff(Qn, U) - (X + 1))
    chk('neg: P at U=R+T-X equals X(T+1)(R-X-1)', Pn.subs(U, Umax) - X*(T + 1)*(R - X - 1))
    chk('neg: dP/dU = -R (<0)', sp.diff(Pn, U) + R)
    # positive branch, all-negative exclusion
    Tq, Rq, Uq = sp.symbols('Tq Rq Uq')
    allneg = {t: -Tq, r: -Rq, u: -Uq}
    Hg22a = ((t - 1)*(x - r)*(-d) - r*t*x*(t - 1)*(x - r - u - 1)).subs(allneg)   # sigma = -1
    chk('pos/all-neg: Hg22 = (T+1)[(x+R)d + RTx(x+R+U-1)]', Hg22a - (Tq + 1)*((x + Rq)*d + Rq*Tq*x*(x + Rq + Uq - 1)))
    chk('pos/all-neg: H = U - T(R+x-1)', H.subs(allneg) - (Uq - Tq*(Rq + x - 1)))
    chk('pos/all-neg: D0 = RT - Ux', D0.subs(allneg) - (Rq*Tq - Uq*x))
    chk('pos/all-neg: D0 at U=T(R+x-1) is -T(x-1)(R+x)', (Rq*Tq - Uq*x).subs(Uq, Tq*(Rq + x - 1)) + Tq*(x - 1)*(Rq + x))
    chk('pos/all-neg: d = (RU/D0)[xT(T+1)(x-1) - JHE/F]',
        (-(-1)*r*u*(x*t*(t - 1)*(x - 1) - J*H*E/Fd)/D0).subs(allneg) - (Rq*Uq/D0.subs(allneg))*(x*Tq*(Tq + 1)*(x - 1) - (J*H).subs(allneg)*E.subs(allneg)/Fd.subs(allneg)))
    # positive branch (9) with sigma = -1, h = 1-t, y = x-1
    h, y = sp.symbols('h y')
    pos = {t: 1 - h, x: y + 1}
    Hp = H.subs(pos); Wp = W.subs(pos); Lp = Ls.subs(pos)
    chk('pos: Hg13 = t[uxhW - Ld]', ((t*Ls*(-d) - t*u*x*(t - 1)*W)).subs(pos) - ((1 - h)*(u*(y + 1)*h*Wp - Lp*d)))
    chk('pos: Ha = h[(x-r)d + rtx(y-r-u)]', (((t - 1)*(x - r)*(-d) - r*t*x*(t - 1)*(x - r - u - 1))).subs(pos) - h*((y + 1 - r)*d + r*(1 - h)*(y + 1)*(y - r - u)))
    chk('pos: g17 = tuxy - d', (-d + t*u*x*(x - 1)).subs(pos) - ((1 - h)*u*(y + 1)*y - d))
    chk('pos: g7 = tx(r+u)', (t*x*(r + u)).subs(pos) - (1 - h)*(y + 1)*(r + u))
    chk('pos: H = tL - h(t+u)', Hp - ((1 - h)*Lp - h*(1 - h + u)))
    chk('pos: L=0 => W = (t-1)(r+t+u)', Wp.subs(u, sp.solve(Lp, u)[0]) - (-h)*(r + 1 - h + sp.solve(Lp, u)[0]))
    tt = 1 - h; xx = y + 1
    Bp = tt*u*xx*y; Cp = tt*xx*(r + u)
    dap = r*tt*xx*(r + u - y)/(xx - r)
    chk('pos: sign(da - B) ~ (r - y): da - B = t x (r+ux)(r-y)/(x-r)', dap - Bp - tt*xx*(r + u*xx)*(r - y)/(xx - r))
    Psip = lambda dd: r*tt*y*(tt/Hp)*(u*xx*h*Wp - Lp*dd) - dd**2
    Pp = r*(r + u) - y*(r*tt + h*xx); Qp = y*Hp + r*h
    chk('pos H<0: Psi(da) = -r t^2 x (r+ux) P/(x-r)^2', Psip(dap) + r*tt**2*xx*(r + u*xx)*Pp/(xx - r)**2)
    chk('pos H<0: Psi(B) = -t^2 u x (r+ux) y Q/H', Psip(Bp) + tt**2*u*xx*(r + u*xx)*y*Qp/Hp)
    chk('pos H<0: P at u = y+h-r equals -h x (y-r)', Pp.subs(u, y + h - r) + h*xx*(y - r))
    chk('pos H<0: dP/du = r', sp.diff(Pp, u) - r)
    k = sp.symbols('k')
    chk('pos H<0: Q = y(h-u-tk) + kh with r = y+k', Qp.subs(r, y + k) - (y*(h - u - tt*k) + k*h))
    Kp = Bp*Cp/(Bp + Cp)
    chk('pos H>0: K = tuxy(r+u)/(r+ux)', Kp - tt*u*xx*y*(r + u)/(r + u*xx))
    aa = lambda dd: (h/Hp)*((xx - r)*dd + r*tt*xx*(y - r - u))
    chi = lambda dd: aa(dd)*(Kp - dd) - dd**2
    chk('pos H>0: a(B-d)(C-d) - d^2(a + B + C) = (B+C) chi(d)', aa(d)*(Bp - d)*(Cp - d) - d**2*(aa(d) + (Bp - d) + Cp + d) - (Bp + Cp)*chi(d))
    chk('pos H>0: chi quadratic coeff = -L/H', sp.Poly(sp.numer(sp.together(chi(d))), d).coeff_monomial(d**2) / sp.denom(sp.together(chi(d))) + Lp/Hp)
    dstar = u*xx*h*Wp/Lp
    Wsym = sp.symbols('Wsym')
    chival = chi(dstar)
    Pstar = (r*Wsym**2 + (r + u)*(r*(2 - tt) + u*(1 - tt))*Wsym + (r + u)**2*(1 - tt)*(r + tt + u)) / tt
    # substitute x = 1 + (r+u+W)/t, i.e. y = (r+u+W)/t
    ysub = (r + u + Wsym) / tt
    chk('pos H>0: W(y=(r+u+W)/t) = W', Wp.subs(y, ysub) - Wsym)
    chk('pos H>0: chi(d*) = u x^2 (r+u) h P*/((r+ux) L) with t P* as displayed',
        (chival - u*xx**2*(r + u)*h*Pstar/((r + u*xx)*Lp)).subs(y, ysub))
    chk('pos H>0: y - r - u = ((r+u)(1-t) + W)/t', (y - r - u).subs(y, ysub) - ((r + u)*(1 - tt) + Wsym)/tt)
    print('time', round(time.time() - T0, 1))
    if FAIL:
        print('FAILED', FAIL); sys.exit(1)
    print('ALL T5 CHECKS PASSED')

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
        for trial in range(3):
            X, Y = rng.uniform(0.2, 4, 2)
            if X * Y <= 1:
                X, Y = 1 / X + 0.3, 1 / Y + 0.3
            a, b, d = X, Y, 1.0
            c = (X + 1) * (Y + 1) / (X * Y - 1)
            val, _ = L.min_phi_over_SO4(fvec({0: a, 6: b, 10: c, 54: d}), restarts=4, seed=trial)
            worst.append(val)
        print('NUMERIC capacity (2) on boundary ab(c-d)=d^2(a+b+c+d), 3 coefficient sets x 4 restarts: min Phi =', min(worst))
        worst = []
        for trial in range(3):
            s = float(np.exp(rng.uniform(-2, 2)))
            val, _ = L.min_phi_over_SO4(fvec({5: s, 2: 1 / s, 54: 1.0}), restarts=4, seed=10 + trial)
            worst.append(val)
        print('NUMERIC SOC bound (T5 (1)) on boundary ab=d^2, 3 coefficient sets x 4 restarts: min Phi =', min(worst))
        # just outside the boundary, the inequality should be able to fail (sanity that the test has teeth)
        val, _ = L.min_phi_over_SO4(fvec({0: 1.0, 6: 1.0, 10: 1.0, 54: 1.0}), restarts=4, seed=99)
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
    print('=== stage_linear_T5'); sys.stdout.flush()
    stage_linear_T5()
    print('=== stage_chart_T5'); sys.stdout.flush()
    stage_chart_T5()
    print('=== stage_capacity'); sys.stdout.flush()
    stage_capacity()
    if '--numeric' in sys.argv:
        print('=== numeric (NUMERIC evidence only)')
        stage_numeric_inner('t5', 2)
    if '--numeric-cap' in sys.argv:
        stage_numeric_inner('cap', 3)
    print('ALL EXACT STAGES PASSED; total time %.0f s' % (time.time() - t0))
