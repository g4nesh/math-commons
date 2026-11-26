"""Own replay for blind_moc4_inventory52_native_harmonic_conditional_20261023.
Chart A = (y, z, w, x), C-points (t, 0, 1, inf); F_sigma = prod_i (A_i - Cpt_sigma(i)) with the
factor at the infinite point dropped (a positive normalization).  Permutations in the library's
lexicographic order.  Everything is rebuilt from definitions with sympy; nothing from the
supplement is read or executed."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import moc4lib as L
import sympy as sp

fails = 0
def check(name, ok, info=""):
    global fails
    print(("PASS " if ok else "FAIL ") + name + (" :: " + str(info) if info != "" else ""))
    if not ok:
        fails += 1

y, zz, w, x, t, r, s, u, h, J, W, eps, T, K = sp.symbols('y z w x t r s u h J W epsilon T K')
Apts = [y, zz, w, x]
Cpts = [t, 0, 1, None]
def F(sig, Ap=Apts, Cp=Cpts):
    p = sp.Integer(1)
    for i in range(4):
        if Cp[sig[i]] is not None:
            p *= (Ap[i] - Cp[sig[i]])
    return sp.expand(p)
Fs = [F(sg) for sg in L.PERMS]
E = Fs[0] - Fs[6]
check("E = F0 - F6 = t(y - z)(w - 1)", sp.expand(E - t * (y - zz) * (w - 1)) == 0)
# solve the tree equations with (3)
V = r * (zz - w) + s * zz
Dt = -V - u * zz * (w - 1) + zz
Dy = V + (w - zz) * (1 - u * (w - 1))
sol = {x: (V + s * zz * (w - 1)) / V, t: w * zz / Dt, y: w * (V + w - zz) / Dy}
sub = lambda e: sp.together(e.subs(sol))
Es = sub(E)
eqs = {"r": Fs[0] - Fs[1] - r * E, "s": Fs[1] - Fs[4] - s * E, "u": Fs[6] - Fs[12] - u * E}
for k, e in eqs.items():
    check("(3) solves the tree equation for %s" % k, sp.simplify(sub(e)) == 0)
Hs = [sp.together((Fs[k] - Fs[0]) / E) for k in range(24)]
# J chart
Jexpr = V / (zz * (w - 1))
Tq = J * (W - 1)
ZD = r + s - Tq
# z = r W / ZD with J = V/(z(w-1)), W = w
zsol = sp.solve(sp.Eq(J, Jexpr.subs(w, W)), zz)
check("J = V/(z(w-1)) inverts to z = r W/(r + s - T)", len(zsol) == 1 and sp.simplify(zsol[0] - r * W / ZD) == 0, zsol)
chart = {zz: r * W / ZD, w: W}
def inJW(k):
    return sp.factor(sp.simplify(sub(Hs[k] * 1).subs(chart)))
amp = {}
H2 = inJW(2)
check("H_2 - J is a constant in r, s, u (so d = -eta Im J)", sp.simplify(H2 - J).free_symbols <= {r, s, u}, sp.simplify(H2 - J))
H7 = inJW(7)
check("f9 amplitude: H_7 - (-r - s/T) is constant in r,s,u", sp.simplify(H7 - (-r - s / Tq)).free_symbols <= {r, s, u}, sp.simplify(H7 - (-r - s / Tq)))
H18 = inJW(18)
check("f0 amplitude: H_18 - (-r - s/T + s u/J) constant", sp.simplify(H18 - (-r - s / Tq + s * u / J)).free_symbols <= {r, s, u}, sp.simplify(H18 - (-r - s / Tq + s * u / J)))
H5 = inJW(5)
check("f21 amplitude: H_5 - (-J - r + (J+r+s)/W) constant", sp.simplify(H5 - (-J - r + (J + r + s) / W)).free_symbols <= {r, s, u}, sp.simplify(H5 - (-J - r + (J + r + s) / W)))
L0 = (r - 1) * Tq + s
Lq = L0 + r * u * (W - 1)
H14 = inJW(14)
check("f7 amplitude: H_14 - (-u((J+r)(W-1) - s)/L) constant", sp.simplify(H14 - (-u * ((J + r) * (W - 1) - s) / Lq)).free_symbols <= {r, s, u}, sp.simplify(H14 - (-u * ((J + r) * (W - 1) - s) / Lq)))
H15 = inJW(15)
# H9 = N/D, curve C = N - h D
H9 = inJW(9)
Nn, Dd = sp.fraction(sp.together(H9))
Cc = sp.expand(Nn - h * Dd)
PJ = sp.Poly(Cc, J, W)
print("INFO H9 =", H9)
print("INFO bidegree of C in (J, W):", sp.degree(Cc, J), sp.degree(Cc, W))
# normalize so that C(J,1) = -s^2
c1 = sp.factor(Cc.subs(W, 1))
lam = sp.simplify(c1 / (-s ** 2))
check("C(J,1) is a nonzero constant multiple of -s^2 (normalization factor found)", lam.free_symbols <= {r, s, u, h} and lam != 0, lam)
Cn = sp.expand(sp.cancel(Cc / lam))
check("bidegree (3,3)", sp.degree(Cn, J) == 3 and sp.degree(Cn, W) == 3)
check("C(0,W) = -s(rW - r - s)(uW - u - 1)", sp.expand(Cn.subs(J, 0) + s * (r * W - r - s) * (u * W - u - 1)) == 0)
check("C(J,0) = (J+r+s)(J+u+1)((r-1)J - s)", sp.expand(Cn.subs(W, 0) - (J + r + s) * (J + u + 1) * ((r - 1) * J - s)) == 0)
check("C(J,1) = -s^2", sp.expand(Cn.subs(W, 1) + s ** 2) == 0)
check("[J^3]C = (r-1)(W-1)^2", sp.expand(sp.Poly(Cn, J).coeff_monomial(J ** 3) - (r - 1) * (W - 1) ** 2) == 0)
check("[W^3]C = -J(h+r)((r-1)J + r u)", sp.expand(sp.Poly(Cn, W).coeff_monomial(W ** 3) + J * (h + r) * ((r - 1) * J + r * u)) == 0)
resL = sp.factor(sp.resultant(Cn, sp.expand(Lq), W))
check("Res_W(C, L) = r s^2 u (J+r+s)(J+u+1)((r-1)J + r u) up to sign", sp.simplify(resL / (r * s ** 2 * u * (J + r + s) * (J + u + 1) * ((r - 1) * J + r * u))) in (1, -1), resL)
Hq = s - Tq
resH = sp.factor(sp.resultant(Cn, sp.expand(Hq), W))
QH = -J ** 3 + (h - s - u) * J ** 2 + (h * (s + u) + r * s - r - s * u + u) * J + s * u * (h + r + 1)
check("Res_W(C, H) = -J r s^2 Q_H(J) up to sign", sp.simplify(resH / (J * r * s ** 2 * QH)) in (1, -1), resH)
# J = infinity branches
Eps = sp.expand(eps ** 3 * J ** 2 * Cn.subs(W, 1 + T / J)).subs(J, 1 / eps)
Eps = sp.expand(sp.simplify(Eps))
check("E(eps,T) polynomial and E(0,T) = T((r-1)T + s)", sp.expand(Eps.subs(eps, 0) - T * ((r - 1) * T + s)) == 0 and sp.Poly(Eps, eps, T) is not None)
cor = sp.expand(sp.simplify(eps ** 5 * Cn.subs({J: 1 / eps, W: K / eps})))
check("[eps^5 C(1/eps, K/eps)]_{eps=0} = -K^2 (r-1)(K(h+r) - 1)", sp.expand(cor.subs(eps, 0) + K ** 2 * (r - 1) * (K * (h + r) - 1)) == 0)
# discriminant example
QH0 = QH.subs({r: 2, s: 3, u: 4, h: 5})
check("Q_H at (2,3,4,5) = -J^3 - 2J^2 + 31J + 96, disc -15616", sp.expand(QH0 - (-J ** 3 - 2 * J ** 2 + 31 * J + 96)) == 0 and sp.discriminant(QH0, J) == -15616)
# inverse divisor E = -J W^2 (W-1) H L in the chart (with the doc's homogeneous inverse we only check the
# ratio F0/E = -rT/H in our normalization)
F0E = sp.factor(sp.simplify(sub(Fs[0] / E).subs(chart)))
check("F0/E = -r T / H", sp.simplify(F0E - (-r * Tq / Hq)) == 0, F0E)
print("RESULT", "all EXACT checks passed" if fails == 0 else "%d failures" % fails)
sys.exit(1 if fails else 0)
