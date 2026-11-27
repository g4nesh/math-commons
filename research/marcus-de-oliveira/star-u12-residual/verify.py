"""Exact replay for the U12=0 residual inequality: cone reduction and selector obstruction.

Rebuilds F0, Fc, F1 exactly as in star-u12-boundary/verify.py (PR 132), then checks with SymPy
(exact rationals only):
  1. F0 = t^2 (l^T X l) (X : Nn), l = (uz, -w), Nn = u^2 q^2 I + t m m^T, m = (z, -uw);
  2. the generic cone-certificate identity in frame coordinates Y (X = E Y E^T);
  3. at two exact rational points with det M < 0: three explicit choices of N fail, and an
     explicit rational N passes (so the certificate family is not empty there).
Usage: python -I -B research/marcus-de-oliveira/star-u12-residual/verify.py
"""
import sys
import sympy as sp
from sympy import Rational as Q

u, w, z, t, s = sp.symbols('u w z t s', positive=True)
X11, X22, X12 = sp.symbols('X11 X22 X12', real=True)
ok_all = True


def check(name, ok):
    global ok_all
    ok_all &= bool(ok)
    print(('ok   ' if ok else 'FAIL ') + name)


def red(e):
    return sp.expand(sp.rem(sp.expand(e), u**2 - (1 - w**2 - z**2), u))


# Construction copied verbatim from star-u12-boundary/verify.py section 2.
q2 = w**2 + z**2
one = X11 + X22
R0 = (u**2*w**2*X11 + z**2*X22 + 2*u*w*z*X12)/q2
H1 = s + (1 - s)*one
A = u**2*s
k = A + H1 - q2*(1 - s)*X11 - (1 - s)*R0
h = (1 - 2*w**2 - z**2)*H1 + q2*(1 - s)*X11 - (u**2*w**2 + z**2)*(1 - s)*one/q2 + (1 - s)*R0
d = (u**2 + z**2)*H1 - (1 - s)*R0
e = s*((u**2 + z**2)**2 - z**2)*H1
F = A**2 + A*(2*k + h)*t + (k**2 + k*h + A*(d - w**2*H1))*t**2 + (k*d + e)*t**3
Fs = sp.Poly(sp.expand(F*q2**2), s)
c0, c1, c2 = (Fs.coeff_monomial(s**i) for i in range(3))
F0, Fc, F1 = red(c0), red(c1 + 2*c0), red(c2 + c1 + c0)
check('Bernstein identity (as in PR 132)', red(F*q2**2 - ((1 - s)**2*F0 + s*(1 - s)*Fc + s**2*F1)) == 0)

Xm = sp.Matrix([[X11, X12], [X12, X22]])
l = sp.Matrix([u*z, -w])
mv = sp.Matrix([z, -u*w])
Nn = u**2*q2*sp.eye(2) + t*mv*mv.T
check('F0 = t^2 (l^T X l)(X : Nn) with trX = X11 + X22',
      red(F0 - t**2*(l.T*Xm*l)[0]*sum(Nn.multiply_elementwise(Xm))) == 0)
M = sp.Matrix([[Fc.coeff(X11), Fc.coeff(X12)/2], [Fc.coeff(X12)/2, Fc.coeff(X22)]])
check('Fc = X : M', sp.expand(Fc - (M[0, 0]*X11 + M[1, 1]*X22 + 2*M[0, 1]*X12)) == 0)

# 2. Generic identity. Frame E = [l, l_perp], l_perp = (w, uz), E^T E = l2 I, X = E Y E^T.
# With tau = sigma >= 0 and trX = 1:  sigma^2 F1 + sigma Fc + F0 = Phi(tau, Y),
#   Phi = F1 tau^2 + tau (Y : Mp) + t^2 l2^2 y1 (Y : Np),  Mp = E^T M E,  Np = E^T Nn E.
f1, tt, L2, n11, mm, tau, y1, y2, y3 = sp.symbols('f1 t l2 n11 m tau y1 y2 y3', real=True)
a11, a12, a22, b11, b12, b22 = sp.symbols('a11 a12 a22 b11 b12 b22', real=True)
Yi = lambda P: P[0]*y1 + P[2]*y2 + 2*P[1]*y3                 # P = (P11, P12, P22)
Phi = f1*tau**2 + tau*Yi((a11, a12, a22)) + tt**2*L2**2*y1*Yi((b11, b12, b22))
cert = (f1*tau**2 + 2*tau*((a11 - n11)/2*y1 + (a12 - mm)*y3)
        + tt**2*L2**2*(b11*y1**2 + 2*b12*y1*y3 + b22*y3**2)    # quadratic form of P3 below
        + tt**2*L2**2*b22*(y1*y2 - y3**2)                      # mu det Y, mu = t^2 l2^2 Np22 >= 0
        + tau*Yi((n11, mm, a22)))                               # tau (Y : N), N = [[n11, m], [m, Mp22]]
check('cone certificate identity Phi = P3 form + mu det Y + tau (Y:N)', sp.expand(Phi - cert) == 0)
Esym = sp.Matrix([[u*z, w], [-w, u*z]])
l2 = sp.expand(u**2*z**2 + w**2)
check('frame: E^T E = l2 I and E^T l = (l2, 0)',
      sp.expand(Esym.T*Esym - l2*sp.eye(2)) == sp.zeros(2) and sp.expand(Esym.T*l - sp.Matrix([l2, 0])) == sp.zeros(2, 1))
Y = sp.Matrix([[y1, y3], [y3, y2]])
XY = Esym*Y*Esym.T
check('frame: l^T X l = l2^2 y1, X : Nn = Y : Np, X : M = Y : Mp',
      all(sp.expand(a - b) == 0 for a, b in [
          ((l.T*XY*l)[0], l2**2*y1),
          (sum(Nn.multiply_elementwise(XY)), sum((Esym.T*Nn*Esym).multiply_elementwise(Y))),
          (sum(M.multiply_elementwise(XY)), sum((Esym.T*M*Esym).multiply_elementwise(Y)))]))


# 3. Exact points. P3 = [[F1, h1, h3], [h1, k11, k12], [h3, k12, k22]], k = t^2 l2^2 Np.
def p3(pt, n11v, mv_):
    sub = dict(zip((u, w, z, t), pt))
    Mp = (Esym.T*M*Esym).subs(sub)
    Np = (Esym.T*Nn*Esym).subs(sub)
    L = l2.subs(sub)
    K = pt[3]**2*L**2*Np
    P = sp.Matrix([[F1.subs(sub), (Mp[0, 0] - n11v)/2, Mp[0, 1] - mv_],
                   [(Mp[0, 0] - n11v)/2, K[0, 0], K[0, 1]],
                   [Mp[0, 1] - mv_, K[0, 1], K[1, 1]]])
    return P, Mp, Np


def pd(P):
    return all(P[:i, :i].det() > 0 for i in (1, 2, 3))


def psd(P):
    idx = [(0,), (1,), (2,), (0, 1), (0, 2), (1, 2), (0, 1, 2)]
    return all(P.extract(list(i), list(i)).det() >= 0 for i in idx)


def sphere(a, b):
    D = 1 + a*a + b*b
    return (2*a/D, (1 - a*a - b*b)/D, 2*b/D)


points = {
    'W1': (*sphere(Q(1, 491), Q(1, 887)), Q(813)),
    'W2': (*sphere(Q(1, 200000), Q(1, 344)), Q(37951, 10)),
}
PASS_M = {'W1': Q(5877, 1000), 'W2': Q(17447, 100)}                   # rational m* in normalised frame units
for name, pt in points.items():
    sub = dict(zip((u, w, z, t), pt))
    check(f'{name}: on the sphere, t > 0', pt[0]**2 + pt[1]**2 + pt[2]**2 == 1 and pt[3] > 0)
    check(f'{name}: det M < 0 (so Fc < 0 somewhere and M is not PSD)', M.det().subs(sub) < 0)
    _, Mp, Np = p3(pt, 0, 0)
    L = l2.subs(sub)
    a22 = Mp[1, 1]
    check(f'{name}: Mp22 = Fc(l_perp) > 0', a22 > 0)
    sel = {'N = Fc(e2) e2e2^T (m = 0, n11 = 0)': (0, 0),
           'N rank one with h3 = 0 (m = Mp12)': (Mp[0, 1]**2/a22, Mp[0, 1])}
    a, b, c = Np[1, 1], -Np[0, 1], Np[0, 0]
    mN = -a22*(b*Mp[0, 0]/2 + c*Mp[0, 1])/(a*Mp[0, 0]/2 + b*Mp[0, 1])
    sel['N rank one, m = one Newton step from 0'] = (mN**2/a22, mN)
    for lab, (nv, mv_) in sel.items():
        check(f'{name}: selector fails: {lab}', not psd(p3(pt, nv, mv_)[0]))
    m_u = PASS_M[name]*L                                          # frame Y uses unnormalised E
    P, _, _ = p3(pt, m_u**2/a22, m_u)
    check(f'{name}: explicit N = [[m^2/Mp22, m], [m, Mp22]] passes (P3 positive definite)', pd(P))

print('ALL CHECKS PASSED' if ok_all else 'SOME CHECKS FAILED')
sys.exit(0 if ok_all else 1)
