"""Exact replay for the U12=0 reduction of the unitary-minor inequality.

Checks symbolic identities with SymPy (exact rationals, no floats decide anything):
  1. the U12=0 frame turns the two-row k, h, d, e into the closed forms of README section 2;
  2. the Bernstein identity q^4 F = (1-s)^2 F0 + s(1-s) Fc + s^2 F1 and the factorisation of F0;
  3. nonnegative homogenised coefficients of F1, T(s=1), M11, M22;
  4. an exact rational point where Fc < 0 and det M < 0 although F > 0 for every s.
Usage: python -B research/marcus-de-oliveira/star-u12-boundary/verify.py [--write-report]
"""
import json, pathlib, sys
import sympy as sp

u, w, z, q, t, s = sp.symbols('u w z q t s', positive=True)
X11, X22, X12 = sp.symbols('X11 X22 X12', real=True)
U, W, Z = sp.symbols('U W Z', nonnegative=True)
results = {}


def check(name, ok):
    results[name] = bool(ok)
    print(('ok   ' if ok else 'FAIL ') + name)


def red(expr):
    """Reduce modulo u^2 + w^2 + z^2 = 1."""
    return sp.expand(sp.rem(sp.expand(expr), u**2 - (1 - w**2 - z**2), u))


# 1. Reduction. Two-row section 2 formulas for real x and complex a with x.a = 0, |a| = 1.
ar, ai, br, bi, gr, gi = sp.symbols('ar ai br bi gr gi', real=True)
al, be, ga = ar + sp.I*ai, br + sp.I*bi, gr + sp.I*gi
x = [u, 0, w, z]
a = [q*ga, al, (z*be - u*w*ga)/q, (-w*be - u*z*ga)/q]          # alpha e2 + beta G + gamma J
ab2 = lambda c: sp.expand(c*sp.conjugate(c))
m = lambda i, j: x[i-1]*a[j-1] - x[j-1]*a[i-1]
gen = dict(
    A=ab2(m(1, 2)),
    k=ab2(m(1, 2)) + u**2 + w**2 + z**2 - ab2(m(3, 4)) - ab2(m(1, 4)),
    h=1 - 2*w**2 - ab2(a[2]) - z**2 + ab2(m(3, 4)),
    d=u**2 + z**2 - ab2(m(1, 4)),
    e=ab2((1 - w**2)*a[1]) - z**2*ab2(a[1]))
sv, bv, rv = ab2(al), ab2(be), ab2(ga)
Rv = ab2(u*w*be + z*ga)/q**2
cl = dict(
    A=u**2*sv,
    k=u**2*sv + 1 - q**2*bv - Rv,
    h=1 - 2*w**2 - z**2 + q**2*bv - (u**2*w**2 + z**2)*(1 - sv)/q**2 + Rv,
    d=u**2 + z**2 - Rv,
    e=sv*((u**2 + z**2)**2 - z**2))
rel = [u**2 + q**2 - 1, w**2 + z**2 - q**2, sv + bv + rv - 1]
gb = sp.groebner([sp.expand(r) for r in rel], ar, br, gr, u, w, ai, bi, gi, z, q, order='lex')
for key in gen:
    diff = sp.together(sp.expand(gen[key] - cl[key]))
    num = sp.expand(sp.numer(diff))
    check(f'reduction {key}', gb.reduce(num)[1] == 0)
xa = sp.expand(sp.numer(sp.together(sum(xi*ai_ for xi, ai_ in zip(x, a)))))
check('frame: x.a = 0', all(gb.reduce(sp.expand(part))[1] == 0 for part in xa.as_real_imag()))
check('frame: |a|^2 = |alpha|^2+|beta|^2+|gamma|^2',
      gb.reduce(sp.expand(sp.numer(sp.together(sum(ab2(c) for c in a) - (sv + bv + rv)))))[1] == 0)

# 2. Bernstein decomposition. beta = sqrt(1-s) B, gamma = sqrt(1-s) G, |B|^2 + |G|^2 = 1.
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
T = 2*A + (2*k + h)*t + d*t**2
Fs = sp.Poly(sp.expand(F*q2**2), s)
c0, c1, c2 = (Fs.coeff_monomial(s**i) for i in range(3))
F0, Fc = c0, c1 + 2*c0
F1 = red(c2 - F0 + Fc)
check('q^4 F is quadratic in s', Fs.degree() <= 2)
check('Bernstein identity', red(F*q2**2 - ((1 - s)**2*F0 + s*(1 - s)*Fc + s**2*F1)) == 0)
p = (u*z)**2*X11 - 2*u*w*z*X12 + w**2*X22          # q^2 k0 = |uz B - w G|^2
n2 = z**2*X11 - 2*u*w*z*X12 + (u*w)**2*X22         # q^2 d0 = |z B - uw G|^2
check('F0 = q^2 t^2 |uzB-wG|^2 (u^2 + t|zB-uwG|^2/q^2)', red(F0*q2 - q2**2*t**2*p*(u**2*q2*one + t*n2)/q2) == 0)
check('F1 free of (B, G)', not (F1.free_symbols & {X11, X22, X12}))
Fc = red(Fc)
M11, M22, M12 = Fc.coeff(X11), Fc.coeff(X22), Fc.coeff(X12)/2
check('Fc is a Hermitian form in (B, G)', sp.expand(Fc - M11*X11 - M22*X22 - 2*M12*X12) == 0)


def homog(expr):
    """Write a polynomial in t, w, z (even in w, z) in U, W, Z, homogenising with U + W + Z = 1 per power of t."""
    expr = sp.expand(expr.subs({w: sp.sqrt(W), z: sp.sqrt(Z)}))
    out = 0
    for (j,), coef in sp.Poly(expr, t).terms():
        P = sp.Poly(coef, W, Z)
        deg = P.total_degree()
        out += t**j*sum(cf*W**i*Z**l*(U + W + Z)**(deg - i - l) for (i, l), cf in P.terms())
    return sp.Poly(sp.expand(out), t, U, W, Z)


def nonneg(expr):
    return all(cf >= 0 for _, cf in homog(expr).terms())


T1 = red(sp.Poly(sp.expand(T*q2), s).eval(1))
check('F1/q^4 has nonnegative coefficients', nonneg(sp.cancel(F1/q2**2)))
check('T(s=1)/q^2 has nonnegative coefficients', nonneg(sp.cancel(T1/q2)))
T0 = red(sp.Poly(sp.expand(T*q2), s).eval(0))
check('T(s=0) = t(p + u^2 q^2) + t^2 n2 (both PSD forms)', red(T0 - (t*(p + u**2*q2*one) + t**2*n2)) == 0)
check('M11/(t q^2) has nonnegative coefficients', nonneg(sp.cancel(M11/(t*q2))))
check('M22/(t q^2) has nonnegative coefficients', nonneg(sp.cancel(M22/(t*q2))))
detM = homog(sp.cancel(red(M11*M22 - M12**2)/(t**2*q2**4)))
negative = sorted(str(sp.Mul(cf, *[g**e for g, e in zip((t, U, W, Z), mon)])) for mon, cf in detM.terms() if cf < 0)
check('det M/(t^2 q^8) has exactly the negative terms -U W^2 Z t^4, -U W Z^2 t^4',
      negative == sorted(['-U*W**2*Z*t**4', '-U*W*Z**2*t**4']))

# 4. Exact obstruction: Fc < 0 is attained, but F > 0 for every s at that point.
pt = {u: sp.Rational(1, 9), w: sp.Rational(8, 9), z: sp.Rational(4, 9), t: 20}
B, G = sp.Rational(620, 947), sp.Rational(353, 467)
nB = B**2 + G**2
Xs = {X11: B**2/nB, X22: G**2/nB, X12: B*G/nB}
fc, f0, f1 = Fc.subs(pt).subs(Xs), F0.subs(pt).subs(Xs), F1.subs(pt)
dm = (M11*M22 - M12**2).subs(pt)
check('witness on the sphere', sum(v**2 for k_, v in pt.items() if k_ != t) == 1)
check('witness: Fc < 0', fc < 0)
check('witness: det M < 0', dm < 0)
check('witness: F1 > 0 and Fc^2 < 4 F0 F1, so F > 0 for all s', f1 > 0 and f0 >= 0 and fc**2 < 4*f0*f1)
witness = dict(u='1/9', w='8/9', z='4/9', t='20', B='620/947', G='353/467', Fc=str(fc), detM=str(dm),
               discriminant=str(fc**2 - 4*f0*f1))

ok = all(results.values())
print('ALL CHECKS PASSED' if ok else 'SOME CHECKS FAILED')
if '--write-report' in sys.argv:
    out = pathlib.Path(__file__).with_name('report.json')
    out.write_text(json.dumps({'sympy': sp.__version__, 'checks': results, 'witness': witness}, indent=2) + '\n', newline='\n')
sys.exit(0 if ok else 1)
