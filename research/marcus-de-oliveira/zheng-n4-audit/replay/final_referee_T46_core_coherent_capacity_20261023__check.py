"""Exact replay for final_referee_T46_core_coherent_capacity_20261023 (historical T46 face; used for paper T_47,
and through the transport for historical T47 = paper T_65).

Own code only; nothing from the supplement is read or executed. The face's Gram completion relations
themselves are CITED from foundations/so4-atoms-part2b.md (G3: Gram G0 checked against the full linear
system). Here: the group transports, the Gram algebra of Thm 1, the template and projection of Thm 2,
the optimisation formulas, the numeric example and Cor. 3. Exit nonzero iff an EXACT check fails.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'lib'))
import moc4lib as L
import sympy as sp
FAIL = []
def check(name, ok, level='EXACT'):
    print(f"[{level}] {'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok and level == 'EXACT': FAIL.append(name)
FIDX = {tuple(L.F_VALUES[j]): j for j in range(82)}
def g_on_f(g, jdx): return FIDX[tuple(int(v) for v in L.g_act_values(g, L.F_VALUES[jdx]))]
def realisers(src, dst): return [g for g in L.G_ELEMENTS if all(g_on_f(g, s) == t for s, t in zip(src, dst))]
z = lambda e: sp.simplify(e) == 0

def transports():
    for name, src, dst, lit in (('T24 roles -> T46 roles', [59,53,55,67,54], [1,59,6,27,28], ((0,3,1,2),(3,1,0,2),1)),
                                ('T46 roles -> T47 roles', [1,59,6,27,28], [45,5,7,0,35], ((0,2,3,1),(3,1,0,2),0))):
        gs = realisers(src, dst)
        check('%s %s -> %s realised in order by some g in G (#g=%d)' % (name, src, dst, len(gs)), len(gs) > 0)
        check('%s: literal reading of the printed action in moc4lib convention' % name, all(g_on_f(lit, s) == t for s, t in zip(src, dst)), level='INFO')
    check('core roles (p,j,q,b,d) = (f1,f59,f27,f6,f28): four minor, f28 non-minor; all generators of T46 = paper T_47',
          [L.F_IS_MINOR[x] for x in (1,59,27,6)] == [True]*4 and not L.F_IS_MINOR[28] and set([1,59,27,6,28]) <= set(L.generators(frozenset([1,3,5,6,7,8])))
          and L.find_g(L.CONTACT_SETS[47], frozenset([1,3,5,6,7,8])) is not None)
    check('T47 images (45,5,7,0,35) = (c,e,a,b,d) of historical T47 face q f2 + p f9 + j f1 + c f45 + e f5 + d f35 + b f0 + a f7',
          set([45,5,7,0,35]) <= set(L.generators(frozenset([0,4,5,6,9,12]))))

def thm1_2():
    q, b, j, be, t = sp.symbols('q b j beta t', real=True)
    G0 = sp.Matrix([[q, be, -t, 0], [be, b, -t, 0], [-t, -t, j - t**2, 0], [0, 0, 0, 1]])  # Q, B, J0, D
    Qv, Bv, J0, Dv = [sp.Matrix([1 if i == k else 0 for i in range(4)]) for k in range(4)]
    ip = lambda x, y: (x.T*G0*y)[0]
    s = q + b + 2*be; h = 1 + be; Dl = s*j - (s + 4)*t**2; l = q*b - be**2
    lam = (-h*(j - t**2) + 2*t*(t - 1))/Dl; mu = (-2*t*h + s*(t - 1))/Dl
    P = lam*(Qv + Bv) + mu*J0; J = J0 + t*Dv
    thr = ((j - t**2)*h**2 - 4*t*h*(t - 1) + s*(t - 1)**2)/Dl
    check('Thm 1: Gram of (Q+B, J0) is [[s,-2t],[-2t,j-t^2]] with determinant Delta', z(ip(Qv+Bv, Qv+Bv) - s) and z(ip(Qv+Bv, J0) + 2*t) and z(ip(J0, J0) - (j - t**2)) and z(s*(j - t**2) - 4*t**2 - Dl))
    check('Thm 1: P.(Q+B) = -h, P.J0 = t - 1, |P|^2 = threshold (1)', z(ip(P, Qv+Bv) + h) and z(ip(P, J0) - (t - 1)) and z(ip(P, P) - thr))
    rel = [ip(Qv, P) + ip(Qv, Bv) + ip(P, Bv) + 1, ip(Qv, J) + ip(J, Dv), ip(P, J) + ip(J, Bv) + 1, ip(J, Dv) + ip(J, Bv), ip(Dv, Qv), ip(Dv, P), ip(Dv, Bv)]
    check('Thm 1 / Thm 2: the seven displayed completion relations hold (|D|^2 = 1)', all(z(r) for r in rel))
    M3 = G0[:3, :3]
    check('Thm 1: leaf Schur complement: det(3x3)/l = (j - t^2) - t^2 (q + b - 2 beta)/l, i.e. j > (1 + (q+b-2beta)/l) t^2',
          z(M3.det()/l - (j - t**2 - t**2*(q + b - 2*be)/l)))
    Tm = sp.Matrix([[1, sp.Rational(-1,2), 0, 0, 0], [sp.Rational(-1,2), 2, -1, 0, sp.Rational(-1,2)], [0, -1, 1, 0, 0], [0, 0, 0, 1, 0], [0, sp.Rational(-1,2), 0, 0, 1]])
    # order of the template: Q, P, J, D, B (read from the page image)
    Qi, Pi, Ji, Di, Bi = range(5); g = lambda x, y: Tm[x, y]
    rel_t = [g(Qi,Pi) + g(Qi,Bi) + g(Pi,Bi) + g(Di,Di), g(Qi,Ji) + g(Ji,Di), g(Qi,Di), g(Pi,Ji) + g(Ji,Bi) + g(Di,Di), g(Pi,Di), g(Ji,Di) + g(Ji,Bi), g(Di,Bi)]
    check('Thm 2: template (order Q,P,J,D,B) is positive definite, |D|^2 = 1 and satisfies the seven equations',
          Tm.is_positive_definite and Tm[3,3] == 1 and all(r == 0 for r in rel_t))
    u = (q + b - 2*be)/l
    check('u - 4/s = (q - b)^2/(s l) with u = (q+b-2beta)/l', z(u - 4/s - (q - b)**2/(s*l)))
    k = 2 - q - b
    form2 = h**2/s + (s + k*t)**2/(s*(s*j - (s + 4)*t**2))
    check('form (2): threshold = h^2/s + (s + k t)^2/(s{sj - (s+4)t^2})', z(thr - form2))
    check('derivative of (2) = 2(s + kt){kj + (s+4)t}/{sj - (s+4)t^2}^2', z(sp.diff(form2, t) - 2*(s + k*t)*(k*j + (s + 4)*t)/(s*j - (s + 4)*t**2)**2))
    Aq, Bq, Cq, Dq, Eq = s - h**2 - 4*h, 4*h - 2*s, j*h**2 + s, s*j, s + 4
    fr = (Aq*t**2 + Bq*t + Cq)/(Dq - Eq*t**2)
    check('stationary quadratic: d/dt (At^2+Bt+C)/(D-Et^2) = (B E t^2 + 2(AD + EC) t + B D)/(D-Et^2)^2 with the displayed A..E',
          z(sp.diff(fr, t) - (Bq*Eq*t**2 + 2*(Aq*Dq + Eq*Cq)*t + Bq*Dq)/(Dq - Eq*t**2)**2))
    check('(1) numerator equals A t^2 + B t + C with the displayed A, B, C', z(((j - t**2)*h**2 - 4*t*h*(t - 1) + s*(t - 1)**2) - ((s - h**2 - 4*h)*t**2 + (4*h - 2*s)*t + (j*h**2 + s))))
    ts = -k*j/(s + 4)
    check('value at t* = -kj/(s+4) = 1/j + (s + 2k)/4 + k^2/(4(s+4))', z(thr.subs(t, ts) - (1/j + (s + 2*k)/4 + k**2/(4*(s + 4)))))
    check('beta = 1-q-b (s = k): value = 1/j + k(k+3)/(k+4); the t=0 value is 1/j + k', z((1/j + (s + 2*k)/4 + k**2/(4*(s + 4))).subs(be, 1 - q - b) - (1/j + k*(k + 3)/(k + 4))) and z(thr.subs({t: 0, be: 1 - q - b}) - (1/j + k)))
    check('t = 0 threshold = 1/j + (1+beta)^2/s; d/dt at 0 = 2(2-q-b)/(s j)', z(thr.subs(t, 0) - (1/j + h**2/s)) and z(sp.diff(thr, t).subs(t, 0) - 2*(2 - q - b)/(s*j)))
    S = q + b
    check('Cor 3: at t = 0, beta = 1 - S: s = h = 2 - S and threshold = 1/j + 2 - S', z(s.subs(be, 1 - S) - (2 - S)) and z(thr.subs({t: 0, be: 1 - S}) - (1/j + 2 - S)))
    check('Cor 3: d/dbeta (1+beta)^2/(S+2beta) = 2(1+beta)(S+beta-1)/(S+2beta)^2', z(sp.diff(h**2/(S + 2*be), be) - 2*(1 + be)*(S + be - 1)/(S + 2*be)**2))
    qq, bb = sp.symbols('qq bb', positive=True)
    check('Cor 3: discriminant in q of (q + b - 1)^2 - q b is b(4 - 3b)', z(sp.discriminant(sp.expand((qq + bb - 1)**2 - qq*bb), qq) + bb*(4 - 3*bb)) or z(sp.discriminant(sp.expand((qq + bb - 1)**2 - qq*bb), qq) - bb*(4 - 3*bb)))
    # numeric example
    qv, bv, jv, bev, tv = sp.Rational(1, 8), 4, 4, sp.Rational(-1, 2), sp.Rational(1, 4)
    sub = {q: qv, b: bv, j: jv, be: bev, t: tv}
    check('example q=1/8, b=4, j=4, beta=-1/2, t=1/4: l>0, Schur condition, threshold = 399/1543 < 13/50 < 1/4 + (1-1/sqrt2)^2/(33/8-sqrt2)',
          l.subs(sub) > 0 and (j - (1 + (q + b - 2*be)/l)*t**2).subs(sub) > 0 and thr.subs(sub) == sp.Rational(399, 1543)
          and sp.Rational(399, 1543) < sp.Rational(13, 50) < sp.Rational(1, 4) + (1 - 1/sp.sqrt(2))**2/(sp.Rational(33, 8) - sp.sqrt(2)))
    check('the "preceding parallel" value equals the t = 0 optimum over beta at these q, b, j (1/j + min_beta (1+beta)^2/(q+b+2beta))',
          sp.simplify(sp.Rational(1,4) + (1 - 1/sp.sqrt(2))**2/(sp.Rational(33,8) - sp.sqrt(2)) - (1/jv + (h**2/s).subs({q: qv, b: bv, be: -1/sp.sqrt(2)}))) == 0, level='INFO')

if __name__ == '__main__':
    transports(); thm1_2()
    print('summary: %d EXACT failures' % len(FAIL)); sys.exit(1 if FAIL else 0)
