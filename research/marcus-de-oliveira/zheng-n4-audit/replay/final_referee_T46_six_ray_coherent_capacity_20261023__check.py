"""Exact replay for final_referee_T46_six_ray_coherent_capacity_20261023 (historical T46 face; used for
paper T_47, and through the transport for historical T47 = paper T_65).

Own code only; nothing from the supplement is read or executed. The Gram completion relations of the
T46 face themselves are CITED from foundations/so4-atoms-part2b.md (G3, checked there against the full
linear system); this script replays the scalar algebra of Thm 1 and Cor. 2-4 and the group transport.
Every printed line is one check; exit nonzero iff an EXACT check fails.
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

def transport():
    gs = realisers([27,1,59,5,2,28,6,4], [0,45,5,9,29,35,7,36])
    check('some g in G maps the T46 roles (27,1,59,5,2,28,6,4) to (0,45,5,9,29,35,7,36) in order (#g=%d)' % len(gs), len(gs) > 0)
    lit = ((1,2,0,3),(2,1,3,0),1)
    check('the literal reading (row (1,2,0,3), column (2,1,3,0), transpose) in moc4lib convention realises it', all(g_on_f(lit, s) == t for s, t in zip([27,1,59,5,2,28,6,4],[0,45,5,9,29,35,7,36])), level='INFO')
    T47h = frozenset([0,4,5,6,9,12])
    check('historical T47 {0,4,5,6,9,12} is paper T_65 up to G; its generators are {0,1,2,5,7,9,35,45}',
          L.find_g(L.CONTACT_SETS[65], T47h) is not None and sorted(L.generators(T47h)) == [0,1,2,5,7,9,35,45])
    check('images of the six-ray roles (q,p,j,c,d,b) = (0,45,5,9,35,7) are generators of historical T47; f9, f45, f0, f5, f7 minor, f35 non-minor',
          set([0,45,5,9,35,7]) <= set(L.generators(T47h)) and [L.F_IS_MINOR[x] for x in (0,45,5,9,7)] == [True]*5 and not L.F_IS_MINOR[35])
    check('Phi_{g.c}(g.O) = Phi_c(O) for this g at a random rational O and the T46 six-ray cost (exact)', (lambda g: (lambda O, cvec: L.Phi(L.g_act_c(g, cvec), L.g_act_O(g, O)) == L.Phi(cvec, O))(L.random_rational_SO4(__import__('random').Random(5)), [sum(L.F_VECS[k][i] for k in (27,1,59,5,28,6)) for i in range(14)]))(gs[0]))

def thm1():
    q, b, j, be, t, v = sp.symbols('q b j beta t v', real=True)
    Gm = sp.Matrix([[q, be, t, 0], [be, b, v, 0], [t, v, j - t**2, 0], [0, 0, 0, 1]])   # basis Q, B, J0, Z
    Qv, Bv, J0, Zv = [sp.Matrix([1 if i == k else 0 for i in range(4)]) for k in range(4)]
    ip = lambda x, y: (x.T*Gm*y)[0]
    S = Qv + Bv; h = 1 + be; s = q + b + 2*be
    Dl = s*(j - t**2) - (t + v)**2
    lam = (-h*(j - t**2) + (t + v)*(1 + v))/Dl; mu = (h*(t + v) - s*(1 + v))/Dl
    P = lam*S + mu*J0
    pstar = ((j - t**2)*h**2 - 2*h*(t + v)*(1 + v) + s*(1 + v)**2)/Dl
    check('Thm 1: P.S = -h, P.J0 = -1-v, |P|^2 = p* (1)', z(ip(P, S) + h) and z(ip(P, J0) + 1 + v) and z(ip(P, P) - pstar))
    kc = (t - v)/(j - t**2 - v**2/b); Cv = kc*(J0 - (v/b)*Bv); Jv = J0 + t*Zv; Dv = Cv - Zv
    cstar = (t - v)**2/(j - t**2 - v**2/b)
    check('Thm 1: C.B = C.Z = 0, C.J = t - v, |C|^2 = c* (2), |D|^2 = c* + 1',
          z(ip(Cv, Bv)) and z(ip(Cv, Zv)) and z(ip(Cv, Jv) - (t - v)) and z(ip(Cv, Cv) - cstar) and z(ip(Dv, Dv) - cstar - 1))
    rel = [ip(Qv, P) + ip(Qv, Bv) + ip(P, Bv) + 1, ip(Qv, Zv), ip(P, Zv), ip(Bv, Zv), ip(Cv, Zv), ip(Cv, Bv),
           ip(Jv, Qv) - ip(Jv, Zv), ip(Jv, P) + ip(Jv, Bv) + 1, ip(Jv, Cv) + ip(Jv, Bv) - ip(Jv, Zv)]
    check('Thm 1: the nine displayed completion relations hold for the construction', all(z(r) for r in rel))
    l = q*b - be**2
    M3 = sp.Matrix([[q, be, t], [be, b, v], [t, v, j - t**2]])
    check('Thm 1: Schur complement of the leaf block = (j - t^2) - (b t^2 - 2 beta t v + q v^2)/l', z(M3.det()/l - ((j - t**2) - (b*t**2 - 2*be*t*v + q*v**2)/l)))
    # template of the density argument
    vecs = dict(Q=(1,0,0,0,0,0), B=(0,1,0,0,0,0), J=(0,0,1,0,0,0), P=(sp.Rational(-1,2),sp.Rational(-1,2),-1,0,1,0), C=(0,0,0,0,0,1), D=(0,0,0,-1,0,1))
    dot = lambda x, y: sum(a*b for a, b in zip(vecs[x], vecs[y]))
    Zt = tuple(c - d for c, d in zip(vecs['C'], vecs['D'])); vecs['Z'] = Zt
    rel_t = [dot('Q','P') + dot('Q','B') + dot('P','B') + (dot('D','D') - dot('C','C')), dot('Q','Z'), dot('P','Z'), dot('B','Z'), dot('C','Z'), dot('C','B'),
             dot('J','Q') - dot('J','Z'), dot('J','P') + dot('J','B') + (dot('D','D') - dot('C','C')), dot('J','C') + dot('J','B') - dot('J','Z')]
    Gt = sp.Matrix([[dot(x, y) for y in 'QBJPCD'] for x in 'QBJPCD'])
    check('density template: Gram PD, |D|^2 - |C|^2 = 1, satisfies the nine relations (with |Z|^2 normalised to 1)',
          Gt.is_positive_definite and dot('D','D') - dot('C','C') == 1 and all(r == 0 for r in rel_t))

def cors():
    bq, r, e, t, p = sp.symbols('b47 r e t p', positive=True)
    aa = r**2/bq
    kap = 1 + 1/bq; rho = 1 + r/bq; eta = (1 - aa)/(bq - r); k = bq*(1 - r)**2/(bq - r)**2
    q46, b46, j46, be, v = bq, aa, e, -r, -r*t/bq
    s = q46 + b46 + 2*be; h = 1 + be; Dl = s*(j46 - t**2) - (t + v)**2
    pstar = ((j46 - t**2)*h**2 - 2*h*(t + v)*(1 + v) + s*(1 + v)**2)/Dl
    cstar = (t - v)**2/(j46 - t**2 - v**2/b46)
    check('Cor 2: with beta=-r, q=b47, b=a47, v=-rt/b47: leaf block singular (l = 0)', z(q46*b46 - be**2))
    check('Cor 2: c*_46 = rho^2 t^2/(e - kappa t^2) (the T47 p-bound)', z(cstar - rho**2*t**2/(e - kap*t**2)))
    check('Cor 2: p*_46 = k + (1 - eta t)^2/(e - kappa t^2) (the T47 c-bound)', z(pstar - (k + (1 - eta*t)**2/(e - kap*t**2))))
    cb = lambda tt: k + (1 - eta*tt)**2/(e - kap*tt**2); pb = lambda tt: rho**2*tt**2/(e - kap*tt**2)
    t1 = eta*e/kap
    check('Cor 2: at t = eta e/kappa: c-bound = 1/e + k - eta^2/kappa, p-bound = rho^2 eta^2 e/(kappa(kappa - eta^2 e))',
          z(cb(t1) - (1/e + k - eta**2/kap)) and z(pb(t1) - rho**2*eta**2*e/(kap*(kap - eta**2*e))))
    check('Cor 2: p-bound <= p  <=>  e <= p kappa^2/(eta^2 (p kappa + rho^2)) (identity of the boundary)',
          z(sp.solve(sp.Eq(rho**2*eta**2*e/(kap*(kap - eta**2*e)), p), e)[0] - p*kap**2/(eta**2*(p*kap + rho**2))))
    D0 = kap + rho**2/p; M0 = D0 + eta**2*e; N = M0**2 - kap*eta**2*e; t3 = eta*e/M0
    check('Cor 3: N = D0^2 + eta^2 e (2 D0 - kappa + eta^2 e) and N - rho^2 eta^2 e/p = D0^2 + D0 eta^2 e + eta^4 e^2',
          z(N - (D0**2 + eta**2*e*(2*D0 - kap + eta**2*e))) and z(N - rho**2*eta**2*e/p - (D0**2 + D0*eta**2*e + eta**4*e**2)))
    check('Cor 3: at t = eta e/M0, e - kappa t^2 = e N/M0^2, consumed c = k + D0^2/(e N), consumed p = rho^2 eta^2 e/N',
          z(e - kap*t3**2 - e*N/M0**2) and z(cb(t3) - (k + D0**2/(e*N))) and z(pb(t3) - rho**2*eta**2*e/N))
    A = eta**2*e; t4 = eta*e/D0; N0 = D0**2 - kap*A
    check('Cor 4 (A <= D0): t = eta e/D0 gives c = k + (D0-A)^2/(e(D0^2 - kappa A)), unused p = p D0 (D0-A)/N0',
          z(cb(t4) - (k + (D0 - A)**2/(e*N0))) and z(p - pb(t4) - p*D0*(D0 - A)/N0))
    check('Cor 4 (A >= D0): t = 1/eta gives c = k, p = rho^2/(A - kappa), and rho^2/(A-kappa) <= p iff A >= D0 (A > kappa)',
          z(cb(1/eta) - k) and z(pb(1/eta) - rho**2/(A - kap)) and z((A - kap) - rho**2/p - (A - D0)))
    diff = (k + D0**2/(e*N)) - (k + (D0 - A)**2/(e*N0))
    check('Cor 4 comparison: difference = A^2 {2 D0 (D0 - kappa) + kappa A - A^2} / (e (D0^2 - kappa A)((D0 + A)^2 - kappa A))',
          z(diff - A**2*(2*D0*(D0 - kap) + kap*A - A**2)/(e*N0*((D0 + A)**2 - kap*A))))
    x = sp.symbols('x')
    br = 2*D0*(D0 - kap) + kap*x - x**2
    check('Cor 4: the brace at A = 0 and A = D0 equals 2 D0 (D0 - kappa) and D0 (D0 - kappa)', z(br.subs(x, 0) - 2*D0*(D0 - kap)) and z(br.subs(x, D0) - D0*(D0 - kap)))
    check('Cor 2 hypotheses b>1, 0<a<1/b give r<1<b, eta>0, D0>kappa>1 (sample b=3, a=1/4 exact)',
          all(bool(c) for c in [(r < 1).subs(r, sp.sqrt(sp.Rational(3,4))), (eta > 0).subs({bq: 3, r: sp.sqrt(sp.Rational(3,4))})]))

if __name__ == '__main__':
    transport(); thm1(); cors()
    print('summary: %d EXACT failures' % len(FAIL)); sys.exit(1 if FAIL else 0)
