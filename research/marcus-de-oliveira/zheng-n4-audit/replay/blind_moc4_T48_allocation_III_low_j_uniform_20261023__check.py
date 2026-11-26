# Auditor's own replay code (Zheng n=4 audit). Nothing from the supplement is executed or imported;
# supplement JSON polynomial data is read with json.load only. Supplement root: env MOC4_SUPP or argv[1].
import sys, os, json, random, itertools
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get('AUDLIB', os.path.join(HERE, '..', 'lib')))
import moc4lib as L
import sympy as sp
SUPP = os.environ.get('MOC4_SUPP') or (sys.argv[1] if len(sys.argv) > 1 else None)
PDIR = os.path.join(SUPP, 'proofs', 'manin_transfer_20261021') if SUPP else None
FAIL = []
def chk(name, ok, info=''):
    print(('PASS ' if ok else 'FAIL ') + name + ((' :: ' + str(info)) if info else ''), flush=True)
    if not ok: FAIL.append(name)
def info(*a): print('INFO', *a, flush=True)
MS = sp.symbols('m0:14'); CUB = sp.expand(L.cubic_C(MS))
HESS_LIN = [[[int(sp.diff(CUB, MS[i], MS[j], MS[k])) for k in range(14)] for j in range(14)] for i in range(14)]
def hess_at(m):
    return [[sum(HESS_LIN[i][j][k] * m[k] for k in range(14) if HESS_LIN[i][j][k]) for j in range(14)] for i in range(14)]
def source_vec(roles, coefs, h=0):
    w = [Fr(0)] * 14
    for r, c in zip(roles, coefs):
        for k in range(14): w[k] += c * L.F_VECS[r][k]
    return [w[k] + h * L.KAPPA[k] for k in range(14)]
def sym_source(roles, syms, hsym):
    w = [0] * 14
    for r, c in zip(roles, syms):
        for k in range(14): w[k] += c * sp.Rational(L.F_VECS[r][k])
    return [sp.expand(w[k] + hsym * sp.Rational(L.KAPPA[k])) for k in range(14)]
def load_radial(name):
    d = json.load(open(os.path.join(PDIR, name)))
    vs = sp.symbols(' '.join(d['variables']))
    poly = lambda terms: sp.Add(*[sp.Integer(int(c)) * sp.Mul(*[v**e for v, e in zip(vs, ex)]) for ex, c in terms])
    facs = [(poly(f['terms']), f['multiplicity'], f['degree'], len(f['terms'])) for f in d['radial_factors']]
    return vs, poly(d['C_terms']), poly(d['Q_terms']), facs, int(d['constant'])
def evalP(terms, vals):
    tot = Fr(0)
    for mon, cf in terms:
        t_ = Fr(int(cf))
        for v, ex in zip(vals, mon):
            if ex: t_ *= v**ex
        tot += t_
    return tot
def finish():
    print('SUMMARY: %d EXACT check(s) failed' % len(FAIL) if FAIL else 'SUMMARY: all EXACT checks passed', flush=True)
    sys.exit(1 if FAIL else 0)
q,p,j,c,e,d,b,a,h = sp.symbols('q p j c e d b a h'); syms = [q,p,j,c,e,d,b,a]
roles = [7,5,29,39,2,40,41,0]; names = 'q p j c e d b a'.split()
FV = L.F_VALUES
T = L.CONTACT_SETS[49]
STAGE6 = {}
for st, ks in {1:[1,2,3,4,5,6,7,15,18,19,20,21,22,23,25,27,28,29,31,34,35,36,37,39,40,42,43,44,72,88],2:[10,26,46,90],3:[9,11,12,47],
               4:[8,13,17,49],5:[14,16,24,33,38,70],6:[52,57,58,65,68,71,76,78],7:[45,51,66,69,75,81],8:[48,50,53,59,62,67,73,79,86],
               9:[30,41,54,55,60,61,63,74,77,80,83,84],10:[56,64,82,85,87],11:[32],12:[89]}.items():
    for k in ks: STAGE6[k] = st
def classify6(S6):
    S6 = frozenset(S6)
    if all(L.F_IS_MINOR[g] for g in L.generators(S6)): return ('minor', 0)
    for k, v in L.CONTACT_SETS.items():
        if L.find_g(v, S6) is not None: return ('T%d' % k, STAGE6[k])
    return ('?', None)
chk('tuple (0,3,5,6,9,12) is literally paper T_49', T == frozenset({0,3,5,6,9,12}))
chk('generators(T_49) = {0,1,2,5,7,9,29,39,40,41}, only f40 non-minor', L.generators(T) == [0,1,2,5,7,9,29,39,40,41] and [k for k in L.generators(T) if not L.F_IS_MINOR[k]] == [40])
add = lambda *ks: tuple(sum(FV[k][s] for k in ks) for s in range(24))
chk('f40 + f1 = f7 + f39 + f0 (all 24 values; P injective)', add(40, 1) == add(7, 39, 0))
chk('f40 + f9 = f7 + f2 + f41', add(40, 9) == add(7, 2, 41))
chk('canonical eight profiles (7,5,29,39,2,40,41,0) are independent', L.rank([list(L.F_VECS[k]) for k in roles]) == 8)
for i, k in enumerate(roles):
    if k == 40: continue
    rest = [kk for kk in roles if kk != k]
    Zs = set.intersection(*[set(L.F_ZERO[kk]) for kk in rest])
    subs = sorted(set(classify6(S6) for S6 in itertools.combinations(sorted(Zs), 6) if frozenset(S6) != T))
    early = [sx for sx in subs if sx[1] is not None and sx[1] < 4]
    info('zero coefficient %s (f%d): extra contacts %s; six-subsets %s' % (names[i], k, sorted(Zs - set(T)), subs))
    if k != 39: chk('zero face %s reduces to an earlier-stage six-set (stage < 4)' % names[i], bool(early), early)
chk('c = 0 face <-> extra contact 17, and T_49+{17} is Table 5 set 23', sorted(set.intersection(*[set(L.F_ZERO[kk]) for kk in roles if kk != 39]) - set(T)) == [17] and frozenset(T | {17}) == L.AUX7[23])
w = sym_source(roles, syms, h)
F = sp.Poly(sp.expand(L.cubic_C(w)), h)
chk('F(g + h kappa) affine in h', F.degree() == 1)
C = sp.expand(F.coeff_monomial(1)); Q = sp.expand(F.coeff_monomial(h))
info('own Q =', Q)
vs, Cj, Qj, facs, const = load_radial('blind_moc4_T48_native_full_radial_20261023.json')
sub = dict(zip(vs[:8], syms))
chk('supplement JSON C,Q equal own C,Q', sp.expand(Cj.subs(sub) - C) == 0 and sp.expand(Qj.subs(sub) - Q) == 0)
G = [f for f, m, deg, n in facs if deg == 9][0].subs(sub)
GP = sp.Poly(G, *syms); GT = GP.terms()
chk('G (supplement data) has 1772 terms, degree 9', len(GT) == 1772 and GP.total_degree() == 9)
rnd = random.Random(5); ok = True; nn = 0
fC = sp.lambdify(syms, C); fQ = sp.lambdify(syms, Q)
for it in range(20):
    vals = [Fr(rnd.randint(-20, 30), rnd.randint(1, 7)) for _ in range(8)]
    if it % 3 == 0: vals[3] = Fr(0)
    Cv, Qv = Fr(fC(*vals)), Fr(fQ(*vals))
    if Qv == 0: continue
    dv = L.det(hess_at(source_vec(roles, vals, -Cv / Qv))); nn += 1
    ok &= (dv * Qv**2 == -4 * evalP(GT, vals)**2)
chk('det H(g - C kappa/Q) * Q^2 == -4 G^2 at %d random rational points (own Hessian; pointwise exact)' % nn, ok)
Q0 = sp.expand(Q.subs(c, 0)); C0 = sp.expand(C.subs(c, 0))
# ===== blind_moc4_T48_allocation_III_low_j_uniform_20261023 =====
Y42 = sp.Poly(sp.expand((e+j)*C - e*j*Q), *syms)
chk('(e+j) C - e j Q has only positive coefficients (forced-height exclusion)', all(cf > 0 for _, cf in Y42.terms()), '%d terms' % len(Y42.terms()))
Qs, js, es, as_ = sp.symbols('Q j e a', positive=True)
u = (1 + 2/es) / (2*(1/Qs + 1/js + 1/es)); v = 1 - u
chk('stationarity: u(1/Q+1/j) - v/e = 1/2', sp.simplify(u*(1/Qs + 1/js) - v/es - sp.Rational(1, 2)) == 0)
pst = u**2*(1/Qs + 1/js) + v - as_/4 + v**2/es
chk('p* = u^2(1/Q+1/j) + v - a/4 + v^2/e = 1 - u/2 + v/e - a/4', sp.simplify(pst - (1 - u/2 + v/es - as_/4)) == 0)
A0s, Ls = sp.symbols('A0 L', positive=True)
uL = Ls/(1+Ls); vL = 1/(1+Ls); QA = 2*uL*(1+A0s)
jform = 2*uL*es*(1+A0s) / (A0s*es + 2*vL*(1+A0s))
chk('stationarity with Q = 2u(1+A0) <=> j = 2ue(1+A0)/(A0 e + 2v(1+A0))', sp.simplify(uL*(1/QA + 1/jform) - vL/es - sp.Rational(1, 2)) == 0)
X_, V_, Aa = sp.symbols('X V A', positive=True)
e2 = sp.solve(sp.Eq(jform.subs({Ls: 1+X_, A0s: Aa}), V_/(1+V_)), es)
chk('chart 2: e = 2(1+A)V / (2(1+X)(1+A) + V[2(1+X)+AX])', len(e2) == 1 and sp.simplify(e2[0] - 2*(1+Aa)*V_/(2*(1+X_)*(1+Aa) + V_*(2*(1+X_) + Aa*X_))) == 0)
e3 = sp.solve(sp.Eq(jform.subs({Ls: X_/(1+X_), A0s: 2*X_*Aa/(1+Aa)}), V_/(1+V_)), es)
R_ = 1 + Aa*(1 + 2*X_)
chk('chart 3: e = (1+X) R V / (X(R+V))', len(e3) == 1 and sp.simplify(e3[0] - (1+X_)*R_*V_/(X_*(R_+V_))) == 0)
nn_, dn_ = sp.fraction(sp.together(jform.subs({Ls: X_/(1+X_), A0s: 2*X_+Aa, es: V_})))
dmn = sp.factor(sp.expand(dn_ - nn_)); info('chart 1: denominator - numerator of j =', dmn)
chk('chart 1: j < 1 (denominator - numerator is a positive-coefficient polynomial)', all(cf > 0 for _, cf in sp.Poly(sp.expand(dn_ - nn_), X_, Aa, V_).terms()))
Kx = lambda Lv, A0v: 2*Lv + A0v*(Lv - 1)
chk('chart 1 bijects onto K <= 0: K(X/(1+X), 2X+A) = -A/(1+X)', sp.simplify(Kx(X_/(1+X_), 2*X_+Aa) + Aa/(1+X_)) == 0)
chk('chart 3 lies in K > 0: K(X/(1+X), 2XA/(1+A)) = 2X/((1+X)(1+A))', sp.simplify(Kx(X_/(1+X_), 2*X_*Aa/(1+Aa)) - 2*X_/((1+X_)*(1+Aa))) == 0)
# atoms used by the pointwise bound (orbit check, roles as printed)
fidx = {FV[k]: k for k in range(82)}
def gf(gm, k):
    vals = [0]*24
    for s_ in range(24): vals[gm[s_]] = FV[k][s_]
    return fidx[tuple(vals)]
chainO = set(); oldH = set()
for gm in L._G_MAPS:
    cc_, l1, l2, kk = [gf(gm, k) for k in (2, 0, 12, 70)]; chainO.add((cc_,) + tuple(sorted((l1, l2))) + (kk,))
    oldH.add(tuple(gf(gm, k) for k in (3, 15, 26)))
chk('chain orbit 576, old-H orbit 288', len(chainO) == 576 and len(oldH) == 288)
for t_ in [(5, 7, 29, 40), (5, 29, 41, 40), (5, 0, 2, 40)]:
    chk('chain (centre f%d; leaves f%d,f%d; non-minor f%d) in the proved chain orbit' % t_, t_ in chainO)
chk('old-H (aux f2 = e, main f5 = p; f40) in the proved old-H orbit', (2, 5, 40) in oldH)
chk('reversed old-H roles (aux f5, main f2) NOT in orbit (roles meaningful)', (5, 2, 40) not in oldH)
# pointwise allocation algebra: e shared in ratio h : a/2, (q,j),(b,j) chains with currents uq/Q, ub/Q
qv, bv, jv, ev, av, uv = sp.symbols('q b j e a u', positive=True); Qv = qv + bv; vv = 1 - uv; hv = vv - av/2
chainQ = (uv*qv/Qv)**2/qv + (uv*qv/Qv)**2/(jv*qv/Qv) + (uv*bv/Qv)**2/bv + (uv*bv/Qv)**2/(jv*bv/Qv)
side_e = hv + hv**2/(ev*hv/vv) + (av/2)**2/av + (av/2)**2/(ev*(av/2)/vv)
chk('allocation algebra: chain requirement = u^2(1/Q+1/j); old-H + (a,e) chain = v - a/4 + v^2/e', sp.simplify(chainQ - uv**2*(1/Qv + 1/jv)) == 0 and sp.simplify(side_e - (vv - av/4 + vv**2/ev)) == 0)
# NUMERIC/pointwise: sign of G on random points of the three charts (d = 1), exact rational evaluation
def point(chart, X0, A, B, C1, V, P, cval, WV=Fr(0)):
    if chart == 1: Lv = X0/(1+X0); A0 = 2*X0 + A
    elif chart == 2: Lv = 1 + X0; A0 = A
    else: Lv = X0/(1+X0); A0 = 2*X0*A/(1+A)
    u_ = Lv/(1+Lv); v_ = 1/(1+Lv); Q_ = 2*u_*(1+A0)
    if chart == 1: e_ = V; j_ = 2*u_*e_*(1+A0)/(A0*e_ + 2*v_*(1+A0))
    else:
        j_ = V/(V+WV)
        e_ = j_*2*v_*(1+A0)/(2*u_*(1+A0) - A0*j_)
    a_ = 2*v_*B/(1+B); q_ = Q_*C1/(1+C1); b_ = Q_/(1+C1)
    ps = 1 - u_/2 + v_/e_ - a_/4; p_ = ps*P/(1+P)
    return [q_, p_, j_, cval, e_, Fr(1), b_, a_], (u_, v_, Q_, j_, e_, a_, ps)
rnd = random.Random(23); cnt = {1: [0, 0], 2: [0, 0], 3: [0, 0]}; worst = None; badhyp = 0
rq = lambda: Fr(rnd.randint(1, 400), rnd.choice([10, 100, 1000]))
for it in range(1800):
    ch = 1 + it % 3
    WV = Fr(0) if (ch > 1 and rnd.random() < 0.15) else Fr(1)
    V = rq(); WVv = Fr(0) if WV == 0 else rq()
    X0 = rq(); A = rq() if rnd.random() > 0.15 else Fr(0)
    cval = rq() if rnd.random() > 0.2 else Fr(0)
    vals, (u_, v_, Q_, j_, e_, a_, ps) = point(ch, X0, A, rq(), rq(), V, rq(), cval, WVv if ch > 1 else Fr(0))
    hyp = (0 < u_ < 1) and Q_ >= 2*u_ and a_ < 2*v_ and j_ <= 1 and e_ > 0 and vals[1] < ps
    stat = u_*(1/Q_ + 1/j_) - v_/e_ == Fr(1, 2)
    if not (hyp and stat): badhyp += 1; continue
    g = evalP(GT, vals); cnt[ch][0] += 1; cnt[ch][1] += (g < 0)
    if worst is None or g > worst[0]: worst = (g, ch)
chk('chart points satisfy the Theorem-1 hypotheses and stationarity (sample)', badhyp == 0, badhyp)
info('pointwise exact evaluation (seed 23): chart -> [points, #G<0] = %s; max G = %.3e (chart %d)' % (cnt, float(worst[0]), worst[1]))
finish()
