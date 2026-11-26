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
chk('Q|c=0 = Q0 = d(d+e+j) - p(e+j)', sp.expand(Q0 - (d*(d+e+j) - p*(e+j))) == 0)
X = sp.expand(sp.cancel((e+j) * C0.subs(p, d*(d+e+j)/(e+j))))
PX = sp.Poly(X, *syms)
chk('(e+j) C0 at Q0 = 0 has only positive coefficients, with a d^3 j term', all(cf > 0 for _, cf in PX.terms()) and any(m == (0,0,1,0,0,3,0,0) for m, _ in PX.terms()), '%d terms' % len(PX.terms()))
G0 = sp.Poly(G.subs(c, 0), *syms)
chk('G|c=0 has 759 terms', len(G0.terms()) == 759, len(G0.terms()))
# Lemma 2 (radial necessity), symbolic sanity: F(x+y)=0 and H(x+y)(x-y)=0 for x,y in X (random Gaussian-rational spectra)
def mv(aa, bb): return [complex(v) for v in L.matching_vector_np(aa, bb, 1.0)]
import numpy as np
rng = np.random.default_rng(2); worst = 0
HLn = np.array(HESS_LIN, float)
for it in range(50):
    x1 = np.array(mv(rng.normal(size=4)+1j*rng.normal(size=4), rng.normal(size=4)+1j*rng.normal(size=4)))*np.exp(1j*rng.uniform(0,6.3))
    wv = x1.real; y1 = x1.conj()/2; x2 = x1/2
    r1 = abs(L.cubic_C(list(wv))); r2 = np.max(np.abs((HLn @ wv) @ (x2 - y1)))
    worst = max(worst, r1 / (1 + np.max(np.abs(wv)))**3, r2 / (1 + np.max(np.abs(wv)))**2)
info('NUMERIC Lemma 2 sanity (seed 2, 50 random realizable w = Re(alpha m)): max scaled |F(w)|, |H(w)(x-y)| = %.1e' % worst)
finish()
