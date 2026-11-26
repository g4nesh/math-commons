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
# ===== root_moc4_T47_zero_q_rational_capacity_20261023 =====
q,p,j,c,e,d,b,a = sp.symbols('q p j c e d b a'); syms = [q,p,j,c,e,d,b,a]
rr, pp, ee, bb, aa = sp.symbols('r p e b a', positive=True)
k_ = bb*(1-rr)**2/(bb-rr)**2; kap = 1 + 1/bb; rho = 1 + rr/bb; eta = (1-aa)/(bb-rr)
D0 = kap + rho**2/pp; N = D0**2 + eta**2*ee*(2*D0 - kap + eta**2*ee); Crat = k_ + D0**2/(ee*N)
tt = eta*ee/(D0 + eta**2*ee)
cst = k_ + (1-eta*tt)**2/(ee - kap*tt**2); pst = rho**2*tt**2/(ee - kap*tt**2)
chk('c*(t) = Crat at t = eta e/(D0 + eta^2 e)', sp.simplify(cst - Crat) == 0)
chk('e - kappa t^2 = e N/(D0+eta^2 e)^2 (positive denominator)', sp.simplify(ee - kap*tt**2 - ee*N/(D0+eta**2*ee)**2) == 0)
chk('N - (rho^2/p) eta^2 e = D0^2 + D0 eta^2 e + eta^4 e^2', sp.simplify(N - (rho**2/pp)*eta**2*ee - (D0**2 + D0*eta**2*ee + eta**4*ee**2)) == 0)
chk('p - p* = p (N - rho^2 eta^2 e/p)/N  (so p* < p)', sp.simplify(pp - pst - (N - rho**2*eta**2*ee/pp)*pp/N) == 0)
chk('N - D0^2 = eta^2 e (kappa + 2 rho^2/p + eta^2 e) > 0', sp.simplify(N - D0**2 - eta**2*ee*(kap + 2*rho**2/pp + eta**2*ee)) == 0)
vs, Cj, Qj, facs, const = load_radial('root_moc4_T47_full_radial_20261022.json')
sub = dict(zip(vs[:8], syms))
G = [f for f, m, deg, n in facs if deg == 8][0].subs(sub)
GT = sp.Poly(G, *syms).terms()
G0 = sp.Poly(G.subs({q: 0, d: 1}), *syms)
info('G(q=0,d=1) has %d terms (document: 431 "original terms")' % len(G0.terms()))
# the identity det H * Q^2 = -4 d^2 G^2 that ties this G to the cubic is replayed in root_moc4_T47_negative_Phi_20261022__check.py
def crat(b_, r_, p_, e_, a_):
    k = b_*(1-r_)**2/(b_-r_)**2; ka = 1 + 1/b_; rh = 1 + r_/b_; et = (1-a_)/(b_-r_)
    D0v = ka + rh**2/p_; Nv = D0v**2 + et**2*e_*(2*D0v - ka + et**2*e_); return k + D0v**2/(e_*Nv)
rnd = random.Random(11); worst = None; n = walls = 0
for it in range(2500):
    b_ = 1 + Fr(rnd.randint(1, 400), rnd.choice([10, 50, 200])) if rnd.random() > 0.1 else 1 + Fr(1, 10**rnd.randint(2, 5))
    r_ = Fr(rnd.randint(1, 999), 1000) if rnd.random() > 0.1 else Fr(rnd.choice([1, 999999]), 10**6)
    a_ = r_**2/b_
    p_ = b_ + (Fr(rnd.randint(0, 300), rnd.choice([10, 100])) if rnd.random() > 0.2 else 0)
    j_ = Fr(rnd.randint(1, 1000), 1000) if rnd.random() > 0.15 else Fr(1)
    if j_ >= b_: continue
    e0 = b_*j_/(b_-j_); e_ = e0 + (Fr(rnd.randint(0, 500), rnd.choice([10, 100])) if rnd.random() > 0.2 else 0)
    C_ = crat(b_, r_, p_, e_, a_)
    c_ = C_*Fr(rnd.randint(0, 999), 1000) if rnd.random() > 0.15 else Fr(0)
    if not (0 <= c_ < C_): continue
    val = evalP(GT, [Fr(0), p_, j_, c_, e_, Fr(1), b_, a_]); n += 1
    if p_ == b_ or j_ == 1 or e_ == e0 or c_ == 0: walls += 1
    if worst is None or val < worst[0]: worst = (val, [str(x) for x in (p_, j_, c_, e_, b_, a_)])
info('pointwise exact evaluation (seed 11): %d points of the Theorem-1 region (%d on at least one wall p=b, j=1, e=bj/(b-j), c=0); min G(0,p,j,c,e,1,b,a) = %.4e; all positive: %s'
     % (n, walls, float(worst[0]), worst[0] > 0))
# probe just above Crat (outside the theorem): sign of G there, to see that the bound is not vacuous
rnd = random.Random(12); negs = 0; m = 0
for it in range(400):
    b_ = 1 + Fr(rnd.randint(1, 300), 50); r_ = Fr(rnd.randint(1, 999), 1000); a_ = r_**2/b_
    p_ = b_ + Fr(rnd.randint(0, 300), 50); j_ = Fr(rnd.randint(1, 1000), 1000)
    e_ = b_*j_/(b_-j_) + Fr(rnd.randint(0, 300), 50); C_ = crat(b_, r_, p_, e_, a_)
    for fac in (Fr(101, 100), Fr(2), Fr(10)):
        m += 1; negs += evalP(GT, [Fr(0), p_, j_, C_*fac, e_, Fr(1), b_, a_]) < 0
info('outside the region (c = 1.01, 2, 10 times Crat; seed 12): %d of %d points have G < 0' % (negs, m))
finish()
