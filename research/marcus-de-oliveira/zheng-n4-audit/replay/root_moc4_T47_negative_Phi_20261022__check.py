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
# ===== root_moc4_T47_negative_Phi_20261022 =====
import mpmath as mp
import numpy as np
q,p,j,c,e,d,b,a,h = sp.symbols('q p j c e d b a h'); syms = [q,p,j,c,e,d,b,a]
roles = [2,9,1,45,5,35,0,7]
Z = sorted(set.intersection(*[set(L.F_ZERO[k]) for k in roles]))
chk('source roles f2,f9,f1,f45,f5,f35,f0,f7 are exactly the generators of their common zero set', L.generators(Z) == sorted(roles), Z)
chk('common zero set {0,4,5,6,9,12} is literally paper T_65', frozenset(Z) == L.CONTACT_SETS[65])
chk('only f35 is non-minor', [k for k in roles if not L.F_IS_MINOR[k]] == [35])
w = sym_source(roles, syms, h)
F = sp.Poly(sp.expand(L.cubic_C(w)), h)
chk('F(w(g,h)) is affine in h', F.degree() == 1)
C = sp.expand(F.coeff_monomial(1)); Q = sp.expand(F.coeff_monomial(h))
chk('C = -d{ab+ad+aq+bd+bp+d^2+de+dp+dq+ep+eq+pq}', sp.expand(C + d*(a*b+a*d+a*q+b*d+b*p+d**2+d*e+d*p+d*q+e*p+e*q+p*q)) == 0)
chk('Q = ce-cj-d^2-de-dp-dq-ep-eq+jp-pq', sp.expand(Q - (c*e-c*j-d**2-d*e-d*p-d*q-e*p-e*q+j*p-p*q)) == 0)
vs, Cj, Qj, facs, const = load_radial('root_moc4_T47_full_radial_20261022.json')
sub = dict(zip(vs[:8], syms))
chk('supplement JSON C,Q equal own C,Q', sp.expand(Cj.subs(sub) - C) == 0 and sp.expand(Qj.subs(sub) - Q) == 0)
G = [f for f, m, deg, n in facs if deg == 8][0].subs(sub)
GT = sp.Poly(G, *syms).terms()
chk('G (supplement data) has 965 terms, total degree 8', len(GT) == 965 and sp.Poly(G, *syms).total_degree() == 8)
rnd = random.Random(1); ok = True; n = 0
fC = sp.lambdify(syms, C); fQ = sp.lambdify(syms, Q)
for t in range(25):
    vals = [Fr(rnd.randint(-20, 30), rnd.randint(1, 7)) for _ in range(8)]
    Cv, Qv = Fr(fC(*vals)), Fr(fQ(*vals))
    if Qv == 0 or vals[5] == 0: continue
    dv = L.det(hess_at(source_vec(roles, vals, -Cv / Qv))); n += 1
    ok &= (dv * Qv**2 == -4 * vals[5]**2 * evalP(GT, vals)**2)
chk('det HF(w(g,-C/Q)) * Q^2 == -4 d^2 G^2 at %d random rational points (own Hessian; pointwise exact, identity only with high probability)' % n, ok)
H0 = hess_at(source_vec(roles, [Fr(1)]*8, 0)); H1 = hess_at(source_vec(roles, [Fr(1)]*8, Fr(7, 3)))
keep = [k for k in range(14) if k not in (0, 5)]
chk('principal block deleting coordinates 0,5 is independent of the height', all(H1[i][k] == H0[i][k] for i in keep for k in keep))
Gq = sp.Poly(G, q)
chk('G has degree 2 in q and 2 in e', Gq.degree() == 2 and sp.Poly(G, e).degree() == 2)
A_, B_, C0 = [sp.expand(Gq.coeff_monomial(q**k)) for k in (2, 1, 0)]
chk('q-content gcd(A,B,C0) = 1', sp.gcd(sp.gcd(A_, B_), C0) == 1)
Phi = ((p-c)*j + (b-c)*(p+d+a) - p*a)**2 - 4*a*c*(p-b)*(p+d+a)
qt, rm = sp.div(sp.Poly(sp.expand(B_**2 - 4*A_*C0), p,j,c,e,d,b,a), sp.Poly(sp.expand(e**2*Phi), p,j,c,e,d,b,a))
chk('B^2 - 4 A C0 = e^2 Phi Psi with Psi a polynomial (162 terms)', rm.is_zero and len(qt.terms()) == 162, len(qt.terms()))
s_ = q + d + b
K, r0 = sp.div(sp.expand(G.subs(e, 0)), sp.expand(s_**2), q)
B1, r1 = sp.div(sp.expand(sp.diff(G, e).subs(e, 0)), s_, q)
G2 = sp.expand(sp.diff(G, e, 2).subs(e, 0) / 2)
chk('G = e^2 A2(q) + e(q+d+b)B1(q) + (q+d+b)^2 K with K free of q,e; deg A2 = 2, deg B1 = 1',
    sp.expand(r0) == 0 and sp.expand(r1) == 0 and not K.has(q) and sp.Poly(B1, q).degree() == 1 and sp.Poly(G2, q).degree() == 2
    and sp.expand(G - (e**2*G2 + e*s_*B1 + s_**2*K)) == 0)
A22, A21, A20 = [sp.Poly(G2, q).coeff_monomial(q**k) for k in (2, 1, 0)]
B11, B10 = [sp.Poly(B1, q).coeff_monomial(q**k) for k in (1, 0)]
Lq = p*j + p*b + d*b + b*a
Hm = sp.Matrix([[2*A22, A21, B11], [A21, 2*A20, B10], [B11, B10, 2*K]])
chk('conic: det(2H) = -2 a c d^2 L^2 Phi^2', sp.expand(Hm.det() + 2*a*c*d**2*Lq**2*Phi**2) == 0)
chk('Phi at b = j = 0 is the square [c(p+d+a) - pa]^2', sp.expand(Phi.subs({b: 0, j: 0}) - (c*(p+d+a) - p*a)**2) == 0)
# PROOF-READ helper: Phi<0 with nonnegative coefficients forces a,c>0 and p>b (checked symbolically: Phi = X^2 - 4ac(p-b)(p+d+a))
# ---- NUMERIC 1: sign of the height-free principal minor P on real roots of G = 0 with Phi < 0
mp.mp.dps = 40
fA, fB, fC0 = [sp.lambdify(syms[1:], X_, 'mpmath') for X_ in (A_, B_, C0)]
def Pminor(vals):
    wv = [sum(vals[i] * float(L.F_VECS[roles[i]][k]) for i in range(8)) for k in range(14)]
    return mp.det(mp.matrix([[sum(HESS_LIN[i2][j2][k] * wv[k] for k in range(14)) for j2 in keep] for i2 in keep]))
rnd = random.Random(7); nphi = nroot = npos = 0; worst = None
for it in range(1500):
    pp_, jj, cc, ee, dd, bb, aa = [mp.mpf(rnd.uniform(0.01, 3)) for _ in range(7)]
    if rnd.random() < 0.2: jj = mp.mpf(0)
    if rnd.random() < 0.2: bb = mp.mpf(0)
    X_ = (pp_-cc)*jj + (bb-cc)*(pp_+dd+aa) - pp_*aa
    if X_**2 - 4*aa*cc*(pp_-bb)*(pp_+dd+aa) >= 0: continue
    nphi += 1
    Av, Bv, Cv = fA(pp_,jj,cc,ee,dd,bb,aa), fB(pp_,jj,cc,ee,dd,bb,aa), fC0(pp_,jj,cc,ee,dd,bb,aa)
    disc = Bv**2 - 4*Av*Cv
    if disc < 0 or Av == 0: continue
    for sg in (1, -1):
        qq = (-Bv + sg*mp.sqrt(disc)) / (2*Av); nroot += 1
        Pv = Pminor([qq, pp_, jj, cc, ee, dd, bb, aa])
        if worst is None or Pv > worst[0]: worst = (Pv, [float(x) for x in (qq, pp_, jj, cc, ee, dd, bb, aa)])
        if Pv > 0: npos += 1
info('NUMERIC P-sign test (seed 7, 1500 draws, dps 40): Phi<0 draws %d, real roots q of G=0 tested %d, #P>0 = %d, max P = %s at %s'
     % (nphi, nroot, npos, mp.nstr(worst[0], 6) if worst else None, worst[1] if worst else None))
# ---- NUMERIC 2: cited global facts on random realizable real matching responses w = Re(alpha m)
HLn = np.array(HESS_LIN, float)
rng = np.random.default_rng(3); bad_rank = 0; neg_minor = 0; ns = 200; mins = []
for it in range(ns):
    aa_ = rng.normal(size=4) + 1j*rng.normal(size=4); bb_ = rng.normal(size=4) + 1j*rng.normal(size=4)
    al = np.exp(1j*rng.uniform(0, 2*np.pi))
    wv = L.matching_vector_np(aa_, bb_, al).real
    Hn = HLn @ wv
    sv = np.linalg.svd(Hn, compute_uv=False)
    if sv[12] / sv[0] > 1e-9: bad_rank += 1
    m12 = [np.linalg.det(Hn[np.ix_(Kk, Kk)]) for Kk in itertools.combinations(range(14), 12)]
    mx = max(abs(x) for x in m12); rel = min(m12) / mx; mins.append(rel)
    if rel < -1e-8: neg_minor += 1
info('NUMERIC cited global facts (seed 3, %d random spectra): #rank>12 = %d, #negative order-12 principal minor (rel < -1e-8) = %d, min rel = %.2e'
     % (ns, bad_rank, neg_minor, min(mins)))
finish()
