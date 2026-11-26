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
# ===== root_true_inventory47_native_uniform_20261023 =====
x,z,w,v,t,r,s,u,h,Y,eps=sp.symbols('x z w v t r s u h Y epsilon')
A=[0,1,z,w]; Cc=[None,v,x,t]
def Fs(sig, inf_factor):
    pr=1
    for i in range(4):
        pr*= inf_factor if sig[i]==0 else (A[i]-Cc[sig[i]])
    return pr
for inf in (-1,1):
    F=[sp.expand(Fs(L.PERMS[k],inf)) for k in range(24)]
    E=sp.factor(F[12]-F[2]); info('point at infinity factor',inf,' F12-F2 =',E, ' F12/E =', sp.cancel(F[12]/E))
F=[sp.expand(Fs(L.PERMS[k],-1)) for k in range(24)]
E=(z-v)*(w-t)
chk('E = F12-F2 = (z-v)(w-t)', sp.expand(F[12]-F[2]-E)==0)
chk('F12/E = x', sp.cancel(F[12]/E-x)==0)
chk('F13-F12 = x(t-v)(z-w)', sp.expand(F[13]-F[12]-x*(t-v)*(z-w))==0)
chk('F0-F2 = (v-x)(z-1)(w-t)', sp.expand(F[0]-F[2]-(v-x)*(z-1)*(w-t))==0)
H=[sp.cancel((F[k]-F[12])/E) for k in range(24)]
chk('H2 = -1, H12 = 0', H[2]==-1 and H[12]==0)
V=r+x*(1+u-r); TD=z*(x+r)+x*(s-r); TN=z*(x**2+x*(s+1)+r)-x*(x+r)
Nv=z*(x+h+1)-x
sub={t:TN/TD, v:Nv/(z+h)}
sub[w]=V*sub[v]/(x*(x+u))
Wd=(z+h)*x*(x+u); Wn=V*Nv; Vd=(z+h); Tn=TN; Td=TD
def clear(expr):
    """numerator after substituting w=Wn/Wd, v=Nv/Vd, t=Tn/Td into a polynomial multilinear in w,v,t, times Wd*Vd*Td"""
    P=sp.Poly(sp.expand(expr),w,v,t); out=0
    for (i,jv,k),cf in P.terms():
        assert i<=1 and jv<=1 and k<=1
        out+=cf*(Wn if i else Wd)*(Nv if jv else Vd)*(Tn if k else Td)
    return sp.expand(out)
def Hnum(k, target):  # numerator of H_k - target
    return clear(F[k]-F[12]-target*E)
def modC(n, name):
    n=sp.expand(n)
    if n==0: chk(name+' (identically)',True); return
    pr=sp.prem(sp.Poly(n,z), sp.Poly(CURVE[0],z))
    chk(name+' (mod C: pseudo-remainder in z is zero)', pr.is_zero)
def numer(e_):
    nn,dd=sp.fraction(sp.together(e_)); return nn
for k,tg,nm in ((5,s,'s'),(9,u,'u'),(0,h,'h')):
    chk('inverse (3) gives H%d = %s identically'%(k,nm), Hnum(k,tg)==0) if k==0 else info('H%d = %s tested modulo C below'%(k,nm))
print('progress',flush=True)
chk('TD - TN = x(x+s)(1-z)', sp.expand(TD-TN-x*(x+s)*(1-z))==0)
chk('Nv - (x+h+1)(z+h) = -(h+1)(x+h)', sp.expand(Nv-(x+h+1)*(z+h)+(h+1)*(x+h))==0)
num=Hnum(13,r)
fl=sp.factor_list(num)
print('H13 - r cleared numerator factors:', [(f if sp.Poly(f,x,z).total_degree()<=1 else 'deg%d'%sp.Poly(f,x,z).total_degree(), m) for f,m in fl[1]], fl[0], flush=True)
Ccurve=[f for f,m in fl[1] if f.has(x) and f.has(z) and sp.Poly(f,x,z).total_degree()>1]
assert len(Ccurve)==1
Cr=sp.expand(Ccurve[0])
c0=sp.expand(Cr.subs(x,0)); ratio=sp.cancel(c0/(r**2*(h+1)*z**2))
print('C(0,z)/(r^2(h+1)z^2) =', ratio)
Cr=sp.expand(Cr/ratio)
CURVE=[Cr]
PC=sp.Poly(Cr,x,z)
print('curve C: bidegree', (sp.degree(Cr,x), sp.degree(Cr,z)), ' terms in (x,z,r,s,u,h):', len(sp.Poly(Cr,x,z,r,s,u,h).terms()), flush=True)
chk('(x+r) divides cleared H13-r numerator', any(sp.expand(f-(x+r))==0 for f,m in fl[1]))
chk('C(0,z) = r^2(h+1) z^2', sp.expand(Cr.subs(x,0)-r**2*(h+1)*z**2)==0)
chk('C(-r,z) = -r^2(r-u)(r-s)(z-1)(z+r)', sp.expand(Cr.subs(x,-r)+r**2*(r-u)*(r-s)*(z-1)*(z+r))==0)
Aa=(u+1)*(r*(h+1)+h-s); Bb=(h-s)*(h+1)*(r+1)
chk('C(1,z) = A(z-1)^2', sp.expand(Cr.subs(x,1)-Aa*(z-1)**2)==0)
chk('Cx(1,1) = B', sp.expand(sp.diff(Cr,x).subs({x:1,z:1})-Bb)==0)
chk('Cz(1,1) = 0 and Czz(1,1)/2 = A', sp.expand(sp.diff(Cr,z).subs({x:1,z:1}))==0 and sp.expand(sp.diff(Cr,z,2).subs({x:1,z:1})/2-Aa)==0)
mindeg=min(sum(m) for m in PC.monoms()); chk('every monomial of C has (x,z)-degree >= 2', mindeg==2, mindeg)
init=sp.expand(sp.expand(Cr.subs({x:eps,z:Y*eps})).coeff(eps,2))
Pdoc=r*((h+1)*r*Y**2+(-h**2-h*r+h*s-h*u-h-2*r+s)*Y+h*(u+1)+r-s)
chk('initial form P(Y) as displayed', sp.expand(init-Pdoc)==0)
chk('P(1) = h r (s-h)', sp.expand(Pdoc.subs(Y,1)-h*r*(s-h))==0)
chk('z-leading coefficient of C vanishes at none of x = 0, 1, -r generically', all(sp.expand(sp.Poly(Cr,z).LC().subs(x,x0))!=0 for x0 in (0,1,-r)), [sp.factor(sp.Poly(Cr,z).LC().subs(x,x0)) for x0 in (0,1,-r)])
info('z-leading coefficient of C:', sp.factor(sp.Poly(Cr,z).LC()), ' x-leading coefficient:', sp.factor(sp.Poly(Cr,x).LC()))
for k5,tg,nm in ((5,s,'H5 = s'),(9,u,'H9 = u')):
    nn5=Hnum(k5,tg)
    modC(nn5, 'inverse (3) gives %s modulo C'%nm)
    info('%s cleared numerator factors:'%nm, [(sp.factor(f) if sp.Poly(f,x,z).total_degree()<=1 else 'deg%d'%sp.Poly(f,x,z).total_degree(), m) for f,m in sp.factor_list(nn5)[1]])
# H_k after substitution = Hnum-style: H_k - T = clear(F_k - F_12 - T*E)/(E*Wd*Vd*Td) with T rational
def Hident(k, T, name):
    Tn_,Td_=sp.fraction(sp.together(T))
    n=clear(sp.expand(Td_*(F[k]-F[12])-Tn_*E))
    modC(n,name)
Hident(11, x*(x+s)*(x+u)/((x-1)*(x+r))-x, 'H11 identity (2)')
n14=sp.expand((x-z)*(F[14]-F[12])-x*(h+z)*E)
# before imposing C means: holds after substituting the inverse (3) only
chk('H14 = x(h+z)/(x-z) after inverse (3), without C', clear(n14)==0)
Hident(19, (z*(x*(s+1)+r)-x*r)/(x*(z-1)), 'H19 identity (4)')
vv=Nv/(z+h); ww=V*Nv/((z+h)*x*(x+u))
Hident(1, sp.together(-(1-vv)*(ww-x)*(x+u)/(vv*(1-x))-x), 'H1 identity (5)')
chk('H19 rewrite', sp.cancel((z*(x*(s+1)+r)-x*r)/(x*(z-1)) - ((s+1)*z/(z-1)+r*z/(x*(z-1))-r/(z-1)))==0)
chk('x F0 = (x-z) F14 (raw)', sp.expand(x*F[0]-(x-z)*F[14])==0)
print('v,w at (x,z)=(-r,1):', sp.simplify(vv.subs({x:-r,z:1})), sp.simplify(ww.subs({x:-r,z:1})))
print('v,w at (x,z)=(-r,-r):', sp.simplify(vv.subs({x:-r,z:-r})), sp.simplify(ww.subs({x:-r,z:-r})))
print('v,w at (1,1):', sp.simplify(vv.subs({x:1,z:1})), sp.simplify(ww.subs({x:1,z:1})))
# limits for infinity analysis
print('H11 at x->oo:', sp.limit(x*(x+s)*(x+u)/((x-1)*(x+r))-x, x, sp.oo), ' H14 z->oo:', sp.limit(x*(h+z)/(x-z),z,sp.oo), ' H14 x->oo:', sp.limit(x*(h+z)/(x-z),x,sp.oo))
# generic parameter check r,s,u,h = 2,3,5,7
g={r:2,s:3,u:5,h:7}
print('at r,s,u,h=2,3,5,7: A,B =', Aa.subs(g), Bb.subs(g), ' disc P(Y) =', sp.discriminant(Pdoc.subs(g),Y), ' r(r+1)(r-u)(h-r)(h+1) =', (r*(r+1)*(r-u)*(h-r)*(h+1)).subs(g))
# gap bound inequality (1): profile values
FV=L.F_VALUES
print('profile values f40 at 1,11,14,19:', [FV[40][k] for k in (1,11,14,19)], ' f43:', [FV[43][k] for k in (1,11,14,19)])
gens=L.generators(L.CONTACT_SETS[48])
print('other generators at 1,11,14,19 (all >=0 by facet property):', {gj:[FV[gj][k] for k in (1,11,14,19)] for gj in gens})
T = L.CONTACT_SETS[48]
chk('document contact set {0,2,5,9,12,13} is literally paper T_48', T == frozenset({0,2,5,9,12,13}))
chk('generators(T_48) = {0,1,5,7,12,39,40,41,43,45}, non-minor exactly {40,43}', L.generators(T) == [0,1,5,7,12,39,40,41,43,45] and [k for k in L.generators(T) if not L.F_IS_MINOR[k]] == [40,43])
chk('inequality (1): f40 and f43 are >= 1 at 1, at 11 or 14, and at 11 or 19 (all other generators >= 0)',
    all(FV[k][1] >= 1 and FV[k][11] + FV[k][14] >= 1 and FV[k][11] + FV[k][19] >= 1 for k in (40, 43)))
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
kill = [jj for jj in range(24) if jj not in T and FV[40][jj] != 0 and FV[43][jj] != 0]
chk('ties forcing both non-minor coefficients to vanish are exactly {1,4,8,10,20,22}', kill == [1,4,8,10,20,22], kill)
for jj in [3,6,7,11,14,15,16,17,18,19,21,23]:
    subs = sorted(set(classify6(S6) for S6 in itertools.combinations(sorted(T | {jj}), 6) if frozenset(S6) != T))
    early = [sx for sx in subs if sx[1] is not None and sx[1] < 8]
    chk('tie %d: T_48+{%d} contains a six-set of an earlier stage (Appendix G stage < 8)' % (jj, jj), bool(early), early)
finish()
