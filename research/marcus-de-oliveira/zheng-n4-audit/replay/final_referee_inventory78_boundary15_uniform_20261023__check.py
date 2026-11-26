# Own-code replay for the Zheng n=4 audit (batch C40-T77-T79). Executes nothing from the supplement.
import os, sys, itertools
LIB = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib')
sys.path.insert(0, LIB)
import moc4lib as L
import moc4lib as Lb
STAGE = {}
for st, ks in {1:[1,2,3,4,5,6,7,15,18,19,20,21,22,23,25,27,28,29,31,34,35,36,37,39,40,42,43,44,72,88],2:[10,26,46,90],3:[9,11,12,47],
               4:[8,13,17,49],5:[14,16,24,33,38,70],6:[52,57,58,65,68,71,76,78],7:[45,51,66,69,75,81],8:[48,50,53,59,62,67,73,79,86],
               9:[30,41,54,55,60,61,63,74,77,80,83,84],10:[56,64,82,85,87],11:[32],12:[89]}.items():
    for k in ks: STAGE['T%d' % k] = st
for st, ks in {1:[1,2,5,7,8,9,10,11,12,13,14,15,16,17,21,22,24],2:[3,4,23,26],4:[6,18],5:[20,27],7:[19,25]}.items():
    for k in ks: STAGE['A%d' % k] = st
REP_CASE = {L.orbit_rep(T): k for k, T in L.CONTACT_SETS.items()}
REP_AUX = {L.orbit_rep(T): k for k, T in L.AUX7.items()}
COMB_OK = True
def cchk(name, cond):
    global COMB_OK
    print(('PASS ' if cond else 'FAIL ') + name, flush=True); COMB_OK = COMB_OK and bool(cond)
def allminor(T): return all(L.F_IS_MINOR[j] for j in L.generators(T))
def earlier_cover(Tbig, stage_limit):
    """paper cases T_k (6-subsets) and Table-5 sets (7-subsets) contained, up to G, in Tbig with stage < stage_limit"""
    out = []
    for size, table, pre in ((6, REP_CASE, 'T'), (7, REP_AUX, 'A')):
        for sub in itertools.combinations(sorted(Tbig), size):
            k = table.get(L.orbit_rep(sub))
            if k is not None and STAGE['%s%d' % (pre, k)] < stage_limit: out.append(('%s%d' % (pre, k), sub))
    return out
def tie_report(T, stage_of_T, cited=None):
    cited = cited or {}
    for s in range(24):
        if s in T: continue
        T7 = set(T) | {s}
        if allminor(T7):
            print('   +%d: all generators minor' % s); continue
        cov = earlier_cover(T7, stage_of_T)
        names = sorted(set(c[0] for c in cov), key=lambda n: (n[0], int(n[1:])))
        msg = '   +%d: earlier cases contained: %s' % (s, names)
        if s in cited: msg += ' ; document cites %s' % cited[s]
        print(msg)
        cchk('added tie %d is covered by an earlier-stage case%s' % (s, '' if s not in cited else ' and the cited one (%s) is among them' % cited[s]),
             len(cov) > 0 and (s not in cited or cited[s] in names))
print('== combinatorics (A25 = {1,3,4,8,11,12,15}, stage 7; document labels are paper index minus one)')
S7=[1,3,4,8,11,12,15]
cchk('S is Table-5 set A25', frozenset(S7)==L.AUX7[25])
tie_report(S7, 7, {20:'T40',17:'T21',10:'T71',9:'T70',14:'T76',2:'T13'})
import sys, random, sympy as sp
from fractions import Fraction as Fr
ok=True
def chk(name,cond):
    global ok
    print(('PASS ' if cond else 'FAIL ')+name, flush=True); ok=ok and bool(cond)
S=[1,3,4,8,11,12,15]
q,j,c,e,d,b,a,zz=sp.symbols('q j c e d b a z')
coef=[q,j,c,e,d,b,a]; rays=[0,1,9,10,56,53,59]; exc=[20,17,10,9,6,14,2]
chk('generators(S) = {0,1,9,10,56,53,59}', sorted(L.generators(S))==sorted(rays))
chk('only f56 nonminor', [j_ for j_ in rays if not L.F_IS_MINOR[j_]]==[56])
# exclusive gaps
for jj,ee in zip(rays,exc):
    chk('ray f%d has exclusive gap %d'%(jj,ee), L.F_VALUES[jj][ee]==1 and all(L.F_VALUES[k][ee]==0 for k in rays if k!=jj))
ms,Cexpr=L.cubic_sympy()
w0=[sum(cf*sp.Rational(L.F_VECS[jj][i].numerator,L.F_VECS[jj][i].denominator) for cf,jj in zip(coef,rays)) for i in range(14)]
wz=[w0[i]+(zz if i==13 else 0) for i in range(14)]
Fz=sp.Poly(sp.expand(Cexpr.subs(dict(zip(ms,wz)))),zz)
print('degree in z of C(w0+z kappa):',Fz.degree())
C0=sp.expand(Fz.coeff_monomial(1)); Q0=sp.expand(Fz.coeff_monomial(zz))
Qd=a*b+a*d+a*e+a*j+b*c+b*d+b*j+c*d+c*e+c*j+d**2+d*e+d*j+e*j
Cd=(a*b*c+a*b*d+a*b*e+a*c*d+a*c*e+a*c*j+a*d**2+a*d*e+a*d*j+a*e*j+b*c*d+b*c*e+b*c*j+b*c*q+b*d**2+b*d*e+b*d*j+b*d*q+b*e*j
    +c*d**2+c*d*e+c*d*j+c*d*q+d**3+d**2*e+d**2*j+d**2*q+d*e*j)
sgnC=1 if sp.expand(C0-Cd)==0 else (-1 if sp.expand(C0+Cd)==0 else 0)
sgnQ=1 if sp.expand(Q0-Qd)==0 else (-1 if sp.expand(Q0+Qd)==0 else 0)
print('paper cubic vs doc: C sign',sgnC,' Q sign',sgnQ)
chk('C(w0+z kappa) = eps (Cdoc + z Qdoc), eps=+-1, linear in z', Fz.degree()==1 and sgnC!=0 and sgnC==sgnQ)
# Hessian of paper cubic
Hs=sp.hessian(Cexpr,ms)
def Hnum(vals):
    # vals: dict coef->Fraction ; returns 14x14 Fraction Hessian at w0 - C/Q kappa, and Q
    Cv=Fr(str(Cd.subs(vals))); Qv=Fr(str(Qd.subs(vals)))
    wv=[Fr(str(w0[i].subs(vals))) for i in range(14)]; wv[13]-=Cv/Qv
    sub=dict(zip(ms,[sp.Rational(x.numerator,x.denominator) for x in wv]))
    return [[Fr(str(Hs[i,k].subs(sub))) for k in range(14)] for i in range(14)], Qv, wv
def ldl_inertia(M):
    M=[row[:] for row in M]; n=len(M); pos=neg=0; idx=list(range(n)); piv=[]
    # symmetric elimination with diagonal pivoting, 2x2 handled by fallback
    remaining=list(range(n))
    while True:
        p=next((i for i in remaining if M[i][i]!=0),None)
        if p is None: break
        pv=M[p][p]; piv.append(pv)
        if pv>0: pos+=1
        else: neg+=1
        rem=[i for i in remaining if i!=p]
        for i in rem:
            f=M[i][p]/pv
            for k in rem: M[i][k]-=f*M[p][k]
        remaining=rem
    restzero=all(M[i][k]==0 for i in remaining for k in remaining)
    return pos,neg,len(remaining),restzero,piv,[[M[i][k] for k in remaining] for i in remaining]
pt={q:Fr(0),j:Fr(1),c:Fr(1),e:Fr(1),d:Fr(1,2),b:Fr(1),a:Fr(57,70)}
H3,Q3,w3=Hnum(pt)
print('point (3): Q =',Q3,' restored vector (paper coords) =',[str(x) for x in w3])
chk('Q = 101/10 at (3)', Q3==Fr(101,10))
pos,neg,nz,rz,piv,rest=ldl_inertia(H3)
print('elimination at (3): pos',pos,'neg',neg,'remaining',nz,'remainder zero?',rz,' remainder',rest)
chk('inertia at (3) is (7,5,2) or (5,7,2) and rank 12', rz and nz==2 and {pos,neg}=={7,5})
A=[row[2:] for row in H3[2:]]
dA=L.det(A); print('det of principal 12x12 on coords 2..13 at (3):',dA)
dH=L.det(H3); print('det H at (3):',dH)
# random exact tests of det H * Q^2 = -k * G^2 (perfect square)
from math import isqrt
def is_sq(fr):
    if fr<0: return False
    n,dn=fr.numerator,fr.denominator
    return isqrt(n)**2==n and isqrt(dn)**2==dn
rng=random.Random(20261125)
const=None; good=True
for it in range(12):
    vals={v:Fr(rng.randint(1,9),rng.randint(1,5)) for v in coef}
    Hm,Qv,_=Hnum(vals); D=L.det(Hm)
    X_=-D*Qv*Qv
    if const is None:
        # find squarefree part
        import sympy
        sf=sympy.Rational(X_.numerator,X_.denominator)
        print('  sample -detH*Q^2 =',X_)
    for k in [4,1]:
        pass
    good = good and (is_sq(X_/4) or is_sq(X_))
    if not (is_sq(X_/4) or is_sq(X_)): print('  NOT square at',vals, X_)
chk('-detH*Q^2 is a rational square at 12 random positive rational points (exact arithmetic; identity test only)', good)
ok = ok and COMB_OK
print('exact part:', 'ALL EXACT CHECKS PASS' if ok else 'SOME CHECK FAILED')
ok = ok and COMB_OK
import sys, random, sympy as sp
from fractions import Fraction as Fr
rays=[0,1,9,10,56,53,59]
ms,Cexpr=L.cubic_sympy(); Hs=sp.hessian(Cexpr,ms)
Hl=[[sp.lambdify(ms,Hs[i,k]) for k in range(14)] for i in range(14)]
def Cq(q,j,c,e,d,b,a):
    return (a*b*c+a*b*d+a*b*e+a*c*d+a*c*e+a*c*j+a*d**2+a*d*e+a*d*j+a*e*j+b*c*d+b*c*e+b*c*j+b*c*q+b*d**2+b*d*e+b*d*j+b*d*q+b*e*j
    +c*d**2+c*d*e+c*d*j+c*d*q+d**3+d**2*e+d**2*j+d**2*q+d*e*j)
def Qq(q,j,c,e,d,b,a): return a*b+a*d+a*e+a*j+b*c+b*d+b*j+c*d+c*e+c*j+d**2+d*e+d*j+e*j
def detQ2(vals):
    Cv=Cq(*vals); Qv=Qq(*vals)
    w=[sum(cf*L.F_VECS[jj][i] for cf,jj in zip(vals,rays)) for i in range(14)]; w[13]-=Cv/Qv
    M=[[Fr(Hl[i][k](*w)) for k in range(14)] for i in range(14)]
    return -L.det(M)*Qv*Qv/4
qs=sp.Symbol('q')
rng=random.Random(7)
allok=True; nsamp=25
for it in range(nsamp):
    o=[Fr(rng.randint(1,12),rng.randint(1,6)) for _ in range(6)]  # j,c,e,d,b,a
    pts=[(Fr(k),detQ2([Fr(k)]+o)) for k in range(16)]
    poly=sp.interpolate([(sp.Rational(p.numerator,p.denominator),sp.Rational(v.numerator,v.denominator)) for p,v in pts],qs)
    Pq=sp.Poly(sp.expand(poly),qs)
    # check one more point
    chkv=detQ2([Fr(37,3)]+o); assert Pq.eval(sp.Rational(37,3))==sp.Rational(chkv.numerator,chkv.denominator)
    cst,fl=Pq.factor_list()
    G=None
    if cst>0 and all(m%2==0 for _,m in fl) and sp.sqrt(cst).is_rational:
        G=sp.Poly(sp.sqrt(cst),qs)
        for f_,m in fl: G=G*f_**(m//2)
        assert G**2==Pq
    if G is None: print('not a square polynomial in q at',o); allok=False; continue
    co=G.all_coeffs()[::-1]  # co[k] coefficient of q^k
    if co[-1]<0: co=[-x for x in co]
    pos_ok=all(x>0 for x in co[1:])
    j,c,e,d,b,a=[sp.Rational(x.numerator,x.denominator) for x in o]
    s1=a+b+c+e+j; s2=sum(x*y for i,x in enumerate([a,b,c,e,j]) for y in [a,b,c,e,j][i+1:])
    s3=a*b*c+a*b*e+a*c*e+a*c*j+a*e*j+b*c*e+b*c*j+b*e*j
    Pd=d**5+s1*d**4+s2*d**3+s3*d**2-4*a*b*c*e*j
    fac_ok=(abs(co[0])==abs(d*Cq(0,j,c,e,d,b,a)*Pd))
    sign_ok= (co[0]==d*Cq(0,j,c,e,d,b,a)*Pd)
    npos=len([rt for rt in sp.Poly(G.as_expr(),qs).real_roots() if rt>0])
    exp_roots = 1 if Pd<0 else 0
    ok= pos_ok and fac_ok and sign_ok and npos==exp_roots
    allok=allok and ok
    print('sample',it,'degq(G)=',G.degree(),'q^k>=1 coeffs positive:',pos_ok,' G(0)=dC P(d):',sign_ok,' P(d)<0:',Pd<0,' #positive q-roots:',npos, 'OK' if ok else 'MISMATCH', flush=True)
print('SUMMARY univariate-in-q tests (exact arithmetic at %d random rational points of the other six variables):'%nsamp, 'all consistent' if allok else 'MISMATCH FOUND')
# ---- Lemma 3 / 4 on actual sheet points: bordered minors divisible by G, inertia at the positive root
import mpmath
mpmath.mp.dps=60
def Hmat(vals):
    Cv=Cq(*vals); Qv=Qq(*vals)
    w=[sum(cf*L.F_VECS[jj][i] for cf,jj in zip(vals,rays)) for i in range(14)]; w[13]-=Cv/Qv
    return [[Fr(Hl[i][k](*w)) for k in range(14)] for i in range(14)],Qv
rng=random.Random(11); done=0; allok2=True
while done<6:
    o=[Fr(rng.randint(1,12),rng.randint(1,6)) for _ in range(6)]
    j,c,e,d,b,a=o
    s1=a+b+c+e+j; s2=sum(x*y for i,x in enumerate([a,b,c,e,j]) for y in [a,b,c,e,j][i+1:])
    s3=a*b*c+a*b*e+a*c*e+a*c*j+a*e*j+b*c*e+b*c*j+b*e*j
    if d**5+s1*d**4+s2*d**3+s3*d**2-4*a*b*c*e*j>=0: continue
    done+=1
    pts=[];m0=[];m1=[]
    for k in range(16):
        H_,Qv=Hmat([Fr(k)]+o)
        pts.append(-L.det(H_)*Qv*Qv/4)
        m0.append(L.det([row[1:] for row in H_[1:]])*Qv*Qv)
        m1.append(L.det([[H_[i][kk] for kk in range(14) if kk!=1] for i in range(14) if i!=1])*Qv*Qv)
    R_=lambda v: sp.Rational(v.numerator,v.denominator)
    Pq=sp.Poly(sp.interpolate([(k,R_(v)) for k,v in enumerate(pts)],qs),qs)
    cst,fl=Pq.factor_list(); G=sp.Poly(sp.sqrt(cst),qs)
    for f_,m in fl: G=G*f_**(m//2)
    M0=sp.Poly(sp.interpolate([(k,R_(v)) for k,v in enumerate(m0)],qs),qs)
    M1=sp.Poly(sp.interpolate([(k,R_(v)) for k,v in enumerate(m1)],qs),qs)
    div0=M0.rem(G).is_zero; div1=M1.rem(G).is_zero
    root=[rt for rt in G.real_roots() if rt>0][0]
    qn=mpmath.mpf(sp.N(root,70))
    # numeric Hessian at root
    vals=[qn]+[mpmath.mpf(x.numerator)/x.denominator for x in o]
    Cv=Cq(*vals); Qv=Qq(*vals)
    w=[sum(cf*mpmath.mpf(L.F_VECS[jj][i].numerator)/L.F_VECS[jj][i].denominator for cf,jj in zip(vals,rays)) for i in range(14)]; w[13]-=Cv/Qv
    Hn=mpmath.matrix([[Hl[i][kk](*w) for kk in range(14)] for i in range(14)])
    ev=mpmath.eigsy(Hn)[0]
    evs=sorted([float(x) for x in ev])
    npos=sum(1 for x in ev if x>mpmath.mpf(10)**-30); nneg=sum(1 for x in ev if x< -mpmath.mpf(10)**-30)
    ok2=div0 and div1 and (npos,nneg)==(7,5)
    allok2=allok2 and ok2
    print('sheet point',done,'q*=%.6g'%float(qn),' G | Q^2 det H_00:',div0,' G | Q^2 det H_11:',div1,' inertia (num, 60 digits):',(npos,nneg,14-npos-nneg),' smallest |eig|=%.2e'%min(abs(x) for x in evs), 'OK' if ok2 else 'MISMATCH', flush=True)
print('SUMMARY sheet-point tests:', 'all consistent' if allok2 else 'MISMATCH FOUND')
print('EXACT checks:', 'ALL PASS' if ok else 'SOME FAILED'); sys.exit(0 if ok else 1)
