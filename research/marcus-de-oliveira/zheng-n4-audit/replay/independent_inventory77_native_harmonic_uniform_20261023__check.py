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
# Shared chart for T77 / T78 documents: a=(y,z,w,x), c=(t,0,1,inf), minus brackets.
import itertools, sympy as sp
PERMS=[tuple(p) for p in itertools.permutations(range(4))]
r,s,u,h,J,X,tau,Y=sp.symbols('r s u h J X tau Y')
y,z,w,x,t=sp.symbols('y z w x t')
def br(p,q):  # homogeneous points (p0,p1),(q0,q1): p0*q1-q0*p1 ; finite (a,1)-(c,1)=a-c, (a,1),(1,0) -> -1
    return p[0]*q[1]-q[0]*p[1]
def Fs(A,B):
    return [sp.expand(sp.prod([br(A[i],B[sg[i]]) for i in range(4)])) for sg in PERMS]
# forward chart
Af=[(y,1),(z,1),(w,1),(x,1)]; Bf=[(t,1),(0,1),(1,1),(1,0)]
Ff=Fs(Af,Bf)
Ef=sp.factor(Ff[1]-Ff[0])
Hf=[sp.cancel((F-Ff[0])/(Ff[1]-Ff[0])) for F in Ff]
# inverse chart (3)
BN=(u*(r-J)+r*s+J)*X+u-r*(s+1)-1
AN=s*J*X**2+(u*(r-J)-s*J-s+J)*X+u+s-1-r
YD=(s+J-r)*BN-(J-r)*AN
Ai=[(s*AN,YD),(1,J-r),(J*X-1,r),(X,1)]; Bi=[(AN,BN),(0,1),(1,1),(1,0)]
Fi=Fs(Ai,Bi)
Ei=sp.factor(Fi[1]-Fi[0])
Hi=[sp.cancel((F-Fi[0])/(Fi[1]-Fi[0])) for F in Fi]
def numden(e):
    n,d=sp.fraction(sp.cancel(sp.together(e)))
    return sp.expand(n),sp.expand(d)
print('== combinatorics (paper T_78 = {1,4,5,6,9,12}; working set {0,1,3,10,11,23}; stage 6)')
Rr=(0,3,2,1)
img={L.perm_index(L.compose(Rr,L.inverse(L.PERMS[j]))) for j in [0,1,3,10,11,23]}
cchk('sigma -> R o sigma^{-1} maps working set to {1,4,5,6,9,12} = T_78', img==set(L.CONTACT_SETS[78]))
cchk('literal source = [0,1,2,5,7,9,35,59]', L.generators([1,4,5,6,9,12])==[0,1,2,5,7,9,35,59])
W=[0,1,3,10,11,23]
cchk('working source = [2,3,4,8,12,13,15,32]', L.generators(W)==[2,3,4,8,12,13,15,32])
cchk('zero set of f32 outside W = {5,7,9,16,20,21,22}', sorted(L.F_ZERO[32]-set(W))==[5,7,9,16,20,21,22])
cchk('W+16 is a G-image of Table-5 set A18', L.find_g(set(W)|{16}, L.AUX7[18]) is not None)
cchk('W+22 is a G-image of Table-5 set A20', L.find_g(set(W)|{22}, L.AUX7[20]) is not None)
tie_report(W, 6, {16:'A18', 22:'A20'})
import sys, sympy as sp
ok=True
def chk(name,cond):
    global ok
    print(('PASS ' if cond else 'FAIL ')+name, flush=True); ok=ok and bool(cond)
# ratios
rf=Hf[3]-1; sf=Hf[11]-Hf[10]
chk('r = -(x-z)/(z(x-w))', sp.cancel(rf+(x-z)/(z*(x-w)))==0)
chk('s = y(t-1)/(z(y-t))', sp.cancel(sf-y*(t-1)/(z*(y-t)))==0)
uf=Hf[10]
chk('(4) w = x + (x-z)/(r z)', sp.cancel(w-(x+(x-z)/(rf*z)))==0)
chk('(4) t formula', sp.cancel(t-(uf*(x-w)+sf*w*(x-1)+w-1)/(uf*(x-w)+sf*(x-1)+w-1))==0)
chk('(4) y formula', sp.cancel(y-sf*z*t/(sf*z-t+1))==0)
chk('J = (1+rw)/x equals H2', sp.cancel(Hf[2]-(1+rf*w)/x)==0)
Kk=(J-r)*X-1; Lk=J*X-r-1
chk('(6) F0/E = (r+1-JX)/K', sp.cancel(Fi[0]/(Fi[1]-Fi[0])-(r+1-J*X)/Kk)==0)
# decomposition on chart: H_sigma = sum over working source
T=[0,1,3,10,11,23]; gens=[2,3,4,8,12,13,15,32]; exc=[21,20,16,9,22,7,5,2]
chk('generators(working T) as listed', sorted(Lb.generators(T))==gens)
def kernel_decomp(T,gens,exc):
    # every real c with (Pc)_T=0 satisfies Pc = sum_j (Pc)_{exc_j} * values(f_j)
    rows=[Lb.P[t_] for t_ in T]
    ns=Lb.nullspace(rows)
    good=len(ns)==14-len(T)
    for cvec in ns:
        vv=Lb.values(cvec)
        good=good and all(vv[sg]==sum(vv[e]*Lb.F_VALUES[j][sg] for j,e in zip(gens,exc)) for sg in range(24))
    return good,len(ns)
kd,dim=kernel_decomp(T,gens,exc)
print('contact kernel dim',dim)
chk('(1) contact kernel (dim 8): Pc = sum_j (Pc)_exc(j) f_j for a basis, hence for all real gap vectors vanishing on T', kd and dim==8)
# denominators of normalized coordinates divide (J-r)(X-1) L AN
dens=set()
for sg in range(24):
    dn=sp.fraction(sp.cancel(Hi[sg]))[1]
    for f_,m in sp.factor_list(dn)[1]: dens.add(sp.factor(f_))
print('distinct denominator factors of H_sigma:',dens)
chk('every denominator factor is one of J-r, X-1, L, AN', all(any(sp.cancel(f_/g).is_constant() for g in [J-r,X-1,Lk,AN]) for f_ in dens))
N23,D23=numden(Hi[23])
Cc=sp.expand(N23-h*D23)
P=sp.Poly(Cc,J,X)
print('C terms',len(sp.Poly(Cc,J,X,r,s,u,h).terms()),'bideg',P.degree(J),P.degree(X), flush=True)
V=s+u-1
cJ4=sp.factor(sp.Poly(Cc,J).coeff_monomial(J**4)); print('[J^4]C raw',cJ4)
cn=sp.cancel(cJ4/(V*(X-1)*X**2)); print('normaliser',cn)
chk('[J^4]C prop V(X-1)X^2', cn.is_constant())
Cn=sp.expand(Cc/cn)
FX=V*J**2-(r*(2*V+1)+s*h+u-1)*J+r*((r+1)*(V+1)+s*h)
chk('[X^3]C = J^2 FX', sp.expand(sp.Poly(Cn,X).coeff_monomial(X**3)-J**2*FX)==0)
# (10) initial polynomial at J=inf, X=tau*Y
e=sp.expand(Cn.subs({J:1/tau,X:tau*Y})*tau**2)
pe=sp.Poly(sp.expand(e*tau**6),tau); lo=min(m[0] for m in pe.monoms())
init=sp.factor(pe.coeff_monomial(tau**lo)); print('lowest tau order (after *tau^2):',lo-6,' initial',init)
chk('(10) initial = (r+1-Y)(VY+r-V), at order 0', lo==6 and sp.cancel(pe.coeff_monomial(tau**lo)/((r+1-Y)*(V*Y+r-V))).is_constant())
# (7) sections (symbolic)
print('    C(J,1) =',sp.factor(Cn.subs(X,1)))
chk('(7) C(J,1) = +(u-1)(J-r-1)^2[(s+u+h-2)J-rh+r-s-u+1]  (document prints a leading minus sign; same zero set; NOTE presentation)', sp.expand(Cn.subs(X,1)-(u-1)*(J-r-1)**2*((s+u+h-2)*J-r*h+r-s-u+1))==0)
RL=sp.factor(sp.expand(Cn.subs(X,(r+1)/J)*J**3)); print('J^3 C(J,(r+1)/J) [= Res_X(C,L) up to sign/power of J] =',RL, flush=True)
chk('Res_X(C,L) prop r^2(s+u)^2 J^2 (J-r-1)^2 (up to power of J)', sp.cancel(RL/(r**2*(s+u)**2*(J-r-1)**2)).as_numer_denom()[1].is_constant() and all(f_ in (J,) or sp.Poly(f_,J).degree()==0 for f_,_ in sp.factor_list(sp.cancel(RL/(r**2*(s+u)**2*(J-r-1)**2)))[1]))
RJr=sp.factor(Cn.subs(J,r)); print('C(r,X) =',RJr)
cV=sp.factor(sp.Poly(sp.expand(Cn*cn).subs(u,1-s),J,X).coeff_monomial(J**3*X)) ; print('V=0: [J^3 X](N23-hD23) =',cV)
chk('V=0 wall: [J^3 X] coefficient = -r (independent of h)', sp.expand(cV+r)==0)
# (8)
lim8=sp.limit(sp.cancel((Hi[23]-Hi[22]-J).subs({r:sp.Rational(2),s:3,u:4})),X,sp.oo); print('(8) limit at (r,s,u)=(2,3,4):',lim8)
chk('(8) H23-H22-J -> -r-1 as X->inf (at r,s,u=2,3,4)', sp.cancel(lim8+3)==0)
# branches at two exact generic parameter points
for gp in [{r:2,s:3,u:4,h:5},{r:sp.Rational(7,3),s:-2,u:sp.Rational(5,2),h:sp.Rational(-3,4)}]:
    Cg=sp.expand(Cn.subs(gp)); Vg=V.subs(gp); rg=gp[r]; sg_=gp[s]; ug=gp[u]; hg=gp[h]
    Hg={k_:sp.cancel(Hi[k_].subs(gp)) for k_ in (5,7,16)}
    y2=sp.Symbol('y2')
    Y1a=rg*(sg_+ug); Y1b=-rg*(sg_+ug)*(rg-Vg)/Vg**2; c1v=-(ug-1)*(hg+sg_+ug-2)/Vg
    brs=[(7,tau*((rg+1)+Y1a*tau+y2*tau**2),4),(16,tau*((Vg-rg)/Vg+Y1b*tau+y2*tau**2),4),(5,1+c1v*tau+y2*tau**2,4)]
    for k_,b,mult in brs:
        ser=sp.Poly(sp.expand(Cg.subs({J:1/tau,X:b})*tau**mult),tau)
        mo=sorted(set(m_[0] for m_ in ser.monoms()))
        low=[sp.expand(ser.coeff_monomial(tau**kk)) for kk in range(0,mo[0]+2)]
        print('   branch for H%d: lowest tau orders of C along trial series: coefficients'%k_,[c_ for c_ in low][:4])
        # first two orders must vanish identically (independent of y2): branch + first correction are right
        o=[sp.expand(ser.coeff_monomial(tau**kk)) for kk in range(0,6)]
        nz=[kk for kk in range(6) if o[kk]!=0]
        firstnz=nz[0] if nz else None
        chk('(11)/(12) trial series solves C through the first-correction order (H%d branch); first nonzero order %s involves y2: %s'%(k_,firstnz, (o[firstnz].has(y2) if firstnz is not None else None)), firstnz is not None and o[firstnz].has(y2))
        ex=sp.cancel((Hg[k_]+J).subs({J:1/tau,X:b}))
        n_,d_=sp.fraction(ex); n_=sp.Poly(sp.expand(n_),tau); d_=sp.Poly(sp.expand(d_),tau)
        ln=min(m_[0] for m_ in n_.monoms()); ld=min(m_[0] for m_ in d_.monoms())
        lim=sp.cancel(n_.coeff_monomial(tau**ln)/d_.coeff_monomial(tau**ld)) if ln==ld else (0 if ln>ld else 'pole')
        print('   H%d+J on its branch: num order %d, den order %d, limit %s'%(k_,ln,ld,lim))
        chk('(13) H%d + J bounded on its branch: orders equal and limit polynomial in the (real) next Taylor coefficient y2 (params)'%k_, ln>=ld and (lim==0 or sp.Poly(lim,y2).degree()>=0))
# ---- symbolic (11)-(13): trial series X = tau*(Y0+Y1*tau+y2*tau^2) (resp. 1+c1*tau+y2*tau^2), J = 1/tau
y2=sp.Symbol('y2')
brs=[(7,tau*((r+1)+r*(s+u)*tau+y2*tau**2),-2),(16,tau*((V-r)/V-r*(s+u)*(r-V)/V**2*tau+y2*tau**2),-2),(5,1-(u-1)*(h+s+u-2)/V*tau+y2*tau**2,-3)]
for k_,b,must in brs:
    n_,d_=sp.fraction(sp.together(Cn.subs({J:1/tau,X:b}))); pn=sp.Poly(sp.expand(n_),tau); pd=sp.Poly(sp.expand(d_),tau)
    lo=min(m[0] for m in pn.monoms())-min(m[0] for m in pd.monoms()); lead=pn.coeff_monomial(tau**min(m[0] for m in pn.monoms()))
    chk('symbolic branch for H%d: C along trial series is O(tau^%d) with leading coefficient involving y2 (all lower orders, down to tau^%d, vanish identically)'%(k_,lo,must), lo>must and lead.has(y2))
    n_,d_=sp.fraction(sp.together((Hi[k_]+J).subs({J:1/tau,X:b}))); pn=sp.Poly(sp.expand(n_),tau); pd=sp.Poly(sp.expand(d_),tau)
    ln=min(m[0] for m in pn.monoms()); ld=min(m[0] for m in pd.monoms()); dl=sp.factor(pd.coeff_monomial(tau**ld))
    print('   H%d+J: num order %d, den order %d, den leading coefficient %s'%(k_,ln,ld,dl))
    chk('symbolic (13): H%d + J bounded on its branch (den leading coefficient free of y2, generically nonzero)'%k_, ln>=ld and not dl.has(y2))
ok = ok and COMB_OK
print('ALL EXACT CHECKS PASS' if ok else 'SOME CHECK FAILED'); sys.exit(0 if ok else 1)
