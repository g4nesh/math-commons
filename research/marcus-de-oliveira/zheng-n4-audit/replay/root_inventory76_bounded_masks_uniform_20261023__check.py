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
print('== combinatorics (T77 = {0,1,3,10,11,16}, Appendix G stage 9)')
T=[0,1,3,10,11,16]
cchk('T is paper T_77', frozenset(T)==L.CONTACT_SETS[77])
cchk('generators = {0,2,3,8,10,12,13,15,27,30,32}', L.generators(T)==[0,2,3,8,10,12,13,15,27,30,32])
cchk('only f30, f32 nonminor', [j for j in L.generators(T) if not L.F_IS_MINOR[j]]==[30,32])
V30,V32=L.F_VALUES[30],L.F_VALUES[32]
cchk('f30(2)=f32(2)=1 and f30(6)+f30(9)=f32(6)+f32(9)=1', V30[2]==V32[2]==1 and V30[6]+V30[9]==1 and V32[6]+V32[9]==1)
cchk('all generator values nonnegative (so (1) holds)', all(v>=0 for j in L.generators(T) for v in L.F_VALUES[j]))
cchk('ties 2,4,8,15,17,18,19 leave only minor generators', all(allminor(set(T)|{s}) for s in [2,4,8,15,17,18,19]))
tie_report(T, 9)
import sys, sympy as sp
ok=True
def chk(name,cond):
    global ok
    print(('PASS ' if cond else 'FAIL ')+name); ok=ok and bool(cond)
print('E forward =',Ef)
chk('E forward = -z(y-t)(x-w)', sp.expand(Ef-(-z*(y-t)*(x-w)))==0)
# forward: r,s,u,J
rf=Hf[3]-1; sf=Hf[11]-Hf[10]; uf=Hf[10]; Jf=Hf[2]
chk('raw identity z(J-r)=1', sp.cancel(z*(Jf-rf)-1)==0)
chk('raw identity w r = J X - 1', sp.cancel(w*rf-(Jf*x-1))==0)
sub={r:rf,s:sf,u:uf,J:Jf,X:x}
BNf=BN.subs(sub); ANf=AN.subs(sub); YDf=YD.subs(sub)
chk('raw identity t*BN = AN', sp.cancel(t*BNf-ANf)==0)
chk('raw identity y*YD = s*AN', sp.cancel(y*YDf-sf*ANf)==0)
chk('AN-BN = s(X-1)(JX-r-1)', sp.expand(AN-BN-s*(X-1)*(J*X-r-1))==0)
print('E inverse =',Ei)
chk('E inverse = s(X-1)(J-r)(JX-r-1)((J-r)X-1)AN (up to sign/const)', sp.cancel(Ei/(s*(X-1)*(J-r)*(J*X-r-1)*((J-r)*X-1)*AN)).is_constant())
print('   ratio', sp.cancel(Ei/(s*(X-1)*(J-r)*(J*X-r-1)*((J-r)*X-1)*AN)))
for k,val in [(0,0),(1,1),(3,r+1),(10,u),(11,s+u),(2,J),(5,r+1-J*X)]:
    chk('H%d = %s'%(k,val), sp.cancel(Hi[k]-val)==0)
N16,D16=numden(Hi[16])
C=sp.Poly(sp.expand(N16-h*D16),J,X)
print('C terms',len(sp.Poly(sp.expand(N16-h*D16),J,X,r,s,u,h).terms()),'degJ',C.degree(J),'degX',C.degree(X))
cJ2=sp.factor(sp.Poly(sp.expand(N16-h*D16),J).coeff_monomial(J**2))
print('[J^2]C =',cJ2)
chk('[J^2]C prop (u-1)X^2', sp.cancel(cJ2/((u-1)*X**2)).is_constant())
const=sp.cancel(cJ2/((u-1)*X**2))
Cn=sp.expand((N16-h*D16)/const)
# Newton at J=inf: J=1/tau, X=tau*Y
e=sp.expand(Cn.subs({J:1/tau,X:tau*Y}))
chk('no negative tau powers in weighted C', all(m[0]>=0 for m in sp.Poly(sp.expand(e*tau**4),tau).monoms() for _ in [0]) and min(m[0] for m in sp.Poly(sp.expand(e*tau**4),tau).monoms())==4)
P0=sp.expand(sp.limit(e,tau,0)) if False else sp.expand(e.subs(tau,0))
V=s+u-1; R=r-V
Pd=(u-1)*Y**2+(h*V+r*(2-u)-s-2*u+2)*Y+R*(h-r-1)
print('initial poly',sp.factor(P0))
chk('P(Y) as in (6)', sp.expand(P0-Pd)==0)
# H6,H9 limits along X=tau*Y, J=1/tau
for k,lim in [(6,s+u),(9,(-Y+r*(s+u)+r+1)/(-Y+r+1))]:
    ex=sp.cancel(Hi[k].subs({J:1/tau,X:tau*Y}))
    n,d=sp.fraction(ex)
    n=sp.Poly(sp.expand(n),tau); d=sp.Poly(sp.expand(d),tau)
    # lowest order terms
    ln=min(m[0] for m in n.monoms()); ld=min(m[0] for m in d.monoms())
    nl=n.coeff_monomial(tau**ln); dl=d.coeff_monomial(tau**ld)
    print('H%d orders num/den'%k,ln,ld,' leading ratio',sp.factor(nl/dl))
    chk('H%d limit (7) and equal orders'%k, ln==ld and sp.cancel(nl/dl-lim)==0)
    if k==9:
        res=sp.factor(sp.resultant(Pd,sp.Poly(sp.numer(sp.factor(dl)),Y).as_expr(),Y))
        print('   resultant(P, weighted den of H9) =',res, ' den=',sp.factor(dl))
# generic witness
gen={r:2,s:3,u:4,h:5}
disc=sp.discriminant(Pd,Y).subs(gen)
print('witness (2,3,4,5): u-1=',3,' disc P=',disc,' const term=',(R*(h-r-1)).subs(gen))
chk('witness: P has two distinct nonzero roots', disc!=0 and (R*(h-r-1)).subs(gen)!=0)
# constant strata
Cu1=sp.factor(Cn.subs(u,1))
print('C|u=1 =',Cu1)
chk('(8) factorization', sp.expand(Cn.subs(u,1)-(X-1)*(J*X*(s*(1-h)-r)+(r-s)*(r+1-h))*sp.cancel(Cn.subs(u,1)/((X-1)*(J*X*(s*(1-h)-r)+(r-s)*(r+1-h)))))==0 and sp.cancel(Cn.subs(u,1)/((X-1)*(J*X*(s*(1-h)-r)+(r-s)*(r+1-h)))).is_constant())
# strata: C|u=1 identically 0 iff s(1-h)-r=0 and (r-s)(r+1-h)=0
sol=sp.solve([s*(1-h)-r,(r-s)*(r+1-h)],[r,s,h],dict=True)
print('constant strata solutions',sol)
H6s2=sp.cancel(Hi[6].subs({u:1,s:-1,h:r+1}))
chk('stratum 2 (u,s)=(1,-1): H6 == 0', H6s2==0)
X0=sp.Symbol('X0')
st1={u:1,r:s,X:X0}
H16s1=sp.cancel(Hi[16].subs(st1))
chk('stratum 1: H16 == 0', H16s1==0)
H6s1=sp.cancel(Hi[6].subs(st1)); H9s1=sp.cancel(Hi[9].subs(st1))
chk('(10) H6', sp.cancel(H6s1-(-J*(s+1)/(s-J)))==0)
chk('(10) H9', sp.cancel(H9s1-(J**2*X0-J*(s+1)**2+s*(s+1)**2)/((s-J)*(-J*X0+s+1)))==0)
# gap-majorant bound check: H2=J nonconstant; H5
ok = ok and COMB_OK
print('ALL EXACT CHECKS PASS' if ok else 'SOME CHECK FAILED'); sys.exit(0 if ok else 1)
