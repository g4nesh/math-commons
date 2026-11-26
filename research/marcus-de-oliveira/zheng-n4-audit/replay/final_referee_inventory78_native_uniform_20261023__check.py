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
print('== combinatorics (T79 = {1,3,4,8,11,12}, stage 8; document labels are paper index minus one)')
T=[1,3,4,8,11,12]
cchk('T is paper T_79', frozenset(T)==L.CONTACT_SETS[79])
cchk('zero set of f56 outside T = {2,9,10,14,15,17,20}', sorted(L.F_ZERO[56]-set(T))==[2,9,10,14,15,17,20])
cchk('T+15 equals Table-5 set A25 literally', set(T)|{15}==set(L.AUX7[25]))
tie_report(T, 8, {2:'T13',9:'T70',10:'T71',14:'T78',17:'T58',20:'T40',15:'A25'})
import sys, itertools, sympy as sp
PERMS=[tuple(p) for p in itertools.permutations(range(4))]
r,s,t,u,x,z,K,tau,Y,al,be=sp.symbols('r s t u x z K tau Y alpha beta')
w,v,k=sp.symbols('w v k')
ok=True
def chk(name,cond):
    global ok
    print(('PASS ' if cond else 'FAIL ')+name, flush=True); ok=ok and bool(cond)
def br(p,q): return p[0]*q[1]-q[0]*p[1]
def Fs(A,B): return [sp.expand(sp.prod([br(A[i],B[sg[i]]) for i in range(4)])) for sg in PERMS]
# forward
Ff=Fs([(w,1),(0,1),(1,1),(z,1)],[(v,1),(1,0),(k*z,1),(x,1)])
Ef=sp.factor(Ff[3]-Ff[1]); print('E fwd',Ef)
chk('E = (w-v)z(1-x)', sp.expand(Ff[3]-Ff[1]-(w-v)*z*(1-x))==0)
Hf=[sp.cancel((F-Ff[1])/(Ff[3]-Ff[1])) for F in Ff]
print('H4 fwd',sp.factor(Hf[4]))
chk('H4 = (k-1)/(x-1), i.e. k = 1+s(x-1)', sp.cancel(Hf[4]-(k-1)/(x-1))==0)
kk=1+s*(x-1)
A0=x-z; B0=kk*z-1; R=r+kk-1; T0=t+kk-1
Lq=x*B0*R-z*kk*A0*T0; n=z*(x*B0*R-kk*A0*T0); Yq=-x*kk*A0*B0*(z-1)
# the two linear equations: check they are equivalent to H8=r, H11=t at k=kk, y=w-v
yv=sp.Symbol('y')
Hs=[sp.together(h.subs({k:kk,w:v+yv})) for h in Hf]
e8=sp.numer(sp.together(Hs[8]-r)); e11=sp.numer(sp.together(Hs[11]-t))
l1=A0*kk*(v-1)+(x-1)*R*yv; l2=x*B0*(v-z)+z*(x-1)*T0*yv
print('H8-r numerator factor:',sp.factor(e8)); print('H11-t numerator factor:',sp.factor(e11))
chk('H8=r <=> linear eq 1 (up to nonzero factor)', sp.cancel(sp.expand(e8)/sp.expand(l1)).is_rational_function() and sp.Poly(sp.cancel(sp.expand(e8)/sp.expand(l1)).as_numer_denom()[0],v,yv).degree()<=0 and sp.rem(sp.expand(e8),sp.expand(l1),v)==0)
chk('H11=t <=> linear eq 2 (up to nonzero factor)', sp.rem(sp.expand(e11),sp.expand(l2),v)==0)
vs=n/Lq; ys=Yq/((x-1)*Lq)
chk('Cramer solution satisfies eq1', sp.cancel(l1.subs({v:vs,yv:ys}))==0)
chk('Cramer solution satisfies eq2', sp.cancel(l2.subs({v:vs,yv:ys}))==0)
# inverse homogeneous
A=[((x-1)*n+Yq,(x-1)*Lq),(0,1),(1,1),(z,1)]; Cc=[(n,Lq),(1,0),(kk*z,1),(x,1)]
Fi=Fs(A,Cc); Ei=Fi[3]-Fi[1]
H=[sp.cancel((F-Fi[1])/Ei) for F in Fi]
chk('F1/E = s(x-1)', sp.cancel(Fi[1]/Ei-s*(x-1))==0)
for kx,val in [(3,1),(4,s),(8,r),(11,t)]:
    chk('H%d = %s'%(kx,val), sp.cancel(H[kx]-val)==0)
chk('(4) H2', sp.cancel(H[2]-(1-kk*(z-1)/(x-1)))==0)
chk('(5) H5', sp.cancel(H[5]-(x*(kk*z-1)/(z*(x-1))-(kk-1)))==0)
chk('(6) H6', sp.cancel(H[6]-((x*(kk*z-1)*R-kk*(x-z)*T0)/(x*kk*(z-1))-(kk-1)))==0)
chk('(7) H10', sp.cancel(H[10]-(s*(x-1)*(x*R/(kk*(x-z))-1)))==0)
# generator decomposition: gap g = sum coef f_j, coefficients = gaps at exclusive positions
T=[1,3,4,8,11,12]; gens=[0,1,9,10,56,53,59,2]; exc=[20,17,10,9,6,14,2,15]
chk('generators(T79) as listed', sorted(L.generators(T))==sorted(gens))
# complete decomposition holds for all 24 positions: g_sigma = sum_j g_{exc_j} f_j(sigma) for gaps vanishing on T
def kernel_decomp(T,gens,exc):
    # every real c with (Pc)_T=0 satisfies Pc = sum_j (Pc)_{exc_j} * values(f_j)
    rows=[L.P[t_] for t_ in T]
    ns=L.nullspace(rows)
    good=len(ns)==14-len(T)
    for cvec in ns:
        vv=L.values(cvec)
        good=good and all(vv[sg]==sum(vv[e]*L.F_VALUES[j][sg] for j,e in zip(gens,exc)) for sg in range(24))
    return good,len(ns)
kd,dim=kernel_decomp(T,gens,exc)
print('contact kernel dim',dim)
chk('(1) contact kernel (dim 8): Pc = sum_j (Pc)_exc(j) f_j for a basis, hence for all real gap vectors vanishing on T', kd and dim==8)
print('f_j(5):',{j:L.F_VALUES[j][5] for j in gens})
# curve
N12,D12=sp.fraction(sp.cancel(H[12]))
Cpoly=sp.expand(N12-u*D12)
P=sp.Poly(Cpoly,x,z)
print('C terms (in x,z,r,s,t,u):',len(sp.Poly(Cpoly,x,z,r,s,t,u).terms()),' bideg',P.degree(x),P.degree(z), flush=True)
# normalise by constant so that [x^5 z^3] = s^3(r+s-u-1)
c53=sp.factor(P.coeff_monomial(x**5*z**3)); print('[x5z3]C raw',c53)
cn=sp.cancel(c53/(s**3*(r+s-u-1))); print('normalising const',cn)
Cn=sp.expand(Cpoly/cn)
Bx=s*x*(r+kk-1)-kk*(t+kk-1)
chk('(8) C(x,1) = (x-1)^3 B(x)^2', sp.expand(Cn.subs(z,1)-(x-1)**3*Bx**2)==0)
chk('(9) C(1,z) = -(r+t)(z-1)^3', sp.expand(Cn.subs(x,1)+(r+t)*(z-1)**3)==0)
chk('(11) [z^3]C at x=0 = -(s-1)^2(s-t)^2', sp.expand(sp.Poly(Cn.subs(x,0),z).coeff_monomial(z**3)+(s-1)**2*(s-t)**2)==0)
print('   C(0,z) =',sp.factor(Cn.subs(x,0)))
chk('(12) [x^5]C at z=1 = s^2(r+s-t-1)^2', sp.expand(sp.Poly(Cn.subs(z,1),x).coeff_monomial(x**5)-s**2*(r+s-t-1)**2)==0)
print('   B linear coefficient', sp.factor(sp.Poly(sp.expand(Bx),x).coeff_monomial(x)))
# Newton at x=z=0, weights x:2, z:3 -> edge 3degx+2degz=6
Pc=sp.Poly(Cn,x,z)
wts=[3*m[0]+2*m[1] for m in Pc.monoms()]
print('min weighted degree',min(wts))
init=sum(c*x**m[0]*z**m[1] for m,c in zip(Pc.monoms(),Pc.coeffs()) if 3*m[0]+2*m[1]==6)
chk('initial form at (0,0)', sp.expand(init-(-(s-1)**2*(s-t)**2*z**3+(s-1)*(s-u)*(r-s)*x**2))==0 and min(wts)==6)
# H5 -> s along x ~ z^{3/2}: param x = c*T^3, z = T^2 -> x/z -> 0; check H5 limit
Tt,cc=sp.symbols('T c')
h5=sp.cancel(H[5].subs({x:cc*Tt**3,z:Tt**2}))
print('H5 along (x,z)=(cT^3,T^2) limit:',sp.limit(h5,Tt,0))
# k=0 section
xs=(s-1+K)/s
CK=sp.expand(s**5*Cn.subs(x,xs))
sec=sp.factor(CK.subs(K,0)); print('K=0 section',sec)
chk('K=0 section = -s^2(r-1)^2(s-1)^2 z', sp.expand(CK.subs(K,0)+s**2*(r-1)**2*(s-1)**2*z)==0)
# residues of H6,H10 at K=0 (k = K there since kk = 1+s(x-1) = K at x=xs)
chk('k = K at x=(s-1+K)/s', sp.expand(kk.subs(x,xs)-K)==0)
S610=sp.factor(sp.cancel((H[6]+H[10]).subs(x,xs)))
print('H6+H10 on k=K chart:',S610)
nS,dS=sp.fraction(sp.cancel((H[6]+H[10]).subs(x,xs)))
print('  den factors', sp.factor_list(dS))
print('  num at K=0,z=0:', sp.expand(nS.subs({K:0,z:0})))
# implicit graph z(K): z = c1 K + ...; check H6+H10 holomorphic: compute series
# solve CK=0 for z as series in K to order 2
c1=sp.Symbol('c1'); c2=sp.Symbol('c2')
zser=c1*K+c2*K**2
ser=sp.series(CK.subs(z,zser),K,0,3).removeO()
sol1=sp.solve(sp.Poly(ser,K).coeff_monomial(K),c1)
print('  z ~ c1 K, c1 =',sol1)
lim=sp.limit(sp.cancel((H[6]+H[10]).subs(x,xs).subs(z,sol1[0]*K)),K,0) if sol1 else None
print('  lim (H6+H10) along graph (first order) =',sp.simplify(lim))
# need second order to be safe: residue depends only on first order of z? compute with c2 too
if sol1:
    sol2=sp.solve(sp.Poly(ser.subs(c1,sol1[0]),K).coeff_monomial(K**2),c2)
    lim2=sp.limit(sp.cancel((H[6]+H[10]).subs(x,xs).subs(z,sol1[0]*K+sol2[0]*K**2)),K,0)
    print('  lim with 2nd order =',sp.simplify(lim2))
    chk('H6+H10 bounded with real limit on k=0 graph (first vs second order agree)', sp.simplify(lim-lim2)==0 and lim.is_finite is not False)
# z=1 slope cubic at (1,1)
Ca=sp.expand(Cn.subs({x:1+al,z:1+be}))
Pa=sp.Poly(Ca,al,be)
md=min(sum(m) for m in Pa.monoms()); print('least total degree at (1,1):',md)
lowf=sum(c*al**m[0]*be**m[1] for m,c in zip(Pa.monoms(),Pa.coeffs()) if sum(m)==md)
slope=sp.factor(sp.expand(lowf.subs(be,Y*al)/al**md))
print('slope cubic',slope)
chk('slope cubic (Yr+Yt+rs-t)(-Y^2+Y(t-u+1)+rs-t) up to const', sp.cancel(sp.expand(lowf.subs(be,Y*al)/al**md)/((Y*r+Y*t+r*s-t)*(-Y**2+Y*(t-u+1)+r*s-t))).is_constant())
# H6 numerator vanish at (1,1)
n6,d6=sp.fraction(sp.cancel(H[6]))
print('H6 num at (1,1):',sp.expand(n6.subs({x:1,z:1})),' den factors',sp.factor_list(d6))
# H6 along x=1+al, z=1+Y al
h6b=sp.cancel(H[6].subs({x:1+al,z:1+Y*al}))
print('H6 on slope branches limit:',sp.factor(sp.limit(h6b,al,0)))
# infinity: z=inf
print('H6 limit z->inf:',sp.factor(sp.limit(H[6],z,sp.oo)))
chk('H6 -> r + (t+s(x-1))/x as z->inf', sp.cancel(sp.limit(H[6],z,sp.oo)-(r+(t+s*(x-1))/x))==0)
print('H6 limit x->inf:',sp.factor(sp.limit(H[6],x,sp.oo)))
chk('H6 -> r+s+(r+s-t-1)/(z-1) as x->inf', sp.cancel(sp.limit(H[6],x,sp.oo)-(r+s+(r+s-t-1)/(z-1)))==0)
# k=0, z=inf double point
CKt=sp.expand(s**5*tau**3*Cn.subs({x:xs,z:1/tau}))
Pk=sp.Poly(CKt,K,tau); mdk=min(sum(m) for m in Pk.monoms())
q2=sum(c*K**m[0]*tau**m[1] for m,c in zip(Pk.monoms(),Pk.coeffs()) if sum(m)==mdk)
print('k=0,z=inf lowest degree',mdk,' quadratic',sp.factor(q2))
qa=sp.Poly(q2,K,tau)
discq=sp.factor(qa.coeff_monomial(K*tau)**2-4*qa.coeff_monomial(K**2)*qa.coeff_monomial(tau**2))
print('   coefficients K^2:',sp.factor(qa.coeff_monomial(K**2)),' tau^2:',sp.factor(qa.coeff_monomial(tau**2)),' disc:',discq)
lam=sp.Symbol('lam')
h6k=sp.cancel(H[6].subs({x:xs,z:1/tau}).subs(K,lam*tau))
print('   H6 limit along K=lam*tau:',sp.factor(sp.limit(h6k,tau,0)))
# constant-C exception
gen={r:2,s:3,t:5,u:7}
print('witness (2,3,5,7):', 'r+s-u-1',(r+s-u-1).subs(gen),' (s-1)(s-t)',((s-1)*(s-t)).subs(gen),' r+s-t-1',(r+s-t-1).subs(gen),' r+t',(r+t).subs(gen),' (s-1)(s-u)(r-s)',((s-1)*(s-u)*(r-s)).subs(gen),' (r-1)',(r-1).subs(gen),' discq',discq.subs(gen))
c1p=Cn.subs({r:0,s:1,t:0,u:0}); print('C at (0,1,0,0):',sp.factor(c1p))
chk('[x^5 z^2]C = 3 at (r,s,t,u)=(0,1,0,0)', sp.Poly(sp.expand(c1p),x,z).coeff_monomial(x**5*z**2)==3)
c2p=sp.expand(Cn.subs({r:1,s:-1,t:-1,u:-1})); chk('C == 0 at (1,-1,-1,-1)', c2p==0)
chk('H6 == 0 at (1,-1,-1,-1)', sp.cancel(H[6].subs({r:1,s:-1,t:-1,u:-1}))==0)
sol=sp.solve([r+t, r+s-t-1, (s-1)*(t-s), u-(r+s-1)],[r,s,t,u],dict=True); print('constant candidates',sol)
ok = ok and COMB_OK
print('ALL EXACT CHECKS PASS' if ok else 'SOME CHECK FAILED'); sys.exit(0 if ok else 1)
