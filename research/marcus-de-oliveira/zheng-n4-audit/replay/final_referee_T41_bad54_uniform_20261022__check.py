# Own-code replay for the Zheng n=4 audit (batch C17-T42). Nothing from the supplement is executed or read.
# Imports only the audit's own library (../lib/moc4lib.py, or the directory in env MOC4LIB).
import sys, os, random, math, time
from fractions import Fraction as Fr
LIB = os.environ.get("MOC4LIB") or os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib")
sys.path.insert(0, LIB)
import moc4lib as M
import sympy as sp
import flint
FAIL = []
def chk(name, ok, extra=""):
    print(("PASS " if ok else "FAIL ") + name + ((" | " + str(extra)) if extra != "" else ""), flush=True)
    if not ok: FAIL.append(name)
def info(name, extra): print("INFO " + name + " | " + str(extra), flush=True)
q,p,j,c,e,d,b,a,h = sp.symbols('q p j c e d b a h')
SRC={"bad17":[(q,3),(p,1),(j,53),(c,6),(e,15),(d,17),(b,16),(a,4)],
     "bad54":[(q,6),(p,2),(j,5),(c,53),(e,0),(d,54),(b,3),(a,15)]}
def src_w(name, dnorm=True):
    w=[sp.Integer(0)]*14
    for coef,k in SRC[name]:
        for i in range(14): w[i]+=coef*sp.Rational(M.F_VECS[k][i].numerator,M.F_VECS[k][i].denominator)
    w[13]+=h
    return [sp.expand(x.subs(d,1)) if dnorm else sp.expand(x) for x in w]
ms=sp.symbols('m0:14')
Cm=sp.expand(M.cubic_C(ms))
HCOEF=[[[int(sp.diff(Cm,ms[i],ms[k],ms[l])) for l in range(14)] for k in range(14)] for i in range(14)]
def hess_at(wv):
    return [[sum(HCOEF[i][k][l]*wv[l] for l in range(14) if HCOEF[i][k][l]) for k in range(14)] for i in range(14)]
def Q_(x): return flint.fmpq(int(sp.numer(x)),int(sp.denom(x)))
def fdet(Mx):
    return flint.fmpq_mat([[x if isinstance(x,flint.fmpq) else flint.fmpq(x.numerator,x.denominator) for x in r] for r in Mx]).det()
def CQ(name):
    w=src_w(name); Cw=sp.Poly(sp.expand(M.cubic_C(w)),h)
    return Cw.degree(), sp.expand(Cw.coeff_monomial(1)), sp.expand(Cw.coeff_monomial(h))
def gen_check(name):
    T=M.CONTACT_SETS[42]; gens=M.generators(T)
    ks=[k for _,k in SRC[name]]
    chk(name+": all 8 source generators are generators of paper T_42 = %s"%sorted(T), set(ks)<=set(gens), "F_T42=%s"%gens)
    nonminor=[k for k in ks if not M.F_IS_MINOR[k]]
    chk(name+": exactly one non-minor generator in the source", len(nonminor)==1, "non-minor: %s"%[(k,M.f_name(k)) for k in nonminor])
    chk(name+": the 8 generators are linearly independent; with kappa rank 9", M.rank([M.F_VECS[k] for k in ks])==8 and M.rank([M.F_VECS[k] for k in ks]+[M.KAPPA])==9)
VAL2IDX={M.F_VALUES[k]:k for k in range(82)}
def gen_map(g):
    t=M.g_index_map(g); out={}
    for k in range(82):
        nv=[None]*24
        for s_ in range(24): nv[t[s_]]=M.F_VALUES[k][s_]
        out[k]=VAL2IDX.get(tuple(nv))
    return out
def find_maps(mapping):
    return [g for g in M.G_ELEMENTS if all(gen_map(g)[k]==v for k,v in mapping.items())]
def fmt(gs): return [("".join(map(str,r1)),"".join(map(str,r2)),eps) for r1,r2,eps in gs]
def det_identity(name, G, pref, npts, seed):
    rnd=random.Random(seed); w=src_w(name); _,Cc,Qc=CQ(name); ratios=set(); used=0
    while used<npts:
        vals={v:sp.Rational(rnd.randint(0,40),rnd.randint(1,9)) for v in (q,p,j,c,e,b,a)}
        Qv=Qc.subs(vals)
        if Qv==0: continue
        used+=1; vals[h]=-Cc.subs(vals)/Qv
        wv=[Q_(sp.Rational(xx.subs(vals))) for xx in w]
        dH=fdet(hess_at(wv)); rhs=Q_(-4*pref.subs(vals)**2*G.subs(vals)**2/Qv**2)
        ratios.add(str(dH/rhs) if rhs!=0 else ("0/0" if dH==0 else "bad"))
    return ratios
def topoly(expr, var):
    P_=sp.Poly(expr,var); cs=P_.all_coeffs()[::-1]
    return flint.fmpq_poly([Q_(x) for x in cs])
def interp(xs, ys):
    res=flint.fmpq_poly([0])
    for i in range(len(xs)):
        num=flint.fmpq_poly([1]); den=flint.fmpq(1)
        for k in range(len(xs)):
            if k!=i: num*=flint.fmpq_poly([-xs[k],1]); den*=(xs[i]-xs[k])
        res+=num*(ys[i]/den)
    return res
def is_sq(x):
    if x<0: return False
    n,dn=int(x.p),int(x.q); return math.isqrt(n)**2==n and math.isqrt(dn)**2==dn
def sqrtq(x): return flint.fmpq(math.isqrt(int(x.p)),math.isqrt(int(x.q)))
I12=list(range(2,14))
def minor_checks(name, var, G, Fp, Fm, mp, mm, factor_pow, npts, seed):
    """At random rational base points (all variables but var), univariate exact checks of the
    height-free principal block on I={2..13}, the resultant/square-class/congruence identities and the bordered minors."""
    rnd=random.Random(seed); w=src_w(name); _,Cc,Qc=CQ(name)
    base=[v for v in (q,p,j,c,e,b,a) if v!=var]; res=[]
    for t in range(npts):
        bv={v:sp.Rational(rnd.randint(1,30),rnd.randint(1,7)) for v in base}
        wsub=[sp.expand(x.subs(bv)) for x in w]
        # each entry is affine in (var, h): exact coefficients
        lin=[(Q_(x.subs({var:0,h:0})), Q_(sp.expand(x).coeff(var)), Q_(sp.expand(x).coeff(h))) for x in wsub]
        assert all(sp.Poly(x,var,h).total_degree()<=1 for x in wsub)
        def wvals(xv,hv): return [c0_+c1_*xv+c2_*hv for c0_,c1_,c2_ in lin]
        Gp=topoly(sp.expand(G.subs(bv)),var)
        xs=[flint.fmpq(s_) for s_ in range(1,17)]; ys=[]
        for xv in xs:
            Hm=hess_at(wvals(xv,flint.fmpq(0))); ys.append(fdet([[Hm[i][k] for k in I12] for i in I12]))
        Bp=interp(xs[:15],ys[:15]); interp_ok=(Bp(xs[15])==ys[15])
        # height-free: same block at a different h
        Hm=hess_at(wvals(xs[3],flint.fmpq(7,3))); hf=(fdet([[Hm[i][k] for k in I12] for i in I12])==ys[3])
        fac=topoly(sp.expand(factor_pow),var); Pq,rem=divmod(Bp,fac)
        al,be,ga=[Gp[i] for i in (2,1,0)]; a2,b1,c0=[Pq[i] for i in (2,1,0)]
        u1=al*b1-a2*be; u0=al*c0-a2*ga; Ucal=flint.fmpq_poly([u0,u1])
        Res=flint.fmpq_mat([[al,be,ga,0],[0,al,be,ga],[a2,b1,c0,0],[0,a2,b1,c0]]).det()
        Tt=2*al*u0-be*u1; found=None
        if is_sq(Res):
            for sgn in (1,-1):
                W=sgn*sqrtq(Res); Sp=Tt+2*al*W; Sm=Tt-2*al*W
                fp=Q_(Fp.subs(bv)); fm=Q_(Fm.subs(bv)); mpv=Q_(mp.subs(bv)); mmv=Q_(mm.subs(bv))
                s1=(fp!=0 and is_sq(Sp/(fp*mpv))); s2=(fm!=0 and is_sq(Sm/(fm*mmv)))
                c1=divmod(Pq*Sp-(Ucal+W)**2, Gp)[1]==0; c2=divmod(Pq*Sm-(Ucal-W)**2, Gp)[1]==0
                if s1 and s2 and c1 and c2: found=sgn
        Ueq=(Ucal-(Pq*al-Gp*a2))==0
        Cp=topoly(sp.expand(Cc.subs(bv)),var); Qp=topoly(sp.expand(Qc.subs(bv)),var); bord=[]
        for kk in (0,1):
            idx=I12+[kk]; xs2=[flint.fmpq(s_) for s_ in range(1,41)]; ys2=[]
            for xv in xs2:
                Qv=Qp(xv); hv=-Cp(xv)/Qv; Hm=hess_at(wvals(xv,hv))
                ys2.append(fdet([[Hm[i][k] for k in idx] for i in idx])*Qv**2)
            Np=interp(xs2[:38],ys2[:38]); okint=(Np(xs2[38])==ys2[38] and Np(xs2[39])==ys2[39])
            bord.append(okint and divmod(Np,Gp)[1]==0)
        res.append(dict(interp=interp_ok, height_free=hf, factor=(rem==0), degP=Pq.degree(), res_square=is_sq(Res),
                        square_classes_and_congruences=(found is not None), U_eq_alphaP_minus_a2G=Ueq, bordered_div_G=all(bord)))
    return res
def rank_test(name, var, G, n, seed):
    import mpmath; mpmath.mp.dps=60
    rnd=random.Random(seed); w=src_w(name); _,Cc,Qc=CQ(name); out=[]; tries=0
    allv=(q,p,j,c,e,b,a)
    Cl=sp.lambdify(allv,Cc,"mpmath"); Ql=sp.lambdify(allv,Qc,"mpmath"); wl=[sp.lambdify(allv+(h,),xx,"mpmath") for xx in w]
    while len(out)<n and tries<300:
        tries+=1
        bv={v:sp.Rational(rnd.randint(1,40),rnd.randint(1,8)) for v in allv if v!=var}
        roots=[rr for rr in sp.Poly(sp.expand(G.subs(bv)),var).nroots(n=50) if abs(sp.im(rr))<1e-40]
        if not roots: continue
        vals={kk:mpmath.mpf(int(sp.numer(vv)))/int(sp.denom(vv)) for kk,vv in bv.items()}; vals[var]=mpmath.mpf(str(sp.re(roots[0])))
        args=[vals[v] for v in allv]; hv=-Cl(*args)/Ql(*args); wv=[f_(*args,hv) for f_ in wl]
        Hm=mpmath.matrix([[sum(HCOEF[i][kk][l]*wv[l] for l in range(14)) for kk in range(14)] for i in range(14)])
        svs=sorted(abs(s_) for s_ in mpmath.svd_r(Hm,compute_uv=False))
        out.append((float(svs[0]),float(svs[1]),float(svs[2])))
    return out
# ======== final_referee_T41_bad54_uniform_20261022 ========
t0=time.time()
gen_check("bad54")
s_=find_maps({2:5,5:2,15:3,3:15,6:6,53:53,0:0,54:54})
chk("bad54 p.1: stabilizer exchanging (p,j),(a,b) and fixing f6,f53,f0,f54 exists and fixes T_42", len(s_)>0 and all(M.g_act_set(g,M.CONTACT_SETS[42])==M.CONTACT_SETS[42] for g in s_), fmt(s_))
deg,Cc,Qc=CQ("bad54")
chk("bad54 p.1: cubic restriction linear in h", deg==1)
chk("bad54 p.1: C = -(q+1)(ab+ac+a+bc+b+c+ce+1+e)", sp.expand(Cc+(q+1)*(a*b+a*c+a+b*c+b+c+c*e+1+e))==0)
chk("bad54 p.1: Q = -ab-ac-a+ap-bc-b+bj-c-1-jp", sp.expand(Qc-(-a*b-a*c-a+a*p-b*c-b+b*j-c-1-j*p))==0)
T5=e+a+b+1; U5=T5+a*b; W5=U5+c*T5; V5=p*j+(p+j+1)*(e+1); H5=(p-j)**2-2*e*(p+j)+e**2
A1=e**2-e*j+e-j*p-j+p**2+p; B1=e**2-e*p+e+j**2-j*p+j-p
K5=V5-e*(e+1); B05=-H5*a*b+p*B1*a+j*A1*b+p*j*K5; C05=T5*(p*j-p*a-j*b)+e*a*b
L=(q+1)*(T5+a*b/(c+1))
G54=sp.expand(sp.cancel((c+1)**2*(V5*L**2-B05*L-p*j*e*C05)))
nt=len(sp.Poly(G54,q,p,j,c,e,b,a).terms())
chk("bad54 p.2 (3): (c+1)^2{VL^2-B0L-pjeC0} is a polynomial of degree 8 with 1011 terms", sp.denom(sp.together(G54))==1 and nt==1011 and sp.Poly(G54,q,p,j,c,e,b,a).total_degree()==8, "terms=%d"%nt)
chk("bad54 p.1 (2): det H_F(w) = -4(q+1)^2 G^2/Q^2 at h=-C/Q (exact, 8 random rational points)", det_identity("bad54",G54,q+1,8,541)=={"1"})
Fa=(j*(p+e+a+1)-a*(p+e))**2-4*p*e*a*(a-j); Fb=(p*(j+e+b+1)-b*(j+e))**2-4*j*e*b*(b-p)
chk("bad54 p.2 (4): B0^2 + 4VpjeC0 = Fa Fb", sp.expand(B05**2+4*V5*p*j*e*C05-Fa*Fb)==0)
Gq=sp.Poly(G54,q); aq,bq,gq=[Gq.coeff_monomial(q**i) for i in (2,1,0)]
chk("bad54 p.2 (4): Disc_q G = (c+1)^2 W^2 Fa Fb and alpha_q = V W^2", Gq.degree()==2 and sp.expand(bq**2-4*aq*gq-(c+1)**2*W5**2*Fa*Fb)==0 and sp.expand(aq-V5*W5**2)==0)
Gc=sp.Poly(G54,c); ac,bc,gc=[Gc.coeff_monomial(c**i) for i in (2,1,0)]
chk("bad54 p.2 (4): Disc_c G = a^2 b^2 (q+1)^2 Fa Fb", Gc.degree()==2 and sp.expand(bc**2-4*ac*gc-a**2*b**2*(q+1)**2*Fa*Fb)==0)
chk("bad54 p.2: Disc_a Fa = 16epj^2V, Disc_b Fb = 16ejp^2V", sp.expand(sp.discriminant(Fa,a)-16*e*p*j**2*V5)==0 and sp.expand(sp.discriminant(Fb,b)-16*e*j*p**2*V5)==0)
LL=sp.symbols('LL'); Phi=V5*LL**2-B05*LL-p*j*e*C05
F5=sp.expand(Phi.subs(LL,T5)); D5=sp.expand(sp.diff(Phi,LL).subs(LL,T5))
chk("bad54 p.3: D T = V T^2 + F + pjeC0", sp.expand(D5*T5-(V5*T5**2+F5+p*j*e*C05))==0)
chk("bad54 p.3: C0 + T(pa+jb) = Tpj + eab (so C0 >= -T(pa+jb))", sp.expand(C05+T5*(p*a+j*b)-(T5*p*j+e*a*b))==0)
chk("bad54 p.3: V - p^2 je = ep(1-pj) + e(j+1) + p + j + 1 + pj", sp.expand(V5-p**2*j*e-(e*p*(1-p*j)+e*(j+1)+p+j+1+p*j))==0)
chk("bad54 p.2: (c+1)(L - T) is a polynomial with nonnegative coefficients", all(cc>0 for cc in sp.Poly(sp.cancel((c+1)*(L-T5)),q,p,j,c,e,b,a).coeffs()))
# Lemma 2 certificate tables (5)
from sympy import field, QQ
Kf, E_,P_,J_,A_,B_ = field("E,P,J,A,B", QQ)
def Rfun(p_,j_,e_,a_,b_):
    T=e_+a_+b_+1; V=p_*j_+(p_+j_+1)*(e_+1); H=(p_-j_)**2-2*e_*(p_+j_)+e_**2
    A1=e_**2-e_*j_+e_-j_*p_-j_+p_**2+p_; B1=e_**2-e_*p_+e_+j_**2-j_*p_+j_-p_
    Kk=V-e_*(e_+1); B0=-H*a_*b_+p_*B1*a_+j_*A1*b_+p_*j_*Kk; C0=T*(p_*j_-p_*a_-j_*b_)+e_*a_*b_
    F=V*T**2-B0*T-p_*j_*e_*C0
    Fa=(j_*(p_+e_+a_+1)-a_*(p_+e_))**2-4*p_*e_*a_*(a_-j_); Fb=(p_*(j_+e_+b_+1)-b_*(j_+e_))**2-4*j_*e_*b_*(b_-p_)
    return F-b_*Fa-a_*Fb
E0=1+E_; N0=1+E0**2*P_
rows=[(E_/(1+E_), P_/(1+P_), J_/(1+J_), A_, B_, 213, False, False),
 (E_/(1+E_), 1+P_, J_/((1+P_)*(1+J_)), A_/(P_*(1+A_)), B_, 473, True, False),
 (E0, P_/(E0*(1+P_)), J_/(E0*(1+J_)), A_, B_, 372, False, False),
 (E0, N0/(E0*(1+P_)), J_/(E0*(1+J_)), E0*(1+P_)*A_/((E0**2-1)*P_*(1+A_)), B_, 957, True, False),
 (E0, E0+P_, J_/((E0+P_)*(1+J_)), E0*A_/((E0*(E0+P_)-1)*(1+A_)), B_, 1161, True, False),
 (E0, N0/(E0*(1+P_)), (N0+E0**2*(1+P_)*J_)/(E0*N0*(1+J_)), E0*(1+P_)*A_/((E0**2-1)*P_*(1+A_)), E0*N0*(1+J_)*B_/((E0**2-1)*J_*(1+B_)), 4135, True, True)]
for i,(ev,pv,jv,av,bv,claimed,abd,bbd) in enumerate(rows,1):
    Rs=Rfun(pv,jv,ev,av,bv); cn=list(Rs.numer.coeffs()); cd=list(Rs.denom.coeffs())
    ok1=(all(x>0 for x in cn) and all(x>0 for x in cd)) or (all(x<0 for x in cn) and all(x<0 for x in cd))
    kap = Kf(1) if i<=2 else 1/ev
    da = (pv-kap) if abd else Kf(0); db = (jv-kap) if bbd else Kf(0)
    lb=(ev+1)**2*(jv+1)*(pv+1)*(1+ev-pv*jv)*(1-da*av)**2*(1-db*bv)**2
    dif=Rs-lb; dn=list(dif.numer.coeffs()); dd=list(dif.denom.coeffs())
    ok2=(not dn) or (all(x>=0 for x in dn) and all(x>0 for x in dd)) or (all(x<=0 for x in dn) and all(x<0 for x in dd))
    chk("bad54 p.2 (5) table row %d: R positive certificate and R - R00(1-da a)^2(1-db b)^2 has one-signed coefficients"%i, ok1 and ok2,
        "R numerator terms %d (claimed %d)"%(len(cn),claimed))
# Lemma 3 identities (variable c) at random base points
res=minor_checks("bad54", c, G54, Fa, Fb, b**2*(q+1)**2, a**2*(q+1)**2, (c+1)**2, 4, 542)
for t_,rr in enumerate(res):
    chk("bad54 p.3 Lemma 3 identities at random base point %d (block I height-free, =(c+1)^2 P(c), P quadratic, Res_c(G,P) square, S+-=Fa b^2(q+1)^2 L^2 / Fb a^2(q+1)^2 M^2, P S+- = (U1+-W)^2 mod G, U1 = alpha P - a2 G; plus bordered I+{0}, I+{1} numerators divisible by G)"%t_,
        all([rr["interp"],rr["height_free"],rr["factor"],rr["degP"]==2,rr["res_square"],rr["square_classes_and_congruences"],rr["U_eq_alphaP_minus_a2G"],rr["bordered_div_G"]]), rr)
rk=rank_test("bad54", c, G54, 5, 543)
info("NUMERIC rank of H_F at real c-roots of G (mpmath 60 digits), 3 smallest singular values", ["%.1e %.1e %.1e"%x for x in rk])
import numpy as np
rng=np.random.default_rng(54)
def ls(n,lo,hi): return np.exp(rng.uniform(np.log(lo),np.log(hi),n))
N=1_000_000; wst=np.inf; tot=0
for rep in range(4):
    ev=ls(N,1e-4,1e4)*(rng.random(N)<0.95); pv=ls(N,1e-4,1e4); jv=ls(N,1e-4,1e4); av=ls(N,1e-4,1e4)*(rng.random(N)<0.95); bv=ls(N,1e-4,1e4)*(rng.random(N)<0.95)
    fail=(pv*jv<1)&(av*(pv-1)<1)&(bv*(jv-1)<1)&~((av>0)&(ev>0)&(pv*av*ev>=av+ev))&~((bv>0)&(ev>0)&(jv*bv*ev>=bv+ev)); tot+=fail.sum()
    Rv=Rfun(pv,jv,ev,av,bv)
    kap=np.minimum(1,1/np.where(ev>0,ev,1e-300)); da=np.maximum(pv-kap,0); db=np.maximum(jv-kap,0)
    lbv=(ev+1)**2*(jv+1)*(pv+1)*(1+ev-pv*jv)*(1-da*av)**2*(1-db*bv)**2
    wst=min(wst,((Rv-lbv)/(np.abs(Rv)+np.abs(lbv)))[fail].min())
info("NUMERIC (5) under failure of all capacities (1), seed 54", "samples=%d, min (R-lower)/scale=%.2e (roundoff level)"%(tot,wst))
info("bad54 total time", "%.1fs"%(time.time()-t0))

print("SUMMARY: %d checks failed"%len(FAIL), FAIL)
sys.exit(1 if FAIL else 0)
