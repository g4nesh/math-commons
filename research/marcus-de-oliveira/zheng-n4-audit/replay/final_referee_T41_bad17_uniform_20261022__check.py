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
# ---- bad17 objects (built from the displayed definitions of final_referee_T41_bad17_uniform, pp.1-2)
S=a+b; A=a*b; R=(a+1)*(b+1); X=a+q+1; Y=b+e+1; T=e+q+1; V=A+(S+1)*(p+1)
x_,y_=sp.symbols('x y')
f=(p-a+b*x_)**2-4*b*p/(b+1)*((S+1)*x_-a)
k=(p-b+a*y_)**2-4*a*p/(a+1)*((S+1)*y_-b)
r00=p**2+(2-S)*p+2-A
r=r00+a*(a+p+2)*y_+b*(b+p+2)*x_-(2*V-A)*x_*y_
R0=sp.factor(sp.together(X*Y*r.subs({x_:e/Y,y_:q/X})))
Fe=sp.expand(sp.cancel(Y**2*f.subs(x_,e/Y))); Fq=sp.expand(sp.cancel(X**2*k.subs(y_,q/X)))
_num=sp.expand(R0**2-Fe*Fq); _qq,_rr=sp.div(sp.Poly(_num,q,p,j,c,e,b,a), sp.Poly(sp.expand(4*V),q,p,j,c,e,b,a))
A0=_qq.as_expr()
D=X*Y*(c+1)+T*c*j
G17=sp.expand(V*D**2+R0*j*D+A0*j**2)
def g17_basic():
    chk("bad17 p.2: R0 = XY r(e/Y, q/X) is a polynomial", sp.denom(sp.together(R0))==1)
    chk("bad17 p.2: 4V divides R0^2 - Fe Fq, so A0 is a polynomial", _rr.is_zero)
    chk("bad17 p.1: A0 and R0 depend only on a,b,p,q,e", not (A0.has(j) or A0.has(c) or sp.sympify(R0).has(j) or sp.sympify(R0).has(c)))
    info("bad17 G = V D^2 + R0 j D + A0 j^2", "terms=%d total degree=%d"%(len(sp.Poly(G17,q,p,j,c,e,b,a).terms()), sp.Poly(G17,q,p,j,c,e,b,a).total_degree()))
# ======== final_referee_T41_bad17_uniform_20261022 ========
t0=time.time()
gen_check("bad17")
deg,Cc,Qc=CQ("bad17")
chk("bad17 p.1: cubic restricted to w = g + h kappa is linear in h", deg==1)
chk("bad17 p.1: C = -(c+1)(ab+a+ae+b+bq+1+e+j+q+ej+eq+jq)", sp.expand(Cc+(c+1)*(a*b+a+a*e+b+b*q+1+e+j+q+e*j+e*q+j*q))==0, "Q = %s"%Qc)
g17_basic()
chk("bad17 p.1: det H_F(w) = -4(c+1)^2 G^2/Q^2 at h=-C/Q (exact, 8 random rational points)", det_identity("bad17",G17,c+1,8,171)=={"1"})
Gc0=sp.expand(G17.subs(c,0))
JJ=sp.symbols('JJ')
chk("bad17 p.1: radial transformation G(c,j) = (D/(XY))^2 G(0, jXY/D)", sp.simplify(G17 - (D/(X*Y))**2*Gc0.subs(j, j*X*Y/D))==0)
chk("bad17 p.1: G(0,0)=V (XY)^2 > 0 (constant term)", sp.expand(G17.subs(j,0)-V*D.subs(j,0)**2)==0)
# symmetry claims
s_=find_maps({4:16,16:4,3:15,15:3,1:1,53:53,6:6,17:17})
chk("bad17 p.2: relabeling stabilizer swapping (a,b),(q,e) exists in G and fixes T_42", len(s_)>0 and all(M.g_act_set(g,M.CONTACT_SETS[42])==M.CONTACT_SETS[42] for g in s_), fmt(s_))
s_=find_maps({59:1,53:53,55:4,67:16,54:17})
chk("bad17 p.2: a group element maps generators [59,53,55,67,54] to [1,53,4,16,17]", len(s_)>0, fmt(s_))
# Section 2 identities (Lemma 1)
Wl=(A-1)*(p**2+1)+(S-A-1)*p
chk("bad17 p.3: r00^2 - f(0)k(0) = -4VW/R", sp.simplify(r00**2-f.subs(x_,0)*k.subs(y_,0)+4*V*Wl/R)==0)
Hh=(a-b)**2-2*p*S+p**2
chk("bad17 p.3: H = f(1) = k(1) = r(1,1)", sp.simplify(f.subs(x_,1)-Hh)==0 and sp.simplify(k.subs(y_,1)-Hh)==0 and sp.expand(r.subs({x_:1,y_:1})-Hh)==0)
fe=sp.expand(f); chk("bad17 p.3: vertex of f at [a(b+1)+p(2a+b+1)]/(b(b+1))", sp.simplify(-fe.coeff(x_,1)/(2*fe.coeff(x_,2)) - (a*(b+1)+p*(2*a+b+1))/(b*(b+1)))==0)
chk("bad17 p.3: disc_x f = 16ab^2pV/(b+1)^2 > 0 and disc_y k = 16a^2bpV/(a+1)^2 (a,b,p>0)",
    sp.simplify(sp.discriminant(sp.expand(f*(b+1)),x_)/(b+1)**2-16*a*b**2*p*V/(b+1)**2)==0 and sp.simplify(sp.discriminant(sp.expand(k*(a+1)),y_)/(a+1)**2-16*a**2*b*p*V/(a+1)**2)==0)
r1=sp.expand(r.subs(y_,1)); C1=r1.subs(x_,0); D1=sp.expand(r1.coeff(x_)); xr=-C1/D1
chk("bad17 p.4: f(x_r) = -4(Ap-ap-a-p-1)HV/((b+1)D1^2)", sp.simplify(f.subs(x_,xr)+4*(A*p-a*p-a-p-1)*Hh*V/((b+1)*D1**2))==0)
hh=sp.symbols('hh'); re_=sp.expand(r+hh)
Fm=sp.Matrix([[fe.coeff(x_,0), fe.coeff(x_,1)/2],[fe.coeff(x_,1)/2, fe.coeff(x_,2)]])
ke=sp.expand(k); Km=sp.Matrix([[ke.coeff(y_,0), ke.coeff(y_,1)/2],[ke.coeff(y_,1)/2, ke.coeff(y_,2)]])
Mm=sp.Matrix([[re_.subs({x_:0,y_:0}), sp.expand(re_.subs(x_,0)).coeff(y_,1)],[sp.expand(re_.subs(y_,0)).coeff(x_,1), re_.coeff(x_,1).coeff(y_,1)]])
Jm=Fm.det()*Km - Mm.T*Fm.adjugate()*Mm
chk("bad17 p.4 (6): J22 = -4[a-(b+1)p]^2 V^2/(b+1)^2 (independent of the level h)", sp.simplify(Jm[1,1]+4*(a-(b+1)*p)**2*V**2/(b+1)**2)==0)
xv=sp.Matrix([1,x_]); yv=sp.Matrix([1,y_]); Nm=(Mm*yv)*(Mm*yv).T - ke*Fm
chk("bad17 p.4: det of the x-quadratic matrix of E equals k(y) y^T J y (auditor check of the 'determinant' step)", sp.cancel(sp.together(Nm.det() - ke*(yv.T*Jm*yv)[0]))==0)
d1=sp.diff(f,x_).subs(x_,1)/2 - sp.diff(r,x_).subs({x_:1,y_:1}); d2=sp.diff(k,y_).subs(y_,1)/2 - sp.diff(r,y_).subs({x_:1,y_:1})
chk("bad17 p.4: edge derivatives at (1,1): 2[(a+1)+p(S+1)/(b+1)] and 2[(b+1)+p(S+1)/(a+1)]", sp.simplify(d1-2*((a+1)+p*(S+1)/(b+1)))==0 and sp.simplify(d2-2*((b+1)+p*(S+1)/(a+1)))==0)
s1,t1=sp.symbols('s t',positive=True)
chk("bad17 p.4 (8): (2b+1)(sa-sb)^2-(a-b) = 2 sb (sa-sb)(sa sb - b - 1)", sp.expand((2*t1**2+1)*(s1-t1)**2-(s1**2-t1**2)-2*t1*(s1-t1)*(s1*t1-t1**2-1))==0)
# (7),(8) stationary values: compare roots of disc_y[(hh+r(x0,y))^2 - f(x0)k(y)] with the claimed closed forms at sample points
for x0,claim,vals in ((0, lambda sg: 2*V*(-(p+1)+sg*2*sp.sqrt(b*p/(b+1)))/(a+1), {a:sp.Rational(3,7),b:sp.Rational(2,9),p:sp.Rational(5,3)}),
                      (1, lambda sg: 2*V*(-a+2*b*p+b+p+sg*2*sp.sqrt(b*p*((b+1)*(p+1)-a)))/(a*(a+1)), {a:sp.Rational(9,10),b:sp.Rational(1,10),p:sp.Rational(1,5)})):
    Ex=sp.expand(((hh+r.subs(x_,x0))**2-f.subs(x_,x0)*k).subs(vals))
    sols=sorted(float(sp.re(s2)) for s2 in sp.Poly(sp.Poly(Ex,y_).discriminant(),hh).nroots() if abs(sp.im(s2))<1e-20)
    cl=sorted(float(sp.N(claim(sg).subs(vals))) for sg in (1,-1))
    chk("bad17 p.4 (%s): stationary edge values at one rational sample (float compare of exact roots)"%("7" if x0==0 else "8"), all(abs(u-v)<1e-9*(1+abs(v)) for u,v in zip(sols,cl)), (sols,cl))
B0=(a-b)**2*e*q+(a+1)*(b**2+b+1-A)*e+(b+1)*(a**2+a+1-A)*q+R*(1-A)
chk("bad17 p.5: G(0,j) at p=0 equals (XY+jT)(RXY+jB0)", sp.expand(G17.subs({c:0,p:0})-(X*Y+j*T)*(R*X*Y+j*B0))==0)
# small-product parametrization and certificate (Section 1, p.3)
from sympy import field, QQ
Kf, U_,Z_,P_,xi,eta = field("U,Z,P,xi,eta", QQ)
rho=1+U_; Dr=rho**2-rho+1; z=(rho+Dr*Z_)/(Dr*(1+Z_)); aa=z*rho; bb=z/rho
k0=rho*U_**2/(Dr*(1+Z_)*(rho+Dr*Z_))
chk("bad17 p.3: 1-z = U^2/[Drho(1+Z)], S-1-z = Z U^2/[rho(1+Z)]", 1-z==U_**2/(Dr*(1+Z_)) and aa+bb-1-z==Z_*U_**2/(rho*(1+Z_)))
chk("bad17 p.2-3: kappa0 = (1-z)^2/(S-2z) = a(1-z)^2/(a-z)^2 = rho U^2/[Drho(1+Z)(rho+Drho Z)]", (1-z)**2/(aa+bb-2*z)==k0 and aa*(1-z)**2/(aa-z)**2==k0)
chk("bad17 p.3: inverse Z = (z Drho - rho)/[Drho(1-z)]", (z*Dr-rho)/(Dr*(1-z))==Z_)
xx=xi/(1+xi); yy=eta/(1+eta); SS=aa+bb; AA=aa*bb; RR=(aa+1)*(bb+1)
def build(pv, layer):
    VV=AA+(SS+1)*(pv+1)
    ff=(pv-aa+bb*xx)**2-4*bb*pv/(bb+1)*((SS+1)*xx-aa)
    kk=(pv-bb+aa*yy)**2-4*aa*pv/(aa+1)*((SS+1)*yy-bb)
    rr0=pv**2+(2-SS)*pv+2-AA
    rr=rr0+aa*(aa+pv+2)*yy+bb*(bb+pv+2)*xx-(2*VV-AA)*xx*yy
    hv=2*VV*P_ if layer==1 else Kf(0)
    HH=RR*(hv+rr); EE=HH**2-RR*((bb+1)*ff)*((aa+1)*kk)
    return HH,EE
for layer,pv,claims in ((1,k0+P_,(1360,14373)),(2,k0*P_/(1+P_),(1330,10831))):
    for nm,ex,cl in zip(("H=R(h+r)","E=H^2-R fN kN"), build(pv,layer), claims):
        cn=list(ex.numer.coeffs()); cd=list(ex.denom.coeffs())
        ok=(all(x>0 for x in cn) and all(x>0 for x in cd)) or (all(x<0 for x in cn) and all(x<0 for x in cd))
        chk("bad17 p.3 small-product layer %d %s: numerator and denominator one-signed with equal sign"%(layer,nm), ok, "numerator terms %d (claimed %d), denominator terms %d"%(len(cn),cl,len(cd)))
# numeric falsification
import numpy as np
rng=np.random.default_rng(2026)
def fkr(a_,b_,p_,x,y):
    S_=a_+b_; A_=a_*b_; V_=A_+(S_+1)*(p_+1)
    f_=(p_-a_+b_*x)**2-4*b_*p_/(b_+1)*((S_+1)*x-a_); k_=(p_-b_+a_*y)**2-4*a_*p_/(a_+1)*((S_+1)*y-b_)
    r_=p_**2+(2-S_)*p_+2-A_+a_*(a_+p_+2)*y+b_*(b_+p_+2)*x-(2*V_-A_)*x*y
    return f_,k_,r_,V_
def ls(n,lo,hi): return np.exp(rng.uniform(np.log(lo),np.log(hi),n))
N=1_000_000; worst=-np.inf; cnt=0
for rep in range(4):
    s_=rng.uniform(0,2,N)*rng.uniform(0,1,N)**0.3; t_=rng.uniform(0,1,N); av=s_*t_; bv=s_*(1-t_)
    m_=rng.random(N); av=np.where(m_<0.05,0,av); bv=np.where((m_>0.05)&(m_<0.1),0,bv)
    pv=ls(N,1e-4,1e3); x=rng.uniform(0,1,N); y=rng.uniform(0,1,N)
    mm=rng.random(N); x=np.where(mm<0.1,1.0,np.where(mm<0.2,0.0,x)); y=np.where((mm>0.3)&(mm<0.4),1.0,np.where((mm>0.4)&(mm<0.5),0.0,y))
    f_,k_,r_,V_=fkr(av,bv,pv,x,y); ok=(f_>=0)&(k_>=0); cnt+=ok.sum()
    worst=max(worst,((np.sqrt(np.clip(f_*k_,0,None))-r_)/(1+np.abs(r_)+np.sqrt(np.abs(f_*k_))))[ok].max())
info("NUMERIC Lemma 1 (5), seed 2026, a+b<2, p in [1e-4,1e3], (x,y) in D incl. edges", "samples in D=%d, max normalized sqrt(fk)-r = %.3e (<= ~1e-8 means no violation)"%(cnt,worst))
w1=w2=np.inf
for rep in range(4):
    bv=ls(N,1e-3,1e3); Av=1+ls(N,1e-6,1e3)*(rng.random(N)<0.9); av=Av/bv; pv=ls(N,1e-5,1e3)*(rng.random(N)<0.95)
    x=rng.uniform(0,1,N); y=rng.uniform(0,1,N); f_,k_,r_,V_=fkr(av,bv,pv,x,y); hv=2*V_*(pv+(Av-1)/((av+1)*(bv+1))); sc=1+np.abs(hv)+np.abs(r_)
    w1=min(w1,((hv+r_)/sc).min()); w2=min(w2,(((hv+r_)**2-f_*k_)/sc**2).min())
info("NUMERIC large-product claim (cited doc), A>=1", "min (h+r)/scale=%.3e, min ((h+r)^2-fk)/scale^2=%.3e"%(w1,w2))
info("bad17 total time", "%.1fs"%(time.time()-t0))

print("SUMMARY: %d checks failed"%len(FAIL), FAIL)
sys.exit(1 if FAIL else 0)
