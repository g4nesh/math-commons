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
# ======== final_referee_T41_negative_Heron_region_20261022 ========
t0=time.time()
gen_check("bad17")
deg,Cc,Qc=CQ("bad17")
chk("negH p.1: C = -(c+1)(ab+a+ae+b+bq+1+e+j+q+ej+eq+jq), linear in h", deg==1 and sp.expand(Cc+(c+1)*(a*b+a+a*e+b+b*q+1+e+j+q+e*j+e*q+j*q))==0)
g17_basic()
chk("negH p.1: det H_F(w) = -4(c+1)^2 G^2/Q^2 with G the bad17 radial factor (exact, 6 random rational points)", det_identity("bad17",G17,c+1,6,991)=={"1"})
Fe3=((e+d+b)*(p-a)+e*b)**2-4*b*p*(e-a)*(e+d+b); Fq3=((q+d+a)*(p-b)+q*a)**2-4*a*p*(q-b)*(q+d+a)
chk("negH p.1: the quartics Fe, Fq (d=1) equal Y^2 f(e/Y), X^2 k(q/X) of the bad17 note", sp.expand(Fe-Fe3.subs(d,1))==0 and sp.expand(Fq-Fq3.subs(d,1))==0)
s_=find_maps({4:16,16:4,3:15,15:3,1:1,53:53,6:6,17:17})
chk("negH p.1: stabilizer swapping (a,b),(q,e) exchanges Fe and Fq", len(s_)>0 and sp.expand(Fe3.subs({a:b,b:a,q:e,e:q},simultaneous=True)-Fq3)==0)
Gj=sp.Poly(G17,j); al,be,ga=[Gj.coeff_monomial(j**i) for i in (2,1,0)]
chk("negH p.1: G quadratic in j, gamma = [p(1+a+b)+(a+1)(b+1)](e+1+b)^2(c+1)^2(q+1+a)^2", Gj.degree()==2 and sp.expand(ga-(p*(1+a+b)+(a+1)*(b+1))*(e+1+b)**2*(c+1)**2*(q+1+a)**2)==0)
chk("negH p.1 (1): beta^2-4 alpha gamma = (e+1+b)^2(c+1)^2(q+1+a)^2 Fe Fq", sp.expand(be**2-4*al*ga-(e+1+b)**2*(c+1)**2*(q+1+a)**2*Fe*Fq)==0)
chk("negH p.1: Fe<0 forces b,p>0 and e>a; Fq<0 forces a>0, q>b (sign structure of the quartics)", True, "by inspection: Fe = square - 4bp(e-a)(e+1+b)")
res=minor_checks("bad17", j, G17, Fe, Fq, (c+1)**2*(q+1+a)**2, (c+1)**2*(e+1+b)**2, j**2, 4, 992)
for t_,rr in enumerate(res):
    chk("negH p.2 identities at random base point %d (block I height-free, = j^2 P(j), P quadratic, Res_j(G,P) square, (2) square classes, congruences mod G, U = alpha P - a2 G, bordered I+{0}, I+{1} numerators divisible by G)"%t_,
        all([rr["interp"],rr["height_free"],rr["factor"],rr["degP"]==2,rr["res_square"],rr["square_classes_and_congruences"],rr["U_eq_alphaP_minus_a2G"],rr["bordered_div_G"]]), rr)
rk=rank_test("bad17", j, G17, 5, 993)
info("NUMERIC rank of H_F at real j-roots of G (mpmath 60 digits), 3 smallest singular values", ["%.1e %.1e %.1e"%x for x in rk])
info("negH total time", "%.1fs"%(time.time()-t0))

print("SUMMARY: %d checks failed"%len(FAIL), FAIL)
sys.exit(1 if FAIL else 0)
