#!/usr/bin/env python3
"""Own-code replay for document final_referee_T25_uniform_20261022 (Zheng n=4 audit). Run: python -I -B this.py
Reads the supplement location from MOC4_SUPP (only JSON data files are loaded, as data).
Prints one line per check; exits nonzero if an EXACT check fails."""
import os, sys
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE,'..','lib'))
SUPP=os.environ.get('MOC4_SUPP') or (sys.argv[1] if len(sys.argv)>1 else None)
import sys, json, functools, random, time
from fractions import Fraction as Fr
import sympy as sp
import moc4lib as L
print=functools.partial(print,flush=True)
def gimage(g,k):
    """index of f_j equal to g.f_k (values transported by g)."""
    mp=L.g_index_map(g); v=L.F_VALUES[k]; w=[None]*24
    for s in range(24): w[mp[s]]=v[s]
    w=tuple(w)
    return L.F_VALUES.index(w) if w in L.F_VALUES else None
def find_sym(T, fixed_pairs):
    T=frozenset(T); out=[]
    for g in L.G_ELEMENTS:
        if L.g_act_set(g,T)!=T: continue
        if all(gimage(g,k)==m for k,m in fixed_pairs): out.append(g)
    return out
MS,CUB=L.cubic_sympy()
HESS=sp.hessian(CUB,MS)
PIVINV=sp.Matrix([[sp.Rational(x.numerator,x.denominator) for x in row] for row in L._P_PIV_INV])
def cvec_sym(vals):
    return list(PIVINV*sp.Matrix([vals[r] for r in L.P_PIVOT_ROWS]))
def check_inP(vals,c):
    Pm=sp.Matrix(L.P); return all(sp.expand(e)==0 for e in (Pm*sp.Matrix(c)-sp.Matrix(vals)))
def load_poly(path, names, deg=None):
    d=json.load(open(path))
    var=[sp.Symbol(n) for n in d['variables']]
    out={}
    for f in d['radial_factors']:
        out[f['degree']]=sp.Poly(sum(sp.Integer(int(co))*sp.prod([var[i]**e for i,e in enumerate(ex)]) for ex,co in f['terms']),*var)
    return d,var,out
def exact_det(M):
    return L.det([[Fr(int(sp.numer(x)),int(sp.denom(x))) for x in row] for row in M.tolist()])

ok=True
from flint import fmpz_mpoly_ctx
supp=SUPP
T={2,3,5,7,8,12}
ok=True
def chk(name,cond):
    global ok; ok&=bool(cond); print(('PASS ' if cond else 'FAIL ')+name)
chk('generators(T25-doc) = {0,1,2,5,6,39,67,70}', L.generators(T)==[0,1,2,5,6,39,67,70])
order=[0,1,2,5,6,39,67,70]
chk('rank 8', sp.Matrix([L.F_VALUES[k] for k in order]).rank()==8)
a,e,p,j,b,c,q,d,h=sp.symbols('a e p j b c q d z')
coef=[a,e,p,j,b,c,q,d]
vals=[sum(cf*L.F_VALUES[k][i] for cf,k in zip(coef,order)) for i in range(24)]
cv=cvec_sym(vals); chk('symbolic preimage', check_inP(vals,cv))
y=[cv[i]+(h if i==13 else 0) for i in range(14)]
Cy=sp.expand(CUB.subs(dict(zip(MS,y))).subs(d,1))
Cst=(1+q)*(a+1)*(b+1)+c*(a+1+e)+e*(1+q); Q=c*(p-1-e)-(1+e)-p*j
lam=sp.cancel(Cy/(-Cst+h*Q)); chk('(2) C(g+h kappa)|_{d=1} = lambda(-C*+hQ), lambda=%s'%lam, lam.is_number and lam!=0)
dj,var,polys=load_poly(supp+'/proofs/manin_transfer_20261021/final_referee_T25_full_radial_20261022.json',None)
G=polys[9]; print('   JSON vars',var,'terms',len(G.terms()),'deg',G.total_degree())
chk('JSON G has 1201 terms, homogeneous deg 9', len(G.terms())==1201 and all(sum(m)==9 for m in G.monoms()))
S={n:sp.Symbol(n) for n in ['q','p','j','c','e','d','b','a','z']}
Gx=G.as_expr().subs(S['z'],0)
G1=sp.expand(Gx.subs(S['d'],1))
rnd=random.Random(77); ratios=set()
for trial in range(6):
    pt={x:sp.Rational(rnd.randint(0,9),rnd.randint(1,5)) for x in [a,e,p,j,b,c,q]}; pt[d]=1
    Qv=Q.subs(pt)
    if Qv==0: continue
    pt[h]=Cst.subs(pt)/Qv
    yv=[sp.nsimplify(yy.subs(pt)) for yy in y]
    dt=exact_det(HESS.subs(dict(zip(MS,yv))))
    Gv=G1.subs({S[str(k)]:v for k,v in pt.items() if str(k) in 'aepjbcq'})
    ratios.add(sp.Rational(dt.numerator,dt.denominator)*Qv**4/(-4*Qv**2*Gv**2))
print('   ratios',ratios)
chk('(3) det Hess*Q^4 = -4Q^2G^2 at 6 random rational points', ratios=={1})
Gp=sp.Poly(G1,S['p']); k0=Gp.coeff_monomial(1)
G6=sp.cancel(k0/Cst.subs({a:S['a'],b:S['b'],c:S['c'],e:S['e'],q:S['q']}))
chk('(13) C* divides k0', sp.denom(sp.together(G6))==1)
G6=sp.expand(G6); P6=sp.Poly(G6,S['e'],S['c'],S['j'])
chk('G6 multiaffine in e,c,j; total degree %d'%sp.Poly(G6).total_degree(), all(max(m)<=1 for m in P6.monoms()))
chk('G6(e=c=j=0) = (a+1)^2(b+1)^2(q+1)^2', sp.expand(G6.subs({S['e']:0,S['c']:0,S['j']:0})-(S['a']+1)**2*(S['b']+1)**2*(S['q']+1)**2)==0)
A_,B_,Q_,J_,E_,C_=S['a'],S['b'],S['q'],S['j'],S['e'],S['c']
chk('G6(b=0) = (e+1+a)(q+c+1)(qj+q+qa+j+ja+1+a)', sp.expand(G6.subs(B_,0)-(E_+1+A_)*(Q_+C_+1)*(Q_*J_+Q_+Q_*A_+J_+J_*A_+1+A_))==0)
Lc=(A_+1)*(B_+1); dl=Lc-E_*(A_*B_-1)
chk('(15) G(j=c=0) = (a+1)(e+L)(q+1)^2{(e+L)[pq+(p+q)(b+1)]+(b+1)delta}', sp.expand(G1.subs({J_:0,C_:0})-(A_+1)*(E_+Lc)*(Q_+1)**2*((E_+Lc)*(S['p']*Q_+(S['p']+Q_)*(B_+1))+(B_+1)*dl))==0)
tt=sp.Symbol('t'); et=tt**2+tt**2*(A_+B_+2*tt)/(A_*B_-tt**2)
chk("e'(t) formula", sp.cancel(sp.diff(et,tt)-2*tt*((A_*B_-tt**2)**2+A_*B_*(A_+B_+3*tt)-tt**3)/(A_*B_-tt**2)**2)==0)
chk('e(1) = E1', sp.cancel(et.subs(tt,1)-(A_+1)*(B_+1)/(A_*B_-1))==0)
chk('q h_b = bq(1-t^2/(ab))', sp.expand(Q_*(B_-tt**2/A_)-B_*Q_*(1-tt**2/(A_*B_)))==0)
# Lemma 1 Gram algebra: base Gram of (A,B,Q,C,D) and E,P
zz=sp.Symbol('zz'); cc=sp.Symbol('cc',positive=True); r=1-tt
Gb=sp.Matrix([[A_,tt,0,0,0],[tt,B_,zz,0,0],[0,zz,Q_,0,0],[0,0,0,cc,0],[0,0,0,0,1]])  # A,B,Q,C,D
Ev=sp.Matrix([tt*(B_+tt)/(A_*B_-tt**2),-tt*(A_+tt)/(A_*B_-tt**2),0,0,tt])
hb=B_-tt**2/A_
Vv=sp.Matrix([-tt/A_,1,-1,0,0]); Pv=(r-zz)/(hb+Q_-2*zz)*Vv+sp.Matrix([0,0,0,r/cc,0])
ip=lambda u,v:(u.T*Gb*v)[0]
basis=[sp.eye(5)[:,i] for i in range(5)]
prods=[sp.cancel(ip(Ev,x)) for x in basis]
chk('E.A=t, E.B=-t, E.D=t (as stated); also E.C=0', prods[0]==tt and prods[1]==-tt and prods[4]==tt and prods[3]==0)
print('   NOTE E.Q =',sp.factor(prods[2]),'(not stated in the document; its admissibility depends on the uncited ten Gram relations)')
chk('|E|^2 = e(t)', sp.cancel(ip(Ev,Ev)-et)==0)
chk('|V|^2 = h_b+q-2z', sp.cancel(ip(Vv,Vv)-(hb+Q_-2*zz))==0)
chk('|P|^2 = r^2/c + (r-z)^2/(h_b+q-2z)', sp.cancel(ip(Pv,Pv)-(r**2/cc+(r-zz)**2/(hb+Q_-2*zz)))==0)
print('   P.A,P.B,P.Q,P.C,P.D,P.E =',[sp.factor(sp.cancel(ip(Pv,x))) for x in basis+[Ev]])
chk('Gram(A,B,Q) PSD iff q h_b >= z^2 (det = a(q h_b - z^2))', sp.cancel(sp.Matrix([[A_,tt,0],[tt,B_,zz],[0,zz,Q_]]).det()-A_*(Q_*hb-zz**2))==0)
w=sp.Symbol('w'); qq=w**2*r**2/hb
chk('(10) at z=wr,q=w^2r^2/h_b: h_b+q-2z=(h_b-wr)^2/h_b and |P|^2 = r^2/c + h_b r^2(1-w)^2/(h_b-wr)^2', sp.cancel(hb+qq-2*w*r-(hb-w*r)**2/hb)==0 and sp.cancel((r**2/cc+(r-w*r)**2/(hb+qq-2*w*r))-(r**2/cc+hb*r**2*(1-w)**2/(hb-w*r)**2))==0)
# G6 sign tables (14) with flint
ctx=fmpz_mpoly_ctx.get(['B','X','E','Y','C','J']); Bv,X,E,Y,C,J=ctx.gens(); one=ctx.constant(1)
src_names=['a','b','e','q','c','j']
P6f=sp.Poly(G6,*[S[n] for n in src_names])
degs=[max(m[i] for m in P6f.monoms()) for i in range(6)]
def subst(nums,dens):
    out=ctx.constant(0)
    pn=[[n**i for i in range(m+1)] for n,m in zip(nums,degs)]; pd=[[dd**i for i in range(m+1)] for dd,m in zip(dens,degs)]
    for ex,co in zip(P6f.monoms(),P6f.coeffs()):
        v=ctx.constant(int(co))
        for i,k in enumerate(ex): v*=pn[i][k]*pd[i][degs[i]-k]
        out+=v
    return out
E1n=(one+X+Bv)*(Bv+one); E1d=Bv*X; E2n=(one+Y+Bv)*(Bv+one); E2d=Bv*Y
cases={('below','below'):([X,Bv,E,Y,C,J],[Bv*(one+X),one,one,Bv*(one+Y),one,one]),
       ('above','below'):([one+X,Bv,E1n*E,Y,C,J],[Bv,one,E1d*(one+E),Bv*(one+Y),one,one]),
       ('below','above'):([X,Bv,E,one+Y,E2n*C,J],[Bv*(one+X),one,one,Bv,E2d*(one+C),one]),
       ('above','above'):([one+X,Bv,E1n*E,one+Y,E2n*C,J],[Bv,one,E1d*(one+E),Bv,E2d*(one+C),one])}
for key,(nu,de) in cases.items():
    pol=subst(nu,de); cs=list(pol.to_dict().values())
    chk('(14) G6 table %s: %d terms, all coefficients >= 0'%(str(key),len(cs)), all(x>=0 for x in cs))



from flint import fmpz_mpoly_ctx

d=json.load(open(supp+'/proofs/manin_transfer_20261021/final_referee_T25_full_radial_20261022.json'))
terms=next(f['terms'] for f in d['radial_factors'] if f['degree']==9)
# vars q p j c e d b a z ; d=1 ; split by power of p ; source order (a,e,j,b,c,q)
src=fmpz_mpoly_ctx.get(['a','e','j','b','c','q'])
kk=[src.constant(0) for _ in range(3)]
for ex,co in terms:
    assert ex[8]==0 and ex[1]<=2
    kk[ex[1]]+=src.from_dict({(ex[7],ex[4],ex[2],ex[6],ex[3],ex[0]):int(co)})
degs=[max(k.degrees()[i] for k in kk) for i in range(6)]
ctx=fmpz_mpoly_ctx.get(['A','T','V','U','R','E']); A,T,V,U,R,E=ctx.gens(); one=ctx.constant(1)
def subst(poly,nums,dens):
    pn=[[n**i for i in range(m+1)] for n,m in zip(nums,degs)]; pd=[[x**i for i in range(m+1)] for x,m in zip(dens,degs)]
    out=ctx.constant(0)
    for ex,co in poly.to_dict().items():
        v=ctx.constant(int(co))
        for i,k in enumerate(ex): v*=pn[i][k]*pd[i][degs[i]-k]
        out+=v
    return out
Lnum=(A+one)*(one+T+A); Lden=A*T
charts={
 # name: (q,b,c,j) as (num,den), threshold (N,D) for P*=N/D
 'I-':([(A,one),(T,A*(one+T)),(V,one),(V+U,one)],(one,V+U)),
 'I+':([(A,one),(one+T,A),(Lnum*V,Lden*(one+V)),(Lnum*V+U*Lden*(one+V),Lden*(one+V))],(Lden*(one+V),Lnum*V+U*Lden*(one+V))),
}
for name,(qbcj,(N,D)) in charts.items():
    (qn,qd),(bn,bd),(cn,cd),(jn,jd)=qbcj
    for prod in ['ab<1','ab>1']:
        if prod=='ab<1': an,ad,en,ed=R*bd,bn*(one+R),E,one
        else: an,ad,en,ed=(one+R)*bd,bn,((one+R)*bd+bn)*(bn+bd)*E,bn*bd*R*(one+E)
        nums=[an,en,jn,bn,cn,qn]; dens=[ad,ed,jd,bd,cd,qd]
        sk=[subst(k,nums,dens) for k in kk]
        for lab,pol in [('middle',2*D*sk[0]+N*sk[1]),('endpoint',D*D*sk[0]+N*D*sk[1]+N*N*sk[2])]:
            cs=list(pol.to_dict().values()); neg=sum(1 for x in cs if x<0)
            ok&=(neg==0); print('T25 baseline %s %s %s: %d terms, %d negative -> %s'%(name,prod,lab,len(cs),neg,('PASS' if neg==0 else 'FAIL')))

import numpy as np

Pn=np.array(L.P,dtype=float)
F={k:np.array(L.F_VALUES[k],dtype=float) for k in range(82)}
def minphi(coefs,restarts=10,seed=0):
    vals=sum(c*F[k] for k,c in coefs.items())
    cv=np.linalg.lstsq(Pn,vals,rcond=None)[0]
    return L.min_phi_over_SO4(list(cv),restarts=restarts,seed=seed)[0]
rng=np.random.default_rng(11)
print('NUMERIC tests below: seed 11, BFGS 10 restarts per sample (lib min_phi_over_SO4); not proofs')
def run(label,gen,n=6):
    worst=np.inf
    for i in range(n):
        co=gen()
        if co is None: continue
        worst=min(worst,minphi(co,seed=i))
    print('NUMERIC %-70s worst min Phi = %.3e'%(label,worst))
U=lambda lo=0.05,hi=3: rng.uniform(lo,hi)
# ---- doc3: g = a f0 + e f1 + p f2 + j f5 + b f6 + c f39 + q f67 + d f70, d=1
def d3(a=0,e=0,p=0,j=0,b=0,c=0,q=0): return {0:a,1:e,2:p,5:j,6:b,39:c,67:q,70:1.0}
run('doc3 (4) SOC pj=1',lambda:(lambda p:d3(p=p,j=1/p))(U()))
run('doc3 (4) old-H c(p-1)=1',lambda:(lambda c:d3(c=c,p=1+1/c))(U()))
def f5a():
    a,b=U(0.5,3),U(0.5,3)
    if a*b<=1.05: return None
    return d3(a=a,b=b,e=(a+1)*(b+1)/(a*b-1))
run('doc3 (5) ab(e-1)=a+b+e+1',f5a)
def f5b():
    b,q=U(0.5,3),U(0.5,3)
    if b*q<=1.05: return None
    return d3(b=b,q=q,c=(q+1)*(b+1)/(b*q-1))
run('doc3 (5) bq(c-1)=b+q+c+1',f5b)
def f6():
    b,q=U(0.5,3),U(0.5,3)
    if b*q<=1.05: return None
    E2=(q+1)*(b+1)/(b*q-1); c=U(0.05,E2); return d3(b=b,q=q,c=c,p=1/c-1/E2)
run('doc3 (6) p=1/c-1/E2',f6)
def f7():
    b=U(1,3); w=U(0.05,0.95); q=w*w/b; c=U(); return d3(b=b,q=q,c=c,p=1/c+b*(1-w)**2/(b-w)**2)
run('doc3 (7) boundary',f7)
def lem9():
    t=U(0,0.95); a=U(0.2,3); b=U(max(t*t/a,0)+0.05,3)
    if a*b<=t*t+1e-3: return None
    hb=b-t*t/a; r=1-t; q=r*r/hb*U(1,2); c=U()
    e=t*t+t*t*(a+b+2*t)/(a*b-t*t); return d3(a=a,b=b,q=q,c=c,e=e,p=r*r/c)
run('doc3 Lemma1 (9) boundary',lem9,8)
def lem10():
    t=U(0,0.95); a=U(0.2,3); b=U(0.2,3)
    if a*b<=t*t+1e-3: return None
    hb=b-t*t/a; r=1-t; w=U(0,1)
    if hb<=w*r: return None
    q=w*w*r*r/hb; c=U(); e=t*t+t*t*(a+b+2*t)/(a*b-t*t)
    return d3(a=a,b=b,q=q,c=c,e=e,p=r*r/c+hb*r*r*(1-w)**2/(hb-w*r)**2)
run('doc3 Lemma1 (10) boundary',lem10,8)
def lem9c():
    co=lem9()
    if co: co[2]*=0.6
    return co
run('doc3 Lemma1 (9) control p=0.6*r^2/c',lem9c,4)


import itertools
kk_I=None
d=json.load(open(supp+'/proofs/manin_transfer_20261021/final_referee_T25_full_radial_20261022.json'))
terms=next(f['terms'] for f in d['radial_factors'] if f['degree']==9)
src=fmpz_mpoly_ctx.get(['a','e','j','b','c','q'])
kk=[src.constant(0) for _ in range(3)]
for ex,co in terms:
    kk[ex[1]]+=src.from_dict({(ex[7],ex[4],ex[2],ex[6],ex[3],ex[0]):int(co)})
degs=[max(k.degrees()[i] for k in kk) for i in range(6)]
NV=sp.symbols('v0:6'); ctx=fmpz_mpoly_ctx.get([str(v) for v in NV])
def fl(expr):
    P=sp.Poly(sp.expand(expr),*NV); return ctx.from_dict({m:int(c) for m,c in P.as_dict().items()})
def nd(expr):
    n,dd=sp.fraction(sp.cancel(sp.together(expr)))
    # make denominator have positive leading behaviour: keep as is, sign checked by caller
    return fl(n),fl(dd)
def subst(poly,nums,dens):
    pn=[[n**i for i in range(m+1)] for n,m in zip(nums,degs)]; pd=[[x**i for i in range(m+1)] for x,m in zip(dens,degs)]
    out=ctx.constant(0)
    for ex,co in poly.to_dict().items():
        v=ctx.constant(int(co))
        for i,k in enumerate(ex): v*=pn[i][k]*pd[i][degs[i]-k]
        out+=v
    return out
def negs(p): return [(m,c) for m,c in p.to_dict().items() if c<0]
def table(label,subs,thr,allow_squares=False):
    """subs: sympy exprs for (a,e,j,b,c,q); thr: sympy threshold. All denominators must be positive polynomials."""
    nums,dens=[],[]
    for ex in subs:
        n,dd=nd(ex)
        if any(c<0 for c in dd.to_dict().values()):  # flip sign so denominator has nonneg coefficients
            n,dd=-n,-dd
        assert all(c>=0 for c in dd.to_dict().values()), (label,'denominator sign')
        nums.append(n); dens.append(dd)
    N,D=nd(thr)
    if any(c<0 for c in D.to_dict().values()): N,D=-N,-D
    assert all(c>=0 for c in D.to_dict().values()) and all(c>=0 for c in N.to_dict().values())
    sk=[subst(k,nums,dens) for k in kk]
    res=[]
    for lab,pol in [('middle',2*D*sk[0]+N*sk[1]),('endpoint',D*D*sk[0]+N*D*sk[1]+N*N*sk[2])]:
        ng=negs(pol); ok=not ng
        extra=''
        if ng and allow_squares:
            # auditor's own AM-GM repair: -k X^m <= (k/2)(X^{m+dl}+X^{m-dl}) using available coefficients
            cf={mm_:sp.Integer(int(v_)) for mm_,v_ in pol.to_dict().items()}; repaired=True
            for m,c in ng:
                k=sp.Rational(-int(c)); found=False
                for dl in itertools.product([-2,-1,0,1,2],repeat=6):
                    if not any(dl): continue
                    mp=tuple(a+b for a,b in zip(m,dl)); mm=tuple(a-b for a,b in zip(m,dl))
                    if min(mm)<0 or min(mp)<0: continue
                    if cf.get(mp,0)>=k/2 and cf.get(mm,0)>=k/2:
                        cf[mp]-=sp.Rational(k,2); cf[mm]-=sp.Rational(k,2); cf[m]=0; found=True; break
                repaired&=found
            ok=repaired and all(v>=0 for v in cf.values()); extra=' ; %d negative coefficients, auditor AM-GM repair %s'%(len(ng),'succeeds' if ok else 'FAILS')
        print('T25 %s %s: %d terms, %d negative -> %s%s'%(label,lab,len(pol),len(ng),'PASS' if ok else 'FAIL',extra))
        res.append(ok)
    return all(res)
A,T,V,U,R,E=NV
allok=True
def ae(bexpr,prod):
    if prod=='ab<1': return R/(bexpr*(1+R)), E
    a=(1+R)/bexpr; E1=(a+1)*(bexpr+1)/(a*bexpr-1); return a, E1*E/(1+E)
L_=lambda x,y:(x+1)*(y+1)/(x*y-1)
base={
 'II':  (lambda: (A,(1+T)/A,L_(A,(1+T)/A)*V/(1+V),L_(A,(1+T)/A)*V*U/((1+V)*(1+U))), lambda x,y,z,s: 1/z-1/L_(x,y)),
 'IIIa-':(lambda: (T*V/((1+T)*(1+V)),T/(1+T),A,A*U/((1+A)*(1+U))), lambda x,y,z,s: 1+1/z),
 'IIIa+':(lambda: (T*V/((1+T)*(1+V)),T/(1+T),A,A/(1+A)+A**2*U/((1+A)*(1+U))), lambda x,y,z,s: 1/s),
 'IIIb+':(lambda: (V/((1+T)*(1+V)),1+T,A,A*(1+T)/(A+1+T)+A**2*U/((A+1+T)*(1+U))), lambda x,y,z,s: 1/s),
 'IIIb-':(lambda: ((V/(1+V))**2/(1+T),1+T,A,A*(1+T)*U/((A+1+T)*(1+U))), lambda x,y,z,s: 1/z+y*(1-V/(1+V))**2/(y-V/(1+V))**2),
}
jobs=[('II',1,'ab<1')]+[(c,o,p) for c in ['IIIa-','IIIa+','IIIb+'] for o in (1,2) for p in ('ab<1','ab>1')]+[('IIIb-',1,'ab<1'),('IIIb-',2,'ab<1'),('IIIb-',2,'ab>1')]
for ch,orient,prod in jobs:
    x,y,z,s=base[ch][0](); thr=base[ch][1](x,y,z,s)
    if orient==1: q,b=x,y
    else: b,q=x,y
    a,e=ae(b,prod)
    allok&=table('%s orient%d %s'%(ch,orient,prod),(a,e,s,b,z,q),thr)
# missing region II: b=B(=A), ab=1+X(=T), bq=1+Y(=V), t=Z/(1+Z) (Z=U), c=E2 Vc/(1+Vc), j=cU'/(1+U') (U'=E), split by W (=R)
B,X,Y,Z,W,Uu=A,T,V,U,R,E
t=Z/(1+Z); a=(1+X)/B; b=B; q=(1+Y)/B
et=t**2+t**2*(a+b+2*t)/(a*b-t**2); E2=(q+1)*(b+1)/(b*q-1)
for lab,Vc,thrf in [('regionII P1',Z*(Z+2)*W/(1+W),lambda c:(1-t)**2/c),('regionII P0',Z*(Z+2)+W,lambda c:1/c-1/E2)]:
    c=E2*Vc/(1+Vc); j=c*Uu/(1+Uu)
    allok&=table(lab,(a,et,j,b,c,q),thrf(c),allow_squares=True)
# missing region IIIb-: b=1+B, ab=1+X, t=Z/(1+Z), c=C(=V), j=CbU/((C+b)(1+U)), q by (20)/(21) with W
Bb,X,Z,C,Uu,W=A,T,U,V,E,R
b=1+Bb; a=(1+X)/b; t=Z/(1+Z); et=t**2+t**2*(a+b+2*t)/(a*b-t**2); Dd=(1+X)*(1+Z)**2-Z**2
hb=b*Dd/((1+X)*(1+Z)**2); r=1/(1+Z); j=C*b*Uu/((C+b)*(1+Uu))
q20=(1+X+Dd*W)/(b*Dd*(1+W)); allok&=table('IIIb- missing (20)',(a,et,j,b,C,q20),r**2/C)
w=W/(1+W); q21=W**2*(1+X)/((1+W)**2*b*Dd); allok&=table('IIIb- missing (21)',(a,et,j,b,C,q21),r**2/C+hb*r**2*(1-w)**2/(hb-w*r)**2)
ok&=allok; print('remaining 32 baseline tables and 8 new-region tables:', 'PASS' if allok else 'FAIL')

print('ALL EXACT CHECKS PASS' if ok else 'SOME EXACT CHECK FAILED'); sys.exit(0 if ok else 1)
