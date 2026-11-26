#!/usr/bin/env python3
"""Own-code replay for document root_moc4_T24_uniform_20261022 (Zheng n=4 audit). Run: python -I -B this.py
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
supp=SUPP
T={3,4,7,8,10,12}
ok=True
def chk(name,cond):
    global ok; ok&=bool(cond); print(('PASS ' if cond else 'FAIL ')+name)
chk('generators(T24-doc) = {0,1,2,6,53,54,55,59,67,70}', L.generators(T)==[0,1,2,6,53,54,55,59,67,70])
V=lambda k: sp.Matrix(L.F_VALUES[k])
chk('f54+f70 = f0+f6+f55+f67', V(54)+V(70)==V(0)+V(6)+V(55)+V(67))
chk('f54+f1 = f0+f6+f53', V(54)+V(1)==V(0)+V(6)+V(53))
order=[0,2,59,6,53,54,55,67]
chk('8 canonical profiles rank 8', sp.Matrix([L.F_VALUES[k] for k in order]).rank()==8)
s1=find_sym(T,[(54,70),(70,54),(2,2),(59,59)])
chk('exists g in G stabilizing T, f54<->f70, fixing f2,f59 (count %d)'%len(s1), len(s1)>0)
print('   e.g.',s1[:2], ' images of all generators under first:', {k:gimage(s1[0],k) for k in L.generators(T)})
s2=find_sym(T,[(67,55),(55,67),(0,6),(6,0),(2,2),(59,59),(53,53),(54,54)])
chk('exists g stabilizing T realizing (a,b,q,c)->(b,a,c,q) with p,j,e,d fixed (count %d)'%len(s2), len(s2)>0)
q,p,j,c,e,d,b,a,h=sp.symbols('q p j c e d b a z')
vals=[q*L.F_VALUES[0][i]+p*L.F_VALUES[2][i]+j*L.F_VALUES[59][i]+c*L.F_VALUES[6][i]+e*L.F_VALUES[53][i]+d*L.F_VALUES[54][i]+b*L.F_VALUES[55][i]+a*L.F_VALUES[67][i] for i in range(24)]
cv=cvec_sym(vals); chk('symbolic preimage c with Pc = g', check_inP(vals,cv))
y=[cv[i]+(h if i==13 else 0) for i in range(14)]
Cy=sp.expand(CUB.subs(dict(zip(MS,y))))
I_=c+d+q; U=(a+d)*(b+d); Vv=U+d*e; S=a+b+d+e; al=I_*U+c*q*(a+b+d); Cst=al+e*(c+d)*(d+q); K=a*b*(e-d)-d**2*(a+b+e+d)
Q=e*j-d*(d+e)-p*j
lam=sp.cancel(Cy/(-Cst+h*Q)); chk('(2) C(g+h kappa) = lambda*(-C*+hQ), lambda = %s'%lam, lam.is_number and lam!=0)
# load G from JSON (data only)
dj,var,polys=load_poly(supp+'/proofs/manin_transfer_20261021/root_moc4_T24_full_radial_20261022.json',None)
G=polys[9]; print('   JSON G: terms',len(G.terms()),'total degree',G.total_degree(), 'vars',var)
Gx=G.as_expr()
chk('JSON G has 2095 terms', len(G.terms())==2095)
# Hessian identity at random rational points
rnd=random.Random(2024); ratios=set()
for trial in range(6):
    pt={q:sp.Rational(rnd.randint(0,9),rnd.randint(1,5)),p:sp.Rational(rnd.randint(0,9),rnd.randint(1,5)),j:sp.Rational(rnd.randint(0,9),rnd.randint(1,5)),
        c:sp.Rational(rnd.randint(0,9),rnd.randint(1,5)),e:sp.Rational(rnd.randint(0,9),rnd.randint(1,5)),d:sp.Rational(rnd.randint(1,9),rnd.randint(1,5)),
        b:sp.Rational(rnd.randint(0,9),rnd.randint(1,5)),a:sp.Rational(rnd.randint(0,9),rnd.randint(1,5))}
    Qv=Q.subs(pt)
    if Qv==0: continue
    hv=Cst.subs(pt)/Qv; pt2=dict(pt); pt2[h]=hv
    yv=[sp.nsimplify(yy.subs(pt2)) for yy in y]
    Hm=HESS.subs(dict(zip(MS,yv)))
    dt=exact_det(Hm)
    Gv=Gx.subs({sp.Symbol('q'):pt[q],sp.Symbol('p'):pt[p],sp.Symbol('j'):pt[j],sp.Symbol('c'):pt[c],sp.Symbol('e'):pt[e],sp.Symbol('d'):pt[d],sp.Symbol('b'):pt[b],sp.Symbol('a'):pt[a],sp.Symbol('z'):0})
    rhs=-4*Qv**2*Gv**2
    rr=sp.Rational(dt.numerator,dt.denominator)*Qv**4/rhs if rhs!=0 else ('rhs0',dt)
    ratios.add(rr); print('   point',trial,'det*Q^4 / (-4Q^2G^2) =',rr)
chk('(3) det Hess(g+h kappa)*Q^4 = const*(-4 Q^2 G^2) at %d random rational points (const %s)'%(len(ratios) and 6,ratios), len(ratios)==1)
# k0 identity (7)
Gp=sp.Poly(Gx.subs(sp.Symbol('z'),0),sp.Symbol('j'))
syms={n:sp.Symbol(n) for n in 'qpjcedba'}
k=[Gp.coeff_monomial(sp.Symbol('j')**i) for i in range(3)]
print('   G degree in j:',Gp.degree())
sub={q:syms['q'],p:syms['p'],c:syms['c'],e:syms['e'],d:syms['d'],b:syms['b'],a:syms['a']}
beta=-K+(p+c+q)*Vv+c*q*S; Wc=p*(b+d); Wq=p*(a+d)
Hh=Vv*(c*(c+d)*Wc+q*(q+d)*Wq)+c*q*(c*q*p*S+c*(S*Wc+p*Vv)+q*(S*Wq+p*Vv)+p*(a+b+2*d)*Vv+a*b*(d+e)*p)
rhs7=sp.expand((Cst*(al*beta+c*q*d**2*p*S+Hh)).subs(sub))
diff7=sp.expand(k[0]-rhs7)
chk('(7) k0 = C*{alpha beta + cq d^2 p S + H} exactly (as polynomials)', diff7==0)
if diff7!=0:
    r7=sp.cancel(k[0]/rhs7); print('   ratio k0/rhs7 =', r7)
# p=e=0: all coefficients positive
G00=sp.Poly(Gx.subs({sp.Symbol('p'):0,sp.Symbol('e'):0,sp.Symbol('z'):0}),*[syms[n] for n in 'qjcdba'])
chk('at p=e=0 every coefficient of G is positive (%d terms), contains d^9: %s'%(len(G00.terms()), G00.coeff_monomial(syms['d']**9)), all(co>0 for co in G00.coeffs()) and G00.coeff_monomial(syms['d']**9)>0)
# invariance of G under (a,b,q,c)->(b,a,c,q)
Gs=sp.expand(Gx.subs({syms['a']:syms['b'],syms['b']:syms['a'],syms['q']:syms['c'],syms['c']:syms['q']},simultaneous=True))
chk('G invariant under (a,b,q,c)->(b,a,c,q)', sp.expand(Gs-Gx)==0)
sp_k=[sp.expand(kk) for kk in k]



from flint import fmpz_mpoly_ctx

which='Im,IIIam,IIIap,Ip,IIIbp,II,IIIbm'.split(',')
d=json.load(open(supp+'/proofs/manin_transfer_20261021/root_moc4_T24_full_radial_20261022.json'))
terms=next(f['terms'] for f in d['radial_factors'] if f['degree']==9)
# variables q p j c e d b a z ; d=1 ; split by power of j
src=fmpz_mpoly_ctx.get(['a','b','e','p','q','c'])
kk=[src.constant(0) for _ in range(3)]
for ex,co in terms:
    assert ex[8]==0
    kk[ex[2]]+=src.from_dict({(ex[7],ex[6],ex[4],ex[1],ex[0],ex[3]):int(co)})
ctx=fmpz_mpoly_ctx.get(['A','T','V','U','Q','C'])
A,T,V,U,Q,C=ctx.gens(); one=ctx.constant(1)
def chart(name):
    # returns (nums, dens) for (a,b,e,p,q,c) and (N,D) of threshold
    if name=='Im':
        return [A,T,V,V+U,Q,C],[one,A*(one+T),one,one,one,one],(one,V+U)
    Lnum=(A+one)*(one+T+A); Lden=A*T   # L for a=A, b=(1+T)/A
    if name=='Ip':
        # e = L V/(1+V); p = e + U
        return [A,one+T,Lnum*V,Lnum*V+U*Lden*(one+V),Q,C],[one,A,Lden*(one+V),Lden*(one+V),one,one],(Lden*(one+V),Lnum*V+U*Lden*(one+V))
    if name=='II':
        # p = e U/(1+U); J* = 1/e - 1/L = 1/(L V)
        return [A,one+T,Lnum*V,Lnum*V*U,Q,C],[one,A,Lden*(one+V),Lden*(one+V)*(one+U),one,one],(Lden,Lnum*V)
    if name=='IIIam':
        # b=T/(1+T), a=bV/(1+V), e=A, p=AU/((1+A)(1+U)), J*=1+1/A
        return [T*V,T,A,A*U,Q,C],[(one+T)*(one+V),one+T,one,(one+A)*(one+U),one,one],(A+one,A)
    if name=='IIIap':
        # p = A/(1+A) + A^2U/((1+A)(1+U)) = A(1+U+AU)/((1+A)(1+U)); J*=1/p
        pn=A*(one+U+A*U); pd=(one+A)*(one+U)
        return [T*V,T,A,pn,Q,C],[(one+T)*(one+V),one+T,one,pd,one,one],(pd,pn)
    if name=='IIIbp':
        b=one+T; pn=A*b*(one+U)+A*A*U; pd=(A+b)*(one+U)
        return [V,b,A,pn,Q,C],[b*(one+V),one,one,pd,one,one],(pd,pn)
    if name=='IIIbm':
        b=one+T; F=(one+T)*(one+V)-V
        return [V*V,b,A,A*b*U,Q,C],[(one+V)**2*b,one,one,(A+b)*(one+U),one,one],(F*F+A*b,A*F*F)
def subst(poly,nums,dens,degs):
    pw_n=[[n**i for i in range(m+1)] for n,m in zip(nums,degs)]
    pw_d=[[dd**i for i in range(m+1)] for dd,m in zip(dens,degs)]
    out=ctx.constant(0)
    for ex,co in poly.to_dict().items():
        v=ctx.constant(int(co))
        for i,k in enumerate(ex):
            v*=pw_n[i][k]*pw_d[i][degs[i]-k]
        out+=v
    return out
degs=[max(k.degrees()[i] for k in kk) for i in range(6)]
for name in which:
    t0=time.time()
    nums,dens,(N,D)=chart(name)
    sk=[subst(k,nums,dens,degs) for k in kk]
    mid=2*D*sk[0]+N*sk[1]
    end=D*D*sk[0]+N*D*sk[1]+N*N*sk[2]
    for lab,pol in [('middle',mid),('endpoint',end)]:
        cs=list(pol.to_dict().values()); neg=sum(1 for x in cs if x<0)
        ok&=(neg==0); print('%s %s: %d terms, %d negative coefficients -> %s  (%.1fs)'%(name,lab,len(cs),neg,('PASS' if neg==0 else 'FAIL'),time.time()-t0))
    k0=sk[0]; print('%s k0 numerator: %d terms, negative %d'%(name,len(k0),sum(1 for x in k0.to_dict().values() if x<0)))

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
# ---- doc2: canonical g = q f0 + p f2 + j f59 + c f6 + e f53 + d f54 + b f55 + a f67, d=1
def d2(q=0,p=0,j=0,c=0,e=0,b=0,a=0): return {0:q,2:p,59:j,6:c,53:e,54:1.0,55:b,67:a}
def soc():
    p=U(); return d2(p=p,j=1/p)
run('doc2 (4) SOC pj=1',soc)
def oldH():
    e=U(); return d2(e=e,j=1+1/e)
run('doc2 (4) old-H e(j-1)=1',oldH)
def four():
    a,b=U(0.5,3),U(0.5,3)
    if a*b<=1.05: return None
    e=(a+1)*(b+1)/(a*b-1); return d2(a=a,b=b,e=e)
run('doc2 (4) K=0 (e=L, j=0)',four)
def cap5():
    a,b=U(0.5,3),U(0.5,3)
    if a*b<=1.05: return None
    Lr=(a+1)*(b+1)/(a*b-1); e=U(0.05,Lr); j=1/e-1/Lr; return d2(a=a,b=b,e=e,j=j)
run('doc2 Lemma1 (5) boundary j=1/e-1/L',cap5)
def cap5c():
    co=cap5()
    if co: co[59]*=0.6
    return co
run('doc2 Lemma1 (5) control j=0.6*(threshold)',cap5c,4)
def cap6(k=1.0):
    b=U(1,3); w=U(0.05,0.95); a=w*w/b; e=U(); j=1/e+b*(1-w)**2/(b-w)**2; return d2(a=a,b=b,e=e,j=k*j)
run('doc2 Lemma1 (6) boundary',cap6)
run('doc2 Lemma1 (6) control 0.6*threshold',lambda:cap6(0.6),4)
def dsoc():
    a=U(0.3,3); b=1/a; e=U(); return d2(a=a,b=b,e=e,j=1/e)
run('doc2 double-SOC ab=1, ej=1',dsoc)

print('ALL EXACT CHECKS PASS' if ok else 'SOME EXACT CHECK FAILED'); sys.exit(0 if ok else 1)
