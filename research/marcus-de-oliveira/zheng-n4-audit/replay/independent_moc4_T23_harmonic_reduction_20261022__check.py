#!/usr/bin/env python3
"""Own-code replay for document independent_moc4_T23_harmonic_reduction_20261022 (Zheng n=4 audit). Run: python -I -B this.py
Reads the supplement location from MOC4_SUPP (only JSON data files are loaded, as data).
Prints one line per check; exits nonzero if an EXACT check fails."""
import os, sys
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE,'..','lib'))
SUPP=os.environ.get('MOC4_SUPP') or (sys.argv[1] if len(sys.argv)>1 else None)
import moc4lib as L, sympy as sp, functools, numpy as np
print=functools.partial(print,flush=True)
OK=True
def chk(name,cond):
    global OK; OK&=bool(cond); print(('PASS ' if cond else 'FAIL ')+name)
T={2,3,4,6,7,8}
gens=L.generators(T)
chk('generators(T) = [0,1,2,3,4,5,6,53,54,67] (document lists 8 of them; f1,f4 eliminated)', gens==[0,1,2,3,4,5,6,53,54,67])
chk('T is paper T_24 up to G', L.find_g(L.CONTACT_SETS[24],frozenset(T)) is not None)
def dec(v):
    d=v[5]-v[16]; out=[0]*24
    for coef,k in [(v[22],0),(v[13],2),(v[12],3),(v[10],5),(v[17],6),(v[16],53),(d,54),(v[1],67)]:
        for i in range(24): out[i]+=coef*L.F_VALUES[k][i]
    return out
for k in gens:
    v=list(L.F_VALUES[k])
    chk('decomposition reproduces f%d and g5-g16=g9-g17 (delta=%d)'%(k,v[5]-v[16]), dec(v)==v and v[5]-v[16]==v[9]-v[17])
chk('rank of the 10 generators = 8', sp.Matrix([L.F_VALUES[k] for k in gens]).rank()==8)
V=lambda k: sp.Matrix(L.F_VALUES[k])
chk('f54+f1 = f0+f6+f53 (so delta<=0 gives a minor decomposition)', V(54)+V(1)==V(0)+V(6)+V(53))
cf4=[L.F_VALUES[4][i] for i in (22,13,12,10,17,16,1)]
chk('f54+f4 = nonnegative minor combination (coefficients %s on f0,f2,f3,f5,f6,f53,f67)'%cf4, all(x>=0 for x in cf4) and sp.Matrix(dec(list(L.F_VALUES[4])))+V(54)==V(4)+V(54))
x,s,t,r,u=sp.symbols('x s t r u'); X=1-x; q=r-u
d=1-t+u*X; l=t*x+s*X; h=1-u*X/t; n=t+q*X
z=1-r*X/d; w=x*n/l; c=w/h; v=1-s*X*n/l
a=[z,0,1,x]; cc=[c,v,w,None]
def Fp(aa,c2,j):
    pr=1;p=L.PERMS[j]
    for i in range(4):
        if c2[p[i]] is None: continue
        pr*=(aa[i]-c2[p[i]])
    return pr
F=[Fp(a,cc,j) for j in range(24)]
norm=(z-c)*w*X
H=[sp.cancel(sp.together((F[j]-F[2])/norm)) for j in range(24)]
for j,val in {2:0,3:1,4:s,6:r-u,7:r+t-u,8:r,9:r+t-u*X}.items(): chk('H%d = %s'%(j,val), sp.cancel(H[j]-val)==0)
J=q*s*u*X**2+(-q*s*t+q*t+s*t*u-s*u+t*u)*X-t*(r+s*(t-1))
chk('z-c = XJ/((t-uX) l d)', sp.cancel(z-c - X*J/((t-u*X)*l*d))==0)
chk('z-v = t h (z-c)', sp.cancel(z-v - t*h*(z-c))==0)
A1=s-q-t
chk('H1 = s + A1/(q(x-1)-t)', sp.cancel(H[1]-(s+A1/(q*(x-1)-t)))==0)
k10=sp.cancel(H[10]+t*r/q*H[1]); chk('H10 + (tr/q) H1 is a real constant (=> g10 = -(tr/q) g1): %s'%k10, sp.diff(k10,x)==0)
n5,d5=sp.fraction(H[5]); chk('H5 -> 1 as x -> oo', sp.Poly(n5,x).degree()==sp.Poly(d5,x).degree() and sp.cancel(sp.Poly(n5,x).LC()/sp.Poly(d5,x).LC())==1)
lin=sp.expand(q*(x-1)-t); Jx=sp.expand(J)
allowed=True
for j in range(24):
    for fac,m in sp.factor_list(sp.fraction(H[j])[1])[1]:
        if sp.Poly(fac,x).degree()==0: continue
        if sp.cancel(fac/x).is_number or sp.cancel(fac/lin).is_number or sp.cancel(fac/Jx).is_number: continue
        allowed=False; print('   unexpected denominator factor in H%d: %s'%(j,fac))
chk('every denominator factor of every H_j (j=0..23) is x, q(x-1)-t or J (so nonreal finite poles are zeros of J)', allowed)
for j in [13,16]:
    dn=sp.fraction(H[j])[1]
    chk('J divides the denominator of H%d exactly once'%j, sp.Poly(dn,x).rem(sp.Poly(Jx,x)).is_zero and not sp.Poly(sp.cancel(dn/Jx),x).rem(sp.Poly(Jx,x)).is_zero)
n13,d13=sp.fraction(H[13]); n16,d16=sp.fraction(H[16])
ex=sp.together(sp.cancel((n13*d16)/(n16*d13))*(z-1)/z-1)
chk('res H13/res H16 = z/(z-1) at zeros of J (difference of ratios vanishes mod J)', sp.Poly(sp.numer(ex),x).rem(sp.Poly(Jx,x)).is_zero)
zz,xx,cv,vv,ww=sp.symbols('zz xx cc vv ww'); aa=[zz,0,1,xx]; c2=[cv,vv,ww,None]
Fs=lambda j: Fp(aa,c2,j)
chk('raw F16-F2 at c=v=z equals (w-z)(x-z)(z-1)', sp.expand((Fs(16)-Fs(2)).subs({cv:zz,vv:zz})-(ww-zz)*(xx-zz)*(zz-1))==0)
chk('raw F13-F2 at c=v=z equals z(w-z)(x-z)', sp.expand((Fs(13)-Fs(2)).subs({cv:zz,vv:zz})-zz*(ww-zz)*(xx-zz))==0)
chk('numerator of z-x is (x-1)(uX-(r+t-1))*(-1) i.e. x=z iff X=0 or d=r', sp.expand(sp.numer(sp.together(z-x))-(x-1)*(r+t+u*x-u-1))==0)
chk('z = (1-t-qX)/(1-t+uX), Mobius determinant -(1-t)r', sp.cancel(z-(1-t-q*X)/(1-t+u*X))==0 and sp.expand((-q)*(1-t)-(1-t)*u+(1-t)*r)==0)
D3=Fs(3)-Fs(2)
chk('F3-F2 = w(z-c)X (normalization)', sp.expand(D3-ww*(zz-cv)*(1-xx))==0)
S_=sp.cancel((Fs(4)-Fs(2))/D3); R_=sp.cancel((Fs(8)-Fs(2))/D3); U_=sp.cancel(R_-(Fs(6)-Fs(2))/D3); T_=sp.cancel((Fs(7)-Fs(6))/D3)
print('   ratios in point coordinates: s=%s, r=%s, u=%s, t=%s'%(sp.factor(S_),sp.factor(R_),sp.factor(U_),sp.factor(T_)))
sub={s:S_,t:T_,r:R_,u:U_,x:xx}
for name,e_,tg in [('z',z,zz),('w',w,ww),('c',c,cv),('v',v,vv)]: chk('chart is exhaustive: ratios of a configuration give back %s'%name, sp.cancel(e_.subs(sub)-tg)==0)
for name,e_ in [('d',d),('l',l),('h',h),('n',n),('t-uX',t-u*X)]: print('   %s in point coordinates = %s'%(name,sp.factor(sp.cancel(e_.subs(sub)))))
# numeric: SOC citation for (f5,f2,f54)
Pn=np.array(L.P,dtype=float); Fv={k:np.array(L.F_VALUES[k],dtype=float) for k in range(82)}
def minphi(co,seed):
    vals=sum(cf*Fv[k] for k,cf in co.items()); cvv=np.linalg.lstsq(Pn,vals,rcond=None)[0]
    return L.min_phi_over_SO4(list(cvv),restarts=10,seed=seed)[0]
rng=np.random.default_rng(11); worst=np.inf
for i in range(6):
    A_,B_=rng.uniform(0.05,3,2); worst=min(worst,minphi({5:A_,2:B_,54:np.sqrt(A_*B_)},i))
print('NUMERIC SOC boundary delta^2=g10*g13 (cost g10 f5+g13 f2+delta f54), 6 samples seed 11, 10 BFGS restarts: worst min Phi = %.3e'%worst)
A_,B_=1.3,0.7; print('NUMERIC control delta=3*sqrt(g10 g13): min Phi = %.3e (negative expected)'%minphi({5:A_,2:B_,54:3*np.sqrt(A_*B_)},0))
print('ALL EXACT CHECKS PASS' if OK else 'SOME EXACT CHECK FAILED'); sys.exit(0 if OK else 1)
