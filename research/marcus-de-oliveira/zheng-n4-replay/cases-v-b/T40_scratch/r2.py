import sys, time, itertools
sys.path.insert(0,'C:/Users/nrvig/Projects/math-commons/research/marcus-de-oliveira/zheng-n4-replay/cases-v-b')
import replay_T38_common as c
from flint import fmpz_mpoly_ctx, Ordering, fmpz
names=['q','p','j','c','e','d','b','a','z']
ctx=fmpz_mpoly_ctx.get(names, Ordering.lex)
q,p,j,cc,e,d,b,a,z=ctx.gens()
src={5:q,1:p,12:j,9:cc,2:e,35:d,33:b,29:a}
u=[ctx.constant(0)]*14
for f,v in src.items():
    for k in range(14):
        x=c.F_COORD[f][k]
        if x!=0: u[k]=u[k]+int(x)*v
u[13]=u[13]+z
C=c.cubic_C(u)
nt=lambda P: len(P.to_dict())
def coef_z(P,k):
    D=P.to_dict(); out=ctx.constant(0)
    zi=8
    for ex,co in D.items():
        if ex[zi]==k:
            ex2=list(ex); ex2[zi]=0
            out=out+ctx.from_dict({tuple(ex2):co})
    return out
C0=coef_z(C,0); Q=coef_z(C,1)
print('degz C',max(ex[8] for ex in C.to_dict()), nt(C0), nt(Q))
A_=1+a+b; B_=1+q+e; N=(q+1)*(e+1)*A_+B_*a*b
# doc at d=1 : C0 = -(AB c + N - j e a)
d1={'d':1}
print(C0.subs(d1) if False else '')
def sub1(P): return P.subs({'d':1}) 
one=ctx.constant(1)
def d1(P): return P.subs({'d':1})
Ad=1+a+b; Bd=1+q+e; Nd=(q+1)*(e+1)*Ad+Bd*a*b; Nj=Nd-j*e*a
print('C0 = -(ABc+Nj)?', d1(C0)== -(Ad*Bd*cc+Nj))
Qdoc_low = -( (q+e+a+b+(1-p))*cc + q*(e+a+1)+b*(e+a+1)+(1-j)*(e+a)+1+p*j )
print('Q = low-pj formula?', d1(Q)==Qdoc_low)
Qh=-((q+e+1+a+b-p)*cc+(q+b+1)*(e+a+1)-j*(e+a)+p*j)
print('Q = high-p formula?', d1(Q)==Qh, d1(Q)==Qdoc_low)
S=a+e;J=j-1
Qj=-((q+b+1+S-p)*cc+1-J*S+(q+b)*(1+S)+p*j)
print('Q = high-j formula (-Q=..)?', d1(Q)==Qj)
print(Q)
