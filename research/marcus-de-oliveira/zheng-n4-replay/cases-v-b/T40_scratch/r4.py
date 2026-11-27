from dets import *
import pickle
from flint import fmpz_mpoly_ctx, Ordering
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
D=ctx.from_dict(pickle.load(open('D.pkl','rb')))
def coef_z(P,k):
    out={}
    for ex,co in P.to_dict().items():
        if ex[8]==k:
            e2=list(ex); e2[8]=0; out[tuple(e2)]=co
    return ctx.from_dict(out)
C=c.cubic_C(u); C0=coef_z(C,0); Q=coef_z(C,1)
Dk=[coef_z(D,k) for k in range(8)]
print([len(x.to_dict()) for x in Dk])
N=ctx.constant(0)
for k in range(5):
    N=N+Dk[k]*(-C0)**k*Q**(4-k)
print('N terms',len(N.to_dict()))
q2,r=divmod(N,Q*Q)
print('rem zero',r.is_zero(), len(q2.to_dict()))
q3,r=divmod(q2,ctx.constant(-4)); print(r.is_zero())
G=q3.sqrt()
print('G terms',len(G.to_dict()), G*G==q3)
pickle.dump((G.to_dict(),C0.to_dict(),Q.to_dict()),open('G.pkl','wb'))
# G at d=1 terms; alpha beta gamma
def collect_c(G):
    out={}
    for ex,co in G.to_dict().items():
        k=ex[3]; e2=list(ex); e2[3]=0
        out.setdefault(k,{})[tuple(e2)]=co
    return {k:ctx.from_dict(v) for k,v in out.items()}
cg=collect_c(G)
print({k:len(v.to_dict()) for k,v in cg.items()})
