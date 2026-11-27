import pickle, time
from flint import fmpz_mpoly_ctx, Ordering
names=['q','p','j','c','e','d','b','a','z']
ctx=fmpz_mpoly_ctx.get(names, Ordering.lex)
G,C0,Q=[ctx.from_dict(x) for x in pickle.load(open('G.pkl','rb'))]
def coef_z(P,k):
    out={}
    for ex,co in P.to_dict().items():
        if ex[8]==k:
            e2=list(ex); e2[8]=0; out[tuple(e2)]=co
    return ctx.from_dict(out)
for k in (0,1):
    Bk=ctx.from_dict(pickle.load(open('Bord%d.pkl'%k,'rb')))
    B0,B1,B2=[coef_z(Bk,i) for i in range(3)]
    num=B2*C0*C0-B1*C0*Q+B0*Q*Q
    qq,rr=divmod(num,G)
    print(k,'numer terms',len(num.to_dict()),'rem zero',rr.is_zero(),'quot terms',len(qq.to_dict()))
