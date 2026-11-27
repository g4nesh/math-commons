import pickle, time
from flint import fmpz_mpoly_ctx, Ordering
names=['q','p','j','c','e','d','b','a','z']
ctx=fmpz_mpoly_ctx.get(names, Ordering.lex)
G,C0,Q=[ctx.from_dict(x) for x in pickle.load(open('G.pkl','rb'))]
D=ctx.from_dict(pickle.load(open('D.pkl','rb')))
def coef_z(P,k):
    out={}
    for ex,co in P.to_dict().items():
        if ex[8]==k:
            e2=list(ex); e2[8]=0; out[tuple(e2)]=co
    return ctx.from_dict(out)
Dk=[coef_z(D,k) for k in range(5)]
# Q^3 D'(z0) = sum k D_k (-C0)^(k-1) Q^(4-k)
S=ctx.constant(0)
for k in range(1,5):
    S=S+k*Dk[k]*(-C0)**(k-1)*Q**(4-k)
print('S terms',len(S.to_dict()))
# claim: D'(z0) = 4 G P6 / Q  => Q^3 D'(z0) = 4 G P6 Q^2
qq,rr=divmod(S,G*Q*Q*4)
print('rem zero',rr.is_zero(),'P6 terms',len(qq.to_dict()),'deg',qq.total_degree())
# sign: G sign convention (gamma<0). check also with -G
