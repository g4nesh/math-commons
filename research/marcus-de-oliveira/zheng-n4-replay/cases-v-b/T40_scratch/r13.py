import pickle, time
from flint import fmpz_mpoly_ctx, Ordering
names=['q','p','j','c','e','d','b','a','z']
ctx=fmpz_mpoly_ctx.get(names, Ordering.lex)
G=ctx.from_dict(pickle.load(open('G.pkl','rb'))[0])
r1,r0=[ctx.from_dict(x) for x in pickle.load(open('r10.pkl','rb'))]
dd,facs=pickle.load(open('disc.pkl','rb'))
Fp,Fm=[ctx.from_dict(f) for f,m in facs]
cont,lf=pickle.load(open('r1fac.pkl','rb'))
Lm=[ctx.from_dict(x) for x,m in lf]
def collect_c(P):
    out={}
    for ex,co in P.to_dict().items():
        k=ex[3]; e2=list(ex); e2[3]=0
        out.setdefault(k,{})[tuple(e2)]=co
    return {k:ctx.from_dict(v) for k,v in out.items()}
cg=collect_c(G); al,be,ga=cg[2],cg[1],cg[0]
T=2*al*r0-be*r1
t=time.time()
print('T terms',len(T.to_dict()),time.time()-t,flush=True)
X,Y=Lm
print('product r1=XY',cont*X*Y==r1)
for (L,M,nm) in [(X,Y,'L=f1,M=f2'),(Y,X,'L=f2,M=f1')]:
    for (FP,FM,nn) in [(Fp,Fm,'F+,F- as returned'),(Fm,Fp,'swapped')]:
        lhs=2*T; rhs=FP*M*M+FM*L*L
        print(nm,nn,lhs==rhs,flush=True)
