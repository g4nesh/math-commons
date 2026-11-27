import pickle, time
from flint import fmpz_mpoly_ctx, Ordering
names=['q','p','j','c','e','d','b','a','z']
ctx=fmpz_mpoly_ctx.get(names, Ordering.lex)
G,C0,Q=[ctx.from_dict(x) for x in pickle.load(open('G.pkl','rb'))]
Pst=ctx.from_dict(pickle.load(open('Pst.pkl','rb')))
def collect_c(P):
    out={}
    for ex,co in P.to_dict().items():
        k=ex[3]; e2=list(ex); e2[3]=0
        out.setdefault(k,{})[tuple(e2)]=co
    return {k:ctx.from_dict(v) for k,v in out.items()}
cc=ctx.gens()[3]
cg=collect_c(G); al,be,ga=cg[2],cg[1],cg[0]
t=time.time()
R=Pst
steps=0
while R.degrees()[3]>=2:
    dg=R.degrees()[3]
    lead=collect_c(R)[dg]
    R=al*R-lead*cc**(dg-2)*G
    steps+=1
print('steps',steps,R.degrees()[3],time.time()-t)
pc=collect_c(R)
r1=pc.get(1,ctx.constant(0)); r0=pc.get(0,ctx.constant(0))
print('r1 terms',len(r1.to_dict()),'r0 terms',len(r0.to_dict()),'deg r1',r1.total_degree(),time.time()-t)
pickle.dump((r1.to_dict(),r0.to_dict()),open('r10.pkl','wb'))
