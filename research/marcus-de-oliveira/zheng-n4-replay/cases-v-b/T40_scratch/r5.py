import pickle, time
from flint import fmpz_mpoly_ctx, Ordering
names=['q','p','j','c','e','d','b','a','z']
ctx=fmpz_mpoly_ctx.get(names, Ordering.lex)
G=ctx.from_dict(pickle.load(open('G.pkl','rb'))[0])
def collect_c(G):
    out={}
    for ex,co in G.to_dict().items():
        k=ex[3]; e2=list(ex); e2[3]=0
        out.setdefault(k,{})[tuple(e2)]=co
    return {k:ctx.from_dict(v) for k,v in out.items()}
cg=collect_c(G)
al,be,ga=cg[2],cg[1],cg[0]
import random
pt={'q':3,'p':2,'j':5,'e':7,'d':1,'b':11,'a':13}
def ev(P): return P.subs(pt)
print('const terms:', [int(x) for x in [ev(al).to_dict().get((0,)*9,0), ev(be).to_dict().get((0,)*9,0), ev(ga).to_dict().get((0,)*9,0)]])
# gamma at d=1 coefficient: sign
g1=ga.subs({'d':1})
print('gamma const', g1.to_dict().get((0,)*9,0))
disc=be*be-4*al*ga
t=time.time()
cont,fac=disc.factor()
print(cont,[(len(f.to_dict()),m,f.total_degree()) for f,m in fac], time.time()-t)
pickle.dump((disc.to_dict(),[(f.to_dict(),m) for f,m in fac]),open('disc.pkl','wb'))
