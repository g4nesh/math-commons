import pickle, time
from flint import fmpz_mpoly_ctx, Ordering
names=['q','p','j','c','e','d','b','a','z']
ctx=fmpz_mpoly_ctx.get(names, Ordering.lex)
r1=ctx.from_dict(pickle.load(open('r10.pkl','rb'))[0])
t=time.time()
sq=r1.factor_squarefree()
print('squarefree',[(len(f.to_dict()),m,f.total_degree()) for f,m in sq[1]],sq[0],time.time()-t,flush=True)
f=r1.factor()
print('factor',f[0],[(len(x.to_dict()),m,x.total_degree()) for x,m in f[1]],time.time()-t,flush=True)
pickle.dump((f[0],[(x.to_dict(),m) for x,m in f[1]]),open('r1fac.pkl','wb'))
