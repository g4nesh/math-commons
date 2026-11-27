from dets import *
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
H=hessian_matrix(u,ctx)
t=time.time()
D=det_dp(H)
print('det terms',len(D.to_dict()),time.time()-t)
import pickle
pickle.dump(D.to_dict(),open('D.pkl','wb'))
