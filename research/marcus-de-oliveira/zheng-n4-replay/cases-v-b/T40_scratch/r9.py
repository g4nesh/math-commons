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
H=hessian_matrix(u,ctx)
t=time.time()
sub=[[H[i][k] for k in range(2,14)] for i in range(2,14)]
Pst=det_dp(sub)
dct=Pst.to_dict()
print('P* terms',len(dct),'degc',Pst.degrees()[3],'zdeg',Pst.degrees()[8],time.time()-t)
pickle.dump(dct,open('Pst.pkl','wb'))
# bordered minors I u {k}
for k in (0,1):
    idx=[k]+list(range(2,14))
    sub=[[H[i][m] for m in idx] for i in idx]
    Bk=det_dp(sub)
    print('bordered',k,len(Bk.to_dict()),'zdeg',Bk.degrees()[8],time.time()-t)
    pickle.dump(Bk.to_dict(),open('Bord%d.pkl'%k,'wb'))
