import pickle, time
from flint import fmpz_mpoly_ctx, Ordering
names=['q','p','j','c','e','d','b','a','z']
ctx=fmpz_mpoly_ctx.get(names, Ordering.lex)
G=ctx.from_dict(pickle.load(open('G.pkl','rb'))[0])
dd,facs=pickle.load(open('disc.pkl','rb'))
Fp,Fm=[ctx.from_dict(f) for f,m in facs]
def collect_c(G):
    out={}
    for ex,co in G.to_dict().items():
        k=ex[3]; e2=list(ex); e2[3]=0
        out.setdefault(k,{})[tuple(e2)]=co
    return {k:ctx.from_dict(v) for k,v in out.items()}
cg=collect_c(G)
al,be,ga=[x.subs({'d':1}) for x in (cg[2],cg[1],cg[0])]
Fp,Fm=[x.subs({'d':1}) for x in (Fp,Fm)]
q,p,j,cc,e,d,b,a,z=ctx.gens()
A=1+a+b;B=1+q+e;N=(q+1)*(e+1)*A+B*a*b;Nj=N-j*e*a
R1=ga; R2=A*B*be-Nj*al; R3=A*B*Fp+N*al; R4=A*B*Fm+N*al
# check which sign for F+/F-: determine ordering
def clear_pj(R, mp, mj, subp, subj):
    pass
# generic: substitute p -> P/(1+P) style given by lambda on exponent
def sub_chart(R, specs):
    """specs: dict var -> ('frac', newname) meaning x = X/(1+X); returns poly in ctx with X replacing x
       after multiplying by (1+X)^deg."""
    degs=R.degrees()
    idx={n:i for i,n in enumerate(names)}
    out=ctx.constant(0)
    X={n:ctx.gens()[idx[n]] for n in names}
    cache={}
    for ex,co in R.to_dict().items():
        t=ctx.constant(int(co))
        for i,n in enumerate(names):
            k=ex[i]
            if k==0 and not (n in specs): continue
            if n in specs:
                dg=degs[i]
                t=t*X[n]**k*(1+X[n])**(dg-k)
            else:
                t=t*X[n]**k
        out=out+t
    return out
for nm,R in [('gamma',R1),('ABb-Nja',R2),('ABF++Na',R3),('ABF-+Na',R4)]:
    S=sub_chart(R,{'p':1,'j':1})
    dct=S.to_dict()
    cs=list(dct.values())
    print(nm,len(cs),'all neg',all(x<0 for x in cs),'const',dct.get((0,)*9))
