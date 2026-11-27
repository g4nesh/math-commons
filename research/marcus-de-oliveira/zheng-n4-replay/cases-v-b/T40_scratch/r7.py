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
Rs={'gamma':ga,'ABb-Nja':A*B*be-Nj*al,'ABF++Na':A*B*Fp+N*al,'ABF-+Na':A*B*Fm+N*al}
pickle.dump({k:v.to_dict() for k,v in Rs.items()},open('Rs.pkl','wb'))
ctx2=fmpz_mpoly_ctx.get(['q','b','J','P','V','L','M'],Ordering.lex)
q2,b2,J2,P2,V2,L2,M2=ctx2.gens()
one=ctx2.constant(1)
def highj(R):
    degs=R.degrees()
    D=max(sum((ex[4],ex[7])) for ex in R.to_dict())
    dp=degs[1]
    assert dp==2,dp
    j_=1+J2
    groups={}
    for ex,co in R.to_dict().items():
        kq,kp,kj,_,ke,_,kb,ka,_=ex
        groups.setdefault((kp,kj,ka,ke),{})[(kq,kb)]=co
    out=ctx2.constant(0)
    pw={}
    for (kp,kj,ka,ke),qb in groups.items():
        base=ctx2.constant(0)
        for (kq,kb),co in qb.items():
            base=base+int(co)*q2**kq*b2**kb
        t=P2**kp*j_**(2-kp)*(1+P2)**(2-kp)*j_**kj
        s=ka+ke
        t=t*V2**s*(J2*(1+V2))**(D-s)*(L2+M2)**(D-s)*L2**ka*M2**ke
        out=out+base*t
    return out,D
t0=time.time()
for nm,R in Rs.items():
    S,D=highj(R)
    dct=S.to_dict(); cs=list(dct.values())
    print(nm,'D',D,'terms',len(cs),'allneg',all(x<0 for x in cs), time.time()-t0)
    # part with q=b=V=0
    part=S.subs({'q':0,'b':0,'V':0})
    j_=1+J2
    expect=-j_*(j_+1)*(j_*(1+P2)+P2)*J2**D*(L2+M2)**D
    print('   q=b=V=0 part equals -j(j+1)[j(1+P)+P]J^D(L+M)^D:', part==expect)
