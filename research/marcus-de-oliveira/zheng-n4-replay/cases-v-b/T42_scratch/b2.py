exec(open("b1.py").read().split("T=e+a+b+1")[0])
import time
Fa=(j*(p+e+a+1)-a*(p+e))**2-4*p*e*a*(a-j)
Fb=(p*(j+e+b+1)-b*(j+e))**2-4*j*e*b*(b-p)
det = from_ms(pickle.load(open("det12_54.pkl","rb")))
print("det12", len(det))
Pc, rem = divmod(det, (c+1)**2); print("det/(c+1)^2 exact", rem==0, len(Pc), Pc.total_degree(), Pc.degrees())
def coeffs_v(P, idx):
    out={}
    for ex,cf in P.to_dict().items():
        k=ex[idx]; ex2=list(ex); ex2[idx]=0
        out.setdefault(k,{})[tuple(ex2)]=cf
    return {k: ctx.from_dict(x) for k,x in out.items()}
ic=names.index("c")
cG=coeffs_v(Gd,ic); cP=coeffs_v(Pc,ic)
al,be,ga=cG[2],cG[1],cG[0]; a2,b1,c0=cP[2],cP[1],cP[0]
print(sorted(cG),sorted(cP))
res=Gd.resultant(Pc,"c"); print("res terms",len(res))
W=res.sqrt(); print("W", len(W), W.total_degree())
W14,rem=divmod(W,(q+1)**2); print("W div (q+1)^2", rem==0, len(W14), W14.total_degree())
u1=al*b1-a2*be; u0=al*c0-a2*ga; T_=2*al*u0-be*u1; U1=u1*c+u0
for s,F_,mul in ((1,Fa,b),(-1,Fb,a)):
    lhs=T_+s*2*al*W
    L,rem=divmod(lhs, F_*mul**2*(q+1)**2)
    print("S",s,"=F mult^2 (q+1)^2 L^2:", rem==0, end=" ")
    if rem==0:
        try: L7=L.sqrt(); print("L7 terms",len(L7), L7.total_degree())
        except Exception as ex: print("not square")
    num=Pc*lhs-(U1+s*W)**2
    qq,rem=divmod(num,Gd); print("   congruence", rem==0)
