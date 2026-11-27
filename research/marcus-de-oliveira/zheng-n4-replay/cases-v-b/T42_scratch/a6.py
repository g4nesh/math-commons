from h import *
import time
src="17"
names=("q","p","j","c","e","b","a")
ctx,v=ctxs(names)
q,p,j,c,e,b,a=[v[n] for n in names]
def from_ms(poly):   # multiset dict: {tuple of var idx: coeff} for vars q,p,j,c,e,d,b,a (index 0..7)
    nm8=["q","p","j","c","e","d","b","a"]
    out={}
    for m,cf in poly.items():
        ex=[0]*7
        for i in m:
            n=nm8[i]
            if n=="d": continue
            ex[names.index(n)]+=1
        out[tuple(ex)]=out.get(tuple(ex),0)+cf
    return ctx.from_dict({k:x for k,x in out.items() if x})
def from_packed(R, pl):
    out={}
    for k,cf in pl.items():
        ex=R.exps(k)
        nm=R.names
        e7=[0]*7
        for n,x in zip(nm,ex):
            if n in ("d","h"): 
                assert n!="h" or x==0
                continue
            e7[names.index(n)]+=x
        out[tuple(e7)]=out.get(tuple(e7),0)+cf
    return ctx.from_dict({k:x for k,x in out.items() if x})
rd = pickle.load(open("rd17.pkl","rb")); R=rd["R"]
G = from_packed(R, rd["G"]); 
Gd, rem = divmod(G, c+1); print("Gd rem", rem==0, len(Gd))
det = from_ms(pickle.load(open("det12_17.pkl","rb")))
print("det12", len(det), det.degrees())
Pj, rem = divmod(det, j*j); print("det/j^2 exact", rem==0, len(Pj), Pj.degrees(), Pj.total_degree())
# coefficients in j
def coeffs_j(P):
    d = P.to_dict()
    out={}
    for ex,cf in d.items():
        k=ex[2]; ex2=list(ex); ex2[2]=0
        out.setdefault(k,{})[tuple(ex2)]=cf
    return {k: ctx.from_dict(x) for k,x in out.items()}
cG=coeffs_j(Gd); cP=coeffs_j(Pj)
al,be,ga=cG[2],cG[1],cG[0]; a2,b1,c0=cP[2],cP[1],cP[0]
print(len(al),len(be),len(ga),len(a2),len(b1),len(c0), sorted(cG), sorted(cP))
t=time.time()
res = Gd.resultant(Pj, "j") if False else None
t=time.time()
res = Gd.resultant(Pj, "j")
print("resultant terms", len(res), res.total_degree(), time.time()-t)
W2 = res
try:
    W = W2.sqrt(); print("sqrt ok", len(W), W.total_degree())
except Exception as ex:
    print("sqrt fail", ex)
    W=None
if W is not None:
    print("W divisible by (c+1)^2", divmod(W,(c+1)**2)[1]==0)
    W14 = divmod(W,(c+1)**2)[0]; print(len(W14), W14.total_degree())
    u1 = al*b1 - a2*be; u0 = al*c0 - a2*ga
    T_ = 2*al*u0 - be*u1
    Fe = ((e+1+b)*(p-a)+e*b)**2 - 4*b*p*(e-a)*(e+1+b)
    Fq = ((q+1+a)*(p-b)+q*a)**2 - 4*a*p*(q-b)*(q+1+a)
    for sgn,F_,mult in ((1,Fe,(q+1+a)),(-1,Fq,(e+1+b))):
        lhs = T_ + sgn*2*al*W
        rhs0 = F_*(c+1)**2*mult**2
        L, rem = divmod(lhs, rhs0)
        print("sign",sgn,"divisible by F(c+1)^2mult^2:", rem==0, end=" ")
        if rem==0:
            try:
                L7=L.sqrt(); print("L sqrt ok", len(L7), L7.total_degree())
            except Exception as ex: print("not square", ex)
        else: print()
