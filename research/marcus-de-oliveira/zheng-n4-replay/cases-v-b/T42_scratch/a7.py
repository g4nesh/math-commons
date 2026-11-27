exec(open("a6.py").read().split("t=time.time()\nres =")[0])
import time
# (a) congruences
res = Gd.resultant(Pj, "j"); W = res.sqrt()
u1 = al*b1 - a2*be; u0 = al*c0 - a2*ga
T_ = 2*al*u0 - be*u1
U_ = u1*j+u0
for s in (1,-1):
    num = Pj*(T_+s*2*al*W) - (U_+s*W)**2
    qq, rem = divmod(num, Gd)
    print("congruence sign",s, rem==0, len(qq) if rem==0 else None)
# (b) bordered
order=[3,1,53,6,15,17,16,4]
Hh = M.hessian_linear_forms(order)
rd = pickle.load(open("rd17.pkl","rb")); R=rd["R"]
Cg = from_packed(R, rd["C0"]); Qg = from_packed(R, rd["Q"])
print("C0 terms", len(Cg), "Q", Qg)
nm9=["q","p","j","c","e","d","b","a","h"]
ctx9,v9=ctxs(("q","p","j","c","e","b","a","h"))
def from_ms9(poly):
    out={}
    for m,cf in poly.items():
        ex=[0]*8
        for i in m:
            n=nm9[i]
            if n=="d": continue
            ex[("q","p","j","c","e","b","a","h").index(n)]+=1
        out[tuple(ex)]=out.get(tuple(ex),0)+cf
    return ctx9.from_dict({k:x for k,x in out.items() if x})
for k in (0,1):
    idx=[k]+list(range(2,14))
    t=time.time()
    sub=[[Hh[i][kk] for kk in idx] for i in idx]
    det,l1=M.det_linear_forms(sub,9)
    Dp=from_ms9(det)
    print("k",k,"det terms",len(Dp), "h-degree", Dp.degrees()[7], time.time()-t)
    # coefficients in h
    cf={}
    for ex,c_ in Dp.to_dict().items():
        ex2=list(ex); hh=ex2[7]; ex2[7]=0
        cf.setdefault(hh,{})[tuple(ex2)]=c_
    cfp={hh: ctx.from_dict({kk:x for kk,x in d_.items()}) if False else None for hh,d_ in cf.items()}
    # build in ctx (7 vars): ex2 has 8 entries with last 0 -> trim
    cfp={hh: ctx.from_dict({kk[:7]:x for kk,x in d_.items()}) for hh,d_ in cf.items()}
    Q=Qg; C0=Cg
    Nk = sum(cfp[m]*(-C0)**m*Q**(2-m) for m in cfp)
    print(" numerator terms", len(Nk))
    # divisible by Q^? and G
    qq,rem=divmod(Nk, Gd*(c+1)); print(" divisible by G_mine:", rem==0)
    qq2,rem2=divmod(Nk, Gd); print(" divisible by G_doc:", rem2==0, len(qq2) if rem2==0 else None)
    Qq,remq=divmod(Nk,Q); print(" divisible by Q:", remq==0, len(Qq) if remq==0 else None)
