from a3 import *
names = ("U","Z","P","xi","eta")
ctx2, v2 = ctxs(names)
import a3
a3.ctx=ctx2; a3.one = v2["U"]**0
a3.RF.__init__.__globals__['one']=a3.one
U_,Z_,P_,xi2,eta2 = [v2[n] for n in names]
one=a3.one
RF=a3.RF
def run(layer):
    rho = 1+RF(U_)
    D = rho*rho-rho+1
    Z = RF(Z_)
    z = (rho + D*Z)/(D*(1+Z))
    a = z*rho; b = z/rho
    k0 = rho*RF(U_)**2/(D*(1+Z)*(rho+D*Z))
    PP = RF(P_)
    if layer==1:
        p = k0 + PP
    else:
        p = k0*PP/(1+PP)
    x = RF(xi2,1+xi2); y = RF(eta2,1+eta2)
    S=a+b; Aa=a*b; R=(a+1)*(b+1)
    V = Aa+(S+1)*(p+1)
    h = 2*V*PP if layer==1 else RF(one*0)
    r00 = p*p+(2-S)*p+2-Aa
    r = r00 + a*(a+p+2)*y + b*(b+p+2)*x - (2*V-Aa)*x*y
    f = (p-a+b*x)**2 - 4*b*p/(b+1)*((S+1)*x - a)
    kk = (p-b+a*y)**2 - 4*a*p/(a+1)*((S+1)*y - b)
    H = R*(h+r)
    E = H*H - R*((b+1)*f)*((a+1)*kk)
    return H,E
import time
for layer in (1,2):
    t=time.time()
    H,E=run(layer)
    print(layer, nn(H.n), nn(H.d), nn(E.n), nn(E.d), time.time()-t)
