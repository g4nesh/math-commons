from h import *
import time
names = ("A","B","p","xi","eta")
ctx, v = ctxs(names)
A_,B_,p,xi,eta = [v[n] for n in names]
one = p**0

class RF:
    def __init__(s, n, d=None):
        d = one if d is None else d
        g = n.gcd(d)
        if not g.is_one():
            n, d = n//g, d//g
        if d.coefficient(0) < 0: n, d = -n, -d
        s.n, s.d = n, d
    def _c(s,o): return o if isinstance(o,RF) else RF(one*o)
    def __add__(s,o): o=s._c(o); return RF(s.n*o.d+o.n*s.d, s.d*o.d)
    __radd__=__add__
    def __neg__(s): return RF(-s.n,s.d)
    def __sub__(s,o): return s+(-s._c(o))
    def __rsub__(s,o): return s._c(o)-s
    def __mul__(s,o): o=s._c(o); return RF(s.n*o.n, s.d*o.d)
    __rmul__=__mul__
    def __truediv__(s,o): o=s._c(o); return RF(s.n*o.d, s.d*o.n)
    def __pow__(s,k):
        r=RF(one)
        for _ in range(k): r=r*s
        return r
def nn(P): 
    return len(P), sum(1 for i in range(len(P)) if P.coefficient(i)<0)
def chart(bR, aR):
    a=aR; b=bR
    x = RF(xi, 1+xi); y = RF(eta, 1+eta)
    S=a+b; Aa=a*b; R=(a+1)*(b+1); X=a+0; 
    V = Aa+(S+1)*(p+1)
    L = (Aa-1)/R
    h = 2*V*(p+L)
    r00 = -Aa - p*S + p*p + 2*p + 2
    r = r00 + a*(a+p+2)*y + b*(b+p+2)*x - (2*V-Aa)*x*y
    f = (p-a+b*x)**2 - 4*b*p/(b+1)*((a+b+1)*x - a)
    g = (p-b+a*y)**2 - 4*a*p/(a+1)*((a+b+1)*y - b)
    H = R*(h+r)
    E = H*H - (a+1)*(b+1)*((b+1)*f)*((a+1)*g)
    return H, E
t=time.time()
B=RF(B_)
H,E = chart(1+B, 1+B+RF(A_))
print("b>=1", nn(H.n), nn(H.d), nn(E.n), nn(E.d), time.time()-t)
bb = RF(one)/(1+B)
H,E = chart(bb, 1+B+RF(A_))
print("b<=1", nn(H.n), nn(H.d), nn(E.n), nn(E.d), time.time()-t)
