from h import *
import time
names=("q","p","j","c","e","b","a")
ctx,v=ctxs(names)
q,p,j,c,e,b,a=[v[n] for n in names]
nm8=["q","p","j","c","e","d","b","a"]
def from_ms(poly):
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
        e7=[0]*7
        for n,x in zip(R.names,ex):
            if n in ("d","h"):
                assert n!="h" or x==0
                continue
            e7[names.index(n)]+=x
        out[tuple(e7)]=out.get(tuple(e7),0)+cf
    return ctx.from_dict({k:x for k,x in out.items() if x})
rd = pickle.load(open("rd54.pkl","rb")); R=rd["R"]
G = from_packed(R, rd["G"]); C0=from_packed(R, rd["C0"]); Q=from_packed(R, rd["Q"])
Gd, rem = divmod(G, q+1); print("G/(q+1) exact", rem==0, len(Gd), Gd.total_degree())
print("C0", C0); print("Q", Q)
Cdoc = -(q+1)*(a*b+a*c+a+b*c+b+c+c*e+1+e)
Qdoc = -a*b-a*c-a+a*p-b*c-b+b*j-c-1-j*p
print("C0==doc", C0==Cdoc, "Q==doc", Q==Qdoc)
T=e+a+b+1; U=T+a*b; W=U+c*T; V=p*j+(p+j+1)*(e+1)
H=(p-j)**2-2*e*(p+j)+e*e
A1=e*e-e*j+e-j*p-j+p*p+p; B1=e*e-e*p+e+j*j-j*p+j-p
K=V-e*(e+1); B0=-H*a*b+p*B1*a+j*A1*b+p*j*K; C0b=T*(p*j-p*a-j*b)+e*a*b
Fa=(j*(p+e+a+1)-a*(p+e))**2-4*p*e*a*(a-j)
Fb=(p*(j+e+b+1)-b*(j+e))**2-4*j*e*b*(b-p)
# nodal identity: G = (c+1)^2 (V L^2 - B0 L - pje C0), L=(q+1)(T+ab/(c+1))
Lc = (q+1)*((c+1)*T+a*b)     # = (c+1) L
nod = V*Lc**2 - B0*Lc*(c+1) - p*j*e*C0b*(c+1)**2
print("nodal (3):", nod==Gd)
print("(4a) B0^2+4VpjeC0 = FaFb:", B0**2+4*V*p*j*e*C0b == Fa*Fb)
def disc(P,var):
    return P.discriminant(var)
print("Disc_q G:", Gd.discriminant("q") == (c+1)**2*W**2*Fa*Fb)
print("Disc_c G:", Gd.discriminant("c") == a**2*b**2*(q+1)**2*Fa*Fb)
print("Disc_a Fa:", Fa.discriminant("a")==16*e*p*j*j*V, "Disc_b Fb:", Fb.discriminant("b")==16*e*j*p*p*V)
