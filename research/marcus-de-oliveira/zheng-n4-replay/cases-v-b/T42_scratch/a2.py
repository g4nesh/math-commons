from h import *
rd = pickle.load(open("rd17.pkl","rb"))
R = rd["R"]
ctx,v = ctxs(R.names)
q,p,j,c,e,d,b,a,hh = [v[n] for n in R.names]
G = fromR(R, rd["G"], ctx); C0 = fromR(R, rd["C0"], ctx); Q = fromR(R, rd["Q"], ctx)
Gd, rem = divmod(G, c+d); assert rem==0
one = d
S=a+b; A=a*b; Rr=(a+d)*(b+d); X=a+q+d; Y=b+e+d; T=e+q+d
V = A + (S+d)*(p+d)
D = X*Y*(c+d) + T*c*j
r00 = p*p+(2*d-S)*p+2*d*d-A  # homogeneous deg 2 using d
# homogenise: r00 = p^2 + (2d - S)p + 2d^2 - A
R0 = X*Y*r00 + a*(a+p+2*d)*q*Y + b*(b+p+2*d)*e*X - (2*V-A)*e*q
Fe = ((e+d+b)*(p-a)+e*b)**2 - 4*b*p*(e-a)*(e+d+b)
Fq = ((q+d+a)*(p-b)+q*a)**2 - 4*a*p*(q-b)*(q+d+a)
num = R0*R0 - Fe*Fq
A0, rem = divmod(num, 4*V); print("A0 exact", rem==0, len(A0), len(R0))
print("identity", Gd - (V*D*D + R0*j*D + A0*j*j) == 0)
d1 = {d:1}
# at d=1 counts
sub = lambda P: P.subs({'d':1}) if False else None
print(R0.degrees(), A0.degrees())
