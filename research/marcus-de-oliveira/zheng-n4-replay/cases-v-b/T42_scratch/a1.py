from h import *
rd = pickle.load(open("rd17.pkl","rb"))
R = rd["R"]
ctx,v = ctxs(R.names)
q,p,j,c,e,d,b,a,hh = [v[n] for n in R.names]
G = fromR(R, rd["G"], ctx); C0 = fromR(R, rd["C0"], ctx); Q = fromR(R, rd["Q"], ctx)
print(len(G), len(C0), len(Q))
print("Q", Q)
print("C0", C0)
Gd = G/(c+d)  if False else None
try:
    q_,r_ = divmod(G, c+d); print("div by (c+d)", r_==0, len(q_))
except Exception as ex: print(ex)
print(G.degrees())
