import numpy as np
rng=np.random.default_rng(1)
def fk(a,b,p,x,y):
    S=a+b; A=a*b
    V=A+(S+1)*(p+1)
    f=(p-a+b*x)**2-4*b*p/(b+1)*((S+1)*x-a)
    k=(p-b+a*y)**2-4*a*p/(a+1)*((S+1)*y-b)
    r00=p*p+(2-S)*p+2-A
    r=r00+a*(a+p+2)*y+b*(b+p+2)*x-(2*V-A)*x*y
    return f,k,r
N=4_000_000
worst=-1e9
for it in range(5):
    u=rng.random(N)<0.9; S=np.where(u,rng.uniform(0,2,N),rng.choice([0.0,1.999999,0.001,1.0],N))
    t=rng.random(N); 
    t=np.where(rng.random(N)<0.05,0.0,t); t=np.where(rng.random(N)<0.05,1.0,t)
    a=S*t; b=S-a
    p=np.exp(rng.normal(0,2,N)); p=np.where(rng.random(N)<0.05, 1e-9, p)
    x=rng.random(N); y=rng.random(N)
    x=np.where(rng.random(N)<0.15,rng.choice([0.0,1.0],N),x); y=np.where(rng.random(N)<0.15,rng.choice([0.0,1.0],N),y)
    f,k,r=fk(a,b,p,x,y)
    m=(f>=0)&(k>=0)
    B=np.sqrt(np.where(m,f*k,0))-r
    B=np.where(m,B,-1e9)
    i=np.argmax(B)
    if B[i]>worst: worst=B[i]; arg=(a[i],b[i],p[i],x[i],y[i])
    print(it,m.mean(),B[i], (a[i],b[i],p[i],x[i],y[i]))
