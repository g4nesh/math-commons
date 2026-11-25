import sys; sys.path.insert(0,'moc')
import numpy as np
from bsearch import *
rng=np.random.default_rng(5)
for n in [3,4,5]:
    pr=Problem(n); N=3
    U=haar(N,n,rng); a=rng.standard_normal((N,n))+1j*rng.standard_normal((N,n)); b=rng.standard_normal((N,n))+1j*rng.standard_normal((N,n)); th=rng.uniform(0,6,N)
    tau=0.05
    J,(gU,gA,gB,gth)=pr.evaluate(U,a,b,th,tau)
    eps=1e-6
    # directional checks
    da=rng.standard_normal((N,n))+1j*rng.standard_normal((N,n))
    Jp,_=pr.evaluate(U,a+eps*da,b,th,tau,False); Jm,_=pr.evaluate(U,a-eps*da,b,th,tau,False)
    print(n,'a', np.max(np.abs((Jp-Jm)/(2*eps) - (np.conj(gA)*da).real.sum(1)*0 - (gA.real*da.real+gA.imag*da.imag).sum(1))))
    db=rng.standard_normal((N,n))+1j*rng.standard_normal((N,n))
    Jp,_=pr.evaluate(U,a,b+eps*db,th,tau,False); Jm,_=pr.evaluate(U,a,b-eps*db,th,tau,False)
    print(n,'b', np.max(np.abs((Jp-Jm)/(2*eps) - (gB.real*db.real+gB.imag*db.imag).sum(1))))
    H=rng.standard_normal((N,n,n))+1j*rng.standard_normal((N,n,n)); H=(H-np.conj(np.swapaxes(H,1,2)))/2
    Jp,_=pr.evaluate(U@expm_skew(eps*H),a,b,th,tau,False); Jm,_=pr.evaluate(U@expm_skew(-eps*H),a,b,th,tau,False)
    print(n,'U', np.max(np.abs((Jp-Jm)/(2*eps) - (np.conj(gU)*H).real.sum((1,2)))))
    Jp,_=pr.evaluate(U,a,b,th+eps,tau,False); Jm,_=pr.evaluate(U,a,b,th-eps,tau,False)
    print(n,'th', np.max(np.abs((Jp-Jm)/(2*eps) - gth)))
