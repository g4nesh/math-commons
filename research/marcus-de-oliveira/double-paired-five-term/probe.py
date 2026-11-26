"""Bounded numerical falsification probe; floating output is not proof."""
import numpy as np,itertools,time
from scipy.linalg import expm
rng=np.random.default_rng(71);subs=list(itertools.combinations(range(6),3))
def hist(s):return tuple(sum(i//2==j for i in s) for j in range(3))
terms=[]
for alpha,beta,c in [((1,0,2),(1,1,1),-.25),((1,1,1),(1,2,0),.25),((0,1,2),(2,1,0),.5),((0,2,1),(2,0,1),.5)]:
 for r in subs:
  for s in subs:
   if hist(r)==alpha and hist(s)==beta:terms.append((r,s,c))
rr=np.array([r for r,s,c in terms]);ss=np.array([s for r,s,c in terms]);weights=np.array([c for r,s,c in terms]);best=(1,None)
def fg(U,grad=True):
 M=U[rr[:,:,None],ss[:,None,:]];a,b,c=M[:,0,:],M[:,1,:],M[:,2,:]
 cof=np.stack([np.cross(b,c),np.cross(c,a),np.cross(a,b)],axis=1);det=np.sum(a*cof[:,0,:],axis=1)
 f=np.sum(weights*abs(det)**2).real
 if not grad:return f
 G=np.zeros((6,6),complex)
 values=weights[:,None,None]*det[:,None,None]*cof.conj()
 for i in range(3):
  for j in range(3):np.add.at(G,(rr[:,i],ss[:,j]),values[:,i,j])
 return f,G
start=time.time()
for run in range(12):
 U,R=np.linalg.qr(rng.normal(size=(6,6))+1j*rng.normal(size=(6,6)))
 for k in range(120):
  f,G=fg(U);K=G@U.conj().T-U@G.conj().T;norm=np.linalg.norm(K)**2
  if norm<1e-20:break
  eta=4.
  while eta>1e-8:
   trial=expm(-eta*K)@U;new=fg(trial,False)
   if new<f-eta*norm*.1:break
   eta*=.5
  U=trial
  if f<best[0]:best=(f,U.copy())
 print(run,f,flush=True)
 if f < -1e-8:break
print('best',best[0],'seconds',time.time()-start)
# No numerical matrix is promoted to a certified witness.
