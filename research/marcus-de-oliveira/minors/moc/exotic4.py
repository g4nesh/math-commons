import sys; sys.path.insert(0,'moc')
import numpy as np, itertools
from scipy.optimize import linprog
from reduced4 import PERMS, perm_matrix
from moclib import compound_vector, haar_unitary, subsets
F = np.load('moc/facets4.npy')
from reduced4 import VERTS
onv = np.array([(np.abs(VERTS@f[:-1]+f[-1])<1e-9) for f in F])
ex = [i for i in range(len(F)) if onv[i].sum()==13]
f = F[ex[0]]
vals = VERTS@f[:-1]+f[-1]
w = np.round(-vals/np.abs(vals[np.abs(vals)>1e-9]).min()).astype(int)  # weights 0,1,2
print("weights w_sigma:", dict(zip(PERMS,w)))
# labels for full coordinates: const, M1 (16), M2 (36)
labels = ['1'] + ['m%d%d'%(i,j) for i in range(4) for j in range(4)] + \
   ['D%s|%s'%(''.join(map(str,S)),''.join(map(str,T))) for S in subsets(4,2) for T in subsets(4,2)]
X = np.array([np.concatenate([[1.0], compound_vector(perm_matrix(s),[1,2])]) for s in PERMS])
d = X.shape[1]
# min L1 coefficients c with X c = w : c = cp - cm
A_eq = np.hstack([X, -X]); res = linprog(np.ones(2*d), A_eq=A_eq, b_eq=w, bounds=[(0,None)]*(2*d), method='highs')
c = res.x[:d]-res.x[d:]
print("L1 norm", res.fun)
print("g(x) = ", ' '.join('%+g*%s'%(round(ci,6),l) for ci,l in zip(c,labels) if abs(ci)>1e-9))
np.save('moc/exotic_g.npy', c)
# check on random unitaries
rng=np.random.default_rng(2); mins=[]
for _ in range(20000):
    U=haar_unitary(4,rng); x=np.concatenate([[1.0],compound_vector(U,[1,2])]); mins.append(x@c)
mins=np.array(mins); print("min g over Haar samples %.4f, frac<0 %.4f"%(mins.min(), (mins<0).mean()))
