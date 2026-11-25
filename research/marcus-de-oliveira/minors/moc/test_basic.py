import numpy as np, sys
sys.path.insert(0, 'moc')
from moclib import *
rng = np.random.default_rng(0)
# 1. Cauchy-Binet expansion check
for n in [2,3,4,5]:
    errs=[]
    for _ in range(20):
        a = rng.standard_normal(n)+1j*rng.standard_normal(n)
        b = rng.standard_normal(n)+1j*rng.standard_normal(n)
        U = haar_unitary(n, rng)
        errs.append(abs(det_sum(a,b,U)-cauchy_binet(a,b,U))/abs(det_sum(a,b,U)))
    print("CB expansion n=%d max rel err %.2e"%(n,max(errs)))
# 2. compound Birkhoff LP
for n in [3,4,5]:
    perms, V = compound_perm_vectors(n)
    vals=[]
    for _ in range(200 if n<5 else 50):
        U = haar_unitary(n, rng)
        t,_ = compound_birkhoff_lp(U, V)
        vals.append(t)
    vals=np.array(vals)
    print("n=%d compound-Birkhoff LP: max infeasibility %.3e, frac infeasible(>1e-9) %.3f"%(n, vals.max(), (vals>1e-9).mean()))
