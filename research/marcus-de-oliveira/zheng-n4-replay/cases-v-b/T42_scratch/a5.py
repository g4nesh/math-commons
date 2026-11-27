from h import *
import time
src = sys.argv[1]
rd = pickle.load(open("rd%s.pkl"%src,"rb"))
R = rd["R"]
order = {"17":[3,1,53,6,15,17,16,4],"54":[6,2,5,53,0,54,3,15]}[src]
t=time.time()
Hh = M.hessian_linear_forms(order)   # 14x14 of dict var->coef, variables 0..8
I = list(range(2,14))
sub = [[Hh[i][k] for k in I] for i in I]
assert all(set(d_)<=set(range(8)) for row in sub for d_ in row)
det, l1 = M.det_linear_forms(sub, 8)
print("det terms", len(det), time.time()-t)
pickle.dump(det, open("det12_%s.pkl"%src,"wb"))
