import sys, numpy as np
from sos_torus import *
L = Level(3, 3000)
bw = W([(1,'D',(0,2),(1,2))])                      # D[13,23]: the second member of the family
for label, bvec in [("D[13,14]", b_), ("D[13,23]", bw)]:
    for t in [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.85, 1.0, 1.5, 10.0, 100.0]:
        st, gam = L.solve([L.X @ (t * g0 + (1 + t) * a_ + t * t * bvec)], maximize_gamma=True)
        print("b = %s  t = %-6g  degree-3 certified lower bound, divided by t: %+.3e" % (label, t, gam / t), flush=True)
