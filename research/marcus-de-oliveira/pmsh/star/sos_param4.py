"""Degree-4 certificate for the whole range t in (0,1]: with t = 1/(1+tau), tau >= 0,
   (1+tau)^2 * [t g0 + (1+t) a + t^2 b]  =  (g0 + 2a + b) + tau (g0 + 3a) + tau^2 a ."""
import sys, time, numpy as np
from sos2 import *
from sos_param import solve_param
t0 = time.time(); L = Level(4, 5000)
print("level built: rank %d  [%.0fs]" % (L.rank, time.time() - t0), flush=True)
for label, bv in [("D[13,14]", b_), ("D[13,23]", bw)]:
    st, gam = solve_param(L, L.X @ (g0 + 2 * a_ + bv), L.X @ (g0 + 3 * a_), L.X @ a_)
    print("degree 4, b = %s: one certificate for ALL t in (0,1]: slack gamma = %+.3e (%s)  [%.0fs]" % (label, gam, st, time.time() - t0), flush=True)
