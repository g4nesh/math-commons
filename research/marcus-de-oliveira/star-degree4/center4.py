"""(a) Scan of the reduced degree-4 problem in t.  (b) Central solutions (zero objective) at chosen t: eigenvalue
histogram of the reduced blocks, to see how many further directions are forced to zero."""
import sys, time, pickle, numpy as np
from reduced4 import *
def blocks_from_x(R, x):
    out = []
    for s in R.slices:
        if s is None: out.append(None); continue
        lo, hi, r = s; S = np.zeros((r, r)); iu = np.triu_indices(r); S[iu] = x[lo:hi]; S = S + S.T - np.diag(np.diag(S)); out.append(S)
    return out
if __name__ == '__main__':
    t0 = time.time(); R = Reduced(load_comp())
    scan = [float(v) for v in sys.argv[1].split(',')] if sys.argv[1] != '-' else []
    for t in scan:
        lam, st, out, x = R.solve(R.f(t)); print("scan t = %-5g best lambda = %s (%s)  [%.0fs]" % (t, None if lam is None else "%+.3e" % lam, st, time.time() - t0), flush=True)
    for t in [float(v) for v in sys.argv[2].split(',')]:
        lam, st, out, x = R.solve(R.f(t), objective='center')
        if x is None: print("centre t = %g: %s" % (t, st)); continue
        Ss = blocks_from_x(R, x); eig = np.concatenate([np.linalg.eigvalsh(S) for S in Ss if S is not None]); lg = np.log10(np.maximum(np.abs(eig), 1e-300))
        print("centre t = %g (%s): %d eigenvalues, min %.2e ; by decade:" % (t, st, len(eig), eig.min()), flush=True)
        for lo in range(-14, 3):
            c = ((lg >= lo) & (lg < lo + 1)).sum()
            if c: print("   1e%d .. 1e%d : %d" % (lo, lo + 1, c))
        pickle.dump({'t': t, 'S': Ss}, open('center4_t%g.pkl' % t, 'wb'))
