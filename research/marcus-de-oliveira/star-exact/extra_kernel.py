"""Second facial-reduction step: the kernel of the central reduced solution, block by block, in monomial coordinates.
Prints the cubic polynomials  c . m_b  that span the extra forced kernel (those that must also be dropped)."""
import sys, pickle, numpy as np
from fractions import Fraction as Fr
from exact_space import MONS, BLOCK_LIST
from reduced_sdp import load_comp
from reduce import rref
t = float(sys.argv[1]); comp = load_comp(1.0); S = pickle.load(open('center_t%g.pkl' % t, 'rb'))['S']
nm = lambda v: "u%d%d" % (v // 4 + 1, v % 4 + 1)
tot = 0; dens = set(); worst = 0
for bi, (b, c, Sb) in enumerate(zip(BLOCK_LIST, comp, S)):
    if Sb is None: continue
    w, V = np.linalg.eigh(Sb); k = int((w < 1e-6).sum())
    if k == 0: continue
    tot += k; C = np.array([[float(q) for q in row] for row in c])
    # kernel vectors of S_b in reduced coordinates; the corresponding polynomials are  (V_k^T C) m
    R, piv = rref(V[:, :k].T)
    Rq = np.array([[float(Fr(float(x)).limit_denominator(200)) for x in row] for row in R]); worst = max(worst, np.abs(Rq - R).max())
    dens |= {Fr(float(x)).limit_denominator(200).denominator for x in R.ravel()}
    if tot <= 12:
        print("block %d (size %d, complement %d): %d extra kernel direction(s) in reduced coordinates:" % (bi, len(b), len(c), k))
        for row in R: print("     coefficients on the complement basis:", np.round(row, 6))
        for row in c: print("     complement basis poly: ", " ".join("%+d*%s" % (int(q), "*".join(nm(v) for v in MONS[b[i]])) for i, q in enumerate(row) if q != 0))
print("extra kernel dimension in total: %d ; max distance of RREF entries from rationals with denominator <= 200: %.1e ; denominators seen: %s" % (tot, worst, sorted(dens)))
