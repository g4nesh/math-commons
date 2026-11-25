"""Second facial-reduction step at a fixed rational t.
Input : complement_t1.json (cubics vanishing on the zero set) and center_t<t>.pkl (a central reduced solution).
Output: amplitudes_t<t>.json -- for each weight block, integer/rational rows A_b such that the certificate is sought as
        f_t = sum_b (A_b m_b)^* T_b (A_b m_b)  with T_b positive DEFINITE.
The extra kernel directions of the central solution are rationalised (complete-pivot row reduction, small denominators)."""
import sys, json, pickle, numpy as np, sympy as sp
from fractions import Fraction as Fr
from exact_space import BLOCK_LIST
from reduced_sdp import load_comp
def rref_cp(K, tol=1e-3):
    """row reduction with complete pivoting (entries stay bounded by 1 in pivot columns)"""
    K = K.copy(); rows = []; used = []
    for r in range(K.shape[0]):
        sub = np.abs(K[r:, :]).copy(); sub[:, used] = 0
        i, j = np.unravel_index(np.argmax(sub), sub.shape); i += r
        if abs(K[i, j]) < tol: break
        K[[r, i]] = K[[i, r]]; K[r] /= K[r, j]
        for q in range(K.shape[0]):
            if q != r: K[q] -= K[q, j] * K[r]
        used.append(j); rows.append(r)
    return K[:len(rows)]
if __name__ == '__main__':
    t = float(sys.argv[1]); comp = load_comp(1.0); S = pickle.load(open('center_t%g.pkl' % t, 'rb'))['S']
    out = []; worst = 0; dropped = 0; dens = set()
    for b, c, Sb in zip(BLOCK_LIST, comp, S):
        if not c: out.append([]); continue
        Cq = sp.Matrix([[sp.Rational(q.numerator, q.denominator) for q in row] for row in c]); r = len(c)
        w, V = np.linalg.eigh(Sb); k = int((w < 1e-6).sum())
        if k == 0: A = Cq
        else:
            R = rref_cp(V[:, :k].T)
            Rq = [[Fr(float(x)).limit_denominator(12) if abs(x) > 2e-3 else Fr(0) for x in row] for row in R]
            worst = max(worst, max(abs(float(q) - x) for rq, rr in zip(Rq, R) for q, x in zip(rq, rr)))
            dens |= {q.denominator for row in Rq for q in row}; dropped += k
            Kq = sp.Matrix([[sp.Rational(q.numerator, q.denominator) for q in row] for row in Rq])
            ns = Kq.nullspace()                                   # vectors w with K w = 0 : allowed amplitude combinations
            if not ns: out.append([]); continue
            Wm = sp.Matrix.hstack(*ns).T
            Wm = sp.Matrix([list(row * sp.ilcm(*[sp.fraction(x)[1] for x in row])) for row in [Wm.row(i) for i in range(Wm.rows)]])
            A = Wm * Cq
        out.append([[str(x) for x in A.row(i)] for i in range(A.rows)])
    tot = sum(len(a) for a in out)
    print("t = %g: dropped %d directions ; amplitudes kept: %d ; max rounding move %.1e ; denominators used %s" % (t, dropped, tot, worst, sorted(dens)))
    json.dump({'t': t, 'amplitudes': out}, open('amplitudes_t%g.json' % t, 'w'))
