"""Second facial-reduction step at degree 4 for a fixed t.
From the central solution (center4_t<t>.pkl): rationalise the null directions of each reduced block and remove them.
Output amplitudes4_t<t>.json: rows A_b (rational) with the certificate sought as f_t = sum_b (A_b m_b)^* T_b (A_b m_b), T_b > 0.
Then test strict feasibility on these amplitude spaces."""
import sys, json, pickle, time, numpy as np, sympy as sp
from fractions import Fraction as Fr
from space4 import *
from reduced4 import Reduced
def rref_cp(K, tol=1e-3):
    K = K.copy(); used = []; nr = 0
    for r in range(K.shape[0]):
        sub = np.abs(K[r:, :]).copy(); sub[:, used] = 0
        i, j = np.unravel_index(np.argmax(sub), sub.shape); i += r
        if abs(K[i, j]) < tol: break
        K[[r, i]] = K[[i, r]]; K[r] /= K[r, j]
        for q in range(K.shape[0]):
            if q != r: K[q] -= K[q, j] * K[r]
        used.append(j); nr += 1
    return K[:nr]
if __name__ == '__main__':
    tq = Fr(sys.argv[1]); t = float(tq); maxden = int(sys.argv[2]) if len(sys.argv) > 2 else 12; t0 = time.time()
    comp = json.load(open('complement4.json'))['complement']; S = pickle.load(open('center4_t%g.pkl' % t, 'rb'))['S']
    out = []; worst = 0; dropped = 0; vals = {}
    for c, Sb in zip(comp, S):
        if not c: out.append([]); continue
        Cq = sp.Matrix([[sp.Rational(q) for q in row] for row in c])
        w, V = np.linalg.eigh(Sb); k = int((w < 1e-6).sum())
        if k == 0: out.append([[str(x) for x in Cq.row(i)] for i in range(Cq.rows)]); continue
        R = rref_cp(V[:, :k].T); assert R.shape[0] == k
        Rq = [[Fr(float(x)).limit_denominator(maxden) if abs(x) > 2e-3 else Fr(0) for x in row] for row in R]
        worst = max(worst, max(abs(float(q) - x) for rq, rr in zip(Rq, R) for q, x in zip(rq, rr))); dropped += k
        for row in Rq:
            for q in row: vals[q] = vals.get(q, 0) + 1
        ns = sp.Matrix([[sp.Rational(q.numerator, q.denominator) for q in row] for row in Rq]).nullspace()
        if not ns: out.append([]); continue
        rowsW = []
        for v in ns:
            den = sp.ilcm(*[sp.fraction(x)[1] for x in v]); rowsW.append(list((v * den).T))
        A = sp.Matrix(rowsW) * Cq
        out.append([[str(x) for x in A.row(i)] for i in range(A.rows)])
    print("t = %s: dropped %d directions ; amplitudes kept %d ; max rounding move %.1e" % (tq, dropped, sum(len(a) for a in out), worst))
    print("   distinct entries of the null directions:", {str(k): v for k, v in sorted(vals.items())})
    json.dump({'t': str(tq), 'amplitudes': out}, open('amplitudes4_t%g.json' % t, 'w'))
    amp = [np.array([[float(Fr(q)) for q in row] for row in a]) if a else None for a in out]
    R2 = Reduced(amp); lam, st, outside, x = R2.solve(R2.f(t))
    print("strict feasibility on the amplitude spaces: unknowns %d, conditions %d, f outside the range %.1e ; best lambda = %+.4e (%s)  [%.0fs]" % (R2.N, R2.rank, outside, lam, st, time.time() - t0))
    pickle.dump({'t': t, 'x': x, 'lam': lam, 'slices': R2.slices}, open('strict4_t%g.pkl' % t, 'wb'))
