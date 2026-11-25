"""Exact certificate at a rational t:  f_t = t*g0 + (1+t)*D[12,12] + t^2*D[13,14] = sum_b (A_b m_b)^* T_b (A_b m_b) on U(4).
1. Exact linear system for the entries of the T_b from evaluations at rational unitaries.
2. Free variables := rounded values of a numerical strictly feasible point; pivot variables solved exactly.
3. Exact check that every T_b is positive definite (leading principal minors).
Writes certificate_t<t>.json.  The independent checker is verify_cert.py."""
import sys, time, json, pickle, random
from fractions import Fraction as Fr
import flint
from exact_space import *
from numeric import pieces
def amp_values(U, amp):
    mv = monomial_values(U); out = []
    for b, a in zip(BLOCK_LIST, amp):
        qs = []
        for row in a:
            z = CQ(0)
            for coef, mi in zip(row, b):
                if coef: z = z + CQ(coef) * mv[mi]
            qs.append(z)
        out.append(qs)
    return out
def var_index(amp):
    idx = {}; n = 0
    for bi, a in enumerate(amp):
        for i in range(len(a)):
            for j in range(i, len(a)): idx[(bi, i, j)] = n; n += 1
    return idx, n
if __name__ == '__main__':
    tq = Fr(sys.argv[1]); t = float(tq); t0 = time.time()
    amp = [[[Fr(q) for q in row] for row in a] for a in json.load(open('amplitudes_t%g.json' % t))['amplitudes']]
    Tnum = pickle.load(open('reduced2_t%g.pkl' % t, 'rb'))['T']
    vidx, N = var_index(amp)
    rng = random.Random(20261125); nU = 300; Us = [rational_unitary(rng) for _ in range(nU)]
    rowsA = []; rhs = []
    for U in Us:
        av = amp_values(U, amp); row = [Fr(0)] * N
        for bi, qs in enumerate(av):
            for i in range(len(qs)):
                for j in range(i, len(qs)):
                    v = qs[i].re * qs[j].re + qs[i].im * qs[j].im
                    row[vidx[(bi, i, j)]] = v if i == j else 2 * v
        g, a, b = pieces(U); rowsA.append(row); rhs.append(tq * g + (1 + tq) * a + tq * tq * b)
    print("exact system assembled: %d equations, %d unknowns  [%.0fs]" % (nU, N, time.time() - t0), flush=True)
    M = flint.fmpq_mat(nU, N + 1, [flint.fmpq(x.numerator, x.denominator) for row, r in zip(rowsA, rhs) for x in row + [r]])
    R, rank = M.rref()
    piv = []
    for r in range(rank):
        c = next(c for c in range(N + 1) if R[r, c] != 0); piv.append(c)
    assert N not in piv, "the exact system is inconsistent"
    free = [c for c in range(N) if c not in piv]
    print("exact rank %d ; free unknowns %d  [%.0fs]" % (rank, len(free), time.time() - t0), flush=True)
    xnum = [0.0] * N
    for (bi, i, j), k in vidx.items(): xnum[k] = float(Tnum[bi][i, j])
    x = [None] * N
    for c in free: x[c] = Fr(xnum[c]).limit_denominator(10**6)
    for r, c in enumerate(piv):
        val = Fr(int(R[r, N].p), int(R[r, N].q))
        for fc in free:
            e = R[r, fc]
            if e != 0: val -= Fr(int(e.p), int(e.q)) * x[fc]
        x[c] = val
    dev = max(abs(float(x[k]) - xnum[k]) for k in range(N))
    print("exact solution found ; max deviation from the numerical point %.2e  [%.0fs]" % (dev, time.time() - t0), flush=True)
    # exact positive definiteness
    T = []; minpm = None; ok = True
    for bi, a in enumerate(amp):
        r = len(a)
        if r == 0: T.append([]); continue
        Tb = [[x[vidx[(bi, min(i, j), max(i, j))]] for j in range(r)] for i in range(r)]; T.append(Tb)
        for k in range(1, r + 1):
            d = flint.fmpq_mat(k, k, [flint.fmpq(Tb[i][j].numerator, Tb[i][j].denominator) for i in range(k) for j in range(k)]).det()
            d = Fr(int(d.p), int(d.q))
            if d <= 0: ok = False
            minpm = d if minpm is None else min(minpm, d)
    print("all leading principal minors positive: %s ; smallest one %.3e" % (ok, float(minpm)))
    json.dump({'t': str(tq), 'statement': 't*g0 + (1+t)*D[12,12] + t^2*D[13,14] = sum_b (A_b m_b)^* T_b (A_b m_b) on U(4)',
               'amplitudes': [[[str(q) for q in row] for row in a] for a in amp], 'T': [[[str(q) for q in row] for row in Tb] for Tb in T]},
              open('certificate_t%g.json' % t, 'w'))
    print("certificate written: certificate_t%g.json  (positive definite: %s)" % (t, ok))
