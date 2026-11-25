"""Numerical search for a certificate valid on a whole interval of s = t/(1+t).

Target:  F_s = (1-s) D[12,12] + s(1-s) g0 + s^2 D[13,14]  ( = (1-s)^2 f_t ),  for every s in [a, b].
Ansatz:  F_s = sum_b (A_b(s) m_b)^* T_b(s) (A_b(s) m_b)  on U(4),
         A_b(s) = A0 + s A1 from family<deg>.json,   T_b(s) = sum_k B_k(s) T_{b,k}  (Bernstein basis of degree d on [a, b]).
If every T_{b,k} >= lambda I with lambda > 0 then T_b(s) >= lambda I on [a, b].
Both sides are polynomials of degree <= d+2 in s, so the identity is imposed at d+3 nodes (on a Haar sample of unitaries).

Usage: python3 interval_sdp.py <deg> <a> <b> <d>      Output: interval<deg>_<a>_<b>_d<d>.pkl (floating point)."""
import sys, json, time, pickle, itertools
from fractions import Fraction as Fr
from math import comb
import numpy as np, cvxpy as cp, scipy.sparse as sps
def setup(deg):
    MONS = list(itertools.combinations_with_replacement(range(16), deg)); MARR = np.array(MONS); blocks = {}
    for i, m in enumerate(MONS):
        r = [0] * 4; c = [0] * 4
        for v in m: r[v // 4] += 1; c[v % 4] += 1
        blocks.setdefault((tuple(r), tuple(c)), []).append(i)
    return MARR, list(blocks.values())
def haar(rng, K):
    z = (rng.standard_normal((K, 4, 4)) + 1j * rng.standard_normal((K, 4, 4))) / np.sqrt(2)
    q, r = np.linalg.qr(z); d = np.diagonal(r, axis1=1, axis2=2); return q * (d / np.abs(d))[:, None, :]
def monomials(Us, MARR):
    flat = Us.reshape(len(Us), 16); out = flat[:, MARR[:, 0]].copy()
    for k in range(1, MARR.shape[1]): out *= flat[:, MARR[:, k]]
    return out
def parts(Us):
    B = np.abs(Us)**2
    def D(S, T): return np.abs(Us[:, S[0], T[0]] * Us[:, S[1], T[1]] - Us[:, S[0], T[1]] * Us[:, S[1], T[0]])**2
    return B[:, 0, 0] + B[:, 2, 2] + D((0, 3), (2, 3)) - D((0, 1), (0, 3)), D((0, 1), (0, 1)), D((0, 2), (0, 3))
def F(s, g0, a, b): return (1 - s) * a + s * (1 - s) * g0 + s * s * b
def bern(d, k, tau): return comb(d, k) * tau**k * (1 - tau)**(d - k)
class Problem:
    def __init__(self, deg, K=None, seed=0):
        fam = json.load(open('family%d.json' % deg)); self.A0 = [np.array(a, dtype=float) if a else None for a in fam['A0']]; self.A1 = [np.array(a, dtype=float) if a else None for a in fam['A1']]
        self.MARR, self.blocks = setup(deg); rng = np.random.default_rng(seed); self.Us = haar(rng, K or (2500 if deg == 3 else 3000))
        self.M = monomials(self.Us, self.MARR); self.g0, self.a, self.b = parts(self.Us); self.slices = []; pos = 0
        for A in self.A0:
            if A is None: self.slices.append(None); continue
            r = A.shape[0]; self.slices.append((pos, pos + r * (r + 1) // 2, r)); pos += r * (r + 1) // 2
        self.N = pos
    def conditions(self, s):
        """(Qr A(s), Qr F_s): the identity at s on the sample, reduced to independent conditions"""
        cols = []
        for blk, A0, A1 in zip(self.blocks, self.A0, self.A1):
            if A0 is None: continue
            Q = self.M[:, blk] @ (A0 + s * A1).T; r = A0.shape[0]; iu = np.triu_indices(r)
            v = np.real(np.conj(Q[:, iu[0]]) * Q[:, iu[1]]); v[:, iu[0] != iu[1]] *= 2; cols.append(v)
        A = np.hstack(cols); u, sv, _ = np.linalg.svd(A, full_matrices=False); keep = sv > 1e-9 * sv[0]; Qr = u[:, keep].T
        f = F(s, self.g0, self.a, self.b); rhs = Qr @ f
        return Qr @ A, rhs, float(np.linalg.norm(f - Qr.T @ rhs) / np.linalg.norm(f))
    def solve(self, a, b, d, verbose=False):
        nodes = [a + (b - a) * (0.5 - 0.5 * np.cos(np.pi * (j + 0.5) / (d + 3))) for j in range(d + 3)]      # Chebyshev points of [a, b]
        rows = []; rhss = []; outside = 0
        for sj in nodes:
            QA, rhs, out = self.conditions(sj); outside = max(outside, out); tau = (sj - a) / (b - a) if b > a else 0.0
            rows.append(np.hstack([bern(d, k, tau) * QA for k in range(d + 1)])); rhss.append(rhs)
        # conditions at different nodes are partly redundant: keep an independent set
        Abig = np.vstack(rows); rbig = np.concatenate(rhss); u, sv, vt = np.linalg.svd(Abig, full_matrices=False); keep = sv > 1e-9 * sv[0]
        Ared = (u[:, keep].T @ Abig); rred = u[:, keep].T @ rbig; incons = float(np.linalg.norm(rbig - u[:, keep] @ rred) / np.linalg.norm(rbig))
        X = cp.Variable((d + 1) * self.N); lam = cp.Variable(); cons = [Ared @ X == rred]
        ones = [s[0] for s in self.slices if s is not None and s[2] == 1]
        for k in range(d + 1):
            off = k * self.N
            if ones: cons.append(X[[off + o for o in ones]] >= lam)
            for s in self.slices:
                if s is None or s[2] == 1: continue
                lo, hi, r = s; S = sps.lil_matrix((r * r, hi - lo))
                for q, (i, j) in enumerate(zip(*np.triu_indices(r))): S[i * r + j, q] = 1; S[j * r + i, q] = 1
                cons.append(cp.reshape(sps.csr_matrix(S) @ X[off + lo:off + hi], (r, r), order='C') - lam * np.eye(r) >> 0)
        prob = cp.Problem(cp.Maximize(lam), cons); prob.solve(solver='CLARABEL', verbose=verbose)
        xs = [None] * (d + 1) if X.value is None else [X.value[k * self.N:(k + 1) * self.N] for k in range(d + 1)]
        return (None if lam.value is None else float(lam.value)), prob.status, xs, max(outside, incons), int(keep.sum())
if __name__ == '__main__':
    deg = int(sys.argv[1]); a, b = Fr(sys.argv[2]), Fr(sys.argv[3]); d = int(sys.argv[4]); t0 = time.time()
    P = Problem(deg); lam, st, xs, outside, ncond = P.solve(float(a), float(b), d)
    print("degree %d, s in [%s, %s] (t in [%s, %s]), Bernstein degree %d: unknowns %d, conditions %d, F outside the range %.1e ; best common margin lambda = %s (%s)  [%.0fs]" %
          (deg, a, b, a / (1 - a), b / (1 - b), d, (d + 1) * P.N, ncond, outside, None if lam is None else "%+.4e" % lam, st, time.time() - t0), flush=True)
    if xs[0] is not None: pickle.dump({'deg': deg, 'a': str(a), 'b': str(b), 'd': d, 'x': xs, 'lam': lam, 'slices': P.slices}, open('interval%d_%s_%s_d%d.pkl' % (deg, str(a).replace('/', '-'), str(b).replace('/', '-'), d), 'wb'))
