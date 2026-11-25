"""The end s -> 1 (t -> infinity).

At s = 1 the target is the single square F_1 = D[13,14] = |det U[13,14]|^2.  On U(4)
    D[13,14] = sum_{i,j} c_ij |det U[13,14] u_ij|^2     whenever  c_ij = alpha_i + beta_j,  sum alpha + sum beta = 1,
because rows and columns of U have norm 1.  Every product p_ij = det U[13,14] * u_ij lies in the amplitude space at s = 1
(checked exactly here), with coordinate vector w_ij.  So T_inf := sum c_ij w_ij w_ij^T is positive semidefinite by
construction when all c_ij >= 0, and represents F_1.

Certificate on [a, 1]:  T_b(s) = sum_{k<d} B_k(s) T_{b,k} + B_d(s) T_inf   with T_{b,k} positive definite for k < d.
For s in [a, 1) the weight B_0(s) is positive, so T_b(s) is positive definite; at s = 1 the statement is F_1 >= 0.

This module: exact coordinates w_ij (rational), and the numerical search (maximise the margin of T_{b,k}, k < d, over
the T_{b,k} and over alpha, beta).   Usage: python3 end1.py <deg> <a> <d>      Output: interval<deg>_<a>_1_d<d>.pkl"""
import sys, json, time, pickle
from fractions import Fraction as Fr
import numpy as np, cvxpy as cp, scipy.sparse as sps, flint
from exactn import space
from interval_sdp import Problem, bern
def fm(rows): return flint.fmpq_mat(len(rows), len(rows[0]), [flint.fmpq(Fr(x).numerator, Fr(x).denominator) for r in rows for x in r])
def square_coordinates(deg):
    """for each entry v = 4i+j: (block index, coordinates w of det U[13,14]*u_v in the rows of A_b(1) = A0 + A1), exact"""
    assert deg == 3
    fam = json.load(open('family%d.json' % deg)); MONS, blocks = space(deg); MIDX = {m: i for i, m in enumerate(MONS)}; where = {}; out = []
    for bi, blk in enumerate(blocks):
        for pos, mi in enumerate(blk): where[mi] = (bi, pos)
    for v in range(16):
        bi, p1 = where[MIDX[tuple(sorted((0, 11, v)))]]; bj, p2 = where[MIDX[tuple(sorted((3, 8, v)))]]; assert bi == bj      # u11 u34 u_v - u14 u31 u_v
        vec = [0] * len(blocks[bi]); vec[p1] += 1; vec[p2] -= 1
        A = [[x0 + x1 for x0, x1 in zip(r0, r1)] for r0, r1 in zip(fam['A0'][bi], fam['A1'][bi])]; r = len(A); n = len(vec)
        R, rk = fm([[A[i][j] for i in range(r)] + [vec[j]] for j in range(n)]).rref()      # solve w A = vec
        assert rk == r and all(R[i, i] == 1 for i in range(r)), "the square does not lie in the amplitude space at s = 1"
        w = [Fr(int(R[i, r].p), int(R[i, r].q)) for i in range(r)]
        assert all(sum(w[i] * A[i][j] for i in range(r)) == vec[j] for j in range(n))
        out.append((bi, w))
    return out
def gram_vector(P, bi, w):
    """the vector x (upper triangles of all blocks) of the rank-one block w w^T in block bi"""
    lo, hi, r = P.slices[bi]; x = np.zeros(P.N); iu = np.triu_indices(r); wf = np.array([float(q) for q in w]); x[lo:hi] = wf[iu[0]] * wf[iu[1]]; return x
if __name__ == '__main__':
    deg = int(sys.argv[1]); a = Fr(sys.argv[2]); d = int(sys.argv[3]); t0 = time.time(); af = float(a)
    P = Problem(deg); sq = square_coordinates(deg); Gv = np.array([gram_vector(P, bi, w) for bi, w in sq]).T      # N x 16
    nodes = [af + (1 - af) * (0.5 - 0.5 * np.cos(np.pi * (j + 0.5) / (d + 3))) for j in range(d + 3)]
    Lmat = np.zeros((16, 8))
    for v in range(16): Lmat[v, v // 4] = 1; Lmat[v, 4 + v % 4] = 1      # c = Lmat (alpha, beta)
    rows = []; rhss = []
    for sj in nodes:
        QA, rhs, out = P.conditions(sj); tau = (sj - af) / (1 - af)
        rows.append(np.hstack([bern(d, k, tau) * QA for k in range(d)] + [bern(d, d, tau) * (QA @ Gv @ Lmat)])); rhss.append(rhs)      # unknowns: X_0..X_{d-1}, alpha, beta
    rows.append(np.hstack([np.zeros((1, d * P.N)), np.ones((1, 8))])); rhss.append(np.array([1.0]))      # sum alpha + sum beta = 1
    Abig = np.vstack(rows); rbig = np.concatenate(rhss); u, sv, vt = np.linalg.svd(Abig, full_matrices=False); keep = sv > 1e-9 * sv[0]
    Ared = u[:, keep].T @ Abig; rred = u[:, keep].T @ rbig; incons = float(np.linalg.norm(rbig - u[:, keep] @ rred) / np.linalg.norm(rbig))
    Z = cp.Variable(d * P.N + 8); X = Z[:d * P.N]; ab = Z[d * P.N:]; lam = cp.Variable(); cons = [Ared @ Z == rred, Lmat @ ab >= 0]
    ones = [s[0] for s in P.slices if s is not None and s[2] == 1]
    for k in range(d):
        off = k * P.N
        if ones: cons.append(X[[off + o for o in ones]] >= lam)
        for s in P.slices:
            if s is None or s[2] == 1: continue
            lo, hi, r = s; S = sps.lil_matrix((r * r, hi - lo))
            for q, (i, j) in enumerate(zip(*np.triu_indices(r))): S[i * r + j, q] = 1; S[j * r + i, q] = 1
            cons.append(cp.reshape(sps.csr_matrix(S) @ X[off + lo:off + hi], (r, r), order='C') - lam * np.eye(r) >> 0)
    prob = cp.Problem(cp.Maximize(lam), cons); prob.solve(solver='CLARABEL')
    print("degree %d, s in [%s, 1] (t >= %s), Bernstein degree %d with the last coefficient fixed to the squares: conditions %d, inconsistency %.1e ; best margin lambda = %s (%s)  [%.0fs]" %
          (deg, a, a / (1 - a), d, int(keep.sum()), incons, None if lam.value is None else "%+.4e" % lam.value, prob.status, time.time() - t0), flush=True)
    if Z.value is not None:
        Xv = Z.value[:d * P.N]; abv = Z.value[d * P.N:]
        print("   alpha = %s ; beta = %s" % (np.round(abv[:4], 6).tolist(), np.round(abv[4:], 6).tolist()))
        pickle.dump({'deg': deg, 'a': str(a), 'd': d, 'x': [Xv[k * P.N:(k + 1) * P.N] for k in range(d)], 'alpha': abv[:4], 'beta': abv[4:], 'lam': float(lam.value)}, open('interval%d_%s_1_d%d.pkl' % (deg, str(a).replace('/', '-'), d), 'wb'))
