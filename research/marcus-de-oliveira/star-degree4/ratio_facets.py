"""Every facet of the projected polytope at ANY given rational ratios, each with an exact certificate that reduces it
to nonnegativity of coordinates or to one member of the family
        L_t = g0 + (1 + 1/t) D[12,12] + t D[13,14]        (t > 0;  t L_t = f_t).

Same method as harmonic_facets.py (which treats the harmonic ratios, where only t = 1 occurs).  For each facet the affine
functional h is computed exactly and one of the following exact identities on the 24 vertices is found:
      'coordinates':  h o y = sum_i lam_i x_i ,  lam_i >= 0 rational, x_i among B[s,j], D[S,T];
      'family'     :  h o y = c * L_t o gamma ,  c > 0, t > 0 rational, gamma a row/column permutation, possibly with
                      transposition.
Printed per ratio: the numbers of both kinds and the set of t that occur.  PMSH(4) for a spectrum with that ratio follows
once L_t >= 0 on U(4) is proved for each t in the set.

Usage: python3 ratio_facets.py 1/3 2/3 ...     Output: ratio_facets.json (certificates for the ratios given).
Needs ../pmsh/star/chiro_coeffs.npy (run `python3 chiro.py` in ../pmsh/star first; exactness: chiro_exact.py)."""
import itertools, json, sys, time, numpy as np
from fractions import Fraction as Fr
from math import comb
from scipy.optimize import linprog
sys.path.insert(0, '../pmsh/star')
from base import *
t0 = time.time(); C = np.load('../pmsh/star/chiro_coeffs.npy')
# ---- facet enumeration (as in ratio_half.py) ----
subs = np.array(list(itertools.combinations(range(24), 12)), dtype=np.int64)
binom = np.array([[comb(v, k) for k in range(13)] for v in range(24)], dtype=np.int64)
A = np.stack([binom[subs[:, k], k + 1] for k in range(11)], 1); Bm = np.stack([binom[subs[:, k + 1], k + 1] for k in range(11)], 1)
cA = np.hstack([np.zeros((len(subs), 1), dtype=np.int64), np.cumsum(A, 1)]); cB = np.hstack([np.zeros((len(subs), 1), dtype=np.int64), np.cumsum(Bm, 1)])
Trank = np.stack([cA[:, i] + (cB[:, 11] - cB[:, i]) for i in range(12)], 1); N11 = comb(24, 11)
T_all = np.array(list(itertools.combinations(range(24), 11)), dtype=np.int64)
trk = sum(binom[T_all[:, k], k + 1] for k in range(11)); tmask = np.zeros(N11, dtype=np.int64); tmask[trk] = (1 << T_all).sum(1)
vbit = (1 << subs).astype(np.float64); sgn_pos = np.array([(-1) ** i for i in range(12)])
def facets_at(rho):
    a_, b_ = rho.numerator, rho.denominator
    s = np.sign(C[:, 0] * b_ * b_ + C[:, 1] * a_ * b_ + C[:, 2] * a_ * a_)
    pos = np.zeros(N11); neg = np.zeros(N11); zm = np.zeros(N11)
    for i in range(12):
        sg = s * sgn_pos[i]
        pos += np.bincount(Trank[:, i], weights=(sg > 0), minlength=N11); neg += np.bincount(Trank[:, i], weights=(sg < 0), minlength=N11)
        zm += np.bincount(Trank[:, i], weights=(sg == 0) * vbit[:, i], minlength=N11)
    one = (pos > 0) ^ (neg > 0)
    return sorted(set((tmask[one] | zm[one].astype(np.int64)).tolist()))
# ---- exact linear algebra over Q ----
def rref(M):
    M = [row[:] for row in M]; piv = []; r = 0
    for c in range(len(M[0])):
        pr = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if pr is None: continue
        M[r], M[pr] = M[pr], M[r]; d = M[r][c]; M[r] = [v / d for v in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]; M[i] = [u - f * v for u, v in zip(M[i], M[r])]
        piv.append(c); r += 1
        if r == len(M): break
    return M, piv
def kernel(M):
    R, piv = rref(M); n = len(M[0]); free = [c for c in range(n) if c not in piv]; out = []
    for f in free:
        v = [Fr(0)] * n; v[f] = Fr(1)
        for i, c in enumerate(piv): v[c] = -R[i][f]
        out.append(v)
    return out
def solve_unique(cols, rhs):
    """exact solution of sum_k z_k cols[k] = rhs if it exists and is unique, else None"""
    m = len(rhs); k = len(cols); R, piv = rref([[cols[j][i] for j in range(k)] + [rhs[i]] for i in range(m)])
    if piv != list(range(k)): return None
    return [R[i][k] for i in range(k)]
# ---- data ----
NAMES = []; COLS = []
for s in range(4):
    for j in range(4): NAMES.append("B[%d,%d]" % (s + 1, j + 1)); COLS.append([Fr(int(v[iB(s, j)])) for v in VX])
for S in P2:
    if 0 not in S: continue
    for T in P2: NAMES.append("D[%d%d,%d%d]" % (S[0] + 1, S[1] + 1, T[0] + 1, T[1] + 1)); COLS.append([Fr(int(v[iD(S, T)])) for v in VX])
VXf = np.array([[float(c[v]) for c in COLS] for v in range(24)])
g0 = np.zeros(NX, dtype=int); g0[iB(0, 0)] += 1; g0[iB(2, 2)] += 1; g0[iD((0, 3), (2, 3))] += 1; g0[iD((0, 1), (0, 3))] -= 1
ea = np.zeros(NX, dtype=int); ea[iD((0, 1), (0, 1))] = 1; eb = np.zeros(NX, dtype=int); eb[iD((0, 2), (0, 3))] = 1      # D[12,12], D[13,14] (1-based)
SYM = {}
for al in PERMS:
    for be in PERMS:
        for tr in (False, True):
            yg, ya, yb = [tuple(int(c) for c in VX @ act(w, al, be, tr)) for w in (g0, ea, eb)]
            SYM.setdefault(sum(1 << v for v in range(24) if yg[v] == 0 and ya[v] == 0 and yb[v] == 0), []).append((al, be, tr, yg, ya, yb))
Nn = [[[int(v[iD(S, T)] + v[iD(S, Tc)]) for (T, Tc) in MATCH] for (S, Sc) in MATCH] for v in VX]
rk = len(rref([[Fr(int(c)) for c in v] + [Fr(1)] for v in VX])[1])
print("affine dimension of the 24 vertices x(P_pi): %d (exact)" % (rk - 1), flush=True)
out = {}; allok = True
for spec in sys.argv[1:]:
    rho = Fr(spec); Y = [[Fr(1)] + [Fr(int(VX[v][iB(s, j)])) for s in range(3) for j in range(3)] + [Nn[v][0][0] + rho * Nn[v][0][1], Nn[v][1][0] + rho * Nn[v][1][1]] for v in range(24)]
    assert len(rref(Y)[1]) == 12       # P_rho is full-dimensional in R^11
    facets = facets_at(rho); recs = []; ncoord = 0; fam = {}; fail = 0
    for Tm in facets:
        verts = [v for v in range(24) if Tm >> v & 1]; ker = kernel([Y[v] for v in verts]); assert len(ker) == 1
        h = ker[0]; H = [sum(a * b for a, b in zip(Y[v], h)) for v in range(24)]
        if any(z < 0 for z in H): h = [-z for z in h]; H = [-z for z in H]
        assert all(z >= 0 for z in H) and [v for v in range(24) if H[v] == 0] == verts      # h is the facet functional
        rec = {'vertices': verts, 'h': [str(z) for z in h]}; done = False
        for al, be, tr, yg, ya, yb in SYM.get(Tm, []):
            sol = solve_unique([[Fr(z) for z in yg], [Fr(z) for z in ya], [Fr(z) for z in yb]], H)      # H = c1 yg + c2 ya + c3 yb
            if sol is None: continue
            c1, c2, c3 = sol
            if c1 > 0 and c3 > 0 and c2 * c3 == c1 * (c1 + c3):      # c2/c1 = 1 + 1/t with t = c3/c1
                t = c3 / c1; rec.update(kind='family', t=str(t), alpha=list(al), beta=list(be), transpose=tr, scale=str(c1)); fam[t] = fam.get(t, 0) + 1; done = True; break
        if not done:
            res = linprog(np.ones(len(COLS)), A_eq=VXf, b_eq=[float(z) for z in H], bounds=(0, None), method='highs-ds')
            if res.status == 0:
                supp = [i for i in range(len(COLS)) if res.x[i] > 1e-9]; lam = solve_unique([COLS[i] for i in supp], H)
                if lam is not None and all(z >= 0 for z in lam) and all(sum(l * COLS[i][v] for l, i in zip(lam, supp)) == H[v] for v in range(24)):
                    rec.update(kind='coordinates', combination={NAMES[i]: str(l) for l, i in zip(lam, supp) if l != 0}); ncoord += 1; done = True
        if not done: fail += 1; rec['kind'] = 'UNCERTIFIED'
        recs.append(rec)
    allok &= (fail == 0); out[str(rho)] = recs
    print("rho = %s: %d facets ; nonnegative coordinates: %d ; family L_t: %d with t in %s ; uncertified: %d  [%.0fs]" %
          (rho, len(facets), ncoord, sum(fam.values()), {str(k): v for k, v in sorted(fam.items())}, fail, time.time() - t0), flush=True)
json.dump(out, open('ratio_facets.json', 'w'), separators=(',', ':'))
print("EVERY FACET AT THESE RATIOS IS REDUCED EXACTLY TO COORDINATES OR TO THE FAMILY" if allok else "SOME FACET IS NOT CERTIFIED")
