"""Every facet of the projected polytope at a harmonic ratio, each with its own exact certificate of validity on U(4).

Setting (../pmsh/equivalence.md).  Projected point of x = (B, D):  y(x) = (B[s,j] for s, j <= 3 ;  n_11 + rho n_12 ;
n_21 + rho n_22)  in R^11, with n_ab = D[S_a,T_b] + D[S_a,T_b^c].  P_rho = conv of the 24 points y(x(P_pi)); it is
full-dimensional.  For each ratio (default 1/2, 2, -1) this script

 1. enumerates the facets of P_rho from the orientation quadratics (shown to be exact by chiro_exact.py):
    an 11-subset of the 24 points spans a facet iff the other 13 lie weakly on one side, not all on the hyperplane;
 2. computes each facet's affine functional h exactly (h = 0 on the facet, h >= 0 on all 24 points);
 3. certifies  h(y(x(U))) >= 0 for every unitary U  in one of two ways, checked in exact rational arithmetic on the
    24 vertices (two affine functions of x that agree on the 24 vertices agree on aff(Q), which contains every x(U)):
      'coordinates':  h o y = sum_i lam_i x_i  with lam_i >= 0 and x_i among B[s,j], D[S,T]      (then use B, D >= 0);
      'theorem S4' :  h o y = c * (g0 + 2 D[12,12] + D[13,14]) o gamma,  c > 0, gamma a row/column permutation, possibly
                      with transposition                                 (then use Theorem S4 at the permuted unitary).

This does not use the facet list of the compound polytope Q.  Output: harmonic_facets.json.
Needs ../pmsh/star/chiro_coeffs.npy (run `python3 chiro.py` in ../pmsh/star first)."""
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
f1 = np.zeros(NX, dtype=int); f1[iB(0, 0)] += 1; f1[iB(2, 2)] += 1; f1[iD((0, 3), (2, 3))] += 1; f1[iD((0, 1), (0, 3))] -= 1
f1[iD((0, 1), (0, 1))] += 2; f1[iD((0, 2), (0, 3))] += 1      # g0 + 2 D[12,12] + D[13,14], 1-based
SYM = {}
for al in PERMS:
    for be in PERMS:
        for tr in (False, True):
            yv = tuple(int(c) for c in VX @ act(f1, al, be, tr)); SYM.setdefault(sum(1 << v for v in range(24) if yv[v] == 0), []).append((al, be, tr, yv))
Nn = [[[int(v[iD(S, T)] + v[iD(S, Tc)]) for (T, Tc) in MATCH] for (S, Sc) in MATCH] for v in VX]
rk = len(rref([[Fr(int(c)) for c in v] + [Fr(1)] for v in VX])[1])
print("affine dimension of the 24 vertices x(P_pi): %d (exact)" % (rk - 1), flush=True)
# ---- sanity of the symmetry action (floating point, not part of the proof) ----
rng = np.random.default_rng(3); worst = 0
for _ in range(20):
    U = haar(rng); x = x_U(U)
    for al, be, tr in [(PERMS[5], PERMS[17], False), (PERMS[11], PERMS[2], True), (PERMS[23], PERMS[9], True)]:
        V = np.array([[U[be[j], al[s]] if tr else U[al[s], be[j]] for j in range(4)] for s in range(4)])
        worst = max(worst, abs(act(f1, al, be, tr) @ x - f1 @ x_U(V)))
print("action check: (f o gamma)(x(U)) = f(x(U')) with U'[s,j] = U[alpha s, beta j] (or U[beta j, alpha s] with transposition), max error %.1e" % worst, flush=True)
out = {}; allok = True
for spec in (sys.argv[1:] or ['1/2', '2', '-1']):
    rho = Fr(spec); Y = [[Fr(1)] + [Fr(int(VX[v][iB(s, j)])) for s in range(3) for j in range(3)] + [Nn[v][0][0] + rho * Nn[v][0][1], Nn[v][1][0] + rho * Nn[v][1][1]] for v in range(24)]
    assert len(rref(Y)[1]) == 12       # P_rho is full-dimensional in R^11
    facets = facets_at(rho); recs = []; counts = {'coordinates': 0, 'theorem S4': 0}; fail = 0
    for Tm in facets:
        verts = [v for v in range(24) if Tm >> v & 1]; ker = kernel([Y[v] for v in verts]); assert len(ker) == 1
        h = ker[0]; H = [sum(a * b for a, b in zip(Y[v], h)) for v in range(24)]
        if any(z < 0 for z in H): h = [-z for z in h]; H = [-z for z in H]
        assert all(z >= 0 for z in H) and [v for v in range(24) if H[v] == 0] == verts      # h is the facet functional
        rec = {'vertices': verts, 'h': [str(z) for z in h]}; done = False
        for al, be, tr, yv in SYM.get(Tm, []):
            v0 = next(v for v in range(24) if H[v] != 0); c = H[v0] / yv[v0]
            if c > 0 and all(H[v] == c * yv[v] for v in range(24)):
                rec.update(kind='theorem S4', alpha=list(al), beta=list(be), transpose=tr, scale=str(c)); done = True; break
        if not done:
            res = linprog(np.ones(len(COLS)), A_eq=VXf, b_eq=[float(z) for z in H], bounds=(0, None), method='highs-ds')
            if res.status == 0:
                supp = [i for i in range(len(COLS)) if res.x[i] > 1e-9]; lam = solve_unique([COLS[i] for i in supp], H)
                if lam is not None and all(z >= 0 for z in lam) and all(sum(l * COLS[i][v] for l, i in zip(lam, supp)) == H[v] for v in range(24)):
                    rec.update(kind='coordinates', combination={NAMES[i]: str(l) for l, i in zip(lam, supp) if l != 0}); done = True
        if done: counts[rec['kind']] += 1
        else: fail += 1; rec['kind'] = 'UNCERTIFIED'
        recs.append(rec)
    allok &= (fail == 0); out[str(rho)] = recs
    print("rho = %s: %d facets ; certified by nonnegative coordinates: %d ; by Theorem S4 at a permuted unitary: %d ; uncertified: %d  [%.0fs]" %
          (rho, len(facets), counts['coordinates'], counts['theorem S4'], fail, time.time() - t0), flush=True)
json.dump(out, open('harmonic_facets.json', 'w'), separators=(',', ':'))
print("EVERY FACET AT THESE RATIOS IS CERTIFIED EXACTLY" if allok else "SOME FACET IS NOT CERTIFIED")
