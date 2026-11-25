"""Two remaining exact checks.
(1) For each of the 288 triples (F, a, b) used as hard facets: solve the multipliers symbolically in rho and verify
    (alpha - 1) * alpha' == 1 identically, and that alpha' has no zero or pole strictly inside an interval where the triple is used.
(2) The compound polytope Q itself: enumerate its facets from the orientation determinants of all 14-subsets of its
    24 vertices (13 affine coordinates), and compare with 16 (B >= 0) + 18 (D >= 0) + 48 exotic."""
import itertools, pickle, json, time, numpy as np, sympy as sp
from math import comb
from base import *
t0 = time.time()
rho = sp.symbols('rho', real=True); q = [rho, -1, 1 - rho]
def shift_sym(k):
    eta = [0, 0, 0]; eta[k] = 1; eta[2] -= 1; dx = [sp.Integer(0)] * NX
    for a, (S, Sc) in enumerate(MATCH):
        for b, (T, Tc) in enumerate(MATCH):
            for S_ in (S, Sc):
                for T_ in (T, Tc): dx[iD(S_, T_)] += sp.Rational(eta[a], 2) * q[b]
    return dx
sh = [shift_sym(0), shift_sym(1)]
EX = list(exotic_facets().values())
triples = pickle.load(open('triples_seen.pkl', 'rb')); info = json.load(open('facets_exact.json'))
crit = [sp.sympify(c) for c in info['critical_values']]
bad = 0; forms = {}
for (F, Da, Db), ivals in triples.items():
    w = EX[F]; ia, ib = iD(*Da), iD(*Db)
    gs = sp.Matrix([sum(int(w[k]) * sh[i][k] for k in range(NX)) for i in range(2)])
    A = sp.Matrix([[sh[0][ia], sh[0][ib]], [sh[1][ia], sh[1][ib]]]); al = A.solve(-gs)
    al = [sp.cancel(v) for v in al]
    ok = sp.cancel((al[0] - 1) * al[1] - 1) == 0
    num, den = sp.fraction(sp.together(al[1])); special = set(sp.real_roots(sp.Poly(num, rho))) | set(sp.real_roots(sp.Poly(den, rho))) if (sp.Poly(num, rho).degree() > 0 or sp.Poly(den, rho).degree() > 0) else set()
    for I in ivals:
        lo = -sp.oo if I == 0 else crit[I - 1]; hi = sp.oo if I == len(crit) else crit[I]
        if any(lo < s < hi for s in special): ok = False
    bad += not ok; forms[(str(sp.factor(al[0])), str(sp.factor(al[1])))] = forms.get((str(sp.factor(al[0])), str(sp.factor(al[1]))), 0) + 1
print("(1) triples checked: %d ; failing the identity (alpha-1)*alpha' = 1 or with a sign change inside a used interval: %d  [%.0fs]" % (len(triples), bad, time.time() - t0))
print("    distinct multiplier pairs (alpha, alpha') and multiplicities:")
for k, v in sorted(forms.items(), key=lambda kv: -kv[1]): print("      alpha = %-22s alpha' = %-22s x %d" % (k[0], k[1], v))
# (2) facets of Q
Y = np.array([np.concatenate([[1], [v[iB(s, j)] for s in range(3) for j in range(3)], [n_matrix(v.astype(float))[0][a, b] for a in range(2) for b in range(2)]]) for v in VX]).astype(float)
subs = np.array(list(itertools.combinations(range(24), 14)), dtype=np.int64)
d = np.concatenate([np.linalg.det(Y[subs[lo:lo + 150000]]) for lo in range(0, len(subs), 150000)]); dr = np.rint(d)
print("(2) 14-subsets: %d ; max |float det - integer| = %.1e ; max |det| = %d" % (len(subs), np.abs(d - dr).max(), np.abs(dr).max()))
s = np.sign(dr); k = 13
binom = np.array([[comb(v, j) for j in range(15)] for v in range(24)], dtype=np.int64)
A = np.stack([binom[subs[:, j], j + 1] for j in range(k)], 1); Bm = np.stack([binom[subs[:, j + 1], j + 1] for j in range(k)], 1)
cA = np.hstack([np.zeros((len(subs), 1), dtype=np.int64), np.cumsum(A, 1)]); cB = np.hstack([np.zeros((len(subs), 1), dtype=np.int64), np.cumsum(Bm, 1)])
N13 = comb(24, 13); pos = np.zeros(N13); neg = np.zeros(N13); zm = np.zeros(N13)
for i in range(14):
    rk = cA[:, i] + (cB[:, k] - cB[:, i]); sg = s * (-1) ** i
    pos += np.bincount(rk, weights=(sg > 0), minlength=N13); neg += np.bincount(rk, weights=(sg < 0), minlength=N13)
    zm += np.bincount(rk, weights=(sg == 0) * (1 << subs[:, i]).astype(float), minlength=N13)
T_all = np.array(list(itertools.combinations(range(24), 13)), dtype=np.int64)
trk = sum(binom[T_all[:, j], j + 1] for j in range(13)); tmask = np.zeros(N13, dtype=np.int64); tmask[trk] = (1 << T_all).sum(1)
one = (pos > 0) ^ (neg > 0); facets = set((tmask[one] | zm[one].astype(np.int64)).tolist())
zmask_of = lambda w: int(sum(1 << v for v in range(24) if VX[v] @ w == 0))
ZF = {zmask_of(w) for w in EX}; ZB = {zmask_of(np.eye(NX, dtype=int)[iB(s_, j)]) for s_ in range(4) for j in range(4)}
ZD = {zmask_of(np.eye(NX, dtype=int)[iD(S, T)]) for S in P2 for T in P2}
print("    facets of Q found: %d ; equal to the union of %d (B=0) + %d (D=0) + %d (exotic): %s  [%.0fs]" % (len(facets), len(ZB), len(ZD), len(ZF), facets == (ZF | ZB | ZD), time.time() - t0))
