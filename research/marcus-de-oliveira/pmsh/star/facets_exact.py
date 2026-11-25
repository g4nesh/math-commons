"""Exact facet enumeration of the projected polytope P_rho on every interval between critical values of rho.
Input: chiro_coeffs.npy (integer quadratics det_S(rho) for all 12-subsets S of the 24 projected vertices).
On an interval between consecutive critical values every det_S has constant sign, so the facets are constant.
Facet criterion: an 11-subset T spans a facet hyperplane iff the other 13 vertices lie weakly on one side.
Each facet is then classified:
  trivial : its vertex set lies in no exotic facet of Q  => it is a nonnegative combination of B >= 0, D >= 0 (LP duality);
  hard    : it lies in exactly one exotic facet F and equals Z(F) & Z(D_a) & Z(D_b); its inequality is
            g_F + alpha D_a + alpha' D_b >= 0, and we check (alpha-1)*alpha' = 1, alpha > 1, with (F,a,b) in the orbit of the
            base triple, so that it follows from  g_F + D_a + 2 sqrt(D_a D_b) >= 0."""
import itertools, sys, time, json, numpy as np, sympy as sp
from fractions import Fraction as Fr
from math import comb
from base import *
t0 = time.time()
C = np.load('chiro_coeffs.npy')
subs = np.array(list(itertools.combinations(range(24), 12)), dtype=np.int64)
binom = np.array([[comb(v, k) for k in range(13)] for v in range(24)], dtype=np.int64)
A = np.stack([binom[subs[:, k], k + 1] for k in range(11)], 1)            # element S_k sitting at T-position k
Bm = np.stack([binom[subs[:, k + 1], k + 1] for k in range(11)], 1)       # element S_{k+1} sitting at T-position k
cA = np.hstack([np.zeros((len(subs), 1), dtype=np.int64), np.cumsum(A, 1)]); cB = np.hstack([np.zeros((len(subs), 1), dtype=np.int64), np.cumsum(Bm, 1)])
Trank = np.stack([cA[:, i] + (cB[:, 11] - cB[:, i]) for i in range(12)], 1)   # rank of S minus its i-th element
N11 = comb(24, 11)
T_all = np.array(list(itertools.combinations(range(24), 11)), dtype=np.int64)
trk = sum(binom[T_all[:, k], k + 1] for k in range(11)); tmask = np.zeros(N11, dtype=np.int64); tmask[trk] = (1 << T_all).sum(1)
vbit = (1 << subs).astype(np.float64)
sgn_pos = np.array([(-1) ** i for i in range(12)])
print("tables built  [%.0fs]" % (time.time() - t0), flush=True)
# critical values
x = sp.symbols('rho'); prim = set()
for row in np.unique(C[np.abs(C).sum(1) > 0], axis=0):
    g = int(np.gcd.reduce(np.abs(row))); prim.add(tuple(int(v) // g for v in row))
crit = set()
for r in prim:
    poly = sp.Poly(r[0] + r[1] * x + r[2] * x**2, x)
    if poly.degree() > 0: crit |= set(sp.real_roots(poly))
crit = sorted(crit, key=lambda v: float(v)); cf = [float(v) for v in crit]
samples = [Fr(int(np.floor(cf[0])) - 1)] + [Fr((cf[i] + cf[i + 1]) / 2).limit_denominator(200) for i in range(len(cf) - 1)] + [Fr(int(np.ceil(cf[-1])) + 1)]
for i, s in enumerate(samples):                                             # exact check that samples separate the critical values
    if i > 0: assert sp.Rational(s.numerator, s.denominator) > crit[i - 1]
    if i < len(crit): assert sp.Rational(s.numerator, s.denominator) < crit[i]
print("%d critical values, %d intervals" % (len(crit), len(samples)), flush=True)
# zero sets
def zmask_of(w, w0=0): return int(sum(1 << v for v in range(24) if VX[v] @ w + w0 == 0))
EXd = exotic_facets(); EX = list(EXd.values()); ZF = [zmask_of(w) for w in EX]
ZB = {(s, j): zmask_of(np.eye(NX, dtype=int)[iB(s, j)]) for s in range(4) for j in range(4)}
ZD = {}
for S in P2:
    for T in P2:
        if 0 in S: ZD[(S, T)] = zmask_of(np.eye(NX, dtype=int)[iD(S, T)])   # Jacobi representative with row 1
# orbit of the base triple (g0; D[12,12], D[13,14]) as (Z(F), Z(a), Z(b)) with a = the minor carrying alpha > 1
g0 = np.zeros(NX, dtype=int); g0[iB(0,0)] += 1; g0[iB(2,2)] += 1; g0[iD((0,3),(2,3))] += 1; g0[iD((0,1),(0,3))] -= 1
ea = np.eye(NX, dtype=int)[iD((0,1),(0,1))]; eb = np.eye(NX, dtype=int)[iD((0,2),(0,3))]
ORB = set()
for al in PERMS:
    for be in PERMS:
        for tr in (False, True): ORB.add((zmask_of(act(g0, al, be, tr)), zmask_of(act(ea, al, be, tr)), zmask_of(act(eb, al, be, tr))))
print("orbit of the base inequality: %d members" % len(ORB), flush=True)
def shift_rat(k, q):
    eta = [Fr(0)] * 3; eta[k] = Fr(1); eta[2] -= 1; dx = [Fr(0)] * NX
    for a, (S, Sc) in enumerate(MATCH):
        for b, (T, Tc) in enumerate(MATCH):
            for S_ in (S, Sc):
                for T_ in (T, Tc): dx[iD(S_, T_)] += eta[a] * q[b] / 2
    return dx
popcount = lambda m: bin(m).count("1")
report = []; allok = True; triples_seen = {}
for I, rho in enumerate(samples):
    a_, b_ = rho.numerator, rho.denominator
    s = np.sign(C[:, 0] * b_ * b_ + C[:, 1] * a_ * b_ + C[:, 2] * a_ * a_)
    pos = np.zeros(N11); neg = np.zeros(N11); zm = np.zeros(N11)
    for i in range(12):
        sg = s * sgn_pos[i]
        pos += np.bincount(Trank[:, i], weights=(sg > 0), minlength=N11); neg += np.bincount(Trank[:, i], weights=(sg < 0), minlength=N11)
        zm += np.bincount(Trank[:, i], weights=(sg == 0) * vbit[:, i], minlength=N11)
    one = (pos > 0) ^ (neg > 0)
    facets = sorted(set((tmask[one] | zm[one].astype(np.int64)).tolist()))
    q = [rho, Fr(-1), 1 - rho]; sh = [shift_rat(0, q), shift_rat(1, q)]
    ntriv = nhard = 0; bad = []
    for Tm in facets:
        inF = [k for k in range(48) if Tm & ZF[k] == Tm]
        if not inF: ntriv += 1; continue
        Ds = [key for key, z in ZD.items() if Tm & z == Tm]; Bs = [key for key, z in ZB.items() if Tm & z == Tm]
        if len(inF) != 1 or len(Ds) != 2 or Bs or popcount(Tm) != 11 or (ZF[inF[0]] & ZD[Ds[0]] & ZD[Ds[1]]) != Tm:
            bad.append((Tm, len(inF), len(Ds), len(Bs))); continue
        w = EX[inF[0]]; ia, ib = iD(*Ds[0]), iD(*Ds[1])
        gs = [sum(int(w[k]) * sh[i][k] for k in range(NX)) for i in range(2)]
        M = [[sh[0][ia], sh[0][ib]], [sh[1][ia], sh[1][ib]]]; det = M[0][0] * M[1][1] - M[0][1] * M[1][0]
        if det == 0: bad.append((Tm, 'singular')); continue
        al = (-gs[0] * M[1][1] + gs[1] * M[0][1]) / det; alp = (-gs[1] * M[0][0] + gs[0] * M[1][0]) / det
        if al > 1 and alp > 0 and (al - 1) * alp == 1 and (ZF[inF[0]], ZD[Ds[0]], ZD[Ds[1]]) in ORB: key = (inF[0], Ds[0], Ds[1])
        elif alp > 1 and al > 0 and (alp - 1) * al == 1 and (ZF[inF[0]], ZD[Ds[1]], ZD[Ds[0]]) in ORB: key = (inF[0], Ds[1], Ds[0])
        else: bad.append((Tm, str(al), str(alp))); continue
        nhard += 1; triples_seen.setdefault(key, []).append(I)
    ok = (not bad) and len(facets) == 184 and nhard == 96
    allok &= ok
    lo = "-inf" if I == 0 else str(crit[I - 1]); hi = "+inf" if I == len(crit) else str(crit[I])
    report.append({'interval': [lo, hi], 'sample': str(rho), 'facets': len(facets), 'trivial': ntriv, 'hard': nhard, 'unexplained': len(bad)})
    print("interval %2d (%s, %s) sample %s: facets %d = %d trivial + %d hard (all from the single inequality) ; unexplained %d  [%.0fs]" %
          (I, lo, hi, rho, len(facets), ntriv, nhard, len(bad), time.time() - t0), flush=True)
print("ALL INTERVALS OK:", allok, "; distinct (F,a,b) triples used:", len(triples_seen))
json.dump({'critical_values': [str(c) for c in crit], 'intervals': report, 'all_ok': bool(allok), 'distinct_triples': len(triples_seen)}, open('facets_exact.json', 'w'), indent=1)
import pickle; pickle.dump(triples_seen, open('triples_seen.pkl', 'wb'))
