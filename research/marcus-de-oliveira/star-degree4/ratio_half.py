"""Exact facets of the projected polytope P_rho AT given rational ratios (default 1/2, 2, -1: the harmonic ratios).
These are critical values of rho, where the general enumeration of ../pmsh/star/facets_exact.py (done at interior
points of the 42 intervals) does not apply directly.  Same method, evaluated at the exact ratio:
an 11-subset of the 24 projected vertices spans a facet iff the other 13 lie weakly on one side.
Each facet is classified:
  trivial: contained in no exotic facet of Q  => nonnegative combination of B >= 0, D >= 0 (LP duality);
  hard   : contained in exactly one exotic facet F, equal to Z(F) & Z(D_a) & Z(D_b), with inequality
           g_F + alpha D_a + alpha' D_b >= 0; we print the multiplier pairs that occur.
Needs ../pmsh/star/chiro_coeffs.npy (run `python3 chiro.py` in ../pmsh/star first; 8 seconds)."""
import itertools, sys, time, numpy as np
from fractions import Fraction as Fr
from math import comb
sys.path.insert(0, '../pmsh/star')
from base import *
t0 = time.time(); C = np.load('../pmsh/star/chiro_coeffs.npy')
subs = np.array(list(itertools.combinations(range(24), 12)), dtype=np.int64)
binom = np.array([[comb(v, k) for k in range(13)] for v in range(24)], dtype=np.int64)
A = np.stack([binom[subs[:, k], k + 1] for k in range(11)], 1); Bm = np.stack([binom[subs[:, k + 1], k + 1] for k in range(11)], 1)
cA = np.hstack([np.zeros((len(subs), 1), dtype=np.int64), np.cumsum(A, 1)]); cB = np.hstack([np.zeros((len(subs), 1), dtype=np.int64), np.cumsum(Bm, 1)])
Trank = np.stack([cA[:, i] + (cB[:, 11] - cB[:, i]) for i in range(12)], 1); N11 = comb(24, 11)
T_all = np.array(list(itertools.combinations(range(24), 11)), dtype=np.int64)
trk = sum(binom[T_all[:, k], k + 1] for k in range(11)); tmask = np.zeros(N11, dtype=np.int64); tmask[trk] = (1 << T_all).sum(1)
vbit = (1 << subs).astype(np.float64); sgn_pos = np.array([(-1) ** i for i in range(12)])
def zmask_of(w): return int(sum(1 << v for v in range(24) if VX[v] @ w == 0))
EX = list(exotic_facets().values()); ZF = [zmask_of(w) for w in EX]
ZB = {(s, j): zmask_of(np.eye(NX, dtype=int)[iB(s, j)]) for s in range(4) for j in range(4)}
ZD = {(S, T): zmask_of(np.eye(NX, dtype=int)[iD(S, T)]) for S in P2 for T in P2 if 0 in S}
g0 = np.zeros(NX, dtype=int); g0[iB(0,0)] += 1; g0[iB(2,2)] += 1; g0[iD((0,3),(2,3))] += 1; g0[iD((0,1),(0,3))] -= 1
ea = np.eye(NX, dtype=int)[iD((0,1),(0,1))]; eb = np.eye(NX, dtype=int)[iD((0,2),(0,3))]; ORB = set()
for al in PERMS:
    for be in PERMS:
        for tr in (False, True): ORB.add((zmask_of(act(g0, al, be, tr)), zmask_of(act(ea, al, be, tr)), zmask_of(act(eb, al, be, tr))))
def shift_rat(k, q):
    eta = [Fr(0)] * 3; eta[k] = Fr(1); eta[2] -= 1; dx = [Fr(0)] * NX
    for a, (S, Sc) in enumerate(MATCH):
        for b, (T, Tc) in enumerate(MATCH):
            for S_ in (S, Sc):
                for T_ in (T, Tc): dx[iD(S_, T_)] += eta[a] * q[b] / 2
    return dx
popcount = lambda m: bin(m).count("1")
for spec in (sys.argv[1:] or ['1/2', '2', '-1']):
    rho = Fr(spec); a_, b_ = rho.numerator, rho.denominator
    s = np.sign(C[:, 0] * b_ * b_ + C[:, 1] * a_ * b_ + C[:, 2] * a_ * a_)
    pos = np.zeros(N11); neg = np.zeros(N11); zm = np.zeros(N11)
    for i in range(12):
        sg = s * sgn_pos[i]
        pos += np.bincount(Trank[:, i], weights=(sg > 0), minlength=N11); neg += np.bincount(Trank[:, i], weights=(sg < 0), minlength=N11)
        zm += np.bincount(Trank[:, i], weights=(sg == 0) * vbit[:, i], minlength=N11)
    one = (pos > 0) ^ (neg > 0); facets = sorted(set((tmask[one] | zm[one].astype(np.int64)).tolist()))
    q = [rho, Fr(-1), 1 - rho]; sh = [shift_rat(0, q), shift_rat(1, q)]
    ntriv = 0; hard = {}; bad = []; sizes = {}
    for Tm in facets:
        sizes[popcount(Tm)] = sizes.get(popcount(Tm), 0) + 1
        inF = [k for k in range(48) if Tm & ZF[k] == Tm]
        if not inF: ntriv += 1; continue
        Ds = [key for key, z in ZD.items() if Tm & z == Tm]; Bs = [key for key, z in ZB.items() if Tm & z == Tm]
        if len(inF) != 1 or len(Ds) != 2 or Bs or (ZF[inF[0]] & ZD[Ds[0]] & ZD[Ds[1]]) != Tm: bad.append((popcount(Tm), len(inF), len(Ds), len(Bs))); continue
        w = EX[inF[0]]; ia, ib = iD(*Ds[0]), iD(*Ds[1])
        gs = [sum(int(w[k]) * sh[i][k] for k in range(NX)) for i in range(2)]
        M = [[sh[0][ia], sh[0][ib]], [sh[1][ia], sh[1][ib]]]; det = M[0][0] * M[1][1] - M[0][1] * M[1][0]
        if det == 0: bad.append((popcount(Tm), 'singular')); continue
        al = (-gs[0] * M[1][1] + gs[1] * M[0][1]) / det; alp = (-gs[1] * M[0][0] + gs[0] * M[1][0]) / det
        if (ZF[inF[0]], ZD[Ds[0]], ZD[Ds[1]]) in ORB and al > 1: key = (str(al), str(alp))
        elif (ZF[inF[0]], ZD[Ds[1]], ZD[Ds[0]]) in ORB and alp > 1: key = (str(alp), str(al))
        else: bad.append((popcount(Tm), str(al), str(alp))); continue
        hard[key] = hard.get(key, 0) + 1
    print("rho = %s: %d facets (vertex counts %s) = %d trivial + %d hard ; hard multiplier pairs (alpha on D_a, alpha' on D_b): %s ; unexplained: %s  [%.0fs]" %
          (rho, len(facets), dict(sorted(sizes.items())), ntriv, sum(hard.values()), hard, bad[:5] if bad else 0, time.time() - t0), flush=True)
