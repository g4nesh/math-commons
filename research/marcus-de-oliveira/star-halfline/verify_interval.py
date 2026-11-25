"""Independent exact check of an interval certificate interval<deg>_<a>_<b>.json.

Claim.  For every s in [a, b] and every U in U(4):
    F_s(U) := (1-s) D[12,12] + s(1-s) g0 + s^2 D[13,14]  =  sum_b (A_b(s) m_b(U))^* T_b(s) (A_b(s) m_b(U))   and   T_b(s) > 0,
hence F_s >= 0 on U(4).   Here A_b(s) = A0_b + s A1_b (integer matrices), T_b(s) = sum_k B_k(s) T_{b,k} with B_k the
Bernstein polynomials of degree d on [a, b], g0 = B11 + B33 + D[14,34] - D[12,14] and D[S,T] = |det U[S,T]|^2 (1-based).
With t = s/(1-s):  F_s = (1-s)^2 (t g0 + (1+t) D[12,12] + t^2 D[13,14]).

Checks, all in exact arithmetic:
 (P) every T_{b,k} is symmetric with positive leading principal minors.  Since B_k >= 0 on [a, b] and sum_k B_k = 1,
     T_b(s) is positive definite for every s in [a, b].
 (I) the identity holds at d+3 nodes s chosen HERE (from the seed; not the nodes used in the construction) and, at each
     node, at N fresh rational unitaries.
 (E) the evaluation functionals at those unitaries have rank r_E on the space W of conjugation-symmetric torus-invariant
     forms of bidegree (deg, deg)  (computed modulo a prime: a lower bound over Q), r_E = 261 for deg 3, 954 for deg 4.
Why this proves the claim.  The relation space K (unitarity relations times monomial pairs) lies in the kernel of every
evaluation and has rank >= dim W - r_E (deg 3: exact rank 1227 of 1488, ../star-exact/exact_space.py; deg 4: rank >= 10254
of 11208, ../star-degree4/lemma4.py).  So an element of W that vanishes at the unitaries of (E) vanishes on U(4).  For a
fixed node s, both sides of the identity restrict on U(4) to elements of W; by (I) their difference vanishes on U(4).
For a fixed U in U(4) both sides are polynomials in s of degree <= d+2 that agree at d+3 nodes; so they agree for all s."""
import sys, json, random, time
from fractions import Fraction as Fr
import numpy as np, flint
from exactn import *
t0 = time.time(); path = sys.argv[1]; seed = int(sys.argv[2]) if len(sys.argv) > 2 else 777
cert = json.load(open(path)); deg = cert['degree']; a, b = Fr(cert['a']), Fr(cert['b']); d = cert['d']; A0, A1 = cert['A0'], cert['A1']
T = [[[[Fr(q) for q in row] for row in Tb] for Tb in Tk] for Tk in cert['T']]; MONS, blocks = space(deg)
nU = int(sys.argv[3]) if len(sys.argv) > 3 else (300 if deg == 3 else 1100); rE_expected = {3: 261, 4: 954}[deg]
assert a < b and len(T) == d + 1 and len(A0) == len(A1) == len(blocks) and all(len(Tk) == len(blocks) for Tk in T)
for a0, a1, blk in zip(A0, A1, blocks):
    assert len(a0) == len(a1) and all(len(r) == len(blk) and all(isinstance(c, int) for c in r) for r in a0 + a1)
# (P)
okP = True; minpm = None; namp = sum(len(a0) for a0 in A0)
for Tk in T:
    for a0, Tb in zip(A0, Tk):
        r = len(a0); assert len(Tb) == r and all(len(row) == r for row in Tb)
        okP &= all(Tb[i][j] == Tb[j][i] for i in range(r) for j in range(r))
        for m in range(1, r + 1):
            dt = flint.fmpq_mat(m, m, [flint.fmpq(Tb[i][j].numerator, Tb[i][j].denominator) for i in range(m) for j in range(m)]).det(); dt = Fr(int(dt.p), int(dt.q))
            okP &= dt > 0; minpm = dt if minpm is None else min(minpm, dt)
print("(P) degree %d, s in [%s, %s], i.e. t in [%s, %s] ; %d amplitudes ; %d Bernstein coefficients per block ; all symmetric with positive leading principal minors: %s (smallest %.3e)  [%.0fs]"
      % (deg, a, b, a / (1 - a), b / (1 - b), namp, d + 1, okP, float(minpm), time.time() - t0), flush=True)
# (I)
rng = random.Random(seed); Us = [rational_unitary(rng) for _ in range(nU)]; assert all(is_unitary(U) for U in Us)
nodes = set()
while len(nodes) < d + 3: nodes.add(Fr(rng.randint(1, 96), 97))
nodes = sorted(nodes); ints = [integer_form(U) for U in Us]; mons = [monomials_int(X, Y, MONS) for X, Y, L in ints]; bad = 0
for s in nodes:
    tau = (s - a) / (b - a); w = [bernstein(d, k, tau) for k in range(d + 1)]; assert sum(w) == 1; sn, sd = s.numerator, s.denominator
    amp = [[[sd * x0 + sn * x1 for x0, x1 in zip(r0, r1)] for r0, r1 in zip(a0, a1)] for a0, a1 in zip(A0, A1)]      # sd * A(s), integer
    Ts = []
    for bi in range(len(blocks)):
        r = len(A0[bi])
        if r == 0: Ts.append(None); continue
        M = [[sum(w[k] * T[k][bi][i][j] for k in range(d + 1)) for j in range(r)] for i in range(r)]; den = 1
        for row in M:
            for q in row: den = lcm(den, q.denominator)
        Ts.append((den, [[int(q * den) for q in row] for row in M]))
    for U, (X, Y, L), (re, im) in zip(Us, ints, mons):
        total = Fr(0)
        for blk, ab, Tb in zip(blocks, amp, Ts):
            if Tb is None: continue
            den, Mi = Tb; qs = []
            for row in ab:
                x = y = 0
                for c, mi in zip(row, blk):
                    if c: x += c * re[mi]; y += c * im[mi]
                qs.append((x, y))
            acc = 0
            for i, (xi, yi) in enumerate(qs):
                Mrow = Mi[i]
                for j, (xj, yj) in enumerate(qs):
                    if Mrow[j]: acc += Mrow[j] * (xi * xj + yi * yj)
            total += Fr(acc, den)
        bad += (total != F_exact(U, s) * L**(2 * deg) * sd * sd)
okI = bad == 0
print("(I) nodes chosen here: %s ; identity fails at %d of %d (node, unitary) pairs (seed %d): exact identity %s  [%.0fs]" % (", ".join(map(str, nodes)), bad, len(nodes) * nU, seed, "HOLDS" if okI else "FAILS", time.time() - t0), flush=True)
# (E)
p = 32749; R_, I_ = entries_mod(Us, p); ma, mb = monomials_mod(R_, I_, p, MONS)
P = np.array([bl[i] for bl in blocks for i in range(len(bl)) for j in range(i, len(bl))]); Q = np.array([bl[j] for bl in blocks for i in range(len(bl)) for j in range(i, len(bl))])
E = (ma[:, P] * ma[:, Q] + mb[:, P] * mb[:, Q]) % p; rE = flint.nmod_mat(nU, E.shape[1], E.reshape(-1).tolist(), p).rank()
print("(E) dim W = %d ; rank of the evaluation functionals at these unitaries modulo %d: %d (needed: %d)  [%.0fs]" % (E.shape[1], p, rE, rE_expected, time.time() - t0))
print("CHECKS (P), (I), (E) PASSED: F_s >= 0 on U(4) for every s in [%s, %s]" % (a, b) if (okP and okI and rE == rE_expected) else "A CHECK FAILED")
