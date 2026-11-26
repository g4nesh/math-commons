"""Independent exact check of a degree-4 certificate (certificate4_t<t>.json).

Claim:  for every U in U(4),  t*g0 + (1+t)*D[12,12] + t^2*D[13,14] = sum_b (A_b m_b)^* T_b (A_b m_b)  with T_b positive
definite; hence the left side is >= 0 on U(4).

 (2) every T_b is symmetric with positive leading principal minors            (exact rational arithmetic)
 (3) the identity holds exactly at N fresh rational unitaries                 (exact integer arithmetic)
 (E) the evaluation functionals at those unitaries have rank 954 on the 11208-dimensional space W4 of
     conjugation-symmetric torus-invariant (4,4)-forms                         (rank modulo a prime: a lower bound over Q)
Together with rank K4 >= 10254 (lemma4.py), (E) shows that the restriction of W4 to U(4) has dimension exactly 954 and that
an element of W4 vanishing at these unitaries vanishes on U(4); then (3) gives the identity on all of U(4)."""
import sys, json, random, time
from fractions import Fraction as Fr
import numpy as np, flint
from exact4 import *
t0 = time.time(); path = sys.argv[1]; seed = int(sys.argv[2]) if len(sys.argv) > 2 else 777; nU = int(sys.argv[3]) if len(sys.argv) > 3 else 1100
cert = json.load(open(path)); tq = Fr(cert['t']); amp = cert['amplitudes']; T = [[[Fr(q) for q in row] for row in Tb] for Tb in cert['T']]
assert len(amp) == len(BLOCK_LIST) == len(T)
# (2)
ok2 = True; minpm = None; namp = 0
for a, Tb in zip(amp, T):
    r = len(a); assert len(Tb) == r; namp += r
    for row in a: assert len(row) == len(BLOCK_LIST[amp.index(a)]) if False else True
    for i in range(r):
        assert len(Tb[i]) == r
        for j in range(r): ok2 &= (Tb[i][j] == Tb[j][i])
    for k in range(1, r + 1):
        d = flint.fmpq_mat(k, k, [flint.fmpq(Tb[i][j].numerator, Tb[i][j].denominator) for i in range(k) for j in range(k)]).det(); d = Fr(int(d.p), int(d.q))
        ok2 &= d > 0; minpm = d if minpm is None else min(minpm, d)
for a, b in zip(amp, BLOCK_LIST):
    for row in a: assert len(row) == len(b) and all(isinstance(c, int) for c in row)
print("(2) %d amplitudes in %d nonempty blocks ; all T_b symmetric with positive leading principal minors: %s (smallest %.3e)  [%.0fs]" % (namp, sum(1 for a in amp if a), ok2, float(minpm), time.time() - t0), flush=True)
# (3)
rng = random.Random(seed); Us = [rational_unitary(rng) for _ in range(nU)]; assert all(is_unitary(U) for U in Us)
blocks = []
for Tb in T:
    r = len(Tb)
    if r == 0: blocks.append(None); continue
    den = 1
    for row in Tb:
        for q in row: den = lcm(den, q.denominator)
    blocks.append((den, [[int(q * den) for q in row] for row in Tb]))
bad = 0
for U in Us:
    X, Y, L = integer_form(U); re, im = monomials_int(X, Y); qs = amplitudes_int(re, im, amp); total = Fr(0)
    for q, blk in zip(qs, blocks):
        if blk is None: continue
        den, Nn = blk; s = 0; r = len(q)
        for i in range(r):
            xi, yi = q[i]; Ni = Nn[i]
            for j in range(r):
                if Ni[j]: s += Ni[j] * (xi * q[j][0] + yi * q[j][1])
        total += Fr(s, den)
    bad += (total != f_exact(U, tq) * L**8)
ok3 = bad == 0
print("(3) identity fails at %d of %d fresh rational unitaries (seed %d): exact identity %s  [%.0fs]" % (bad, nU, seed, "HOLDS" if ok3 else "FAILS", time.time() - t0), flush=True)
# (E)
p = 32749; R_, I_ = entries_mod(Us, p); a, b = monomials_mod(R_, I_, p)
P = np.array([bl[i] for bl in BLOCK_LIST for i in range(len(bl)) for j in range(i, len(bl))]); Q = np.array([bl[j] for bl in BLOCK_LIST for i in range(len(bl)) for j in range(i, len(bl))])
E = (a[:, P] * a[:, Q] + b[:, P] * b[:, Q]) % p
rE = flint.nmod_mat(nU, E.shape[1], E.reshape(-1).tolist(), p).rank()
print("(E) dim W4 = %d ; rank of the evaluation functionals at these unitaries modulo %d: %d (expected 954)  [%.0fs]" % (E.shape[1], p, rE, time.time() - t0))
print("CHECKS (2), (3), (E) PASSED" if (ok2 and ok3 and rE == 954) else "A CHECK FAILED")
