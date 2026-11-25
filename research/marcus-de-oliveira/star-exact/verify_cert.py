"""Independent exact check of certificate_t<t>.json.  Uses only exact rational arithmetic.

Claim checked:  for every U in U(4),
    t*g0 + (1+t)*D[12,12] + t^2*D[13,14]  =  sum_b (A_b m_b(U))^* T_b (A_b m_b(U)),      T_b positive definite,
hence the left side is >= 0 on U(4).   (g0 = B11 + B33 + D[14,34] - D[12,14], D[S,T] = |det U[S,T]|^2, 1-based.)

Steps:
 (1) Lemma 1: rank K = 1227 and rank of evaluations at FRESH rational unitaries = 261, with 1227 + 261 = 1488 = dim W.
     Then a function in W vanishes on U(4) iff it vanishes at those unitaries.
 (2) Every T_b is symmetric with positive leading principal minors.
 (3) The identity holds exactly at each of the fresh rational unitaries (they were not used to construct the certificate)."""
import sys, json, random, time
from fractions import Fraction as Fr
import flint
from exact_space import (MONS, BLOCK_LIST, NW, PAIRS, kernel_generators, rational_unitary, is_unitary, monomial_values, evaluation_row, det2, CQ)
t0 = time.time(); path = sys.argv[1]; seed = int(sys.argv[2]) if len(sys.argv) > 2 else 777
cert = json.load(open(path)); tq = Fr(cert['t'])
amp = [[[Fr(q) for q in row] for row in a] for a in cert['amplitudes']]; T = [[[Fr(q) for q in row] for row in Tb] for Tb in cert['T']]
assert len(amp) == len(BLOCK_LIST) == len(T)
# (1)
gens = kernel_generators(); M = flint.fmpz_mat(len(gens), NW)
for r, row in enumerate(gens):
    for k, v in row.items(): M[r, k] = v
rK = M.rank()
rng = random.Random(seed); Us = [rational_unitary(rng) for _ in range(300)]; assert all(is_unitary(U) for U in Us)
E = flint.fmpq_mat(300, NW, [flint.fmpq(x.numerator, x.denominator) for U in Us for x in evaluation_row(U)]); rE = E.rank()
ok1 = (rK + rE == NW)
print("(1) rank K = %d ; rank of evaluations at 300 fresh rational unitaries (seed %d) = %d ; sum = %d = dim W: %s  [%.0fs]" % (rK, seed, rE, rK + rE, ok1, time.time() - t0), flush=True)
# (2)
ok2 = True; nblk = 0; minpm = None
for a, Tb in zip(amp, T):
    r = len(a); assert len(Tb) == r
    if r == 0: continue
    nblk += 1
    for i in range(r):
        assert len(a[i]) > 0 and len(Tb[i]) == r
        for j in range(r): ok2 &= (Tb[i][j] == Tb[j][i])
    for k in range(1, r + 1):
        d = flint.fmpq_mat(k, k, [flint.fmpq(Tb[i][j].numerator, Tb[i][j].denominator) for i in range(k) for j in range(k)]).det(); d = Fr(int(d.p), int(d.q))
        ok2 &= (d > 0); minpm = d if minpm is None else min(minpm, d)
print("(2) %d nonempty blocks, %d amplitudes ; all T_b symmetric with positive leading principal minors: %s (smallest minor %.3e)" % (nblk, sum(len(a) for a in amp), ok2, float(minpm)), flush=True)
# (3)
def f_exact(U):
    B = lambda s, j: U[s][j].abs2(); D = lambda S, Tt: det2(U, S, Tt).abs2()
    g0 = B(0, 0) + B(2, 2) + D((0, 3), (2, 3)) - D((0, 1), (0, 3))
    return tq * g0 + (1 + tq) * D((0, 1), (0, 1)) + tq * tq * D((0, 2), (0, 3))
bad = 0
for U in Us:
    mv = monomial_values(U); total = Fr(0)
    for b, a, Tb in zip(BLOCK_LIST, amp, T):
        if not a: continue
        qs = []
        for row in a:
            z = CQ(0)
            for coef, mi in zip(row, b):
                if coef: z = z + CQ(coef) * mv[mi]
            qs.append(z)
        for i in range(len(qs)):
            for j in range(len(qs)): total += Tb[i][j] * (qs[i].re * qs[j].re + qs[i].im * qs[j].im)
    bad += (total != f_exact(U))
ok3 = (bad == 0)
print("(3) identity fails at %d of 300 fresh rational unitaries: exact identity %s  [%.0fs]" % (bad, "HOLDS" if ok3 else "FAILS", time.time() - t0))
print("CERTIFICATE VALID" if (ok1 and ok2 and ok3) else "CERTIFICATE NOT VALID")
