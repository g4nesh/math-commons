"""Exact checks of the identities used in the proof of Theorem F (README.md), on Gaussian-rational 4x4 unitaries.

 (J)  Jacobi for unitary U:   det U[I^c, T^c] = (-1)^(sum I + sum T) * det(U) * conj(det U[I, T])        (1-based sums)
 (P)  Pluecker for rows and for columns:  C[I,12]C[I,34] - C[I,13]C[I,24] + C[I,14]C[I,23] = 0,  and the same with the
      roles of rows and columns exchanged
 (T)  for each column pair T:  sum_M s_M * C[S_M,T] * conj(C[S_M,T^c]) = 0  with signs s_M = (-1)^(M+1+sum S_M+sum T^c);
      hence the moduli |C[S_M,T] C[S_M,T^c]|, M = 1,2,3, satisfy the triangle inequality (checked in squared form)
 (F)  for S = U^T diag(e) U with rational unit e:   det S[T,T^c] = sum_M ( e_{S_M} c_M + e1 e2 e3 e4 conj(e_{S_M}) d^2 conj(c_M) ),
      c_M = C[S_M,T] C[S_M,T^c], d = det U; and  det S[12,34] - det S[13,24] + det S[14,23] = 0
 (U)  the algebra of Step F: for S = Lambda + x x^T (symbolic check of the off-diagonal entry of S S^*)."""
import itertools, random, sys, time
from fractions import Fraction as Fr
sys.path.insert(0, '../star-exact')
from exact_space import rational_unitary, is_unitary, CQ
t0 = time.time(); rng = random.Random(7)
P2 = list(itertools.combinations(range(4), 2)); MATCH = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
def mul(a, b): return a * b
def det2(U, S, T): return U[S[0]][T[0]] * U[S[1]][T[1]] - U[S[0]][T[1]] * U[S[1]][T[0]]
def det4(U):
    tot = CQ(0)
    for p in itertools.permutations(range(4)):
        sg = 1
        for i in range(4):
            for j in range(i + 1, 4): sg *= -1 if p[i] > p[j] else 1
        term = CQ(sg)
        for i in range(4): term = term * U[i][p[i]]
        tot = tot + term
    return tot
def eq(a, b): return a.re == b.re and a.im == b.im
def comp(S): return tuple(sorted(set(range(4)) - set(S)))
def sgn(S, T): return -1 if (sum(S) + sum(T)) % 2 else 1        # (-1)^(sum S + sum T); the parity is the same for 0- and 1-based labels of 2-sets
TRIPLES = [(3, 4, 5), (5, 12, 13), (8, 15, 17), (7, 24, 25), (20, 21, 29)]
def unit(): a, b, h = rng.choice(TRIPLES); return CQ(Fr(a, h), Fr(b, h)) if rng.random() < 0.5 else CQ(Fr(-b, h), Fr(a, h))
bad = dict(J=0, P=0, T=0, F=0); n = 0; tri_min = None
for trial in range(300):
    U = rational_unitary(rng)
    if not is_unitary(U):
        raise ValueError(f"Trial {trial} did not produce a unitary matrix")
    d = det4(U); C = {(S, T): det2(U, S, T) for S in P2 for T in P2}; n += 1
    for S in P2:
        for T in P2:
            bad['J'] += not eq(C[(comp(S), comp(T))], CQ(sgn(S, T)) * d * C[(S, T)].conj())
    for I in P2:
        bad['P'] += not eq(C[(I, (0, 1))] * C[(I, (2, 3))] - C[(I, (0, 2))] * C[(I, (1, 3))] + C[(I, (0, 3))] * C[(I, (1, 2))], CQ(0))
        bad['P'] += not eq(C[((0, 1), I)] * C[((2, 3), I)] - C[((0, 2), I)] * C[((1, 3), I)] + C[((0, 3), I)] * C[((1, 2), I)], CQ(0))
    for T, Tc in MATCH:
        tot = CQ(0); mods = []
        for M, (S, Sc) in enumerate(MATCH):
            s = (-1) ** M * sgn(S, Tc); tot = tot + CQ(s) * C[(S, T)] * C[(S, Tc)].conj(); mods.append(C[(S, T)].abs2() * C[(S, Tc)].abs2())      # squared moduli
        bad['T'] += not eq(tot, CQ(0))
        # triangle inequality a <= b + c for moduli, in squared form: (a^2 - b^2 - c^2)^2 <= 4 b^2 c^2  or  a^2 <= b^2 + c^2
        for M in range(3):
            a2 = mods[M]; b2, c2 = [mods[k] for k in range(3) if k != M]; ok = a2 <= b2 + c2 or (a2 - b2 - c2) ** 2 <= 4 * b2 * c2
            bad['T'] += not ok
    e = [unit() for _ in range(4)]; eps0 = e[0] * e[1] * e[2] * e[3]
    S_ = [[CQ(0)] * 4 for _ in range(4)]
    for j in range(4):
        for k in range(4):
            tot = CQ(0)
            for i in range(4): tot = tot + e[i] * U[i][j] * U[i][k]
            S_[j][k] = tot
    Fs = []
    for T, Tc in MATCH:
        lhs = det2(S_, T, Tc); rhs = CQ(0)
        for S, Sc in MATCH:
            eI = e[S[0]] * e[S[1]]; c = C[(S, T)] * C[(S, Tc)]; rhs = rhs + eI * c + eps0 * eI.conj() * d * d * c.conj()
        bad['F'] += not eq(lhs, rhs); Fs.append(lhs)
    bad['F'] += not eq(Fs[0] - Fs[1] + Fs[2], CQ(0))
print("(J)(P)(T)(F) on %d Gaussian-rational unitaries: failures %s  [%.0fs]" % (n, bad, time.time() - t0), flush=True)
import sympy as sp
L1, L2, x1, x2, x3, x4 = sp.symbols('L1 L2 x1 x2 x3 x4'); L1c, L2c, x1c, x2c, x3c, x4c = sp.symbols('L1c L2c x1c x2c x3c x4c')     # conjugates as independent symbols
xs = [x1, x2, x3, x4]; xc = [x1c, x2c, x3c, x4c]; n2 = sum(a * b for a, b in zip(xs, xc))
# entry (1,2) of S S^* for S = diag(L) + x x^T:  sum_l S_1l conj(S_2l)
S1 = [L1 + x1 * x1, x1 * x2, x1 * x3, x1 * x4]; S2c = [x2c * x1c, L2c + x2c * x2c, x2c * x3c, x2c * x4c]
offd = sp.expand(sum(a * b for a, b in zip(S1, S2c)) - (L1 * x1c * x2c + L2c * x1 * x2 + n2 * x1 * x2c))
diag = sp.expand(sum(a * b for a, b in zip(S1, [L1c + x1c * x1c, x1c * x2c, x1c * x3c, x1c * x4c])) - (L1 * L1c + L1 * x1c**2 + L1c * x1**2 + n2 * x1 * x1c))
print("(U) entries of S S^* for S = Lambda + x x^T: off-diagonal formula %s ; diagonal formula %s" % (offd == 0, diag == 0))
if any(bad.values()) or offd != 0 or diag != 0:
    raise SystemExit("A CHECK FAILED")
print("ALL CHECKS PASSED")
