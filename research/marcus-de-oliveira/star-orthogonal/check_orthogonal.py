"""Exact checks of every identity used in the proof of Theorem O (README.md).  The proof itself is written out in the
README; this script only confirms its algebra, in exact arithmetic.

 (1) Rational orthogonal 4x4 matrices (Cayley transforms; det +1 and, after a row sign change, det -1):
     the second compound in the self-dual / anti-self-dual basis is A (+) B with A, B in SO(3);
     4 D[S_M, S_N] = (A_MN + B_MN)^2 and 4 D[S_M, S_N^c] = (A_MN - B_MN)^2;
     4 (g0 + D[12,12]) = M22^2 + M33^2 + M23^2 + N23^2 + M11^2 - N11^2   (M = A + B, N = A - B);
     the lower 2x2 blocks of A and B satisfy |zeta| + |omega| = 1 and |zeta|^2 - |omega|^2 = A11 (checked in squared form);
     the inequality (star) itself, in the square-root-free form  w >= 0  or  w^2 <= 4 D[12,12] D[13,14],  w = g0 + D[12,12].
 (2) Symbolic identities (sympy): the expansion of E/2, the shift identity for Qf, Qf(1,1,1,1) = 2k^2, the identity
     behind the reversed Cauchy-Schwarz inequality, and the identities of facts (F1)-(F4) in the proof of the lemma.
 (3) The lemma  max(Qf(w), Qf(1-w)) >= 0  on random rational data, including vertices and faces of the cube."""
import itertools, random, sys, time
from fractions import Fraction as Fr
import sympy as sp
t0 = time.time(); rng = random.Random(20261125)
P2 = list(itertools.combinations(range(4), 2)); IDX = {S: k for k, S in enumerate(P2)}
MATCH = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
# rows of the basis change, times sqrt(2): self-dual e12+e34, e13-e24, e14+e23 ; anti-self-dual e12-e34, e13+e24, e14-e23
BASIS = [{(0, 1): 1, (2, 3): 1}, {(0, 2): 1, (1, 3): -1}, {(0, 3): 1, (1, 2): 1}, {(0, 1): 1, (2, 3): -1}, {(0, 2): 1, (1, 3): 1}, {(0, 3): 1, (1, 2): -1}]
def matmul(X, Y): return [[sum(X[i][k] * Y[k][j] for k in range(len(Y))) for j in range(len(Y[0]))] for i in range(len(X))]
def inv(Mx):
    n = len(Mx); A = [row[:] + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(Mx)]
    for c in range(n):
        pr = next(r for r in range(c, n) if A[r][c] != 0); A[c], A[pr] = A[pr], A[c]; d = A[c][c]; A[c] = [v / d for v in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]; A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [row[n:] for row in A]
def cayley(n=4):
    K = [[Fr(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            v = Fr(rng.randint(-6, 6), rng.choice([1, 2, 3, 5])); K[i][j] = v; K[j][i] = -v
    I = [[Fr(int(i == j)) for j in range(n)] for i in range(n)]
    return matmul([[I[i][j] - K[i][j] for j in range(n)] for i in range(n)], inv([[I[i][j] + K[i][j] for j in range(n)] for i in range(n)]))
def det2(O, S, T): return O[S[0]][T[0]] * O[S[1]][T[1]] - O[S[0]][T[1]] * O[S[1]][T[0]]
def det3(A): return (A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]))
def star_holds(O):
    D = lambda S, T: det2(O, S, T)**2
    g0 = O[0][0]**2 + O[2][2]**2 + D((0, 3), (2, 3)) - D((0, 1), (0, 3)); D12 = D((0, 1), (0, 1)); D13 = D((0, 2), (0, 3)); w = g0 + D12
    alt = D((0, 2), (0, 2)) + D((0, 3), (0, 3)) + D((0, 2), (0, 3)) + D((0, 2), (1, 2)) - D((0, 1), (2, 3))      # the six-minor form of g0
    return (w >= 0 or w * w <= 4 * D12 * D13), g0 == alt, w, D12, D13
n1 = 0; bad = 0; nneg = 0
for trial in range(400):
    O = cayley(); assert matmul(O, [list(r) for r in zip(*O)]) == [[Fr(int(i == j)) for j in range(4)] for i in range(4)]
    C2 = [[det2(O, S, T) for T in P2] for S in P2]
    R = [[sum(Fr(ci * cj, 2) * C2[IDX[S]][IDX[T]] for S, ci in BASIS[i].items() for T, cj in BASIS[j].items()) for j in range(6)] for i in range(6)]
    A = [row[:3] for row in R[:3]]; B = [row[3:] for row in R[3:]]
    ok = all(R[i][j] == 0 and R[j][i] == 0 for i in range(3) for j in range(3, 6))
    for X in (A, B): ok &= matmul(X, [list(r) for r in zip(*X)]) == [[Fr(int(i == j)) for j in range(3)] for i in range(3)] and det3(X) == 1
    for m in range(3):
        for n in range(3):
            S = MATCH[m][0]; T, Tc = MATCH[n]
            ok &= 4 * det2(O, S, T)**2 == (A[m][n] + B[m][n])**2 and 4 * det2(O, S, Tc)**2 == (A[m][n] - B[m][n])**2
    M = [[A[i][j] + B[i][j] for j in range(3)] for i in range(3)]; N = [[A[i][j] - B[i][j] for j in range(3)] for i in range(3)]
    holds, same, w, D12, D13 = star_holds(O)
    ok &= same and 4 * w == M[1][1]**2 + M[2][2]**2 + M[1][2]**2 + N[1][2]**2 + M[0][0]**2 - N[0][0]**2
    for X in (A, B):       # lower 2x2 block: zeta, omega as (re, im) pairs
        z = (Fr(1, 2) * (X[1][1] + X[2][2]), Fr(1, 2) * (X[2][1] - X[1][2])); om = (Fr(1, 2) * (X[1][1] - X[2][2]), Fr(1, 2) * (X[2][1] + X[1][2]))
        nz = z[0]**2 + z[1]**2; no = om[0]**2 + om[1]**2            # |zeta|^2, |omega|^2
        ok &= nz - no == X[0][0] and 4 * nz == (1 + X[0][0])**2 and 4 * no == (1 - X[0][0])**2      # i.e. |zeta| = (1+A11)/2, |omega| = (1-A11)/2
    Om = [[-v for v in O[0]]] + O[1:]                                  # det -1: same squared minors
    ok &= all(det2(Om, S, T)**2 == det2(O, S, T)**2 for S in P2 for T in P2) and all(Om[i][j]**2 == O[i][j]**2 for i in range(4) for j in range(4))
    bad += (not ok) or (not holds); n1 += 1; nneg += (w < 0)
print("(1) %d rational matrices in SO(4), each also with one row negated: structure, formulas or (star) failing at %d ; cases with g0 + D[12,12] < 0: %d  [%.0fs]" % (n1, bad, nneg, time.time() - t0), flush=True)
ok1 = bad == 0
# (1b) the hard region.  Two rational base points and small rational perturbations of them:
#      O7 = (I + 4S)/7 on Drury's family (g0 = -158/2401 at one of its symmetric images), with random symmetric images;
#      Ob = Cayley transform of K0, a rational point near the minimiser of g0 + D[12,12] over SO(4) (about -0.0255).
Sk = [[0, 1, 1, -1], [-1, 0, -1, -1], [-1, 1, 0, 1], [1, 1, -1, 0]]; O7 = [[Fr(int(i == j) + 4 * Sk[i][j], 7) for j in range(4)] for i in range(4)]
def cayley_of(vals):
    K = [[Fr(0)] * 4 for _ in range(4)]; it = iter(vals)
    for i in range(4):
        for j in range(i + 1, 4):
            v = next(it); K[i][j] = v; K[j][i] = -v
    I = [[Fr(int(i == j)) for j in range(4)] for i in range(4)]
    return matmul([[I[i][j] - K[i][j] for j in range(4)] for i in range(4)], inv([[I[i][j] + K[i][j] for j in range(4)] for i in range(4)]))
Ob = cayley_of([Fr(-5, 17), Fr(-10, 17), Fr(-10, 17), Fr(10, 17), Fr(-10, 17), Fr(5, 17)])
nb = 0; badb = 0; negb = 0; tight = None; perms = list(itertools.permutations(range(4)))
for trial in range(3000):
    pert = cayley_of([Fr(rng.randint(-3, 3), rng.choice([20, 40, 80, 200])) for _ in range(6)]) if trial % 3 else [[Fr(int(i == j)) for j in range(4)] for i in range(4)]
    if trial % 2 == 0:
        O = matmul(Ob, pert)
    else:
        base = matmul(O7, pert); al = rng.choice(perms); be = rng.choice(perms); O = [[base[al[i]][be[j]] for j in range(4)] for i in range(4)]
        if rng.random() < 0.5: O = [list(r_) for r_ in zip(*O)]
    holds, same, w, D12, D13 = star_holds(O); nb += 1; badb += (not holds) or (not same)
    if w < 0:
        negb += 1; slack = 4 * D12 * D13 - w * w; tight = slack if tight is None else min(tight, slack)
print("(1b) %d rational orthogonal matrices in and near the hard region: (star) failing at %d ; cases with g0 + D[12,12] < 0, where the square-root term is needed: %d ; smallest 4 D12 D13 - w^2 among them: about %.3e  [%.0fs]"
      % (nb, badb, negb, float(tight) if tight is not None else float('nan'), time.time() - t0), flush=True)
ok1 = ok1 and badb == 0 and negb > 0
# (2) symbolic identities
p, r, k, th = sp.symbols('p r k theta', real=True); q = 1 - p; u = 1 - r
x1, x2, x3, x4 = sp.symbols('x1:5', real=True); ca, cb, cg, cd, sa, sb, sg, sd = sp.symbols('ca cb cg cd sa sb sg sd', real=True)
def Qf(w, p_, q_, r_, u_): return p_ * r_ * (w[0]**2 + w[2]**2) + q_ * u_ * (w[1]**2 + w[3]**2) - 2 * p_ * q_ * w[0] * w[1] - 2 * r_ * u_ * w[2] * w[3]
kk = p - u; X = p * x1 - q * x2 + r * x3 - u * x4
idA = sp.expand(Qf((x1 + th, x2 + th, x3 + th, x4 + th), p, q, r, u) - Qf((x1, x2, x3, x4), p, q, r, u) - 2 * th * kk * X - 2 * th**2 * kk**2)
idB = sp.expand(Qf((1, 1, 1, 1), p, q, r, u) - 2 * kk**2)
half_E = (p * ca - r * cg)**2 + (q * cb - u * cd)**2 + (p * sa - q * sb)**2 + (r * sg - u * sd)**2 + 2 * (p - q) * (r - u)
target = p * r * (ca - cg)**2 + q * u * (cb - cd)**2 + Qf((sa, sb, sg, sd), p, q, r, u) + 2 * kk**2
idC = sp.expand(half_E - target)
for c_, s_ in ((ca, sa), (cb, sb), (cg, sg), (cd, sd)): idC = sp.expand(idC.subs(c_**2, 1 - s_**2))
X1, X2, Y1, Y2 = sp.symbols('X1 X2 Y1 Y2', real=True)
idD = sp.expand((X1 * X2 - Y1 * Y2)**2 - (X1**2 - Y1**2) * (X2**2 - Y2**2) - (X1 * Y2 - X2 * Y1)**2)
m11, m12, m21, m22 = sp.symbols('m11 m12 m21 m22', real=True); zr, zi, wr, wi = (m11 + m22) / 2, (m21 - m12) / 2, (m11 - m22) / 2, (m21 + m12) / 2
idE = [sp.expand((zr**2 + zi**2) - (wr**2 + wi**2) - (m11 * m22 - m12 * m21)), sp.expand(2 * (zr**2 + zi**2 + wr**2 + wi**2) - (m11**2 + m12**2 + m21**2 + m22**2)),
       sp.expand(zr + wr - m11), sp.expand(zr - wr - m22), sp.expand(wi - zi - m12), sp.expand(zi + wi - m21)]
# facts (F1), (F2), (F3) and the two identities used in (F4), for general p, q, r, u (only p + q = r + u is used in F3)
pp, qq, rr, uu = sp.symbols('p q r u', real=True); z1, z2, y1, y2, v1, v2, t1, t2 = sp.symbols('z1 z2 y1 y2 v1 v2 t1 t2', real=True)
Q1 = lambda a_, b_: pp * rr * a_**2 - 2 * pp * qq * a_ * b_ + qq * uu * b_**2; Q2 = lambda a_, b_: pp * rr * a_**2 - 2 * rr * uu * a_ * b_ + qq * uu * b_**2
B1 = lambda a_, b_, c_, d_: pp * rr * a_ * c_ - pp * qq * (a_ * d_ + b_ * c_) + qq * uu * b_ * d_; B2 = lambda a_, b_, c_, d_: pp * rr * a_ * c_ - rr * uu * (a_ * d_ + b_ * c_) + qq * uu * b_ * d_
F1 = sp.expand(Q1(z1, z2) * Q1(y1, y2) - B1(z1, z2, y1, y2)**2 - pp * qq * (rr * uu - pp * qq) * (z1 * y2 - z2 * y1)**2)
F2 = sp.expand(B2(v1, v2, t1, t2)**2 - Q2(v1, v2) * Q2(t1, t2) - rr * uu * (rr * uu - pp * qq) * (v1 * t2 - v2 * t1)**2)
F3a = sp.expand(Q1(z1 + y1, z2 + y2) - Q1(z1, z2) - Q1(y1, y2) - 2 * B1(z1, z2, y1, y2)); F3b = sp.expand(Q2(v1 + t1, v2 + t2) - Q2(v1, v2) - Q2(t1, t2) - 2 * B2(v1, v2, t1, t2))
F3c = sp.expand((Q1(1, 1) + Q2(1, 1) - 2 * (pp - uu) * (rr - qq)))
F4a = sp.expand((pp * rr * v1**2 + qq * uu * v2**2) * (pp * rr * t1**2 + qq * uu * t2**2) - (pp * rr * v1 * t1 + qq * uu * v2 * t2)**2 - pp * rr * qq * uu * (v1 * t2 - v2 * t1)**2)
F4b = sp.expand((v1 * t2 + v2 * t1)**2 - 4 * v1 * v2 * t1 * t2 - (v1 * t2 - v2 * t1)**2)
F4c = sp.expand(-B2(v1, v2, t1, t2) - (rr * uu * (v1 * t2 + v2 * t1) - (pp * rr * v1 * t1 + qq * uu * v2 * t2)))
okF = all(e == 0 for e in (F1, F2, F3a, F3b, F3c, F4a, F4b, F4c))
ok2 = idA == 0 and idB == 0 and idC == 0 and idD == 0 and all(e == 0 for e in idE) and okF
print("(2) symbolic identities: shift identity %s ; Qf(1,1,1,1) = 2k^2 %s ; expansion of E/2 %s ; reversed Cauchy-Schwarz identity %s ; 2x2 block identities %s ; facts F1-F4 of the lemma %s" % (idA == 0, idB == 0, idC == 0, idD == 0, all(e == 0 for e in idE), okF), flush=True)
# (3) the lemma on random rational data
def Qn(w, p_, r_):
    q_ = 1 - p_; u_ = 1 - r_; return p_ * r_ * (w[0]**2 + w[2]**2) + q_ * u_ * (w[1]**2 + w[3]**2) - 2 * p_ * q_ * w[0] * w[1] - 2 * r_ * u_ * w[2] * w[3]
viol = 0; n3 = 0; smallest = None
for trial in range(200000):
    p_ = Fr(rng.randint(0, 40), 40); r_ = Fr(rng.randint(0, 40), 40); w = [Fr(rng.randint(0, 24), 24) for _ in range(4)]
    if trial % 3 == 0: w = [Fr(rng.randint(0, 1)) if rng.random() < 0.5 else v for v in w]
    m = max(Qn(w, p_, r_), Qn([1 - v for v in w], p_, r_)); viol += m < 0; n3 += 1; smallest = m if smallest is None else min(smallest, m)
ok3 = viol == 0
print("(3) lemma: max(Qf(w), Qf(1-w)) < 0 in %d of %d random rational cases ; smallest value %s  [%.0fs]" % (viol, n3, smallest, time.time() - t0))
print("ALL CHECKS PASSED" if (ok1 and ok2 and ok3) else "A CHECK FAILED")
