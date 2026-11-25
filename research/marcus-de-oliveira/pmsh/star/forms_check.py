"""Numerical check of two reformulations of (*) used in equivalence.md.
(1) Complex coordinates of R in SO(6): each 2x2 block R_MN acts on C as z -> zeta_MN z + omega_MN conj(z), and
    D[S_M, T_N] = |zeta_MN|^2,  D[S_M, T_N^c] = |omega_MN|^2   (S_M, T_N = the members of the matchings containing index 1).
(2) Block form: with A = R[blocks 2,3 ; blocks 2,3] (4x4 real), singular values (1, 1, s1, s2), and det A < 0,
    margin of (*)  =  |z11|^2 + |z22|^2 + (1/2)|A12|_F^2 + (s1 - s2)|z12| + det A ,
    where z11, z22, z12 are the complex-linear parts of the 2x2 blocks of A."""
import numpy as np
from base import *
from so6 import R_of
def parts(Bk): return ((Bk[0, 0] + Bk[1, 1]) + 1j * (Bk[1, 0] - Bk[0, 1])) / 2, ((Bk[0, 0] - Bk[1, 1]) + 1j * (Bk[0, 1] + Bk[1, 0])) / 2
rng = np.random.default_rng(3); e1 = e2 = e3 = 0; nneg = 0; worst = np.inf
for _ in range(4000):
    U = haar(rng); x = x_U(U); R = R_of(U)
    z = np.zeros((3, 3), complex); w = np.zeros((3, 3), complex)
    for M in range(3):
        for N in range(3): z[M, N], w[M, N] = parts(R[2 * M:2 * M + 2, 2 * N:2 * N + 2])
    for M, (S, Sc) in enumerate(MATCH):
        for N, (T, Tc) in enumerate(MATCH):
            e1 = max(e1, abs(abs(z[M, N])**2 - x[iD(S, T)]), abs(abs(w[M, N])**2 - x[iD(S, Tc)]))
    D = lambda S, T: x[iD(S, T)]
    margin = D((0,2),(0,2)) + D((0,3),(0,3)) + D((0,2),(1,2)) + (np.sqrt(D((0,1),(0,1))) + np.sqrt(D((0,2),(0,3))))**2 - D((0,1),(2,3))
    worst = min(worst, margin)
    A = R[2:6, 2:6]; sv = np.linalg.svd(A, compute_uv=False); dA = np.linalg.det(A)
    e2 = max(e2, abs(sv[0] - 1), abs(sv[1] - 1), abs(dA - np.linalg.det(R[0:2, 0:2])))
    if dA < 0:
        nneg += 1
        blockform = abs(z[1, 1])**2 + abs(z[2, 2])**2 + 0.5 * np.linalg.norm(R[2:4, 4:6])**2 + (sv[2] - sv[3]) * abs(z[1, 2]) + dA
        e3 = max(e3, abs(blockform - margin))
print("(1) |zeta|^2, |omega|^2 reproduce all 18 squared minors: max error %.1e" % e1)
print("(2) A has two unit singular values and det A = det R_11: max error %.1e ; cases with det A < 0: %d of 4000" % (e2, nneg))
print("    block form equals the margin of (*) on those cases: max error %.1e" % e3)
print("smallest margin of (*) over 4000 Haar unitaries: %.4f" % worst)
