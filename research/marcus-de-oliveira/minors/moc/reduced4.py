"""Reduced coordinates for n=4: p_U is determined by (M1, N) where
M1 = |U|^2 (4x4 unistochastic) and N is the 3x3 'pairing matrix'
N[pi,pi'] = |det U[S,T]|^2 + |det U[S,T^c]|^2,  pi={S,S^c}, pi'={T,T^c}."""
import itertools, numpy as np
PAIRS = [((0,1),(2,3)), ((0,2),(1,3)), ((0,3),(1,2))]

def red_coords(U):
    M1 = np.abs(U)**2
    N = np.zeros((3,3))
    for p,(S,Sc) in enumerate(PAIRS):
        for q,(T,Tc) in enumerate(PAIRS):
            N[p,q] = abs(np.linalg.det(U[np.ix_(S,T)]))**2 + abs(np.linalg.det(U[np.ix_(S,Tc)]))**2
    return M1, N

def perm_matrix(s):
    P = np.zeros((len(s),len(s)))
    for i,j in enumerate(s): P[i,j]=1
    return P

PERMS = list(itertools.permutations(range(4)))
def vec(M1,N):
    # affine coords: M1[0:3,0:3] (9) and N[0:2,0:2] (4)
    return np.concatenate([M1[:3,:3].ravel(), N[:2,:2].ravel()])

VERTS = np.array([vec(*red_coords(perm_matrix(s))) for s in PERMS])
