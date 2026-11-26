"""Is Zheng's rational orthogonal O* on the family O(y) = y I + sqrt((1-y^2)/3) S ?   (review finding 2)
Exact integer check: 7*O(1/7) = I + 4S is NOT equal to 7*O*, but it equals D1 P_alpha V P_beta D2 with V in {O*, O*^T},
permutation matrices P and diagonal sign matrices D.  Such maps preserve every squared minor up to relabelling, so the
value of the exotic functional over the symmetry orbit is unchanged: g(O(1/7)) = -158/2401 is attained by an image of O*."""
import itertools
import numpy as np
S = np.array([[0, 1, 1, -1], [-1, 0, -1, -1], [-1, 1, 0, 1], [1, 1, -1, 0]]); OY = np.eye(4, dtype=int) + 4 * S        # 7 * O(1/7)
OZ = np.array([[1, 4, 4, 4], [-4, 1, 4, -4], [-4, -4, 1, 4], [-4, 4, -4, 1]])                                          # 7 * O*
assert (OY @ OY.T == 49 * np.eye(4, dtype=int)).all() and (OZ @ OZ.T == 49 * np.eye(4, dtype=int)).all()
print("O(1/7) == O*:", bool((OY == OZ).all()))
count = 0; example = None
for tr in (False, True):
    V = OZ.T if tr else OZ
    for al in itertools.permutations(range(4)):
        for be in itertools.permutations(range(4)):
            W = V[np.ix_(al, be)]
            if not (np.abs(W) == np.abs(OY)).all(): continue
            # signs: OY[i,j] = d1[i] * W[i,j] * d2[j]; fix d1[0] = 1
            d2 = OY[0] // W[0]; d1 = OY[:, 0] // (W[:, 0] * d2[0])
            if (d1[:, None] * W * d2[None, :] == OY).all():
                count += 1; example = example or (tr, al, be, d1.tolist(), d2.tolist())
print("representations O(1/7) = D1 P_alpha V P_beta D2 (D1 normalised to start with +1): %d ; one of them: transposed=%s, alpha=%s, beta=%s, D1=%s, D2=%s" % ((count,) + example))
def g_num(U7):
    """numerator of g(U) over 2401 for U = U7 / 7"""
    D = lambda A, B: (U7[A[0], B[0]] * U7[A[1], B[1]] - U7[A[0], B[1]] * U7[A[1], B[0]])**2
    return 49 * (U7[0, 0]**2 + U7[2, 2]**2) + D((0, 3), (2, 3)) - D((0, 1), (0, 3))
print("g(O(1/7)) = %d/2401 (exact)" % g_num(OY))
imgs = [g_num(V[np.ix_(al, be)]) for V in (OZ, OZ.T) for al in itertools.permutations(range(4)) for be in itertools.permutations(range(4))]
print("min of g over the 1152 images P_alpha O* P_beta and transposes: %d/2401 ; attained %d times" % (min(imgs), imgs.count(min(imgs))))
