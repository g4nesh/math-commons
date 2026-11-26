"""Exact replay for "A shared spectral value in order four: a Birkhoff construction"
(final_referee_shared_zero_birkhoff_20261022).  Own code only; nothing from the supplement is
read or executed.  One printed line per check; exit status nonzero iff an EXACT check fails.

Conventions (the document's): indices 0..3; sigma acts as P_sigma with (P_sigma)_{i,sigma(i)} = 1;
det(diag(a) - U diag(c) U^*) = sum_sigma p_sigma prod_i (a_i - c_sigma(i)), with a_0 = c_0 = d.
In the paper's plus-sign notation c = -b.
"""
import itertools
import random
import sys
from fractions import Fraction as Fr

import sympy as sp

FAIL = []
PERMS = list(itertools.permutations(range(4)))


def check(name, ok, level="EXACT"):
    print(f"[{level}] {'PASS' if ok else 'FAIL'}  {name}")
    if not ok and level == "EXACT":
        FAIL.append(name)


# ---------------------------------------------------------------- exact unitaries
def gaussian_unitary(rng, n=4, size=3):
    """Cayley transform (I - S)(I + S)^{-1} of a Gaussian-rational skew-Hermitian S: exactly unitary."""
    S = sp.zeros(n, n)
    for i in range(n):
        S[i, i] = sp.I * sp.Rational(rng.randint(-size, size), rng.randint(1, size))
        for j in range(i + 1, n):
            v = sp.Rational(rng.randint(-size, size), rng.randint(1, size)) + \
                sp.I * sp.Rational(rng.randint(-size, size), rng.randint(1, size))
            S[i, j] = v
            S[j, i] = -sp.conjugate(v)
    Id = sp.eye(n)
    U = ((Id - S) * (Id + S).inv(method="LU")).applyfunc(lambda e: sp.expand(sp.radsimp(e)))
    return U


def real_rational_SO4(rng, size=3):
    S = sp.zeros(4, 4)
    for i in range(4):
        for j in range(i + 1, 4):
            v = sp.Rational(rng.randint(-size, size), rng.randint(1, size))
            S[i, j] = v
            S[j, i] = -v
    Id = sp.eye(4)
    return (Id - S) * (Id + S).inv()


def is_unitary(U):
    n = U.shape[0]
    return all(sp.expand(sp.radsimp(sp.expand(e))) == 0 for e in (U * U.H - sp.eye(n)))


def abs2(z):
    z = sp.expand(sp.radsimp(sp.expand(z)))
    re_, im_ = z.as_real_imag()
    r = sp.Rational(sp.expand(re_**2 + im_**2))
    return Fr(int(r.p), int(r.q))


def sub(U, I, J):
    return U.extract(list(I), list(J))


# ---------------------------------------------------------------- the construction
def birkhoff_exact(W):
    """Greedy exact Birkhoff decomposition of a doubly stochastic Fraction matrix."""
    W = [row[:] for row in W]
    out = {}
    while True:
        tot = sum(W[0])
        if tot == 0:
            break
        best = None
        for p in PERMS:
            if all(W[i][p[i]] > 0 for i in range(4)):
                mn = min(W[i][p[i]] for i in range(4))
                best = (p, mn)
                break
        if best is None:
            raise RuntimeError("no permutation in support: not doubly stochastic")
        p, mn = best
        out[p] = out.get(p, 0) + mn
        for i in range(4):
            W[i][p[i]] -= mn
    return out


def phi(pi):
    r = pi.index(0)
    k, l = [i for i in range(1, 4) if i != r]
    s = list(pi)
    s[k], s[l] = pi[l], pi[k]
    return tuple(s)


N = (1, 2, 3)


def m_of_sigma(sig):
    return {(i, j): int({sig[0], sig[i]} == set(N) - {j}) for i in N for j in N}


def shared_zero_law(U):
    """The document's law for the shared pair (row index 0, column index 0)."""
    u00 = abs2(U[0, 0])
    s = 1 - u00
    m = {(i, j): abs2(sub(U, (0, i), [k for k in N if k != j]).det()) for i in N for j in N}
    if s == 0:
        return {(0, 1, 2, 3): Fr(1)}, m, s
    W = [[Fr(0)] + [(abs2(U[0, j])) / s for j in N]]
    for i in N:
        W.append([(abs2(U[i, 0])) / s] + [Fr(m[(i, j)]) / s for j in N])
    dec = birkhoff_exact(W)
    law = {}
    for pi, wt in dec.items():
        assert pi[0] != 0
        sg = phi(pi)
        law[sg] = law.get(sg, 0) + s * wt
    if s != 1:
        law[(0, 1, 2, 3)] = law.get((0, 1, 2, 3), 0) + 1 - s
    return law, m, s


def law_for_cross_pair(U, i0, j0):
    """Law for the shared pair a_{i0} = c_{j0}, by relabeling: U'[i,j] = U[tr(i), tc(j)],
    sigma = tc o sigma' o tr."""
    tr = list(range(4)); tr[0], tr[i0] = tr[i0], tr[0]
    tc = list(range(4)); tc[0], tc[j0] = tc[j0], tc[0]
    Up = sp.Matrix(4, 4, lambda i, j: U[tr[i], tc[j]])
    lawp, _, _ = shared_zero_law(Up)
    law = {}
    for sp_, w in lawp.items():
        sg = tuple(tc[sp_[tr[k]]] for k in range(4))
        law[sg] = law.get(sg, 0) + w
    return law


def det_poly(U, a, c):
    n = U.shape[0]
    return sp.expand((sp.diag(*a) - U * sp.diag(*c) * U.H).det(method="berkowitz"))


def law_poly(law, a, c):
    return sp.expand(sum(sp.Rational(w.numerator, w.denominator) *
                         sp.prod([a[i] - c[sg[i]] for i in range(4)]) for sg, w in law.items()))


# ---------------------------------------------------------------- checks
def main():
    rng = random.Random(20261125)
    d, a1, a2, a3, c1, c2, c3 = sp.symbols("d a1 a2 a3 c1 c2 c3")

    # (6)-(5) as a generic identity given the marginal relations (2)
    msym = {(i, j): sp.Symbol(f"m{i}{j}") for i in N for j in N}
    A = {1: a1, 2: a2, 3: a3}
    Cc = {1: c1, 2: c2, 3: c3}
    astar, cstar = a1*a2*a3, c1*c2*c3
    ahat = {i: sp.prod([A[k] for k in N if k != i]) for i in N}
    chat = {j: sp.prod([Cc[l] for l in N if l != j]) for j in N}
    sS = sum(msym.values()) / 2
    col = {j: sS - sum(msym[(i, j)] for i in N) for j in N}   # |u0j|^2 from (2)
    row = {i: sS - sum(msym[(i, j)] for j in N) for i in N}   # |ui0|^2 from (2)
    rhs5 = -astar*sum(Cc[j]*col[j] for j in N) + sum(ahat[i]*chat[j]*msym[(i, j)] for i in N for j in N) \
        - cstar*sum(A[i]*row[i] for i in N)
    H = {(i, j): ahat[i]*chat[j] + astar*(Cc[j] - sp.Rational(1, 2)*sum(Cc.values()))
         + cstar*(A[i] - sp.Rational(1, 2)*sum(A.values())) for i in N for j in N}
    check("(5)+(2) => (6): -a*sum c_j|u0j|^2 + sum ahat chat m - c*sum a_i|ui0|^2 = sum H_ij m_ij (generic m)",
          sp.expand(rhs5 - sum(H[k]*msym[k] for k in msym)) == 0)

    # permutation-matrix identity, all 24 sigma, a0 = c0 = 0
    a = [0, a1, a2, a3]
    c = [0, c1, c2, c3]
    ok = True
    for sg in PERMS:
        ms = m_of_sigma(sg)
        ok &= sp.expand(sp.prod([a[i] - c[sg[i]] for i in range(4)]) - sum(H[k]*ms[k] for k in ms)) == 0
    check("prod_i (a_i - c_sigma(i)) = sum H_ij m_ij(sigma) for all 24 sigma (a0=c0=0)", ok)

    # m(sigma) equals |det P_sigma[{0,i}, N\{j}]|^2
    ok = True
    for sg in PERMS:
        Pm = sp.Matrix(4, 4, lambda i, j: 1 if sg[i] == j else 0)
        ms = m_of_sigma(sg)
        for (i, j), v in ms.items():
            ok &= sub(Pm, (0, i), [k for k in N if k != j]).det()**2 == v
    check("m_ij(sigma) = 1{{sigma(0),sigma(i)} = N\\{j}} equals the squared 2x2 minor of P_sigma", ok)

    # phi: involution on the 18 permutations not fixing 0; m(phi(pi)) = lower-right block of P_pi
    ok_inv = ok_blk = True
    for pi in PERMS:
        if pi[0] == 0:
            continue
        ok_inv &= phi(phi(pi)) == pi and phi(pi)[0] != 0
        ms = m_of_sigma(phi(pi))
        ok_blk &= all(ms[(i, j)] == (1 if pi[i] == j else 0) for i in N for j in N)
    check("phi is an involution on {pi : pi(0) != 0} (18 perms)", ok_inv)
    check("m(phi(pi)) = lower-right 3x3 block of P_pi for all 18 pi with pi(0) != 0", ok_blk)

    # face {W00 = 0} of the Birkhoff polytope: affine dimension 8
    vecs = [sp.Matrix([1 if p[i] == j else 0 for i in range(4) for j in range(4)]) for p in PERMS if p[0] != 0]
    aff = sp.Matrix.hstack(*[v - vecs[0] for v in vecs[1:]])
    check(f"face W00=0 has affine dimension 8 (computed {aff.rank()})", aff.rank() == 8)

    # instances: exact unitaries, law construction, (2), (4), (1) as a polynomial identity
    a_full = [d, a1, a2, a3]
    c_full = [d, c1, c2, c3]
    nterms = []
    Us = [gaussian_unitary(rng) for _ in range(5)]
    Us.append(real_rational_SO4(rng))
    # s = 0 case and s = 1 case
    V = gaussian_unitary(rng, n=3)
    Us.append(sp.diag(1, V))
    Us.append(sp.Matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]) * sp.diag(1, 1, 1, 1))
    for idx, U in enumerate(Us):
        assert is_unitary(U)
        law, m, s = shared_zero_law(U)
        okm = True
        for i in N:
            okm &= sum(m[(i, j)] for j in N) == s - abs2(U[i, 0])
        for j in N:
            okm &= sum(m[(i, j)] for i in N) == s - abs2(U[0, j])
        okm &= sum(m.values()) == 2 * s
        check(f"U#{idx} (s={s}): marginal identities (2)", okm)
        okl = all(w >= 0 for w in law.values()) and sum(law.values()) == 1
        okl &= sum(w for sg, w in law.items() if sg[0] == 0) == (abs2(U[0, 0]))
        okl &= all(sum(w * m_of_sigma(sg)[k] for sg, w in law.items()) == m[k] for k in m)
        check(f"U#{idx}: law is a probability law with sum_(sigma(0)=0) p = |u00|^2 and (4)", okl)
        nterms.append(len([w for w in law.values() if w > 0]))
        lhs = det_poly(U, a_full, c_full)
        check(f"U#{idx}: (1) holds as a polynomial identity in d,a1..a3,c1..c3 ({nterms[-1]} weights)",
              sp.expand(lhs - law_poly(law, a_full, c_full)) == 0)
    print(f"[INFO] number of nonzero weights from the greedy exact Birkhoff decomposition: {nterms} "
          "(the bound <= 10 is from Caratheodory, PROOF-READ; greedy decompositions need not attain it)")

    # cross pairs other than (0,0), via relabeling (used by the single-collision note)
    U = Us[0]
    ok = True
    for i0 in range(4):
        for j0 in range(4):
            law = law_for_cross_pair(U, i0, j0)
            aa = list(sp.symbols("x0:4"))
            cc = list(sp.symbols("y0:4"))
            aa[i0] = d
            cc[j0] = d
            ok &= sp.expand(det_poly(U, aa, cc) - law_poly(law, aa, cc)) == 0
    check("relabeling: a law exists for every shared pair (i0, j0), 16 pairs, U#0", ok)

    # ---------------- Theorem 2 (order >= 5 failure) -------------------------------
    O = sp.Matrix([[1, 4, 4, 4], [-4, 1, -4, 4], [-4, 4, 1, -4], [-4, -4, 4, 1]]) / 7
    check("Thm 2: the displayed O is orthogonal", O * O.T == sp.eye(4))
    print(f"[INFO] det O = {O.det()}")
    B = (0, 1, 2, 3)

    def dIJ(M, I, J):
        return sub(M, I, J).det()**2

    s1 = dIJ(O, (0,), (0,)) + sum(dIJ(O, tuple(x for x in B if x != i), tuple(x for x in B if x != i)) for i in (1, 2, 3))
    vals = {"d00+sum": s1, "d01,02": dIJ(O, (0, 1), (0, 2)), "d01,13": dIJ(O, (0, 1), (1, 3)),
            "d02,01": dIJ(O, (0, 2), (0, 1)), "d02,23": dIJ(O, (0, 2), (2, 3)), "dBB": dIJ(O, B, B)}
    check(f"Thm 2: minor values {vals} match 4/49, 144/2401, 144/2401, 400/2401, 400/2401, 1",
          [vals[k] for k in vals] == [sp.Rational(4, 49), sp.Rational(144, 2401), sp.Rational(144, 2401),
                                      sp.Rational(400, 2401), sp.Rational(400, 2401), 1])
    GU = sp.Rational(1, 2)*(s1 + vals["d01,02"] + vals["d01,13"] - vals["d02,01"] - vals["d02,23"]) \
        + sp.Rational(1, 2)*(1 - vals["dBB"])
    check(f"Thm 2: G(U) = {GU} = -158/2401", GU == sp.Rational(-158, 2401))

    def G_perm(sg):
        def dd(I, J):
            return 1 if set(sg[i] for i in I) == set(J) else 0
        F = sp.Rational(1, 2)*(dd((0,), (0,)) + sum(dd(tuple(x for x in B if x != i), tuple(x for x in B if x != i)) for i in (1, 2, 3))
                               + dd((0, 1), (0, 2)) + dd((0, 1), (1, 3)) - dd((0, 2), (0, 1)) - dd((0, 2), (2, 3)))
        return F + sp.Rational(1, 2)*(1 - dd(B, B))
    for n in (4, 5, 6, 7):
        mn = min(G_perm(sg) for sg in itertools.permutations(range(n)))
        check(f"Thm 2: min over S_{n} of G(sigma) = {mn} >= 0", mn >= 0)
    # the surviving-monomial claim: 0 in I and 4 not in J for every minor used
    used = [((0,), (0,)), ((0, 2, 3), (0, 2, 3)), ((0, 1, 3), (0, 1, 3)), ((0, 1, 2), (0, 1, 2)),
            ((0, 1), (0, 2)), ((0, 1), (1, 3)), ((0, 2), (0, 1)), ((0, 2), (2, 3)), (B, B)]
    check("Thm 2: every minor used has 0 in I and 4 not in J", all(0 in I and 4 not in J for I, J in used))
    # n = 5: coefficient of prod_{i notin I} a_i prod_{j in J} c_j in det(diag(a) - U diag(c) U^*)
    # at a0 = c4 = 0 is (-1)^|J| |det U[I,J]|^2 -- replay for U = O + I_1 on the used minors
    U5 = sp.diag(O, 1)
    a5 = list(sp.symbols("p0:5"))
    c5 = list(sp.symbols("q0:5"))
    a5[0] = 0
    c5[4] = 0
    D5 = sp.Poly(sp.expand((sp.diag(*a5) - U5 * sp.diag(*c5) * U5.T).det(method="berkowitz")),
                 *[v for v in a5 + c5 if v != 0])
    ok = True
    for I, J in used:
        mono = sp.prod([a5[i] for i in range(5) if i not in I]) * sp.prod([c5[j] for j in J])
        coeff = D5.coeff_monomial(mono)
        ok &= coeff == (-1)**len(J) * dIJ(U5, I, J)
    check("Thm 2 (n=5): each used minor is (-1)^|J| times the coefficient of its surviving monomial", ok)

    print(f"EXACT failures: {len(FAIL)}")
    return 1 if FAIL else 0


def numeric_checks(seed=7):
    """Haar-random complex unitaries, float law, float evaluation of (1) at random labels."""
    import numpy as np
    rng = np.random.default_rng(seed)
    worst = 0.0
    for _ in range(300):
        Z = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
        Q, R = np.linalg.qr(Z)
        U = Q * (np.diag(R) / abs(np.diag(R)))
        s = 1 - abs(U[0, 0])**2
        m = np.zeros((3, 3))
        for i in N:
            for j in N:
                cols = [k for k in N if k != j]
                m[i-1, j-1] = abs(np.linalg.det(U[np.ix_([0, i], cols)]))**2
        W = np.zeros((4, 4))
        W[0, 1:] = abs(U[0, 1:])**2 / s
        W[1:, 0] = abs(U[1:, 0])**2 / s
        W[1:, 1:] = m / s
        law = {}
        Wr = W.copy()
        for _k in range(30):
            if Wr.sum() < 1e-13:
                break
            cands = [p for p in PERMS if all(Wr[i, p[i]] > 1e-14 for i in range(4))]
            p = max(cands, key=lambda p: min(Wr[i, p[i]] for i in range(4)))
            mn = min(Wr[i, p[i]] for i in range(4))
            law[phi(p)] = law.get(phi(p), 0) + s * mn
            for i in range(4):
                Wr[i, p[i]] -= mn
        law[(0, 1, 2, 3)] = law.get((0, 1, 2, 3), 0) + 1 - s
        for _r in range(5):
            dd = rng.normal() + 1j*rng.normal()
            a = np.r_[dd, rng.normal(size=3) + 1j*rng.normal(size=3)]
            c = np.r_[dd, rng.normal(size=3) + 1j*rng.normal(size=3)]
            lhs = np.linalg.det(np.diag(a) - U @ np.diag(c) @ U.conj().T)
            rhs = sum(w * np.prod([a[i] - c[sg[i]] for i in range(4)]) for sg, w in law.items())
            worst = max(worst, abs(lhs - rhs) / max(1, abs(lhs)))
    check(f"(1) on 300 Haar unitaries x 5 random label sets (seed {seed}): max rel error {worst:.1e}",
          worst < 1e-8, level="NUMERIC")


if __name__ == "__main__":
    rc = main()
    numeric_checks()
    sys.exit(rc)
