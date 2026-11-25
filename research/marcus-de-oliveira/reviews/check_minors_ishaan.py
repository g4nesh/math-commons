"""Independent checks for the mdo-minors review (reviewer ishaan-air).

Run from the repository root:
    python3 research/marcus-de-oliveira/reviews/check_minors_ishaan.py
Needs only the standard library and numpy. Deterministic (seed 2026). Runs in a few seconds.

Part A uses exact rational arithmetic (fractions.Fraction). Part B is floating point.
"""
import itertools, json, math, os
from fractions import Fraction as Fr
import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'minors')
PERMS4 = list(itertools.permutations(range(4)))


def det(M):  # exact Gaussian elimination over Fractions (or floats)
    M = [list(r) for r in M]; k = len(M); d = Fr(1)
    for i in range(k):
        p = next((r for r in range(i, k) if M[r][i] != 0), None)
        if p is None:
            return Fr(0)
        if p != i:
            M[i], M[p] = M[p], M[i]; d = -d
        d *= M[i][i]
        for r in range(i + 1, k):
            f = M[r][i] / M[i][i]
            for c in range(i, k):
                M[r][c] -= f * M[i][c]
    return d


def sub(M, I, J): return [[M[i][j] for j in J] for i in I]
def mul(A, B): return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def tr(A): return [list(r) for r in zip(*A)]
def pmat(p): return [[Fr(int(p[i] == j)) for j in range(4)] for i in range(4)]
def eye(n): return [[Fr(int(i == j)) for j in range(n)] for i in range(n)]


def g(U):  # exotic facet functional of report 5.1 / drury-constant.md (0-based indices), real U
    return U[0][0] ** 2 + U[2][2] ** 2 + det(sub(U, [0, 3], [2, 3])) ** 2 - det(sub(U, [0, 1], [0, 3])) ** 2


print("== A. Exact checks (fractions.Fraction)")
# A1. g >= 0 on all 24 permutation matrices
vals = [g(pmat(p)) for p in PERMS4]
print("A1 g on the 24 permutation matrices: min", min(vals), "| counts", {int(v): vals.count(v) for v in sorted(set(vals))})
# A2. Zheng's O* and the skew-conference family of drury-constant.md
Ostar = [[Fr(v, 7) for v in r] for r in [[1, 4, 4, 4], [-4, 1, 4, -4], [-4, -4, 1, 4], [-4, 4, -4, 1]]]
S = [[0, 1, 1, -1], [-1, 0, -1, -1], [-1, 1, 0, 1], [1, 1, -1, 0]]
S = [[Fr(v) for v in r] for r in S]
O17 = [[(Fr(int(i == j)) + 4 * S[i][j]) / 7 for j in range(4)] for i in range(4)]  # O(1/7) = I/7 + (4/7) S
print("A2 S^T = -S:", tr(S) == [[-x for x in r] for r in S], "| S^2 = -3I:", mul(S, S) == [[-3 * x for x in r] for r in eye(4)])
print("   O* orthogonal:", mul(tr(Ostar), Ostar) == eye(4), "| O(1/7) orthogonal:", mul(tr(O17), O17) == eye(4))
print("   g(O(1/7)) =", g(O17), "| O(1/7) == O*:", O17 == Ostar)
orbit = [g(mul(mul(pmat(a), V), pmat(b))) for V in (Ostar, tr(Ostar)) for a in PERMS4 for b in PERMS4]
print("   min of g over the 1152 images P_a O* P_b, P_a O*^T P_b:", min(orbit))
def sign_equiv(M, N):  # M = D1 N D2 for diagonal sign matrices D1, D2 ?
    if any(abs(M[i][j]) != abs(N[i][j]) for i in range(4) for j in range(4)):
        return False
    s = [[1 if M[i][j] * N[i][j] > 0 else -1 if M[i][j] != 0 else 0 for j in range(4)] for i in range(4)]
    return all(s[i][j] * s[0][0] == s[i][0] * s[0][j] for i in range(4) for j in range(4))
hits = sum(sign_equiv(mul(mul(pmat(a), V), pmat(b)), O17) for V in (Ostar, tr(Ostar)) for a in PERMS4 for b in PERMS4)
print("   number of (V, a, b) with O(1/7) = D1 P_a V P_b D2, V in {O*, O*^T}, D1, D2 sign matrices:", hits)

# A3. Proposition H5: replay the committed PMSH certificates at Zheng's O*
data = json.load(open(os.path.join(ROOT, 'moc', 'runs', 'pmsh_exact_zheng.json')))
assert data['O_star_times_7'] == [[1, 4, 4, 4], [-4, 1, 4, -4], [-4, -4, 1, 4], [-4, 4, -4, 1]]
Ts = [T for k in range(1, 5) for T in itertools.combinations(range(4), k)]  # all 15 nonempty T
for cert in data['certificates']:
    c = [Fr(x) for x in cert['c']]
    H = mul(mul(Ostar, [[c[i] if i == j else Fr(0) for j in range(4)] for i in range(4)]), tr(Ostar))
    lam = {tuple(int(ch) for ch in k): Fr(v) for k, v in cert['lambda'].items()}
    ok_minors = all(det(sub(H, T, T)) == sum(l * math.prod(c[p[t]] for t in T) for p, l in lam.items()) for T in Ts)
    # Proposition H3: marginals P(pi(t) = j) must equal O*_{tj}^2 (spectrum distinct)
    ok_marg = all(sum(l for p, l in lam.items() if p[t] == j) == Ostar[t][j] ** 2 for t in range(4) for j in range(4))
    # Theorem H1 identity at two rational a: det(D_a + H) = sum_pi lam_pi prod_i (a_i + c_pi(i))
    ok_h1 = True
    for a in ([Fr(1, 2), Fr(-3), Fr(5, 7), Fr(2)], [Fr(-1), Fr(1, 3), Fr(4), Fr(-2, 5)]):
        lhs = det([[H[i][j] + (a[i] if i == j else 0) for j in range(4)] for i in range(4)])
        ok_h1 &= lhs == sum(l * math.prod(a[i] + c[p[i]] for i in range(4)) for p, l in lam.items())
    print("A3 c =", cert['c'], "| sum 1:", sum(lam.values()) == 1, "| min weight >= 0:", min(lam.values()) >= 0,
          "| 15 principal minors:", ok_minors, "| H3 marginals = |O*|^2:", ok_marg, "| H1 identity:", ok_h1)

print("== B. Floating-point checks (numpy, seed 2026)")
rng = np.random.default_rng(2026)


def haar(n):
    Z = (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z); return Q * (np.diag(R) / abs(np.diag(R)))


def cplx(k): return rng.standard_normal(k) + 1j * rng.standard_normal(k)
def subsets(n, k): return list(itertools.combinations(range(n), k))


# B1. Proposition 1 (Cauchy-Binet expansion) and Lemma 3
e1 = e3 = 0.0
for n in range(2, 7):
    U = haar(n); a, b = cplx(n), cplx(n)
    lhs = np.linalg.det(np.diag(a) + U @ np.diag(b) @ U.conj().T)
    rhs = sum(abs(np.linalg.det(U[np.ix_(S, T)])) ** 2 * np.prod(a[[i for i in range(n) if i not in S]]) * np.prod(b[list(T)])
              for k in range(n + 1) for S in subsets(n, k) for T in subsets(n, k))
    e1 = max(e1, abs(lhs - rhs) / abs(lhs))
    for r in range(n + 1):
        Y = U[:, :r]
        for q in range(n + 1):
            for R in subsets(n, q):
                L = sum(abs(np.linalg.det(Y[list(X), :])) ** 2 for X in subsets(n, r) if set(R) <= set(X))
                e3 = max(e3, abs(L - np.linalg.det(Y[list(R), :] @ Y[list(R), :].conj().T).real))
print("B1 Prop. 1 expansion max rel err %.1e (n=2..6) | Lemma 3 max abs err %.1e (all r, R, n=2..6)" % (e1, e3))


def maxflow(cap, s, t):  # Edmonds-Karp on a dense float capacity matrix; returns flow matrix
    N = len(cap); F = np.zeros_like(cap)
    while True:
        par = [-1] * N; par[s] = s; q = [s]
        while q and par[t] < 0:
            u = q.pop(0)
            for v in range(N):
                if par[v] < 0 and cap[u, v] - F[u, v] > 1e-15:
                    par[v] = u; q.append(v)
        if par[t] < 0:
            return F
        v, d = t, np.inf
        while v != s:
            d = min(d, cap[par[v], v] - F[par[v], v]); v = par[v]
        v = t
        while v != s:
            F[par[v], v] += d; F[v, par[v]] -= d; v = par[v]


# B2. Theorem D end to end: monotone coupling X ~ P^V, Y = X + {j} ~ P^V' found by max-flow, then E z = d(U)
worst, minflow, cases = 0.0, 1.0, 0
for n in range(2, 7):
    for k in range(0, n):
        for trial in range(2):
            U = haar(n); J = list(range(k)); m = k
            Xs, Ys = subsets(n, k), subsets(n, k + 1)
            pX = [abs(np.linalg.det(U[np.ix_(X, J)])) ** 2 if k else 1.0 for X in Xs]
            pY = [abs(np.linalg.det(U[np.ix_(Y, J + [m])])) ** 2 for Y in Ys]
            N = 2 + len(Xs) + len(Ys); cap = np.zeros((N, N))
            for i, X in enumerate(Xs):
                cap[0, 1 + i] = pX[i]
                for j, Y in enumerate(Ys):
                    if set(X) <= set(Y):
                        cap[1 + i, 1 + len(Xs) + j] = 10.0
            for j in range(len(Ys)):
                cap[1 + len(Xs) + j, N - 1] = pY[j]
            F = maxflow(cap, 0, N - 1); minflow = min(minflow, F[0].sum())
            a = cplx(n); beta, gamma, delta = cplx(3)
            if trial == 1:
                gamma = beta  # coincident values are allowed
            b = np.full(n, beta); b[J] = gamma; b[m] = delta
            d = np.linalg.det(np.diag(a) + U @ np.diag(b) @ U.conj().T)
            Ez = 0
            for i, X in enumerate(Xs):
                for j, Y in enumerate(Ys):
                    f = F[1 + i, 1 + len(Xs) + j]
                    if f > 0:
                        jj = (set(Y) - set(X)).pop()
                        Ez += f * np.prod([a[t] + (gamma if t in X else delta if t == jj else beta) for t in range(n)])
            worst = max(worst, abs(d - Ez) / abs(d)); cases += 1
print("B2 Theorem D: %d cases (n=2..6, every k=|J| from 0 to n-1, incl. gamma=beta): min coupling mass %.15f, max rel err %.1e"
      % (cases, minflow, worst))


# B3. drury-constant.md: closed form on the family, minimum, cubic, Cardano form
def gf(U):
    m = lambda I, J: abs(np.linalg.det(U[np.ix_(I, J)])) ** 2
    return abs(U[0, 0]) ** 2 + abs(U[2, 2]) ** 2 + m([0, 3], [2, 3]) - m([0, 1], [0, 3])


Sf = np.array(S, dtype=float)
Of = lambda y: y * np.eye(4) + np.sqrt((1 - y * y) / 3) * Sf
closed = lambda y: 2 * y * y - 4 * np.sqrt(3) / 9 * y * (1 - y * y) ** 1.5
ys = np.linspace(-1, 1, 2001)
print("B3 closed form vs direct g on 2001 grid points: max err %.1e" % max(abs(gf(Of(y)) - closed(y)) for y in ys))
lo, hi = 0.0, 0.5
for _ in range(200):  # golden-section search for the minimum on [0, 0.5]
    x1, x2 = hi - 0.618033988749895 * (hi - lo), lo + 0.618033988749895 * (hi - lo)
    lo, hi = (lo, x2) if closed(x1) < closed(x2) else (x1, hi)
ystar = (lo + hi) / 2; gmin = closed(ystar)
realroot = lambda co: min(np.roots(co), key=lambda r: abs(r.imag)).real
Droot = realroot([4, 0, 2, -7])
card = np.cbrt(7 / 8 + 11 * np.sqrt(33) / 72) + np.cbrt(7 / 8 - 11 * np.sqrt(33) / 72)
print("   y* = %.12f, y*^2 = %.12f (root of 16u^3-24u^2+36u-1: %.12f)" % (ystar, ystar ** 2, realroot([16, -24, 36, -1])))
print("   g_min = %.13f, root of 4G^3-12G^2+14G+1 = %.13f, 1 - g_min = %.13f, root of 4D^3+2D-7 = %.13f, Cardano = %.13f"
      % (gmin, realroot([4, -12, 14, 1]), 1 - gmin, Droot, card))
print("   other half of exp(phi S/sqrt3) (negative square root): grid min %.10f >= family grid min %.10f"
      % (min(gf(y * np.eye(4) - np.sqrt((1 - y * y) / 3) * Sf) for y in ys), min(closed(y) for y in ys)))
