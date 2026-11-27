"""EXACT replay (Gaussian rationals, Python fractions only) of every computational claim in
the README: the rank-2 relation at n = 3, an exact infinitesimally rigid 4-contact spectrum
(rectangle type), and the support-gap identity for exact Cayley unitaries.

Linear trick: with v_3 = (q1, q2) and the chart for the other 2n-1 vectors,
x_sigma = prod_{i<3}(a_i + b_sigma(i)) * (q1 + q2 b_sigma(3)) is C-linear in (q1, q2), so
'Re x_sigma equal for sigma in T' is a real-linear system in 4 real unknowns.

Usage: python -I -B replay_exact.py
"""
import itertools, random
from fractions import Fraction as F

n = 3
PERMS = list(itertools.permutations(range(n)))
fails = 0


class G:  # Gaussian rational
    __slots__ = ('r', 'i')
    def __init__(s, r, i=0): s.r, s.i = F(r), F(i)
    def __add__(s, o): o = o if isinstance(o, G) else G(o); return G(s.r + o.r, s.i + o.i)
    __radd__ = __add__
    def __neg__(s): return G(-s.r, -s.i)
    def __sub__(s, o): return s + (-(o if isinstance(o, G) else G(o)))
    def __mul__(s, o): o = o if isinstance(o, G) else G(o); return G(s.r*o.r - s.i*o.i, s.r*o.i + s.i*o.r)
    __rmul__ = __mul__
    def conj(s): return G(s.r, -s.i)
    def __truediv__(s, o):
        o = o if isinstance(o, G) else G(o); d = o.r*o.r + o.i*o.i; t = s * o.conj(); return G(t.r/d, t.i/d)
    def __eq__(s, o): o = o if isinstance(o, G) else G(o); return s.r == o.r and s.i == o.i
    def __hash__(s): return hash((s.r, s.i))
    def __repr__(s): return f'({s.r})+({s.i})i'


def prod(xs):
    p = G(1)
    for x in xs: p = p * x
    return p


def rank(M):
    M = [row[:] for row in M]; r = 0
    for c in range(len(M[0])):
        piv = next((k for k in range(r, len(M)) if not M[k][c] == 0), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        for k in range(len(M)):
            if k != r and not M[k][c] == 0:
                f = M[k][c] / M[r][c]; M[k] = [x - f*y for x, y in zip(M[k], M[r])]
        r += 1
    return r


def nullspace_real(M):  # rational real matrix, returns basis of kernel
    M = [row[:] for row in M]; ncol = len(M[0]); pivs = []; r = 0
    for c in range(ncol):
        piv = next((k for k in range(r, len(M)) if M[k][c] != 0), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]; M[r] = [x / M[r][c] for x in M[r]]
        for k in range(len(M)):
            if k != r and M[k][c] != 0:
                f = M[k][c]; M[k] = [x - f*y for x, y in zip(M[k], M[r])]
        pivs.append(c); r += 1
    basis = []
    for fc in [c for c in range(ncol) if c not in pivs]:
        v = [F(0)] * ncol; v[fc] = F(1)
        for k, pc in enumerate(pivs): v[pc] = -M[k][fc]
        basis.append(v)
    return basis


def check(ok, msg):
    global fails
    fails += (not ok)
    print(('PASS ' if ok else 'FAIL ') + msg, flush=True)


gi = lambda r: G(random.randint(-r, r), random.randint(-r, r))
collinear = lambda p: all(((x - p[0]) * (p[1] - p[0]).conj()).i == 0 for x in p[2:])


def rect_type(T):
    """T = {sigma : sigma(k) != l} for some k, l."""
    return any(set(T) == {s for s in PERMS if s[k] != l} for k in range(n) for l in range(n))


def attempt(T):
    a12 = [gi(3), gi(3)]; b = [gi(3) for _ in range(n)]
    # x_sigma = c1 q1 + c2 q2 ; Re(c q) = Re c Re q - Im c Im q
    coef = {s: ((a12[0] + b[s[0]]) * (a12[1] + b[s[1]]), (a12[0] + b[s[0]]) * (a12[1] + b[s[1]]) * b[s[2]]) for s in PERMS}
    rows = []
    for s in T[1:]:
        c1, c2 = coef[s][0] - coef[T[0]][0], coef[s][1] - coef[T[0]][1]
        rows.append([c1.r, -c1.i, c2.r, -c2.i])
    ns = nullspace_real(rows)
    if len(ns) != 1: return None
    v = ns[0]; Q1, Q2 = G(v[0], v[1]), G(v[2], v[3])
    if Q2 == 0: return None
    a = a12 + [Q1 / Q2]
    z = {s: prod([a[i] + b[s[i]] for i in range(n)]) for s in PERMS}
    others = [(Q2 * (z[s] - z[T[0]])).r for s in PERMS if s not in T]
    if any(o == 0 for o in others) or not (all(o > 0 for o in others) or all(o < 0 for o in others)): return None
    theta = Q2 if others[0] > 0 else -Q2
    if len(set(a)) < n or len(set(b)) < n or collinear(a) or collinear(b): return None
    if any(ai + bj == 0 for ai in a for bj in b): return None
    return a, b, theta


def verify(name, a, b, theta):
    print(f'== {name}: a = {a}, b = {b}, theta = {theta}', flush=True)
    z = {s: prod([a[i] + b[s[i]] for i in range(n)]) for s in PERMS}
    re = {s: (theta * z[s]).r for s in PERMS}
    m = min(re.values()); T = [s for s in PERMS if re[s] == m]
    check(len(T) == 2 * n - 2, f'{name}: exactly 2n-2 = {2*n-2} permutations on the supporting line Re(theta w) = {m}: {T}')
    check(len(set(a)) == n and len(set(b)) == n and not collinear(a) and not collinear(b),
          f'{name}: spectra distinct, neither collinear')
    # Jacobian of (lam, a, b) -> lam z(a, b) at (theta, a, b); holomorphic, exact
    J = []
    for s in PERMS:
        f = [prod([a[i] + b[s[i]] for i in range(n)])]                       # d/d lam
        f += [theta * prod([a[i] + b[s[i]] for i in range(n) if i != k]) for k in range(n)]   # d/d a_k
        f += [theta * prod([a[i] + b[s[i]] for i in range(n) if s[i] != k]) for k in range(n)]  # d/d b_k
        J.append(f)
    pi = PERMS.index(T[0]); idx = [PERMS.index(s) for s in T[1:]]
    CJ = [[J[k][c] - J[pi][c] for c in range(2 * n + 1)] for k in idx]
    CJ.append([sum((J[k][c] for k in range(len(PERMS))), G(0)) - len(PERMS) * J[pi][c] for c in range(2 * n + 1)])
    rJ, rCJ = rank(J), rank(CJ)
    check(rJ == 2 * n - 2 and rCJ == rJ, f'{name}: rank dPhi = {rJ} = 2n-2 and the {2*n-2} Lambda-conditions have rank {rCJ} on it (infinitesimally rigid)')
    check(rect_type(T) == (name == 'rect'), f'{name}: contact type {"rectangle {sigma(k) != l}" if rect_type(T) else "3-cycle complement"}')
    fl = lambda t: [float(t.r), float(t.i)]
    return dict(a=[fl(t) for t in a], b=[fl(t) for t in b], theta=fl(theta),
                exact=dict(a=[repr(t) for t in a], b=[repr(t) for t in b], theta=repr(theta)), T=[list(s) for s in T])


def sgn(s):
    return (-1) ** sum(1 for i in range(n) for j in range(i + 1, n) if s[i] > s[j])


def matmul(X, Y):
    return [[sum((X[i][k] * Y[k][j] for k in range(n)), G(0)) for j in range(n)] for i in range(n)]


def det3(M):
    return sum((sgn(s) * prod([M[i][s[i]] for i in range(n)]) for s in PERMS), G(0))


def inv3(M):
    d = det3(M)
    cof = lambda i, j: (M[(i+1) % 3][(j+1) % 3] * M[(i+2) % 3][(j+2) % 3] - M[(i+1) % 3][(j+2) % 3] * M[(i+2) % 3][(j+1) % 3])
    return [[cof(j, i) / d for j in range(n)] for i in range(n)]


def cayley(rng):
    """Exact unitary (I - iH)(I + iH)^{-1}, H Hermitian with Gaussian-integer entries."""
    H = [[G(0)] * n for _ in range(n)]
    for i in range(n):
        H[i][i] = G(rng.randint(-3, 3))
        for j in range(i + 1, n):
            H[i][j] = G(rng.randint(-3, 3), rng.randint(-3, 3)); H[j][i] = H[i][j].conj()
    iH = [[G(0, 1) * H[i][j] for j in range(n)] for i in range(n)]
    Im = [[G(int(i == j)) for j in range(n)] for i in range(n)]
    P = [[Im[i][j] - iH[i][j] for j in range(n)] for i in range(n)]
    Q = [[Im[i][j] + iH[i][j] for j in range(n)] for i in range(n)]
    return matmul(P, inv3(Q))


# (1) The rank-2 relation: sum_sigma sgn(sigma) z_sigma = det(a_i + b_j) = 0 (n = 3), at random exact spectra.
random.seed(7)
ok = True
for _ in range(20):
    a, b = [gi(9) for _ in range(n)], [gi(9) for _ in range(n)]
    ok &= sum((sgn(s) * prod([a[i] + b[s[i]] for i in range(n)]) for s in PERMS), G(0)) == 0
check(ok, 'sum_sigma sgn(sigma) z_sigma(a, b) = 0 at 20 random Gaussian-integer spectra (n = 3)')

# (2) An exact rigid rectangle-type configuration (the only possible type, by (1)).
random.seed(1)
cfgs = {}
for _ in range(20000):
    T = random.sample(PERMS, 2 * n - 2)
    if not rect_type(T): continue
    r = attempt(T)
    if r:
        a, b, theta = r; cfgs['rect'] = verify('rect', a, b, theta); break
else:
    check(False, 'no rectangle configuration found')

# (3) The support-gap identity Re(theta d_U) = m + |U_kl|^2 (Re(theta z_sigma') - m) for exact unitaries.
z = {s: prod([a[i] + b[s[i]] for i in range(n)]) for s in PERMS}
m = min((theta * z[s]).r for s in PERMS)
T = [s for s in PERMS if (theta * z[s]).r == m]
k, l = next((k, l) for k in range(n) for l in range(n) if set(T) == {s for s in PERMS if s[k] != l})
s1, s2 = [s for s in PERMS if s not in T]
check((theta * z[s1]).r == (theta * z[s2]).r, f'the two non-contacts {s1}, {s2} (sigma({k+1}) = {l+1}) have equal real gaps')
rng = random.Random(3)
ok = True
for _ in range(5):
    U = cayley(rng)
    Ustar = [[U[j][i].conj() for j in range(n)] for i in range(n)]
    UU = matmul(U, Ustar)
    ok &= all(UU[i][j] == G(int(i == j)) for i in range(n) for j in range(n))
    Bp = matmul(matmul(U, [[b[i] if i == j else G(0) for j in range(n)] for i in range(n)]), Ustar)
    d = det3([[Bp[i][j] + (a[i] if i == j else G(0)) for j in range(n)] for i in range(n)])
    w = U[k][l]
    ok &= (theta * d).r == m + (w * w.conj()).r * ((theta * z[s1]).r - m)
check(ok, f'Re(theta d_U) = m + |U_{k+1}{l+1}|^2 (Re(theta z_sigma\') - m) exactly for 5 Cayley unitaries')

print(f'{fails} failures', flush=True)
raise SystemExit(1 if fails else 0)
