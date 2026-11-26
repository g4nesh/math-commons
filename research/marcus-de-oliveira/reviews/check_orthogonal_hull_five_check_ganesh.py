"""Independent exact re-check of ../../orthogonal-hull-five-check/bound/certificate.json  (reviewer: owner g4nesh).

Written from the statement in that folder's bound/README.md; it imports nothing from that folder and uses the standard
library only.  With x_ij the entries of a real 5x5 matrix it expands, in exact rational arithmetic,

    R(x) = c - Phi(x) - sum_b v_b(x)^T G_b v_b(x) - sum_k h_k(x) g_k(x),

checks that every G_b is positive definite, and reports B = c + (sum of |coefficients of R|).  At an orthogonal matrix the
relations g_k vanish and every monomial has absolute value at most 1, so Phi <= B there.
It also checks Phi(O) = 1 - 6 det(O) E_S(O) exactly at rational orthogonal matrices of both determinants."""
import itertools, json, os, sys
from fractions import Fraction as Fr
here = os.path.dirname(os.path.abspath(__file__)); cert = json.load(open(os.path.join(here, "..", "orthogonal-hull-five-check", "bound", "certificate.json")))
DEN = cert["denominator"]; num, den = cert["c"].split("/"); c = Fr(int(num), int(den))
V = lambda i, j: 5 * (i % 5) + (j % 5)
def mono(pairs): return tuple(sorted(V(i, j) for i, j in pairs))
def add(p, q, k=1):
    for m, w in q.items():
        nv = p.get(m, 0) + k * w
        if nv: p[m] = nv
        else: p.pop(m, None)
    return p
def mul(p, q):
    out = {}
    for a, u in p.items():
        for b, w in q.items():
            m = tuple(sorted(a + b)); nv = out.get(m, 0) + u * w
            if nv: out[m] = nv
            else: out.pop(m, None)
    return out
# Phi from the definition (indices mod 5)
terms = [((v,), (v,), -1) for v in range(5)]
for v in range(5):
    a = tuple(sorted(((v + 1) % 5, (v - 1) % 5))); b = tuple(sorted(((v + 2) % 5, (v - 2) % 5))); terms += [(a, b, 1), (b, a, 1)]
    for d in (1, 2):
        s = tuple(sorted((v, (v + d) % 5))); t = tuple(sorted((v, (v - d) % 5))); terms += [(s, t, -1), (t, s, -1)]
PHI = {}
for S, T, k in terms:
    if len(S) == 1: d = {(V(S[0], T[0]),): Fr(1)}
    else: d = add({tuple(sorted((V(S[0], T[0]), V(S[1], T[1])))): Fr(1)}, {tuple(sorted((V(S[0], T[1]), V(S[1], T[0])))): Fr(1)}, -1)
    add(PHI, mul(d, d), k)
def relation(kind, i, k):
    r = {}
    for j in range(5): add(r, {tuple(sorted((V(i, j), V(k, j)) if kind == "row" else (V(j, i), V(j, k)))): Fr(1)})
    if i == k: add(r, {(): Fr(-1)})
    return r
def positive_definite(M):
    n = len(M); A = [row[:] for row in M]
    for k in range(n):
        if A[k][k] <= 0: return False
        for i in range(k + 1, n):
            if A[i][k]:
                f = A[i][k] / A[k][k]
                for j in range(k, n): A[i][j] -= f * A[k][j]
    return True
R = add({(): c}, PHI, -1); pd = True; sizes = {}
for blk in cert["blocks"]:
    ms = [mono(m) for m in blk["monomials"]]; n = len(ms); up = blk["gram_upper"]; assert len(up) == n and all(len(up[a]) == n - a for a in range(n))
    G = [[None] * n for _ in range(n)]
    for a in range(n):
        for b in range(a, n): G[a][b] = G[b][a] = Fr(up[a][b - a], DEN)
    pd &= positive_definite(G); sizes[n] = sizes.get(n, 0) + 1
    for a in range(n):
        for b in range(a, n):
            if G[a][b]: add(R, {tuple(sorted(ms[a] + ms[b])): G[a][b] * (1 if a == b else 2)}, -1)
seen = set()
for mlt in cert["multipliers"]:
    kind, i, k = mlt["relation"]; assert kind in ("row", "col") and 0 <= i <= k <= 4 and (kind, i, k) not in seen; seen.add((kind, i, k))
    h = {}
    for m, numer in mlt["terms"]: add(h, {mono(m): Fr(numer, DEN)})
    add(R, mul(h, relation(kind, i, k)), -1)
l1 = sum(abs(w) for w in R.values()); B = c + l1
print("blocks (size: count):", dict(sorted(sizes.items())), "; all positive definite:", pd)
print("multipliers for %d distinct relations ; residual R: %d coefficients, sum of absolute values %.3e, largest %.2e" % (len(seen), len(R), float(l1), float(max(abs(w) for w in R.values()))))
print("c = %s ; B = c + ||R||_1 = %.10f ; B < 11371/10000: %s ; B < 37/25: %s" % (c, float(B), B < Fr(11371, 10000), B < Fr(37, 25)))
# the identity Phi = 1 - 6 det(O) E_S(O) at exact orthogonal matrices
def ev(p, M):
    tot = Fr(0)
    for m, w in p.items():
        t = w
        for v in m: t *= M[v // 5][v % 5]
        tot += t
    return tot
def sgn(p): return (-1) ** sum(1 for i in range(5) for j in range(i) if p[j] > p[i])
PERMS = list(itertools.permutations(range(5)))
def pmat(p): return [[Fr(1 if p[i] == j else 0) for j in range(5)] for i in range(5)]
S = [p for p in PERMS if ev(PHI, pmat(p)) == -5]; vals = {int(ev(PHI, pmat(p))) for p in PERMS}
def det5(M): return sum(sgn(p) * M[0][p[0]] * M[1][p[1]] * M[2][p[2]] * M[3][p[3]] * M[4][p[4]] for p in PERMS)
def ES(M): return sum(sgn(p) * M[0][p[0]] * M[1][p[1]] * M[2][p[2]] * M[3][p[3]] * M[4][p[4]] for p in S)
def giv(a, b, x, y, z):
    M = [[Fr(int(i == j)) for j in range(5)] for i in range(5)]; M[a][a] = M[b][b] = Fr(x, z); M[a][b] = Fr(-y, z); M[b][a] = Fr(y, z); return M
def mm(A, Bm): return [[sum(A[i][k] * Bm[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
tests = []; triples = [(3, 4, 5), (5, 12, 13), (8, 15, 17), (7, 24, 25), (20, 21, 29)]
import random; rnd = random.Random(5)
for t in range(12):
    M = pmat(PERMS[rnd.randrange(120)])
    for _ in range(6):
        a, b = rnd.sample(range(5), 2); M = mm(M, giv(a, b, *triples[rnd.randrange(5)]))
    if t % 2: M = [[-M[i][j] if j == 0 else M[i][j] for j in range(5)] for i in range(5)]
    orth = all(sum(M[i][k] * M[j][k] for k in range(5)) == (i == j) for i in range(5) for j in range(5))
    tests.append(orth and ev(PHI, M) == 1 - 6 * det5(M) * ES(M) and ev(PHI, M) <= B)
print("Phi on permutation matrices takes the values %s ; |S| = %d ; Phi = 1 - 6 det(O) E_S(O) and Phi <= B at 12 rational orthogonal matrices (6 of each determinant): %s" % (sorted(vals), len(S), all(tests)))
ok = pd and len(seen) == 30 and B < Fr(11371, 10000) and all(tests) and vals == {1, -5} and len(S) == 60
print("RE-CHECK PASSED" if ok else "RE-CHECK FAILED"); sys.exit(0 if ok else 1)
