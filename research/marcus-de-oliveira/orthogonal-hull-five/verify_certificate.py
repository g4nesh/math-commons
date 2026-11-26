"""Exact verification of certificate.json.  Self-contained: standard library only, no import from the builder.

Claim checked:  for every real orthogonal 5x5 matrix O,   Phi(O) <= beta   (beta = 13/10 in the committed certificate).

The certificate gives rational Gram matrices G_b, one for each of the 225 sign-character blocks of cubic monomials, and
rational multipliers mu_k for products (monomial m_k) x (orthogonality relation r_k).  The script checks, in exact arithmetic:

 (A) Phi, expanded as a polynomial in the 25 entries, agrees with its definition by minors on rational test matrices
     (both are coded from the list phi_terms(), which is the formula of README.md);
 (B) every G_b is symmetric positive definite (all pivots of Gaussian elimination are positive);
 (C) the polynomial identity   beta - Phi(x) - sum_b v_b(x)^T G_b v_b(x)  =  (1/40) sum_{g in GAMMA} g.( sum_k mu_k m_k r_k )
     holds coefficient by coefficient, as polynomials in 25 variables;
 (D) every element of GAMMA maps each orthogonality relation to an orthogonality relation.

By (C) and (D) the right-hand side vanishes at every orthogonal matrix, so there beta - Phi = sum_b v_b^T G_b v_b, which is
>= 0 by (B)."""
import itertools, json, random, sys
from fractions import Fraction as Fr
cert = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "certificate.json")); beta = Fr(*cert["beta"])
X = lambda i, j: 5 * (i % 5) + (j % 5)                      # variable index of the entry (i, j)
def mul(p, q):
    out = {}
    for a, u in p.items():
        for b, w in q.items():
            m = tuple(sorted(a + b)); out[m] = out.get(m, 0) + u * w
    return {m: c for m, c in out.items() if c}
def add(p, q, c=1):
    out = dict(p)
    for m, w in q.items(): out[m] = out.get(m, 0) + c * w
    return {m: c_ for m, c_ in out.items() if c_}
def ev(p, M):
    tot = Fr(0)
    for m, c in p.items():
        t = c
        for v in m: t *= M[v // 5][v % 5]
        tot += t
    return tot
# Phi from its definition
def phi_terms():
    out = [((v,), (v,), -1) for v in range(5)]
    for v in range(5):
        a = tuple(sorted(((v + 1) % 5, (v - 1) % 5))); b = tuple(sorted(((v + 2) % 5, (v - 2) % 5))); out += [(a, b, 1), (b, a, 1)]
        for d in (1, 2):
            s = tuple(sorted((v, (v + d) % 5))); t = tuple(sorted((v, (v - d) % 5))); out += [(s, t, -1), (t, s, -1)]
    return out
def phi_direct(M):
    tot = Fr(0)
    for S, T, c in phi_terms():
        d = M[S[0]][T[0]] if len(S) == 1 else M[S[0]][T[0]] * M[S[1]][T[1]] - M[S[0]][T[1]] * M[S[1]][T[0]]; tot += c * d * d
    return tot
PHI = {}
for S, T, c in phi_terms():
    if len(S) == 1: d = {(X(S[0], T[0]),): Fr(1)}
    else: d = add({tuple(sorted((X(S[0], T[0]), X(S[1], T[1])))): Fr(1)}, {tuple(sorted((X(S[0], T[1]), X(S[1], T[0])))): Fr(1)}, -1)
    PHI = add(PHI, mul(d, d), c)
rnd = random.Random(1); okA = True
for _ in range(20):
    M = [[Fr(rnd.randint(-9, 9), rnd.randint(1, 5)) for _ in range(5)] for _ in range(5)]; okA &= ev(PHI, M) == phi_direct(M)
print("(A) polynomial form of Phi agrees with the definition on 20 rational matrices:", okA)
# symmetry group and relations
AFF = [[(a * i + b) % 5 for i in range(5)] for a in (1, 2, 3, 4) for b in range(5)]
GAMMA = [[5 * g[v // 5] + g[v % 5] for v in range(25)] for g in AFF] + [[5 * g[v % 5] + g[v // 5] for v in range(25)] for g in AFF]
def act(g, p): return {tuple(sorted(g[v] for v in m)): c for m, c in p.items()}
RELS = []
for i in range(5):
    for k in range(i, 5):
        r = {}
        for j in range(5): r = add(r, {tuple(sorted((X(i, j), X(k, j)))): Fr(1)})
        RELS.append(add(r, {(): Fr(-1)}) if i == k else r)
for j in range(5):
    for l in range(j, 5):
        r = {}
        for i in range(5): r = add(r, {tuple(sorted((X(i, j), X(i, l)))): Fr(1)})
        RELS.append(add(r, {(): Fr(-1)}) if j == l else r)
okD = all(any(act(g, r) == r2 for r2 in RELS) for g in GAMMA for r in RELS)
okO = True                                                       # the relations do vanish on orthogonal matrices (sanity, exact)
for (a, b, c) in ((3, 4, 5), (5, 12, 13)):
    Mx = [[Fr(0)] * 5 for _ in range(5)]
    for i in range(5): Mx[i][i] = Fr(1)
    Mx[0][0] = Mx[1][1] = Fr(a, c); Mx[0][1] = Fr(-b, c); Mx[1][0] = Fr(b, c); Mx[2][2] = Mx[4][4] = Fr(b, c); Mx[2][4] = Fr(a, c); Mx[4][2] = Fr(-a, c)
    okO &= all(ev(r, Mx) == 0 for r in RELS)
print("(D) the 40 symmetries permute the 30 orthogonality relations:", okD, "; relations vanish on rational orthogonal test matrices:", okO)
# Gram blocks
def char(m):
    r = [0] * 5; c = [0] * 5
    for v in m: r[v // 5] ^= 1; c[v % 5] ^= 1
    return tuple(r), tuple(c)
blocks = {}
for m in itertools.combinations_with_replacement(range(25), 3): blocks.setdefault(char(m), []).append(m)
theta = {(tuple(k[0]), tuple(k[1])): Fr(*t) for k, t in zip(cert["orbit_keys"], cert["theta"])}
if len(cert["orbit_keys"]) != len(cert["theta"]) or len(theta) != len(cert["theta"]): raise SystemExit("VERIFICATION FAILED: malformed certificate (orbit keys)")
def orbit_key(m1, m2):
    best = None
    for g in GAMMA:
        u = tuple(sorted(g[v] for v in m1)); w = tuple(sorted(g[v] for v in m2)); pr = (u, w) if u <= w else (w, u)
        if best is None or pr < best: best = pr
    return best
def positive_definite(Mx):
    n = len(Mx); A = [row[:] for row in Mx]
    for k in range(n):
        if A[k][k] <= 0: return False
        for i in range(k + 1, n):
            if A[i][k]:
                f = A[i][k] / A[k][k]
                for j in range(k, n): A[i][j] -= f * A[k][j]
    return True
SOS = {}; okB = True; sizes = {}
for ch, ms in blocks.items():
    s = len(ms); G = [[None] * s for _ in range(s)]
    for a in range(s):
        for b in range(a, s):
            kk = orbit_key(ms[a], ms[b])
            if kk not in theta: raise SystemExit("VERIFICATION FAILED: the certificate has no entry for a pair of monomials")
            G[a][b] = G[b][a] = theta[kk]
    okB &= positive_definite(G); sizes[s] = sizes.get(s, 0) + 1
    for a in range(s):
        for b in range(a, s):
            if G[a][b]:
                m = tuple(sorted(ms[a] + ms[b])); SOS[m] = SOS.get(m, 0) + (G[a][b] if a == b else 2 * G[a][b])
print("(B) Gram blocks (size: count) %s all positive definite:" % sizes, okB)
# identity
lhs = add(add({(): beta}, PHI, -1), {m: c for m, c in SOS.items() if c}, -1)
T = {}
for ri, m, num, den in cert["multipliers"]:
    if not (0 <= ri < len(RELS)) or not all(0 <= x < 25 for x in m): raise SystemExit("VERIFICATION FAILED: malformed multiplier")
    T = add(T, mul({tuple(sorted(m)): Fr(1)}, RELS[ri]), Fr(num, den))
rhs = {}
for g in GAMMA: rhs = add(rhs, act(g, T), Fr(1, 40))
okC = lhs == rhs
print("(C) beta - Phi - sum of squares equals the averaged combination of %d multiples of relations: %s (polynomial with %d monomials)" % (len(cert["multipliers"]), okC, len(lhs)))
ok = okA and okB and okC and okD and okO
print("VERIFIED: Phi(O) <= %s for every real orthogonal 5x5 matrix O" % beta if ok else "VERIFICATION FAILED")
if not ok: raise SystemExit(1)
