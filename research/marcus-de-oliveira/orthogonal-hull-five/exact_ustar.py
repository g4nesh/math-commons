"""Exact check: the matrix V below is unitary and Phi(V) = 37/25.  Arithmetic in Q(sqrt6, i sqrt15) with Fractions.
Standard library only.
An element is (a, b, c, d) = a + b*s + c*t + d*s*t with s = sqrt(6), t = i*sqrt(15); s^2 = 6, t^2 = -15; conj(t) = -t."""
import itertools
from fractions import Fraction as Fr
def terms():                 # the definition of Phi (README.md); rebuilt independently in verify_certificate.py and phi.py
    out = [((v,), (v,), -1) for v in range(5)]
    for v in range(5):
        a = tuple(sorted(((v + 1) % 5, (v - 1) % 5))); b = tuple(sorted(((v + 2) % 5, (v - 2) % 5))); out += [(a, b, 1), (b, a, 1)]
        for d in (1, 2):
            s_ = tuple(sorted((v, (v + d) % 5))); t_ = tuple(sorted((v, (v - d) % 5))); out += [(s_, t_, -1), (t_, s_, -1)]
    return out
TERMS = terms()
class K:
    __slots__ = ('v',)
    def __init__(self, a=0, b=0, c=0, d=0): self.v = (Fr(a), Fr(b), Fr(c), Fr(d))
    def __add__(self, o): return K(*[x + y for x, y in zip(self.v, o.v)])
    def __sub__(self, o): return K(*[x - y for x, y in zip(self.v, o.v)])
    def __mul__(self, o):
        a, b, c, d = self.v; e, f, g, h = o.v
        return K(a*e + 6*b*f - 15*c*g - 90*d*h, a*f + b*e - 15*c*h - 15*d*g, a*g + c*e + 6*b*h + 6*d*f, a*h + d*e + b*g + c*f)
    def conj(self): a, b, c, d = self.v; return K(a, b, -c, -d)
    def __eq__(self, o): return self.v == o.v
    def __repr__(self): return "K%s" % (tuple(str(x) for x in self.v),)
def E(a, c=0): return K(Fr(a, 20), 0, Fr(c, 20), 0)          # (a + c i sqrt15) / 20
S6 = K(0, Fr(4, 20), 0, 0)                                   # 4 sqrt6 / 20
V = [[E(4), S6, S6, S6, S6],
     [S6, E(-1, -1), E(-6, -2), E(-6, 2), E(9, 1)],
     [S6, E(-6, -2), E(-1, 1), E(9, -1), E(-6, 2)],
     [S6, E(-6, 2), E(9, -1), E(-1, 1), E(-6, -2)],
     [S6, E(9, 1), E(-6, 2), E(-6, -2), E(-1, -1)]]
n = 5; ZERO = K(); ONE = K(1)
ok_unitary = all(sum((V[i][k] * V[j][k].conj() for k in range(n)), ZERO) == (ONE if i == j else ZERO) for i in range(n) for j in range(n))
def det(S, T):
    if len(S) == 1: return V[S[0]][T[0]]
    return V[S[0]][T[0]] * V[S[1]][T[1]] - V[S[0]][T[1]] * V[S[1]][T[0]]
total = ZERO; vals = {}
for S, T, c in TERMS:
    d = det(S, T); m = d * d.conj(); assert m.v[2] == 0 and m.v[3] == 0; total = total + K(c) * m
    if len(S) == 2: vals[m.v[:2]] = vals.get(m.v[:2], 0) + 1
print("V V^* = I exactly:", ok_unitary)
print("Phi(V) =", total.v[0], "+", total.v[1], "* sqrt6   (expected 37/25 + 0)")
sq = sorted({(V[i][j] * V[i][j].conj()).v[:2] for i in range(n) for j in range(n)}); print("squared moduli of entries:", [str(a) for a, b in sq], "(sqrt6 parts:", [str(b) for a, b in sq], ")")
allm = {}
for S in itertools.combinations(range(n), 2):
    for T in itertools.combinations(range(n), 2):
        d = det(S, T); m = (d * d.conj()).v; assert m[1] == 0 and m[2] == 0 and m[3] == 0; allm[str(m[0])] = allm.get(str(m[0]), 0) + 1
print("squared 2x2 minors (value: count):", allm)
# Phi on the 120 permutation matrices, exactly
cnt = {}
for p in itertools.permutations(range(n)):
    val = 0
    for S, T, c in TERMS:
        dd = (1 if p[S[0]] == T[0] else 0) if len(S) == 1 else (1 if (p[S[0]], p[S[1]]) == T else 0) - (1 if (p[S[0]], p[S[1]]) == (T[1], T[0]) else 0)
        val += c * dd * dd
    cnt[val] = cnt.get(val, 0) + 1
print("Phi on the 120 permutation matrices (value: count):", cnt)
# spectral structure: V is symmetric, (V - I)(V + I)(V^2 + V + I) = 0 and tr V = 0, so the eigenvalues are 1, 1, -1, w, w^2 (w^3 = 1)
def mmul(A, B): return [[sum((A[i][k] * B[k][j] for k in range(n)), ZERO) for j in range(n)] for i in range(n)]
def madd(A, B, c=1): return [[A[i][j] + K(c) * B[i][j] for j in range(n)] for i in range(n)]
I5 = [[ONE if i == j else ZERO for j in range(n)] for i in range(n)]; V2 = mmul(V, V)
Pm = mmul(mmul(madd(V, I5, -1), madd(V, I5)), madd(madd(V2, V), I5))
ok_sym = all(V[i][j] == V[j][i] for i in range(n) for j in range(n)); ok_min = all(Pm[i][j] == ZERO for i in range(n) for j in range(n)); ok_tr = sum((V[i][i] for i in range(n)), ZERO) == ZERO
print("V symmetric: %s ; (V - I)(V + I)(V^2 + V + I) = 0: %s ; trace 0: %s" % (ok_sym, ok_min, ok_tr))
print("EXACT CHECK PASSED" if ok_unitary and ok_sym and ok_min and ok_tr and total.v == (Fr(37, 25), 0, 0, 0) else "EXACT CHECK FAILED")
