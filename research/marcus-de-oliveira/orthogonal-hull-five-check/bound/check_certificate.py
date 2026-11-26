"""Independent exact checker for certificate.json (Python standard library only).

It shares no code with generate_certificate.py.  It

  1. builds Phi from its definition by 1x1 and 2x2 minors, and the 30
     orthogonality relations, from scratch;
  2. reads c, the Gram blocks (v_b, G_b) and the multipliers h_k;
  3. checks that every G_b is positive definite (exact LDL^T, all pivots > 0);
  4. computes R = c - Phi - sum_b v_b^T G_b v_b - sum_k h_k g_k exactly;
  5. prints  B = c + sum |coefficients of R|.

Why B bounds Phi on O(5): at an orthogonal O every g_k vanishes, every entry
satisfies |O_ij| <= 1 (rows are unit vectors), so every monomial has absolute
value <= 1 and |R(O)| <= sum |coefficients of R|; and v^T G v >= 0.  Hence
Phi(O) = c - sum_b v_b(O)^T G_b v_b(O) - R(O) <= B.

Usage: python check_certificate.py [certificate.json]
"""
import json, os, sys, time
from fractions import Fraction

t_start = time.time()
NV = 25                       # variables x[i][j], position 5*i + j in an exponent vector
ZERO = (0,) * NV

def require(cond, msg):
    """Explicit check (not assert, so it also runs under python -O)."""
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)

def mono(pairs):
    e = [0] * NV
    for i, j in pairs:
        require(0 <= i < 5 and 0 <= j < 5, "bad variable index")
        e[5 * i + j] += 1
    return tuple(e)

def times(e, f):
    return tuple(a + b for a, b in zip(e, f))

class Poly:
    def __init__(self, terms=None):
        self.t = dict(terms or {})
    def addterm(self, e, c):
        v = self.t.get(e, 0) + c
        if v:
            self.t[e] = v
        else:
            self.t.pop(e, None)
    def __add__(self, o):
        r = Poly(self.t)
        for e, c in o.t.items():
            r.addterm(e, c)
        return r
    def scale(self, s):
        return Poly({e: s * c for e, c in self.t.items()})
    def __sub__(self, o):
        return self + o.scale(-1)
    def __mul__(self, o):
        r = Poly()
        for e, c in self.t.items():
            for f, d in o.t.items():
                r.addterm(times(e, f), c * d)
        return r

def X(i, j):
    return Poly({mono([(i, j)]): 1})
ONE = Poly({ZERO: 1})

def det_poly(rows, cols):
    """Determinant of the submatrix x[rows, cols] (len 1 or 2) by the Leibniz formula."""
    rows, cols = sorted(rows), sorted(cols)
    if len(rows) == 1:
        return X(rows[0], cols[0])
    (a, b), (c, d) = rows, cols
    return X(a, c) * X(b, d) - X(a, d) * X(b, c)

# ---- Phi, literally from the definition (indices mod 5) ----
def Phi():
    terms = []   # (sign, rows, cols)
    for v in range(5):
        terms.append((-1, {v}, {v}))
        S1, S2 = {(v + 1) % 5, (v - 1) % 5}, {(v + 2) % 5, (v - 2) % 5}
        terms += [(+1, S1, S2), (+1, S2, S1)]
        for d in (1, 2):
            terms += [(-1, {v, (v + d) % 5}, {v, (v - d) % 5}),
                      (-1, {v, (v - d) % 5}, {v, (v + d) % 5})]
    keys = [(frozenset(I), frozenset(J)) for _, I, J in terms]
    require(len(terms) == 35 and len(set(keys)) == 35, "the 35 terms must be distinct minors")
    P = Poly()
    for s, I, J in terms:
        m = det_poly(I, J)
        P = P + (m * m).scale(s)
    return P

# sanity check of Phi on permutation matrices: values 1 and -5, 60 each
def phi_at_perm(p, P):
    total = 0
    for e, c in P.t.items():
        val = 1
        for k, a in enumerate(e):
            if a and p[k // 5] != k % 5:
                val = 0
                break
        total += c * val
    return total

def relations():
    rel = {}
    for i in range(5):
        for k in range(i, 5):
            g = Poly()
            for j in range(5):
                g = g + X(i, j) * X(k, j)
            rel[("row", i, k)] = g - ONE if i == k else g
    for j in range(5):
        for l in range(j, 5):
            g = Poly()
            for i in range(5):
                g = g + X(i, j) * X(i, l)
            rel[("col", j, l)] = g - ONE if j == l else g
    return rel

def positive_definite(G):
    """Exact LDL^T without pivoting; True iff all pivots are > 0 (Sylvester)."""
    n = len(G)
    A = [row[:] for row in G]
    for k in range(n):
        p = A[k][k]
        if p <= 0:
            return False
        for i in range(k + 1, n):
            if A[i][k]:
                f = A[i][k] / p
                Ai, Ak = A[i], A[k]
                for j in range(k + 1, n):
                    if Ak[j]:
                        Ai[j] -= f * Ak[j]
    return True

def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "certificate.json")
    with open(path) as f:
        cert = json.load(f)
    P = Phi()
    import itertools
    vals = [phi_at_perm(p, P) for p in itertools.permutations(range(5))]
    print("Phi on the 120 permutation matrices: value 1 x", vals.count(1), ", value -5 x", vals.count(-5),
          ", other", 120 - vals.count(1) - vals.count(-5))
    c = Fraction(cert["c"])
    DEN = int(cert["denominator"])
    require(DEN > 0, "denominator must be positive")
    lhs = ONE.scale(c) - P

    # sum of squares part
    sos = {}
    npd = 0
    sizes = []
    for blk in cert["blocks"]:
        v = [mono(m) for m in blk["monomials"]]
        n = len(v)
        U = blk["gram_upper"]
        require(len(U) == n and all(len(U[a]) == n - a for a in range(n)), "Gram block shape")
        G = [[None] * n for _ in range(n)]
        for a in range(n):
            for b in range(a, n):
                G[a][b] = G[b][a] = Fraction(int(U[a][b - a]), DEN)
        if positive_definite(G):
            npd += 1
        sizes.append(n)
        for a in range(n):
            for b in range(n):
                if G[a][b]:
                    e = times(v[a], v[b])
                    sos[e] = sos.get(e, 0) + G[a][b]
    print("Gram blocks:", len(sizes), "; positive definite (exact):", npd, "; sizes", sorted(set(sizes)))
    require(npd == len(sizes), "a Gram block is not positive definite")

    # multiplier part
    rel = relations()
    ideal = {}
    nterms = 0
    for m in cert["multipliers"]:
        require(tuple(m["relation"]) in rel, "unknown relation label")
        g = rel[tuple(m["relation"])]
        for pairs, coef in m["terms"]:
            q = Fraction(int(coef), DEN)
            e = mono(pairs)
            nterms += 1
            for f, d in g.t.items():
                k = times(e, f)
                ideal[k] = ideal.get(k, 0) + q * d
    print("multiplier terms:", nterms)

    # remainder
    R = dict(lhs.t)
    for e, q in sos.items():
        R[e] = R.get(e, 0) - q
    for e, q in ideal.items():
        R[e] = R.get(e, 0) - q
    R = {e: q for e, q in R.items() if q}
    L1 = sum(abs(q) for q in R.values())
    B = c + L1
    print("c =", c, "=", float(c))
    print("remainder R: nonzero coefficients", len(R), "; sum of |coefficients| =", float(L1),
          "; largest |coefficient| =", float(max(map(abs, R.values()))) if R else 0)
    print("bound B = c + |R|_1 =", float(B))
    print("B < 37/25 :", B < Fraction(37, 25))
    print("B <= 13/10:", B <= Fraction(13, 10))
    for k in (Fraction(11371, 10000),):
        print(f"B <= {k} :", B <= k)
    print(f"checked in {time.time() - t_start:.1f}s")

if __name__ == "__main__":
    main()
