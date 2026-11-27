"""Checker A: all principal minors of (zeta^(ij))_{0<=i,j<14} over F_1331, both primes above 11."""
import functools, itertools, json, os, random, sys
P = 11
random.seed(12345)

# --- F_1331 = F_11[x]/(x^3 + a x^2 + b x + c); elements are tuples (e0,e1,e2) ---
def find_cubic():
    for a, b, c in itertools.product(range(P), repeat=3):
        # a cubic over a field is irreducible iff it has no root
        if all((v**3 + a*v*v + b*v + c) % P for v in range(P)):
            return (c, b, a)  # (c, b, a) = low-to-high non-leading coefficients
CUB = find_cubic()
assert all((v**3 + CUB[2]*v*v + CUB[1]*v + CUB[0]) % P for v in range(P)), "cubic has a root"
ZERO, ONE = (0, 0, 0), (1, 0, 0)

def add(u, v): return tuple((a + b) % P for a, b in zip(u, v))
def sub(u, v): return tuple((a - b) % P for a, b in zip(u, v))
def mul(u, v):
    r = [0] * 5
    for i in range(3):
        if u[i]:
            for j in range(3):
                r[i + j] += u[i] * v[j]
    for d in (4, 3):  # x^d = -(c + b x + a x^2) x^(d-3)
        t = r[d] % P; r[d] = 0
        if t:
            for k in range(3):
                r[d - 3 + k] -= t * CUB[k]
    return tuple(x % P for x in r[:3])
def fpow(u, e):
    r = ONE
    while e:
        if e & 1: r = mul(r, u)
        u = mul(u, u); e >>= 1
    return r
def inv(u):
    assert u != ZERO
    return fpow(u, 1329)
def neg(u): return tuple((-a) % P for a in u)
def const(n): return (n % P, 0, 0)

def det(M):
    n = len(M)
    if n == 0: return ONE
    A = [row[:] for row in M]; d = ONE
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c] != ZERO), None)
        if p is None: return ZERO
        if p != c: A[c], A[p] = A[p], A[c]; d = neg(d)
        d = mul(d, A[c][c]); iv = inv(A[c][c])
        for r in range(c + 1, n):
            if A[r][c] != ZERO:
                f = mul(A[r][c], iv)
                A[r] = [sub(x, mul(f, y)) for x, y in zip(A[r], A[c])]
    return d

def det_cofactor(M):  # independent second method (Laplace along first row), small sizes only
    n = len(M)
    if n == 0: return ONE
    if n == 1: return M[0][0]
    t = ZERO
    for j in range(n):
        minor = [r[:j] + r[j+1:] for r in M[1:]]
        term = mul(M[0][j], det_cofactor(minor))
        t = add(t, term) if j % 2 == 0 else sub(t, term)
    return t

# --- find zeta of order exactly 14 ---
def find_g():
    for g in itertools.product(range(P), repeat=3):
        if g != ZERO and all(fpow(g, 1330 // q) != ONE for q in (2, 5, 7, 19)):
            return g
G = find_g()
assert fpow(G, 1330) == ONE
ZETA = fpow(G, 95)
assert fpow(ZETA, 14) == ONE and fpow(ZETA, 7) != ONE and fpow(ZETA, 2) != ONE
assert next(k for k in range(1, 15) if fpow(ZETA, k) == ONE) == 14
ZINV = fpow(ZETA, 13); assert mul(ZETA, ZINV) == ONE

def powers(z): return [fpow(z, k) for k in range(14)]
def sub_matrix(z, S, T=None):
    pw = powers(z); T = S if T is None else T
    return [[pw[(i * j) % 14] for j in T] for i in S]
def full(z): return sub_matrix(z, list(range(14)))
def subset(mask): return [i for i in range(14) if mask >> i & 1]
def rnd(): return tuple(random.randrange(P) for _ in range(3))
def matmul(A, B):
    return [[functools.reduce(add, (mul(A[i][k], B[k][j]) for k in range(len(B))), ZERO)
             for j in range(len(B[0]))] for i in range(len(A))]

def test(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + (" :: " + detail if detail else "")); assert ok, name

def self_tests():
    print(f"field: F_11[x]/(x^3 + {CUB[2]}x^2 + {CUB[1]}x + {CUB[0]}), irreducible (no root in F_11)")
    print(f"generator g = {G} (order 1330), zeta = g^95 = {ZETA} (order 14), zeta^-1 = {ZINV}")
    # (1) field axioms
    for _ in range(300):
        a, b, c = rnd(), rnd(), rnd()
        assert mul(a, b) == mul(b, a) and mul(mul(a, b), c) == mul(a, mul(b, c))
        assert mul(a, add(b, c)) == add(mul(a, b), mul(a, c))
        if a != ZERO: assert mul(a, inv(a)) == ONE
        assert fpow(a, 1331) == a  # Frobenius
    test("(1) commutativity/assoc/distrib/inverse + x^1331==x on 300 random triples", True)
    # (2)
    test("(2) empty det = 1", det([]) == ONE)
    sing = {(z, i): det(sub_matrix(zz, [i])) for z, zz in enumerate((ZETA, ZINV)) for i in range(14)}
    test("(2) singletons: det = zeta^(i*i), nonzero; equals 1 iff i*i=0 mod 14 (only i=0)",
         all(s != ZERO for s in sing.values()) and all((s == ONE) == (i * i % 14 == 0) for (z, i), s in sing.items()))
    # (3)
    M, Mi = full(ZETA), full(ZINV)
    prod = matmul(M, Mi)
    test("(3) M(zeta)*M(zeta^-1) == 14*I", all(prod[i][j] == (const(14) if i == j else ZERO) for i in range(14) for j in range(14)))
    test("(3) det M(zeta) != 0 and det M(zeta^-1) != 0", det(M) != ZERO and det(Mi) != ZERO)
    # (4) zero detection
    rep = [r[:] for r in full(ZETA)[:6]]; rep[3] = rep[1][:]
    test("(4) repeated rows -> det 0 (Gaussian and cofactor)", det(rep) == ZERO and det_cofactor([r[:4] for r in rep[:4]]) == ZERO)
    # F_1331* has order 1330, not divisible by 4, so the spec's N=4 example (primitive 4th root, S={0,2}) is not available.
    # Analogues in F_1331: root -1 (order 2), S={0,2}: entries (-1)^(ij) are all 1 since ij is even.
    test("(4) analogue of N=4,S={0,2}: root -1 of order 2, S={0,2} -> [[1,1],[1,1]], det 0",
         det([[fpow(const(-1), i * j) for j in (0, 2)] for i in (0, 2)]) == ZERO)
    # order-7 root w=zeta^2 used as a root on Z_14 (non-squarefree-type degeneracy): S={0,7} gives all ones.
    w = fpow(ZETA, 2)
    test("(4) order-7 root w on Z_14, S={0,7}: w^(ij)=1 for all i,j in S -> det 0",
         det([[fpow(w, i * j) for j in (0, 7)] for i in (0, 7)]) == ZERO)
    ok = True
    for _ in range(60):
        S = random.sample(range(14), random.randint(0, 6)); z = random.choice((ZETA, ZINV))
        ok &= det(sub_matrix(z, S)) == det_cofactor(sub_matrix(z, S))
    test("(4) Gaussian == cofactor on 60 random S (|S|<=6)", ok)
    # (5) Jacobi: det M_z[S^c] = eps * det M_z * 14^-|S| * det M_{z^-1}[S]
    inv14 = inv(const(14)); epss = set()
    for _ in range(300):
        S = subset(random.randrange(1 << 14)); Sc = [i for i in range(14) if i not in S]
        for z, zi in ((ZETA, ZINV), (ZINV, ZETA)):
            lhs = det(sub_matrix(z, Sc))
            rhs = mul(mul(det(full(z)), fpow(inv14, len(S))), det(sub_matrix(zi, S)))
            if lhs == rhs: epss.add(1)
            elif lhs == neg(rhs): epss.add(-1)
            else: epss.add(None)
    test("(5) Jacobi complement identity on 600 (S,root) samples", None not in epss, f"eps values seen: {sorted(e for e in epss if e)}")
    # (6) translation: det_{S+t} = det_S * prod_{s in S} z^(2ts + t^2)  (diagonal scalings)
    ok = True
    for _ in range(300):
        S = subset(random.randrange(1 << 14)); t = random.randrange(14); z = random.choice((ZETA, ZINV))
        a, b = det(sub_matrix(z, S)), det(sub_matrix(z, [(s + t) % 14 for s in S]))
        f = ONE
        for s in S: f = mul(f, fpow(z, (2 * t * s + t * t) % 14))
        ok &= mul(a, f) == b and (a == ZERO) == (b == ZERO)
    test("(6) translation: det_{S+t} = det_S * prod z^(2ts+t^2) exactly, 300 samples", ok)

def main():
    print("python", sys.version.split()[0])
    self_tests()
    counts, zeros = {}, []
    for name, z in (("zeta", ZETA), ("zeta_inv", ZINV)):
        nz = [0] * 15; zz = [0] * 15
        for mask in range(1 << 14):
            S = subset(mask); d = det(sub_matrix(z, S))
            if d == ZERO:
                zz[len(S)] += 1
                rec = {"root": name, "S": S, "zeta": z, "cubic_low_to_high_nonleading_cba": CUB}
                if len(S) <= 8: rec["cofactor_det"] = det_cofactor(sub_matrix(z, S))
                zeros.append(rec)
            else: nz[len(S)] += 1
        counts[name] = {"nonzero_by_size": nz, "zero_by_size": zz, "total": sum(nz) + sum(zz)}
        print(f"\nroot {name}: total {sum(nz)+sum(zz)}, zeros {sum(zz)}")
        print("size  total  nonzero  zero")
        for k in range(15): print(f"{k:4d} {nz[k]+zz[k]:6d} {nz[k]:7d} {zz[k]:5d}")
    out = {"field": f"F_11[x]/(x^3+{CUB[2]}x^2+{CUB[1]}x+{CUB[0]})", "generator": G, "zeta": ZETA, "zeta_inv": ZINV,
           "counts": counts, "zeros": zeros, "all_nonzero": not zeros}
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "counts.json"), "w"), indent=1)
    print("\nVERDICT:", "ALL 32768 determinants nonzero" if not zeros else f"{len(zeros)} ZEROS: {zeros}")

main()
