#!/usr/bin/env python3
"""Independent exact check of the higher-order lifting report (README.md in this folder).

Written by a referee pass from the statements of the report only.  It shares no code
with ho_lib.py, validate.py, window.c, apply*.py, blind105.py or with so_lib.py and
exact_check.py of the second-order folder (none of them was opened).

What is tested
  (a)  Lemma 1':   varpi, the Galois action, Gamma_m = g_m varpi^m, in (Z/q^r)[zeta].
  (b1) Theorem A (i)-(iii): brute force over all pairs of fillings on tiny index sets.
  (b2) Theorem A (ii),(iv): window polynomial from its definition against v(D_K)
       computed q-adically; also lambda_0 = A_0 and lambda_1 = Lambda.
  (c0) Lemma C1 and the statements on the excess of a fibre, combinatorially, by brute
       force over all pairs of fillings.
  (c)  Theorem C (i)-(iii), including orders
       beyond the window; base 10 above 11, base 14 above 29, and base 4 above 5 (where
       the singular sets have d_S = 0 exactly).
  (d)  Theorem B: brute force over fillings with prescribed level sizes.
  (e)  Propositions E, F, Theorem D and the dual form modulo a prime l = 1 (mod N);
       the valuation statements of Proposition F and Theorem D q-adically.

Conventions.  All q-adic work is done for pairs (M, q) with q = 1 (mod M), so that the
residue field of a chart is F_q, W = Z_q, and eta is an integer M-th root of unity
modulo q^r.  O = Z_q[zeta] is represented modulo q^r, first as (Z/q^r)[t]/Phi_q(1+t)
(t = pi = zeta - 1) and then, after varpi has been constructed and verified there, as
(Z/q^r)[w]/(w^(q-1) + q) (w = varpi), which is the same ring in another basis.

Python 3.9, standard library only.  Usage:  python3 exact_check.py > exact_check.txt
"""
import itertools
import random
import struct
import sys
import time
from math import comb

FAILS = []


def fail(tag, info):
    FAILS.append((tag, info))
    print("  FAIL [%s] %s" % (tag, info))
    sys.stdout.flush()


# =============================================================== integer helpers
def vq_int(a, q):
    """q-adic valuation of a nonzero integer."""
    v = 0
    while a % q == 0:
        a //= q
        v += 1
    return v


def order_mod(x, p):
    k, y = 1, x % p
    while y != 1:
        y = y * x % p
        k += 1
    return k


def det_int(A):
    """Exact determinant of an integer matrix (fraction-free elimination)."""
    n = len(A)
    if n == 0:
        return 1
    A = [row[:] for row in A]
    sign, prev = 1, 1
    for k in range(n - 1):
        if A[k][k] == 0:
            sw = None
            for i in range(k + 1, n):
                if A[i][k] != 0:
                    sw = i
                    break
            if sw is None:
                return 0
            A[k], A[sw] = A[sw], A[k]
            sign = -sign
        akk, Ak = A[k][k], A[k]
        for i in range(k + 1, n):
            Ai = A[i]
            aik = Ai[k]
            for j in range(k + 1, n):
                Ai[j] = (Ai[j] * akk - aik * Ak[j]) // prev
        prev = akk
    return sign * A[n - 1][n - 1]


def det_mod(A, p):
    """Determinant modulo a prime p."""
    n = len(A)
    if n == 0:
        return 1
    A = [[x % p for x in row] for row in A]
    d = 1
    for c in range(n):
        piv = None
        for i in range(c, n):
            if A[i][c]:
                piv = i
                break
        if piv is None:
            return 0
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            d = -d
        Ac = A[c]
        pc = Ac[c]
        d = d * pc % p
        inv = pow(pc, p - 2, p)
        for i in range(c + 1, n):
            Ai = A[i]
            f = Ai[c] * inv % p
            if f:
                for j in range(c + 1, n):
                    Ai[j] = (Ai[j] - f * Ac[j]) % p
    return d % p


def inv_mat_mod(A, m):
    """Inverse of a square matrix modulo m (a prime or prime power); raises if singular."""
    n = len(A)
    B = [[x % m for x in row] + [1 if i == j else 0 for j in range(n)] for i, row in enumerate(A)]
    for c in range(n):
        piv = None
        inv = None
        for i in range(c, n):
            try:
                inv = pow(B[i][c], -1, m)
                piv = i
                break
            except ValueError:
                continue
        if piv is None:
            raise ZeroDivisionError("singular matrix")
        B[c], B[piv] = B[piv], B[c]
        Bc = [x * inv % m for x in B[c]]
        B[c] = Bc
        for i in range(n):
            if i != c and B[i][c]:
                f = B[i][c]
                B[i] = [(x - f * y) % m for x, y in zip(B[i], Bc)]
    return [row[n:] for row in B]


def vdm(xs, mod):
    """V(x) = prod_{i<j} (x_j - x_i) modulo mod."""
    p = 1
    for j in range(len(xs)):
        for i in range(j):
            p = p * (xs[j] - xs[i]) % mod
    return p


# ============================================ the ring (Z/q^r)[t]/Phi_q(1+t), t = pi
class TRing:
    def __init__(self, q, r):
        self.q, self.r, self.Q, self.n = q, r, q ** r, q - 1
        # Phi_q(1+t) = sum_{k=1}^{q} C(q,k) t^(k-1), so t^(q-1) = -sum_{k=1}^{q-1} C(q,k) t^(k-1)
        self.red = [(-comb(q, k)) % self.Q for k in range(1, q)]
        self.one = [1] + [0] * (self.n - 1)
        self.t = [0, 1] + [0] * (self.n - 2)

    def add(self, a, b):
        Q = self.Q
        return [(x + y) % Q for x, y in zip(a, b)]

    def sub(self, a, b):
        Q = self.Q
        return [(x - y) % Q for x, y in zip(a, b)]

    def smul(self, s, a):
        Q = self.Q
        return [s * x % Q for x in a]

    def mul(self, a, b):
        n, Q = self.n, self.Q
        c = [0] * (2 * n - 1)
        for i, ai in enumerate(a):
            if ai:
                for j, bj in enumerate(b):
                    c[i + j] += ai * bj
        red = self.red
        for k in range(2 * n - 2, n - 1, -1):
            ck = c[k] % Q
            if ck:
                base = k - n
                for j in range(n):
                    c[base + j] += ck * red[j]
        return [x % Q for x in c[:n]]

    def pw(self, a, e):
        res, base = self.one, a
        while e:
            if e & 1:
                res = self.mul(res, base)
            base = self.mul(base, base)
            e >>= 1
        return res

    def val(self, a):
        """pi-adic valuation: 1, pi, ..., pi^(q-2) is a Z_q-basis, the terms a_k pi^k have
        valuations k + (q-1) v_q(a_k), distinct modulo q-1.  None for 0."""
        best = None
        for k, x in enumerate(a):
            if x:
                v = k + self.n * vq_int(x, self.q)
                if best is None or v < best:
                    best = v
        return best

    def subst(self, a, s):
        """The image of a under t -> s."""
        res = [0] * self.n
        p = self.one
        for x in a:
            if x:
                res = self.add(res, self.smul(x, p))
            p = self.mul(p, s)
        return res


# ============================================ the ring (Z/q^r)[w]/(w^(q-1) + q), w = varpi
class WRing:
    def __init__(self, q, r):
        self.q, self.r, self.Q, self.n = q, r, q ** r, q - 1
        self.one = [1] + [0] * (self.n - 1)
        self.zero = [0] * self.n
        self.prec = r * self.n          # everything is known modulo w^prec

    def mul(self, a, b):
        n, Q, q = self.n, self.Q, self.q
        c = [0] * (2 * n)
        for i in range(n):
            ai = a[i]
            if ai:
                for j in range(n):
                    c[i + j] += ai * b[j]
        return [(c[k] - q * c[k + n]) % Q for k in range(n)]

    def submul(self, x, m, p):
        """x - m*p"""
        n, Q, q = self.n, self.Q, self.q
        c = [0] * (2 * n)
        for i in range(n):
            mi = m[i]
            if mi:
                for j in range(n):
                    c[i + j] += mi * p[j]
        return [(x[k] - c[k] + q * c[k + n]) % Q for k in range(n)]

    def add(self, a, b):
        Q = self.Q
        return [(x + y) % Q for x, y in zip(a, b)]

    def sub(self, a, b):
        Q = self.Q
        return [(x - y) % Q for x, y in zip(a, b)]

    def smul(self, s, a):
        Q = self.Q
        return [s * x % Q for x in a]

    def val(self, a):
        best = None
        for k, x in enumerate(a):
            if x:
                if best is not None and k >= best:
                    break
                v = k + self.n * vq_int(x, self.q)
                if best is None or v < best:
                    best = v
        return best

    def shift(self, a, v):
        """a / w^v, for an element whose stored valuation is at least v."""
        q, Q, n = self.q, self.Q, self.n
        d, s = divmod(v, n)
        if d:                                   # w^(n d) = (-q)^d
            qd = q ** d
            sg = -1 if d & 1 else 1
            for x in a:
                assert x % qd == 0
            a = [(sg * (x // qd)) % Q for x in a]
        else:
            a = a[:]
        for _ in range(s):                      # a_0/w = -(a_0/q) w^(n-1)
            a0 = a[0]
            assert a0 % q == 0
            a = a[1:] + [(-(a0 // q)) % Q]
        return a

    def inv(self, u):
        """Inverse of a unit, by Newton iteration."""
        Q = self.Q
        y = [pow(u[0], -1, Q)] + [0] * (self.n - 1)
        k = 1
        while k < self.prec:
            uy = self.mul(u, y)
            two = [(-x) % Q for x in uy]
            two[0] = (two[0] + 2) % Q
            y = self.mul(y, two)
            k *= 2
        return y


def padic_det_val(rows, R):
    """v(det) for a matrix over R = WRing, by elimination with pivots of least valuation.
    The entries are known modulo w^A; a pivot of valuation v < A is known exactly as far as
    its valuation goes, and costs v units of absolute precision in the rest of the matrix.
    Returns (v, True) when exact, or (lower bound, False) when the precision ran out."""
    N = len(rows)
    rows = [row[:] for row in rows]
    A, total = R.prec, 0
    val, mul, sm, shift = R.val, R.mul, R.submul, R.shift
    for c in range(N):
        best, bv = None, None
        for i in range(c, N):
            v = val(rows[i][c])
            if v is not None and v < A and (bv is None or v < bv):
                best, bv = i, v
                if v == 0:
                    break
        if best is None:
            return (total + A, False)
        if best != c:
            rows[c], rows[best] = rows[best], rows[c]
        v = bv
        total += v
        rp = rows[c]
        pinv = R.inv(shift(rp[c], v))
        A -= v
        for i in range(c + 1, N):
            ri = rows[i]
            x = ri[c]
            if not any(x):
                continue
            m = mul(shift(x, v) if v else x, pinv)
            for j in range(c + 1, N):
                ri[j] = sm(ri[j], m, rp[j])
    return (total, True)


def det_laplace(A, mul, addsub, zero, one):
    """Division-free determinant over a commutative ring (row expansion, memo on column
    sets).  addsub(x, y, s) returns x + s*y for s = +1 or -1."""
    n = len(A)
    memo = {0: one}

    def rec(mask):
        got = memo.get(mask)
        if got is not None:
            return got
        k = bin(mask).count("1")
        row = k - 1                      # rows 0..k-1 against the columns in mask
        acc, pos = zero, 0
        for c in range(n):
            if mask >> c & 1:
                sgn = 1 if (k - 1 - pos) % 2 == 0 else -1
                acc = addsub(acc, mul(A[row][c], rec(mask & ~(1 << c))), sgn)
                pos += 1
        memo[mask] = acc
        return acc

    return rec((1 << n) - 1)


# ================================================== Z_q[zeta] modulo q^r, from scratch
class Cyclo:
    def __init__(self, q, r):
        self.q, self.r = q, r
        Q = self.Q = q ** r
        T = self.T = TRing(q, r)
        n = q - 1
        self.teich = [pow(x, q ** (r - 1), Q) for x in range(q)]
        # Gamma_m from the defining system zeta^x = sum_m Gamma_m [x]^m   ([x]^0 = 1 also for x = 0)
        V = [[1] + [pow(self.teich[x], m, Q) for m in range(1, q)] for x in range(q)]
        Vinv = inv_mat_mod(V, Q)
        self.zeta_t = T.add(T.one, T.t)
        zp = [T.one]
        for _ in range(1, q):
            zp.append(T.mul(zp[-1], self.zeta_t))
        self.zeta_pows_t = zp
        self.Gamma = []
        for m in range(q):
            g = [0] * n
            for x in range(q):
                g = T.add(g, T.smul(Vinv[m][x], zp[x]))
            self.Gamma.append(g)
        # varpi = pi*u with u^(q-1) * c = 1, where pi^(q-1) = -q*c, c = sum_k (C(q,k)/q) t^(k-1)
        self.c = [(comb(q, k) // q) % Q for k in range(1, q)]
        self.u = self.newton_u(T.one)
        self.varpi_t = T.mul(T.t, self.u)
        cols = [T.one]
        for _ in range(1, n):
            cols.append(T.mul(cols[-1], self.varpi_t))
        self.varpi_pows_t = cols
        B = [[cols[k][i] for k in range(n)] for i in range(n)]      # varpi-basis -> t-basis
        self.Binv = inv_mat_mod(B, Q)
        # the same ring in the varpi-basis
        W = self.W = WRing(q, r)
        self.zeta_w = self.to_w(self.zeta_t)
        zw = [W.one]
        for _ in range(1, q):
            zw.append(W.mul(zw[-1], self.zeta_w))
        self.zeta_pows_w = zw
        # g_m: Gamma_m = g_m varpi^m  (filled in and verified by check (a) / g_values)
        self.g = None

    def newton_u(self, start):
        T, q, Q = self.T, self.q, self.Q
        u = start
        inv = pow(q - 1, -1, Q)
        for _ in range(200):
            d = T.sub(T.mul(T.pw(u, q - 1), self.c), T.one)
            if not any(d):
                return u
            u = T.mul(u, T.sub(T.one, T.smul(inv, d)))
        raise RuntimeError("Newton iteration for varpi did not converge")

    def to_w(self, a):
        Q = self.Q
        return [sum(bi * x for bi, x in zip(row, a)) % Q for row in self.Binv]

    def g_values(self):
        """g_m read off from the varpi-coordinates of Gamma_m; raises if Gamma_m is not a
        Z_q-multiple of varpi^m.  g_(q-1) is determined modulo q^(r-1) only."""
        if self.g is None:
            q, Q = self.q, self.Q
            g = []
            for m in range(q):
                co = self.to_w(self.Gamma[m])
                idx = m if m < q - 1 else 0
                if any(x for k, x in enumerate(co) if k != idx):
                    raise ArithmeticError("Gamma_%d is not a multiple of varpi^%d" % (m, m))
                if m < q - 1:
                    g.append(co[m])
                else:
                    if co[0] % q:
                        raise ArithmeticError("Gamma_(q-1) not divisible by q")
                    g.append((-(co[0] // q)) % Q)
            self.g = g
        return self.g


_CYCLO = {}


def get_cyclo(q, r):
    C = _CYCLO.get((q, r))
    if C is None:
        C = Cyclo(q, r)
        W, T, Q = C.W, C.T, C.Q
        # the varpi-ring is the t-ring: varpi^(q-1) = -q there, and zeta is a primitive q-th root
        if T.pw(C.varpi_t, q - 1) != [(-q) % Q] + [0] * (q - 2):
            raise ArithmeticError("varpi^(q-1) != -q for q=%d r=%d" % (q, r))
        s = W.zero
        for z in C.zeta_pows_w:
            s = W.add(s, z)
        if any(s) or C.zeta_w == W.one or W.mul(C.zeta_pows_w[-1], C.zeta_w) != W.one:
            raise ArithmeticError("zeta is not a primitive q-th root in the varpi-ring")
        _CYCLO[(q, r)] = C
    return C


# ================================================================ index sets
def fibres_of(K, M):
    fib = [[] for _ in range(M)]
    for (a, b) in K:
        fib[a].append(b)
    for f in fib:
        f.sort()
    return fib


def points_of(fib):
    return [(a, b) for a, f in enumerate(fib) for b in f]


def level_sets(fib, q):
    """T_m = {a : c_a > m} for m = 0..q-1 (sorted tuples)."""
    return [tuple(a for a, f in enumerate(fib) if len(f) > m) for m in range(q)]


def e_of(fib):
    return sum(comb(len(f), 2) for f in fib)


class Chart:
    """(M, q) with q = 1 (mod M) and eta an element of order M in F_q."""

    def __init__(self, M, q, eta):
        assert (q - 1) % M == 0 and order_mod(eta, q) == M
        self.M, self.q, self.eta = M, q, eta
        self.ep = [pow(eta, k, q) for k in range(M)]
        self._d = {}
        self._minor = {}
        self._sing = None

    def minor(self, rowsT, colsT):
        """det F_M[I, J] in F_q, rows and columns in the given (increasing) order."""
        key = (rowsT, colsT)
        v = self._minor.get(key)
        if v is None:
            if len(rowsT) != len(colsT):
                v = 0
            else:
                ep, M = self.ep, self.M
                v = det_mod([[ep[a * b % M] for b in colsT] for a in rowsT], self.q)
            self._minor[key] = v
        return v

    def d(self, T):
        T = tuple(T)
        return self.minor(T, T)

    def singular_sets(self):
        if self._sing is None:
            out = []
            for mask in range(1, 1 << self.M):
                T = tuple(a for a in range(self.M) if mask >> a & 1)
                if self.d(T) == 0:
                    out.append(T)
            self._sing = out
        return self._sing


def charts(M, q):
    return [Chart(M, q, x) for x in range(1, q) if order_mod(x, q) == M]


# ============================== window polynomial, straight from its definition
def E_coeffs(z, q, L):
    """coefficients of E(tau z) = sum_{m=0}^{q-1} z^m tau^m / m! over F_q (first L of them)."""
    out, zp, f = [], 1, 1
    for m in range(min(q, L)):
        if m:
            f = f * m % q
        out.append(zp * pow(f, q - 2, q) % q)
        zp = zp * z % q
    return out


class Ser:
    """F_q[[tau]]/tau^L; a series is one Python int with 32-bit limbs (Kronecker packing)."""

    def __init__(self, q, L):
        assert (L + 1) * q * q < 2 ** 32
        self.q, self.L = q, L
        self.fmt = "<%dI" % L
        self.nb = 4 * L
        self.mask = (1 << (32 * L)) - 1
        self.QQ = int.from_bytes(struct.pack(self.fmt, *([q] * L)), "little")

    def pack(self, coeffs):
        c = list(coeffs)[: self.L]
        c += [0] * (self.L - len(c))
        return int.from_bytes(struct.pack(self.fmt, *c), "little")

    def unpack(self, z):
        return list(struct.unpack(self.fmt, z.to_bytes(self.nb, "little")))

    def red(self, z):
        q = self.q
        t = struct.unpack(self.fmt, (z & self.mask).to_bytes(self.nb, "little"))
        return int.from_bytes(struct.pack(self.fmt, *[c % q for c in t]), "little")

    @staticmethod
    def val(z):
        if z == 0:
            return None
        return ((z & -z).bit_length() - 1) // 32

    def inv(self, u):
        q = self.q
        y = pow(u & 0xFFFFFFFF, q - 2, q)
        k, red = 1, self.red
        while k < self.L:
            uy = red(u * y)
            y = red(y * (self.QQ - uy + 2))
            k *= 2
        return y


def series_det(rows, S):
    """det of a matrix over F_q[[tau]]/tau^L as its list of L coefficients, exact modulo
    tau^L.  Pivots of least valuation; a pivot of valuation v costs v terms of absolute
    precision in the rest, and the product of the pivots is still known modulo tau^L."""
    n, L, q = len(rows), S.L, S.q
    rows = [row[:] for row in rows]
    A, total, sign, unit = L, 0, 1, 1
    red, val, QQ = S.red, S.val, S.QQ
    for c in range(n):
        best, bv = None, None
        for i in range(c, n):
            v = val(rows[i][c])
            if v is not None and v < A and (bv is None or v < bv):
                best, bv = i, v
                if v == 0:
                    break
        if best is None:
            return [0] * L
        if best != c:
            rows[c], rows[best] = rows[best], rows[c]
            sign = -sign
        v = bv
        total += v
        if total >= L:
            return [0] * L
        rp = rows[c]
        sh = 32 * v
        pu = rp[c] >> sh
        pinv = S.inv(pu)
        unit = red(unit * pu)
        A -= v
        for i in range(c + 1, n):
            ri = rows[i]
            x = ri[c]
            if x == 0:
                continue
            mneg = QQ - red((x >> sh) * pinv)
            for j in range(c + 1, n):
                ri[j] = red(ri[j] + mneg * rp[j])
    u = S.unpack(unit)
    return [0] * total + [(sign * x) % q for x in u[: L - total]]


def window_coeffs(K, ch, L):
    """Coefficients of tau^0 .. tau^(L-1) of P_K(tau) = det( eta^(aa') E(tau b b') )."""
    M, q = ch.M, ch.q
    S = Ser(q, L)
    Ec = [E_coeffs(z, q, L) for z in range(q)]
    packed = {}
    rows = []
    for (a, b) in K:
        row = []
        for (a2, b2) in K:
            key = (a * a2 % M, b * b2 % q)
            z = packed.get(key)
            if z is None:
                h = ch.ep[key[0]]
                z = S.pack([h * c % q for c in Ec[key[1]]])
                packed[key] = z
            row.append(z)
        rows.append(row)
    return series_det(rows, S)


def pmul(a, b, q):
    if not a or not b:
        return []
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                c[i + j] += x * y
    return [x % q for x in c]


def paddsub(a, b, s, q):
    n = max(len(a), len(b))
    c = [0] * n
    for i, x in enumerate(a):
        c[i] = x
    for i, y in enumerate(b):
        c[i] = (c[i] + s * y) % q
    return c


def window_exact(K, ch):
    """The whole polynomial P_K(tau) for a tiny K, by Laplace expansion over F_q[tau]."""
    M, q = ch.M, ch.q
    Ec = [E_coeffs(z, q, q) for z in range(q)]
    A = [[[ch.ep[a * a2 % M] * c % q for c in Ec[b * b2 % q]] for (a2, b2) in K] for (a, b) in K]
    P = det_laplace(A, lambda x, y: pmul(x, y, q), lambda x, y, s: paddsub(x, y, s, q), [], [1])
    while P and P[-1] == 0:
        P.pop()
    return P


def lambdas(coeffs, fib, q):
    """lambda_j = coeff of tau^(e+j) divided by prod_a V(K_a)^2 prod_m (m!)^(-|T_m|)."""
    e = e_of(fib)
    den = 1
    for f in fib:
        den = den * vdm(f, q) ** 2 % q
    T = level_sets(fib, q)
    fact = 1
    for m in range(q):
        if m:
            fact = fact * m % q
        den = den * pow(pow(fact, q - 2, q), len(T[m]), q) % q
    dinv = pow(den, q - 2, q)
    return [c * dinv % q for c in coeffs[e:]]


# ============================================= q-adic principal minor D_K
def DK_matrix(K, M, q, C, eta_Q):
    W, Q = C.W, C.Q
    zp = C.zeta_pows_w
    ep = [pow(eta_Q, k, Q) for k in range(M)]
    cache = {}
    rows = []
    for (a, b) in K:
        row = []
        for (a2, b2) in K:
            key = (a * a2 % M, b * b2 % q)
            z = cache.get(key)
            if z is None:
                z = W.smul(ep[key[0]], zp[key[1]])
                cache[key] = z
            row.append(z)
        rows.append(row)
    return rows


def vDK(K, ch, need):
    """v(D_K) in the chart, with working precision comfortably above `need`."""
    q = ch.q
    r = need // (q - 1) + 3
    C = get_cyclo(q, r)
    eta_Q = pow(ch.eta, q ** (r - 1), q ** r)        # Teichmueller lift, an M-th root of unity
    return padic_det_val(DK_matrix(K, ch.M, q, C, eta_Q), C.W)


# ============================================= formulas of the second-order report
def Q_level(c, fib, T, ch):
    """Q_c = sum_{a,a' in T_(c-1) minus T_c} sigma_a sigma_a' Cof_(a,a') E_(a,a')  in F_q."""
    q = ch.q
    Tc1, Tc = list(T[c - 1]), list(T[c])
    shell = [a for a in Tc1 if a not in Tc]
    sig = {a: sum(fib[a]) % q for a in shell}
    tot = 0
    for a in shell:
        ra = tuple(x for x in Tc1 if x != a)
        ea = tuple(sorted(Tc + [a]))
        for a2 in shell:
            ca = tuple(x for x in Tc1 if x != a2)
            ea2 = tuple(sorted(Tc + [a2]))
            cof = ch.minor(ra, ca)
            if (Tc1.index(a) + Tc1.index(a2)) & 1:
                cof = -cof
            E = ch.minor(ea, ea2)
            if (ea.index(a) + ea2.index(a2)) & 1:
                E = -E
            tot += sig[a] * sig[a2] * cof * E
    return tot % q


def A0_and_Lambda(fib, ch):
    """A_0 = prod_m d_(T_m) and Lambda of Theorem 3, in F_q."""
    q = ch.q
    T = level_sets(fib, q)
    d = [ch.d(T[m]) for m in range(q)]
    A0 = 1
    for x in d:
        A0 = A0 * x % q
    Lam = 0
    for c in range(1, q):
        if not T[c - 1] or len(T[c - 1]) == len(T[c]):
            continue                                       # empty shell
        rest = 1
        for m in range(q):
            if m != c - 1 and m != c:
                rest = rest * d[m] % q
        if rest:
            Lam += pow(c, q - 2, q) * rest * Q_level(c, fib, T, ch)
    return A0, Lam % q


# ====================================================================== fillings
def eps_sign(filling):
    """(-1)^inv(I), inv = number of pairs (a,m), (a',m') with a < a' and m > m'."""
    inv = 0
    for i in range(len(filling)):
        Ia = filling[i]
        for j in range(i + 1, len(filling)):
            for m2 in filling[j]:
                for m in Ia:
                    if m > m2:
                        inv += 1
    return -1 if inv & 1 else 1


def schur_value(xs, expo, mod):
    """s_expo(xs) = det(x_i^(m_j)) / V(xs) modulo mod."""
    if not xs:
        return 1
    D = det_int([[pow(x, m, mod) for m in expo] for x in xs]) % mod
    return D * pow(vdm(xs, mod), -1, mod) % mod


def filling_record(support, filling, schur, q, mod):
    lev = [[] for _ in range(q)]
    wt = 0
    for a, sub in zip(support, filling):
        for m in sub:
            lev[m].append(a)
            wt += m
    return (tuple(tuple(l) for l in lev), wt, eps_sign(filling), schur % mod, tuple(filling))


def all_fillings(fib, q, xs_by_fibre, mod):
    """All fillings: (level sets, weight, sign, product of Schur values, the filling)."""
    support = [a for a, f in enumerate(fib) if f]
    per = []
    for a in support:
        per.append([(sub, schur_value(xs_by_fibre[a], sub, mod))
                    for sub in itertools.combinations(range(q), len(fib[a]))])
    out = []
    for combo in itertools.product(*per):
        sch = 1
        for c in combo:
            sch = sch * c[1] % mod
        out.append(filling_record(support, [c[0] for c in combo], sch, q, mod))
    return out


def fillings_with_sizes(fib, q, p, xs_by_fibre, mod):
    """All fillings whose vector of level sizes is p (depth-first, pruned)."""
    support = [a for a, f in enumerate(fib) if f]
    allowed = [m for m in range(q) if p[m] > 0]
    per = []
    for a in support:
        per.append([(sub, schur_value(xs_by_fibre[a], sub, mod))
                    for sub in itertools.combinations(allowed, len(fib[a]))])
    out = []
    cnt = [0] * q
    chosen = []

    def rec(i, sch):
        if i == len(support):
            if all(cnt[m] == p[m] for m in range(q)):
                out.append(filling_record(support, [c for c in chosen], sch, q, mod))
            return
        for sub, sv in per[i]:
            ok = True
            for m in sub:
                if cnt[m] >= p[m]:
                    ok = False
                    break
            if not ok:
                continue
            for m in sub:
                cnt[m] += 1
            chosen.append(sub)
            rec(i + 1, sch * sv % mod)
            chosen.pop()
            for m in sub:
                cnt[m] -= 1

    rec(0, 1)
    return out


# =============================================================== check (a)
def check_a():
    print("(a) Lemma 1': varpi, Galois action, Gamma_m = g_m varpi^m")
    total = 0
    for q, r in ((3, 10), (5, 8), (7, 7), (11, 6), (13, 6)):
        bad0 = len(FAILS)
        n_checks = 0
        C = Cyclo(q, r)
        T, Q, n = C.T, C.Q, q - 1
        tag = "a q=%d" % q
        # Teichmueller representatives
        for x in range(q):
            tx = C.teich[x]
            if tx % q != x or (x and pow(tx, q - 1, Q) != 1) or (x == 0 and tx != 0):
                fail(tag, "Teichmueller representative of %d" % x)
            n_checks += 1
        # the defining identity of the Gamma_m (Lemma 1), and their leading terms
        for x in range(q):
            acc = [0] * n
            for m in range(q):
                acc = T.add(acc, T.smul(pow(C.teich[x], m, Q) if m else 1, C.Gamma[m]))
            if acc != C.zeta_pows_t[x]:
                fail(tag, "zeta^x != sum Gamma_m [x]^m at x=%d" % x)
            n_checks += 1
        fact = 1
        for m in range(q):
            if m:
                fact *= m
            tm = T.pw(T.t, m)
            diff = T.sub(T.smul(fact, C.Gamma[m]), tm)
            v = T.val(diff)
            if not (v is None or v >= m + 1):
                fail(tag, "m! Gamma_m != pi^m mod pi^(m+1) at m=%d" % m)
            n_checks += 1
        # existence: varpi^(q-1) = -q and varpi = pi mod pi^2
        minus_q = [(-q) % Q] + [0] * (n - 1)
        if T.pw(C.varpi_t, q - 1) != minus_q:
            fail(tag, "varpi^(q-1) != -q")
        v = T.val(T.sub(C.varpi_t, T.t))
        if not (v is None or v >= 2):
            fail(tag, "varpi != pi mod pi^2")
        n_checks += 2
        # uniqueness: a different Newton start gives the same element, and the other roots
        # [xi] varpi of X^(q-1) = -q are not congruent to pi modulo pi^2
        u2 = C.newton_u(T.add(T.one, T.t))
        if u2 != C.u:
            fail(tag, "Newton iteration from another start gives a different varpi")
        n_checks += 1
        for x in range(2, q):
            other = T.smul(C.teich[x], C.varpi_t)
            if T.pw(other, q - 1) != minus_q:
                fail(tag, "([x] varpi)^(q-1) != -q at x=%d" % x)
            if T.val(T.sub(other, T.t)) != 1:
                fail(tag, "[x] varpi = pi mod pi^2 for x=%d" % x)
            n_checks += 1
        # Galois action: tau_s(varpi) = [s] varpi and tau_s(Gamma_m) = [s]^m Gamma_m
        for s in range(1, q):
            img_t = T.sub(C.zeta_pows_t[s], T.one)            # tau_s(pi) = zeta^s - 1
            if T.subst(C.varpi_t, img_t) != T.smul(C.teich[s], C.varpi_t):
                fail(tag, "tau_s(varpi) != [s] varpi at s=%d" % s)
            n_checks += 1
            for m in range(q):
                if T.subst(C.Gamma[m], img_t) != T.smul(pow(C.teich[s], m, Q), C.Gamma[m]):
                    fail(tag, "tau_s(Gamma_m) != [s]^m Gamma_m at s=%d m=%d" % (s, m))
                n_checks += 1
        # Gamma_m = g_m varpi^m with g_m in Z_q, g_m = 1/m! mod q
        try:
            g = C.g_values()
        except ArithmeticError as ex:
            fail(tag, str(ex))
            g = None
        if g is not None:
            fact = 1
            for m in range(q):
                if m:
                    fact *= m
                if m < q - 1:
                    ok = (C.Gamma[m] == T.smul(g[m], T.pw(C.varpi_t, m)))
                else:
                    ok = (C.Gamma[m] == [(-q * g[m]) % Q] + [0] * (n - 1))
                    # Gamma_(q-1) = -q/(q-1), so g_(q-1) = 1/(q-1) modulo q^(r-1)
                    ok = ok and (g[m] * (q - 1) - 1) % (q ** (r - 1)) == 0
                if not ok:
                    fail(tag, "Gamma_m != g_m varpi^m at m=%d" % m)
                if (g[m] * fact - 1) % q:
                    fail(tag, "g_m != 1/m! mod q at m=%d" % m)
                n_checks += 2
            if g[0] != 1:
                fail(tag, "g_0 != 1")
            n_checks += 1
        total += n_checks
        print("  q=%2d r=%2d (precision q^%d; g_(q-1) to q^%d): %4d assertions, %d failed"
              % (q, r, r, r - 1, n_checks, len(FAILS) - bad0))
    print("  (a) total: %d assertions" % total)


# =============================================================== check (b1)
def lambda_hat(fib, ch, C):
    """Lambda-hat_w in Z/q^r for every w, by brute force over all pairs of fillings."""
    M, q, Q = ch.M, ch.q, C.Q
    eta_Q = pow(ch.eta, q ** (C.r - 1), Q)
    g = C.g_values()
    xs = [[C.teich[b] for b in f] for f in fib]
    fl = all_fillings(fib, q, xs, Q)
    groups = {}
    for f in fl:
        groups.setdefault(tuple(len(l) for l in f[0]), []).append(f)
    Lh = [0] * (max(f[1] for f in fl) + 1)
    memo = {}

    def minor(I, J):
        key = (I, J)
        v = memo.get(key)
        if v is None:
            v = det_int([[pow(eta_Q, a * b, Q) for b in J] for a in I]) % Q
            memo[key] = v
        return v

    npairs = 0
    for sizes, lst in groups.items():
        G = 1
        for m, k in enumerate(sizes):
            if k:
                G = G * pow(g[m], k, Q) % Q
        for (levI, wI, eI, sI, _) in lst:
            acc = 0
            for (levJ, wJ, eJ, sJ, _) in lst:
                p = eJ * sJ
                for m in range(q):
                    if levI[m]:
                        p = p * minor(levI[m], levJ[m]) % Q
                        if p == 0:
                            break
                acc += p
                npairs += 1
            Lh[wI] = (Lh[wI] + eI * sI * G * acc) % Q
    return Lh, npairs, xs


def check_b1():
    print("(b1) Theorem A (i)-(iii): brute force over all pairs of fillings (tiny K)")
    rng = random.Random(20261128)
    configs = [
        (4, 5, [2, 1, 1, 0]), (4, 5, [2, 2, 0, 1]), (4, 5, [2, 0, 2, 0]), (4, 5, [3, 1, 0, 0]),
        (4, 5, [2, 2, 1, 0]), (4, 5, [1, 0, 1, 0]), (4, 5, [3, 0, 3, 0]), (4, 5, [4, 1, 0, 0]),
        (4, 5, [5, 0, 1, 0]), (2, 5, [3, 2]), (2, 3, [2, 2]), (2, 3, [3, 1]), (2, 3, [3, 2]),
        (3, 7, [2, 2, 1]), (3, 7, [3, 1, 1]), (3, 7, [1, 1, 1]), (3, 7, [4, 1, 0]), (3, 7, [6, 0, 0]),
        (3, 7, [0, 5, 1]), (6, 7, [2, 1, 1, 0, 0, 1]), (6, 7, [1, 1, 1, 1, 0, 0]), (6, 7, [0, 3, 0, 0, 2, 0]),
        (6, 7, [2, 0, 2, 0, 1, 0]), (5, 11, [2, 1, 1, 0, 0]), (5, 11, [1, 1, 1, 0, 1]), (5, 11, [3, 0, 1, 0, 0]),
        (10, 11, "S1"), (10, 11, "S2"), (10, 11, "S1"), (10, 11, "S2"),
    ]
    n_cfg = n_pairs = 0
    n_i = n_ii = n_ser = 0
    for (M, q, sizes) in configs:
        chs = charts(M, q)
        ch = rng.choice(chs)
        if isinstance(sizes, str):
            # a singular set of Z_10 as bottom level set: single layer (S1) or one fibre of size 2 (S2)
            small = [S for S in ch.singular_sets() if len(S) == min(len(x) for x in ch.singular_sets())]
            S = rng.choice(small)
            szs = [0] * M
            for a in S:
                szs[a] = 1
            if sizes == "S2":
                szs[rng.choice(S)] = 2
            sizes = szs
        r = 5
        C = get_cyclo(q, r)
        W, Q, n = C.W, C.Q, q - 1
        fib = [sorted(rng.sample(range(q), c)) for c in sizes]
        K = points_of(fib)
        tag = "b1 (%d,%d) eta=%d K=%s" % (M, q, ch.eta, K)
        Lh, npairs, xs = lambda_hat(fib, ch, C)
        n_cfg += 1
        n_pairs += npairs
        e = e_of(fib)
        if any(Lh[:e]):
            fail(tag, "Lambda-hat_w != 0 for some w < e")
        # (i) and (iii): D_K = prod V(x_a)^2 sum_w Lh_w varpi^w = prod V^2 sum_r Y_r varpi^r
        V2 = 1
        for x in xs:
            V2 = V2 * vdm(x, Q) ** 2 % Q
        Y = [0] * n
        for w, val in enumerate(Lh):
            i, k = divmod(w, n)
            Y[k] = (Y[k] + pow(-q, i, Q) * val) % Q
        pred = [V2 * y % Q for y in Y]
        eta_Q = pow(ch.eta, q ** (r - 1), Q)
        D = det_laplace(DK_matrix(K, M, q, C, eta_Q), W.mul,
                        lambda x, y, s: W.add(x, y) if s > 0 else W.sub(x, y), W.zero, W.one)
        if D != pred:
            fail(tag, "D_K != prod V^2 sum_w Lambda-hat_w varpi^w (mod q^%d)" % r)
        else:
            # valuation formula of (iii)
            vv = [k + n * vq_int(y, q) for k, y in enumerate(Y) if y]
            if (W.val(D) is None) != (not vv) or (vv and W.val(D) != min(vv)):
                fail(tag, "v(D_K) != min_r (r + (q-1) v(Y_r))")
        n_i += 1
        # (ii): the whole window polynomial
        P = window_exact(K, ch)
        V2q = 1
        for f in fib:
            V2q = V2q * vdm(f, q) ** 2 % q
        predP = [V2q * (x % q) % q for x in Lh]
        while predP and predP[-1] == 0:
            predP.pop()
        if P != predP:
            fail(tag, "P_K(tau) != prod V(K_a)^2 sum_w (Lambda-hat_w mod q) tau^w")
        n_ii += 1
        # the truncated-series routine used in (b2) and (c) against the exact polynomial
        L = e + 2 * q
        if window_coeffs(K, ch, L) != (P + [0] * L)[:L]:
            fail(tag, "series elimination disagrees with the exact polynomial determinant")
        n_ser += 1
    print("  %d index sets, %d pairs of fillings summed; A(i)/(iii) identity mod q^5: %d checked;"
          % (n_cfg, n_pairs, n_i))
    print("  A(ii) all coefficients of P_K: %d checked; series routine vs exact polynomial: %d checked"
          % (n_ii, n_ser))


# =============================================================== check (b2)
def random_fibres(M, q, sizes, rng):
    return [sorted(rng.sample(range(q), c)) for c in sizes]


def gen_b2(ch, rng, nmax):
    """A stream of index sets of several kinds for the chart ch."""
    M, q = ch.M, ch.q
    sing = [S for S in ch.singular_sets() if len(S) * 1 <= nmax - 1]
    kinds = ["sparse", "sparse", "tall", "tall"]
    if sing:
        kinds += ["level", "level", "level", "ferrers", "block"]
    kind = rng.choice(kinds)
    if kind == "sparse":
        n = rng.randint(2, nmax)
        pts = rng.sample([(a, b) for a in range(M) for b in range(q)], n)
        return kind, fibres_of(pts, M)
    if kind == "tall":
        sizes = [0] * M
        budget = nmax
        for a in rng.sample(range(M), min(M, rng.randint(1, 3))):
            c = rng.randint(2, min(q, max(2, budget - 1), 9))
            sizes[a] = c
            budget -= c
        for a in range(M):
            if sizes[a] == 0 and budget > 0 and rng.random() < 0.3:
                sizes[a] = 1
                budget -= 1
        return kind, random_fibres(M, q, sizes, rng)
    S = rng.choice(sing)
    out = [a for a in range(M) if a not in S]
    if kind == "level":
        # S as a level set of multiplicity delta at level s
        for _ in range(100):
            s = rng.randint(0, 2)
            delta = rng.randint(1, 2)
            if len(S) * (s + delta) + (s > 0) <= nmax:
                break
        else:
            s, delta = 0, 1
        sizes = [0] * M
        room = nmax - len(S) * (s + delta)
        for a in S:
            sizes[a] = s + delta
            if room > 0 and rng.random() < 0.4 and s + delta < q:
                sizes[a] += 1
                room -= 1
        if s > 0:
            for a in rng.sample(out, min(len(out), rng.randint(1, 3))):
                c = rng.randint(1, s)
                if room >= c:
                    sizes[a] = c
                    room -= c
        return kind, random_fibres(M, q, sizes, rng)
    if kind == "ferrers":
        # nested fibres with S among the level sets: window-blind by Proposition F
        perm = rng.sample(range(q), q)
        sizes = [0] * M
        top = rng.randint(1, max(1, min(3, nmax // len(S))))
        for a in S:
            sizes[a] = top
        room = nmax - top * len(S)
        if top > 1:
            for a in rng.sample(out, min(len(out), 2)):
                c = rng.randint(1, top - 1)
                if room >= c:
                    sizes[a] = c
                    room -= c
        return kind, [sorted(perm[:c]) for c in sizes]
    # block: S x B plus a few further points
    bmax = max(1, min(3, (nmax - 1) // len(S)))
    B = rng.sample(range(q), rng.randint(1, bmax))
    pts = set((a, b) for a in S for b in B)
    room = nmax - len(pts)
    rest = [(a, b) for a in range(M) for b in range(q) if (a, b) not in pts]
    pts |= set(rng.sample(rest, rng.randint(0, min(3, room))))
    return kind, fibres_of(sorted(pts), M)


def check_b2():
    print("(b2) Theorem A (ii),(iv): window polynomial from the definition vs v(D_K) (q-adic)")
    rng = random.Random(1128)
    plan = [(6, 7, 70, 22), (3, 7, 40, 16), (10, 11, 110, 20), (5, 11, 40, 18), (7, 29, 24, 12), (14, 29, 40, 12)]
    gt = dict(cases=0, win=0, blind=0, ege=0, lam01=0)
    for (M, q, ncases, nmax) in plan:
        chs = charts(M, q)
        hist = {}
        n_blind = n_ege = n_lam = 0
        kinds_seen = {}
        for _ in range(ncases):
            ch = rng.choice(chs)
            kind, fib = gen_b2(ch, rng, nmax)
            K = points_of(fib)
            if not K:
                continue
            kinds_seen[kind] = kinds_seen.get(kind, 0) + 1
            e = e_of(fib)
            tag = "b2 (%d,%d) eta=%d %s K=%s" % (M, q, ch.eta, kind, K)
            L = e + q - 1
            co = window_coeffs(K, ch, L)
            if any(co[:e]):
                fail(tag, "a coefficient of P_K below tau^e is nonzero")
            w0 = next((w for w in range(e, e + q - 1) if co[w]), None)
            v, exact = vDK(K, ch, e + 3 * (q - 1))
            if w0 is not None:
                if not exact or v != w0:
                    fail(tag, "least nonzero window order %d (e=%d) but v(D_K)=%s%s"
                         % (w0, e, v, "" if exact else " (lower bound)"))
                hist[w0 - e] = hist.get(w0 - e, 0) + 1
            else:
                if v < e + q - 1:
                    fail(tag, "no nonzero window coefficient but v(D_K)=%d < e+q-1=%d" % (v, e + q - 1))
                n_blind += 1
            if e >= q - 1:
                n_ege += 1
            # lambda_0 = A_0 and lambda_1 = Lambda (formulas of the second-order report)
            lam = lambdas(co, fib, q)
            A0, Lam = A0_and_Lambda(fib, ch)
            if lam[0] != A0 or (q > 2 and lam[1] != Lam):
                fail(tag, "lambda_0, lambda_1 = %d, %d but A_0, Lambda = %d, %d" % (lam[0], lam[1], A0, Lam))
            n_lam += 1
            gt["cases"] += 1
        gt["win"] += sum(hist.values())
        gt["blind"] += n_blind
        gt["ege"] += n_ege
        gt["lam01"] += n_lam
        print("  (%2d,%2d): %3d sets %s; v(D_K)-e = least window order: %s; window-blind with v >= e+q-1: %d; e >= q-1: %d"
              % (M, q, sum(kinds_seen.values()), dict(sorted(kinds_seen.items())),
                 dict(sorted(hist.items())), n_blind, n_ege))
    print("  (b2) total: %d sets; %d with a nonzero window coefficient (valuation agrees), %d window-blind"
          " (v >= e+q-1 confirmed), %d with e >= q-1; lambda_0 = A_0 and lambda_1 = Lambda on %d"
          % (gt["cases"], gt["win"], gt["blind"], gt["ege"], gt["lam01"]))


# =============================================================== check (c)
def block_K(M, q, S, s, delta, rng, n_up, n_low):
    """Fibres with T_m = S exactly for s <= m <= s+delta-1; n_up fibres over S of size
    exactly s+delta (if s+delta < q) and n_low fibres outside S of size exactly s (if s > 0)."""
    Sset = set(S)
    sizes = [0] * M
    top = s + delta
    if top == q:
        for a in S:
            sizes[a] = q
    else:
        n_up = max(1, min(len(S), n_up))
        cs = [top] * n_up + [rng.randint(top + 1, min(q, top + 3)) for _ in range(len(S) - n_up)]
        rng.shuffle(cs)
        for a, c in zip(S, cs):
            sizes[a] = c
    if s > 0:
        outside = [a for a in range(M) if a not in Sset]
        rng.shuffle(outside)
        n_low = max(1, min(len(outside), n_low))
        cs = [s] * n_low
        if len(outside) > n_low and s > 1 and rng.random() < 0.3:
            cs.append(rng.randint(1, s - 1))
        for a, c in zip(outside, cs):
            sizes[a] = c
    fib = random_fibres(M, q, sizes, rng)
    T = level_sets(fib, q)
    S = tuple(S)
    assert all(T[m] == S for m in range(s, s + delta))
    assert (s == 0 or T[s - 1] != S) and (s + delta == q or T[s + delta] != S)
    return fib


def block_numbers(pos, delta):
    if pos == "interior":
        m = [min(i + 1, delta - i) for i in range(delta)]
    elif pos == "bottom":
        m = [delta - i for i in range(delta)]
    else:
        m = [i + 1 for i in range(delta)]
    return sum(m), max(m)


def shells_large_enough(pos, delta, nL, nU):
    """Referee's observation (not a statement of the report): a pair of excess exactly w that
    clears every block level uses g single gains by distinct fibres of the lower shell
    T_(s-1) minus S and l single losses by distinct fibres of the upper shell S minus
    T_(s+delta); the level below (above) the block then keeps the factor d_S unless the
    shell has at least g+1 (l+1) fibres.  So lambda_w = 0 unless this returns True."""
    if pos == "bottom":
        splits = [(0, delta)]
    elif pos == "top":
        splits = [(delta, 0)]
    elif delta % 2 == 0:
        splits = [(delta // 2, delta // 2)]
    else:
        splits = [(delta // 2 + 1, delta // 2), (delta // 2, delta // 2 + 1)]
    return any((g == 0 or nL >= g + 1) and (l == 0 or nU >= l + 1) for g, l in splits)


def check_c():
    print("(c) Theorem C: a singular level set of multiplicity delta")
    rng = random.Random(110)
    stats = {}
    tot = dict(n=0, beyond=0, beyond_w=0, formula=0, formula_nz=0, small=0, small_zero=0, large=0, large_nz=0)

    def one(ch, S, pos, s, delta, n_up, n_low):
        M, q = ch.M, ch.q
        fib = block_K(M, q, S, s, delta, rng, n_up, n_low)
        K = points_of(fib)
        e = e_of(fib)
        T = level_sets(fib, q)
        w, mmax = block_numbers(pos, delta)
        assert mmax < q - 1
        if pos == "interior":
            assert s >= 1 and s + delta <= q - 1
            assert w == (delta + 1) ** 2 // 4 and mmax == (delta + 1) // 2
        else:
            assert (s == 0 and delta <= q - 1) if pos == "bottom" else (s >= 1 and s + delta == q)
            assert w == delta * (delta + 1) // 2 and mmax == delta
        jmax = w + q - 2 - mmax
        tag = "c (%d,%d) eta=%d S=%s %s s=%d delta=%d K=%s" % (M, q, ch.eta, S, pos, s, delta, K)
        co = window_coeffs(K, ch, e + jmax + 1)
        lam = lambdas(co, fib, q)                       # lambda_0 .. lambda_jmax
        st = stats.setdefault((M, q, pos, delta), dict(n=0, at_w=0, above=0, none=0, beyond=0))
        st["n"] += 1
        tot["n"] += 1
        # (i) lambda_j = 0 for j < w, and v(D_K) >= e + w
        if any(co[:e]) or any(lam[:w]):
            fail(tag, "lambda_j != 0 for some j < w = %d: %s" % (w, lam[:w]))
        j0 = next((j for j in range(len(lam)) if lam[j]), None)
        v, exact = vDK(K, ch, e + jmax + q)
        if v < e + w:
            fail(tag, "v(D_K) = %d < e + w = %d" % (v, e + w))
        # (ii) v(D_K) = e + j0 whenever j0 <= w + q - 2 - m_max
        if j0 is not None:
            if not exact or v != e + j0:
                fail(tag, "j0 = %d (range up to %d) but v(D_K) - e = %s%s"
                     % (j0, jmax, v - e, "" if exact else " (lower bound)"))
            st["at_w" if j0 == w else "above"] += 1
            if j0 > q - 2:
                st["beyond"] += 1
                tot["beyond"] += 1
                tot["beyond_w"] += 1 if j0 == w else 0
        else:
            # no nonzero lambda_j up to the end of the range.  The theorem asserts v(D_K) > e + w;
            # its proof gives more, v(D_K) >= e + w + q - 1 - m_max.
            if v <= e + w:
                fail(tag, "lambda_w = 0 but v(D_K) = e + w")
            if v < e + jmax + 1:
                fail(tag, "all lambda_j = 0 for j <= %d but v(D_K) - e = %d" % (jmax, v - e))
            st["none"] += 1
        if (exact and v == e + w) != (lam[w] != 0):
            fail(tag, "v(D_K) = e + w is not equivalent to lambda_w != 0")
        # (iii) the product formula for an interior block with delta = 2
        if pos == "interior" and delta == 2:
            rhs = pow(s * (s + 2), q - 2, q)
            for m in range(q):
                if m not in (s - 1, s, s + 1, s + 2):
                    rhs = rhs * ch.d(T[m]) % q
            rhs = rhs * Q_level(s, fib, T, ch) * Q_level(s + 2, fib, T, ch) % q
            if lam[2] != rhs:
                fail(tag, "lambda_2 = %d but the product formula gives %d" % (lam[2], rhs))
            tot["formula"] += 1
            tot["formula_nz"] += 1 if rhs else 0
        # referee's observation on the size of the shells (informational, not a pass/fail item)
        nL = len(T[s - 1]) - len(S) if s > 0 else 0
        nU = len(S) - len(T[s + delta]) if s + delta < q else 0
        if shells_large_enough(pos, delta, nL, nU):
            tot["large"] += 1
            tot["large_nz"] += 1 if lam[w] else 0
        else:
            tot["small"] += 1
            tot["small_zero"] += 1 if lam[w] == 0 else 0

    # ---- (10, 11): singular sets of size 3 and 7 in two charts, 4 and 6 in the other two
    chs = charts(10, 11)

    def pick(ch, least):
        sing = ch.singular_sets()
        k = min(len(S) for S in sing if len(S) >= least)
        return rng.choice([S for S in sing if len(S) == k])

    for ch in chs:
        for delta in (1, 2, 3, 4, 5, 6):
            need = (delta + 1) // 2 + 1
            for s in ((1, 2, 3) if delta <= 2 else (1, 2) if delta <= 4 else (1,)):
                for k in range(3 if delta <= 2 else 2 if delta <= 4 else 1):
                    one(ch, pick(ch, need), "interior", s, delta, 99, need + (k % 2))
                if delta <= 4:
                    one(ch, pick(ch, 3), "interior", s, delta, 1, need)
                    one(ch, pick(ch, 3), "interior", s, delta, 99, 1)
                    one(ch, pick(ch, 3), "interior", s, delta, 1, 1)
        one(ch, pick(ch, 5), "interior", 1, 1, 3, 2)
        one(ch, pick(ch, 5), "interior", 1, 2, 3, 3)
        one(ch, pick(ch, 3), "interior", 1, 5, 99, 1)
        for delta in (1, 2, 3, 4, 5):
            for k in range(2):
                one(ch, pick(ch, delta + 1), "bottom", 0, delta, 99, 0)
            one(ch, pick(ch, delta + 1), "bottom", 0, delta, delta, 0)
            one(ch, pick(ch, 3), "bottom", 0, delta, 1, 0)
    # top blocks contain S x Z_11; a few only
    top_chs = [ch for ch in chs if min(len(S) for S in ch.singular_sets()) == 3]
    for delta in (1, 2, 3, 4):
        one(top_chs[delta % 2], pick(top_chs[delta % 2], 3), "top", 11 - delta, delta, 0, delta + 1)
        one(top_chs[(delta + 1) % 2], pick(top_chs[(delta + 1) % 2], 3), "top", 11 - delta, delta, 0, 1)
    # ---- (4, 5): the sets {0,2}, {1,3} have d_S = 0 exactly, which the hypothesis d_S in q allows
    for ch in charts(4, 5):
        sing = ch.singular_sets()
        for (pos, s, delta) in [("interior", 1, 1), ("interior", 2, 1), ("interior", 3, 1), ("interior", 1, 2),
                                ("interior", 2, 2), ("interior", 1, 3), ("bottom", 0, 1), ("bottom", 0, 2),
                                ("bottom", 0, 3), ("top", 4, 1), ("top", 3, 2), ("top", 2, 3)]:
            for k in range(12):
                one(ch, rng.choice(sing), pos, s, delta, rng.randint(1, 2), rng.randint(1, 2))
    # ---- (14, 29): small multiplicities only (the window is long)
    for ch in charts(14, 29)[:3]:
        small = [S for S in ch.singular_sets() if len(S) == 4]
        for (pos, s, delta) in [("interior", 1, 1), ("interior", 1, 2), ("bottom", 0, 1), ("bottom", 0, 2),
                                ("bottom", 0, 3)]:
            one(ch, rng.choice(small), pos, s, delta, 99, 2)
            one(ch, rng.choice(small), pos, s, delta, 1, 1)
    for key in sorted(stats):
        st = stats[key]
        w, mmax = block_numbers(key[2], key[3])
        print("  (%2d,%2d) %-8s delta=%d w=%2d range<=%2d: %3d sets | j0=w: %3d | w<j0<=range: %3d |"
              " no j0 in range: %3d | j0 > q-2: %3d"
              % (key[0], key[1], key[2], key[3], w, w + key[1] - 2 - mmax, st["n"], st["at_w"], st["above"],
                 st["none"], st["beyond"]))
    n_j0 = sum(st["at_w"] + st["above"] for st in stats.values())
    print("  (c) total: %d sets; (i) lambda_j = 0 for j < w and v(D_K) >= e+w in all; (ii) v(D_K) = e+j0"
          " confirmed in %d," % (tot["n"], n_j0))
    print("      of which %d have j0 > q-2, beyond the window of Theorem A (%d of them with j0 = w);"
          % (tot["beyond"], tot["beyond_w"]))
    print("      no j0 in range in %d, all with v(D_K) >= e+w+q-1-m_max" % (tot["n"] - n_j0))
    print("  (c) Theorem C(iii), product formula for delta = 2 (interior): %d sets agree, product nonzero in %d"
          % (tot["formula"], tot["formula_nz"]))
    print("  (c) observation: shells too small for a term of excess w: %d sets, lambda_w = 0 in %d;"
          " shells large enough: %d sets, lambda_w != 0 in %d"
          % (tot["small"], tot["small_zero"], tot["large"], tot["large_nz"]))


# =============================================================== check (c0)
def check_c0():
    """Lemma C1 is purely combinatorial (S need not be singular): brute force over all fillings."""
    print("(c0) Lemma C1 and the excess of a fibre: brute force over all pairs of fillings")
    configs = [
        # (q, s, delta, sizes of the fibres over S, sizes of the fibres outside S)
        (5, 1, 1, [2, 2], [1, 1]), (5, 1, 2, [3, 3], [1, 1]), (5, 1, 3, [4, 4], [1, 1]), (5, 2, 2, [4, 4], [2, 1]),
        (5, 2, 1, [3, 4], [2, 2]), (5, 3, 1, [4, 5], [3, 1]), (5, 1, 2, [3, 4, 5], [1]), (5, 0, 1, [1, 1, 2], []),
        (5, 0, 2, [2, 2, 3], []), (5, 0, 3, [3, 3, 3], []), (5, 0, 4, [4, 4], []), (5, 4, 1, [5, 5], [4, 4]),
        (5, 3, 2, [5, 5], [3, 3, 2]), (5, 2, 3, [5], [2, 2, 2]), (5, 1, 4, [5, 5], [1, 1, 1]),
        (7, 1, 3, [4, 4], [1, 1]), (7, 2, 2, [4, 5], [2, 1]), (7, 1, 4, [5, 5], [1, 1]), (7, 2, 3, [5, 6], [2, 2]),
        (7, 0, 3, [3, 3, 4], []), (7, 0, 2, [2, 2, 3], []), (7, 5, 2, [7], [5, 5, 4]), (7, 4, 3, [7, 7], [4, 4, 3]),
        (7, 1, 5, [6, 6], [1, 1]), (7, 3, 2, [5, 5], [3]),
    ]
    n_fill = n_pairs = n_classes = n_attained = 0
    for (q, s, delta, S_sizes, out_sizes) in configs:
        tag = "c0 q=%d s=%d delta=%d S-fibres %s other fibres %s" % (q, s, delta, S_sizes, out_sizes)
        assert all(c >= s + delta for c in S_sizes) and all(c <= s for c in out_sizes)
        pos = "bottom" if s == 0 else "top" if s + delta == q else "interior"
        assert pos != "bottom" or not any(out_sizes)
        if pos == "interior":
            m = [min(i + 1, delta - i) for i in range(delta)]
        elif pos == "bottom":
            m = [delta - i for i in range(delta)]
        else:
            m = [i + 1 for i in range(delta)]
        cost = [sum(m[i] for i in range(delta) if P >> i & 1) for P in range(1 << delta)]
        fibs = [(True, c) for c in S_sizes] + [(False, c) for c in out_sizes if c > 0]
        per = []
        for inS, c in fibs:
            lst = []
            for sub in itertools.combinations(range(q), c):
                mask = sum(1 << k for k in sub)
                exc = sum(sub) - comb(c, 2)
                g = {k: k - bin(mask & ((1 << k) - 1)).count("1") for k in sub}          # non-members below k
                h = {l: bin(mask >> (l + 1)).count("1") for l in range(q) if l not in sub}   # members above l
                if sum(g.values()) != exc or sum(h.values()) != exc:
                    fail(tag, "excess of the fibre %s is not sum g = sum h" % (sub,))
                if (exc == 0) != (sub == tuple(range(c))) or (exc == 1) != (sub == tuple(range(c - 1)) + (c,)):
                    fail(tag, "excess 0 / 1 are not I_0 / the single raise: %s" % (sub,))
                P = ev = 0
                for i in range(delta):
                    lvl = s + i
                    if inS and lvl not in sub:              # a loss
                        if h[lvl] < delta - i:
                            fail(tag, "loss at level s+%d of cost %d < delta - i" % (i, h[lvl]))
                        P |= 1 << i
                        ev += delta - i
                    if (not inS) and lvl in sub:            # a gain
                        if g[lvl] < i + 1:
                            fail(tag, "gain at level s+%d of cost %d < i + 1" % (i, g[lvl]))
                        P |= 1 << i
                        ev += i + 1
                if exc < ev:
                    fail(tag, "excess of a fibre below the cost of its gains/losses: %s" % (sub,))
                lst.append((sum(1 << (4 * k) for k in sub), exc, P))
            per.append(lst)
        groups = {}
        for combo in itertools.product(*per):
            key = exc = P = 0
            for (kk, ee, pp) in combo:
                key += kk
                exc += ee
                P |= pp
            grp = groups.get(key)
            if grp is None:
                groups[key] = [exc, {P: 1}]
            else:
                grp[1][P] = grp[1].get(P, 0) + 1
            n_fill += 1
        full = (1 << delta) - 1
        best = None
        for key, (exc, Ps) in groups.items():
            tot = sum(Ps.values())
            n_pairs += tot * tot
            for PI in Ps:
                for PJ in Ps:
                    n_classes += 1
                    if exc < cost[PI | PJ]:
                        fail(tag, "a pair of excess %d clears levels of total cost %d" % (exc, cost[PI | PJ]))
                    if PI | PJ == full and (best is None or exc < best):
                        best = exc
        if best == sum(m):
            n_attained += 1
    print("  %d configurations (interior, bottom, top; delta = 1..5), %d fillings, %d pairs with equal level sizes:"
          % (len(configs), n_fill, n_pairs))
    print("  j >= sum_{i in P} m_i holds for every pair; the bound j >= w for pairs clearing every block level"
          " is attained in %d configurations" % n_attained)


# =============================================================== check (d)
def check_d():
    print("(d) Theorem B: fillings with prescribed level sizes (brute force)")
    rng = random.Random(4)
    configs = [
        (6, 7, [1, 1, 3, 3, 0, 0], (1, 3)), (6, 7, [1, 1, 3, 3, 2, 0], (1, 3)), (6, 7, [3, 1, 3, 1, 2, 4], (1, 3)),
        (6, 7, [1, 1, 3, 3, 5, 5], (1, 3, 5)), (6, 7, [2, 2, 4, 4, 6, 0], (2, 4, 6)), (6, 7, [2, 2, 1, 0, 3, 0], (2,)),
        (6, 7, [1, 4, 4, 1, 2, 0], (1, 4)), (6, 7, [2, 2, 5, 5, 3, 0], (2, 5)), (6, 7, [1, 1, 2, 2, 0, 0], (1, 3)),
        (6, 7, [3, 3, 3, 0, 0, 0], (1, 3)), (6, 7, [2, 2, 2, 4, 4, 4], (2, 4)), (6, 7, [1, 3, 5, 0, 0, 0], (1, 3, 5)),
        (3, 7, [1, 3, 3], (1, 3)), (3, 7, [2, 2, 4], (2, 4)), (3, 7, [1, 1, 1], (1,)), (3, 7, [6, 6, 4], (4, 6)),
        (5, 11, [2, 2, 4, 4, 0], (2, 4)), (5, 11, [1, 1, 3, 3, 6], (1, 3, 6)), (5, 11, [2, 2, 5, 5, 3], (2, 5)),
        (5, 11, [3, 3, 3, 5, 5], (3, 5)), (10, 11, "B", (1, 3)), (10, 11, "B", (1, 3)), (10, 11, "B", (1, 3)),
        (10, 11, "B", (1, 3)), (4, 5, [1, 1, 3, 3], (1, 3)), (4, 5, [2, 2, 4, 4], (2, 4)), (4, 5, [1, 3, 1, 3], (1, 3)),
    ]
    n_cfg = n_fill = n_pairs = n_nonzero = 0
    for (M, q, sizes, levels) in configs:
        ch = rng.choice(charts(M, q))
        if sizes == "B":
            # the situation of Theorem C(iii): a singular S at the levels 1 and 2
            S = rng.choice([S for S in ch.singular_sets() if len(S) <= 4])
            sizes = [0] * M
            for a in S:
                sizes[a] = 3
            for a in rng.sample([a for a in range(M) if a not in S], 2):
                sizes[a] = 1
        fib = random_fibres(M, q, sizes, rng)
        T = level_sets(fib, q)
        tag = "d (%d,%d) eta=%d K=%s levels=%s" % (M, q, ch.eta, points_of(fib), levels)
        assert all(levels[i + 1] - levels[i] >= 2 for i in range(len(levels) - 1))
        assert all(1 <= c <= q - 1 for c in levels)
        p = [len(T[m]) + sum(1 for c in levels if c == m) - sum(1 for c in levels if c - 1 == m)
             for m in range(q)]
        fl = fillings_with_sizes(fib, q, p, fib, q)
        # the claimed description: one single raise in a fibre of size c_i for every i
        support = [a for a, f in enumerate(fib) if f]
        shells = [[a for a in T[c - 1] if a not in T[c]] for c in levels]
        expected = set()
        for choice in itertools.product(*shells):
            filling = []
            for a in support:
                c = len(fib[a])
                if a in choice:
                    filling.append(tuple(range(c - 1)) + (c,))
                else:
                    filling.append(tuple(range(c)))
            expected.add(tuple(filling))
        if set(f[4] for f in fl) != expected:
            fail(tag, "fillings with level sizes p are not exactly the products of single raises "
                      "(%d found, %d expected)" % (len(fl), len(expected)))
        # total contribution to lambda_r
        fact = [1] * q
        for m in range(1, q):
            fact[m] = fact[m - 1] * m % q
        coef = 1
        for m in range(q):
            ex = len(T[m]) - p[m]
            coef = coef * (pow(fact[m], ex, q) if ex >= 0 else pow(pow(fact[m], q - 2, q), -ex, q)) % q
        tot = 0
        for (levI, wI, eI, sI, _) in fl:
            for (levJ, wJ, eJ, sJ, _) in fl:
                t = eI * eJ * sI * sJ
                for m in range(q):
                    if levI[m] or levJ[m]:
                        t = t * ch.minor(levI[m], levJ[m]) % q
                tot += t
                n_pairs += 1
        tot = tot * coef % q
        rhs = 1
        for c in levels:
            rhs = rhs * pow(c, q - 2, q) % q
        for m in range(q):
            if all(m != c and m != c - 1 for c in levels):
                rhs = rhs * ch.d(T[m]) % q
        for c in levels:
            rhs = rhs * Q_level(c, fib, T, ch) % q
        if tot != rhs:
            fail(tag, "contribution %d but the product formula gives %d" % (tot, rhs))
        if fl and any(f[1] != e_of(fib) + len(levels) for f in fl):
            fail(tag, "a filling with level sizes p does not have excess r")
        n_cfg += 1
        n_fill += len(fl)
        n_nonzero += 1 if rhs else 0
    print("  %d configurations (r = 1, 2, 3 raises); %d fillings found, all products of single raises;"
          " %d pairs summed; product formula agrees in all, nonzero in %d" % (n_cfg, n_fill, n_pairs, n_nonzero))


# =============================================================== check (e)
def is_prime(n):
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for p in small:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in small:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def prime_factors(n):
    out, p = [], 2
    while p * p <= n:
        if n % p == 0:
            out.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        out.append(n)
    return out


class ModL:
    """F_N = (omega^(kl)) modulo a prime l = 1 (mod N), N = qM, with the splitting of the Setting."""

    def __init__(self, M, q, start=3 * 10 ** 15):
        N = M * q
        l = start - start % N + 1
        while not is_prime(l):
            l += N
        self.M, self.q, self.N, self.l = M, q, N, l
        g = 2
        while True:
            w = pow(g, (l - 1) // N, l)
            if all(pow(w, N // p, l) != 1 for p in prime_factors(N)):
                break
            g += 1
        self.omega = w
        nu = q * pow(q, -1, M)           # = 1 mod M, = 0 mod q
        nu2 = M * pow(M, -1, q)          # = 0 mod M, = 1 mod q
        self.eta, self.zeta = pow(w, nu, l), pow(w, nu2, l)
        self.ep = [pow(self.eta, k, l) for k in range(M)]
        self.zp = [pow(self.zeta, k, l) for k in range(q)]
        self.crt = lambda a, b: (a * nu + b * nu2) % N

    def D(self, K):
        ep, zp, M, q = self.ep, self.zp, self.M, self.q
        return det_mod([[ep[a * a2 % M] * zp[b * b2 % q] for (a2, b2) in K] for (a, b) in K], self.l)

    def D_omega(self, K, conj=False):
        """The same minor as a minor of F_N in the index k <-> (a, b); conj: omega -> omega^(-1)."""
        N, l = self.N, self.l
        w = pow(self.omega, -1, l) if conj else self.omega
        wp = [pow(w, k, l) for k in range(N)]
        ks = [self.crt(a, b) for (a, b) in K]
        return det_mod([[wp[k * k2 % N] for k2 in ks] for k in ks], l)

    def d(self, T):
        ep, M = self.ep, self.M
        return det_mod([[ep[a * b % M] for b in T] for a in T], self.l)

    def Delta(self, B):
        zp, q = self.zp, self.q
        return det_mod([[zp[b * b2 % q] for b2 in B] for b in B], self.l)


def matmul_mod(A, B, m):
    Bt = list(zip(*B))
    return [[sum(x * y for x, y in zip(row, col)) % m for col in Bt] for row in A]


def prop_E_case(ctx, S, B, K, tag, deep):
    """Proposition E for S x B contained in K.  Returns False if d_S or Delta_B is 0 mod l."""
    M, q, l, ep, zp = ctx.M, ctx.q, ctx.l, ctx.ep, ctx.zp
    dS, DB = ctx.d(S), ctx.Delta(B)
    if dS == 0 or DB == 0:
        return False
    FSi = inv_mat_mod([[ep[a * b % M] for b in S] for a in S], l)
    ZBi = inv_mat_mod([[zp[b * b2 % q] for b2 in B] for b in B], l)
    FaS = [[ep[a * b % M] for b in S] for a in range(M)]
    FSa = [[ep[a * b % M] for b in range(M)] for a in S]
    PS = matmul_mod(matmul_mod(FaS, FSi, l), FSa, l)
    ZbB = [[zp[b * b2 % q] for b2 in B] for b in range(q)]
    ZBb = [[zp[b * b2 % q] for b2 in range(q)] for b in B]
    PB = matmul_mod(matmul_mod(ZbB, ZBi, l), ZBb, l)
    Sset, Bset = set(S), set(B)
    rest = [(a, b) for (a, b) in K if not (a in Sset and b in Bset)]
    mat = []
    for (a, b) in rest:
        row = []
        for (a2, b2) in rest:
            h = ep[a * a2 % M]
            RB = (zp[b * b2 % q] - PB[b][b2]) % l
            RS = (h - PS[a][a2]) % l
            row.append((h * RB + RS * PB[b][b2]) % l)
        mat.append(row)
    rhs = pow(dS, len(B), l) * pow(DB, len(S), l) * det_mod(mat, l) % l
    if ctx.D(K) != rhs:
        fail(tag, "Proposition E fails for S=%s B=%s K=%s" % (S, B, K))
    if deep:
        # the remark after the proof: R_S vanishes on S, and is a bordered minor over d_S otherwise
        inv_dS = pow(dS, -1, l)
        for a in range(M):
            for a2 in range(M):
                RS = (ep[a * a2 % M] - PS[a][a2]) % l
                if a in Sset or a2 in Sset:
                    want = 0
                else:
                    rows, cols = list(S) + [a], list(S) + [a2]
                    want = det_mod([[ep[x * y % M] for y in cols] for x in rows], l) * inv_dS % l
                if RS != want:
                    fail(tag, "R_S(%d,%d) is not the bordered minor / d_S" % (a, a2))
    return True


def theorem_D_identity(ctx, fib, xs, tag, DK=None):
    """Checks the hypotheses of Theorem D as stated (returns None if they fail) and then the identity."""
    M, q, l = ctx.M, ctx.q, ctx.l
    r = len(xs)
    c = [len(f) for f in fib]
    Tl = [tuple(a for a in range(M) if c[a] > m) for m in range(q + 1)]
    for j, x in enumerate(xs):
        if tuple(a for a in range(M) if x in fib[a]) != Tl[j]:
            return None
    dT = [ctx.d(Tl[j]) for j in range(r + 1)]
    if any(x == 0 for x in dT):
        return None
    X = set(xs)
    Tr = Tl[r]
    fib2 = [sorted((set(f) - X) | (X if a in Tr else set())) for a, f in enumerate(fib)]
    # level sets of K-tilde: T_r at the levels m <= r, T_m above
    c2 = [len(f) for f in fib2]
    for m in range(q):
        want = Tr if m <= r else Tl[m]
        if tuple(a for a in range(M) if c2[a] > m) != want:
            fail(tag, "level sets of K-tilde are not as stated (K=%s, x=%s)" % (points_of(fib), xs))
            break
    Dl = [1] + [ctx.Delta(xs[:j]) for j in range(1, r + 1)]
    if any(x == 0 for x in Dl):
        fail(tag, "a leading minor Delta_j vanishes modulo l")
        return None
    lhs = ctx.D(points_of(fib)) if DK is None else DK
    rhs = ctx.D(points_of(fib2))
    for j in range(r):
        ratio = Dl[j + 1] * pow(Dl[j], -1, l) % l
        lhs = lhs * pow(ratio, len(Tr), l) * dT[r] % l
        rhs = rhs * pow(ratio, len(Tl[j]), l) * dT[j] % l
    if lhs != rhs:
        fail(tag, "Theorem D identity fails (K=%s, x=%s)" % (points_of(fib), xs))
    return True


def peelable_fibres(M, q, r, rng, Tr=None, lower=None):
    """Fibres satisfying the hypothesis of Theorem D with the points x_0..x_(r-1)."""
    xs = rng.sample(range(q), r)
    others = [b for b in range(q) if b not in xs]
    if Tr is None:
        Tr = [a for a in range(M) if rng.random() < 0.4]
    fib = []
    for a in range(M):
        if a in Tr:
            extra = rng.sample(others, rng.randint(1, len(others)))
            fib.append(sorted(xs + extra))
        else:
            c = rng.randint(0, r) if lower is None else lower.get(a, 0)
            fib.append(sorted(xs[:c]))
    return fib, xs


def check_e():
    print("(e) Propositions E, F, Theorem D and the dual form")
    rng = random.Random(105110)
    ctxs = {}
    for (M, q) in [(15, 7), (10, 11), (21, 5), (35, 3), (6, 7), (3, 5), (5, 3), (4, 3)]:
        ctxs[(M, q)] = ModL(M, q)
    print("  primes l = 1 (mod N): " + ", ".join("N=%d: %d" % (N, l) for N, l in
                                                    sorted(set((c.N, c.l) for c in ctxs.values()))))
    # the identification omega^(kl) = eta^(aa') zeta^(bb')
    n_id = 0
    for key in [(15, 7), (10, 11), (6, 7)]:
        ctx = ctxs[key]
        pts = [(a, b) for a in range(ctx.M) for b in range(ctx.q)]
        for _ in range(4):
            K = rng.sample(pts, rng.randint(1, 30))
            if ctx.D(K) != ctx.D_omega(K):
                fail("e setting", "D_K as a minor of F_N differs from det(eta^(aa') zeta^(bb')), K=%s" % K)
            n_id += 1
    print("  setting: D_K = minor of F_N under k <-> (a,b): %d sets" % n_id)

    # ---- Proposition E
    n_E = n_E_skip = 0
    for key, reps in [((15, 7), 40), ((10, 11), 40), ((21, 5), 20), ((35, 3), 10)]:
        ctx = ctxs[key]
        M, q = key
        pts = [(a, b) for a in range(M) for b in range(q)]
        for i in range(reps):
            S = sorted(rng.sample(range(M), rng.randint(1, min(M - 1, 8))))
            B = sorted(rng.sample(range(q), rng.randint(1, q - 1)))
            blk = set((a, b) for a in S for b in B)
            rest = [p for p in pts if p not in blk]
            extra = rng.sample(rest, rng.randint(0, min(len(rest), 14)))
            K = sorted(blk | set(extra))
            if prop_E_case(ctx, S, B, K, "e PropE N=%d" % ctx.N, deep=(i < 3)):
                n_E += 1
            else:
                n_E_skip += 1
    print("  Proposition E modulo l: %d sets agree (%d skipped: d_S or Delta_B = 0 mod l)" % (n_E, n_E_skip))

    # ---- Proposition F
    n_F = 0
    for key, reps in [((15, 7), 30), ((10, 11), 30), ((21, 5), 15), ((35, 3), 10)]:
        ctx = ctxs[key]
        M, q, l = ctx.M, ctx.q, ctx.l
        for _ in range(reps):
            perm = rng.sample(range(q), q)
            sizes = [rng.randint(0, q) if rng.random() < 0.6 else 0 for _ in range(M)]
            fib = [sorted(perm[:c]) for c in sizes]
            K = points_of(fib)
            if not K:
                continue
            Dl = [1] + [ctx.Delta(perm[:m]) for m in range(1, q + 1)]
            rhs = 1
            for m in range(q):
                Tm = [a for a in range(M) if sizes[a] > m]
                if Tm:
                    rhs = rhs * pow(Dl[m + 1] * pow(Dl[m], -1, l) % l, len(Tm), l) * ctx.d(Tm) % l
            if ctx.D(K) != rhs:
                fail("e PropF N=%d" % ctx.N, "Proposition F fails for K=%s" % K)
            n_F += 1
    print("  Proposition F modulo l: %d Ferrers sets agree" % n_F)

    # ---- Theorem D, random sets satisfying the hypothesis
    n_D = n_D_skip = 0
    for key, reps in [((15, 7), 30), ((10, 11), 30), ((21, 5), 15), ((35, 3), 8)]:
        ctx = ctxs[key]
        M, q = key
        for _ in range(reps):
            r = rng.randint(1, q - 1)
            fib, xs = peelable_fibres(M, q, r, rng)
            res = theorem_D_identity(ctx, fib, xs, "e ThmD N=%d" % ctx.N)
            if res is None:
                n_D_skip += 1
            else:
                n_D += 1
    print("  Theorem D modulo l, random sets with r nested levels: %d agree (%d skipped: hypothesis not met)"
          % (n_D, n_D_skip))

    # ---- Theorem D, exhaustively: every K and every admissible x_0..x_(r-1), small N
    for key in [(3, 5), (5, 3), (4, 3)]:
        ctx = ctxs[key]
        M, q = key
        pts = [(a, b) for a in range(M) for b in range(q)]
        n_K = n_inst = 0
        for mask in range(1, 1 << (M * q)):
            K = [pts[i] for i in range(M * q) if mask >> i & 1]
            fib = fibres_of(K, M)
            c = [len(f) for f in fib]
            slices = {x: tuple(a for a in range(M) if x in fib[a]) for x in range(q)}
            T0 = tuple(a for a in range(M) if c[a] > 0)
            if T0 not in slices.values():
                continue
            DK = ctx.D(K)
            found = False
            stack = [[]]
            while stack:
                xs = stack.pop()
                j = len(xs)
                Tj = tuple(a for a in range(M) if c[a] > j)
                for x in range(q):
                    if x not in xs and slices[x] == Tj:
                        xs2 = xs + [x]
                        res = theorem_D_identity(ctx, fib, xs2, "e ThmD exhaustive N=%d" % ctx.N, DK)
                        if res:
                            n_inst += 1
                            found = True
                        stack.append(xs2)
            n_K += 1 if found else 0
        print("  Theorem D exhaustively over all subsets of Z_%d x Z_%d: %d sets admit a peeling, "
              "%d (set, x-sequence) instances agree" % (M, q, n_K, n_inst))

    # ---- dual form: D_K = det F_N * N^(-|K^c|) * conj(D_(K^c))
    n_dual = 0
    for key, reps in [((6, 7), 12), ((3, 5), 12), ((15, 7), 3)]:
        ctx = ctxs[key]
        M, q, N, l = ctx.M, ctx.q, ctx.N, ctx.l
        pts = [(a, b) for a in range(M) for b in range(q)]
        detF = ctx.D_omega(pts)
        for _ in range(reps):
            K = sorted(rng.sample(pts, rng.randint(1, N - 1)))
            Kc = [p for p in pts if p not in set(K)]
            rhs = detF * pow(pow(N, len(Kc), l), -1, l) * ctx.D_omega(Kc, conj=True) % l
            if ctx.D_omega(K) != rhs:
                fail("e dual N=%d" % N, "Jacobi identity fails for K=%s" % K)
            # slices and level sets of the complement
            fib, fibc = fibres_of(K, M), fibres_of(Kc, M)
            T, Tc = level_sets(fib, q), level_sets(fibc, q)
            for m in range(q):
                if set(Tc[m]) != set(range(M)) - set(T[q - 1 - m]):
                    fail("e dual N=%d" % N, "level sets of the complement, K=%s" % K)
            n_dual += 1
    print("  dual form (Jacobi) modulo l, and level sets of the complement: %d sets agree" % n_dual)

    # ---- valuation statements, q-adically, base 10 above 11
    chs = charts(10, 11)
    q, M = 11, 10
    r_d = 12

    def v_d(ch, T):
        """v_q(d_T) for the Teichmueller lift of eta (exact if < r_d)."""
        if not T:
            return 0
        Q = q ** r_d
        eta_Q = pow(ch.eta, q ** (r_d - 1), Q)
        d = det_int([[pow(eta_Q, a * b, Q) for b in T] for a in T]) % Q
        assert d != 0
        return vq_int(d, q)

    n_vF = n_vD = n_vF_sing = n_vD_sing = 0
    for i in range(24):
        ch = chs[i % len(chs)]
        sing = ch.singular_sets()
        S = rng.choice([S for S in sing if len(S) <= 4])
        # Ferrers set with S among its level sets (multiplicity 1..3)
        perm = rng.sample(range(q), q)
        mult = rng.randint(1, 3)
        s = rng.randint(0, 2)
        sizes = [0] * M
        for a in S:
            sizes[a] = s + mult + (1 if rng.random() < 0.3 else 0)
        if s:
            for a in rng.sample([a for a in range(M) if a not in S], rng.randint(1, 3)):
                sizes[a] = rng.randint(1, s)
        fib = [sorted(perm[:c]) for c in sizes]
        e = e_of(fib)
        T = level_sets(fib, q)
        pred = e + (q - 1) * sum(v_d(ch, T[m]) for m in range(q))
        v, exact = vDK(points_of(fib), ch, pred + q)
        if not exact or v != pred:
            fail("e PropF valuation", "eta=%d K=%s: v(D_K)=%s, formula %d" % (ch.eta, points_of(fib), v, pred))
        n_vF += 1
        n_vF_sing += 1 if pred > e else 0
    for i in range(24):
        ch = chs[i % len(chs)]
        sing = set(ch.singular_sets())
        for _ in range(1000):
            # nested up to level r, with a singular set among T_0..T_(r-1), arbitrary above
            S = rng.choice([S for S in sing if len(S) <= 4])
            r = rng.randint(1, 3)
            Tr = [a for a in S if rng.random() < 0.6]
            lower = {a: r for a in S if a not in Tr}       # S = T_(r-1) at least
            for a in rng.sample([a for a in range(M) if a not in S], rng.randint(0, 2)):
                lower[a] = rng.randint(1, r)
            xs = rng.sample(range(q), r)
            others = [b for b in range(q) if b not in xs]
            fib = []
            for a in range(M):
                if a in Tr:
                    fib.append(sorted(xs + rng.sample(others, rng.randint(1, 3))))
                else:
                    fib.append(sorted(xs[:lower.get(a, 0)]))
            T = level_sets(fib, q)
            if all(T[m] not in sing for m in range(r, q)) and any(T[j] in sing for j in range(r)):
                break
        else:
            raise RuntimeError("no peelable example found")
        ctx = ctxs[(10, 11)]
        if theorem_D_identity(ctx, fib, xs, "e ThmD valuation (identity mod l)") is None:
            fail("e ThmD valuation", "constructed set does not satisfy the hypothesis: %s" % points_of(fib))
        e = e_of(fib)
        pred = e + (q - 1) * sum(v_d(ch, T[j]) for j in range(r))
        v, exact = vDK(points_of(fib), ch, pred + q)
        if not exact or v != pred:
            fail("e ThmD valuation", "eta=%d K=%s x=%s: v(D_K)=%s, formula %d"
                 % (ch.eta, points_of(fib), xs, v, pred))
        n_vD += 1
        n_vD_sing += 1 if pred > e else 0
    print("  valuation of Proposition F, q-adically at (10,11): %d Ferrers sets agree (%d with a singular level, v > e)"
          % (n_vF, n_vF_sing))
    print("  valuation of Theorem D, q-adically at (10,11): %d peelable sets agree (%d with v > e)"
          % (n_vD, n_vD_sing))


# ====================================================================== main
def main():
    t0 = time.time()
    print("exact_check.py: independent check of the higher-order lifting report")
    print("python %s" % sys.version.split()[0])
    for fn in (check_a, check_b1, check_b2, check_c0, check_c, check_d, check_e):
        t = time.time()
        fn()
        print("  [%.1f s]" % (time.time() - t))
        sys.stdout.flush()
    print("total time %.1f s" % (time.time() - t0))
    if FAILS:
        print("%d FAILURES" % len(FAILS))
        for tag, info in FAILS[:50]:
            print("  [%s] %s" % (tag, info))
        sys.exit(1)
    print("all checks passed")


if __name__ == "__main__":
    main()
