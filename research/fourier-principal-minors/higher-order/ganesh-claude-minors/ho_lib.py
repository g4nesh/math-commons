#!/usr/bin/env python3
"""Higher-order lifting: shared code.

Notation as in README.md and in the second-order report.  M is a positive
integer, q a prime not dividing M, a chart is an irreducible factor h of Phi_M
mod q with residue field FF = F_q[x]/(h) and eta = x.  An index set K in
Z_M x Z_q is given by its fibres: fib[a] = sorted list of the b with (a, b) in K.

The window polynomial of K in a chart is
    P(tau) = det( eta^(a a') E(tau b b') ),   E(z) = sum_{m=0}^{q-1} z^m / m!,
a polynomial over FF.  By Theorem A of README.md
    P(tau) = prod_a V(K_a)^2 * prod_m (m!)^(-|T_m|) * sum_j Lambda_j tau^(e + j),
and the principal minor D_K is nonzero as soon as one Lambda_j with
0 <= j <= q - 2 is nonzero in FF.

lam_coeffs(...) returns Lambda_0 .. Lambda_(L-1).  It is computed from the
divided-difference matrix H of window.c, here in pure Python as a reference;
window_c(...) runs the C program on many sets.

Field elements are in the log representation of so_lib.Field (zero is -1).
Standard library only.
"""
import os
import subprocess
import sys
from math import comb, factorial

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "second-order", "ganesh-claude-minors"))
import so_lib as so  # noqa: E402


def code_table(K):
    """Map the log representation of so_lib.Field to the integer sum_i c_i q^i of
    the element sum_i c_i x^i, the format printed by window.c, and back."""
    p, f = K.p, K.f
    to_code, from_code = {-1: 0}, {0: -1}
    xp = [0]
    for _ in range(f - 1):
        xp.append(K.mul(xp[-1], K.eta))
    for code in range(1, p ** f):
        acc, c = -1, code
        for i in range(f):
            d = c % p
            c //= p
            if d:
                acc = K.add(acc, K.mul(K.from_int(d), xp[i]))
        to_code[acc] = code
        from_code[code] = acc
    return to_code, from_code


def hom_sym(nodes, n, q):
    """h_0 .. h_(n-1) of the given nodes, modulo q."""
    hs = [1] + [0] * (n - 1)
    for x in nodes:
        for i in range(1, n):
            hs[i] = (hs[i] + x * hs[i - 1]) % q
    return hs


def build_H(ch, fib, L):
    """The matrix H of window.c: rows (a, k), columns (a', b'), entries lists of L field elements."""
    K, q, M, ep = ch.K, ch.q, ch.M, ch.ep
    invfact = [pow(factorial(n), q - 2, q) if q > 2 else 1 for n in range(q)]
    rows = [(a, k) for a in sorted(fib) for k in range(len(fib[a]))]
    cols = [(a, b) for a in sorted(fib) for b in fib[a]]
    H = []
    for a, k in rows:
        hs = hom_sym(fib[a][:k + 1], L, q)
        line = []
        for a2, b2 in cols:
            root = ep[(a * a2) % M]
            entry = []
            for i in range(L):
                v = pow(b2, k + i, q) * hs[i] * invfact[k + i] % q if k + i <= q - 1 else 0
                entry.append(K.mul(K.from_int(v), root) if v else -1)
            line.append(entry)
        H.append(line)
    return H


def series_det(K, A, L):
    """Determinant modulo tau^L of a square matrix of power series over the field K
    (lists of L log-represented coefficients), by elimination with pivots of least
    valuation.  Returns the list of L coefficients."""
    n = len(A)
    A = [[e[:] for e in row] for row in A]

    def val(s, cur):
        for i in range(cur):
            if s[i] >= 0:
                return i
        return cur

    cur, V, neg = L, 0, False
    U = [0] + [-1] * (L - 1)
    for c in range(n):
        best, bv = None, cur
        for r in range(c, n):
            v = val(A[r][c], cur)
            if v < bv:
                best, bv = r, v
        if best is None:
            return [-1] * L
        if best != c:
            A[c], A[best] = A[best], A[c]
            neg = not neg
        nc = cur - bv
        u = A[c][c][bv:bv + nc]
        ui = [K.inv(u[0])]
        for i in range(1, nc):
            s = -1
            for j in range(1, i + 1):
                s = K.add(s, K.mul(u[j], ui[i - j]))
            ui.append(K.neg(K.mul(s, ui[0])))
        newU = []
        for i in range(nc):
            s = -1
            for j in range(i + 1):
                s = K.add(s, K.mul(U[j], u[i - j]))
            newU.append(s)
        U = newU + [-1] * (L - nc)
        V += bv
        for r in range(c + 1, n):
            x = A[r][c]
            if val(x, cur) >= cur:
                continue
            mu = []
            for i in range(nc):
                s = -1
                for j in range(i + 1):
                    s = K.add(s, K.mul(x[j + bv], ui[i - j]))
                mu.append(s)
            if all(m < 0 for m in mu):
                continue
            Ar, Ac = A[r], A[c]
            for j in range(c + 1, n):
                p, y = Ac[j], Ar[j]
                for i in range(nc):
                    s = -1
                    for t in range(i + 1):
                        if mu[t] >= 0 and p[i - t] >= 0:
                            s = K.add(s, K.mul(mu[t], p[i - t]))
                    if s >= 0:
                        y[i] = K.sub(y[i], s)
        cur = nc
    out = [-1] * L
    for j in range(V, L):
        if j - V < cur:
            out[j] = K.neg(U[j - V]) if neg else U[j - V]
    return out


def normaliser(ch, fib):
    """prod_a V(K_a) * prod_m (m!)^(-|T_m|) as a field element, V with increasing order."""
    K, q = ch.K, ch.q
    T = so.level_sets(fib, q)
    v = 1
    for b in fib.values():
        for i in range(len(b)):
            for j in range(i + 1, len(b)):
                v = v * (b[j] - b[i]) % q
    for m in range(q):
        v = v * pow(pow(factorial(m), q - 2, q), len(T[m]), q) % q
    return K.from_int(v)


def lam_coeffs(ch, fib, L, detH=None):
    """Lambda_0 .. Lambda_(L-1) of README.md, Theorem A, in the chart ch."""
    K = ch.K
    if detH is None:
        detH = series_det(K, build_H(ch, fib, L), L)
    inv = K.inv(normaliser(ch, fib))
    return [K.mul(x, inv) for x in detH]


def window_c(binary, ch, fibs, L):
    """det H for many sets through the C program; returns lists of L field elements."""
    K = ch.K
    if not hasattr(K, "_codes"):
        K._codes = code_table(K)
    from_code = K._codes[1]
    q, M = ch.q, ch.M
    text = "\n".join(" ".join(str(sum(1 << b for b in fib.get(a, []))) for a in range(M)) for fib in fibs) + "\n"
    args = [binary, str(M), str(q), str(L), str(K.f)] + [str(c) for c in ch.h[:-1]]
    out = subprocess.run(args, input=text, capture_output=True, text=True, check=True).stdout.split("\n")
    return [[from_code[int(x)] for x in line.split()] for line in out if line.strip()]


def q_form(ch, T_low, T_high, shell, sig):
    """Q = sum over a, a' in shell of sig_a sig_a' Cof_(a,a') E_(a,a'), where Cof is
    the cofactor of F_M[T_low] and E the bordered minor of F_M[T_high]; shell is
    T_low minus T_high.  This is the inner sum of Lambda in the second-order report."""
    K = ch.K
    acc = -1
    for a in shell:
        if sig[a] < 0:
            continue
        for a2 in shell:
            if sig[a2] < 0:
                continue
            term = K.mul(K.mul(sig[a], sig[a2]), K.mul(ch.cofactor(T_low, a, a2), ch.bordered(T_high, a, a2)))
            acc = K.add(acc, term)
    return acc


def separated_raises(ch, fib, shells):
    """Theorem B: kappa_p * Lambda^(p) for the profile of one single raise in each
    of the given shells c_1 < c_2 < ... with c_(i+1) - c_i >= 2."""
    K, q = ch.K, ch.q
    assert all(y - x >= 2 for x, y in zip(shells, shells[1:]))
    T = so.level_sets(fib, q)
    sig = {a: K.from_int(sum(v)) for a, v in fib.items()}
    used = {c - 1 for c in shells} | set(shells)
    acc = 0
    for m in range(q):
        if m not in used:
            acc = K.mul(acc, ch.minor(T[m], T[m]))
    for c in shells:
        shell = [a for a in T[c - 1] if a not in T[c]]
        acc = K.mul(acc, K.mul(K.from_int(pow(c, q - 2, q)), q_form(ch, T[c - 1], T[c], shell, sig)))
    return acc


def excess(fib):
    return sum(comb(len(v), 2) for v in fib.values())


# ---- exact arithmetic modulo a large prime l = 1 (mod N), for identities over Z[zeta_N]

def is_prime(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def prime_one_mod(N, start=10 ** 9):
    """Least prime l >= start with l = 1 (mod N)."""
    l = start - start % N + 1
    while not is_prime(l):
        l += N
    return l


def root_mod(N, l):
    """An element of order exactly N in F_l."""
    assert (l - 1) % N == 0
    primes = [r for r in range(2, N + 1) if N % r == 0 and is_prime(r)]
    for g in range(2, 10000):
        w = pow(g, (l - 1) // N, l)
        if all(pow(w, N // r, l) != 1 for r in primes):
            return w
    raise ValueError


def det_mod(A, l):
    n = len(A)
    A = [r[:] for r in A]
    d = 1
    for c in range(n):
        piv = next((r for r in range(c, n) if A[r][c] % l), None)
        if piv is None:
            return 0
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            d = -d
        d = d * A[c][c] % l
        inv = pow(A[c][c], l - 2, l)
        for r in range(c + 1, n):
            f = A[r][c] * inv % l
            if f:
                Ar, Ac = A[r], A[c]
                for j in range(c, n):
                    Ar[j] = (Ar[j] - f * Ac[j]) % l
    return d % l


class ModN:
    """F_N = (omega^(kl)) over F_l with omega of order N = M * Mp, gcd(M, Mp) = 1, and
    the splitting omega^(kl) = eta^(a a') zeta^(b b') under Z_N = Z_M x Z_Mp."""

    def __init__(self, M, Mp, l=None):
        self.M, self.Mp, self.N = M, Mp, M * Mp
        self.l = l or prime_one_mod(self.N)
        w = root_mod(self.N, self.l)
        nu = next(x for x in range(self.N) if x % M == 1 % M and x % Mp == 0)
        nup = next(x for x in range(self.N) if x % M == 0 and x % Mp == 1 % Mp)
        self.eta, self.zeta = pow(w, nu, self.l), pow(w, nup, self.l)
        self.ep = [pow(self.eta, d, self.l) for d in range(M)]
        self.zp = [pow(self.zeta, d, self.l) for d in range(Mp)]

    def minor(self, K):
        """D_K for a list K of pairs (a, b)."""
        l, M, Mp = self.l, self.M, self.Mp
        return det_mod([[self.ep[a * a2 % M] * self.zp[b * b2 % Mp] % l for (a2, b2) in K] for (a, b) in K], l)

    def d(self, T):
        """det F_M[T, T]."""
        return det_mod([[self.ep[a * a2 % self.M] for a2 in T] for a in T], self.l) if T else 1

    def delta(self, X):
        """det (zeta^(b b'))_(b, b' in X)."""
        return det_mod([[self.zp[b * b2 % self.Mp] for b2 in X] for b in X], self.l) if X else 1
