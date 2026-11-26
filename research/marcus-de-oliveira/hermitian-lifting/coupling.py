"""Exact check of Lemma C (coupling of interlacing sequences) and Lemma I (interlacing lemma) of README.md.
Standard library only.

Lemma C.  Let sigma_1 >= ... >= sigma_n and sigma_i >= s_i >= sigma_{i+1} (i < n).  The recursive rule

    slot 1 takes index 1 with probability (s_1 - sigma_2)/(sigma_1 - sigma_2) and index 2 otherwise
    (index 1 with probability 1 if sigma_1 = sigma_2); then the same rule is applied to slots 2..n-1 and the remaining pool

defines a probability distribution on injections tau: {1..n-1} -> {1..n} with
    E prod_{i in T} sigma_{tau(i)} = prod_{i in T} s_i     for every subset T of {1..n-1}.

Lemma I.  Let a_1 <= ... <= a_n, a_i <= alpha_i <= a_{i+1}, beta real with a_1 + beta > 0, c_1..c_{n-1} complex.  Then
    prod_i (a_i + beta) * prod_{i<n} (alpha_i + c_i)/(alpha_i + beta) = sum_tau w_tau (a_{j(tau)} + beta) prod_{i<n} (a_{tau(i)} + c_i),
where j(tau) is the index outside the image of tau and w is the distribution of Lemma C for sigma_i = 1/(a_i + beta),
s_i = 1/(alpha_i + beta).  Each term on the right is one of the n! products of the conjecture for the spectra a and (c, beta).

The script builds the distribution with rational arithmetic and checks both statements on random rational data, including
ties and boundary cases, and Gaussian-rational c."""
import itertools, random, sys
from fractions import Fraction as Fr
def coupling(sigma, s):
    """list of (weight, injection as a tuple of pool indices); sigma decreasing, s interlacing"""
    def rec(pool, slots):                      # pool: list of indices into sigma (sorted by decreasing value); slots: list of s values
        if not slots: return [(Fr(1), ())]
        i1, i2 = pool[0], pool[1]; s1 = slots[0]
        assert sigma[i1] >= s1 >= sigma[i2], "not interlacing"
        out = []
        p1 = Fr(1) if sigma[i1] == sigma[i2] else (s1 - sigma[i2]) / (sigma[i1] - sigma[i2])
        if p1: out += [(p1 * w, (i1,) + t) for w, t in rec(pool[1:], slots[1:])]
        if p1 != 1: out += [((1 - p1) * w, (i2,) + t) for w, t in rec([pool[0]] + pool[2:], slots[1:])]
        return out
    return rec(list(range(len(sigma))), list(s))
class GQ:                                                     # Gaussian rationals
    def __init__(self, re, im=0): self.re = Fr(re); self.im = Fr(im)
    def __add__(self, o): o = o if isinstance(o, GQ) else GQ(o); return GQ(self.re + o.re, self.im + o.im)
    __radd__ = __add__
    def __mul__(self, o): o = o if isinstance(o, GQ) else GQ(o); return GQ(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
    __rmul__ = __mul__
    def __eq__(self, o): return self.re == o.re and self.im == o.im
rnd = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 0); NT = int(sys.argv[2]) if len(sys.argv) > 2 else 300
def rq(lo=-9, hi=9): return Fr(rnd.randint(lo * 4, hi * 4), rnd.randint(1, 6))
okC = okI = True; cases = 0; support = 0; ties = 0
for trial in range(NT):
    n = rnd.randint(2, 7)
    a = sorted(rq() for _ in range(n))
    if trial % 4 == 0: a[rnd.randrange(n - 1) + 1] = a[rnd.randrange(n - 1)]; a.sort()            # ties
    alpha = []
    for i in range(n - 1):
        mode = rnd.randrange(4)
        alpha.append(a[i] if mode == 0 else a[i + 1] if mode == 1 else a[i] + (a[i + 1] - a[i]) * Fr(rnd.randint(0, 8), 8))
    beta = -a[0] + Fr(rnd.randint(1, 30), rnd.randint(1, 7))                                    # a_1 + beta > 0
    sigma = [1 / (x + beta) for x in a]; s = [1 / (x + beta) for x in alpha]
    assert all(sigma[i] >= s[i] >= sigma[i + 1] for i in range(n - 1))
    dist = coupling(sigma, s); cases += 1; support += len(dist); ties += len(set(a)) < n
    okC &= sum(w for w, _ in dist) == 1 and all(w > 0 and len(set(t)) == n - 1 for w, t in dist)
    for k in range(n):
        for T in itertools.combinations(range(n - 1), k):
            lhs = Fr(1);
            for i in T: lhs *= s[i]
            rhs = Fr(0)
            for w, t in dist:
                term = w
                for i in T: term *= sigma[t[i]]
                rhs += term
            okC &= lhs == rhs
    c = [GQ(rq(), rq()) for _ in range(n - 1)]
    left = GQ(1)
    for x in a: left = left * (x + beta)
    for i in range(n - 1): left = left * (GQ(alpha[i]) + c[i]) * (1 / (alpha[i] + beta))
    right = GQ(0)
    for w, t in dist:
        j = (set(range(n)) - set(t)).pop(); term = GQ(w * (a[j] + beta))
        for i in range(n - 1): term = term * (GQ(a[t[i]]) + c[i])
        right = right + term
    okI &= left == right
print("%d random rational instances, n = 2..7 (%d with repeated a_i); %d injections in the supports in total" % (cases, ties, support))
print("Lemma C: weights positive with sum 1, supported on injections, all product moments exact:", okC)
print("Lemma I: the scaled product equals the mixture of products of the conjecture, exactly:", okI)
print("EXACT CHECK PASSED" if okC and okI else "EXACT CHECK FAILED"); sys.exit(0 if okC and okI else 1)
