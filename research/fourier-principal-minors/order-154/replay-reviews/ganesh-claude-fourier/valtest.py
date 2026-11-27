"""Test of the valuation formula in the proof of Theorem 1.4 (review.md, Lemma V).

For N = p N' square-free and S in Z_N, let n_a = |{s in S : s = a mod N'}| and A_r = {a : n_a >= r}. Lemma V says
  v_p(Nm D_S) >= phi(N') * sum_a C(n_a, 2),
with equality iff p does not divide Nm(D'_{A_r}) for every r >= 1, where D'_A is the principal minor of F_{N'} on A.
(Nm over Q(w_N), resp. Q(w_{N'}).) Checked exhaustively over all S for small N and every p | N,
and on random S for N in {30, 42, 154}.

usage: valtest.py SEED NRANDOM [N ...]   (N: orders to test exhaustively; default 6 10 14 15)
"""
import collections
import math
import random
import sys

from normcrt import all_subsets, norms, valuation


def phi(n):
    return sum(1 for c in range(1, n + 1) if math.gcd(c, n) == 1)


_cache = {}


def norm_small(Np, A):
    key = (Np, tuple(sorted(A)))
    if key not in _cache:
        _cache[key] = norms(Np, [sorted(A)])[0] if A else 1
    return _cache[key]


def check(N, subsets_by_size, stats):
    primes = [p for p in range(2, N + 1) if N % p == 0 and all(p % r for r in range(2, p))]
    for m, subs in subsets_by_size.items():
        nm = norms(N, subs)
        for S, x in zip(subs, nm):
            assert x != 0, (N, S)
            for p in primes:
                Np = N // p
                n = collections.Counter(s % Np for s in S)
                base = phi(Np) * sum(math.comb(k, 2) for k in n.values())
                chain = [[a for a in n if n[a] >= r] for r in range(1, max(n.values()) + 1)]
                good = all(norm_small(Np, A) % p for A in chain)
                v = valuation(x, p)
                key = (N, p)
                stats[key]["S"] += 1
                if v < base:
                    stats[key]["below bound"] += 1
                elif v == base:
                    stats[key]["equal, chain good" if good else "EQUAL BUT CHAIN BAD"] += 1
                else:
                    stats[key]["above, chain bad" if not good else "ABOVE BUT CHAIN GOOD"] += 1


def report(stats, keys):
    bad = 0
    for key in keys:
        c = stats[key]
        bad += c["below bound"] + c["EQUAL BUT CHAIN BAD"] + c["ABOVE BUT CHAIN GOOD"]
        print(f"  {key}: " + ", ".join(f"{k} {v}" for k, v in sorted(c.items())), flush=True)
    return bad


def main():
    seed, nrand = int(sys.argv[1]), int(sys.argv[2])
    exhaustive = [int(x) for x in sys.argv[3:]] or [6, 10, 14, 15]
    rng = random.Random(seed)
    stats = collections.defaultdict(collections.Counter)
    bad = 0
    for N in exhaustive:
        check(N, {m: all_subsets(N, m) for m in range(1, N + 1)}, stats)
        print(f"N = {N}: all {2 ** N - 1} nonempty subsets (N, p): counts", flush=True)
        bad += report(stats, [k for k in sorted(stats) if k[0] == N])
    for N, maxm in ((30, 12), (42, 12), (154, 12)) if nrand else ():
        by = collections.defaultdict(list)
        for _ in range(nrand):
            m = rng.randint(2, maxm)
            # bias toward repeated residues mod the small factors, so that the bound is not trivially 0
            base = rng.sample(range(N), rng.randint(1, m))
            S = set(base)
            while len(S) < m:
                b = rng.choice(base)
                S.add((b + rng.choice([d for d in range(2, N) if N % d == 0]) * rng.randrange(N)) % N)
            by[m].append(sorted(S))
        check(N, dict(by), stats)
        print(f"N = {N}: {nrand} random subsets (N, p): counts", flush=True)
        bad += report(stats, [k for k in sorted(stats) if k[0] == N])
    print("\nLemma V: no counterexample" if bad == 0 else f"\nLemma V: {bad} COUNTEREXAMPLES")


if __name__ == "__main__":
    main()
