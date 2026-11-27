#!/usr/bin/env python3
"""Render the JSON written by `frontier.py scan` as Markdown.

Usage: python3 make_table.py scan.json > table.md

Each composite square-free order is tagged with the earliest located source that
already proves it (see ../../order-154/prior-art/ganesh-claude-minors/README.md
for the literature search and its limits):

  CLMP    arXiv 2505.24326v3, Theorem 1.3: 2p, 3p, 5p, 6p, 7p
  GZW     arXiv 2608.17746v1, Theorem 1.1: 70 and 143
  WZ      arXiv 2608.22119v1, Theorem 1.1: 10p for p not in {2, 5, 11}
  EPOCH   arXiv 2610.06986v1, Supplement K.5: 11p for p outside
          {3, 11, 23, 67, 199, 419}
  seed    claimed inside this repository only (order 154 and the 14p seed),
          awaiting independent review; not a literature source
"""
import json
import sys


def factorint(n):
    out, d = [], 2
    while d * d <= n:
        while n % d == 0:
            out.append(d)
            n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def is_prime(n):
    return n > 1 and len(factorint(n)) == 1


def literature(n):
    for k in (2, 3, 5, 6, 7):
        if n % k == 0 and is_prime(n // k):
            return "CLMP"
    if n in (70, 143):
        return "GZW"
    if n % 10 == 0 and is_prime(n // 10) and n // 10 not in (2, 5, 11):
        return "WZ"
    if n % 11 == 0 and is_prime(n // 11) and n // 11 not in (3, 11, 23, 67, 199, 419):
        return "EPOCH"
    if n % 14 == 0 and is_prime(n // 14) and n // 14 not in (2, 7, 29, 43):
        return "seed"
    return ""


def divisors(n):
    return [d for d in range(2, n + 1) if n % d == 0]


def star_obstruction(p, M):
    """Proposition 1 of README.md: pairwise coprime divisors d_1..d_k >= 2 of M with
    p | c, c = (k-1) prod d_j - sum_j prod_{i != j} d_i.  Returns the divisors or None."""
    primes = sorted(set(factorint(M)))
    best = None

    def go(i, blocks):
        nonlocal best
        if i == len(primes):
            ds = [d for d in blocks if d > 1]
            if len(ds) >= 2:
                prod = 1
                for d in ds:
                    prod *= d
                c = (len(ds) - 1) * prod - sum(prod // d for d in ds)
                if c % p == 0 and (best is None or prod < best[0]):
                    best = (prod, sorted(ds))
            return
        q = primes[i]
        go(i + 1, blocks)  # prime unused
        for j in range(len(blocks)):
            if blocks[j] > 1:
                go(i + 1, blocks[:j] + [blocks[j] * q] + blocks[j + 1:])
        go(i + 1, blocks + [q])

    go(0, [])
    return best[1] if best else None


def refine(p, M, outcome, pairs):
    """Apply the propagation rule to a presentation that was not decided directly:
    a vanishing minor of F_M0 modulo a prime above p gives one of F_M for M0 | M,
    and a covering chain at (p, M0) gives a covering chain at (p, M)."""
    if outcome.startswith("ordinary lifting holds") or outcome.startswith("selective lifting holds") or outcome == "ordinary and selective fail":
        return outcome
    own = pairs.get((M, p))
    if own and own["selective"] is False:
        return "ordinary and selective fail (run outside the scan, see special-runs.md)"
    sel = [d for d in divisors(M) if d < M and (d, p) in pairs and pairs[(d, p)]["selective"] is False]
    if sel:
        return "ordinary and selective fail (covering chain propagated from base %d)" % sel[0]
    star = star_obstruction(p, M)
    if star:
        return "ordinary and selective fail (Proposition 1 with subgroup orders %s)" % ", ".join(map(str, star))
    if outcome.startswith("ordinary fails"):
        return outcome
    od = [d for d in divisors(M) if d < M and (d, p) in pairs and pairs[(d, p)]["ordinary"] is False]
    if od:
        return "ordinary fails (vanishing minor propagated from base %d); selective not decided" % od[0]
    return outcome


def main(path):
    data = json.load(open(path))
    status = {int(k): v for k, v in data["status"].items()}
    pairs = {tuple(json.loads(k)): v for k, v in data["pairs"].items()}
    comp = [n for n in sorted(status) if n > 1 and not is_prime(n)]
    settled = [n for n in comp if status[n]["proved"]]
    open_ = [n for n in comp if not status[n]["proved"]]
    new = [n for n in settled if literature(n) in ("", "seed")]
    missing = [n for n in open_ if literature(n) not in ("", "seed")]
    assert not missing, "orders proved in the literature but not settled here: %s" % missing
    print("Bound %d. Composite square-free orders: %d. Settled: %d (%d by ordinary lifting, %d by selective lifting). Not settled: %d." % (
        data["bound"], len(comp), len(settled),
        sum(1 for n in settled if status[n]["how"] == "ordinary"),
        sum(1 for n in settled if status[n]["how"] == "selective"), len(open_)))
    print()
    print("Settled orders with no located literature source (%d): %s" % (len(new), ", ".join(
        "%d%s" % (n, "*" if literature(n) == "seed" else "") for n in new)))
    print()
    print("An asterisk marks orders also claimed by this repository's own seed (154 and 14p).")
    print()
    print("### Orders not settled")
    print()
    print("| N | factors | every presentation fails both tests? | presentations (p, base): outcome |")
    print("|---|---|---|---|")
    for n in open_:
        outs = [(t["p"], t["M"], refine(t["p"], t["M"], t["outcome"], pairs)) for t in status[n]["tried"]]
        closed = all(o.startswith("ordinary and selective fail") for _, _, o in outs)
        print("| %d | %s | %s | %s |" % (n, " · ".join(map(str, factorint(n))), "yes" if closed else "no",
                                     "; ".join("(%d, %d): %s" % o for o in outs)))
    print()
    print("### Pairs that fail the ordinary test (computed directly)")
    print()
    print("All pairs with base at most 22 and prime at most 31 were decided, plus the pairs the scan needed.")
    print("Counts are per chart; \"not counted\" means the run stopped at the first vanishing minor or was a targeted search.")
    print()
    print("| base M | prime p | f | charts | vanishing minors per chart | selective test |")
    print("|---|---|---|---|---|---|")
    for (M, p), r in sorted(pairs.items()):
        if r["ordinary"]:
            continue
        zeros = r.get("zeros_by_chart")
        ztxt = ", ".join(map(str, zeros)) if zeros else "not counted"
        sel = {True: "holds: no covering chain", False: "fails: covering chain exists", None: "not decided"}[r["selective"]]
        print("| %d | %d | %d | %d | %s | %s |" % (M, p, r["f"], r["n_factors"], ztxt, sel))
    print()
    print("### All composite square-free orders")
    print()
    print("| N | factors | settled by | located source |")
    print("|---|---|---|---|")
    for n in comp:
        st = status[n]
        how = "%s lifting, p = %d over base %d" % (st["how"], st["p"], st["M"]) if st["proved"] else "**not settled**"
        print("| %d | %s | %s | %s |" % (n, " · ".join(map(str, factorint(n))), how, literature(n)))


if __name__ == "__main__":
    main(sys.argv[1])
