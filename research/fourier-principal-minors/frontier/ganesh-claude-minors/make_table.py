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


def main(path):
    data = json.load(open(path))
    status = {int(k): v for k, v in data["status"].items()}
    pairs = {tuple(json.loads(k)): v for k, v in data["pairs"].items()}
    comp = [n for n in sorted(status) if n > 1 and not is_prime(n)]
    settled = [n for n in comp if status[n]["proved"]]
    open_ = [n for n in comp if not status[n]["proved"]]
    new = [n for n in settled if literature(n) in ("", "seed")]
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
    print("| N | factors | presentations tested (p, base): outcome |")
    print("|---|---|---|")
    for n in open_:
        tried = "; ".join("(%d, %d): %s" % (t["p"], t["M"], t["outcome"]) for t in status[n]["tried"])
        print("| %d | %s | %s |" % (n, " · ".join(map(str, factorint(n))), tried))
    print()
    print("### Presentations that fail the ordinary test")
    print()
    print("| base M | prime p | f | charts | vanishing minors per chart | selective test |")
    print("|---|---|---|---|---|---|")
    for (M, p), r in sorted(pairs.items()):
        if r["ordinary"]:
            continue
        zeros = r.get("zeros_by_chart")
        ztxt = ", ".join(map(str, zeros)) if zeros else "not counted (stopped at the first)"
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
