#!/usr/bin/env python3
"""Lifting frontier for the square-free Fourier principal-minor conjecture.

For a square-free M and a prime p not dividing M, H(M, p) is the statement that
every principal minor of the Fourier matrix F_M is nonzero modulo every prime
ideal of Z[zeta_M] above p.  The prime ideals above p correspond to the
irreducible factors h of the cyclotomic polynomial Phi_M modulo p, and reduction
modulo the ideal is F_p[x]/(h) with zeta = x.  H(M, p) is decided by the C
program pmscan (see pmscan.c), one run per factor.

By the lifting theorem of Caragea, Lee, Malikiosis and Pfander (arXiv
2505.24326v3, Theorem 1.4), P(M) and H(M, p) imply P(pM).  Starting from
Chebotarev's theorem for primes, this script walks through the square-free N in
increasing order and records, for each composite N, the first presentation
N = p * M with P(M) already established and H(M, p) true.

Usage:
    python3 frontier.py scan BOUND BIN_DIR OUT_JSON [BUDGET]
    python3 frontier.py pairs MAX_BASE MAX_PRIME BIN_DIR OUT_JSON [BUDGET]
    python3 frontier.py pair M p BIN_DIR [list]

BIN_DIR must contain the two builds of pmscan.c (pmscan, pmscan_poly).
Standard library only.
"""
import json
import os
import random
import subprocess
import sys
import time
from math import gcd

TABLE_LIMIT = 16_000_000  # largest field handled with log/Zech tables


# ---------- integers ----------
def factorint(n):
    out, d = {}, 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def is_squarefree(n):
    return all(e == 1 for e in factorint(n).values())


def mult_order(a, n):
    if n == 1:
        return 1
    k, x = 1, a % n
    while x != 1:
        x, k = x * a % n, k + 1
    return k


def phi(n):
    r = n
    for q in factorint(n):
        r = r // q * (q - 1)
    return r


# ---------- polynomials over F_p, low degree first ----------
def trim(a):
    while a and a[-1] == 0:
        a.pop()
    return a


def padd(a, b, p):
    n = max(len(a), len(b))
    return trim([((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p for i in range(n)])


def psub(a, b, p):
    n = max(len(a), len(b))
    return trim([((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)) % p for i in range(n)])


def pmul(a, b, p):
    if not a or not b:
        return []
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] = (r[i + j] + x * y) % p
    return trim(r)


def pdivmod(a, b, p):
    a = a[:]
    inv = pow(b[-1], p - 2, p)
    q = [0] * max(0, len(a) - len(b) + 1)
    while len(a) >= len(b) and a:
        c = a[-1] * inv % p
        s = len(a) - len(b)
        q[s] = c
        for i, y in enumerate(b):
            a[s + i] = (a[s + i] - c * y) % p
        trim(a)
    return trim(q), a


def pgcd(a, b, p):
    while b:
        a, b = b, pdivmod(a, b, p)[1]
    if a:
        inv = pow(a[-1], p - 2, p)
        a = [x * inv % p for x in a]
    return a


def ppowmod(a, e, m, p):
    r, a = [1], pdivmod(a, m, p)[1]
    while e:
        if e & 1:
            r = pdivmod(pmul(r, a, p), m, p)[1]
        a = pdivmod(pmul(a, a, p), m, p)[1]
        e >>= 1
    return r


# ---------- cyclotomic polynomials ----------
_CYC = {}


def cyclotomic(n):
    """Phi_n over Z, low degree first."""
    if n in _CYC:
        return _CYC[n]
    num = [-1] + [0] * (n - 1) + [1]
    for d in range(1, n):
        if n % d == 0:
            den = cyclotomic(d)
            q = [0] * (len(num) - len(den) + 1)
            num = num[:]
            for s in range(len(q) - 1, -1, -1):
                c = num[s + len(den) - 1]  # den is monic
                q[s] = c
                for i, y in enumerate(den):
                    num[s + i] -= c * y
            assert not any(num)
            num = q
    _CYC[n] = num
    return num


def equal_degree_split(h, d, p, rng):
    n = len(h) - 1
    if n == d:
        return [h]
    while True:
        a = trim([rng.randrange(p) for _ in range(n)])
        if len(a) < 2:
            continue
        if p == 2:
            b, t = a[:], a[:]
            for _ in range(d - 1):
                t = pdivmod(pmul(t, t, p), h, p)[1]
                b = padd(b, t, p)
        else:
            b = psub(ppowmod(a, (p ** d - 1) // 2, h, p), [1], p)
        g = pgcd(b, h, p)
        if 1 < len(g) < len(h):
            return equal_degree_split(g, d, p, rng) + equal_degree_split(pdivmod(h, g, p)[0], d, p, rng)


def cyclotomic_factors(M, p):
    """Irreducible factors of Phi_M mod p.  They all have degree ord_M(p)."""
    assert gcd(M, p) == 1
    f = mult_order(p, M)
    phi_m = trim([c % p for c in cyclotomic(M)])
    rng = random.Random(1000003 * M + p)
    facs = sorted(equal_degree_split(phi_m, f, p, rng))
    prod = [1]
    for h in facs:
        assert len(h) == f + 1 and h[-1] == 1
        prod = pmul(prod, h, p)
    assert prod == phi_m and len(facs) * f == phi(M)
    return f, facs


def chart_classes(M, p, f, facs):
    """Group the factors: h and h' are equivalent when a root of h' is a root of h
    raised to a unit square.  det(rho^(u^2 a b))_{a,b in S} = det(rho^(ab))_{a,b in uS},
    so equivalent charts have zero sets that differ by scaling the index set."""
    h1 = facs[0]

    def which(s):
        xs = ppowmod([0, 1], s, h1, p)
        for j, h in enumerate(facs):
            acc = []
            for c in reversed(h):  # Horner in F_p[x]/(h1)
                acc = padd(pdivmod(pmul(acc, xs, p), h1, p)[1], [c], p)
            if not acc:
                return j
        raise AssertionError("root not found")

    units = [s for s in range(1, M) if gcd(s, M) == 1]
    factor_of = {s: which(s) for s in units}
    cls = list(range(len(facs)))

    def find(i):
        while cls[i] != i:
            i = cls[i]
        return i

    for s in units:
        for u in units:
            a, b = find(factor_of[s]), find(factor_of[s * u * u % M])
            if a != b:
                cls[max(a, b)] = min(a, b)
    groups = {}
    for j in range(len(facs)):
        groups.setdefault(find(j), []).append(j)
    return list(groups.values())


# ---------- running the scanner ----------
def run_chart(M, p, f, h, bin_dir, extra=()):
    exe = "pmscan" if p ** f <= TABLE_LIMIT else "pmscan_poly"
    cmd = [os.path.join(bin_dir, exe), str(M), str(p), str(f)] + [str(c) for c in h[:-1]] + list(extra)
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError("pmscan failed: %s" % out.stderr)
    res = {"h": h, "backend": exe, "zero_masks": []}
    for line in out.stdout.splitlines():
        if line.startswith("Z "):
            res["zero_masks"].append(int(line.split()[1], 16))
        elif line.startswith("FIRSTZERO"):
            res["first_zero"] = int(line.split()[1], 16)
            res["zeros"] = None
        elif line.startswith("RESULT"):
            kv = dict(tok.split("=", 1) for tok in line.split()[1:])
            res["subsets"] = int(kv["subsets"])
            res["zeros"] = int(kv["zeros"])
            res["byrank"] = [int(x) for x in kv["byrank"].split(",")]
            assert res["subsets"] == 2 ** M
    return res


def cost(M, p):
    f = mult_order(p, M)
    per = 1 if p ** f <= TABLE_LIMIT else 3 * f * f
    return (2 ** M) * per


def test_pair(M, p, bin_dir, full=True, all_charts=True, list_zeros=False):
    """Decide H(M, p).  Returns a dict with holds = True/False."""
    f, facs = cyclotomic_factors(M, p)
    classes = chart_classes(M, p, f, facs)
    todo = range(len(facs)) if all_charts else [c[0] for c in classes]
    extra = (["list"] if list_zeros else []) + ([] if full else ["first"])
    t0 = time.time()
    charts = []
    for j in todo:
        r = run_chart(M, p, f, facs[j], bin_dir, extra)
        r["factor_index"] = j
        charts.append(r)
        if not full and "first_zero" in r:
            break
    holds = all(c.get("zeros") == 0 for c in charts) and len(charts) == len(list(todo))
    return {"M": M, "p": p, "f": f, "n_factors": len(facs), "classes": classes, "all_charts": all_charts,
            "full": full, "holds": holds, "charts": charts, "seconds": round(time.time() - t0, 2)}


class SearchLimit(Exception):
    pass


def covering_chain(bad, maxlen, work_limit=200_000_000):
    """bad[i] is the list of index sets (bit masks) whose principal minor vanishes in
    chart i.  Return a chain of at most maxlen nested sets that meets every bad[i],
    or None if there is none.  If there is none, every occupancy flag avoids all
    singular carriers in at least one chart, and Lemma 4.1 of Wang-Zhang
    (arXiv 2608.22119v1) gives P(M p).  Raises SearchLimit when the exhaustive
    search would exceed work_limit comparisons."""
    order = sorted(range(len(bad)), key=lambda i: len(bad[i]))
    bad_sets = [set(b) for b in bad]
    work = [0]

    def go(k, chain):
        while k < len(order) and any(c in bad_sets[order[k]] for c in chain):
            k += 1
        if k == len(order):
            return chain
        if len(chain) == maxlen:
            return None
        cand = bad[order[k]]
        work[0] += len(cand) * max(1, len(chain))
        if work[0] > work_limit:
            raise SearchLimit()
        for s in cand:
            if all((s & c) == s or (s & c) == c for c in chain):
                r = go(k + 1, chain + [s])
                if r is not None:
                    return r
        return None

    if any(not b for b in bad):
        return None
    return go(0, [])


def selective(M, p, bin_dir):
    """Selective lifting test for the presentation (p, M).  Returns the dict of
    test_pair with two extra keys: covering_chain (None if no chain covers every
    chart) and selective_holds."""
    r = test_pair(M, p, bin_dir, full=True, all_charts=True, list_zeros=True)
    bad = [c["zero_masks"] for c in r["charts"]]
    try:
        chain = covering_chain(bad, min(p, M))
        r["covering_chain"] = chain
        r["selective_holds"] = chain is None
    except SearchLimit:
        r["covering_chain"] = None
        r["selective_holds"] = None
    return r


LIST_LIMIT = 22  # bases up to this size are enumerated in full straight away
SEL_LIMIT = 30   # larger bases get the selective test only if the ordinary test fails


def decide_pair(M, p, bin_dir, budget):
    """Outcome of the two first-order tests for the presentation N = p * M."""
    def full():
        r = selective(M, p, bin_dir)
        return r, {"mode": "full enumeration, all charts", "ordinary": r["holds"], "selective": r["selective_holds"],
                   "zeros_by_chart": [c["zeros"] for c in r["charts"]],
                   "byrank": [c["byrank"] for c in r["charts"]],
                   "covering_chain": r["covering_chain"]}

    if M <= LIST_LIMIT:
        r, out = full()
    else:
        r = test_pair(M, p, bin_dir, full=False, all_charts=False)
        out = {"mode": "stop at first zero, one chart per square class", "ordinary": r["holds"], "selective": None,
               "first_zero": [c.get("first_zero") for c in r["charts"]]}
        quick_seconds = r["seconds"]
        if not r["holds"] and M <= SEL_LIMIT and cost(M, p) * r["n_factors"] <= budget:
            r, out = full()
            r["seconds"] += quick_seconds
    out.update({"f": r["f"], "n_factors": r["n_factors"], "classes": r["classes"],
                "factors": [c["h"] for c in r["charts"]], "backend": r["charts"][0]["backend"], "seconds": r["seconds"]})
    return out


def scan(bound, bin_dir, out_json, budget):
    """Walk through the square-free N <= bound.  status[N]["proved"] is True when N is
    prime, or some presentation N = p * M passes the ordinary test with P(M) already
    proved, or passes the selective test (which needs no hypothesis on M)."""
    cache = {}
    if os.path.exists(out_json):
        cache = {tuple(json.loads(k)): v for k, v in json.load(open(out_json)).get("pairs", {}).items()}
    status = {1: {"proved": True, "how": "trivial"}}
    for N in range(2, bound + 1):
        if not is_squarefree(N):
            continue
        primes = sorted(factorint(N))
        if len(primes) == 1:
            status[N] = {"proved": True, "how": "prime (Chebotarev)"}
            continue
        tried = []
        entry = {"proved": False, "tried": tried}
        for p in sorted(primes, reverse=True):  # largest prime first: smallest base
            M = N // p
            c = cost(M, p)
            if c > budget:
                tried.append({"p": p, "M": M, "outcome": "not run (cost %.1e over budget)" % c})
                continue
            if (M, p) not in cache:
                cache[(M, p)] = decide_pair(M, p, bin_dir, budget)
                save(out_json, bound, budget, cache, status)
            r = cache[(M, p)]
            if r["ordinary"] and status[M]["proved"]:
                tried.append({"p": p, "M": M, "outcome": "ordinary lifting holds"})
                entry.update({"proved": True, "how": "ordinary", "p": p, "M": M})
                break
            if r["selective"]:
                tried.append({"p": p, "M": M, "outcome": "selective lifting holds"})
                entry.update({"proved": True, "how": "selective", "p": p, "M": M})
                break
            if r["ordinary"]:
                tried.append({"p": p, "M": M, "outcome": "ordinary test holds but P(%d) is not established" % M})
            elif r["selective"] is None:
                tried.append({"p": p, "M": M, "outcome": "ordinary fails; selective not decided"})
            else:
                tried.append({"p": p, "M": M, "outcome": "ordinary and selective fail"})
        status[N] = entry
        if entry["proved"]:
            print("%d = %s: %s lifting with p = %d over base %d" % (N, "*".join(map(str, primes)), entry["how"], entry["p"], entry["M"]), flush=True)
        else:
            print("%d = %s: NOT SETTLED; %s" % (N, "*".join(map(str, primes)), "; ".join("(%d,%d) %s" % (t["p"], t["M"], t["outcome"]) for t in tried)), flush=True)
    save(out_json, bound, budget, cache, status)
    return status


def extra_pairs(max_base, max_prime, bin_dir, out_json, budget):
    """Decide every pair (M, p) with M <= max_base square-free, p <= max_prime prime,
    p not dividing M, and merge the outcomes into the cache of out_json.  These
    small pairs feed the propagation rule: a vanishing minor of F_M modulo a prime
    above p gives one of F_(M m), so a failed test at (p, M) fails at (p, M m)."""
    data = json.load(open(out_json))
    cache = {tuple(json.loads(k)): v for k, v in data.get("pairs", {}).items()}
    status = {int(k): v for k, v in data["status"].items()}
    for M in range(2, max_base + 1):
        if not is_squarefree(M):
            continue
        for q in range(2, max_prime + 1):
            if len(factorint(q)) != 1 or factorint(q)[q] != 1 or M % q == 0:
                continue
            if (M, q) in cache or cost(M, q) > budget:
                continue
            cache[(M, q)] = decide_pair(M, q, bin_dir, budget)
            r = cache[(M, q)]
            print("(M,p)=(%d,%d): ordinary %s, selective %s" % (M, q, r["ordinary"], r["selective"]), flush=True)
            save(out_json, data["bound"], data["budget"], cache, status)
    save(out_json, data["bound"], data["budget"], cache, status)


def save(out_json, bound, budget, cache, status):
    data = {"bound": bound, "budget": budget,
            "pairs": {json.dumps(list(k)): v for k, v in sorted(cache.items())},
            "status": {str(k): v for k, v in sorted(status.items())}}
    tmp = out_json + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(data, fh, indent=1, sort_keys=True)
    os.replace(tmp, out_json)


if __name__ == "__main__":
    if len(sys.argv) >= 5 and sys.argv[1] == "scan":
        scan(int(sys.argv[2]), sys.argv[3], sys.argv[4], float(sys.argv[5]) if len(sys.argv) > 5 else 3e9)
    elif len(sys.argv) >= 6 and sys.argv[1] == "pairs":
        extra_pairs(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5], float(sys.argv[6]) if len(sys.argv) > 6 else 3e9)
    elif len(sys.argv) >= 5 and sys.argv[1] == "pair":
        M, p = int(sys.argv[2]), int(sys.argv[3])
        r = test_pair(M, p, sys.argv[4], list_zeros="list" in sys.argv[5:])
        print("H(%d,%d): %s; f=%d, %d factor(s), classes %s, %.2fs" % (M, p, "holds" if r["holds"] else "FAILS", r["f"], r["n_factors"], r["classes"], r["seconds"]))
        for c in r["charts"]:
            print("  h=%s zeros=%s byrank=%s" % (c["h"], c["zeros"], c["byrank"]))
    else:
        print(__doc__)
        sys.exit(2)
