#!/usr/bin/env python3
"""Exact check of Kovacec's former Conjecture 7 (polygon enclosure) for one witness.

Conjecture 7 (DMUC p2611 v1, withdrawn in v2) claimed  Delta(A,B) is contained in cl(J(A,B)),
where J(A,B) is the bounded region enclosed by the outer boundary of the metric graph
whose vertices are the sigma-points z_sigma and whose edges are the segments
[z_sigma, z_tau] for permutations sigma, tau differing by one transposition.

For a witness (a, b, U) this script computes, in exact Gaussian-rational arithmetic:

  1. unitarity of U and the determinant p = det(D_a + U D_b U^*);
  2. the ORIGINAL Marcus-de Oliveira test: is p in conv{z_sigma}?  (facet slacks);
  3. the STRONGER test: is p outside cl(J)?  Certified by an explicit ray from p
     that meets no graph edge, so p lies in the unbounded component of the
     complement of the graph;
  4. the outer boundary walk of the graph (display aid; neighbour order uses float angles).

The two verdicts are reported separately on purpose: a point outside cl(J) but inside
the convex hull refutes Conjecture 7 and says nothing against the original conjecture.

Standard library only. Written independently of ../certificates/verify.py.
Usage:  python adversarial/polygon_enclosure.py certificates/fixtures/goyal-2026.json [...]
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction as F
from itertools import permutations
from pathlib import Path
import math

MAX_N = 4  # n! vertices and the O(E^2) arrangement stay small; larger n is out of scope.


# ---------- Gaussian rationals as (re, im) tuples of Fractions ----------
def cadd(x, y): return (x[0] + y[0], x[1] + y[1])
def csub(x, y): return (x[0] - y[0], x[1] - y[1])
def cmul(x, y): return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])
def cconj(x): return (x[0], -x[1])
def cdiv(x, y):
    d = y[0] * y[0] + y[1] * y[1]
    if d == 0:
        raise ZeroDivisionError
    n = cmul(x, cconj(y))
    return (n[0] / d, n[1] / d)
ZERO, ONE = (F(0), F(0)), (F(1), F(0))


def parse(v):
    return (F(v[0]), F(v[1]))


def det(m):
    """Gaussian elimination over Q(i) with exact pivots."""
    m = [row[:] for row in m]
    n, sign, acc = len(m), 1, ONE
    for c in range(n):
        piv = next((r for r in range(c, n) if m[r][c] != ZERO), None)
        if piv is None:
            return ZERO
        if piv != c:
            m[c], m[piv] = m[piv], m[c]
            sign = -sign
        acc = cmul(acc, m[c][c])
        for r in range(c + 1, n):
            f = cdiv(m[r][c], m[c][c])
            m[r] = [csub(m[r][k], cmul(f, m[c][k])) for k in range(n)]
    return acc if sign == 1 else (-acc[0], -acc[1])


# ---------- exact planar geometry on (x, y) Fraction points ----------
def cross(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def on_segment(p, a, b):
    return cross(a, b, p) == 0 and min(a[0], b[0]) <= p[0] <= max(a[0], b[0]) \
        and min(a[1], b[1]) <= p[1] <= max(a[1], b[1])


def ray_hits_segment(p, d, a, b):
    """Does the closed ray {p + t d : t >= 0} meet the closed segment [a, b]?"""
    e = (b[0] - a[0], b[1] - a[1])
    den = d[0] * e[1] - d[1] * e[0]
    ap = (a[0] - p[0], a[1] - p[1])
    if den == 0:
        if ap[0] * d[1] - ap[1] * d[0] != 0:
            return False  # parallel, not collinear
        # collinear: hit iff some endpoint has nonnegative projection on d
        return any((q[0] - p[0]) * d[0] + (q[1] - p[1]) * d[1] >= 0 for q in (a, b))
    t = (ap[0] * e[1] - ap[1] * e[0]) / den
    s = (ap[0] * d[1] - ap[1] * d[0]) / den
    return t >= 0 and 0 <= s <= 1


def seg_intersections(a, b, c, d):
    """Points of [a,b] that also lie on [c,d] which matter for splitting [a,b]."""
    out = [q for q in (c, d) if on_segment(q, a, b)]
    r = (b[0] - a[0], b[1] - a[1])
    s = (d[0] - c[0], d[1] - c[1])
    den = r[0] * s[1] - r[1] * s[0]
    if den != 0:
        ca = (c[0] - a[0], c[1] - a[1])
        t = (ca[0] * s[1] - ca[1] * s[0]) / den
        u = (ca[0] * r[1] - ca[1] * r[0]) / den
        if 0 <= t <= 1 and 0 <= u <= 1:
            out.append((a[0] + t * r[0], a[1] + t * r[1]))
    return out


def convex_hull(pts):
    pts = sorted(set(pts))
    if len(pts) <= 2:
        return pts
    lo, hi = [], []
    for p in pts:
        while len(lo) >= 2 and cross(lo[-2], lo[-1], p) <= 0:
            lo.pop()
        lo.append(p)
    for p in reversed(pts):
        while len(hi) >= 2 and cross(hi[-2], hi[-1], p) <= 0:
            hi.pop()
        hi.append(p)
    return lo[:-1] + hi[:-1]


def in_hull(p, hull):
    if len(hull) == 1:
        return p == hull[0], None
    if len(hull) == 2:
        return on_segment(p, hull[0], hull[1]), None
    slacks = [cross(hull[i], hull[(i + 1) % len(hull)], p) for i in range(len(hull))]
    return all(s >= 0 for s in slacks), min(slacks)


def outer_walk(segments):
    """Outer face boundary of the arrangement of the given segments (exact)."""
    pieces = set()
    for i, (a, b) in enumerate(segments):
        cuts = {a, b}
        for j, (c, d) in enumerate(segments):
            if i != j:
                cuts.update(seg_intersections(a, b, c, d))
        key = (lambda q: (q[0] - a[0]) * (b[0] - a[0]) + (q[1] - a[1]) * (b[1] - a[1]))
        ordered = sorted(cuts, key=key)
        for u, v in zip(ordered, ordered[1:]):
            if u != v:
                pieces.add(frozenset((u, v)))
    adj = {}
    for e in pieces:
        u, v = tuple(e)
        adj.setdefault(u, set()).add(v)
        adj.setdefault(v, set()).add(u)

    def next_vertex(prev, cur):
        ref = (prev[0] - cur[0], prev[1] - cur[1])
        cands = [w for w in adj[cur] if w != prev] or [prev]
        # first edge met when rotating counterclockwise from the incoming direction.
        # Float angles are used only to ORDER exact neighbours for this display walk;
        # the outside-cl(J) certificate is the exact ray test, not this walk.
        def angle_key(w):
            v = (w[0] - cur[0], w[1] - cur[1])
            a = (math.atan2(float(v[1]), float(v[0])) - math.atan2(float(ref[1]), float(ref[0]))) % (2 * math.pi)
            return a
        cands.sort(key=angle_key)
        return cands[0]

    start = min(adj)  # leftmost, then lowest: the region x < start.x is exterior
    virtual = (start[0] - 1, start[1])
    walk, prev, cur = [start], virtual, start
    first = None
    for _ in range(10 * len(pieces) + 10):
        nxt = next_vertex(prev, cur)
        if first is None:
            first = (cur, nxt)
        elif (cur, nxt) == first:
            break
        walk.append(nxt)
        prev, cur = cur, nxt
    return walk[:-1], len(pieces)


def winding(p, poly):
    w = 0
    for i in range(len(poly)):
        a, b = poly[i], poly[(i + 1) % len(poly)]
        if a[1] <= p[1] < b[1] and cross(a, b, p) > 0:
            w += 1
        elif b[1] <= p[1] < a[1] and cross(a, b, p) < 0:
            w -= 1
    return w


def fmt(q):
    return f"({q[0]}, {q[1]})"


def check(path):
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    a = [parse(v) for v in data["a"]]
    b = [parse(v) for v in data["b"]]
    U = [[parse(v) for v in row] for row in data["U"]]
    n = len(a)
    if not (1 <= n <= MAX_N and len(b) == n and len(U) == n and all(len(r) == n for r in U)):
        raise SystemExit(f"{path}: need square data with 1 <= n <= {MAX_N}")
    Ustar = [[cconj(U[j][i]) for j in range(n)] for i in range(n)]
    UU = [[ZERO for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            s = ZERO
            for k in range(n):
                s = cadd(s, cmul(Ustar[i][k], U[k][j]))
            UU[i][j] = s
    unitary = all(UU[i][j] == (ONE if i == j else ZERO) for i in range(n) for j in range(n))
    M = [[ZERO for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            s = a[i] if i == j else ZERO
            for k in range(n):
                s = cadd(s, cmul(cmul(U[i][k], b[k]), Ustar[k][j]))
            M[i][j] = s
    p = det(M)

    perms = list(permutations(range(n)))
    z = {}
    for s in perms:
        v = ONE
        for i in range(n):
            v = cmul(v, cadd(a[i], b[s[i]]))
        z[s] = v
    edges = {frozenset((s, t)) for s in perms for t in perms
             if sum(x != y for x, y in zip(s, t)) == 2}
    segs = sorted({tuple(sorted((z[s], z[t]))) for s, t in map(tuple, edges) if z[s] != z[t]})

    hull = convex_hull(list(z.values()))
    inside_moc, min_slack = in_hull(p, hull)

    on_graph = any(on_segment(p, u, v) for u, v in segs)
    ray = None
    if not on_graph:
        for r in range(1, 25):
            for dx in range(-r, r + 1):
                for dy in (-r, r) if abs(dx) != r else range(-r, r + 1):
                    d = (F(dx), F(dy))
                    if d != (0, 0) and not any(ray_hits_segment(p, d, u, v) for u, v in segs):
                        ray = d
                        break
                if ray:
                    break
            if ray:
                break

    walk, npieces = outer_walk(segs)
    wind = winding(p, walk)

    print(f"== {data.get('name', path)}  (n={n}, {len(edges)} transposition edges, {len(segs)} distinct segments)")
    print(f"U unitary exactly: {unitary}")
    print(f"det(D_a + U D_b U*) = {fmt(p)}")
    exp = data.get("expectedDeterminant")
    if exp:
        print(f"matches fixture expectedDeterminant: {p == parse(exp)}")
    print("sigma-points (one-line word sigma(1)..sigma(n), 1-based):")
    for s in perms:
        print(f"   z_{''.join(str(i + 1) for i in s)} = {fmt(z[s])}")
    print(f"ORIGINAL MOC: p in conv(z): {inside_moc}"
          + (f"  (min facet slack {min_slack}; hull has {len(hull)} vertices)" if min_slack is not None else ""))
    print(f"STRONGER (Conjecture 7): p on a graph edge: {on_graph}")
    if ray:
        print(f"   ray p + t*{fmt(ray)}, t>=0, meets no edge -> p is in the unbounded component, i.e. OUTSIDE cl(J)")
    else:
        print("   no avoiding ray found among the searched directions; outside cl(J) NOT certified")
    print(f"   outer boundary walk ({npieces} arrangement pieces, {len(walk)} vertices): winding number of p = {wind}")
    print("   walk: " + " -> ".join(fmt(q) for q in walk))
    refutes_c7 = bool(ray) and inside_moc
    print(f"VERDICT: Conjecture 7 refuted by this witness: {refutes_c7}; "
          f"original MOC refuted by this witness: {not inside_moc}")
    return {"name": data.get("name"), "unitary": unitary, "det": [str(p[0]), str(p[1])],
            "insideHull": inside_moc, "minFacetSlack": None if min_slack is None else str(min_slack),
            "onGraph": on_graph, "avoidingRay": None if not ray else [str(ray[0]), str(ray[1])],
            "outerWalk": [[str(q[0]), str(q[1])] for q in walk], "winding": wind,
            "refutesConjecture7": refutes_c7, "refutesMOC": not inside_moc}


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    results = [check(arg) for arg in sys.argv[1:]]
    out = Path(__file__).with_name("polygon-enclosure-results.json")
    out.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {out.name}")
