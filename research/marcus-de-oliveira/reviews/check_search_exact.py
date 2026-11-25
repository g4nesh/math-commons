#!/usr/bin/env python3
"""Exact re-examination of all 200 saved mdo-search samples (search/report.json @ 67528c1).

Review aid for mdo-search by vn-claude. The search stores binary64 data whose U is only
approximately unitary, so its 86 "outside the stored-float hull" flags cannot be judged
from the floats alone. This script turns every sample into a nearby EXACT instance and
decides membership exactly:

  1. a, b: the stored binary64 values, converted to Fractions without rounding.
  2. U: pick a scalar phase w (|w| = 1, from a 64-point grid) maximizing the distance of
     the spectrum of w*U from -1. Since d(wU) = d(U) (the LUR invariance proved in
     reductions §3 with L = wI, R = I), this changes nothing mathematically. Then take the
     binary64 Cayley parameter K = (I - wU)(I + wU)^-1, convert it to Fractions, take its
     exact skew-Hermitian part, and form U' = (I - K)(I + K)^-1 exactly. By reductions
     Lemma 1, U' is exactly unitary over Q(i), and it agrees with wU to about 1e-15.
  3. d = det(diag(a) + U' diag(b) U'*) and all n! permutation points, exactly, then an
     exact inside/outside certificate. Steps 3 use the core functions of
     certificates/verify.py (independently reviewed: reviews/mdo-certificates.md),
     called directly because Cayley denominators exceed its 64-digit file-input limit.
     Unitarity of U' is re-checked here exactly.

What this does and does not show: each result is an exact statement about the exact
instance (a, b, U'), which lies within about 1e-15 of the saved sample. It is not a
statement about the saved floats themselves, and it is not a search.
Bounds: the 200 saved cases only; no new random samples. numpy is used only to choose
the phase and the binary64 Cayley parameter.
Usage: python -I research/marcus-de-oliveira/reviews/check_search_exact.py
"""
import json
import sys
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "certificates"))
sys.dont_write_bytecode = True
import verify as V  # noqa: E402

Q = V.QComplex
REPORT = HERE.parent / "search" / "report.json"


def qc(z):
    return Q(F(z.real), F(z.imag))


def exact_unitary_near(U):
    n = U.shape[0]
    eig = np.linalg.eigvals(U)
    phases = np.exp(2j * np.pi * np.arange(64) / 64)
    w = max(phases, key=lambda p: np.min(np.abs(p * eig + 1)))
    W = w * U
    I = np.eye(n)
    Kf = (I - W) @ np.linalg.inv(I + W)
    K = [[qc(Kf[i, j]) for j in range(n)] for i in range(n)]
    half = F(1, 2)
    K = [[Q(half * (K[i][j].real - K[j][i].real), half * (K[i][j].imag + K[j][i].imag)) for j in range(n)]
         for i in range(n)]
    one = V.ONE
    IpK = [[(one if i == j else V.ZERO) + K[i][j] for j in range(n)] for i in range(n)]
    ImK = [[(one if i == j else V.ZERO) - K[i][j] for j in range(n)] for i in range(n)]
    M = [row[:] + [one if i == j else V.ZERO for j in range(n)] for i, row in enumerate(IpK)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != V.ZERO)
        M[c], M[p] = M[p], M[c]
        piv = M[c][c]
        M[c] = [x / piv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != V.ZERO:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    inv = [row[n:] for row in M]
    Up = V.multiply(ImK, inv)
    star = V.adjoint(Up)
    assert V.multiply(star, Up) == V.identity(n) and V.multiply(Up, star) == V.identity(n)
    dist = max(abs(complex(float(Up[i][j].real), float(Up[i][j].imag)) - W[i, j]) for i in range(n) for j in range(n))
    return Up, dist


def rel_slack(points, target, cert):
    """Exact inward slack (inside) or outward gap (outside), in Euclidean units / scale."""
    scale = max([1.0] + [abs(complex(float(x), float(y))) for x, y in points])
    hull = cert["hull"]
    if cert["membership"] == "outside":
        nx, ny = cert["normal"]
        return float(cert["gap"]) / float((nx * nx + ny * ny)) ** 0.5 / scale
    if len(hull) < 3:
        return 0.0  # on a point or segment: membership is exact equality / collinearity
    best = None
    for i in range(len(hull)):
        a, b = hull[i], hull[(i + 1) % len(hull)]
        cr = V.cross(a, b, target)
        length = float((b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2) ** 0.5
        s = float(cr) / length
        best = s if best is None else min(best, s)
    return best / scale


report = json.loads(REPORT.read_text(encoding="utf-8"))
rows, counts, outside = [], Counter(), []
max_dist = 0.0
for case in report["cases"]:
    a = [qc(complex(*z)) for z in case["a"]]
    b = [qc(complex(*z)) for z in case["b"]]
    U = np.array([[complex(*z) for z in row] for row in case["U"]])
    Up, dist = exact_unitary_near(U)
    max_dist = max(max_dist, dist)
    n = len(a)
    B = V.multiply(V.multiply(Up, V.diagonal(b)), V.adjoint(Up))
    S = [[(a[i] if i == j else V.ZERO) + B[i][j] for j in range(n)] for i in range(n)]
    d = V.determinant(S)
    pts = {z.pair() for _, z in V.permutation_points(a, b)}
    cert = V.certify_hull(list(pts), d.pair())  # includes the direct exact re-check
    dim = min(2, len(cert["hull"]) - 1)
    slack = rel_slack(list(pts), d.pair(), cert)
    counts[(dim, cert["membership"])] += 1
    float_flag = case["analysis"]["hull"]["outside_stored_float_hull"]
    rows.append({"id": case["id"], "hullDimension": dim, "membership": cert["membership"],
                 "relativeExactSlack": slack, "storedFloatOutside": float_flag,
                 "classification": case["analysis"]["classification"]})
    if cert["membership"] == "outside":
        outside.append(case["id"])

print(f"cases: {len(rows)}; max |U' - wU| entry distance: {max_dist:.2e}")
for (dim, mem), k in sorted(counts.items()):
    print(f"   hull dimension {dim}: {mem:7s} {k}")
flagged = [r for r in rows if r["storedFloatOutside"]]
print(f"stored-float 'outside' flags: {len(flagged)}; exact membership of those instances: "
      f"{dict(Counter(r['membership'] for r in flagged))}")
poly = sorted((r for r in rows if r["hullDimension"] == 2), key=lambda r: r["relativeExactSlack"])
print("tightest polygon instances (smallest exact relative inward slack):")
for r in poly[:8]:
    print(f"   {r['id']:28s} slack {r['relativeExactSlack']:.3e}  search label {r['classification']}"
          f"{'  [float-outside]' if r['storedFloatOutside'] else ''}")
out = HERE / "search-exact-recheck.json"
out.write_text(json.dumps({"source": "search/report.json @ 67528c1", "maxUnitaryDistance": max_dist,
                           "outside": outside, "cases": rows}, indent=2) + "\n", encoding="utf-8", newline="\n")
print(f"exact outside results: {len(outside)} {outside if outside else ''}")
print(f"wrote {out.name}")
