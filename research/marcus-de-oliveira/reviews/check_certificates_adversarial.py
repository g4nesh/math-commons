#!/usr/bin/env python3
"""Adversarial differential tests of certificates/verify.py (commit 67528c1).

Review aid for mdo-certificates by vn-claude. The oracle code below is written
independently of verify.py and uses a different method for every decision:

  * hull membership: brute force over Caratheodory supports (a point is in the convex
    hull of a finite planar set iff it equals a point, lies on a segment, or lies in a
    closed triangle with vertices from the set). No hull is constructed.
  * determinant: Leibniz expansion over all permutations.
  * certificates: weights and separators returned by verify.py are re-checked here.

Parts
  A  geometry: random small rational point sets, biased toward degeneracies
     (single points, duplicates, collinear sets, targets on vertices and edges).
  B  end to end: random rational Cayley unitaries with Gaussian-integer spectra.
  C  malformed witness files: expected exit status 2 and no certificate.

Bounds: seed 20261125; A 3000 cases; B 150 witnesses with n in 2..4; C 14 files.
Usage: python -I research/marcus-de-oliveira/reviews/check_certificates_adversarial.py
"""
import json
import random
import subprocess
import sys
import tempfile
from fractions import Fraction as F
from itertools import combinations, permutations
from pathlib import Path

HERE = Path(__file__).resolve().parent
CERT = HERE.parent / "certificates"
sys.path.insert(0, str(CERT))
import verify as V  # noqa: E402  (the module under review; read in full before use)

SEED = 20261125
rng = random.Random(SEED)


# ---------------- independent oracle ----------------
def orient(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def on_seg(p, a, b):
    return orient(a, b, p) == 0 and min(a[0], b[0]) <= p[0] <= max(a[0], b[0]) \
        and min(a[1], b[1]) <= p[1] <= max(a[1], b[1])


def in_tri(p, a, b, c):
    s1, s2, s3 = orient(a, b, p), orient(b, c, p), orient(c, a, p)
    return (s1 >= 0 and s2 >= 0 and s3 >= 0) or (s1 <= 0 and s2 <= 0 and s3 <= 0)


def oracle_inside(points, p):
    pts = sorted(set(points))
    if p in pts:
        return True
    if any(on_seg(p, a, b) for a, b in combinations(pts, 2)):
        return True
    return any(orient(a, b, c) != 0 and in_tri(p, a, b, c) for a, b, c in combinations(pts, 3))


def recheck(points, p, cert):
    """Independent re-validation of a returned certificate."""
    pts = set(points)
    if cert["membership"] == "inside":
        w = cert["weights"]
        assert all(q in pts and x >= 0 for q, x in w) and sum(x for _, x in w) == 1
        assert (sum(q[0] * x for q, x in w), sum(q[1] * x for q, x in w)) == p
    else:
        nx, ny = cert["normal"]
        bound = max(nx * q[0] + ny * q[1] for q in pts)
        assert (nx, ny) != (0, 0) and nx * p[0] + ny * p[1] > bound


def rand_point(k=4):
    return (F(rng.randint(-k, k), rng.choice([1, 1, 2, 3])), F(rng.randint(-k, k), rng.choice([1, 1, 2, 3])))


def geometry_case():
    kind = rng.random()
    if kind < 0.15:                       # single point, possibly repeated
        q = rand_point()
        pts = [q] * rng.randint(1, 4)
    elif kind < 0.40:                     # collinear set with duplicates
        o, d = rand_point(), rand_point(2)
        pts = [(o[0] + t * d[0], o[1] + t * d[1]) for t in (F(rng.randint(-3, 3)) for _ in range(rng.randint(2, 6)))]
    else:                                 # general, with duplicates
        pts = [rand_point() for _ in range(rng.randint(3, 8))]
        pts += rng.sample(pts, rng.randint(0, 2))
    r = rng.random()
    if r < 0.2:
        p = rng.choice(pts)                # a vertex or duplicate
    elif r < 0.45:
        a, b = rng.choice(pts), rng.choice(pts)
        t = F(rng.randint(0, 6), 6)
        p = (a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]))  # on a chord or edge
    else:
        p = rand_point(5)
    return pts, p


# ---------------- Part A ----------------
a_cases = a_inside = a_outside = 0
mismatch = []
for _ in range(3000):
    pts, p = geometry_case()
    cert = V.certify_hull(pts, p)
    recheck(pts, p, cert)
    got = cert["membership"] == "inside"
    if got != oracle_inside(pts, p):
        mismatch.append((pts, p, cert["membership"]))
    a_cases += 1
    a_inside += got
    a_outside += not got
print(f"A geometry: {a_cases} cases ({a_inside} inside, {a_outside} outside), "
      f"oracle mismatches: {len(mismatch)}")


# ---------------- Part B ----------------
Q = V.QComplex
I1 = Q(F(1))


def gi(k=3):
    return Q(F(rng.randint(-k, k)), F(rng.randint(-k, k)))


def rational_cayley(n):
    K = [[V.ZERO] * n for _ in range(n)]
    for i in range(n):
        K[i][i] = Q(F(0), F(rng.randint(-2, 2)))
        for j in range(i + 1, n):
            K[i][j] = gi(2)
            K[j][i] = -K[i][j].conjugate()
    IpK = [[(I1 if i == j else V.ZERO) + K[i][j] for j in range(n)] for i in range(n)]
    ImK = [[(I1 if i == j else V.ZERO) - K[i][j] for j in range(n)] for i in range(n)]
    # invert I+K by Gauss-Jordan (oracle-side code, not verify.py)
    M = [row[:] + [I1 if i == j else V.ZERO for j in range(n)] for i, row in enumerate(IpK)]
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
    return [[sum((ImK[i][k] * inv[k][j] for k in range(n)), V.ZERO) for j in range(n)] for i in range(n)]


def leibniz(M):
    n, total = len(M), V.ZERO
    for s in permutations(range(n)):
        inv = sum(1 for i in range(n) for j in range(i + 1, n) if s[i] > s[j])
        term = I1
        for i in range(n):
            term = term * M[i][s[i]]
        total = total + (term if inv % 2 == 0 else -term)
    return total


def enc(z):
    return [str(z.real), str(z.imag)]


b_cases = b_outside = b_det_mismatch = b_member_mismatch = b_singular = 0
min_rel_slack = None
for _ in range(150):
    n = rng.randint(2, 4)
    a = [gi() for _ in range(n)]
    b = [gi() for _ in range(n)]
    if rng.random() < 0.2:                 # repeated eigenvalues
        a[-1] = a[0]
    if rng.random() < 0.15:                # force a singular sum direction
        b[0] = -a[0]
    U = rational_cayley(n)
    w = {"schemaVersion": 1, "a": [enc(x) for x in a], "b": [enc(x) for x in b],
         "U": [[enc(x) for x in row] for row in U]}
    res = V.verify_witness(w)
    Ustar = [[U[j][i].conjugate() for j in range(n)] for i in range(n)]
    B = [[sum((U[i][k] * b[k] * Ustar[k][j] for k in range(n)), V.ZERO) for j in range(n)] for i in range(n)]
    S = [[(a[i] if i == j else V.ZERO) + B[i][j] for j in range(n)] for i in range(n)]
    det = leibniz(S)
    if res["determinant"] != enc(det):
        b_det_mismatch += 1
    zs = []
    for s in permutations(range(n)):
        z = I1
        for i in range(n):
            z = z * (a[i] + b[s[i]])
        zs.append((z.real, z.imag))
    if (res["membership"] == "inside") != oracle_inside(zs, (det.real, det.imag)):
        b_member_mismatch += 1
    b_outside += res["membership"] == "outside"
    b_singular += res["singularSum"]
    b_cases += 1
print(f"B end-to-end: {b_cases} exact witnesses (n=2..4), {b_singular} singular sums, "
      f"determinant mismatches vs Leibniz: {b_det_mismatch}, membership mismatches vs oracle: "
      f"{b_member_mismatch}, outside results: {b_outside}")


# ---------------- Part C ----------------
good = json.loads((CERT / "fixtures" / "goyal-2026.json").read_text(encoding="utf-8"))


def variant(fn):
    w = json.loads(json.dumps(good))
    fn(w)
    return json.dumps(w)


cases = {
    "float entry": variant(lambda w: w["a"].__setitem__(0, [0.5, "0"])),
    "decimal string": variant(lambda w: w["a"].__setitem__(0, ["0.5", "0"])),
    "exponent string": variant(lambda w: w["a"].__setitem__(0, ["1e3", "0"])),
    "zero denominator": variant(lambda w: w["a"].__setitem__(0, ["1/0", "0"])),
    "negative denominator": variant(lambda w: w["a"].__setitem__(0, ["1/-2", "0"])),
    "padded string": variant(lambda w: w["a"].__setitem__(0, [" 1", "0"])),
    "65-digit integer": variant(lambda w: w["a"].__setitem__(0, ["1" * 65, "0"])),
    "bool schemaVersion": variant(lambda w: w.__setitem__("schemaVersion", True)),
    "unknown field": variant(lambda w: w.__setitem__("extra", 1)),
    "U off by 1e-30": variant(lambda w: w["U"][0].__setitem__(0, ["-1/3", "1/" + "1" + "0" * 30])),
    "wrong expected determinant": variant(lambda w: w.__setitem__("expectedDeterminant", ["2/3", "2/9 "])),
    "expected determinant mismatch": variant(lambda w: w.__setitem__("expectedDeterminant", ["2/3", "1/9"])),
    "NaN literal": json.dumps(good).replace('"schemaVersion": 1', '"schemaVersion": NaN'),
    "duplicate key": json.dumps(good)[:-1] + ', "a": [["0","0"],["1","0"],["0","1"]]}',
}
c_ok = 0
with tempfile.TemporaryDirectory() as tmp:
    for name, text in cases.items():
        path = Path(tmp) / "w.json"
        path.write_text(text, encoding="utf-8")
        proc = subprocess.run([sys.executable, "-I", str(CERT / "verify.py"), str(path)],
                              capture_output=True, text=True, timeout=60)
        rejected = proc.returncode == 2 and not proc.stdout.strip()
        c_ok += rejected
        if not rejected:
            print(f"   C NOT REJECTED as expected: {name} (exit {proc.returncode})")
    # Robustness probe, not a correctness requirement: deeply nested JSON
    path = Path(tmp) / "deep.json"
    path.write_text("[" * 100000 + "]" * 100000, encoding="utf-8")
    proc = subprocess.run([sys.executable, "-I", str(CERT / "verify.py"), str(path)],
                          capture_output=True, text=True, timeout=60)
    deep = proc.returncode
print(f"C malformed inputs: {c_ok}/{len(cases)} rejected with exit 2 and no certificate on stdout")
print(f"   deep-nesting probe (100000 levels, under 1 MiB): exit {deep} "
      f"({'clean rejection' if deep == 2 else 'uncaught exception, no certificate emitted'})")

ok = not mismatch and b_det_mismatch == 0 and b_member_mismatch == 0 and c_ok == len(cases)
print(f"seed={SEED}: {'no defects found' if ok else 'DEFECTS FOUND'}")
raise SystemExit(0 if ok else 1)
