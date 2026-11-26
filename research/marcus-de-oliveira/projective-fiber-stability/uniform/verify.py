"""Exact full-support angular laws and rank at all ordinary projective parameters."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import permutations
from pathlib import Path
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(ROOT / "four-cycle-spectral"))
import spectral as sp
spec = importlib.util.spec_from_file_location("angular_law", ROOT / "complex-direction-cover/check.py")
cover = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cover)

EPSILON = F(1, 10_000_000)
PERMS = tuple(permutations(range(4)))


def strict_patches(patches, corrections):
    if len(patches) != 8 or len(corrections) != len(patches):
        raise ValueError("Expected eight angular patches and corrections")
    charts = defaultdict(list)
    for p, correction in zip(patches, corrections):
        if len(correction) != len(p["terms"]):
            raise ValueError("Correction and archived support disagree")
        top = max(len(centers) for _, _, centers, _ in p["terms"])
        if top not in (1, 2):
            raise ValueError("Unexpected radial degree")
        terms = defaultdict(F)
        for (j, word, centers, wt), v in zip(p["terms"], correction):
            mass = F(wt) - EPSILON * F(v)
            if mass <= 0:
                raise ValueError("Adjusted base term is not positive")
            terms[j, word, tuple(centers)] += mass
        for j in range(p["degree"] + 1):
            for perm in PERMS:
                word = "".join(str(i + 1) for i in perm)
                for centers in ((), ("0",) * top):
                    terms[j, word, centers] += EPSILON
        augmented = dict(p, terms=[[j, word, list(centers), str(mass)]
                                   for (j, word, centers), mass in terms.items()])
        chart, lo, hi = cover.verify_patch(augmented)
        charts[chart].append((lo, hi))
    edges = [F(-1), F(-1, 2), F(0), F(1, 2), F(1)]
    expected = list(zip(edges, edges[1:]))
    if set(charts) != {"h", "v"} or any(sorted(ranges) != expected for ranges in charts.values()):
        raise ArithmeticError("Incomplete angular coverage")
    print("All eight exact laws strictly support all 24 permutations, including their radial limits.")


def integer(z):
    if z.denominator != 1:
        raise ValueError("Grid arithmetic requires an integer")
    return z.numerator


def vertex_rows(z, reciprocal=False):
    # Columns are listed by the row-to-column word (i,j,l,k) here.
    out = []
    b = sp.B_SPECTRUM
    for i, j, l, k in PERMS:
        middle = (sp.O + z * b[j]) if reciprocal else (z + b[j])
        parts = [b[i] * b[l] * b[k], middle * b[l] * b[k], middle * b[i] * b[k], middle * b[k]]
        out.append([1] + [integer(getattr(v, part)) for part in ("real", "imag") for v in parts])
    return out


def determinant(rows):
    """Fraction-free determinant on the small integer matrices in the rank proof."""
    a = [r.copy() for r in rows]
    previous, sign = 1, 1
    for k in range(len(a) - 1):
        pivot = next((i for i in range(k, len(a)) if a[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        for i in range(k + 1, len(a)):
            for j in range(k + 1, len(a)):
                v = a[k][k] * a[i][j] - a[i][k] * a[k][j]
                a[i][j], rem = divmod(v, previous)
                if rem:
                    raise ArithmeticError("Fraction-free division failed")
        previous = a[k][k]
    return sign * a[-1][-1]


def check_basis(basis):
    if (len(basis) != 9 or any(type(i) is not int or not 0 <= i < 24 for i in basis)
            or len(set(basis)) != 9):
        raise ValueError("Invalid nine-permutation rank certificate")


def rank_identities(data):
    for basis in [*data["bases"], data["infinityBasis"]]:
        check_basis(basis)
    if len(data["bases"]) != 6 or len(data["multipliers"]) != 2:
        raise ValueError("Expected six determinants and two identities")
    identities = []
    for terms in data["multipliers"]:
        decoded = []
        for n, i, j, v in terms:
            if any(type(k) is not int for k in (n, i, j)) or not 0 <= n < 6 or min(i, j) < 0 or i+j > 2:
                raise ValueError("Invalid quadratic multiplier")
            decoded.append((n, i, j, F(v)))
        identities.append(decoded)
    # The products of degree-six determinants and quadratic multipliers have
    # total degree at most eight. The Newton triangular grid is unisolvent.
    for x in range(9):
        for y in range(9 - x):
            rows = vertex_rows(sp.Q(x, y))
            ds = [determinant([rows[i] for i in inds]) for inds in data["bases"]]
            left = [sum(v * x**i * y**j * ds[n] for n, i, j, v in terms) for terms in identities]
            right = [6*x + 5*y**4 - 20*y**2 + 3*y + 12, y**2 * (y-2) * (y-1) * (y+2)]
            if left != right:
                raise ArithmeticError("Polynomial rank identity failed")
    finite = {(-F(5*y**4 - 20*y**2 + 3*y + 12, 6), F(y)) for y in (0, 2, 1, -2)}
    if finite != {(-z.real, -z.imag) for z in sp.B_SPECTRUM}:
        raise ArithmeticError("The exceptional points differ from the four collisions")
    infinity = vertex_rows(sp.Z, reciprocal=True)
    if determinant([infinity[i] for i in data["infinityBasis"]]) != 100_000:
        raise ArithmeticError("Reciprocal chart rank check failed")
    print("Two degree-eight polynomial identities isolate all four rank failures. Rank at infinity is nine.")


if __name__ == "__main__":
    strict_patches(json.loads((ROOT / "complex-direction-cover/patches.json").read_text()),
                   json.loads((HERE / "corrections.json").read_text()))
    rank_identities(json.loads((HERE / "rank.json").read_text()))
