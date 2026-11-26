"""Replay nonnegative radial weights uniformly over closed angular intervals."""
from fractions import Fraction as F
from functools import lru_cache
from math import comb
from pathlib import Path
import importlib.util
import json

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("real_atlas", HERE.parent / "two-free-real-atlas" / "check.py")
atlas = importlib.util.module_from_spec(spec)
spec.loader.exec_module(atlas)


@lru_cache(None)
def gaps():
    target, vertices = atlas.fixed_data()
    return tuple(atlas.poly_add(v, atlas.poly_scale(-1, target)) for v in vertices)


@lru_cache(None)
def square_product(centers):
    h = [F(1)]
    for r in centers:
        out = [F(0)]*(len(h)+2)
        for i, coeff in enumerate(h):
            out[i] += r*r*coeff
            out[i+1] -= 2*r*coeff
            out[i+2] += coeff
        h = out
    return h


@lru_cache(None)
def radial_gap(chart, slope, centers, permutation):
    """Coefficients in (free x, free y, real radius) of the weighted gap."""
    omega = (F(1), slope) if chart == "h" else (slope, F(1))
    result = {}
    for (x, middle, y), coeff in gaps()[permutation].items():
        if middle:
            coeff = atlas.mul(coeff, omega)
        for r, factor in enumerate(square_product(centers)):
            key = (x, y, middle+r)
            result[key] = atlas.add(result.get(key, atlas.ZERO), atlas.scale(factor, coeff))
    return result


def column(chart, lo, hi, degree, control, centers, perm):
    entries = {("mass", control): F(1)}
    for j, k, slope in ((control, degree+1-control, lo), (control+1, control+1, hi)):
        for monomial, value in radial_gap(chart, slope, centers, perm).items():
            for coord in range(2):
                if value[coord]:
                    entries[(j, *monomial, coord)] = k*value[coord]
    return entries


def verify_patch(patch):
    chart = patch["chart"]
    lo, hi = map(F, patch["interval"])
    degree = patch["degree"]
    assert chart in ("h", "v") and -1 <= lo < hi <= 1 and 1 <= degree <= 10
    totals = {}
    roots = {j: None for j in range(degree+1)}
    used = set()
    for control, word, origin, weight in patch["terms"]:
        w = F(weight)
        c = tuple(map(F, origin))
        p = tuple(int(i)-1 for i in word)
        assert w > 0 and 0 <= control <= degree and tuple(sorted(c)) == c
        assert (control, word, c) not in used
        used.add((control, word, c))
        idx = atlas.PERMS.index(p)
        roots[control] = set(c) if roots[control] is None else roots[control].intersection(c)
        for key, value in column(chart, lo, hi, degree, control, c, idx).items():
            totals[key] = totals.get(key, F(0))+w*value
    assert {k: v for k, v in totals.items() if v} == {("mass", j): F(1) for j in roots}
    assert all(r == set() for r in roots.values())
    return chart, lo, hi


def weights_at(real, imag, patches=None):
    """Return column-to-row pairing weights for a Gaussian rational middle eigenvalue."""
    if any(isinstance(v, bool) or not isinstance(v, (int, F)) for v in (real, imag)):
        raise TypeError("Use integers or fractions for the real and imaginary parts")
    a, b = F(real), F(imag)
    if a == b == 0:
        choices = [("h", F(0), F(0))]
    else:
        choices = []
        if a:
            choices.append(("h", b/a, a))
        if b:
            choices.append(("v", a/b, b))
    if patches is None:
        patches = json.loads((HERE / "patches.json").read_text())
    for chart, slope, radius in choices:
        for patch in patches:
            lo, hi = map(F, patch["interval"])
            if patch["chart"] != chart or not lo <= slope <= hi:
                continue
            verify_patch(patch)
            alpha = (slope-lo)/(hi-lo)
            degree = patch["degree"]
            result = {}
            for j, word, origin, mass in patch["terms"]:
                w = F(mass)*comb(degree, j)*alpha**j*(1-alpha)**(degree-j)
                for center in origin:
                    w *= (radius-F(center))**2
                result[word] = result.get(word, F(0))+w
            total = sum(result.values())
            assert total > 0
            return {word: w/total for word, w in result.items() if w}
    raise ValueError("No supplied certificate covers this direction")


def main():
    data = json.loads((HERE / "patches.json").read_text())
    verified = [verify_patch(p) for p in data]
    for chart in ("h", "v"):
        ranges = sorted((a, b) for c, a, b in verified if c == chart)
        assert ranges == [(F(-1), F(-1, 2)), (F(-1, 2), F(0)),
                          (F(0), F(1, 2)), (F(1, 2), F(1))]
        print(chart, " ".join(f"[{a},{b}]" for a, b in ranges))
    print(f"Verified all mixed-polynomial coefficients for {len(data)} patches covering every complex direction.")


if __name__ == "__main__":
    main()
