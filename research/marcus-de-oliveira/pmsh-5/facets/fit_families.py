"""Closed forms of the multiplier families of the hard facets (see multipliers.py).

For each orbit of decompose.json, on a grid of Moebius classes (0, 1, s, t, inf), 1 < s < t, this finds
where the vertex set Z (from (s, t) = (3, 9)) is still the vertex set of a facet with the same Q5 support,
and fits every weight (divided by the first non-minor weight) by a monomial
    C * s^e1 * t^e2 * (s - 1)^e3 * (t - 1)^e4 * (t - s)^e5
with integer exponents (accepted below relative error 1e-4: the grid reaches ill-conditioned spectra; the
formulas are then checked exactly by exact_families.py). Weights that are not monomials are fitted by a
rational function of lowest total degree <= 4 (ratfit). Reports the grid region and, per term, the fit.

usage: fit_families.py ORBIT[,ORBIT...]
"""
import json
import os
import sys
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from multipliers import facet_gap, xvec  # noqa: E402
from decompose import gap_from_zero_set  # noqa: E402

FORMS = ["s", "t", "(s-1)", "(t-1)", "(t-s)"]


def forms(s, t):
    return np.array([s, t, s - 1, t - 1, t - s])


def setup(k):
    rec = json.load(open(os.path.join(HERE, "decompose.json")))["decompositions"][k]
    gap0 = np.array(json.load(open(os.path.join(HERE, "orbits.json")))["orbits"][k]["gap"])
    cols, names = [], []
    for t in rec["terms"]:
        if t["kind"] == "minor":
            cols.append(xvec(*t["pair"]))
            names.append("x" + "".join(map(str, t["pair"][0])) + ">" + "".join(map(str, t["pair"][1])))
        else:
            cols.append(gap_from_zero_set(t["zero_set"]).astype(float))
            names.append(t["kind"])
    first = next(i for i, t in enumerate(rec["terms"]) if t["kind"] != "minor")
    return np.flatnonzero(gap0 == 0), np.array(cols).T, names, first


def weights(Z, A, first, s, t):
    c = np.array([0, 0.5, s / (s + 1), t / (t + 1), 1.0])
    f = facet_gap(c, Z)
    if f is None:
        return None
    w, *_ = np.linalg.lstsq(A, f, rcond=None)
    if np.abs(A @ w - f).max() > 1e-9 or w.min() < -1e-9:
        return None  # same vertex set, different face of the gap cone
    return w / w[first]


def monomial(points, vals):
    L = np.array([np.r_[1.0, np.log(np.abs(forms(s, t)))] for s, t in points])
    coef, *_ = np.linalg.lstsq(L, np.log(vals), rcond=None)
    e = np.round(coef[1:]).astype(int)
    C = np.exp(np.median(np.log(vals) - L[:, 1:] @ e))
    pred = C * np.exp(L[:, 1:] @ e)
    err = np.abs(pred / vals - 1).max()
    Cf = Fraction(C).limit_denominator(1000)
    num = [f"{FORMS[i]}^{e[i]}" if e[i] != 1 else FORMS[i] for i in range(5) if e[i] > 0]
    den = [f"{FORMS[i]}^{-e[i]}" if e[i] != -1 else FORMS[i] for i in range(5) if e[i] < 0]
    txt = (f"{Cf}*" if Cf != 1 else "") + ("*".join(num) or "1") + ("/(" + "*".join(den) + ")" if den else "")
    return txt, err


def ratfit(points, vals, maxdeg=4):
    """Rational function N/D in (s, t) of lowest total degree fitting the values (numerically, then
    rationalized); returns a factored sympy expression and the max relative error, or None."""
    import sympy as sp
    s_, t_ = sp.symbols("s t")
    good = [(p, v) for p, v in zip(points, vals) if 1.2 < p[0] < 12 and p[1] < 60 and p[1] > 1.1 * p[0]]
    P = np.array([p for p, _ in good])
    v = np.array([x for _, x in good])
    for deg in range(1, maxdeg + 1):
        mons = [(i, j) for i in range(deg + 1) for j in range(deg + 1 - i)]
        M = np.array([[s ** i * t ** j for i, j in mons] for s, t in P])
        A = np.c_[M, -v[:, None] * M]
        A = A / np.abs(A).max(1, keepdims=True)
        _, sv, vt = np.linalg.svd(A)
        if sv[-1] > 1e-9 * sv[0]:
            continue
        x = vt[-1] / vt[-1][np.argmax(np.abs(vt[-1]))]
        q = [sp.Rational(str(Fraction(float(a)).limit_denominator(10000))) for a in x]
        N = sum(q[k] * s_ ** i * t_ ** j for k, (i, j) in enumerate(mons))
        D = sum(q[len(mons) + k] * s_ ** i * t_ ** j for k, (i, j) in enumerate(mons))
        expr = sp.factor(sp.cancel(N / D))
        f = sp.lambdify((s_, t_), expr)
        err = max(abs(f(s, t) / val - 1) for (s, t), val in zip(points, vals))
        return expr, err
    return None


def main():
    grid = [(s, t) for s in np.geomspace(1.05, 40, 23) for t in np.geomspace(1.1, 400, 31) if t > s * 1.02]
    for k in (int(x) for x in sys.argv[1].split(",")):
        Z, A, names, first = setup(k)
        pts, W = [], []
        for s, t in grid:
            w = weights(Z, A, first, s, t)
            if w is not None:
                pts.append((s, t))
                W.append(w)
        if not pts:
            print(f"orbit {k}: the support found at (3, 9) gives the facet at no grid point")
            continue
        W = np.array(W)
        S = np.array([p[0] for p in pts])
        T = np.array([p[1] for p in pts])
        print(f"orbit {k}: facet with the same support at {len(pts)} of {len(grid)} grid points; "
              f"s in [{S.min():.3g}, {S.max():.3g}], t in [{T.min():.3g}, {T.max():.3g}], t/s in [{(T / S).min():.3g}, {(T / S).max():.3g}]")
        for j, name in enumerate(names):
            if j == first:
                continue
            txt, err = monomial(pts, W[:, j])
            if err > 1e-4:
                r = ratfit(pts, W[:, j])
                txt, err = (str(r[0]), r[1]) if r else ("no fit up to degree 4", np.nan)
            print(f"   {name}: {txt}  (max rel. err {err:.1e})")


if __name__ == "__main__":
    main()
