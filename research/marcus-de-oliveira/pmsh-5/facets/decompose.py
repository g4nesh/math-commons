"""Write each hard facet of P5(c) as a nonnegative combination of facets of the minor polytope Q5.

Q5 = conv{x(sigma)} with x_{K,J}(sigma) = 1[sigma(K) = J] (squared minors of permutation matrices); see
../../n5-contact-program/strata/README.md (PR 115), whose facets_n5.json lists two orbits of
squared-minor facets x_{K,J} >= 0 and seven certified orbits A-G of non-minor facets under the group of
order 28,800 (sigma -> a o sigma o b and sigma -> a o sigma^-1 o b).

A facet of P(c) pulls back to an affine function on aff(Q5) that is >= 0 on Q5. Its gap vector
f(sigma) = h - l.v_sigma (orbits.json) lies in the gap cone, whose extreme rays are the facets of Q5
(Lemma of the strata README). So f = sum mu_F g_F + sum mu_{K,J} x_{K,J} with mu >= 0, over the facets
of Q5 vanishing on the vertex set Z of the P(c)-facet. "Hard" means some non-minor g_F is needed.

For every orbit representative this solves that LP over the known facets (exact gap vectors of the
seven non-minor orbits, 57,780 facets), minimizing the total weight on non-minor facets, and reports the
non-minor facets used (orbit label and multiplier) and the number of squared-minor terms. If the LP is
infeasible, the known list misses a facet of Q5 through Z. Then f is decomposed exactly into extreme rays
of the gap cone by a Caratheodory peel (peel()); every non-minor ray it uses is certified exactly as a facet
of Q5 (its zero set has affine rank 41 and the vanishing affine function is one-signed), its orbit under the
group of order 28,800 is added, and the LP is solved again. New orbits are labelled N1, N2, ... in the
order (number of vertices, canonical vertex set) at the end; the canonical vertex set is the
lexicographically smallest sorted image.

usage: decompose.py
"""
import itertools
import json
import os
import sys

import numpy as np
import sympy as sp
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
STRATA = os.path.join(HERE, "..", "..", "n5-contact-program", "strata")
n = 5
PERMS = list(itertools.permutations(range(n)))
PIDX = {p: i for i, p in enumerate(PERMS)}
PAIRS = [(K, J) for k in (1, 2) for K in itertools.combinations(range(n), k) for J in itertools.combinations(range(n), k)]
X = np.array([[int(tuple(sorted(p[i] for i in K)) == J) for K, J in PAIRS] for p in PERMS])  # 120 x 125
LABELS = "ABCDEFG"  # order of the non-minor orbits in facets_n5.json, as named in the strata README


def group():
    inv = [PIDX[tuple(sorted(range(n), key=lambda i: s[i]))] for s in PERMS]
    rows = []
    for a in PERMS:
        for b in PERMS:
            row = [PIDX[tuple(a[s[b[i]]] for i in range(n))] for s in PERMS]
            rows.append(row)
            rows.append([row[inv[k]] for k in range(len(PERMS))])
    return np.array(rows, dtype=np.int64)


_B = sp.Matrix(np.c_[X, np.ones(len(PERMS), int)].tolist())
BASIS = _B[:, list(_B.rref()[1])]  # 42 independent columns: affine functions on aff(Q5), as values on vertices


def gap_from_zero_set(Z):
    """Exact gap vector (primitive, nonnegative integers) of the facet of Q5 with vertex set Z."""
    B = BASIS
    ns = B[sorted(Z), :].nullspace()
    assert len(ns) == 1, len(ns)
    v = B * ns[0]
    den = sp.ilcm(*[sp.Rational(x).q for x in v])
    w = [int(x * den) for x in v]
    gcd = sp.igcd(*w)
    w = [x // gcd for x in w]
    if min(w) < 0:
        w = [-x for x in w]
    assert min(w) == 0 and sum(1 for x in w if x == 0) == len(Z)
    return np.array(w, dtype=np.int64)


def canon(Z, G):
    """Canonical form of a vertex set under G: the lexicographically smallest sorted image."""
    return min(tuple(sorted(row[list(Z)])) for row in G)


def peel(f, rng):
    """Exact Caratheodory peel of a gap vector f (integers): repeatedly subtract squared-minor rays as far
    as possible, then walk from the remainder to an extreme ray of its face of the gap cone by exact random
    pushes inside the face, and subtract that. Returns the zero sets of the non-minor extreme rays used."""
    from fractions import Fraction as Fr
    B = BASIS
    g = [Fr(int(x)) for x in f]
    found = []
    while any(g):
        for k in range(X.shape[1]):
            if all(g[s] == 0 for s in range(len(PERMS)) if X[s, k]) or not any(X[:, k]):
                continue
            if all(X[s, k] == 0 for s in range(len(PERMS)) if g[s] == 0):
                t = min(g[s] for s in range(len(PERMS)) if X[s, k])
                g = [g[s] - t * int(X[s, k]) for s in range(len(PERMS))]
        if not any(g):
            break
        r = list(g)
        while True:
            Zr = [s for s in range(len(PERMS)) if r[s] == 0]
            ns = B[Zr, :].nullspace() if Zr else [sp.eye(B.shape[1])[:, j] for j in range(B.shape[1])]
            ns = [v for v in ns if any(x != 0 for x in B * v)]
            if len(ns) <= 1:
                break  # the face of r is a ray
            y = sum((int(rng.integers(-9, 10)) * v for v in ns), sp.zeros(B.shape[1], 1))
            d = [Fr(int(sp.Rational(x).p), int(sp.Rational(x).q)) for x in B * y]
            if all(x >= 0 for x in d):
                d = [-x for x in d]
            if all(x == 0 for x in d):
                continue
            t = min(r[s] / -d[s] for s in range(len(PERMS)) if d[s] < 0)
            r = [r[s] + t * d[s] for s in range(len(PERMS))]
        Zr = tuple(s for s in range(len(PERMS)) if r[s] == 0)
        t = min(g[s] / r[s] for s in range(len(PERMS)) if r[s] > 0)
        g = [g[s] - t * r[s] for s in range(len(PERMS))]
        assert min(g) >= 0
        if not any(all(X[s, k] == 0 for s in Zr) and sum(X[:, k]) == len(PERMS) - len(Zr) for k in range(X.shape[1])):
            found.append(Zr)
    return found


def main():
    fac = json.load(open(os.path.join(STRATA, "facets_n5.json")))
    nonminor = [o for o in fac["orbits"] if not o["minor"]]
    G = group()
    gaps, labels = [], []
    for lab, o in zip(LABELS, nonminor):
        g0 = gap_from_zero_set(o["zero_set"])
        seen = {}
        for row in G:
            img = np.empty(len(PERMS), np.int64)
            img[row] = g0  # value at the image of sigma is g0(sigma)
            seen.setdefault(img.tobytes(), img)
        gaps.extend(seen.values())
        labels.extend([lab] * len(seen))
        print(f"orbit {lab}: {len(seen)} facets, gap values {sorted(set(g0.tolist()))}", flush=True)
    GA = np.array(gaps)
    labels = np.array(labels)
    orbits = json.load(open(os.path.join(HERE, "orbits.json")))["orbits"]
    known = {g.tobytes() for g in GA}
    minor_sets = {tuple(np.flatnonzero(X[:, k] == 0).tolist()) for k in range(X.shape[1])}
    new = {}  # canonical zero set -> (label, exact gap of the canonical representative, orbit size)
    rng = np.random.default_rng(0)
    out = []
    for o in orbits:
        f = np.array(o["gap"], dtype=np.int64)
        Z = f == 0
        for attempt in range(2):  # one round of facet discovery
            ex = np.flatnonzero(~(GA[:, Z] != 0).any(1))
            mi = np.flatnonzero(~(X[Z] != 0).any(0))
            A = np.c_[GA[ex].T, X[:, mi]].astype(float)
            cost = np.r_[np.ones(len(ex)), np.zeros(len(mi))]
            r = linprog(cost, A_eq=A, b_eq=f.astype(float), bounds=(0, None), method="highs")
            if r.status == 0 or attempt == 1:
                break
            for Zr in peel(f, rng):
                gr = gap_from_zero_set(Zr)
                if gr.tobytes() in known:
                    continue
                key = tuple(int(x) for x in canon(Zr, G))
                if key not in new:
                    imgs = {}
                    for row in G:
                        img = np.empty(len(PERMS), np.int64)
                        img[row] = gr
                        imgs.setdefault(img.tobytes(), img)
                    new[key] = (None, gap_from_zero_set(key), len(imgs))
                    GA = np.r_[GA, np.array(list(imgs.values()))]
                    labels = np.r_[labels, [str(key)] * len(imgs)]
                    known |= set(imgs)
                    print(f"  new non-minor facet orbit: {len(key)} vertices, orbit size {len(imgs)}, "
                          f"gap values {sorted(set(gr.tolist()))}", flush=True)
        rec = dict(orbit=o["orbit"], status=o["status"], vertices=int(Z.sum()), candidates_nonminor=len(ex),
                   candidates_minor=len(mi), feasible=r.status == 0)
        if r.status == 0:
            sup = [k for k in range(A.shape[1]) if r.x[k] > 1e-9]
            As = sp.Matrix(A[:, sup].astype(int).tolist())
            w = (As.T * As).LUsolve(As.T * sp.Matrix(f.tolist()))  # exact; the support columns are independent
            assert As * w == sp.Matrix(f.tolist()) and min(w) > 0, "LP support does not give an exact decomposition"
            terms = []
            for k, wk in zip(sup, w):
                if k < len(ex):
                    terms.append(dict(kind=str(labels[ex[k]]), weight=str(wk), zero_set=np.flatnonzero(GA[ex[k]] == 0).tolist()))
                else:
                    K, J = PAIRS[mi[k - len(ex)]]
                    terms.append(dict(kind="minor", weight=str(wk), pair=[list(K), list(J)]))
            rec.update(terms=terms)
        out.append(rec)
    # name the new orbits N1, N2, ... by (number of vertices, canonical vertex set)
    names = {str(k): f"N{i + 1}" for i, k in enumerate(sorted(new, key=lambda k: (len(k), k)))}
    for rec in out:
        for t in rec.get("terms", []):
            t["kind"] = names.get(t["kind"], t["kind"])
        summary = {k: v for k, v in rec.items() if k != "terms"}
        summary["nonminor"] = [(t["kind"], t["weight"]) for t in rec.get("terms", []) if t["kind"] != "minor"]
        summary["minor_terms"] = sum(1 for t in rec.get("terms", []) if t["kind"] == "minor")
        print(json.dumps(summary), flush=True)
    for k in sorted(new, key=lambda k: (len(k), k)):
        print(f"{names[str(k)]}: {len(k)} vertices, orbit size {new[k][2]}, gap values {sorted(set(new[k][1].tolist()))}, "
              f"canonical vertex set {list(k)}")
    json.dump(dict(new_orbits=[dict(label=names[str(k)], zero_set=[int(x) for x in k], gap=new[k][1].tolist(), orbit_size=new[k][2])
                               for k in sorted(new, key=lambda k: (len(k), k))], decompositions=out),
              open(os.path.join(HERE, "decompose.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
