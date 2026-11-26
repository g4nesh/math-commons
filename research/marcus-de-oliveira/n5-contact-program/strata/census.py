#!/usr/bin/env python3
"""Census of the automatic 8-element contact sets at n = 5, by exact cutting planes.

Loop:
  1. An integer program (HiGHS) proposes an 8-set T meeting the support of every known non-minor
     facet of Q_5 and avoiding every orbit of automatic sets already recorded.
  2. Exact test (rational LPs, cddlib/GMP).  For each functional phi vanishing on span N_T (the
     squared-minor vectors vanishing on T), maximise and minimise phi.g over
     F_T = {g in L : g >= 0, g_T = 0, sum g = 1}.
       - If some optimum is nonzero, the optimal vertex is a gap vector outside span N_T; it is an
         extreme ray of C, i.e. a non-minor facet of Q_5.  It is certified exactly and its orbit
         is added to the known facets.  T is not automatic.
       - If all optima are zero, F_T lies in span N_T.  If moreover every vector of N_T has a
         private permutation, the coefficients of any g in F_T are coordinates of g, hence >= 0,
         so F_T = cone N_T and T is automatic (EXACT).  Its orbit is recorded and excluded.
  3. Stop when the integer program is infeasible: then every 8-set either lies in a known non-minor
     facet (not automatic, exactly) or is in a recorded automatic orbit.  The infeasibility claim
     rests on HiGHS (floating point branch and bound); it is not an exact certificate.

Usage:  python3 -I -B census.py ITERATIONS TIME_LIMIT_PER_ILP
Reads and extends facets_n5.json; writes census_n5.json.
"""
import json
import os
import sys
import time
from fractions import Fraction

import cdd.gmp as cdd
import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from automatic_certificate import nullspace_rows  # noqa: E402
from facet_search import canon, certify, group  # noqa: E402
from minor_polytope import setup  # noqa: E402


def exact_test(T, perms, pairs, pts, cols, Lbasis):
    """Return ('automatic', None) / ('not-automatic', zero_set of a new non-minor ray) / ('unresolved', why)."""
    N = len(perms)
    NT = [k for k in range(len(pairs)) if all(pts[s][k] == 0 for s in T)]
    phis = nullspace_rows([[cols[k][s] for s in range(N)] for k in NT] or [[0] * N], N)
    # restrict to functionals nonzero on L (others vanish identically on F_T)
    phis = [ph for ph in phis if any(sum(ph[s] * g[s] for s in range(N)) != 0 for g in Lbasis)]
    H = [[Fraction(0)] + [g[s] for g in Lbasis] for s in range(N)]                      # g >= 0
    for s in T:                                                                          # g_T = 0
        H.append([Fraction(0)] + [-g[s] for g in Lbasis])
    tot = [sum(g[s] for s in range(N)) for g in Lbasis]
    H.append([Fraction(-1)] + tot)
    H.append([Fraction(1)] + [-v for v in tot])                                         # sum g = 1
    for ph in phis:
        obj = [sum(ph[s] * g[s] for s in range(N)) for g in Lbasis]
        for sign in (1, -1):
            mat = cdd.matrix_from_array(H, rep_type=cdd.RepType.INEQUALITY)
            mat.obj_type = cdd.LPObjType.MAX
            mat.obj_func = [Fraction(0)] + [sign * v for v in obj]
            lp = cdd.linprog_from_matrix(mat)
            cdd.linprog_solve(lp)
            if lp.status != cdd.LPStatusType.OPTIMAL:
                return "unresolved", f"LP status {lp.status}"
            if Fraction(lp.obj_value) != 0:
                x = [Fraction(v) for v in lp.primal_solution]
                g = [sum(xi * b[s] for xi, b in zip(x, Lbasis)) for s in range(N)]
                return "not-automatic", sorted(s for s in range(N) if g[s] == 0)
    private = all(any(pts[s][k] == 1 and all(pts[s][j] == 0 for j in NT if j != k) for s in range(N)) for k in NT)
    if private:
        return "automatic", None
    # F_T lies in span N_T: F_T = span(N_T) cap orthant.  Enumerate its extreme rays exactly.
    from flint import fmpq_mat
    rref, rank = fmpq_mat([[pts[s][k] for s in range(N)] for k in NT]).rref()
    basis = [[Fraction(int(rref[i, s].p), int(rref[i, s].q)) for s in range(N)] for i in range(rank)]
    Hs = [[Fraction(0)] + [b[s] for b in basis] for s in range(N) if any(b[s] != 0 for b in basis)]
    gens = cdd.copy_generators(cdd.polyhedron_from_matrix(cdd.matrix_from_array(Hs, rep_type=cdd.RepType.INEQUALITY)))
    minor_zero = {frozenset(s for s in range(N) if pts[s][k] == 0) for k in NT}
    bad = []
    for row in gens.array:
        if row[0] != 0:
            continue
        g = [sum(Fraction(x) * b[s] for x, b in zip(row[1:], basis)) for s in range(N)]
        Z = frozenset(s for s in range(N) if g[s] == 0)
        if Z not in minor_zero:
            bad.append(sorted(Z))
    return ("not-automatic-rays", bad) if bad else ("automatic", None)


def main(iterations, tlimit):
    perms, pairs, pts = setup(5)
    N = len(perms)
    G = group(perms)
    cols = [[pts[s][k] for s in range(N)] for k in range(len(pairs))] + [[1] * N]
    # a basis of L as exact vectors
    from flint import fmpq_mat
    rref, rank = fmpq_mat(cols).rref()
    Lbasis = [[Fraction(int(rref[i, s].p), int(rref[i, s].q)) for s in range(N)] for i in range(rank)]
    fpath = os.path.join(HERE, "facets_n5.json")
    facets = json.load(open(fpath))
    cpath = os.path.join(HERE, "census_n5.json")
    census = json.load(open(cpath)) if os.path.exists(cpath) else {"automatic_orbits": [], "log": []}
    words = lambda S: ["".join(map(str, perms[s])) for s in sorted(S)]

    def supports():
        fam = set()
        for o in facets["orbits"]:
            if o["minor"]:
                continue
            sup = sorted(set(range(N)) - set(o["zero_set"]))
            fam |= {tuple(sorted(G[g][sup])) for g in range(G.shape[0])}
        return sorted(fam)

    fam = supports()
    for it in range(iterations):
        A = np.zeros((len(fam), N))
        for i, s in enumerate(fam):
            A[i, list(s)] = 1
        cons = [LinearConstraint(A, lb=1), LinearConstraint(np.ones((1, N)), lb=8, ub=8)]
        excl = []
        for a in census["automatic_orbits"]:
            T0 = a["T"]
            excl += sorted({tuple(sorted(G[g][T0])) for g in range(G.shape[0])})
        if excl:
            E = np.zeros((len(excl), N))
            for i, s in enumerate(excl):
                E[i, list(s)] = 1
            cons.append(LinearConstraint(E, ub=7))
        t0 = time.time()
        r = milp(np.ones(N), constraints=cons, integrality=np.ones(N), bounds=Bounds(0, 1), options={"time_limit": tlimit})
        dt = time.time() - t0
        if r.x is None:
            msg = f"iteration {it}: ILP {r.message} after {dt:.0f}s with {len(fam)} supports, {len(census['automatic_orbits'])} automatic orbits"
            print(msg, flush=True)
            census["log"].append(msg)
            census["stopped"] = r.message
            break
        T = sorted(i for i in range(N) if r.x[i] > 0.5)
        verdict, info = exact_test(T, perms, pairs, pts, cols, Lbasis)
        msg = f"iteration {it}: T={words(T)} ILP {dt:.0f}s -> {verdict}"
        if verdict == "not-automatic":
            cert = certify(info, pts, 41) is not None
            key = canon(info, G)
            facets["orbits"].append({"hits": 0, "vertices": len(info), "minor": False, "certified": cert,
                                     "zero_set": info, "found_by": "census", "key": key})
            json.dump(facets, open(fpath, "w"), indent=1)
            fam = supports()
            msg += f"; new non-minor facet with {len(info)} vertices, certified={cert}; {len(fam)} supports"
        elif verdict == "not-automatic-rays":
            seen, added = {o.get("key") for o in facets["orbits"]}, 0
            for Z in sorted(info, key=len, reverse=True):
                key = canon(Z, G)
                if key in seen:
                    continue
                seen.add(key)
                cert = certify(Z, pts, 41) is not None
                facets["orbits"].append({"hits": 0, "vertices": len(Z), "minor": False, "certified": cert,
                                         "zero_set": Z, "found_by": "census-rays", "key": key})
                added += 1
                if added == 3:
                    break
            json.dump(facets, open(fpath, "w"), indent=1)
            fam = supports()
            msg += f"; {len(info)} non-minor rays, {added} new orbits added; {len(fam)} supports"
        elif verdict == "automatic":
            census["automatic_orbits"].append({"T": T, "words": words(T)})
        else:
            census.setdefault("unresolved", []).append({"T": T, "why": info})
            census["automatic_orbits"].append({"T": T, "words": words(T), "unresolved": True})
        print(msg, flush=True)
        census["log"].append(msg)
        json.dump(census, open(cpath, "w"), indent=1)


if __name__ == "__main__":
    main(int(sys.argv[1]), float(sys.argv[2]))
