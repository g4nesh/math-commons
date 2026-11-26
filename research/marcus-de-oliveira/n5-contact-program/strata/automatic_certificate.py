#!/usr/bin/env python3
"""Exact certificate that a contact set T at n = 5 is automatic.

Let N_T be the squared-minor vectors vanishing on T, U the union of their supports and Z = S_5 \\ U.
(1) For every s in Z \\ T an exact u in L^perp with u >= 0 off T and u_s > 0 (from an exact LP over
    rationals, cddlib/GMP) shows g_s = 0 for every gap vector g in C with g_T = 0: indeed
    0 = u.g >= u_s g_s.  Hence F_T lies in L_Z = {g in L : g_Z = 0}.
(2) dim L_Z = rank N_T, so L_Z = span N_T.
(3) Every vector of N_T has a private permutation (in its support, outside all the others), so the
    coefficient of g in N_T equals a coordinate of g and is >= 0.
Together: F_T = cone(N_T), i.e. T is automatic.  Usage: python3 -I -B automatic_certificate.py PERM ...
"""
import os
import sys
from fractions import Fraction

import cdd.gmp as cdd
from flint import fmpq, fmpq_mat

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from minor_polytope import setup  # noqa: E402


def nullspace_rows(rows, ncols):
    """Rational basis of {x : R x = 0} for the matrix R given by rows."""
    R = fmpq_mat([[fmpq(x.numerator, x.denominator) if isinstance(x, Fraction) else x for x in r] for r in rows])
    rref, rank = R.rref()
    piv, r = [], 0
    for c in range(ncols):
        if r < rank and rref[r, c] != 0:
            piv.append(c); r += 1
    basis = []
    for f in (c for c in range(ncols) if c not in piv):
        y = [Fraction(0)] * ncols
        y[f] = Fraction(1)
        for i, c in enumerate(piv):
            q = rref[i, f]
            y[c] = -Fraction(int(q.p), int(q.q))
        basis.append(y)
    return basis


def main(words):
    perms, pairs, pts = setup(5)
    N = len(perms)
    idx = {"".join(map(str, p)): i for i, p in enumerate(perms)}
    T = sorted(idx[w] for w in words)
    cols = [[pts[s][k] for s in range(N)] for k in range(len(pairs))] + [[1] * N]
    NT = [k for k in range(len(pairs)) if all(pts[s][k] == 0 for s in T)]
    U = {s for k in NT for s in range(N) if pts[s][k] == 1}
    Z = [s for s in range(N) if s not in U]
    # (1) L^perp = {u : sum_s u_s v_s = 0 for every column v of L}
    perp = nullspace_rows(cols, N)
    ok1 = True
    for s in Z:
        if s in T:
            continue
        # variables: coefficients c of u = sum c_i perp_i ; constraints u_t >= 0 (t not in T), u_s >= 1
        H = []
        for t in range(N):
            if t in T:
                continue
            row = [Fraction(-1 if t == s else 0)] + [b[t] for b in perp]
            H.append(row)
        mat = cdd.matrix_from_array(H, rep_type=cdd.RepType.INEQUALITY)
        mat.obj_type = cdd.LPObjType.MAX
        mat.obj_func = [Fraction(0)] * (len(perp) + 1)
        lp = cdd.linprog_from_matrix(mat)
        cdd.linprog_solve(lp)
        if lp.status != cdd.LPStatusType.OPTIMAL:
            ok1 = False
            print("no certificate for", "".join(map(str, perms[s])))
            break
        x = [Fraction(v) for v in lp.primal_solution]
        u = [sum(xi * b[t] for xi, b in zip(x, perp)) for t in range(N)]
        assert all(sum(u[t] * c[t] for t in range(N)) == 0 for c in cols)
        assert all(u[t] >= 0 for t in range(N) if t not in T) and u[s] >= 1
    # (2)
    dimL = fmpq_mat(cols).rank()
    dimLZ = dimL - fmpq_mat([[c[s] for s in Z] for c in cols]).rank()
    rankNT = fmpq_mat([cols[k] for k in NT]).rank()
    # (3)
    private = all(any(pts[s][k] == 1 and all(pts[s][j] == 0 for j in NT if j != k) for s in range(N)) for k in NT)
    # (2') if dim L_Z = rank N_T + extra, test exactly whether any g >= 0 in L_Z leaves span N_T
    extra_ok = None
    if ok1 and private and dimLZ > rankNT:
        basisLZ = nullspace_rows([[cols[i][s] for i in range(len(cols))] for s in Z], len(cols))   # coefficients a with (sum a_i col_i)_Z = 0
        gens = [[sum(a[i] * cols[i][s] for i in range(len(cols))) for s in range(N)] for a in basisLZ]
        # the quotient of L_Z by span N_T: functionals phi vanishing on N_T but not on L_Z
        Mq = nullspace_rows([[cols[k][s] for s in range(N)] for k in NT], N)   # phi with phi.n_k = 0
        phis = [ph for ph in Mq if any(sum(ph[s] * g[s] for s in range(N)) != 0 for g in gens)]
        extra_ok = True
        for ph in phis[:1]:
            for sign in (1, -1):
                # variables: coefficients x of g = sum x_j gens_j; constraints g >= 0, sum g = 1; objective sign*phi.g
                H = [[Fraction(0)] + [g[s] for g in gens] for s in range(N)]
                tot = [sum(g[s] for s in range(N)) for g in gens]
                H.append([Fraction(-1)] + tot); H.append([Fraction(1)] + [-v for v in tot])
                mat = cdd.matrix_from_array(H, rep_type=cdd.RepType.INEQUALITY)
                mat.obj_type = cdd.LPObjType.MAX
                mat.obj_func = [Fraction(0)] + [sign * sum(ph[s] * g[s] for s in range(N)) for g in gens]
                lp = cdd.linprog_from_matrix(mat); cdd.linprog_solve(lp)
                print("  extra direction, sign", sign, "max", lp.obj_value, lp.status)
                if lp.status != cdd.LPStatusType.OPTIMAL or Fraction(lp.obj_value) != 0:
                    extra_ok = False
        if len(phis) == 0:
            extra_ok = True
    print(f"T={words}: |N_T|={len(NT)}, |Z|={len(Z)}, certificates for Z\\T: {ok1}, "
          f"dim L_Z={dimLZ}, rank N_T={rankNT}, private: {private}")
    auto = ok1 and private and (dimLZ == rankNT or extra_ok)
    print("AUTOMATIC (exact)" if auto else ("NOT automatic: a nonnegative gap vector leaves span N_T" if extra_ok is False else "not certified"))


if __name__ == "__main__":
    main(sys.argv[1:])
