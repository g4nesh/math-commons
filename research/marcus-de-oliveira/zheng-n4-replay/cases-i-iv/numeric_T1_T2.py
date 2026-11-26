#!/usr/bin/env python3
"""NUMERIC falsification attempts for the intermediate claims in the proofs of T_1 and T_2.

Floating point only.  Nothing here is evidence for the proofs; a violation would be evidence against them.
Usage:  python numeric_T1_T2.py T1    or    python numeric_T1_T2.py T2
Exit code 0 iff no violation beyond the stated tolerance was found.  Packages: numpy, scipy.

Notation: spectra A, C in C^4 and lam in C; F_sigma = lam * prod_i (A_i - C_sigma(i)), permutations in
lexicographic order (index 0 = identity).  "Equal real parts on T" means Re F_sigma = Re F_tau for sigma, tau in T.

How points are produced.  The unknown is z = (A_0..A_3, C_0..C_3, lam) in C^9 = R^18.  The equations
Re(F_sigma - F_tau0) = 0 (sigma in T) and any extra complex-linear side conditions are solved by Gauss-Newton
with minimum-norm steps from Gaussian random starts.  No structure from the proofs is used to produce points,
so a component of the solution set on which a claim fails would be hit if it attracts random starts.
For claims that need the inequalities Re F_sigma >= min as well, a penalised least-squares solve is used.
"""
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import functools
import itertools
import numpy as np
from scipy.optimize import least_squares, linprog

import t12_common as c

print = functools.partial(print, flush=True)
PERMS = np.array(c.PERMS)
T1 = sorted(c.PAPER_T[1])
T2 = sorted(c.PAPER_T[2])
MINORS = np.array(sorted(c.ALL_MINOR_PROFILES), dtype=float).T          # 24 x 34


def F_and_jac(z):
    """F (24,) complex and J (24, 9) complex with J[s, k] = dF_s / dz_k, z = (A, C, lam)."""
    A, C, lam = z[0:4], z[4:8], z[8]
    fac = A[None, :] - C[PERMS]                     # (24, 4): A_i - C_sigma(i)
    prod = np.prod(fac, axis=1)
    F = lam * prod
    J = np.zeros((24, 9), dtype=complex)
    for i in range(4):
        others = np.prod(np.delete(fac, i, axis=1), axis=1) * lam
        J[:, i] += others
        for s in range(24):
            J[s, 4 + PERMS[s, i]] -= others[s]
    J[:, 8] = prod
    return F, J


def solve_equalities(T, rng, side=(), tol=1e-13, maxit=200):
    """Gauss-Newton, minimum-norm steps.  side = list of (coefficient vector in C^9, constant): complex-linear
    conditions coeff . z = const.  Returns z or None."""
    x = rng.normal(size=18)
    for _ in range(maxit):
        z = x[:9] + 1j * x[9:]
        F, J = F_and_jac(z)
        res, rows = [], []
        for s in T[1:]:
            res.append((F[s] - F[T[0]]).real)
            dj = J[s] - J[T[0]]
            rows.append(np.concatenate([dj.real, -dj.imag]))
        for coef, const in side:
            v = np.dot(coef, z) - const
            res += [v.real, v.imag]
            rows.append(np.concatenate([coef.real, -coef.imag]))
            rows.append(np.concatenate([coef.imag, coef.real]))
        res = np.array(res)
        scale = max(1.0, np.max(np.abs(F)))
        if np.max(np.abs(res)) < tol * scale:
            return z
        step = np.linalg.lstsq(np.array(rows), -res, rcond=None)[0]
        nrm = np.linalg.norm(step)
        if nrm > 2.0:
            step *= 2.0 / nrm
        x = x + step
        if not np.all(np.isfinite(x)) or np.max(np.abs(x)) > 1e3:
            return None
    return None


def nondegenerate(z, T, dmin):
    """All within-spectrum distances, |lam| and all factors of the products in T are at least dmin, and the
    point is of moderate size."""
    A, C, lam = z[0:4], z[4:8], z[8]
    if max(np.max(np.abs(A)), np.max(np.abs(C)), abs(lam)) > 30 or abs(lam) < dmin:
        return False
    for i, j in itertools.combinations(range(4), 2):
        if abs(A[i] - A[j]) < dmin or abs(C[i] - C[j]) < dmin:
            return False
    for s in T:
        if np.min(np.abs(A - C[PERMS[s]])) < dmin:
            return False
    return True


def cross_ratio(a, b, cc, d):
    return (a - cc) * (b - d) / ((a - d) * (b - cc))


def rel_im(w):
    return abs(w.imag) / max(abs(w), 1e-300)


def cone_residual(g, gens=None):
    """Relative sup-norm distance of g (24-vector) from the cone spanned by the columns of `gens` (default: the
    34 minor profiles): min over a >= 0 of max_i |(gens a - g)_i|, divided by max_i |g_i|.  Solved as a linear
    programme (HiGHS).  (scipy 1.13's nnls raised 'maximum number of iterations' on some inputs, so it is not used.)"""
    M = MINORS if gens is None else gens
    n = M.shape[1]
    cost = np.zeros(n + 1)
    cost[n] = 1.0
    ones = np.ones((24, 1))
    A_ub = np.vstack([np.hstack([M, -ones]), np.hstack([-M, -ones])])
    b_ub = np.concatenate([g, -g])
    res = linprog(cost, A_ub=A_ub, b_ub=b_ub, bounds=[(0, None)] * (n + 1), method="highs")
    if res.status != 0:
        return float("inf")
    return float(res.x[n]) / max(float(np.max(np.abs(g))), 1e-300)


def control(rp):
    """The cone test must be able to say no: f_17 and f_17 + 0.5 f_9 are not in the minor cone, f_17 + f_9 is."""
    f17 = np.array(c.F_PROF[17], dtype=float)
    f9 = np.array(c.F_PROF[9], dtype=float)
    d_out1, d_out2, d_in = cone_residual(f17), cone_residual(f17 + 0.5 * f9), cone_residual(f17 + f9)
    ok = d_out1 > 1e-2 and d_out2 > 1e-2 and d_in < 1e-7
    if not ok:
        rp.bad += 1
    print("[%s] control of the cone test: distance of f_17 = %.4f, of f_17 + f_9/2 = %.4f (both outside the minor "
          "cone), of f_17 + f_9 = %.1e (inside)" % ("ok" if ok else "VIOLATION", d_out1, d_out2, d_in))
    rp.info("the linear programme has feasibility tolerance about 1e-7, so a reported distance 0.000e+00 means "
            "'below about 1e-7 of the largest gap'")


def solve_contact(T, rng, side=(), margin=0.0, tol=1e-11):
    """Penalised least squares for: Re F equal on T, Re F_sigma >= that value + margin elsewhere, side conditions,
    normalised so that the mean of the gaps off T is 1 (this excludes the zero gap vector).
    One start per call; returns z or None."""
    others = [s for s in range(24) if s not in T]

    def res(x):
        z = x[:9] + 1j * x[9:]
        F, _ = F_and_jac(z)
        r = [(F[s] - F[T[0]]).real for s in T[1:]]
        r += [min(0.0, (F[s] - F[T[0]]).real - margin) for s in others]
        for coef, const in side:
            v = np.dot(coef, z) - const
            r += [v.real, v.imag]
        r.append(np.mean([(F[s] - F[T[0]]).real for s in others]) - 1.0)   # normalisation: mean gap off T = 1
        return np.array(r)

    x0 = rng.normal(size=18)
    try:
        sol = least_squares(res, x0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=600)
    except Exception:
        return None
    if np.max(np.abs(res(sol.x))) < tol:
        return sol.x[:9] + 1j * sol.x[9:]
    return None


def unit(k):
    e = np.zeros(9, dtype=complex)
    e[k] = 1.0
    return e


class Report:
    def __init__(self, title):
        self.bad = 0
        print("=" * 100)
        print(title)
        print("=" * 100)

    def line(self, name, n, worst, tol, note=""):
        ok = worst <= tol
        if not ok:
            self.bad += 1
        print("[%s] %s :: samples %d, worst value %.3e, tolerance %.1e%s"
              % ("no violation" if ok else "VIOLATION", name, n, worst, tol, (" :: " + note) if note else ""))

    def info(self, text):
        print("       " + text)

    def finish(self):
        print()
        print("FINAL: %s" % ("NO VIOLATION FOUND (numeric evidence only)" if self.bad == 0
                             else "VIOLATIONS FOUND: %d" % self.bad))
        return 0 if self.bad == 0 else 1


# ----------------------------------------------------------------------------------------------------
def run_T1():
    rp = Report("NUMERIC falsification attempts, case T_1 = {0,1,2,3,6,7}  (floating point; seeds printed)")
    control(rp)
    six = lambda z: [z[1], z[2], z[3], z[5], z[6], z[7]]            # A1 A2 A3 C1 C2 C3

    def claims(z):
        A, C = z[0:4], z[4:8]
        F, _ = F_and_jac(z)
        sc = np.max(np.abs(F[T1]))
        e_first = rel_im(cross_ratio(A[2], A[3], C[2], C[3]))
        # all quadruples of pairwise distinct points among the six (on the collision A1 = C3 two of them coincide)
        e_six = max([rel_im(cross_ratio(*q)) for q in itertools.combinations(six(z), 4)
                     if min(abs(q[i] - q[j]) for i, j in itertools.combinations(range(4), 2)) > 1e-9] + [0.0])
        e_re0 = abs(F[0].real) / sc
        e_45 = max(abs((F[4] - F[0]).real), abs((F[5] - F[0]).real)) / max(sc, abs(F[4]), abs(F[5]))
        return e_first, e_six, e_re0, e_45

    # E1: equalities only, generic
    seed, want, dmin = 101, 2000, 0.05
    rng = np.random.default_rng(seed)
    got, tries, worst = 0, 0, np.zeros(4)
    while got < want and tries < 200000:
        tries += 1
        z = solve_equalities(T1, rng)
        if z is None or not nondegenerate(z, T1, dmin):
            continue
        got += 1
        worst = np.maximum(worst, claims(z))
    rp.info("E1: solutions of 'Re F equal on T_1' only (no inequalities), all distances and factors >= %.2f; "
            "seed %d, %d starts" % (dmin, seed, tries))
    rp.line("E1a first cross-ratio (A2,A3;C2,C3) is real [relative imaginary part]", got, worst[0], 1e-7)
    rp.line("E1b A1,A2,A3,C1,C2,C3 concyclic [max relative imaginary part over all quadruples of distinct points]", got,
            worst[1], 1e-7)
    rp.line("E1c Re F_0 = 0 [relative]", got, worst[2], 1e-7)
    rp.line("E1d Re F_4 = Re F_5 = Re F_0 [relative]", got, worst[3], 1e-7)

    # E2: the exceptional collision A1 = C3
    seed, want = 102, 500
    rng = np.random.default_rng(seed)
    side = [(unit(1) - unit(7), 0.0)]
    got, tries, worst = 0, 0, np.zeros(4)
    while got < want and tries < 100000:
        tries += 1
        z = solve_equalities(T1, rng, side=side)
        if z is None or not nondegenerate(z, T1, dmin):
            continue
        got += 1
        worst = np.maximum(worst, claims(z))
    rp.info("E2: same, on the collision A1 = C3 that the proof treats separately; seed %d, %d starts" % (seed, tries))
    rp.line("E2a first cross-ratio real", got, worst[0], 1e-7)
    rp.line("E2b six points concyclic", got, worst[1], 1e-7)
    rp.line("E2c Re F_0 = 0", got, worst[2], 1e-7)
    rp.line("E2d Re F_4 = Re F_5 = Re F_0", got, worst[3], 1e-7)

    # E3: approaching the excluded loci
    for label, coef, k0 in (("A3 - A2 = eps (towards a repeated eigenvalue)", unit(3) - unit(2), 103),
                            ("A1 - A2 = eps", unit(1) - unit(2), 104),
                            ("C0 - C1 = eps", unit(4) - unit(5), 105),
                            ("A0 - C0 = eps (towards a zero contact product)", unit(0) - unit(4), 106),
                            ("A1 - C3 = eps (towards the exceptional collision)", unit(1) - unit(7), 107)):
        for eps in (1e-2, 1e-4):
            rng = np.random.default_rng(k0)
            got, tries, worst = 0, 0, np.zeros(4)
            while got < 200 and tries < 100000:
                tries += 1
                z = solve_equalities(T1, rng, side=[(coef, eps)], tol=1e-15)
                if z is None:
                    continue
                # nondegenerate apart from the one small quantity
                zz = z.copy()
                A, C = z[0:4], z[4:8]
                small = sum(1 for i, j in itertools.combinations(range(4), 2)
                            if abs(A[i] - A[j]) < 0.05) + sum(1 for i, j in itertools.combinations(range(4), 2)
                                                             if abs(C[i] - C[j]) < 0.05)
                smallx = sum(1 for s in T1 for i in range(4) if abs(A[i] - C[PERMS[s][i]]) < 0.05)
                if "eps (towards a zero" in label:
                    ok = small == 0 and smallx == 4          # the factor A0 - C0 occurs in four products of T_1
                elif "exceptional" in label:
                    ok = small == 0 and smallx == 0
                else:
                    ok = small == 1 and smallx == 0
                if not ok or max(np.max(np.abs(z[:8])), abs(z[8])) > 30 or abs(z[8]) < 0.05:
                    continue
                got += 1
                worst = np.maximum(worst, claims(z))
            rp.line("E3 %s, eps = %.0e: max of the four errors" % (label, eps), got, float(np.max(worst)),
                    1e-6 if eps > 1e-3 else 1e-5,
                    "errors (first, six, Re F_0, F_4/F_5) = %.1e %.1e %.1e %.1e; seed %d" % (*worst, k0))

    # E4: full contact conditions (with inequalities), nondegenerate points.  The penalised solver only finds
    # degenerate points for T_1 (see E4x), so nondegenerate contact points are taken from equality solutions
    # whose 24 real parts happen to lie on one side of the common value (lam is replaced by -lam if needed).
    seed, want = 108, 300
    rng = np.random.default_rng(seed)
    got, tries, w_add, w_cone, exact8 = 0, 0, 0.0, 0.0, 0
    while got < want and tries < 400000:
        tries += 1
        z = solve_equalities(T1, rng)
        if z is None or not nondegenerate(z, T1, dmin):
            continue
        F, _ = F_and_jac(z)
        g = F.real - F[0].real
        sc = np.max(np.abs(g))
        if np.max(g) <= 1e-9 * sc:
            g = -g
        if np.min(g) < -1e-9 * sc:
            continue
        g = g / np.mean(g)
        got += 1
        w_add = max(w_add, abs(g[4]), abs(g[5]))
        w_cone = max(w_cone, cone_residual(np.maximum(g, 0.0)))
        if sorted(np.where(np.abs(g) < 1e-7)[0].tolist()) == list(range(8)):
            exact8 += 1
    rp.info("E4: nondegenerate points with contact set containing T_1 (all 24 real parts on one side; gap vector "
            "normalised to mean 1); seed %d, %d starts; %d of them have contact set exactly {0,..,7}"
            % (seed, tries, exact8))
    rp.line("E4a the gaps at indices 4 and 5 vanish", got, w_add, 1e-7)
    rp.line("E4b gap vector lies in the cone of the 34 minor profiles [relative distance, linear programme]", got, w_cone, 1e-7)

    # E4x: whatever the penalised solver returns for 'contact set contains T_1' (in practice degenerate points)
    seed, want = 118, 200
    rng = np.random.default_rng(seed)
    got, tries, w_cone, kinds = 0, 0, 0.0, {"repeated eigenvalue": 0, "zero contact product": 0, "neither": 0}
    while got < want and tries < 3000:
        tries += 1
        z = solve_contact(T1, rng, margin=0.0)
        if z is None:
            continue
        F, _ = F_and_jac(z)
        g = F.real - F[0].real
        if np.min(g) < -1e-9 or np.max(np.abs(F)) > 1e4:
            continue
        got += 1
        A, C = z[0:4], z[4:8]
        rep = min(min(abs(A[i] - A[j]), abs(C[i] - C[j])) for i, j in itertools.combinations(range(4), 2)) < 1e-6
        zero = min(abs(F[s]) for s in T1) < 1e-6
        kinds["repeated eigenvalue" if rep else ("zero contact product" if zero else "neither")] += 1
        w_cone = max(w_cone, cone_residual(np.maximum(g, 0.0)))
    rp.info("E4x: points returned by the penalised solver for 'contact set contains T_1' (seed %d, %d solves): %s"
            % (seed, tries, ", ".join("%s %d" % kv for kv in kinds.items())))
    rp.line("E4x gap vector lies in the cone of the 34 minor profiles", got, w_cone, 1e-7)

    # E5: the repeated-eigenvalue locus A2 = A3 (delegated by the document to the one-circle theorem)
    seed, want = 109, 300
    rng = np.random.default_rng(seed)
    side = [(unit(2) - unit(3), 0.0)]
    got, tries, w_cone, exactT, w45 = 0, 0, 0.0, 0, 0.0
    while got < want and tries < 5000:
        tries += 1
        z = solve_contact(T1, rng, side=side, margin=float(rng.choice([0.0, 0.05, 0.2])))
        if z is None:
            continue
        F, _ = F_and_jac(z)
        g = F.real - F[0].real
        if np.min(g) < -1e-9 or np.max(np.abs(F)) > 1e4:
            continue
        got += 1
        w_cone = max(w_cone, cone_residual(np.maximum(g, 0.0)))
        w45 = max(w45, min(abs(g[4]), abs(g[5])))
        if sorted(np.where(np.abs(g) < 1e-7)[0].tolist()) == T1:
            exactT += 1
    rp.info("E5: points with A2 = A3 and contact set containing T_1; seed %d, %d solves; %d of them have contact "
            "set exactly T_1; largest gap min(g_4, g_5) seen: %.3f (so the main-branch conclusion fails there, "
            "as expected)" % (seed, tries, exactT, w45))
    rp.line("E5 gap vector lies in the cone of the 34 minor profiles", got, w_cone, 1e-7)

    # E6: zero contact product A0 = C0, and A2 = C2
    for label, coef, k0 in (("A0 = C0", unit(0) - unit(4), 110), ("A2 = C2", unit(2) - unit(6), 111)):
        rng = np.random.default_rng(k0)
        got, tries, w_cone = 0, 0, 0.0
        while got < 200 and tries < 5000:
            tries += 1
            z = solve_contact(T1, rng, side=[(coef, 0.0)], margin=0.0)
            if z is None:
                continue
            F, _ = F_and_jac(z)
            g = F.real - F[0].real
            if np.min(g) < -1e-9 or np.max(np.abs(F)) > 1e4:
                continue
            got += 1
            w_cone = max(w_cone, cone_residual(np.maximum(g, 0.0)))
        rp.line("E6 zero contact product (%s): gap vector in the cone of minor profiles; seed %d" % (label, k0),
                got, w_cone, 1e-7)
    return rp.finish()


# ----------------------------------------------------------------------------------------------------
def run_T2():
    rp = Report("NUMERIC falsification attempts, case T_2 = {0,1,2,3,6,8}  (floating point; seeds printed)")
    J_MINOR = [0, 1, 2, 3, 4, 5, 6, 9, 10, 15, 16]
    MJ = np.array([c.F_PROF[j] for j in J_MINOR], dtype=float).T
    control(rp)

    def d79(z):
        F, _ = F_and_jac(z)
        return abs((F[7] - F[9]).real) / max(np.max(np.abs(F[T2])), abs(F[7]), abs(F[9]))

    def d02(z):
        A, C = z[0:4], z[4:8]
        return abs(z[8] * (A[0] - C[0]) * (A[3] - C[3]) * (A[1] - A[2]) * (C[1] - C[2]))

    # E1: equalities only; the region of the claim is D_02 != 0
    seed, want = 201, 2000
    rng = np.random.default_rng(seed)
    got, tries, worst = 0, 0, 0.0
    while got < want and tries < 200000:
        tries += 1
        z = solve_equalities(T2, rng)
        if z is None or max(np.max(np.abs(z[:8])), abs(z[8])) > 30:
            continue
        A, C = z[0:4], z[4:8]
        if min(abs(z[8]), abs(A[0] - C[0]), abs(A[3] - C[3]), abs(A[1] - A[2]), abs(C[1] - C[2])) < 0.05:
            continue
        got += 1
        worst = max(worst, d79(z))
    rp.info("E1: solutions of 'Re F equal on T_2' only (no inequalities) with the four factors of D_02 and |lam| "
            ">= 0.05 (no other nondegeneracy imposed); seed %d, %d starts" % (seed, tries))
    rp.line("E1 Re F_7 = Re F_9 [relative]", got, worst, 1e-7)

    # E2: approaching D_02 = 0
    for label, coef, k0 in (("A1 - A2 = eps", unit(1) - unit(2), 202), ("C1 - C2 = eps", unit(5) - unit(6), 203),
                            ("A0 - C0 = eps", unit(0) - unit(4), 204), ("A3 - C3 = eps", unit(3) - unit(7), 205)):
        for eps in (1e-2, 1e-4):
            rng = np.random.default_rng(k0)
            got, tries, worst = 0, 0, 0.0
            while got < 200 and tries < 100000:
                tries += 1
                z = solve_equalities(T2, rng, side=[(coef, eps)], tol=1e-15)
                if z is None or max(np.max(np.abs(z[:8])), abs(z[8])) > 30 or abs(z[8]) < 0.05:
                    continue
                A, C = z[0:4], z[4:8]
                fs = sorted([abs(A[0] - C[0]), abs(A[3] - C[3]), abs(A[1] - A[2]), abs(C[1] - C[2])])
                if fs[1] < 0.05:
                    continue
                got += 1
                worst = max(worst, d79(z))
            rp.line("E2 %s, eps = %.0e: Re F_7 = Re F_9" % (label, eps), got, worst,
                    1e-6 if eps > 1e-3 else 1e-5, "seed %d" % k0)

    # E3: ON the excluded loci the equality g_7 = g_9 is not claimed; record what happens there
    for label, coef, k0 in (("A1 = A2", unit(1) - unit(2), 206), ("C1 = C2", unit(5) - unit(6), 207)):
        rng = np.random.default_rng(k0)
        got, tries, big = 0, 0, 0.0
        while got < 200 and tries < 100000:
            tries += 1
            z = solve_equalities(T2, rng, side=[(coef, 0.0)])
            if z is None or max(np.max(np.abs(z[:8])), abs(z[8])) > 30 or abs(z[8]) < 0.05:
                continue
            got += 1
            big = max(big, d79(z))
        rp.info("E3 on the locus %s (equalities only, seed %d, %d samples): largest relative |Re(F_7 - F_9)| = "
                "%.3f, so the equality really needs D_02 != 0" % (label, k0, got, big))
        rng = np.random.default_rng(k0 + 10)
        got, tries, w_cone, w79 = 0, 0, 0.0, 0.0
        while got < 200 and tries < 5000:
            tries += 1
            z = solve_contact(T2, rng, side=[(coef, 0.0)], margin=0.0)
            if z is None:
                continue
            F, _ = F_and_jac(z)
            g = F.real - F[0].real
            if np.min(g) < -1e-9 or np.max(np.abs(F)) > 1e4:
                continue
            got += 1
            w_cone = max(w_cone, cone_residual(np.maximum(g, 0.0)))
            w79 = max(w79, abs(g[7] - g[9]))
        rp.line("E3 locus %s with contact set containing T_2: gap vector in the cone of minor profiles; seed %d"
                % (label, k0 + 10), got, w_cone, 1e-7, "largest |g_7 - g_9| among these points: %.3f" % w79)

    # E4: full contact conditions, generic
    seed, want = 220, 300
    rng = np.random.default_rng(seed)
    got, tries, w79, w_cone, w_cone_all, exactT, g7max = 0, 0, 0.0, 0.0, 0.0, 0, 0.0
    while got < want and tries < 5000:
        tries += 1
        z = solve_contact(T2, rng, margin=float(rng.choice([0.0, 0.05, 0.2])))
        if z is None or d02(z) < 1e-3:
            continue
        F, _ = F_and_jac(z)
        g = F.real - F[0].real
        if np.min(g) < -1e-9 or np.max(np.abs(F)) > 1e4:
            continue
        got += 1
        w79 = max(w79, abs(g[7] - g[9]) / max(1.0, np.max(np.abs(g))))
        w_cone = max(w_cone, cone_residual(np.maximum(g, 0.0), MJ))
        w_cone_all = max(w_cone_all, cone_residual(np.maximum(g, 0.0)))
        g7max = max(g7max, g[7])
        if sorted(np.where(np.abs(g) < 1e-7)[0].tolist()) == T2:
            exactT += 1
    rp.info("E4: points with contact set containing T_2 (inequalities enforced, mean gap off T = 1), |D_02| >= 1e-3; seed %d, "
            "%d solves; %d of them have contact set exactly T_2; largest g_7 seen %.3f" % (seed, tries, exactT, g7max))
    rp.line("E4a g_7 = g_9 [relative to max gap]", got, w79, 1e-7)
    rp.line("E4b gap vector in the cone of the 11 minor generators of N_{T_2} [relative distance, linear programme]", got,
            w_cone, 1e-7)
    rp.line("E4c gap vector in the cone of all 34 minor profiles", got, w_cone_all, 1e-7)
    return rp.finish()


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else ""
    if which == "T1":
        sys.exit(run_T1())
    if which == "T2":
        sys.exit(run_T2())
    print("usage: numeric_T1_T2.py T1|T2")
    sys.exit(2)
