"""Certificate replay (auditor cert-T39) for independent_moc4_T39_high_p_uniform_20261022.

Own exact code. Every polynomial is rebuilt from the paper's definitions through the audit
library (cubic C of (A.2), facet vectors f_j of Appendix D, kappa = e_13) with python-flint.
Nothing from the supplement is executed or imported; its JSON files are loaded as DATA only
(to compare term counts and coefficients), from the directory given by MOC4_SUPP.

This file also holds the shared CORE (Hessian, determinant identity, G, octics, principal
minor, pseudo-division, trace-norm identity, bordered minors) that the low_pj and high_j
replays import from this same directory.

usage: MOC4_SUPP=<supplement root> python -I -B this_file.py
Exit code nonzero if any EXACT check fails.
"""
import json
import os
import sys
import time

import flint

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(os.path.dirname(HERE), "lib")
sys.path.insert(0, LIB)
import moc4lib as L  # noqa: E402

T0 = time.monotonic()
FAILS = []
SUPP_DIR_REL = os.path.join("proofs", "manin_transfer_20261021")


def check(label, ok, extra=""):
    print(("PASS " if ok else "FAIL ") + label + ((" | " + extra) if extra else ""), flush=True)
    if not ok:
        FAILS.append(label)
    return ok


def info(msg):
    print("INFO " + msg + "  [t=%.1fs]" % (time.monotonic() - T0), flush=True)


# ------------------------------------------------------------------------------------
# Base ring: q p j c e b a h   (h = height along kappa = e_13); d normalized to 1.
# ------------------------------------------------------------------------------------
VARS = ["q", "p", "j", "c", "e", "b", "a", "h"]
CTX = flint.fmpz_mpoly_ctx.get(VARS)
q, p, j, c, e, b, a, h = CTX.gens()
IQ, IP, IJ, IC, IE, IB, IA, IH = range(8)
ONE = CTX.constant(1)
SOURCE = [5, 1, 12, 9, 2, 35, 33, 29]          # canonical source of the T_40 chamber papers
BUDGET = [q, p, j, c, e, ONE, b, a]            # coefficient of f_5, f_1, f_12, f_9, f_2, f_35(d=1), f_33, f_29


def coeff(poly, var, n):
    """Coefficient of var^n in poly (as a polynomial with var-exponent 0)."""
    out = {}
    for ex, co in poly.to_dict().items():
        if ex[var] == n:
            ex2 = list(ex)
            ex2[var] = 0
            out[tuple(ex2)] = co
    return CTX.from_dict(out) if out else CTX.constant(0)


def deg(poly, var):
    return int(poly.degrees()[var]) if not poly.is_zero() else -1


def bareiss_det(M):
    """Fraction-free Gaussian elimination with full pivoting (smallest pivot), exact division."""
    M = [row[:] for row in M]
    n = len(M)
    prev = CTX.constant(1)
    sign = 1
    for k in range(n - 1):
        best = None
        for i in range(k, n):
            for jj in range(k, n):
                if not M[i][jj].is_zero():
                    ln = len(M[i][jj])
                    if best is None or ln < best[0]:
                        best = (ln, i, jj)
        if best is None:
            return CTX.constant(0)
        _, i, jj = best
        if i != k:
            M[k], M[i] = M[i], M[k]
            sign = -sign
        if jj != k:
            for row in M:
                row[k], row[jj] = row[jj], row[k]
            sign = -sign
        piv = M[k][k]
        for i in range(k + 1, n):
            for jj in range(k + 1, n):
                num = piv * M[i][jj] - M[i][k] * M[k][jj]
                if k:
                    qt, rm = divmod(num, prev)
                    assert rm.is_zero(), "Bareiss division not exact"
                    M[i][jj] = qt
                else:
                    M[i][jj] = num
            M[i][k] = CTX.constant(0)
        prev = piv
    return sign * M[n - 1][n - 1]


def project_height(Dpoly, C, Q):
    """Q^m * Dpoly(h = -C/Q) with m = deg_h Dpoly; returns (value, m)."""
    m = deg(Dpoly, IH)
    out = CTX.constant(0)
    for i in range(m + 1):
        out += coeff(Dpoly, IH, i) * (-C) ** i * Q ** (m - i)
    return out, m


def core(full=True):
    """Rebuild every shared polynomial.  full=False stops after G, F+-, (for certificates)."""
    R = {}
    # point w = g + h kappa in the paper's matching coordinates m_0..m_13
    w = [CTX.constant(0)] * 14
    for co, jdx in zip(BUDGET, SOURCE):
        vec = L.F_VECS[jdx]
        assert all(v.denominator == 1 for v in vec)
        for k in range(14):
            if vec[k]:
                w[k] = w[k] + int(vec[k]) * co
    w[13] = w[13] + h                                   # kappa = e_13 (paper (3.2))
    check("kappa: P kappa = (1,...,1)", L.values(L.KAPPA) == [1] * 24)
    # source rays vanish exactly on T_40 = {0,1,5,9,14,18}
    common = frozenset(range(24))
    for jdx in SOURCE:
        common = common & L.F_ZERO[jdx]
    check("the eight source rays vanish together exactly on T_40 = {0,1,5,9,14,18}",
          sorted(common) == [0, 1, 5, 9, 14, 18], str(sorted(common)))
    # cubic of (A.2) as a 14-variable integer polynomial, then restricted
    C14 = flint.fmpz_mpoly_ctx.get(["m%d" % i for i in range(14)])
    mg = C14.gens()
    cub = L.cubic_C(list(mg))
    Fw = cub.compose(*w, ctx=CTX)
    check("C(g + h kappa) has degree 1 in h", deg(Fw, IH) == 1)
    Cc, Qc = coeff(Fw, IH, 0), coeff(Fw, IH, 1)
    A_ = 1 + a + b
    B_ = 1 + q + e
    N_ = (q + 1) * (e + 1) * A_ + B_ * a * b
    Nj = N_ - j * e * a
    check("C = -(A B c + N_j)", Cc == -(A_ * B_ * c + Nj))
    Qdisp = -((q + e + 1 + a + b - p) * c + (q + b + 1) * (e + a + 1) - j * (e + a) + p * j)
    check("Q = -{(q+e+1+a+b-p)c + (q+b+1)(e+a+1) - j(e+a) + pj}", Qc == Qdisp)
    R.update(C=Cc, Q=Qc, A=A_, B=B_, N=N_, Nj=Nj)
    # Hessian of the cubic at w (entries are linear in q..a,h)
    H = [[cub.derivative(r).derivative(s).compose(*w, ctx=CTX) for s in range(14)]
         for r in range(14)]
    check("Hessian symmetric", all(H[r][s] == H[s][r] for r in range(14) for s in range(14)))
    hk = [[coeff(H[r][s], IH, 1) for s in range(14)] for r in range(14)]
    nzk = sorted((r, s) for r in range(14) for s in range(14) if not hk[r][s].is_zero())
    check("H(kappa) = height part of H(w) has entries only at (0,6),(1,5) and transposes",
          nzk == [(0, 6), (1, 5), (5, 1), (6, 0)], str(nzk))
    info("Hessian built")
    D = bareiss_det(H)
    info("det H: %d terms, deg_h %d, deg_c %d" % (len(D), deg(D, IH), deg(D, IC)))
    R["detH_terms"] = len(D)
    proj, m = project_height(D, Cc, Qc)
    check("deg_h det H(g + h kappa) = 3", m == 3, "deg_h = %d" % m)
    G2, rem = divmod(proj, -4 * Qc ** (m - 2))
    check("Q^3 det H(-C/Q) divisible by -4 Q", rem.is_zero() and G2 * (-4 * Qc ** (m - 2)) == proj)
    G = G2.sqrt()          # raises if not a perfect square
    check("quotient is a perfect square G^2: identity (1) det H(w(c)) = -4 G^2/Q^2",
          G * G == G2)
    const = G.to_dict().get((0,) * 8, 0)
    check("constant term of G is +-1", abs(int(const)) == 1, "const = %s" % const)
    if int(const) == 1:
        G = -G                                     # normalization: constant term -1
    check("G: 1313 terms, deg_c 2, total degree 9, height free",
          len(G) == 1313 and deg(G, IC) == 2 and G.total_degree() == 9 and deg(G, IH) <= 0,
          "%d terms, deg_c %d, tdeg %d" % (len(G), deg(G, IC), G.total_degree()))
    k0, k1, k2 = (coeff(G, IC, i) for i in range(3))   # gamma, beta, alpha
    R.update(G=G, gamma=k0, beta=k1, alpha=k2)
    info("G rebuilt")
    # derivative identity (2)/(3): d/dh det H at h=-C/Q equals 4 G P6 / Q
    Dd = CTX.constant(0)
    for i in range(1, m + 1):
        Dd += i * coeff(D, IH, i) * (-Cc) ** (i - 1) * Qc ** (m - i)   # = Q^(m-1) D'(-C/Q)
    P6, rem = divmod(Dd, 4 * G * Qc ** (m - 2))
    check("Q^2 D'(-C/Q) = 4 G P6 Q exactly (so D'(-C/Q) = 4 G P6/Q), P6 has 464 terms",
          rem.is_zero() and P6 * 4 * G * Qc ** (m - 2) == Dd and len(P6) == 464,
          "P6 %d terms" % len(P6))
    # discriminant and octics
    disc = k1 ** 2 - 4 * k0 * k2
    unit, facs = disc.factor()
    ok = int(unit) == 1 and len(facs) == 2 and all(int(mu) == 1 and len(f) == 306 for f, mu in facs)
    check("beta^2 - 4 alpha gamma = F+ F- (unit 1, two simple 306-term factors)", ok,
          "unit %s, factors %s" % (unit, [(len(f), int(mu), f.total_degree()) for f, mu in facs]))
    Fa, Fb = facs[0][0], facs[1][0]
    check("octics: product equals the discriminant, each total degree 8",
          Fa * Fb == disc and Fa.total_degree() == 8 and Fb.total_degree() == 8)
    R.update(Fa=Fa, Fb=Fb, disc=disc)
    info("discriminant factored")
    if not full:
        return R
    # principal minor P* deleting indices 0,1
    I = list(range(2, 14))
    Pst = bareiss_det([[H[r][s] for s in I] for r in I])
    check("P* (delete 0,1): height independent, degree 4 in c, 5188 terms",
          deg(Pst, IH) <= 0 and deg(Pst, IC) == 4 and len(Pst) == 5188,
          "%d terms, deg_c %d" % (len(Pst), deg(Pst, IC)))
    info("P* built")
    # pseudo-division alpha^3 P* = Quot*G + r1 c + r0, with Quot tracked
    rr = [coeff(Pst, IC, n) for n in range(5)]
    quot = CTX.constant(0)
    for power in [4, 3, 2]:
        lead = rr[power]
        rr = [x * k2 for x in rr]
        quot = quot * k2 + lead * c ** (power - 2)
        for i in range(3):
            rr[power - 2 + i] -= lead * [k0, k1, k2][i]
        assert rr[power].is_zero()
    r0, r1 = rr[0], rr[1]
    lhs = k2 ** 3 * Pst
    check("pseudo-division: alpha^3 P* = Quot * G + r1 c + r0 (full polynomial identity)",
          lhs == quot * G + r1 * c + r0, "r0 %d terms, r1 %d terms" % (len(r0), len(r1)))
    info("pseudo-division done")
    unit, rf = r1.factor()
    ok = int(unit) == 1 and len(rf) == 2 and all(int(mu) == 1 and len(f) == 13576 for f, mu in rf)
    check("r1 = L M with two simple 13576-term factors (unit 1)", ok,
          "unit %s, factors %s" % (unit, [(len(f), int(mu)) for f, mu in rf]))
    L1, M1 = rf[0][0], rf[1][0]
    check("r1 == L*M exactly", L1 * M1 == r1)
    info("r1 factored")
    T = 2 * k2 * r0 - k1 * r1
    # identify which octic pairs with which square
    if 2 * T == Fa * M1 ** 2 + Fb * L1 ** 2:
        Fp, Fm, Lp, Mp = Fa, Fb, L1, M1
    elif 2 * T == Fa * L1 ** 2 + Fb * M1 ** 2:
        Fp, Fm, Lp, Mp = Fa, Fb, M1, L1
    else:
        Fp = None
    check("trace identity 2T = F+ M^2 + F- L^2, T = 2 alpha r0 - beta r1 (for one labelling)", Fp is not None)
    if Fp is not None:
        Z = Fp * Mp ** 2 - Fm * Lp ** 2
        check("Z = F+ M^2 - F- L^2 is not the zero polynomial", not Z.is_zero(), "%d terms" % len(Z))
        # norm identity as polynomial identity: T^2 - disc r1^2 = Z^2/4
        check("norm identity 4(T^2 - (beta^2-4 alpha gamma) r1^2) = Z^2",
              4 * (T * T - disc * r1 * r1) == Z * Z)
        check("Vieta numerator: T^2 - disc r1^2 = 4 alpha (alpha r0^2 - beta r0 r1 + gamma r1^2)",
              T * T - disc * r1 * r1 == 4 * k2 * (k2 * r0 * r0 - k1 * r0 * r1 + k0 * r1 * r1))
        R.update(Z_terms=len(Z))
    R.update(Pst=Pst, r0=r0, r1=r1, L=L1, M=M1, P6=P6)
    info("trace-norm identity done")
    # bordered 13-minors on I u {0}, I u {1}
    for added, want in [(0, 2021), (1, 1091)]:
        idx = I + [added]
        Bd = bareiss_det([[H[r][s] for s in idx] for r in idx])
        pr, md = project_height(Bd, Cc, Qc)
        qt, rm = divmod(pr, G)
        check("bordered minor on I u {%d}: height degree 2; Q^2 * minor(-C/Q) = G * quotient, quotient %d terms"
              % (added, want),
              md == 2 and rm.is_zero() and qt * G == pr and len(qt) == want,
              "raw %d terms, projected %d, quotient %d" % (len(Bd), len(pr), len(qt)))
    info("bordered minors done")
    return R


# ------------------------------------------------------------------------------------
# substitution helpers
# ------------------------------------------------------------------------------------
def substitute(pol, octx, images, dens, degs=None):
    """pol in CTX (no c, h).  images[var] = numerator poly in octx, dens[var] = denominator
    poly in octx (or None).  Returns pol(images/dens) * prod dens[var]^deg_var(pol)."""
    d = pol.to_dict()
    assert all(ex[IC] == 0 and ex[IH] == 0 for ex in d)
    degsl = degs or [deg(pol, v) for v in range(8)]
    cache = {}

    def pw(v, k):
        key = (v, k)
        if key not in cache:
            num = images[v] ** k
            if dens[v] is not None:
                num = num * dens[v] ** (degsl[v] - k)
            cache[key] = num
        return cache[key]
    out = octx.constant(0)
    for ex, co in d.items():
        t = octx.constant(int(co))
        for v in range(8):
            if v in (IC, IH):
                continue
            if images[v] is None:
                assert ex[v] == 0
                continue
            t = t * pw(v, ex[v])
        out += t
    return out


def all_negative(poly):
    vals = [int(x) for x in poly.to_dict().values()]
    return len(vals) > 0 and all(v < 0 for v in vals)


def all_nonneg(poly):
    return all(int(x) >= 0 for x in poly.to_dict().values())


def load_supp_json(name):
    root = os.environ.get("MOC4_SUPP") or (sys.argv[1] if len(sys.argv) > 1 else None)
    if not root:
        return None
    path = os.path.join(root, SUPP_DIR_REL, name)
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return json.load(fh)


def compare_with_data(label, mine, octx_vars, data_cert):
    """Compare my polynomial with a coefficient table from the supplement JSON (data only)."""
    if data_cert is None:
        print("INFO no supplement data for " + label)
        return
    theirs = {tuple(ex): int(co) for ex, co in data_cert}
    ours = {tuple(int(x) for x in ex): int(co) for ex, co in mine.to_dict().items()}
    check("data comparison: my %s equals the supplement's stored table (%d terms)" % (label, len(theirs)),
          ours == theirs)


def compare_factor_pair(label, mine_pair, data_pair):
    if data_pair[0] is None:
        print("INFO no supplement data for " + label)
        return
    def canon(dd):   # sign-normalized frozen form
        items = sorted(dd.items())
        s = 1 if items[0][1] > 0 else -1
        return frozenset((k, s * v) for k, v in items)
    th = {canon({tuple(ex): int(co) for ex, co in dd}) for dd in data_pair}
    ou = {canon({tuple(int(x) for x in ex): int(co) for ex, co in m.to_dict().items()}) for m in mine_pair}
    match = th == ou
    check("data comparison: my %s equal the supplement's stored pair up to order/sign" % label, match)


# ------------------------------------------------------------------------------------
# high-p chamber: p >= 1, 0 <= pj < 1 (d = 1); chart p = 1+P, j = J/((1+P)(1+J))
# ------------------------------------------------------------------------------------
def main():
    R = core(full=True)
    A_, B_, N_, Nj, k0, k1, k2, Fa, Fb = (R[k] for k in ["A", "B", "N", "Nj", "gamma", "beta", "alpha", "Fa", "Fb"])
    check("N - ea has nonnegative coefficients and constant 1",
          all_nonneg(N_ - e * a) and (N_ - e * a).to_dict().get((0,) * 8) == 1)
    octx = flint.fmpz_mpoly_ctx.get(["q", "P", "J", "e", "b", "a"])
    oq, oP, oJ, oe, ob, oa = octx.gens()
    one = octx.constant(1)
    images = [oq, one + oP, oJ, None, oe, ob, oa, None]
    dens = [None, None, (one + oP) * (one + oJ), None, None, None, None, None]
    data = load_supp_json("independent_moc4_T39_high_p_uniform_exact_20261022.json")
    stored = {cc["name"]: cc["coefficients"] for cc in data["strict_certificates"]} if data else {}
    certs = [("gamma", k0, 832, -2),
             ("AB_beta_minus_Nj_alpha", A_ * B_ * k1 - Nj * k2, 3084, -2),
             ("AB_Fsum_plus_2N_alpha", A_ * B_ * (Fa + Fb) + 2 * N_ * k2, 1821, -4)]
    for name, pol, want, const in certs:
        out = substitute(pol, octx, images, dens)
        cst = int(out.to_dict().get((0,) * 6, 0))
        check("high-p certificate %s: %d terms (claimed %d), every coefficient < 0, constant %d (claimed %d)"
              % (name, len(out), want, cst, const),
              len(out) == want and all_negative(out) and cst == const)
        compare_with_data(name, out, None, stored.get(name))
    # cover of the chamber by the chart
    print("INFO chart inverse: P = p - 1 >= 0, J = pj/(1 - pj) >= 0 for p >= 1, 0 <= pj < 1; "
          "then (1+P) = p and J/((1+P)(1+J)) = (pj)/p = j.  All P, J >= 0 occur (p>=1, j>=0, pj<1).")
    # the octic pair and L, M against stored data
    if data:
        compare_factor_pair("octics F+, F-", [Fa, Fb], [data.get("Fa_terms"), data.get("Fb_terms")])
        compare_factor_pair("factors L, M of r1", [R["L"], R["M"]], [data.get("U16_terms"), data.get("V16_terms")])
    print("SUMMARY %d failures; elapsed %.1f s" % (len(FAILS), time.monotonic() - T0))
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
