"""Exact polynomial tools for the independent replay of the group V cases T_35, T_36, T_37.

All code here is new (nothing copied).  Exact integer polynomial arithmetic is done with
python-flint (class fmpz_mpoly: sparse multivariate polynomials over Z, arbitrary precision).
Only +, -, *, **, exact division, gcd and term access of that class are used.

Contents
  * PR        a polynomial ring Z[names] with conversion from/to {exponent tuple: int} dicts;
  * RF        reduced fractions of such polynomials (for the rational charts of the documents);
  * radial_data(order, names)   the "radial" objects of a group V proof, all rebuilt from the
              definitions: u = sum_i r_i f_{order[i]} + h kappa, the restriction C(u) = C* + h Q
              of the paper's cubic, the full 14 x 14 Hessian of C at u, its determinant D (a
              polynomial in the eight coefficients and h), N = Q^4 D(h = -C*/Q), and the
              polynomial G with N = -4 Q^2 G^2 (exact division and exact square root);
  * cleared(...)  substitution of a rational chart into a polynomial and clearing of the
              positive denominators;  sign statistics of the result;
  * chart bookkeeping (positivity of denominators on the stated parameter domain).
"""
from fractions import Fraction
import math

from flint import fmpz, fmpz_mpoly

import replay_T35_common as M


# ----------------------------------------------------------------------------------------
# polynomial rings over Z
# ----------------------------------------------------------------------------------------


class PR:
    """Z[names].  Elements are flint fmpz_mpoly objects with len(names) variables."""

    def __init__(self, names):
        self.names = tuple(names)
        self.n = len(self.names)
        assert self.n >= 1 and len(set(self.names)) == self.n
        self.gens = tuple(fmpz_mpoly.gens(self.n))
        self.zero = self.gens[0] * 0
        self.one = self.gens[0] ** 0
        self.idx = {nm: i for i, nm in enumerate(self.names)}

    def __getitem__(self, name):
        return self.gens[self.idx[name]]

    def const(self, c):
        return self.one * fmpz(int(c))

    def vars(self, *names):
        return [self[nm] for nm in names]

    @staticmethod
    def terms(p):
        """{exponent tuple: int coefficient}."""
        return {tuple(int(e) for e in p.exponent_tuple(i)): int(p.coefficient(i)) for i in range(len(p))}

    def from_terms(self, d):
        """Polynomial with the given {exponent tuple: coefficient} (Horner recursion)."""
        return self._horner(list(d.items()), 0, [(g, None) for g in self.gens], None)

    def _horner(self, items, v, subs, degs):
        if not items:
            return self.zero
        if v == len(subs):
            tot = 0
            for _, c in items:
                tot += c
            return self.const(tot)
        groups = {}
        for ex, c in items:
            groups.setdefault(ex[v], []).append((ex, c))
        num, den = subs[v]
        out = self.zero
        npow = {0: self.one}
        top = max(groups)
        cur = self.one
        for e in range(1, top + 1):
            cur = cur * num
            npow[e] = cur
        for e, sub in groups.items():
            inner = self._horner(sub, v + 1, subs, degs)
            if len(inner) == 0:
                continue
            t = inner * npow[e] if e else inner
            if den is not None and degs[v] - e > 0:
                t = t * (den ** (degs[v] - e))
            out = out + t
        return out

    def compose(self, d, subs):
        """sum_terms c * prod_v subs[v]^e_v for a term dict d over other variables;
        subs[v] are elements of this ring."""
        return self._horner(list(d.items()), 0, [(s, None) for s in subs], None)

    def compose_cleared(self, d, nums, dens, degs):
        """sum_terms c * prod_v nums[v]^e_v dens[v]^(degs[v]-e_v): the numerator of the
        polynomial with term dict d at x_v = nums[v]/dens[v], over prod_v dens[v]^degs[v]."""
        for ex in d:
            assert all(ex[v] <= degs[v] for v in range(len(degs)))
        return self._horner(list(d.items()), 0, list(zip(nums, dens)), list(degs))


def poly_stats(p):
    """(number of terms, number of negative coefficients, minimal coefficient)."""
    n = len(p)
    neg = 0
    mn = None
    for i in range(n):
        c = p.coefficient(i)
        if c < 0:
            neg += 1
        if mn is None or c < mn:
            mn = c
    return n, neg, (int(mn) if mn is not None else 0)


def nonneg(p):
    return all(p.coefficient(i) >= 0 for i in range(len(p)))


def degree_in(p, v):
    return max((p.exponent_tuple(i)[v] for i in range(len(p))), default=0)


def term_dict_degrees(d, nvars):
    return [max((ex[v] for ex in d), default=0) for v in range(nvars)]


def coeffs_in(d, v):
    """Split a term dict by the exponent of variable v: {e: term dict with exponent of v set to 0}."""
    out = {}
    for ex, c in d.items():
        e = ex[v]
        key = ex[:v] + (0,) + ex[v + 1:]
        out.setdefault(e, {})[key] = c
    return out


def set_var(d, v, value):
    """Term dict after substituting the integer `value` for variable v (exponent set to 0)."""
    out = {}
    for ex, c in d.items():
        key = ex[:v] + (0,) + ex[v + 1:]
        w = out.get(key, 0) + c * value ** ex[v]
        if w:
            out[key] = w
        else:
            out.pop(key, None)
    return out


def eval_terms(d, point):
    """Exact value of a term dict at a point (list of Fractions/ints)."""
    tot = 0
    for ex, c in d.items():
        t = c
        for x, e in zip(point, ex):
            if e:
                t = t * x ** e
        tot += t
    return tot


# ----------------------------------------------------------------------------------------
# reduced fractions of polynomials
# ----------------------------------------------------------------------------------------


class RF:
    """num/den with num, den in a PR ring, reduced by the gcd, den with positive leading
    coefficient.  Only used for the small rational expressions of the charts."""
    __slots__ = ("R", "num", "den")

    def __init__(self, R, num, den=None, reduce=True):
        self.R = R
        if den is None:
            den = R.one
        if isinstance(num, int):
            num = R.const(num)
        if isinstance(den, int):
            den = R.const(den)
        assert len(den) > 0, "zero denominator"
        if reduce and len(num) > 0:
            g = num.gcd(den)
            if not g.is_one():
                num, den = num // g, den // g
        elif len(num) == 0:
            den = R.one
        if den.leading_coefficient() < 0:
            num, den = -num, -den
        self.num, self.den = num, den

    def _c(self, o):
        return o if isinstance(o, RF) else RF(self.R, o)

    def __add__(self, o):
        o = self._c(o)
        return RF(self.R, self.num * o.den + o.num * self.den, self.den * o.den)
    __radd__ = __add__

    def __neg__(self):
        return RF(self.R, -self.num, self.den, reduce=False)

    def __sub__(self, o):
        return self + (-self._c(o))

    def __rsub__(self, o):
        return self._c(o) - self

    def __mul__(self, o):
        o = self._c(o)
        return RF(self.R, self.num * o.num, self.den * o.den)
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = self._c(o)
        return RF(self.R, self.num * o.den, self.den * o.num)

    def __rtruediv__(self, o):
        return self._c(o) / self

    def __pow__(self, k):
        return RF(self.R, self.num ** k, self.den ** k, reduce=False)

    def __eq__(self, o):
        o = self._c(o)
        return self.num * o.den == o.num * self.den

    def is_zero(self):
        return len(self.num) == 0

    def value(self, point):
        n = eval_terms(PR.terms(self.num), point)
        dd = eval_terms(PR.terms(self.den), point)
        return Fraction(n) / Fraction(dd)


def positive_on_domain(R, p, positive_vars):
    """Sufficient exact test that the polynomial p is > 0 whenever all variables are >= 0
    and the variables named in positive_vars are > 0: all coefficients are >= 0 and some
    monomial involves only variables of positive_vars."""
    ok_mono = False
    pv = set(R.idx[nm] for nm in positive_vars)
    for i in range(len(p)):
        if p.coefficient(i) < 0:
            return False
        ex = p.exponent_tuple(i)
        if all((e == 0) or (v in pv) for v, e in enumerate(ex)):
            ok_mono = True
    return ok_mono


def cleared(R, d, src_names, chart, degs=None):
    """Numerator of the polynomial with term dict d (variables src_names) after the
    substitution chart[name] = RF, over the denominator prod_v den_v^deg_v, where deg_v is the
    degree of d in the variable (or the given degs).  Returns (numerator, list of the
    denominators with their powers)."""
    nv = len(src_names)
    if degs is None:
        degs = term_dict_degrees(d, nv)
    nums = [chart[nm].num for nm in src_names]
    dens = [chart[nm].den for nm in src_names]
    num = R.compose_cleared(d, nums, dens, degs)
    return num, [(dens[v], degs[v]) for v in range(nv)]


# ----------------------------------------------------------------------------------------
# radial objects of a group V proof
# ----------------------------------------------------------------------------------------


def hessian_det(H, one):
    """Determinant of a square matrix of fmpz_mpoly entries by Laplace expansion along rows
    with memoisation over column subsets (rows ordered by sparsity)."""
    n = len(H)
    rows = sorted(range(n), key=lambda i: sum(1 for j in range(n) if len(H[i][j])))
    inv = sum(1 for a in range(n) for b in range(a + 1, n) if rows[a] > rows[b])
    dp = {0: one}
    for k, r in enumerate(rows):
        new = {}
        nzcols = [j for j in range(n) if len(H[r][j])]
        for mask, poly in dp.items():
            for j in nzcols:
                if mask >> j & 1:
                    continue
                nm = mask | (1 << j)
                pos = bin(nm & ((1 << j) - 1)).count("1")
                term = poly * H[r][j]
                if (k + pos) % 2:
                    term = -term
                new[nm] = new[nm] + term if nm in new else term
        dp = {m: p for m, p in new.items() if len(p)}
    D = dp.get((1 << n) - 1, one * 0)
    return -D if inv % 2 else D


def det_bareiss_check(H, R, D, points):
    """Independent exact check of a symbolic determinant: evaluate the matrix at integer
    points and compare the integer determinant (Fraction elimination) with D at the point."""
    ok = True
    Dt = PR.terms(D)
    for pt in points:
        A = [[Fraction(eval_terms(PR.terms(H[i][j]), pt)) for j in range(len(H))] for i in range(len(H))]
        n = len(A)
        det = Fraction(1)
        for c in range(n):
            k = next((i for i in range(c, n) if A[i][c] != 0), None)
            if k is None:
                det = Fraction(0)
                break
            if k != c:
                A[c], A[k] = A[k], A[c]
                det = -det
            det *= A[c][c]
            for i in range(c + 1, n):
                if A[i][c] != 0:
                    f = A[i][c] / A[c][c]
                    A[i] = [a - f * b for a, b in zip(A[i], A[c])]
        ok &= det == eval_terms(Dt, pt)
    return ok


def poly_sqrt(R, S, dvar):
    """Exact square root of S in the ring R, assuming S = G^2 with G homogeneous-like in the
    grading "total degree in the variables other than dvar" and G_0 = dvar^m (coefficient +1).
    Returns G or None.  The result is verified by squaring by the caller."""
    v = R.idx[dvar]
    St = PR.terms(S)
    if not St:
        return None
    # grade by total degree in the other variables
    grade = {}
    for ex, c in St.items():
        k = sum(ex) - ex[v]
        grade.setdefault(k, {})[ex] = c
    kmax = max(grade)
    if kmax % 2:
        return None
    half = kmax // 2
    S0 = grade.get(0, {})
    if len(S0) != 1:
        return None
    (ex0, c0), = S0.items()
    if c0 != 1 or ex0[v] % 2:
        return None
    m = ex0[v] // 2
    G = [None] * (half + 1)
    G[0] = R.gens[v] ** m
    Gp = [R.from_terms(grade.get(k, {})) for k in range(kmax + 1)]
    dm = R.gens[v] ** m
    for k in range(1, half + 1):
        acc = Gp[k]
        for i in range(1, k):
            acc = acc - G[i] * G[k - i]
        # acc = 2 * dvar^m * G_k
        t = PR.terms(acc)
        new = {}
        for ex, c in t.items():
            if c % 2 or ex[v] < m:
                return None
            new[ex[:v] + (ex[v] - m,) + ex[v + 1:]] = c // 2
        G[k] = R.from_terms(new)
    tot = R.zero
    for g in G:
        tot = tot + g
    return tot


def radial_data(order, names, dname="d"):
    """order: the eight generator indices of a normal form, names: the eight coefficient
    names in the same order.  Everything is computed from the paper's cubic and from the
    coordinates of the f_j (module replay_T35_common); nothing is read from any file."""
    assert len(order) == 8 and len(names) == 8
    R = PR(list(names) + ["h"])
    hv = 8
    # u_k = sum_i r_i (f_order[i])_k + h kappa_k
    u = []
    for k in range(14):
        co = [M.F_COORD[j][k] for j in order] + [M.KAPPA[k]]
        assert all(Fraction(x).denominator == 1 for x in co)
        e = R.zero
        for c, g in zip(co, R.gens):
            if c:
                e = e + g * int(c)
        u.append(e)
    Cu = M.cubic_C(u)
    parts = coeffs_in(PR.terms(Cu), hv)
    lin_in_h = set(parts) <= {0, 1}
    Cstar = R.from_terms(parts.get(0, {}))
    Q = R.from_terms(parts.get(1, {}))
    # Hessian: entries are linear forms in u
    M.hess_at([Fraction(0)] * 14)
    table = M._CUBIC_CACHE["hess_lin"]
    H = [[R.zero] * 14 for _ in range(14)]
    for i in range(14):
        for j in range(14):
            e = R.zero
            for cf, k in table[i][j]:
                e = e + u[k] * cf
            H[i][j] = e
    D = hessian_det(H, R.one)
    Dparts = coeffs_in(PR.terms(D), hv)
    degh = max(Dparts) if Dparts else 0
    # N = sum_k D_k (-C*)^k Q^(4-k)
    N = R.zero
    if degh <= 4:
        for k, dk in Dparts.items():
            N = N + R.from_terms(dk) * ((-Cstar) ** k) * (Q ** (4 - k))
    return dict(R=R, u=u, C=Cu, lin_in_h=lin_in_h, Cstar=Cstar, Q=Q, H=H, D=D, degh=degh, N=N,
                order=list(order), names=list(names), dname=dname)


def radial_factor(rd):
    """Exact division N / (-4 Q^2) and exact square root.  Adds S, G (normalised so that the
    coefficient of d^9 is +1) and the flags div_ok, sqrt_ok to the dict."""
    R, N, Q = rd["R"], rd["N"], rd["Q"]
    den = Q * Q * (-4)
    S = N // den
    rd["div_ok"] = (S * den == N)
    rd["S"] = S
    G = poly_sqrt(R, S, rd["dname"]) if rd["div_ok"] else None
    rd["sqrt_ok"] = G is not None and (G * G == S)
    rd["G"] = G
    return rd


def drop_h_set_d(rd, p, dvalue=1):
    """Term dict of a polynomial of the radial ring that does not involve h, as a term dict in
    the seven variables other than d (d set to dvalue).  Returns (names7, dict)."""
    R = rd["R"]
    t = PR.terms(p)
    assert all(ex[8] == 0 for ex in t)
    dv = R.idx[rd["dname"]]
    t = set_var(t, dv, dvalue)
    keep = [i for i in range(8) if i != dv]
    out = {}
    for ex, c in t.items():
        out[tuple(ex[i] for i in keep)] = c
    return [R.names[i] for i in keep], out


def load_supplement_G(path, names9, alias=None):
    """Read the degree-nine radial factor stored in a JSON file of the supplement AS DATA.
    alias maps a name of names9 to the name used in the file.  The ninth variable (the height h,
    called 'z' in the files) does not occur in the factor.  Returns a dict or None.  Nothing is
    executed."""
    import json
    import os
    if not os.path.isfile(path):
        return None
    with open(path) as fh:
        J = json.load(fh)
    var = list(J["variables"])
    fac = [f for f in J["radial_factors"] if f.get("degree") == 9]
    if len(fac) != 1:
        return None
    alias = alias or {}
    perm = []
    for nm in names9[:8]:
        nm2 = alias.get(nm, nm)
        if nm2 not in var:
            return None
        perm.append(var.index(nm2))
    rest = [i for i in range(len(var)) if i not in perm]
    out = {}
    for ex, c in fac[0]["terms"]:
        if any(ex[i] for i in rest):
            return None
        out[tuple(ex[i] for i in perm) + (0,)] = int(c)
    return dict(variables=var, terms=out, multiplicity=fac[0].get("multiplicity"),
                constant=J.get("constant"))


# ----------------------------------------------------------------------------------------
# helpers shared by the three replay scripts
# ----------------------------------------------------------------------------------------


def prod_dens(R, dens):
    W = R.one
    for den, k in dens:
        if k:
            W = W * den ** k
    return W


def restrict_zero(R, p, zero_names):
    """The part of p that does not involve the variables in zero_names (p at those variables = 0)."""
    idx = [R.idx[nm] for nm in zero_names]
    t = {ex: c for ex, c in PR.terms(p).items() if all(ex[i] == 0 for i in idx)}
    return R.from_terms(t)


def chart_denominators_positive(R, chart, positive_vars):
    """Every denominator of the chart is > 0 when the variables are >= 0 and those named in
    positive_vars are > 0, and every numerator has nonnegative coefficients (so every chart
    value is a nonnegative real number on the domain)."""
    ok = True
    for nm, rf in chart.items():
        ok &= positive_on_domain(R, rf.den, positive_vars)
        ok &= nonneg(rf.num)
    return ok


def chart_point(R, chart, rng, positive_vars, zero_prob=0.0, big=6):
    """Random rational point of the chart domain and the chart values there."""
    pt = []
    for nm in R.names:
        if nm not in positive_vars and rng.random() < zero_prob:
            pt.append(Fraction(0))
        else:
            pt.append(Fraction(rng.randint(1, 40), rng.randint(1, 40)) ** rng.choice([1, 1, 2])
                      * rng.choice([1, 1, 1, big, Fraction(1, big)]))
    vals = {nm: rf.value(pt) for nm, rf in chart.items()}
    return pt, vals


def bernstein(R, ks, names, chart, thr, degs):
    """ks = (k0, k1, k2): term dicts of the coefficients of a quadratic in its central variable;
    thr = RF threshold N/D.  Returns the three cleared Bernstein numerators
    (K0, 2 D K0 + N K1, D^2 K0 + N D K1 + N^2 K2) and the common positive multiplier W."""
    K = []
    W = None
    for k in ks:
        num, dens = cleared(R, k, names, chart, degs=degs)
        K.append(num)
        W = prod_dens(R, dens)
    Nn, Dd = thr.num, thr.den
    mid = K[0] * Dd * 2 + K[1] * Nn
    end = K[0] * Dd * Dd + K[1] * Nn * Dd + K[2] * Nn * Nn
    return K[0], mid, end, W


def check_radial(ck, rd, label, expect_terms, json_name, alias=None):
    """Generic exact checks of the radial objects; returns (names7, term dict of G at d = 1)."""
    import random
    R = rd["R"]
    ck.check("EXACT", "%s: C(g + h kappa) is of degree <= 1 in h; C* = C(g) has %d terms, Q = dC/dh has %d "
             "terms" % (label, len(rd["Cstar"]), len(rd["Q"])), rd["lin_in_h"])
    rng = random.Random(12345)
    pts = [[rng.randint(-5, 5) for _ in range(9)] for _ in range(3)]
    ck.check("EXACT", "%s: the full 14 x 14 Hessian of C at g + h kappa has determinant D with %d terms, of "
             "degree %d in h (symbolic Laplace expansion; confirmed at 3 random integer points by "
             "rational Gaussian elimination)" % (label, len(rd["D"]), rd["degh"]),
             rd["degh"] <= 4 and det_bareiss_check(rd["H"], R, rd["D"], pts))
    radial_factor(rd)
    G = rd["G"]
    ck.check("EXACT", "%s: N := Q^4 D(h = -C*/Q) = sum_k D_k (-C*)^k Q^(4-k) (%d terms) is exactly divisible by "
             "-4 Q^2" % (label, len(rd["N"])), rd["div_ok"])
    dv = R.idx[rd["dname"]]
    homog = G is not None and all(sum(ex) == 9 for ex in PR.terms(G))
    d9 = None if G is None else PR.terms(G).get(tuple(9 if i == dv else 0 for i in range(9)))
    ck.check("EXACT", "%s: N / (-4 Q^2) = G^2 for a polynomial G (own exact square root, verified by squaring): "
             "identity (2) of the document; G is homogeneous of degree 9 with %s terms (document: %d) and "
             "coefficient of d^9 equal to +1" % (label, None if G is None else len(G), expect_terms),
             rd["sqrt_ok"] and homog and d9 == 1 and len(G) == expect_terms)
    path = M.supplement_path(M.DOC_DIR, json_name)
    if path is None:
        ck.info("MOC4_SUPP not set: G is not compared with the supplement's stored polynomial")
    else:
        names9 = list(R.names)
        J = load_supplement_G(path, names9, alias)
        if J is None:
            ck.check("EXACT", "%s: supplement data file %s readable" % (label, json_name), False)
        else:
            Gj = R.from_terms(J["terms"])
            same = (Gj == G)
            opp = (len(Gj + G) == 0)
            ck.check("EXACT", "%s: the degree-nine factor stored in the supplement's %s (read as data, %d terms) "
                     "equals %s G" % (label, json_name, len(J["terms"]), "+" if same else "-"), same or opp)
    return drop_h_set_d(rd, G)
