"""Gap check (second pass) for independent_moc4_T47_high_j_uniform_20261022, Theorem 2,
last range (e >= L^2/j, pc <= a^2, Phi >= 0): replay of the two-guard identity
    42 * Gtilde = R + M_pc * num(a^2 - pc) + M_Phi * num(Phi)
with own code (python-flint).  Only the supplement's JSON data file holding the polynomial G
is loaded (json.load, as data).  The two multipliers (5 and 6 monomials, weights over the
common denominator 42) are transcribed below from the companion verifier named in the
document's dependency paragraph (IIIb_low_j_stationary_exact); they are data, checked here.

Usage:  MOC4_SUPP=<supplement root> python -I -B this_file.py
Exits nonzero if an EXACT check fails.
"""
import os, sys, json, random
from fractions import Fraction
from flint import fmpz_mpoly_ctx

SUPP = os.environ.get("MOC4_SUPP") or (sys.argv[1] if len(sys.argv) > 1 else None)
if not SUPP:
    sys.exit("set MOC4_SUPP")
GJSON = os.path.join(SUPP, "proofs", "manin_transfer_20261021", "root_moc4_T47_q_squareclass_20261022.json")

VARS = ["A", "P", "T", "U", "V", "S", "Q"]
ctx = fmpz_mpoly_ctx.get(VARS)
A, P, T, U, V, S, Q = ctx.gens()
one = ctx.constant(1)
FAIL = []


def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + (" : " + info if info else ""), flush=True)
    if not ok:
        FAIL.append(name)


class RF:
    """reduced fraction n/d of fmpz_mpoly"""
    def __init__(s, n, d=None):
        if d is None:
            d = one
        g = n.gcd(d)
        s.n, s.d = n // g, d // g

    def __add__(s, o):
        o = o if isinstance(o, RF) else RF(o if not isinstance(o, int) else ctx.constant(o))
        return RF(s.n * o.d + o.n * s.d, s.d * o.d)
    __radd__ = __add__

    def __sub__(s, o):
        o = o if isinstance(o, RF) else RF(o if not isinstance(o, int) else ctx.constant(o))
        return RF(s.n * o.d - o.n * s.d, s.d * o.d)

    def __mul__(s, o):
        o = o if isinstance(o, RF) else RF(o if not isinstance(o, int) else ctx.constant(o))
        return RF(s.n * o.n, s.d * o.d)
    __rmul__ = __mul__

    def __truediv__(s, o):
        o = o if isinstance(o, RF) else RF(o if not isinstance(o, int) else ctx.constant(o))
        return RF(s.n * o.d, s.d * o.n)


def allpos(p):
    return len(p) > 0 and all(int(c) > 0 for c in p.to_dict().values())


def theta(X):
    return RF(X, 1 + X)


# ---------------- chart of the last range of Theorem 2 ----------------
# 1 < j < b, 0 < ab < 1, old-H failure p(j-1) < 1, e >= L^2/j (L = a+b),
# stationary hard part 0 < c < 1/(L+j), q >= 0.
# j = 1+T, b = j+A, a = theta(V)/b, p = theta(P)/T, e = L^2/j + U, c = theta(S)/(L+j), q = Q
# with A, P, T, V, S > 0 and U, Q >= 0.
j = RF(1 + T)
b = j + RF(A)
a = theta(V) / b
L = a + b
p = theta(P) / RF(T)
e = L * L / j + RF(U)
c = theta(S) / (L + j)
q = RF(Q)
dd = RF(one)  # d = 1
vals = [q, p, j, c, e, dd, b, a, RF(one)]   # order q p j c e d b a z

# ---------------- G (data) ----------------
data = json.load(open(GJSON))
assert data["variables"] == ["q", "p", "j", "c", "e", "d", "b", "a", "z"]
Gterms = [(tuple(int(x) for x in ex), int(co)) for ex, co in data["polynomials"]["G"]]
deg = [max(ex[i] for ex, _ in Gterms) for i in range(9)]
print("G terms", len(Gterms), "per-variable degrees", deg, "total degrees", sorted({sum(ex) for ex, _ in Gterms}))

for k, v in enumerate(vals):
    check(f"chart denominator {k} has positive coefficients", allpos(v.d))

# cleared G: Gtilde = sum co prod n_i^h d_i^(m_i-h) = G(chart) * prod d_i^m_i
pw = []
for v, m in zip(vals, deg):
    pw.append([v.n ** h * v.d ** (m - h) for h in range(m + 1)])
Gt = ctx.constant(0)
for ex, co in Gterms:
    t = ctx.constant(co)
    for i, h in enumerate(ex):
        if h:
            t *= pw[i][h]
        elif deg[i]:
            t *= pw[i][0]
    Gt += t
print("Gtilde terms", len(Gt), flush=True)

# guards
fpc = a * a - p * c
rr = (p - c) * j + (b - c) * (p + 1 + a) - p * a
fphi = rr * rr - RF(ctx.constant(4)) * a * c * (p - b) * (p + 1 + a)
check("denominator of a^2-pc has positive coefficients", allpos(fpc.d))
check("denominator of Phi has positive coefficients", allpos(fphi.d))

# multipliers: (weight, guard, exponent in A,P,T,U,V,S,Q); guard 0 = pc, 2 = Phi
CERT = [("64/3", 0, (0, 0, 1, 2, 5, 1, 2)), ("12", 0, (0, 0, 1, 2, 6, 1, 2)),
        ("96/7", 0, (0, 1, 0, 1, 6, 1, 2)), ("96/7", 0, (0, 1, 0, 2, 6, 1, 1)),
        ("96/7", 0, (0, 1, 1, 2, 6, 1, 2)),
        ("13/2", 2, (0, 0, 0, 1, 3, 0, 2)), ("6/7", 2, (0, 0, 0, 1, 4, 0, 2)),
        ("13/2", 2, (0, 0, 0, 2, 3, 0, 1)), ("6/7", 2, (0, 0, 0, 2, 4, 0, 1)),
        ("13/2", 2, (0, 0, 1, 2, 3, 0, 2)), ("6/7", 2, (0, 0, 1, 2, 4, 0, 2))]
D = 42
M = {0: ctx.constant(0), 2: ctx.constant(0)}
for w, i, ex in CERT:
    W = Fraction(w) * D
    assert W > 0 and W.denominator == 1
    M[i] += ctx.from_dict({ex: int(W)})
check("multiplier sizes 5 and 6, positive integer coefficients",
      len(M[0]) == 5 and len(M[2]) == 6 and allpos(M[0]) and allpos(M[2]))

R = D * Gt - M[0] * fpc.n - M[2] * fphi.n
check("remainder R has only positive integer coefficients", allpos(R), f"{len(R)} terms (document: 134706)")
check("remainder term count equals the document's 134706", len(R) == 134706)
# strict marker: a monomial with U = Q = 0 (A,P,T,V,S > 0 on the chart)
marks = [ex for ex in R.to_dict() if ex[3] == 0 and ex[6] == 0]
check("R has a monomial free of the nonnegative variables U, Q", len(marks) > 0, f"{len(marks)} such monomials")

# consistency: Gtilde / prod d^m equals G at a random rational chart point, guards likewise
rnd = random.Random(7)


def poly_eval(poly, pt):
    s = Fraction(0)
    for ex, co in poly.to_dict().items():
        t = Fraction(int(co))
        for x, h in zip(pt, ex):
            if h:
                t *= x ** int(h)
        s += t
    return s


for trial in range(2):
    pt = [Fraction(rnd.randint(1, 30), rnd.randint(1, 30)) for _ in range(7)]
    xv = [poly_eval(v.n, pt) / poly_eval(v.d, pt) for v in vals]
    Gdir = sum(co * eval_prod for co, eval_prod in
               ((co, __import__("math").prod([xv[i] ** ex[i] for i in range(9)])) for ex, co in Gterms))
    den = __import__("math").prod([poly_eval(v.d, pt) ** m for v, m in zip(vals, deg)])
    check(f"Gtilde/prod(den^m) == G(chart point) at random point {trial}", poly_eval(Gt, pt) / den == Gdir)
    av, bv, pv, cv, jv = xv[7], xv[6], xv[1], xv[3], xv[2]
    check(f"num/den of a^2-pc correct at point {trial}", poly_eval(fpc.n, pt) / poly_eval(fpc.d, pt) == av * av - pv * cv)
    rv = (pv - cv) * jv + (bv - cv) * (pv + 1 + av) - pv * av
    phiv = rv * rv - 4 * av * cv * (pv - bv) * (pv + 1 + av)
    check(f"num/den of Phi correct at point {trial}", poly_eval(fphi.n, pt) / poly_eval(fphi.d, pt) == phiv)
    # chart lands in the stated region
    ev_ = xv[4]
    ok = (jv > 1 and bv > jv and 0 < av * bv < 1 and pv * (jv - 1) < 1 and pv > 0
          and ev_ * jv >= (av + bv) ** 2 and 0 < cv < 1 / (av + bv + jv))
    check(f"chart point {trial} lies in the stated region", ok)

print("ALL EXACT CHECKS PASS" if not FAIL else f"FAILURES: {FAIL}")
sys.exit(1 if FAIL else 0)
