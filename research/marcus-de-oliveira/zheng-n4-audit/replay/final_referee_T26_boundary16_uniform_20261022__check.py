"""Replay for final_referee_T26_boundary16_uniform_20261022 (face T_27 u {16}). Auditor code."""
# ---------------------------------------------------------------------------
# Own code (auditor). Rebuilds everything from the paper's definitions via moc4lib:
# generators f_j (Appendix D), the cubic C of (A.2), kappa = e_13.
# Nothing from the supplement is executed or imported.
# ---------------------------------------------------------------------------
import sys, os, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L
import flint

FAIL = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + ((" : " + str(info)) if info != "" else ""))
    if not ok:
        FAIL.append(name)

def fvals_sum(d):
    s = [0] * 24
    for j, c in d.items():
        for t in range(24):
            s[t] += c * L.F_VALUES[j][t]
    return s

def identity(lhs, rhs):
    return fvals_sum({j: 1 for j in lhs}) == fvals_sum({j: 1 for j in rhs})

def f_action(g, js):
    out = {}
    for j in js:
        v = tuple(L.g_act_values(g, L.F_VALUES[j]))
        ks = [k for k in range(82) if tuple(L.F_VALUES[k]) == v]
        out[j] = ks[0] if ks else None
    return out

def stabilizer_has_action(T, claimed):
    """claimed: dict j->k on the generators; True if some g in G fixes T and acts so."""
    G = L.generators(T)
    for g in L.G_ELEMENTS:
        if L.g_act_set(g, T) != frozenset(T):
            continue
        m = f_action(g, G)
        if all(m[j] == k for j, k in claimed.items()):
            return g
    return None

def _split_z(P, ctx):
    out = {}
    for mon, c in zip(P.monoms(), P.coeffs()):
        e = mon[-1]; mm = list(mon); mm[-1] = 0
        out[e] = out.get(e, ctx.from_dict({})) + ctx.from_dict({tuple(mm): c})
    return out

def bareiss_det(M):
    M = [r[:] for r in M]
    n = len(M); sign = 1; prev = None
    for k in range(n - 1):
        piv = None
        for i in range(k, n):
            if not M[i][k].is_zero():
                if piv is None or len(M[i][k].coeffs()) < len(M[piv][k].coeffs()):
                    piv = i
        if piv is None:
            return M[0][0] * 0
        if piv != k:
            M[k], M[piv] = M[piv], M[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                t = M[i][j] * M[k][k] - M[i][k] * M[k][j]
                if prev is not None:
                    q, r = divmod(t, prev)
                    assert r.is_zero()
                    t = q
                M[i][j] = t
            M[i][k] = M[i][k] * 0
        prev = M[k][k]
    return sign * M[n - 1][n - 1]

def radial(weights):
    """weights: list of (f index, variable name). y = sum w_j f_j + z*kappa (matching
    coordinates). Returns ctx, vars, C, Q (with C(y0+z kappa) = C + z Q checked),
    D(z) = det Hess C (y0 + z kappa), N = deg_z D, R = Q^N D(-C/Q)."""
    names = [v for _, v in weights] + ["z"]
    ctx = flint.fmpz_mpoly_ctx.get(tuple(names), 'lex')
    V = dict(zip(names, ctx.gens()))
    y = [ctx.from_dict({}) for _ in range(14)]
    for (j, v) in weights:
        for k in range(14):
            fk = L.F_VECS[j][k]
            assert fk.denominator == 1
            if fk != 0:
                y[k] = y[k] + int(fk) * V[v]
    y[13] = y[13] + V["z"]                       # kappa = e_13 (paper (3.2))
    Cz = L.cubic_C(y)
    parts = _split_z(Cz, ctx)
    assert max(parts) <= 1, "C(y0 + z kappa) is not affine in z"
    Cc = parts.get(0, ctx.from_dict({})); Qc = parts.get(1, ctx.from_dict({}))
    mctx = flint.fmpz_mpoly_ctx.get(tuple("m%d" % i for i in range(14)), 'lex')
    Cm = L.cubic_C(list(mctx.gens()))
    H = []
    for k in range(14):
        dk = Cm.derivative(k); row = []
        for l in range(14):
            dkl = dk.derivative(l); e = ctx.from_dict({})
            for mon, c in zip(dkl.monoms(), dkl.coeffs()):
                idx = [i for i, ex in enumerate(mon) if ex][0]
                e = e + int(c) * y[idx]
            row.append(e)
        H.append(row)
    D = bareiss_det(H)
    Dk = _split_z(D, ctx)
    N = max(Dk)
    R = ctx.from_dict({})
    for k, dk in Dk.items():
        R = R + dk * (-Cc) ** k * Qc ** (N - k)
    return ctx, V, Cc, Qc, D, N, R

def nonneg_coeffs(P):
    return all(int(c) >= 0 for c in P.coeffs())

T = {2, 3, 4, 7, 8, 12}
T16 = T | {16}
G = L.generators(T16)
check("generators of T u {16} = {0,2,5,6,54,55,67,70}", G == [0, 2, 5, 6, 54, 55, 67, 70], G)
check("non-minor generators are exactly 54, 70", [j for j in G if not L.F_IS_MINOR[j]] == [54, 70])
check("identity (1) f54 + f70 = f0 + f6 + f55 + f67", identity([54, 70], [0, 6, 55, 67]))
g = stabilizer_has_action(T16, {54: 70, 70: 54, 2: 5, 5: 2, 55: 67, 67: 55, 0: 0, 6: 6})
check("a group element fixes T u {16} and acts 54<->70, 2<->5, 55<->67, 0, 6 fixed", g is not None, g)

t0 = time.time()
W = [(0, 'q'), (2, 'a'), (5, 'b'), (6, 'c'), (54, 'd'), (55, 'h'), (67, 'w')]
ctx, V, Cc, Qc, D, N, R = radial(W)
q, a, b, c, d, h, w = [V[k] for k in 'qabcdhw']
print("INFO Hessian determinant: %d terms, degree %d in z, %.0f s" % (len(D.coeffs()), N, time.time() - t0))
check("Q = -ab - d^2", Qc == -a * b - d ** 2)
check("C = -(d+h)(d+q)(d+w) - c(d^2 + d(h+q+w) + hq + hw + qw)",
      Cc == -(d + h) * (d + q) * (d + w) - c * (d ** 2 + d * (h + q + w) + h * q + h * w + q * w))
LHS = R * Qc ** (4 - N)
quo, rem = divmod(LHS, -4 * Qc ** 2 * Cc ** 2)
check("Q^4 det H(...) is divisible by -4 Q^2 C^2", rem.is_zero())
G6 = quo.sqrt()
check("quotient is a perfect square G6^2 (3)", G6 is not None)
zero = {'q': 0, 'c': 0, 'h': 0, 'w': 0}
G6r = G6.subs(zero)
target = -d ** 2 * (a + d) * (b + d) * (a * b - d ** 2)
if G6r == -target:
    G6 = -G6; G6r = -G6r
check("(4) G6 at q=c=h=w=0 equals -d^2 (a+d)(b+d)(ab-d^2) (sign of the square root fixed by this)", G6r == target)
check("G6 has 277 monomials and degree 6", len(G6.coeffs()) == 277 and G6.total_degree() == 6, (len(G6.coeffs()), G6.total_degree()))
Ga0 = G6.subs({'a': 0}); Gb0 = G6.subs({'b': 0})
check("G6|a=0 has 136 coefficients, all positive", len(Ga0.coeffs()) == 136 and all(int(x) > 0 for x in Ga0.coeffs()), len(Ga0.coeffs()))
check("G6|b=0 has 136 coefficients, all positive", len(Gb0.coeffs()) == 136 and all(int(x) > 0 for x in Gb0.coeffs()), len(Gb0.coeffs()))
# (5): a^2 (1+u)^2 G6(q, a, d^2/(a(1+u)), c, d, h, w)
ctx2 = flint.fmpz_mpoly_ctx.get(('q', 'a', 'u', 'c', 'd', 'h', 'w'), 'lex')
q2, a2, u2, c2, d2, h2, w2 = ctx2.gens()
degb = G6.degrees()[2]
check("deg_b G6 <= 2", degb <= 2, degb)
tot = ctx2.from_dict({})
for mon, cf in zip(G6.monoms(), G6.coeffs()):
    eq, ea, eb, ec, ed, eh, ew, ez = mon
    term = int(cf) * q2 ** eq * a2 ** ea * c2 ** ec * d2 ** ed * h2 ** eh * w2 ** ew
    term = term * (d2 ** 2) ** eb * (a2 * (1 + u2)) ** (2 - eb)
    tot = tot + term
check("(5) cleared polynomial has 525 monomials, all coefficients positive", len(tot.coeffs()) == 525 and all(int(x) > 0 for x in tot.coeffs()), len(tot.coeffs()))
print("RESULT: %d exact checks failed" % len(FAIL) + ((": " + ", ".join(FAIL)) if FAIL else ""))
sys.exit(1 if FAIL else 0)
