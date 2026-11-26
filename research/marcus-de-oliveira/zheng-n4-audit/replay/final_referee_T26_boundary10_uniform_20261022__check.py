"""Replay for final_referee_T26_boundary10_uniform_20261022 (faces T_27 u {10}, T_27 u {13}). Auditor code."""
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
T10 = T | {10}; T13 = T | {13}
G = L.generators(T10)
check("generators of T u {10} = {0,1,2,6,53,54,55,67,70}", G == [0, 1, 2, 6, 53, 54, 55, 67, 70], G)
check("non-minor generators are exactly 54, 70", [j for j in G if not L.F_IS_MINOR[j]] == [54, 70])
check("identity (2) f54 + f70 = f0 + f6 + f55 + f67", identity([54, 70], [0, 6, 55, 67]))
check("identity (3) f54 + f1 = f0 + f6 + f53", identity([54, 1], [0, 6, 53]))
g = stabilizer_has_action(T10, {0: 67, 67: 0, 1: 53, 53: 1, 6: 55, 55: 6, 54: 70, 70: 54, 2: 2})
check("a group element fixes T u {10} and acts 0<->67, 1<->53, 6<->55, 54<->70, 2 fixed", g is not None, g)
tr = (L.PERMS[0], L.PERMS[0], 1)
check("transposition (sigma -> sigma^-1) maps T u {10} onto T u {13}", L.g_act_set(tr, T10) == frozenset(T13))

t0 = time.time()
W = [(0, 'q'), (2, 'p'), (6, 'c'), (53, 'e'), (54, 'd'), (55, 'b'), (67, 'a')]
ctx, V, Cc, Qc, D, N, R = radial(W)
q, p, c, e, d, b, a = [V[k] for k in 'qpcedba']
print("INFO Hessian determinant: %d terms, degree %d in z, %.0f s" % (len(D.coeffs()), N, time.time() - t0))
Tt = (c + d) * (d + q); A0 = (c + d + q) * (a + d) * (b + d) + c * q * (a + b + d)
check("C = -(A0 + eT)", Cc == -(A0 + e * Tt))
check("Q = -d(d+e)", Qc == -d * (d + e))
G6 = A0 * (A0 - e * (a * b - Tt)) + p * (a + c + d) * (b + d + q) * (A0 + e * Tt)
check("G6 of (6) has 206 monomials", len(G6.coeffs()) == 206, len(G6.coeffs()))
LHS = R * Qc ** (4 - N)
check("(5) Q^4 det H(y0 - (C/Q) kappa) = -4 Q^2 C^2 G6^2", LHS == -4 * Qc ** 2 * Cc ** 2 * G6 ** 2)
K = a * b * (e - d) - d ** 2 * (a + b + e + d)
check("sign identity A0 - e(ab - T) = -K + (c+q)(a+d)(b+d) + cq(a+b+d) + e(d(c+q)+cq)",
      A0 - e * (a * b - Tt) == -K + (c + q) * (a + d) * (b + d) + c * q * (a + b + d) + e * (d * (c + q) + c * q))
check("A0 - d(a+d)(b+d) has nonnegative coefficients", nonneg_coeffs(A0 - d * (a + d) * (b + d)))
check("p(a+c+d)(b+d+q)(A0+eT) has nonnegative coefficients", nonneg_coeffs(p * (a + c + d) * (b + d + q) * (A0 + e * Tt)))
print("RESULT: %d exact checks failed" % len(FAIL) + ((": " + ", ".join(FAIL)) if FAIL else ""))
sys.exit(1 if FAIL else 0)
