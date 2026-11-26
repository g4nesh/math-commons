"""Replay for final_referee_T26_uniform_20261022 (paper case T_27). Auditor code; run from the replay directory:  python -I -B final_referee_T26_uniform_20261022__check.py"""
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

import numpy as np
T = {2, 3, 4, 7, 8, 12}
check("contact set equals paper T_27 exactly", frozenset(T) == L.CONTACT_SETS[27])
G = L.generators(T)
check("generator list F_T equals {0,1,2,5,6,53,54,55,67,70}", G == [0, 1, 2, 5, 6, 53, 54, 55, 67, 70], G)
check("non-minor generators are exactly 54, 70", [j for j in G if not L.F_IS_MINOR[j]] == [54, 70])
check("identity (3a) f54 + f70 = f0 + f6 + f55 + f67", identity([54, 70], [0, 6, 55, 67]))
check("identity (3b) f54 + f1 = f0 + f6 + f53", identity([54, 1], [0, 6, 53]))
g = stabilizer_has_action(T, {0: 67, 67: 0, 1: 53, 53: 1, 6: 55, 55: 6, 54: 70, 70: 54, 2: 2, 5: 5})
check("a group element fixes T and acts 0<->67, 1<->53, 6<->55, 54<->70, 2, 5 fixed", g is not None, g)

t0 = time.time()
W = [(0, 'q'), (2, 'p'), (5, 'j'), (6, 'c'), (53, 'e'), (54, 'd'), (55, 'b'), (67, 'a')]
ctx, V, Cc, Qc, D, N, R = radial(W)
q, p, j, c, e, d, b, a = [V[k] for k in 'qpjcedba']
print("INFO Hessian determinant: %d terms, degree %d in z, %.0f s" % (len(D.coeffs()), N, time.time() - t0))
Tt = (c + d) * (d + q); A0 = (c + d + q) * (a + d) * (b + d) + c * q * (a + b + d)
check("C(y0 + z kappa) = C + zQ with C = -(A0 + eT)", Cc == -(A0 + e * Tt))
check("Q = -d(d+e) - pj", Qc == -d * (d + e) - p * j)
I = c + d + q; Rr = p * j; m = p + j; U = (a + d) * (b + d); Vv = U + d * e; S = a + b + d + e
Wc = j * (a + d) + p * (b + d) + 2 * Rr; Wq = p * (a + d) + j * (b + d) + 2 * Rr
K = a * b * (e - d) - d ** 2 * (a + b + e + d)
al = I * (U - Rr) + c * q * (a + b + d)
be = -K + (m + c + q) * Vv + (Rr + c * q) * S
Hh = Vv * (c * (c + d) * Wc + q * (q + d) * Wq) + c * q * (c * q * m * S + c * (S * Wc + m * Vv) + q * (S * Wq + m * Vv) + m * (a + b + 2 * d) * Vv + 4 * Rr * d * S + a * b * (d + e) * m)
G6 = al * be + c * q * (d ** 2 - Rr) * m * S + Hh
check("compact G6 of (6) has 401 monomials and degree 6", len(G6.coeffs()) == 401 and G6.total_degree() == 6, (len(G6.coeffs()), G6.total_degree()))
LHS = R * Qc ** (4 - N)          # = Q^4 det H(y0 - (C/Q) kappa)
check("(5) Q^4 det H(y0 - (C/Q) kappa) = -4 Q^2 C^2 G6^2 with G6 from (6)", LHS == -4 * Qc ** 2 * Cc ** 2 * G6 ** 2)
check("all coefficients of H in (6) are nonnegative", nonneg_coeffs(Hh))
check("all coefficients of -C are nonnegative and d^3 occurs", nonneg_coeffs(-Cc) and int((-Cc).to_dict().get((0, 0, 0, 0, 0, 3, 0, 0, 0), 0)) == 1)

# NUMERIC: Lemma 1 sign claim on random points with pj < d^2, K < 0 (redundant with the proof)
vars_ = ['q', 'p', 'j', 'c', 'e', 'd', 'b', 'a']
import fractions
def ev(P, x):
    return float(sum(float(cf) * np.prod([x[n] ** int(ex) for n, ex in zip(vars_ + ['z'], mon)]) for mon, cf in zip(P.monoms(), P.coeffs())))
rng = np.random.default_rng(1); worst = np.inf; cnt = 0
while cnt < 2000:
    x = dict(zip(vars_, rng.exponential(1.0, 8) * rng.choice([0, 1, 1, 1], 8))); x['d'] = 1.0; x['z'] = 0.0
    if x['p'] * x['j'] >= 1 or x['a'] * x['b'] * (x['e'] - 1) - (x['a'] + x['b'] + x['e'] + 1) >= 0:
        continue
    cnt += 1
    worst = min(worst, ev(G6, x))
print("NUMERIC Lemma 1: 2000 random points (seed 1) with pj < d^2 = 1 and K < 0: min G6 = %.4g (> 0 expected)" % worst)
print("RESULT: %d exact checks failed" % len(FAIL) + ((": " + ", ".join(FAIL)) if FAIL else ""))
sys.exit(1 if FAIL else 0)
