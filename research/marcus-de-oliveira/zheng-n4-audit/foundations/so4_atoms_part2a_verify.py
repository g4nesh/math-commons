"""so4_atoms_part2a_verify.py -- independent verification (label atoms2a-verify) of the part-2a SO(4) atoms:
(R) rank-two family, (F1) T24 coherent/parallel/double-SOC, (F2) T25 partial-four, (F3) T41 bad17, (F4) inventory56.

Sections:  R  own proof of the rank-two inequality (exact identities, sympy/flint)
           G  own re-verification of the Gram certificates from the supplement's cubic features (JSON read as DATA,
              parsed by a small own tokenizer; nothing from the supplement is executed or imported)
           T  orbits / transports / roles under the group G of order 1152 (shared lib moc4lib)
           N  numeric attack (NUMERIC evidence only; deterministic seeds)
Run:  MOC4_SUPP=<supplement root> python -I -B so4_atoms_part2a_verify.py   (lib located at ../lib; override MOC4_LIB)
"""
import os, sys, json, re, time, math
import numpy as np
from scipy.optimize import minimize
import mpmath as mp
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.environ.get("MOC4_LIB", os.path.join(HERE, "..", "lib"))
sys.path.insert(0, os.path.abspath(LIB))
import moc4lib as L
import flint
import sympy as sp

SUPP = os.environ.get("MOC4_SUPP") or (sys.argv[1] if len(sys.argv) > 1 else None)
BASE = os.path.join(SUPP, "proofs", "manin_transfer_20261021") if SUPP else None

NPASS = NFAIL = 0
def check(name, ok, info=""):
    global NPASS, NFAIL
    if ok: NPASS += 1
    else: NFAIL += 1
    print(("PASS " if ok else "FAIL ") + name + ((" : " + str(info)) if info != "" else ""), flush=True)

def fq(x): return flint.fmpq(x.numerator, x.denominator)

# ---------------------------------------------------------------- quaternions over any ring
def qmul(a, b):
    a0, a1, a2, a3 = a; b0, b1, b2, b3 = b
    return [a0*b0 - a1*b1 - a2*b2 - a3*b3, a0*b1 + a1*b0 + a2*b3 - a3*b2,
            a0*b2 - a1*b3 + a2*b0 + a3*b1, a0*b3 + a1*b2 - a2*b1 + a3*b0]
def qconj(a): return [a[0], -a[1], -a[2], -a[3]]
def omat(p, q, conj_q, unit):
    """4x4 matrix of x -> p x q  (conj_q False) or x -> p x conj(q) (conj_q True); columns = images of 1,i,j,k."""
    qq = qconj(q) if conj_q else q
    cols = []
    for j in range(4):
        e = [unit * 0] * 4; e[j] = unit
        cols.append(qmul(qmul(p, e), qq))
    return [[cols[j][i] for j in range(4)] for i in range(4)]

def Phi_poly(j, O):
    w = L.leibniz_weights(O)
    tot = 0
    for v, x in zip(L.F_VALUES[j], w):
        if v: tot = tot + v * x
    return tot

# ================================================================ section R: own proof of the rank-two family
def section_R():
    print("\n=== R: rank-two family, own proof (identities EXACT) ===")
    ctx = flint.fmpq_mpoly_ctx.get(('u0', 'u1', 'u2', 'u3', 'v0', 'v1', 'v2', 'v3'), 'lex')
    g = ctx.gens(); u = list(g[:4]); v = list(g[4:]); one = ctx.from_dict({(0,)*8: 1})
    c1 = [one, one, 0*one, 0*one]; c2 = [0*one, 0*one, one, one]          # c = 1+i, c' = j+k (unnormalised)
    p = qmul(u, qconj(c1)); q = qmul(v, qconj(c2))                         # |p|^2 = 2|u|^2, |q|^2 = 2|v|^2
    O = omat(p, q, True, one)
    O = [[x * flint.fmpq(1, 2) for x in r] for r in O]                     # Ohat = O(p,q)/2, Ohat^T Ohat = |u|^2|v|^2 I
    nu = sum(x*x for x in u); nv = sum(x*x for x in v)
    OtO = [[sum(O[k][i]*O[k][j] for k in range(4)) for j in range(4)] for i in range(4)]
    check("R1 Ohat^T Ohat = |u|^2|v|^2 I", all(OtO[i][j] == (nu*nv if i == j else 0*one) for i in range(4) for j in range(4)))
    check("R1 det Ohat = (|u|^2|v|^2)^2 (so Ohat in SO(4) on S^3 x S^3)", L.det_leibniz4(O) == (nu*nv)**2)
    print("   names:", {j: L.f_name(j) for j in (2, 12, 29, 35)})
    P35 = Phi_poly(35, O); P12 = Phi_poly(12, O); P2 = Phi_poly(2, O); P29 = Phi_poly(29, O)
    w2 = sum((u[i]*v[i])**2 for i in range(4))
    rhs35 = 2*w2*nu*nv + 4*(u[0]*u[1]*u[2]*u[3]*nv*nv + v[0]*v[1]*v[2]*v[3]*nu*nu)
    check("R1 Phi35 = 2 sum u_i^2 v_i^2 + 4(u0u1u2u3 + v0v1v2v3)  (homogenised)", P35 == rhs35)
    x12 = L.minor_det(O, (1, 4), (2, 4))
    X33 = u[0]**2 - u[1]**2 - u[2]**2 + u[3]**2; Y33 = v[0]**2 - v[1]**2 - v[2]**2 + v[3]**2
    xx = (X33*nv + Y33*nu) * flint.fmpq(1, 2)
    sgn = 1 if x12 == xx else (-1 if x12 == -xx else 0)
    check("R1 x12 := det Ohat[14,24] = +-(X33 + Y33)/2, X33 = u0^2-u1^2-u2^2+u3^2", sgn != 0, "sign %d" % sgn)
    check("R1 Phi12 = x12^2", P12 == x12*x12)
    S = 2*((u[1]*u[2] + u[0]*u[3])**2 * nv*nv + (v[1]*v[2] + v[0]*v[3])**2 * nu*nu)
    check("R1 Phi2 + Phi29 = 2[(u1u2+u0u3)^2 + (v1v2+v0v3)^2]", P2 + P29 == S)
    # --- elementary algebra of the proof (sympy)
    al, alp, c1_, c2_, c1p, c2p, h, t = sp.symbols('alpha alphap c1 c2 c1p c2p h t', real=True)
    be, bep = 1 - al, 1 - alp
    y = al*c1_ - be*c2_; z = alp*c1p - bep*c2p; dl = alp - al; x = al + alp - 1
    Qd = al*alp*(c1_**2 + c1p**2)/2 + be*bep*(c2_**2 + c2p**2)/2 - al*be*c1_*c2_ - alp*bep*c1p*c2p
    K = al*c1_**2 - be*c2_**2 - alp*c1p**2 + bep*c2p**2
    check("R2 Qd = (y^2+z^2)/2 + (delta/2) K", sp.expand(Qd - (y**2 + z**2)/2 - dl*K/2) == 0)
    fu = (y - h)**2 + dl*(al*c1_**2 - be*c2_**2)
    fp = (z + h)**2 - dl*(alp*c1p**2 - bep*c2p**2)
    Lf = Qd + h*h + h*(z - y)          # = Qd + x^2/4 + |x|(z-y)/2 with h = |x|/2
    check("R2 Qd + h^2 + h(z-y) = (f_u + f_p)/2", sp.expand(Lf - (fu + fp)/2) == 0)
    k = be*c2_ + h
    check("R2 f_u = al*alp*(c1 - k/alp)^2 + (delta/alp)(k^2 - alp*be*c2^2)",
          sp.simplify(sp.expand(fu - (al*alp*(c1_ - k/alp)**2 + dl/alp*(k**2 - alp*be*c2_**2)))) == 0)
    sub = {al: 1 - alp, alp: 1 - al, c1_: c2p, c2_: c1p}
    fu_sw = fu.subs(sub, simultaneous=True)
    check("R2 f_p(al,alp,c1',c2') = f_u(1-alp, 1-al, c2', c1')", sp.expand(fp - fu_sw) == 0)
    xs = sp.expand(x.subs(sub, simultaneous=True))
    check("R2 x is invariant under (al,alp) -> (1-alp,1-al) up to sign", sp.expand(xs + x) == 0)
    sA, sB = sp.symbols('sA sB', positive=True)   # sA = sqrt(alp), sB = sqrt(beta)
    check("R2 (alp-beta)/2 - sqrt(beta)(sqrt(alp)-sqrt(beta)) = (sqrt(alp)-sqrt(beta))^2/2",
          sp.expand((sA**2 - sB**2)/2 - sB*(sA - sB) - (sA - sB)**2/2) == 0)
    a0, a3, A0, A3 = sp.symbols('a0 a3 A0 A3')
    check("R2 a0A0 + a3A3 = [(a0+a3)(A0+A3) + (a0-a3)(A0-A3)]/2",
          sp.expand(a0*A0 + a3*A3 - ((a0 + a3)*(A0 + A3) + (a0 - a3)*(A0 - A3))/2) == 0)
    s1, s1p = sp.symbols('s1 s1p')
    check("R2 2 - 2 s s' - c^2 - c'^2 = (s - s')^2 when c^2+s^2 = c'^2+s'^2 = 1",
          sp.expand((2 - 2*s1*s1p - (1 - s1**2) - (1 - s1p**2)) - (s1 - s1p)**2) == 0)
    # infimum over t of t S + X/(4t) = sqrt(S X):  AM-GM, and the coefficient form
    d_, xm, S_, X_ = sp.symbols('d x S X', positive=True)
    tt = sp.sqrt(X_/(4*S_))
    check("R2 at t* = sqrt(X/(4S)): t S + X/(4t) = sqrt(S X)", sp.simplify(tt*S_ + X_/(4*tt) - sp.sqrt(S_*X_)) == 0)
    check("R2 x(S) + (d^2/(4x)) X - d sqrt(S X) = (sqrt(x S) - d sqrt(X)/(2 sqrt x))^2",
          sp.simplify(xm*S_ + d_**2/(4*xm)*X_ - d_*sp.sqrt(S_*X_) - (sp.sqrt(xm*S_) - d_*sp.sqrt(X_)/(2*sp.sqrt(xm)))**2) == 0)

# ================================================================ data parsing (no eval)
TOK = re.compile(r'\s*(\*\*|\*|/|\+|-|[pq][0-3]|\d+)')
def parse_poly(s, ctx):
    gens = dict(zip(['p0', 'p1', 'p2', 'p3', 'q0', 'q1', 'q2', 'q3'], ctx.gens()))
    s = s.strip(); pos = 0; toks = []
    while pos < len(s):
        m = TOK.match(s, pos)
        if not m: raise ValueError("bad token in data: %r" % s[pos:pos+10])
        toks.append(m.group(1)); pos = m.end()
    one = ctx.from_dict({(0,)*8: 1})
    i = 0; total = 0*one
    while i < len(toks):
        sign = 1
        while toks[i] in '+-':
            if toks[i] == '-': sign = -sign
            i += 1
        term = one * sign; num, den = 1, 1
        while True:
            t = toks[i]; i += 1
            if t in gens: base = gens[t]
            elif t.isdigit(): base = one * int(t)
            else: raise ValueError(t)
            if i < len(toks) and toks[i] == '**':
                e = int(toks[i+1]); i += 2; base = base**e
            term = term * base
            if i < len(toks) and toks[i] == '*': i += 1; continue
            if i < len(toks) and toks[i] == '/':
                term = term * flint.fmpq(1, int(toks[i+1])); i += 2
                if i < len(toks) and toks[i] == '*': i += 1; continue
            break
        total = total + term
    return total

def frac(s): return Fr(s)

# ================================================================ section G: Gram certificates
class Face:
    def __init__(self, name, cubfile, budfile, contacts_expected):
        self.name = name
        self.C = json.load(open(os.path.join(BASE, cubfile)))
        self.D = json.load(open(os.path.join(BASE, budfile)))
        self.src = list(self.D['source_order'])
        self.m = len(self.src)
        self.contacts = list(self.C['contacts'])
        assert sorted(self.contacts) == sorted(contacts_expected)

def reduce_N(f, ctx, s_pows):
    """normal form mod N = |p|^2 - |q|^2: p0^k -> p0^(k mod 2) s^(k div 2), s = |q|^2 - p1^2 - p2^2 - p3^2."""
    groups = {}
    for mon, c in f.terms():
        k = mon[0]; rest = (k % 2,) + tuple(mon[1:])
        groups.setdefault(k // 2, {})[rest] = groups.setdefault(k // 2, {}).get(rest, 0) + c
    out = f * 0
    for e, d in groups.items():
        out = out + ctx.from_dict(d) * s_pows(e)
    return out

def setup_face(F):
    ctx = flint.fmpq_mpoly_ctx.get(('p0', 'p1', 'p2', 'p3', 'q0', 'q1', 'q2', 'q3'), 'lex')
    g = ctx.gens(); one = ctx.from_dict({(0,)*8: 1})
    P = list(g[:4]); Q = list(g[4:])
    s = sum(x*x for x in Q) - P[1]**2 - P[2]**2 - P[3]**2
    cache = {0: one}
    def s_pows(e):
        if e not in cache: cache[e] = s_pows(e - 1) * s
        return cache[e]
    red = lambda f: reduce_N(f, ctx, s_pows)
    np2 = sum(x*x for x in P)
    O = omat(P, Q, False, one)                       # x -> p x q (the data's model)
    Phis = [Phi_poly(j, O) for j in F.src]
    feats = [[parse_poly(sx, ctx) for sx in par] for par in F.C['features']]
    maps = F.C['left_four']
    T = [[frac(x) for x in row] for row in F.D['T']]
    return dict(ctx=ctx, one=one, red=red, np2=np2, O=O, Phis=Phis, feats=feats, maps=maps, T=T, g=g)

def gram_M(S, conv):
    """M_ik = sum_parity sum_j g_i(S_j y) g_k(S_j y), g = conv(T) F."""
    ctx, g = S['ctx'], S['g']
    T = S['T']; m = len(T)
    if conv == 'Tinv': A = L.inv(T)
    elif conv == 'T': A = T
    elif conv == 'TinvT': A = L.transpose(L.inv(T))
    else: A = L.transpose(T)
    M = [[0*S['one'] for _ in range(m)] for _ in range(m)]
    for par in S['feats']:
        gf = [sum((fq(A[i][k]) * par[k] for k in range(m) if A[i][k] != 0), 0*S['one']) for i in range(m)]
        for mp in S['maps']:
            imgs = [g[idx] * sg for idx, sg in mp]
            gs = [f.compose(*imgs) for f in gf]
            for i in range(m):
                for k in range(i, m):
                    M[i][k] = M[i][k] + gs[i] * gs[k]
    for i in range(m):
        for k in range(i):
            M[i][k] = M[k][i]
    return M

def binv_cols(F):
    """cost coefficient of source l as a linear form in the Gram diagonal: returns m x m matrix B with c = B diag."""
    m = F.m
    bd = F.D.get('budget_diagonal')
    if bd is None or isinstance(bd[0], list):
        Bm = [[Fr(1) if i == j else Fr(0) for j in range(m)] for i in range(m)]
        if bd is not None:      # explicit matrix diag = bd * c  -> c = bd^{-1} diag
            Bm = L.inv([[frac(x) for x in r] for r in bd])
        return Bm
    names = [str(x) for x in bd]          # T24: (q,p,j,c,e,q+d,b,a)
    syms = {nm: sp.Symbol(nm) for nm in ['q', 'p', 'j', 'c', 'e', 'd', 'b', 'a']}
    exprs = [sp.sympify(nm, locals=syms) for nm in names]
    order = ['q', 'p', 'j', 'c', 'e', 'd', 'b', 'a']
    A = [[Fr(str(sp.diff(ex, syms[o]))) for o in order] for ex in exprs]     # diag = A c
    return L.inv(A)

def delta_matrix(F, S, M, kappa):
    """columns: entries (i<=k); value: normal form of Delta(E_ik) as dict monomial -> coeff."""
    m = F.m; Bm = binv_cols(F); red = S['red']; np2 = S['np2']
    cols = []; ents = [(i, k) for i in range(m) for k in range(i, m)]
    for (i, k) in ents:
        if i == k:
            f = np2 * M[i][i]
            for l in range(m):
                if Bm[l][i] != 0: f = f - fq(kappa * Bm[l][i]) * S['Phis'][l]
        else:
            f = 2 * np2 * M[i][k]
        cols.append(red(f).to_dict())
    return ents, cols

def mat_from_cols(cols):
    mons = sorted(set().union(*[c.keys() for c in cols]))
    idx = {mm: r for r, mm in enumerate(mons)}
    rows = [[Fr(0)] * len(cols) for _ in mons]
    for j, c in enumerate(cols):
        for mm, v in c.items():
            rows[idx[mm]][j] = Fr(int(v.p), int(v.q))
    return rows

def find_kappa(F, S, M, Gtest):
    """Phi_cost proportionality constant: |p|^2 sum G M = kappa * Phi_cost (mod N) for a test Gram."""
    m = F.m; Bm = binv_cols(F); red = S['red']
    diag = [Gtest[i][i] for i in range(m)]
    cost = [sum(Bm[l][i] * diag[i] for i in range(m)) for l in range(m)]
    lhs = red(S['np2'] * sum((fq(Fr(Gtest[i][k])) * M[i][k] for i in range(m) for k in range(m) if Gtest[i][k] != 0), 0*S['one']))
    rhs = red(sum((fq(cost[l]) * S['Phis'][l] for l in range(m) if cost[l] != 0), 0*S['one']))
    dl, dr = lhs.to_dict(), rhs.to_dict()
    if not dr: return None
    mm = next(iter(dr)); kap = Fr(int(dl.get(mm, 0).p), int(dl.get(mm, 0).q)) / Fr(int(dr[mm].p), int(dr[mm].q)) if mm in dl else None
    if kap is None: return None
    ok = (lhs - fq(kap) * rhs) == 0
    return kap if ok else None

def sym_G(m, entries):
    G = sp.zeros(m, m)
    for (i, k), val in entries.items():
        G[i, k] = val; G[k, i] = val
    return G

def kernel_check(Rrows, ents, G):
    vec = [G[i, k] for (i, k) in ents]
    bad = []
    for r in Rrows:
        s = sp.simplify(sum(sp.Rational(c.numerator, c.denominator) * v for c, v in zip(r, vec) if c != 0))
        if s != 0: bad.append(s)
    return bad

def ldl_pivots(G, order):
    """symbolic LDL^T pivots of G[order, order] (no pivoting); returns pivots and residual."""
    A = sp.Matrix([[G[i, k] for k in order] for i in order])
    n = A.shape[0]; piv = []
    for c in range(n):
        d = sp.simplify(A[c, c]); piv.append(d)
        if d == 0:
            if any(sp.simplify(A[c, r]) != 0 for r in range(c + 1, n)): return piv, 'zero pivot with nonzero row'
            continue
        for r in range(c + 1, n):
            for s_ in range(c + 1, n):
                A[r, s_] = sp.simplify(A[r, s_] - A[r, c] * A[c, s_] / d)
    return piv, 'ok'

def e2e(F, S, M, kappa, G, subs):
    """end-to-end: |p|^2 sum G M - kappa Phi_cost == 0 mod N at a rational parameter point."""
    m = F.m; Bm = binv_cols(F); red = S['red']
    Gn = [[Fr(str(sp.nsimplify(G[i, k].subs(subs)))) for k in range(m)] for i in range(m)]
    diag = [Gn[i][i] for i in range(m)]
    cost = [sum(Bm[l][i] * diag[i] for i in range(m)) for l in range(m)]
    lhs = S['np2'] * sum((fq(Gn[i][k]) * M[i][k] for i in range(m) for k in range(m) if Gn[i][k] != 0), 0*S['one'])
    rhs = sum((fq(kappa * cost[l]) * S['Phis'][l] for l in range(m) if cost[l] != 0), 0*S['one'])
    return red(lhs - rhs) == 0, cost

def section_G():
    print("\n=== G: Gram certificates of (R), (F1), (F2), re-verified with own code ===")
    faces = {
        'inv56': Face('inv56', 'final_referee_inventory56_boundary17_cubic_reduction_20261023.json',
                      'final_referee_inventory56_boundary17_budget_PSD7_20261023.json', [0, 1, 5, 7, 9, 12, 17]),
        'T24': Face('T24', 'root_moc4_T24_cubic_reduction_20261022.json', 'root_moc4_T24_budget_PSD8_20261022.json', [3, 4, 7, 8, 10, 12]),
        'T25': Face('T25', 'final_referee_T25_cubic_reduction_20261022.json', 'final_referee_T25_budget_PSD8_20261022.json', [2, 3, 5, 7, 8, 12]),
    }
    printed = {'inv56': [0, 2, 5, 7, 12, 29, 35], 'T24': [0, 2, 59, 6, 53, 54, 55, 67], 'T25': [0, 1, 2, 5, 6, 39, 67, 70]}
    for nm, F in faces.items():
        gens = L.generators(F.contacts)
        nonmin = [j for j in F.src if not L.F_IS_MINOR[j]]
        check("G0 %s: data source order = printed order" % nm, F.src == printed[nm], F.src)
        check("G0 %s: Gram sources are generators of the contact set %s" % (nm, F.contacts), set(F.src) <= set(gens),
              "generators %s" % gens)
        print("   %s: non-minor Gram sources %s (%s); generators not on the Gram: %s" %
              (nm, nonmin, [L.f_name(j) for j in nonmin], sorted(set(gens) - set(F.src))))
    # constructions -------------------------------------------------------------
    x, d, j, e, a, b, r, s_, t, z, c, q, w = sp.symbols('x d j e a b r s t z c q w', positive=True)
    out = {}
    # (R): inv56, positions (Q,J,C,E,B,A,D) = 0..6
    vv = sp.Matrix([0, 2*x, 0, 0, d, -2*x, 0]); ww = sp.Matrix([0, 0, 0, 0, 1, 0, 2])
    G4 = vv * vv.T + d * x * ww * ww.T
    cons = [('inv56', 'R: G4 = v v^T + d x w w^T', G4, {x: sp.Rational(2, 5), d: sp.Rational(3, 2)},
             lambda cost: [sp.simplify(cc) for cc in cost] == [0, 4*x**2, 0, 0, sp.simplify(d*x + d**2), 4*x**2, 4*d*x], None)]
    # (F1) T24 positions (Q,P,J,C,E,D,B,A) = 0..7, d = 1
    J, E, D, B, A = 2, 4, 5, 6, 7
    G5 = sym_G(8, {(J, J): j, (E, E): e, (D, D): 1, (B, B): b, (A, A): a, (B, A): 1, (E, D): x, (E, B): x, (E, A): -x, (E, J): x - 1})
    G5j0 = sym_G(8, {(E, E): e, (D, D): 1, (B, B): b, (A, A): a, (B, A): 1, (E, D): 1, (E, B): 1, (E, A): -1})
    Gds = sym_G(8, {(J, J): j, (E, E): e, (D, D): 1, (B, B): b, (A, A): a, (B, A): 1, (E, J): -1})
    tau = (1 - r*s_) / (r - s_)
    G6 = sym_G(8, {(B, B): r**2, (A, A): s_**2, (B, A): r*s_, (J, B): r*tau, (J, A): s_*tau, (J, E): -1, (D, D): 1, (E, E): e, (J, J): j})
    cons += [('T24', '(5) coherent', G5, {j: sp.Rational(1, 3), e: sp.Rational(5, 2), a: 3, b: 2, x: sp.Rational(2, 7)}, None, [J, D, B, A, E]),
             ('T24', '(5) j = 0 member', G5j0, {e: 7, a: 3, b: 2}, None, [D, B, A, E]),
             ('T24', 'double-SOC (x = 0 member)', Gds, {j: 2, e: 1, a: 3, b: 1}, None, [J, E, D, B, A]),
             ('T24', '(6) parallel', G6, {r: 3, s_: sp.Rational(1, 5), e: 2, j: 7}, None, [E, D, B, J, A])]
    # (F2) T25 positions (A,E,P,J,B,C,Q,D) = 0..7, d = 1
    hb = b - t**2/a; rr = 1 - t
    def T25G(zv, kap):
        # basis (A, B, Q, C, D): Gbase = Gram(A,B,Q) (+) c (+) 1
        Gb = sp.zeros(5, 5)
        Gb[0, 0], Gb[0, 1], Gb[1, 1], Gb[1, 2], Gb[2, 2], Gb[3, 3], Gb[4, 4] = a, t, b, zv, q, c, 1
        Gb[1, 0], Gb[2, 1] = t, zv
        K = sp.zeros(8, 5)
        K[0, 0] = 1                                              # A
        K[1, 0], K[1, 1], K[1, 4] = t*(b + t)/(a*b - t**2), -t*(a + t)/(a*b - t**2), t     # E
        K[2, 0], K[2, 1], K[2, 2], K[2, 3] = -kap*t/a, kap, -kap, rr/c                      # P = kap V + (r/c) C, V = B - (t/a)A - Q
        K[4, 1] = 1; K[5, 3] = 1; K[6, 2] = 1; K[7, 4] = 1        # B, C, Q, D ; J = 0
        return (K * Gb * K.T).applyfunc(sp.simplify), Gb, K
    kap12 = (rr - z) / (hb + q - 2*z)
    G12, Gb12, K12 = T25G(z, kap12)
    G9, Gb9, K9 = T25G(rr, 0)
    q10 = w**2 * rr**2 / hb
    G10, Gb10, K10 = T25G(w*rr, ((rr - w*rr) / (hb + q - 2*w*rr)).subs(q, q10))
    G10 = G10.subs(q, q10).applyfunc(sp.simplify); Gb10 = Gb10.subs(q, q10)
    kk = sp.Symbol('kappa')
    Gk, _, _ = T25G(z, kk)
    pt = {a: 3, b: 2, c: sp.Rational(1, 2), t: sp.Rational(1, 3), q: 5, z: sp.Rational(1, 4), w: sp.Rational(1, 2)}
    cons += [('T25', '(12) general z', G12, pt, None, None),
             ('T25', '(9) z = r, P = (r/c)C', G9, pt, None, None),
             ('T25', '(10) z = w r, q = w^2 r^2/hb', G10, pt, None, None)]
    res = {}
    for nm, F in faces.items():
        t0 = time.time()
        S = setup_face(F)
        # choose the congruence convention by the data's own test: the PDF construction must lie in the kernel
        test = {'inv56': G4.subs({x: sp.Rational(2, 5), d: sp.Rational(3, 2)}),
                'T24': G5.subs({j: sp.Rational(1, 3), e: sp.Rational(5, 2), a: 3, b: 2, x: sp.Rational(2, 7)}),
                'T25': G9.subs(pt)}[nm]
        Gt = [[Fr(str(test[i, k])) for k in range(F.m)] for i in range(F.m)]
        chosen = None
        for conv in ('Tinv', 'T', 'TinvT', 'TT'):
            M = gram_M(S, conv)
            kap = find_kappa(F, S, M, Gt)
            print("   %s: convention g = %s F -> proportionality constant %s" % (nm, conv, kap))
            if kap is not None and kap > 0:
                chosen = (conv, M, kap); break
        check("G1 %s: some congruence convention gives |p|^2 sum G M = kappa Phi_cost (mod N), kappa > 0" % nm, chosen is not None,
              chosen[0] + ", kappa = %s" % chosen[2] if chosen else "")
        conv, M, kap = chosen
        # minors: diagonal of M vs Phi (sanity: M_ii is a square sum, hence >= 0)
        ents, cols = delta_matrix(F, S, M, kap)
        Mat = mat_from_cols(cols)
        Rrows = L.rref(Mat)[0]
        Rrows = [rw for rw in Rrows if any(v != 0 for v in rw)]
        print("   %s: %d normal-form monomial rows, relation rank %d (%.1fs)" % (nm, len(Mat), len(Rrows), time.time() - t0))
        check("G1 %s: relation space has dimension 10" % nm, len(Rrows) == 10)
        # cross-check with the data's R matrix (data consistency only)
        Rd = [[frac(v) for v in row] for row in F.D['R']]
        both = L.rank(Rrows + Rd)
        check("G1 %s: own relation space = row space of the data's R" % nm, both == len(Rrows) == L.rank(Rd))
        res[nm] = (F, S, M, kap, ents, Rrows)
    # printed relations of rank-two p.2 (PDF) for inv56
    F, S, M, kap, ents, Rrows = res['inv56']
    pos = {'Q': 0, 'J': 1, 'C': 2, 'E': 3, 'B': 4, 'A': 5, 'D': 6}
    def rel(terms):
        vec = [Fr(0)] * len(ents)
        for coef, u1, u2 in terms:
            i, k = sorted((pos[u1], pos[u2]))
            vec[ents.index((i, k))] += coef
        return vec
    printed_rel = [rel([(1, 'Q', 'J'), (1, 'Q', 'B'), (1, 'J', 'B'), (-1, 'B', 'D')]),
                   rel([(1, 'Q', 'C'), (-1, 'Q', 'B'), (1, 'C', 'A'), (-1, 'B', 'A'), (1, 'B', 'D'), (-1, 'D', 'D')]),
                   rel([(1, 'Q', 'E'), (-1, 'Q', 'B'), (1, 'E', 'A'), (-1, 'B', 'A'), (1, 'B', 'D'), (-1, 'D', 'D')]),
                   rel([(1, 'Q', 'D')]), rel([(1, 'J', 'D')]), rel([(1, 'E', 'D')]), rel([(1, 'A', 'D')]),
                   rel([(1, 'J', 'C'), (-1, 'J', 'B'), (-1, 'C', 'A'), (-1, 'C', 'D'), (1, 'B', 'A'), (1, 'D', 'D')]),
                   rel([(1, 'C', 'E'), (1, 'C', 'A'), (1, 'C', 'D'), (-1, 'E', 'A')]),
                   rel([(1, 'E', 'B'), (-1, 'E', 'A'), (1, 'B', 'A'), (1, 'B', 'D')])]
    check("G1 inv56: the ten relations printed on rank-two p.2 are independent", L.rank(printed_rel) == 10)
    check("G1 inv56: printed relations span exactly the own relation space", L.rank(printed_rel + Rrows) == 10)
    # constructions
    for (nm, label, G, pt_, costcheck, order) in cons:
        F, S, M, kap, ents, Rrows = res[nm]
        bad = kernel_check(Rrows, ents, G)
        check("G2 %s %s: G satisfies the relations identically in its parameters" % (nm, label), not bad, bad[:2])
        ok, cost = e2e(F, S, M, kap, G, pt_)
        check("G2 %s %s: end-to-end identity |p|^2 sum G M = kappa Phi_cost mod N at a rational point" % (nm, label), ok, "cost %s" % cost)
        Bm = binv_cols(F)
        sc = [sp.simplify(sum(sp.Rational(Bm[l][i].numerator, Bm[l][i].denominator) * G[i, i] for i in range(F.m))) for l in range(F.m)]
        print("   cost coefficients on %s: %s" % (F.src, sc))
        out[(nm, label)] = (G, sc)
        if costcheck is not None:
            check("G2 %s %s: diagonal/cost as printed" % (nm, label), costcheck(sc))
        if order is not None:
            piv, st = ldl_pivots(G, order)
            print("   LDL pivots on order %s: %s (%s)" % (order, [sp.factor(pv) for pv in piv], st))
    # T25: kappa forced
    F, S, M, kap, ents, Rrows = res['T25']
    bad = kernel_check(Rrows, ents, Gk)
    sol = sp.solve([bb for bb in bad if bb != 0], kk, dict=True)
    print("   T25 with P = kappa V + (r/c) C: failing relations", len(bad), "solutions for kappa:", sol)
    check("G3 T25: exactly one relation constrains kappa, and it forces the PDF choice (12)",
          len(set(sp.simplify(bb) for bb in bad)) >= 1 and len(sol) == 1 and sp.simplify(sol[0][kk] - kap12) == 0)
    # ---- PSD region formulas (own symbolic derivation)
    Lc = (a + 1)*(b + 1)/(a*b - 1)
    piv, st = ldl_pivots(G5, [J, D, B, A, E])
    last = sp.simplify(piv[-1].subs(x, 1/(1 + j*Lc)))
    check("G4 (5): last LDL pivot at x = 1/(1+jL) equals e - 1/(j + 1/L)", sp.simplify(last - (e - 1/(j + 1/Lc))) == 0)
    check("G4 (5): last pivot minus [(1-x)^2/j + L x^2] form", sp.simplify(piv[-1] - (e - ((1 - x)**2/j + Lc*x**2))) == 0)
    check("G4 (5): e >= 1/(j+1/L)  <=>  j >= 1/e - 1/L (e, j, L > 0): difference of the two sides has the sign of e(jL+1) - L",
          sp.simplify((e - 1/(j + 1/Lc)) * (j*Lc + 1) - (e*(j*Lc + 1) - Lc)) == 0 and
          sp.simplify((j - (1/e - 1/Lc)) * e*Lc - (e*(j*Lc + 1) - Lc)) == 0)
    dd = sp.Symbol('dd', positive=True)
    pr5 = dd**2/e - dd*(a*b - dd**2)/((a + dd)*(b + dd))
    check("G4 (5): printed homogeneous threshold = d * normalised threshold at (a,b,e)/d",
          sp.simplify(pr5 - dd*(1/(e/dd) - 1/Lc.subs({a: a/dd, b: b/dd}, simultaneous=True))) == 0)
    piv6, st6 = ldl_pivots(G6, [E, D, B, A, J])
    check("G4 (6): LDL on (E,D,B,A,J): pivots e, 1, r^2, 0 (zero residual row), j - 1/e - tau^2",
          st6 == 'ok' and sp.simplify(piv6[3]) == 0 and sp.simplify(piv6[4] - (j - 1/e - tau**2)) == 0, [sp.factor(v) for v in piv6])
    wv = r*s_
    check("G4 (6): tau^2 = b(1-w)^2/(b-w)^2 = a(1-w)^2/(a-w)^2, b = r^2, a = s^2, w = sqrt(ab)",
          sp.simplify(tau**2 - r**2*(1 - wv)**2/(r**2 - wv)**2) == 0 and sp.simplify(tau**2 - s_**2*(1 - wv)**2/(s_**2 - wv)**2) == 0)
    pr6 = dd**2/e + r**2*(dd - wv)**2/(r**2 - wv)**2
    check("G4 (6): printed homogeneous threshold = d * normalised threshold",
          sp.simplify(pr6 - dd*(1/(e/dd) + (tau**2).subs({r: r/sp.sqrt(dd), s_: s_/sp.sqrt(dd)}, simultaneous=True))) == 0)
    pivd, _ = ldl_pivots(Gds, [J, E, D, B, A])
    check("G4 double-SOC: pivots j, e - 1/j, 1, b, a - 1/b  (PSD iff ej >= 1 and ab >= 1)",
          all(sp.simplify(u1 - u2) == 0 for u1, u2 in zip(pivd, [j, e - 1/j, 1, b, a - 1/b])))
    # T25
    pb, _ = ldl_pivots(Gb12, [0, 1, 2, 3, 4])
    check("G4 T25: Gbase pivots a, hb, q - z^2/hb, c, 1",
          all(sp.simplify(u1 - u2) == 0 for u1, u2 in zip(pb, [a, hb, q - z**2/hb, c, 1])))
    e_t = t**2 + t**2*(a + b + 2*t)/(a*b - t**2)
    _, sc12 = out[('T25', '(12) general z')]
    _, sc9 = out[('T25', '(9) z = r, P = (r/c)C')]
    _, sc10 = out[('T25', '(10) z = w r, q = w^2 r^2/hb')]
    check("G4 T25: e-coefficient = e(t) of (8) in (9), (10), (12)", all(sp.simplify(s[1] - e_t) == 0 for s in (sc9, sc10, sc12)))
    check("G4 T25 (12): p-coefficient = r^2/c + (r-z)^2/(hb+q-2z)", sp.simplify(sc12[2] - (rr**2/c + (rr - z)**2/(hb + q - 2*z))) == 0)
    check("G4 T25 (9): p-coefficient = r^2/c", sp.simplify(sc9[2] - rr**2/c) == 0)
    check("G4 T25 (10): p-coefficient = r^2/c + hb r^2 (1-w)^2/(hb - w r)^2, q-coefficient = w^2 r^2/hb",
          sp.simplify(sc10[2] - (rr**2/c + hb*rr**2*(1 - w)**2/(hb - w*rr)**2)) == 0 and sp.simplify(sc10[6] - w**2*rr**2/hb) == 0)
    check("G4 T25 (10): hb + q - 2z = (hb - w r)^2/hb at z = wr, q = w^2r^2/hb", sp.simplify(hb + q10 - 2*w*rr - (hb - w*rr)**2/hb) == 0)
    check("G4 T25 (10): third Gbase pivot q - z^2/hb vanishes at z = wr, q = w^2r^2/hb", sp.simplify(q10 - (w*rr)**2/hb) == 0)
    check("G4 T25 usage, region II: q hb = bq(1 - t^2/(ab)) and 1 - t^2 - r^2 = 2t(1-t)",
          sp.simplify(q*hb - b*q*(1 - t**2/(a*b))) == 0 and sp.expand(1 - t**2 - rr**2 - 2*t*(1 - t)) == 0)
    X_, Z_, W_, B_ = sp.symbols('X Z W B', positive=True)
    bb = 1 + B_; tZ = Z_/(1 + Z_); aa = (1 + X_)/bb; Dd = (1 + X_)*(1 + Z_)**2 - Z_**2
    hbZ = bb - tZ**2/aa; rZ = 1/(1 + Z_); wW = W_/(1 + W_)
    check("G4 T25 usage, IIIb- (21): q = W^2(1+X)/((1+W)^2 b D) equals w^2 r^2/hb with w = W/(1+W)",
          sp.simplify(W_**2*(1 + X_)/((1 + W_)**2*bb*Dd) - wW**2*rZ**2/hbZ) == 0 and sp.simplify(hbZ - bb*Dd/((1 + X_)*(1 + Z_)**2)) == 0)
    # F3: T41 transports of (5), (6)
    A_ = a*b; R_ = (a + 1)*(b + 1)
    check("F3 T41 (3): 1/j <= p + (A-1)/R is T24 (5) with (j_T24, e_T24, b_T24, a_T24) = (p, j, a, b)",
          sp.simplify((1/j - 1/Lc) - (1/j - (A_ - 1)/R_)) == 0)
    zz = sp.Symbol('zz', positive=True)
    check("F3 T41 (4): (1-z)^2/(S-2z) = a(1-z)^2/(a-z)^2 when b = z^2/a",
          sp.simplify(((1 - zz)**2/(a + zz**2/a - 2*zz)) - a*(1 - zz)**2/(a - zz)**2) == 0)
    check("F3 T41 (4): a(1-z)^2/(a-z)^2 is T24 (6)'s tau^2 with b_T24 = a = r^2, a_T24 = b = s^2, z = rs",
          sp.simplify(tau**2 - (r**2*(1 - r*s_)**2/(r**2 - r*s_)**2)) == 0)
    check("F3 T41 (4) at ab = 0: tau^2 = 1/S (s = 0 member of (6))", sp.simplify(tau.subs(s_, 0)**2 - 1/r**2) == 0)
    pr41 = dd**2/j - dd*(a*b - dd**2)/((a + dd)*(b + dd))
    check("F3 T41 large_product (1): d^2/j <= p + d(ab-d^2)/((a+d)(b+d)) is the homogeneous T24 (5)", sp.simplify(pr41 - pr5.subs(e, j)) == 0)
    # F4: inventory56 (1) allocation of d between two chains
    qq, jj, ee, aa2 = sp.symbols('qq jj ee aa2', positive=True)
    al_ = qq*jj/(qq + jj); be_ = ee*aa2/(ee + aa2)
    d1 = dd*al_/(al_ + be_); d2 = dd*be_/(al_ + be_)
    Bch = dd**2*(qq + jj)*(ee + aa2)/(qq*jj*(ee + aa2) + ee*aa2*(qq + jj))
    check("F4 inv56 (1): d1 + d2 = d and d1^2(1/q+1/j) + d2^2(1/e+1/a) = Bch",
          sp.simplify(d1 + d2 - dd) == 0 and sp.simplify(d1**2*(1/qq + 1/jj) + d2**2*(1/ee + 1/aa2) - Bch) == 0)
    check("F4 inv56 large_product (1): b(qj/(q+j) + ea/(e+a)) >= d^2  <=>  b >= Bch", sp.simplify(dd**2/(al_ + be_) - Bch) == 0)
    return out

# ================================================================ section T: orbits, transports, conic combinations
FMAP = {tuple(v): j for j, v in enumerate(L.F_VALUES)}
def gperm(g):
    """facet permutation j -> index of g.f_j."""
    return [FMAP[tuple(L.g_act_values(g, list(L.F_VALUES[jj])))] for jj in range(82)]
GPERMS = None
def gperms():
    global GPERMS
    if GPERMS is None:
        GPERMS = [gperm(g) for g in L.G_ELEMENTS]
    return GPERMS

def in_orbit(rep, target, unordered=()):
    """rep, target: tuples of facet indices with roles; unordered: list of index groups whose order is free."""
    for gp in gperms():
        img = [gp[k] for k in rep]
        ok = True
        used = set()
        for grp in unordered:
            if sorted(img[i] for i in grp) != sorted(target[i] for i in grp): ok = False
            used |= set(grp)
        if ok and all(img[i] == target[i] for i in range(len(rep)) if i not in used):
            return True
    return False

def section_T():
    print("\n=== T: orbits of the representatives, transports, conic combinations ===")
    gp = gperms()
    check("T0 G acts on the 82 facets by 1152 permutations, preserving minor/non-minor",
          all(sorted(p) == list(range(82)) for p in gp) and all(L.F_IS_MINOR[p[k]] == L.F_IS_MINOR[k] for p in gp for k in range(82)))
    # rank-two representative (leaves 2, 29; centre 12; k 35): orbit size of the ordered tuple with {2,29} unordered
    orb = set()
    for p in gp:
        orb.add((frozenset((p[2], p[29])), p[12], p[35]))
    print("   rank-two orbit size (leaf pair unordered):", len(orb))
    check("T1 rank-two orbit has 144 members", len(orb) == 144)
    # T24 core transports printed in the documents: (59,53,55,67,54) with {55,67} unordered
    core = (59, 53, 55, 67, 54)
    for doc, tgt in [("T41_bad17 p.2", (1, 53, 4, 16, 17)), ("T27 p.2", (5, 2, 0, 67, 40)), ("T38 p.1", (5, 2, 67, 0, 40))]:
        check("T2 T24 core -> %s %s (b,a pair unordered)" % (doc, tgt), in_orbit(core, tgt, unordered=[(2, 3)]))
        check("T2 T24 core -> %s %s with the printed order of (b, a)" % (doc, tgt), in_orbit(core, tgt))
    # T25 (6), (7): j -> 2, e -> 39, d -> 70, {b,a} -> {6,67}
    check("T2 T24 core (59,53,{55,67},54) -> T25 (2,39,{6,67},70)", in_orbit(core, (2, 39, 6, 67, 70), unordered=[(2, 3)]))
    for o in [(2, 39, 6, 67, 70), (2, 39, 67, 6, 70)]:
        print("   ordered T24 core -> %s: %s" % (o, in_orbit(core, o)))
    # first-pass atoms: representatives
    SOC = (1, 9, 22); OLDH = (3, 15, 26); CHAIN = (2, 0, 12, 70); FOUR = (0, 6, 10, 54)
    tests = [("SOC (2,59;54) T24 (4)", SOC, (2, 59, 54), [(0, 1)]),
             ("old-H (53,59;54) aux,main T24 (4)", OLDH, (53, 59, 54), []),
             ("four (67,55;53;54) T24 (4)", FOUR, (67, 55, 53, 54), [(0, 1)]),
             ("chain centre 12 leaves (0,2), k 35  inv56 (1)", CHAIN, (12, 0, 2, 35), [(1, 2)]),
             ("chain centre 12 leaves (7,29), k 35  inv56 (1)", CHAIN, (12, 7, 29, 35), [(1, 2)]),
             ("four (7,0;5;35)  inv56 (3)", FOUR, (7, 0, 5, 35), [(0, 1)]),
             ("old-H (29,12;35) rank-two p.3", OLDH, (29, 12, 35), []),
             ("old-H (2,12;35) rank-two p.3", OLDH, (2, 12, 35), []),
             ("SOC (2,5;70) T25 (4)", SOC, (2, 5, 70), [(0, 1)]),
             ("old-H (39,2;70) T25 (4)", OLDH, (39, 2, 70), []),
             ("four (0,6;1;70) T25 (5)", FOUR, (0, 6, 1, 70), [(0, 1)]),
             ("four (6,67;39;70) T25 (5)", FOUR, (6, 67, 39, 70), [(0, 1)])]
    for nm, rep, tgt, un in tests:
        check("T3 " + nm, in_orbit(rep, tgt, unordered=un))
    # does the T25 partial-four Gram sit in the same face as T24 under G?  (informational)
    return

# ================================================================ section N: numeric attack (NUMERIC only)
FV = np.array(L.F_VALUES, dtype=float)          # 82 x 24
SIG = np.array(L.PERMS); SGNA = np.array(L.SGN, dtype=float)

def omat_num(y, lib=np):
    """O for x -> p x conj(q), (p, q) = y[:4], y[4:] normalised."""
    if lib is np:
        p = y[:4] / np.linalg.norm(y[:4]); q = y[4:] / np.linalg.norm(y[4:])
    else:
        npn = mp.sqrt(sum(t*t for t in y[:4])); nqn = mp.sqrt(sum(t*t for t in y[4:]))
        p = [t / npn for t in y[:4]]; q = [t / nqn for t in y[4:]]
    qc = [q[0], -q[1], -q[2], -q[3]]
    cols = []
    for j in range(4):
        e = [0, 0, 0, 0]; e[j] = 1
        cols.append(qmul(qmul(p, e), qc))
    if lib is np:
        return np.array(cols).T
    return [[cols[j][i] for j in range(4)] for i in range(4)]

def weights(O):
    return SGNA * O[0, SIG[:, 0]] * O[1, SIG[:, 1]] * O[2, SIG[:, 2]] * O[3, SIG[:, 3]]

def weights_mp(O):
    out = []
    for s, sg in zip(L.PERMS, L.SGN):
        t = O[0][s[0]] * O[1][s[1]] * O[2][s[2]] * O[3][s[3]]
        out.append(t if sg > 0 else -t)
    return out

def cost_fn(coef):
    """coef: dict facet -> coefficient; returns f(y) = Phi_cost(O(y)) and an mpmath evaluator."""
    vals = sum(c * FV[j] for j, c in coef.items())
    f = lambda y: float(vals @ weights(omat_num(y)))
    def fmp(y):
        O = omat_num([mp.mpf(t) for t in y], lib=mp)
        w = weights_mp(O)
        return sum(mp.mpf(sum(c * L.F_VALUES[j][k] for j, c in coef.items())) * w[k] for k in range(24))
    return f, fmp

def minor(O, I, J):
    return O[I[0]-1][J[0]-1] * O[I[1]-1][J[1]-1] - O[I[0]-1][J[1]-1] * O[I[1]-1][J[0]-1]

def attack(f, nstart, rng):
    best = (np.inf, None)
    for _ in range(nstart):
        y0 = rng.normal(size=8)
        r = minimize(f, y0, method='BFGS', options={'gtol': 1e-10, 'maxiter': 400})
        if r.fun < best[0]: best = (float(r.fun), r.x)
    return best

def run_family(name, samples, nper, seed, weaken=None):
    """samples: list of coefficient dicts on the threshold; weaken: function coef -> weakened coef."""
    rng = np.random.default_rng(seed)
    t0 = time.time(); worst = (np.inf, None, None)
    for coef in samples:
        f, fmp = cost_fn(coef)
        v, y = attack(f, nper, rng)
        if v < worst[0]: worst = (v, y, coef)
    v, y, coef = worst
    _, fmp = cost_fn(coef)
    mp.mp.dps = 50
    vm = fmp(list(y))
    line = "%-34s %3d coeff points x %2d starts = %4d starts: min %.3e (mpmath 50 digits %.3e) [%.0fs]" % (
        name, len(samples), nper, len(samples) * nper, v, float(vm), time.time() - t0)
    if weaken is not None:
        wv = np.inf
        rng2 = np.random.default_rng(seed + 100)
        for coef in samples[:6]:
            f2, _ = cost_fn(weaken(coef))
            wv = min(wv, attack(f2, 8, rng2)[0])
        line += "; weakened variant min %.3e" % wv
    print(line, flush=True)
    return v, vm

def section_N():
    print("\n=== N: numeric attack (multistart BFGS over (p,q) in R^8, normalised; NUMERIC only) ===")
    rng = np.random.default_rng(20261125)
    res = {}
    def lu(lo, hi): return float(np.exp(rng.uniform(np.log(lo), np.log(hi))))
    # (R) coefficient form on the threshold b = d/4 + d^2/(4 min(a,j)), d = 1
    S = []
    for i in range(25):
        a, j = lu(0.02, 50), lu(0.02, 50)
        if i % 4 == 0: j = a
        x = min(a, j)
        S.append({12: 0.25 + 1/(4*x), 2: j, 29: a, 35: 1.0})
    res['R coeff'] = run_family("(R) b = d/4 + d^2/(4 min(a,j))", S, 20, 11,
                                weaken=lambda c: {**c, 12: 0.0})
    # (R) pointwise form
    def fR(y):
        O = omat_num(y); w = weights(O)
        return float(FV[35] @ w + (FV[12] @ w)/4 + abs(minor(O, (1, 4), (2, 4))) * math.sqrt(max(FV[2] @ w + FV[29] @ w, 0)))
    rngR = np.random.default_rng(12)
    v, y = attack(fR, 500, rngR)
    mp.mp.dps = 50
    O = omat_num([mp.mpf(t) for t in y], lib=mp); w = weights_mp(O)
    ph = lambda j: sum(L.F_VALUES[j][k] * w[k] for k in range(24))
    vm = ph(35) + ph(12)/4 + abs(minor(O, (1, 4), (2, 4))) * mp.sqrt(ph(2) + ph(29))
    def fRw(y):
        O = omat_num(y); w = weights(O)
        return float(FV[35] @ w + (FV[12] @ w)/4 + 0.5 * abs(minor(O, (1, 4), (2, 4))) * math.sqrt(max(FV[2] @ w + FV[29] @ w, 0)))
    vw = attack(fRw, 40, np.random.default_rng(112))[0]
    def fRw2(y):
        O = omat_num(y); w = weights(O)
        return float(FV[35] @ w + (FV[12] @ w)/4 + abs(minor(O, (1, 4), (2, 4))) * abs(minor(O, (1, 3), (1, 2))))
    vw2 = attack(fRw2, 60, np.random.default_rng(113))[0]
    print("%-34s 500 starts: min %.3e (mpmath %.3e); weakened (sqrt term x 1/2) min %.3e; "
          "variant with |x29| in place of sqrt(Phi2+Phi29) min %.3e" % ("(R) pointwise (2)", v, float(vm), vw, vw2), flush=True)
    # (F1) T24 face, d = 1, q = p = c = 0; facets j:59 e:53 d:54 b:55 a:67
    def T24(j, e, b, a): return {59: j, 53: e, 54: 1.0, 55: b, 67: a}
    S5 = []
    for i in range(25):
        a = lu(0.05, 20); b = lu(1.0/a * 1.0001, 1.0/a * 50) if i % 5 else 1.0/a * 1.02
        Lc = (a + 1)*(b + 1)/(a*b - 1)
        e = lu(0.05, 3*Lc) if i % 6 else Lc
        S5.append(T24(max(0.0, 1/e - 1/Lc), e, b, a))
    res['F1 (5)'] = run_family("(F1) (5) coherent boundary", S5, 20, 21, weaken=lambda c: {**c, 59: 0.5*c[59], 53: 0.5*c[53]})
    S6 = []
    for i in range(25):
        r = lu(0.05, 10); s = lu(0.05, 10)
        if i % 5 == 0: s = 0.0
        if i % 7 == 1: s = 1/r * 0.7        # ab < 1 chart of the PDF
        if abs(r - s) < 1e-3: s = r * 1.5
        e = lu(0.05, 20); tau = (1 - r*s)/(r - s)
        S6.append(T24(1/e + tau*tau, e, r*r, s*s))
    res['F1 (6)'] = run_family("(F1) (6) parallel, all a != b", S6, 20, 22, weaken=lambda c: {**c, 59: 0.5*c[59]})
    Sd = []
    for i in range(25):
        a = lu(0.05, 20); b = (1.0/a) * (1 if i % 3 == 0 else lu(1, 30)); e = lu(0.05, 20)
        Sd.append(T24(1/e, e, b, a))
    res['F1 dSOC'] = run_family("(F1) double-SOC, ab >= 1, ej = 1", Sd, 20, 23, weaken=lambda c: {**c, 59: 0.5*c[59]})
    # (F2) T25 face: a f0 + e f1 + p f2 + j f5 + b f6 + c f39 + q f67 + f70
    def T25(a, e, p, b, c, q): return {0: a, 1: e, 2: p, 6: b, 39: c, 67: q, 70: 1.0}
    S9 = []; S10 = []
    for i in range(25):
        a = lu(0.05, 20); t = 0.0 if i % 6 == 0 else float(rng.uniform(0, 0.99))
        b = (t*t/a) * (1.0001 if i % 5 == 1 else 1) + lu(0.01, 20)
        c = lu(0.05, 20); r = 1 - t; hb = b - t*t/a
        et = t*t + t*t*(a + b + 2*t)/(a*b - t*t)
        S9.append(T25(a, et, r*r/c, b, c, r*r/hb))
        w = [0.0, 1.0][i % 2] if i % 5 == 0 else float(rng.uniform(0, 1))
        if hb <= w*r: w = 0.5 * hb / r
        S10.append(T25(a, et, r*r/c + hb*r*r*(1 - w)**2/(hb - w*r)**2, b, c, w*w*r*r/hb))
    res['F2 (9)'] = run_family("(F2) (9) boundary q hb = r^2", S9, 20, 31, weaken=lambda c: {**c, 2: 0.0})
    res['F2 (10)'] = run_family("(F2) (10) boundary, w in [0,1]", S10, 20, 32, weaken=lambda c: {**c, 2: 0.0})
    # (F3) T41 bad17: q f3 + p f1 + j f53 + c f6 + e f15 + f17 + b f16 + a f4
    def T41(p, j, b, a): return {1: p, 53: j, 17: 1.0, 16: b, 4: a}
    S3 = []; S4 = []
    for i in range(25):
        a = lu(0.05, 20); b = (1/a) * (1 if i % 4 == 0 else lu(1, 30)); j = lu(0.05, 20)
        R = (a + 1)*(b + 1)
        S3.append(T41(max(0.0, 1/j - (a*b - 1)/R), j, b, a))
        z = float(rng.uniform(0.02, 0.98)); a2 = z * lu(1.0001, 30); b2 = z*z/a2
        if i % 6 == 0: b2 = 0.0; a2 = lu(1, 20)
        k0 = a2*(1 - z)**2/(a2 - z)**2 if b2 > 0 else 1/a2
        S4.append(T41(1/j + k0, j, b2, a2))
    res['F3 (3)'] = run_family("(F3) T41 (3) boundary", S3, 20, 41, weaken=lambda c: {**c, 1: 0.5*c[1], 53: 0.5*c[53]})
    res['F3 (4)'] = run_family("(F3) T41 (4) boundary, incl. ab = 0", S4, 20, 42, weaken=lambda c: {**c, 1: 0.0})
    # (F4) inventory56: q f0 + j f2 + c f5 + e f7 + b f12 + a f29 + f35
    def I56(q, j, c, e, b, a): return {0: q, 2: j, 5: c, 7: e, 12: b, 29: a, 35: 1.0}
    S1 = []; S3b = []
    for i in range(25):
        q, j, e, a = lu(0.05, 20), lu(0.05, 20), lu(0.05, 20), lu(0.05, 20)
        S1.append(I56(q, j, 0.0, e, (q + j)*(e + a)/(q*j*(e + a) + e*a*(q + j)), a))
        q2 = lu(0.05, 20); e2 = (1/q2) * lu(1.0001, 40)
        S3b.append(I56(q2, 0.0, (q2*e2 + e2 + q2 + 1)/(q2*e2 - 1), e2, 0.0, 0.0))
    res['F4 (1)'] = run_family("(F4) inv56 (1) b = Bch", S1, 20, 51, weaken=lambda c: {**c, 12: 0.0})
    res['F4 (3)'] = run_family("(F4) inv56 (3) four boundary", S3b, 20, 52, weaken=lambda c: {**c, 5: 0.0})
    worst = min(min(v[0], float(v[1])) for v in res.values())
    print("NUMERIC: global minimum over all in-region families %.3e (pointwise (R) %.3e)" % (worst, v))
    return res

if __name__ == "__main__":
    t0 = time.time()
    section_R()
    if BASE: section_G()
    else: print("MOC4_SUPP not set: section G skipped")
    section_T()
    print("\nEXACT PART: %d passed, %d failed (%.0f s)" % (NPASS, NFAIL, time.time() - t0))
    section_N()
    print("total %.0f s" % (time.time() - t0))
