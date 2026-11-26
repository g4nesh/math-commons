"""Independent verification of the four SO(4) atom families (SOC, old-H, chain, four).

Sections (run one at a time, deterministic seeds):
  exact     polynomial identities behind the hand proofs (bihomogenised, exact over Q),
            chain Gram identity re-checked from certificate JSON data (read as data only)
  orbits    covariance under G (1152 elements), det -1 handling, orbit sizes, role symmetries,
            membership of the cited tuples
  attack F  numeric attack on family F in {SOC, oldH, chain, four}: multistarts on the
            representative and on 30 further orbit members, mpmath re-evaluation
  variants  wrong variants (roles swapped, factor 1.8, one absolute value dropped), the
            half-compensation chain, and an extra multiplier certificate for chain; negative
            values are confirmed EXACTLY at rational points of SO(4)

Usage:  MOC4_SUPP=<supplement root> python -I -B so4_atoms_verify.py SECTION [FAMILY]
Only the chain check in `exact` reads the supplement, with json.load; nothing there is executed.
"""
import os, sys, json, itertools, time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lib'))
import moc4lib as L
import flint
import sympy as sp
import numpy as np
from scipy.optimize import minimize
import mpmath

SUPP = os.environ.get('MOC4_SUPP')
ARGS = sys.argv[1:]
NFAIL = 0


def check(name, ok, info=''):
    global NFAIL
    if not ok:
        NFAIL += 1
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({info})" if info else ''), flush=True)


REP = {'SOC': (1, 9, 22), 'oldH': (12, 2, 70), 'chain': (2, 0, 12, 70), 'four': (0, 6, 10, 54)}
MINOR_IJ = {j: ([int(c) - 1 for c in L.F_TABLE[j][1]], [int(c) - 1 for c in L.F_TABLE[j][2]])
            for j in range(82) if L.F_IS_MINOR[j]}


def minor(O, j):
    I, J = MINOR_IJ[j]
    if len(I) == 1:
        return O[I[0]][J[0]]
    return O[I[0]][J[0]] * O[I[1]][J[1]] - O[I[0]][J[1]] * O[I[1]][J[0]]


def _comp(j):
    I, J = MINOR_IJ[j]
    if len(I) == 1:
        return lambda O: minor(O, j)
    Ic = [i for i in range(4) if i not in I]
    Jc = [i for i in range(4) if i not in J]
    return lambda O: O[Ic[0]][Jc[0]] * O[Ic[1]][Jc[1]] - O[Ic[0]][Jc[1]] * O[Ic[1]][Jc[0]]


COMP = {j: _comp(j) for j in MINOR_IJ}


def weights(O):
    out = []
    for s, sg in zip(L.PERMS, L.SGN):
        t = O[0][s[0]] * O[1][s[1]] * O[2][s[2]] * O[3][s[3]]
        out.append(t if sg > 0 else -t)
    return out


def Phi(j, W, zero):
    return sum((v * w for v, w in zip(L.F_VALUES[j], W) if v), zero)


def Lmat(p):
    p0, p1, p2, p3 = p
    return [[p0, -p1, -p2, -p3], [p1, p0, -p3, p2], [p2, p3, p0, -p1], [p3, -p2, p1, p0]]


def Rmat(q):
    q0, q1, q2, q3 = q
    return [[q0, -q1, -q2, -q3], [q1, q0, q3, -q2], [q2, -q3, q0, q1], [q3, q2, -q1, q0]]


def mm(A, B, zero):
    return [[sum((A[i][k] * B[k][j] for k in range(4)), zero) for j in range(4)] for i in range(4)]


PI_OLDH = [[1, 0, 0, 0], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0]]
PI_ID = [[1 if i == j else 0 for j in range(4)] for i in range(4)]


def Omodel(p, q, Pi, zero):
    qb = (q[0], -q[1], -q[2], -q[3])
    return mm(mm(Lmat(p), Rmat(qb), zero), [[zero + c for c in r] for r in Pi], zero)


# =============================================================================================
def section_exact():
    names = ('p0', 'p1', 'p2', 'p3', 'q0', 'q1', 'q2', 'q3')
    C = flint.fmpq_mpoly_ctx.get(names, 'lex')
    Y = C.gens()
    Z = C.from_dict({})
    p, q = Y[:4], Y[4:]
    NP = sum((x * x for x in p), Z)
    NQ = sum((x * x for x in q), Z)
    P = [x * x for x in p]
    Q = [x * x for x in q]
    Dh = [P[i] * NQ - Q[i] * NP for i in range(4)]          # homogenised D_i = P_i - Q_i, (2,2)

    def signed_eq(a, b):
        return 1 if a == b else (-1 if a == -b else 0)

    print('Reduction method: every identity is bihomogenised with |p|^2 and |q|^2 to bidegree (4,4) '
          '(or (2,2) for 2x2 minors) and compared as polynomials over Q (python-flint). '
          'On S^3 x S^3 the homogenised and original forms agree.')
    for nm, Pi in (('Pi=I', PI_ID), ('Pi=oldH', PI_OLDH)):
        O = Omodel(p, q, Pi, Z)
        ok = all(sum((O[k][i] * O[k][j] for k in range(4)), Z) == (NP * NQ if i == j else Z)
                 for i in range(4) for j in range(4))
        check(f'E0 [{nm}] O^T O = |p|^2|q|^2 I', ok)
        check(f'E0 [{nm}] det O = +(|p|^2|q|^2)^2', sum(weights(O), Z) == (NP * NQ) ** 2)

    # ---------------- SOC, representative (f1, f9; f22), Pi = I
    O = Omodel(p, q, PI_ID, Z)
    W = weights(O)
    x1, x9 = minor(O, 1), minor(O, 9)
    s1 = signed_eq(x1, Dh[0] + Dh[1])
    s9 = signed_eq(x9, Dh[0] + Dh[3])
    check('E1 SOC x1 = +-(D0+D1)', s1 != 0, f'sign {s1}')
    check('E1 SOC x9 = +-(D0+D3)', s9 != 0, f'sign {s9}')
    Kh = (p[0] * p[2] - p[1] * p[3]) ** 2 * NQ * NQ + (q[0] * q[2] - q[1] * q[3]) ** 2 * NP * NP
    check('E1 SOC Phi22 = 2K - 2 D0 D2 - 2 D1 D3', Phi(22, W, Z) == 2 * Kh - 2 * Dh[0] * Dh[2] - 2 * Dh[1] * Dh[3])
    check('E1 SOC sum D = 0 (homogenised, exact)', sum(Dh, Z) == Z)
    check('E1 SOC x1 x9 = +-(D1 D3 - D0 D2)', signed_eq(x1 * x9, Dh[1] * Dh[3] - Dh[0] * Dh[2]) != 0)
    check('E1 SOC Phi22 + 2 sgn x1 x9 = 2K - 4 D0 D2', Phi(22, W, Z) + 2 * s1 * s9 * x1 * x9 == 2 * Kh - 4 * Dh[0] * Dh[2])
    check('E1 SOC Phi22 - 2 sgn x1 x9 = 2K - 4 D1 D3', Phi(22, W, Z) - 2 * s1 * s9 * x1 * x9 == 2 * Kh - 4 * Dh[1] * Dh[3])
    check('E1 SOC Phi1 = x1^2, Phi9 = x9^2', Phi(1, W, Z) == x1 * x1 and Phi(9, W, Z) == x9 * x9)

    # ---------------- old-H, representative (aux f12, main f2; f70), Pi = PI_OLDH
    O = Omodel(p, q, PI_OLDH, Z)
    W = weights(O)
    x2, x12 = minor(O, 2), minor(O, 12)
    s2 = signed_eq(x2, Dh[0] + Dh[2])
    s12 = signed_eq(x12, Dh[0] + Dh[3])
    check('E2 oldH x2 = +-(D0+D2)', s2 != 0, f'sign {s2}')
    check('E2 oldH x12 = +-(D0+D3)', s12 != 0, f'sign {s12}')
    Kp = ((p[0] * q[3] - p[2] * q[1]) ** 2 + (p[1] * q[2] - p[3] * q[0]) ** 2) * NP * NQ
    Phis = Phi(2, W, Z) + Phi(70, W, Z)
    e = s2 * s12
    check('E2 oldH Phi2 + Phi70 + 2 sgn x2 x12 = 2K\' - 4 D0 D1', Phis + 2 * e * x2 * x12 == 2 * Kp - 4 * Dh[0] * Dh[1])
    check('E2 oldH Phi2 + Phi70 - 2 sgn x2 x12 = 2K\' - 4 D2 D3', Phis - 2 * e * x2 * x12 == 2 * Kp - 4 * Dh[2] * Dh[3])
    check('E2 oldH Phi2 = x2^2', Phi(2, W, Z) == x2 * x2)
    r0, r3, t0, t3 = sp.symbols('r0 r3 t0 t3')
    check('E2 Aczel identity (r0 t3 - t0 r3)^2 - (r0^2 - t0^2)(t3^2 - r3^2) = (r0 r3 - t0 t3)^2',
          sp.expand((r0 * t3 - t0 * r3) ** 2 - (r0 ** 2 - t0 ** 2) * (t3 ** 2 - r3 ** 2) - (r0 * r3 - t0 * t3) ** 2) == 0)
    # K' symmetric under p <-> q
    sw = [Y[4], Y[5], Y[6], Y[7], Y[0], Y[1], Y[2], Y[3]]
    check('E2 K\' invariant under p <-> q', Kp.compose(*sw) == Kp)
    # chain reformulation used in the report: Phi70 = 2K' - (D1 - D3)(D0 - D2) on the spheres
    check('E2 Phi70 = 2K\' - (D1-D3)(D0-D2) (homogenised)', Phi(70, W, Z) == 2 * Kp - (Dh[1] - Dh[3]) * (Dh[0] - Dh[2]))

    # ---------------- four, representative (f0, f6; f10; f54), Pi = I
    O = Omodel(p, q, PI_ID, Z)
    W = weights(O)
    x0, x6, x10 = minor(O, 0), minor(O, 6), minor(O, 10)
    m03 = p[0] * q[3] - p[3] * q[0]
    m12 = p[1] * q[2] - p[2] * q[1]
    a0 = signed_eq(x0, m03 - m12)
    a6 = signed_eq(x6, m03 + m12)
    a10 = signed_eq(x10, Dh[0] + Dh[2])
    check('E3 four x0 = +-(m03 - m12)', a0 != 0, f'sign {a0}')
    check('E3 four x6 = +-(m03 + m12)', a6 != 0, f'sign {a6}')
    check('E3 four x10 = +-(D0 + D2)', a10 != 0, f'sign {a10}')
    S2 = NP * NQ
    Xh = (P[0] - P[3]) * NQ - (Q[1] - Q[2]) * NP
    Yh = (P[1] - P[2]) * NQ - (Q[0] - Q[3]) * NP
    check('E3 four Phi54 = |p|^2|q|^2 (x0^2 + x6^2) - X Y', Phi(54, W, Z) == S2 * (x0 * x0 + x6 * x6) - Xh * Yh)
    check('E3 four Phi10 = x10^2', Phi(10, W, Z) == x10 * x10)
    check('E3 four Phi0 = |p|^2|q|^2 x0^2, Phi6 = |p|^2|q|^2 x6^2',
          Phi(0, W, Z) == S2 * x0 * x0 and Phi(6, W, Z) == S2 * x6 * x6)
    Ps = sp.symbols('P0:4'); Qs = sp.symbols('Q0:4')
    lam0 = (Ps[0] - Qs[0]) - (Ps[3] - Qs[3]); lam1 = (Ps[1] - Qs[1]) - (Ps[2] - Qs[2])
    eps0 = Ps[0] + Qs[0] - Ps[3] - Qs[3]; eps1 = Ps[1] + Qs[1] - Ps[2] - Qs[2]
    X = (Ps[0] - Ps[3]) - (Qs[1] - Qs[2]); Yv = (Ps[1] - Ps[2]) - (Qs[0] - Qs[3])
    check('E3 four X Y = (lam0+lam1)^2/4 - sigma^2, sigma = (eps0-eps1)/2',
          sp.expand(X * Yv - ((lam0 + lam1) ** 2 / 4 - (eps0 - eps1) ** 2 / 4)) == 0)
    Dl = [Ps[i] - Qs[i] for i in range(4)]
    check('E3 four D0 + D2 - (lam0 - lam1)/2 = (sum D)/2', sp.expand(Dl[0] + Dl[2] - (lam0 - lam1) / 2 - sum(Dl) / 2) == 0)
    a, b, c, d = sp.symbols('a b c d', real=True)
    z, w = a + sp.I * b, c + sp.I * d
    lhs = sp.expand((z ** 2 - w ** 2) * sp.conjugate(z ** 2 - w ** 2))
    rhs = sp.expand((a * a + b * b - c * c - d * d) ** 2 + 4 * (a * d - b * c) ** 2)
    check('E3 four pair fact |z^2 - w^2|^2 = (|z|^2 - |w|^2)^2 + 4 Im(conj(z) w)^2', sp.simplify(lhs - rhs) == 0)
    check('E3 four Re(z^2 - w^2) = lam-type, Im(conj z w) = a d - b c',
          sp.expand(sp.re(z ** 2 - w ** 2) - ((a * a - b * b) - (c * c - d * d))) == 0
          and sp.expand(sp.im(sp.conjugate(z) * w) - (a * d - b * c)) == 0)

    # ---------------- chain: re-check of the supplement (PB) Gram identity from JSON DATA
    if not SUPP:
        check('E4 chain Gram identity (needs MOC4_SUPP)', False, 'MOC4_SUPP not set')
        return
    base = os.path.join(SUPP, 'proofs', 'manin_transfer_20261021')
    with open(os.path.join(base, 'blind_moc4_next_singlebad_cubic_reduction_20261022.json')) as fh:
        CUB = json.load(fh)
    with open(os.path.join(base, 'blind_moc4_next_singlebad_budget_PSD8_20261022.json')) as fh:
        GD = json.load(fh)
    SY = sp.symbols(names)

    def parse(s):
        e = sp.Poly(sp.sympify(s, locals=dict(zip(names, SY))), *SY)
        return C.from_dict({m: flint.fmpq(int(cf.p), int(cf.q)) for m, cf in zip(e.monoms(), e.coeffs())})

    Fe = [[parse(s) for s in par] for par in CUB['features']]
    T = sp.Matrix([[sp.Rational(v) for v in row] for row in GD['T']])
    Ti = T.inv()
    src = GD['source_order']
    check('E4 data: source order lists f0, f2, f12, f70', all(f in src for f in (0, 2, 12, 70)), str(src))
    g = [[sum((flint.fmpq(int(Ti[i, k].p), int(Ti[i, k].q)) * Fe[par][k] for k in range(8)), Z)
          for i in range(8)] for par in range(2)]
    gs = [[[gi.compose(*[sgn * Y[ix] for ix, sgn in sm]) for gi in g[par]] for sm in CUB['left_four']]
          for par in range(2)]

    def ip(i, k):
        return sum((gs[par][j][i] * gs[par][j][k] for par in range(2) for j in range(4)), Z)

    def qmul(x, y):
        a0, a1, a2, a3 = x
        b0, b1, b2, b3 = y
        return (a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3, a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
                a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1, a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0)
    one = Z + 1
    cols = [qmul(qmul(p, tuple(one if r == cc else Z for r in range(4))), q) for cc in range(4)]
    Op = [[cols[cc][r] for cc in range(4)] for r in range(4)]
    check('E4 model x -> p x q: det O = +(|p|^2|q|^2)^2', sum(weights(Op), Z) == (NP * NQ) ** 2)
    Wp = weights(Op)
    SR = NQ - p[1] ** 2 - p[2] ** 2 - p[3] ** 2

    def redN(f):   # normal form modulo N = |p|^2 - |q|^2: p0^2 -> |q|^2 - p1^2 - p2^2 - p3^2
        out = Z
        for mon, cf in f.to_dict().items():
            k = mon[0]
            out += C.from_dict({(0,) + tuple(mon[1:]): cf}) * (p[0] ** (k % 2)) * (SR ** (k // 2))
        return out
    ix = {f: i for i, f in enumerate(src)}
    for f in (0, 2, 12):
        xf = minor(Op, f)
        Pf = Phi(f, Wp, Z)
        ent = len(MINOR_IJ[f][0]) == 1
        check(f'E4 chain Phi{f} = {"|p|^2|q|^2 " if ent else ""}x{f}^2', Pf == ((NP * NQ) if ent else one) * xf * xf)
        check(f'E4 chain 2 Phi{f} = |p|^2 |u{f}|^2 mod N', redN(2 * Pf - NP * ip(ix[f], ix[f])) == Z)
    pb = ip(ix[70], ix[70]) - 2 * ip(ix[0], ix[2]) - 2 * ip(ix[2], ix[12])
    check('E4 chain (PB) 2 Phi70 = |p|^2 (|u70|^2 - 2<u0,u2> - 2<u2,u12>) mod N', redN(2 * Phi(70, Wp, Z) - NP * pb) == Z)
    check('E4 chain |u70|^2 is a sum of squares of 8 real cubics (by construction)', True,
          'u70 = (g_70(S_j y))_{parity,j}, 8 coordinates')


# =============================================================================================
VAL2IDX = {L.F_VALUES[j]: j for j in range(82)}


def gmap(g):
    return [VAL2IDX[tuple(L.g_act_values(g, list(L.F_VALUES[j])))] for j in range(82)]


def orbit(rep):
    return sorted(set(tuple(m[j] for j in rep) for m in GM))


def section_orbits():
    global GM
    GM = [gmap(g) for g in L.G_ELEMENTS]
    check('O0 |G| = 1152 and every g permutes the 82 facets', len(GM) == 1152 and all(sorted(m) == list(range(82)) for m in GM))
    check('O0 G maps minors to minors', all(L.F_IS_MINOR[m[j]] == L.F_IS_MINOR[j] for m in GM for j in range(82)))
    C16 = flint.fmpz_mpoly_ctx.get(tuple(f'o{i}{j}' for i in range(4) for j in range(4)), 'lex')
    G16 = C16.gens()
    Z = C16.from_dict({})
    O = [[G16[4 * i + j] for j in range(4)] for i in range(4)]
    WO = weights(O)
    detO = sum(WO, Z)
    okw = okm = True
    neps = {1: 0, -1: 0}
    for g, m in zip(L.G_ELEMENTS, GM):
        gO = L.g_act_O(g, O)
        Wg = weights(gO)
        t = L.g_index_map(g)
        eps = 1 if sum(Wg, Z) == detO else (-1 if sum(Wg, Z) == -detO else 0)
        neps[eps] = neps.get(eps, 0) + 1
        if eps == 0 or any(Wg[t[s]] != eps * WO[s] for s in range(24)):
            okw = False
        for j in MINOR_IJ:
            a = minor(gO, m[j])
            if a not in (minor(O, j), -minor(O, j), COMP[j](O), -COMP[j](O)):
                okm = False
    check('O1 w_{g.sigma}(g.O) = eps(g) w_sigma(O), eps = det(g.O)/det(O) = +-1 (generic O, all g)', okw, f'eps counts {neps}')
    check('O1 hence Phi_{g.j}(g.O) = eps(g) Phi_j(O) and Phi^O4 = Phi/det is G-invariant', okw)
    check('O1 minor_{g.j}(g.O) = +- minor_j(O) or +- its complementary minor, all 34 minors and all g', okm)
    # Jacobi on SO(4): the complementary 2x2 minor equals +- the minor (exact in the quaternion
    # model, which covers SO(4)); for det -1 apply it to O J. Hence |x_{g.j}(g.O)| = |x_j(O)| on O(4).
    Cq = flint.fmpq_mpoly_ctx.get(('p0', 'p1', 'p2', 'p3', 'q0', 'q1', 'q2', 'q3'), 'lex')
    Yq = Cq.gens()
    Zq = Cq.from_dict({})
    Om = Omodel(Yq[:4], Yq[4:], PI_ID, Zq)
    okj = all(COMP[j](Om) in (minor(Om, j), -minor(Om, j)) for j in MINOR_IJ if len(MINOR_IJ[j][0]) == 2)
    check('O1 Jacobi: complementary 2x2 minor = +- minor on SO(4) (quaternion model, exact)', okj)
    J = [[Z + (1 if (i == j and i < 3) else (-1 if i == j else 0)) for j in range(4)] for i in range(4)]
    OJ = [[sum((O[i][k] * J[k][j] for k in range(4)), Z) for j in range(4)] for i in range(4)]
    WJ = weights(OJ)
    check('O2 det -1: Phi_c(O J) = -Phi_c(O) for J = diag(1,1,1,-1), all 82 c',
          all(Phi(c, WJ, Z) == -Phi(c, WO, Z) for c in range(82)))
    check('O2 det -1: minor_j(O J) = +- minor_j(O)', all(minor(OJ, j) in (minor(O, j), -minor(O, j)) for j in MINOR_IJ))

    orbs = {k: orbit(v) for k, v in REP.items()}
    soc_u = set((frozenset(t[:2]), t[2]) for t in orbs['SOC'])
    check('O3 SOC ordered orbit 288, unordered pairs 144; pair swap in orbit',
          len(orbs['SOC']) == 288 and len(soc_u) == 144 and all((t[1], t[0], t[2]) in set(orbs['SOC']) for t in orbs['SOC']),
          f'{len(orbs["SOC"])} / {len(soc_u)}')
    oh = set(orbs['oldH'])
    check('O3 old-H ordered orbit 288; reversed order never in orbit',
          len(oh) == 288 and all((t[1], t[0], t[2]) not in oh for t in oh), str(len(oh)))
    check('O3 old-H (3,15;26) [their representative] lies in the orbit of (12,2;70)', (3, 15, 26) in oh)
    ch = set(orbs['chain'])
    ch_u = set((t[0], frozenset(t[1:3]), t[3]) for t in ch)
    swap_in = all((t[0], t[2], t[1], t[3]) in ch for t in ch)
    check('O3 chain: 576 distinct (centre, {leaves}, k)', len(ch_u) == 576,
          f'ordered {len(ch)}, unordered-leaves {len(ch_u)}, leaf swap in orbit: {swap_in}')
    fo = set(orbs['four'])
    fo_u = set((frozenset(t[:2]), t[2], t[3]) for t in fo)
    swap4 = all((t[1], t[0], t[2], t[3]) in fo for t in fo)
    check('O3 four: ordered 576, ({a,b}, c, d) 288, a<->b swap in orbit', len(fo) == 576 and len(fo_u) == 288 and swap4,
          f'{len(fo)} / {len(fo_u)} / {swap4}')

    cited = {
        'SOC': [(1, 11, 17), (1, 9, 22), (2, 12, 26), (2, 15, 26), (12, 15, 26), (5, 2, 54), (2, 5, 70), (1, 9, 18),
                (2, 39, 46), (2, 8, 69), (2, 15, 56), (1, 12, 35)],
        'oldH': [(3, 15, 26), (23, 15, 26), (12, 2, 70), (39, 2, 70), (3, 5, 54), (23, 5, 54), (9, 2, 69), (10, 8, 69),
                 (53, 8, 69), (29, 12, 35), (2, 12, 35), (9, 5, 68), (11, 12, 26), (15, 10, 30), (8, 10, 30),
                 (33, 2, 46), (5, 2, 46), (9, 2, 56), (10, 15, 56), (53, 15, 56)],
        'chain': [(2, 0, 12, 70), (2, 39, 67, 70), (15, 0, 3, 26), (5, 0, 3, 54), (5, 23, 6, 54), (5, 7, 9, 68),
                  (2, 0, 15, 54), (15, 0, 10, 56), (2, 4, 33, 46), (2, 5, 6, 46), (59, 53, 55, 54), (59, 53, 67, 54)],
        'four': [(0, 6, 10, 54), (6, 0, 5, 34), (0, 13, 3, 26), (0, 7, 12, 40), (7, 0, 1, 68), (7, 67, 53, 68),
                 (0, 67, 2, 40), (67, 41, 39, 40), (7, 21, 9, 35), (4, 6, 2, 28), (4, 7, 3, 22), (0, 6, 1, 70),
                 (6, 67, 39, 70), (67, 55, 53, 54)],
    }
    for fam, lst in cited.items():
        S = set(orbs[fam])
        bad = [t for t in lst if t not in S]
        extra = ''
        if fam == 'chain':   # the chain form is symmetric in its two leaves
            sw = [t for t in bad if (t[0], t[2], t[1], t[3]) in S]
            bad = [t for t in bad if t not in sw]
            extra = f'; in orbit only with the two leaves listed in the other order: {sw}'
        if fam == 'oldH':
            rev = [t for t in lst if (t[1], t[0], t[2]) in S]
            extra = f'; reversed-order hits {rev}'
        check(f'O4 cited {fam} tuples ({len(lst)}) lie in the orbit with the printed roles', not bad, f'misses {bad}{extra}')
    return orbs


# =============================================================================================
SGN_NP = np.array(L.SGN, dtype=float)
PER_NP = np.array(L.PERMS)
FV_NP = np.array(L.F_VALUES, dtype=float)
ROWS = np.arange(4)[None, :]


def O_np(z):
    p = z[:4] / np.linalg.norm(z[:4])
    q = z[4:] / np.linalg.norm(z[4:])
    qb = q * np.array([1, -1, -1, -1])
    return np.array(Lmat(p)) @ np.array(Rmat(qb))


def eval_terms(O, lib=np):
    w = SGN_NP * np.prod(O[ROWS, PER_NP], axis=1)
    return FV_NP @ w


def x_np(O, j):
    I, J = MINOR_IJ[j]
    if len(I) == 1:
        return O[I[0], J[0]]
    return O[I[0], J[0]] * O[I[1], J[1]] - O[I[0], J[1]] * O[I[1], J[0]]


def form(fam, t, Ph, X, lam=2.0, drop=False, xabs=abs):
    if fam == 'SOC':
        i, j, k = t
        prod = X(i) * X(j)
        return Ph(k) + lam * (prod if drop else xabs(prod))
    if fam == 'oldH':
        a, c, k = t
        prod = X(a) * X(c)
        return Ph(c) + Ph(k) + lam * (prod if drop else xabs(prod))
    if fam == 'chain':
        c, l, m, k = t
        first = X(c) * X(l) if drop else xabs(X(c)) * xabs(X(l))
        return Ph(k) + lam * (first + xabs(X(c)) * xabs(X(m)))
    if fam == 'four':
        a, b, c, d = t
        ab = X(a) * X(b) if drop else xabs(X(a) * X(b))
        return Ph(d) + Ph(c) + lam * (ab + xabs(X(a) * X(c)) + xabs(X(b) * X(c)))


def f_np(fam, t, **kw):
    def f(z):
        O = O_np(z)
        ph = eval_terms(O)
        return form(fam, t, lambda j: ph[j], lambda j: x_np(O, j), **kw)
    return f


def mp_value(fam, t, z, dps=40, **kw):
    mpmath.mp.dps = dps
    zz = [mpmath.mpf(float(v)) for v in z]
    np_ = mpmath.sqrt(sum(v * v for v in zz[:4]))
    nq = mpmath.sqrt(sum(v * v for v in zz[4:]))
    p = [v / np_ for v in zz[:4]]
    q = [v / nq for v in zz[4:]]
    O = mm(Lmat(p), Rmat((q[0], -q[1], -q[2], -q[3])), mpmath.mpf(0))
    W = weights(O)
    return form(fam, t, lambda j: Phi(j, W, mpmath.mpf(0)), lambda j: minor(O, j), xabs=abs, **kw)


def multistart(f, n, seed, maxiter=700):
    rng = np.random.default_rng(seed)
    best = (np.inf, None)
    vals = []
    for _ in range(n):
        r = minimize(f, rng.normal(size=8), method='Nelder-Mead',
                     options={'maxiter': maxiter, 'xatol': 1e-10, 'fatol': 1e-14})
        vals.append(r.fun)
        if r.fun < best[0]:
            best = (r.fun, r.x)
    return best, np.array(vals)


def section_attack(fam):
    global GM
    GM = [gmap(g) for g in L.G_ELEMENTS]
    orb = [t for t in orbit(REP[fam]) if t != REP[fam]]
    rng = np.random.default_rng({'SOC': 11, 'oldH': 12, 'chain': 13, 'four': 14}[fam])
    members = [orb[i] for i in rng.choice(len(orb), 30, replace=False)]
    t0 = time.time()
    (bv, bz), vals = multistart(f_np(fam, REP[fam]), 2000, seed=100)
    mpv = mp_value(fam, REP[fam], bz)
    check(f'A {fam} representative {REP[fam]}: 2000 multistarts, min >= -1e-12', bv >= -1e-12,
          f'float min {bv:.3e}, mpmath(40 digits) at minimiser {mpmath.nstr(mpv, 5)}')
    worst = (np.inf, None, None)
    for k, t in enumerate(members):
        (v, z), _ = multistart(f_np(fam, t), 70, seed=200 + k)
        if v < worst[0]:
            worst = (v, t, z)
    mpw = mp_value(fam, worst[1], worst[2])
    check(f'A {fam} 30 other orbit members x 70 multistarts, min >= -1e-12', worst[0] >= -1e-12,
          f'worst member {worst[1]}, float {worst[0]:.3e}, mpmath {mpmath.nstr(mpw, 5)}; members {members}')
    print(f'   time {time.time() - t0:.0f}s', flush=True)


# ---------------------------------------------------------------------------------------------
def rational_unit_quat(v, den=10 ** 6):
    """Exact rational unit quaternion near the unit vector v (inverse stereographic projection)."""
    v = np.asarray(v, float) / np.linalg.norm(v)
    if v[0] < -0.9:
        v = -v
    u = [Fraction(round(x / (1 + v[0]) * den), den) for x in v[1:]]
    n = sum(x * x for x in u)
    return [(1 - n) / (1 + n)] + [2 * x / (1 + n) for x in u]


def exact_value(fam, t, z, **kw):
    p = rational_unit_quat(z[:4])
    q = rational_unit_quat(z[4:])
    O = mm(Lmat(p), Rmat((q[0], -q[1], -q[2], -q[3])), Fraction(0))
    assert all(sum(O[k][i] * O[k][j] for k in range(4)) == (1 if i == j else 0) for i in range(4) for j in range(4))
    W = weights(O)
    return form(fam, t, lambda j: Phi(j, W, Fraction(0)), lambda j: minor(O, j), **kw)


def section_variants():
    global GM
    GM = [gmap(g) for g in L.G_ELEMENTS]
    soc_orb = set(orbit(REP['SOC']))
    # SOC is symmetric in its pair, so "roles swapped" = the pair (1, 9) paired with a wrong facet k
    k_wrong = None
    for k in range(82):
        if L.F_IS_MINOR[k] or (1, 9, k) in soc_orb:
            continue
        (v, z), _ = multistart(f_np('SOC', (1, 9, k)), 20, seed=400 + k)
        if v < -1e-3:
            k_wrong = k
            break
    variants = [
        ('SOC', (1, 9, k_wrong), {}, f'wrong facet: pair (1,9) with k = {k_wrong} (first non-orbit k found negative)'),
        ('SOC', (1, 9, 22), {'lam': 1.8}, 'factor 1.8'),
        ('SOC', (1, 9, 22), {'drop': True}, 'absolute value dropped'),
        ('oldH', (2, 12, 70), {}, 'roles swapped (aux <-> main)'),
        ('oldH', (12, 2, 70), {'lam': 1.8}, 'factor 1.8'),
        ('oldH', (12, 2, 70), {'drop': True}, 'absolute value dropped'),
        ('chain', (0, 2, 12, 70), {}, 'roles swapped (centre <-> leaf)'),
        ('chain', (2, 0, 12, 70), {'lam': 1.8}, 'factor 1.8'),
        ('chain', (2, 0, 12, 70), {'drop': True}, 'absolute value dropped on |x2||x0|'),
        ('four', (0, 10, 6, 54), {}, 'roles swapped (b <-> c)'),
        ('four', (0, 6, 10, 54), {'lam': 1.8}, 'factor 1.8'),
        ('four', (0, 6, 10, 54), {'drop': True}, 'absolute value dropped on |x0 x6|'),
        ('chain', (2, 0, 12, 70), {'lam': 1.0}, 'half compensation (factor 1 instead of 2)'),
    ]
    for k, (fam, t, kw, desc) in enumerate(variants):
        (v, z), _ = multistart(f_np(fam, t, **kw), 300, seed=500 + k, maxiter=1500)
        if fam == 'chain' and kw.get('lam') == 1.8:
            check(f'V chain {t} factor 1.8: still nonnegative (1.8 is NOT a failing variant for chain)', v >= -1e-12,
                  f'float min {v:.3e}')
            continue
        kwx = dict(kw)
        if 'lam' in kwx:
            kwx['lam'] = Fraction(kwx['lam']).limit_denominator(10)
        ex = exact_value(fam, t, z, **kwx)
        check(f'V {fam} {t} {desc}: negative minimum found', v < -1e-6,
              f'float {v:.4e}; EXACT value at nearby rational SO(4) point {float(ex):.4e} (< 0: {ex < 0})')

    # extra multiplier certificate for chain (NUMERIC): in the old-H coordinates
    #   F = Phi70 - 2 O33 x0' x2' - (A11 + B11) x2' x12' >= 0, |O33| <= 1, |A11 + B11| <= 2
    def Fcert(z):
        p = z[:4] / np.linalg.norm(z[:4]); q = z[4:] / np.linalg.norm(z[4:])
        P, Q = p * p, q * q
        D = P - Q
        Kp = (p[0] * q[3] - p[2] * q[1]) ** 2 + (p[1] * q[2] - p[3] * q[0]) ** 2
        phi70 = 2 * Kp - (D[1] - D[3]) * (D[0] - D[2])
        x2 = D[0] + D[2]; x12 = D[0] + D[3]
        x0 = p[0] * q[1] - p[1] * q[0] - p[2] * q[3] + p[3] * q[2]
        O33 = -p[0] * q[1] - p[1] * q[0] + p[2] * q[3] + p[3] * q[2]
        A11 = P[0] + P[1] - P[2] - P[3]; B11 = Q[0] + Q[1] - Q[2] - Q[3]
        return phi70 - 2 * O33 * x0 * x2 - (A11 + B11) * x2 * x12
    (v, z), _ = multistart(Fcert, 1000, seed=900, maxiter=1500)
    check('C chain multiplier form F = Phi70 - 2 O33 x0 x2 - (A11+B11) x2 x12 >= 0 (1000 multistarts)', v >= -1e-12,
          f'min {v:.3e}  [NUMERIC only; implies the chain since |O33| <= 1, |A11+B11| <= 2]')
    # chain compensation threshold: bisection on lam in [1, 1.8] (150 multistarts per step)
    lo, hi = 1.0, 1.8
    for step in range(5):
        mid = (lo + hi) / 2
        (v, z), _ = multistart(f_np('chain', REP['chain'], lam=mid), 150, seed=950 + step, maxiter=1500)
        if v < -1e-9:
            lo = mid
        else:
            hi = mid
    check('C chain: Phi70 + lam |x2|(|x0|+|x12|) has a negative minimum for lam <= lo, none found for lam >= hi',
          True, f'lo = {lo:.3f}, hi = {hi:.3f} (NUMERIC)')
    # nondegenerate equality of the chain atom: penalise |x2| < 0.3
    fc = f_np('chain', REP['chain'])

    def pen(z):
        O = O_np(z)
        return fc(z) + 100 * max(0.0, 0.3 - abs(x_np(O, 2))) ** 2
    (v, z), _ = multistart(pen, 200, seed=901, maxiter=3000)
    O = O_np(z)
    val = fc(z)
    check('C chain attains 0 at a point with |x2| >= 0.29 (equality is not only at degenerate points)',
          abs(val) < 1e-10 and abs(x_np(O, 2)) >= 0.29,
          f'value {val:.2e}, mpmath {mpmath.nstr(mp_value("chain", REP["chain"], z), 4)}, x2 {x_np(O, 2):.4f}, '
          f'x0 {x_np(O, 0):.2e}, x12 {x_np(O, 12):.2e}, Phi70 {eval_terms(O)[70]:.2e}')


if __name__ == '__main__':
    sec = ARGS[0] if ARGS else 'exact'
    print(f'== section {sec} {" ".join(ARGS[1:])}', flush=True)
    if sec == 'exact':
        section_exact()
    elif sec == 'orbits':
        section_orbits()
    elif sec == 'attack':
        section_attack(ARGS[1])
    elif sec == 'variants':
        section_variants()
    print(f'== section {sec} done: {NFAIL} failures', flush=True)
