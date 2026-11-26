"""Independent re-verification of so4-atoms-part2b (label atoms2b-verify).

Own code only. Uses the audit library ../lib/moc4lib.py (paper Appendix D vectors f_j, the
group G of order 1152) and reads supplement JSON files as DATA only (cubic feature basis,
four sign maps, congruence matrices T and V). No supplement code is executed or imported.

Run:  MOC4_SUPP=<supplement root> python -I -B so4_atoms_part2b_verify.py [SECTIONS]
SECTIONS (default all): ROLES SRC GRAM44 GRAM46 GRAM47 ALG NUM

Conventions (re-derived here, not taken from the other auditor):
  * O(p,q) is the matrix of x -> p x q on the quaternions (basis 1,i,j,k), column j = image
    of e_j.  This map S^3 x S^3 -> SO(4) is onto.
  * Gram identity checked:  (|p|^2+|q|^2) * sum_{parity} sum_{4 sign maps S}
        F(Sy)^T M_F F(Sy)  ==  4 * Phi_c(O(p,q))   modulo N = |p|^2 - |q|^2,
    both sides homogeneous of degree 8 in y = (p,q).  On |p| = |q| = 1 this reads
    SOS = 2 Phi_c, so M_F PSD  =>  Phi_c >= 0 on SO(4).  Normal form mod N by division
    (deglex, leading term p0^2; {N} is a Groebner basis).
  * A document Gram M in coordinates g = W F (W from data) corresponds to M_F = W^T M W.
"""
import os, sys, json, re, math, time, itertools, random
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lib'))
import moc4lib as L                      # noqa: E402
import flint                             # noqa: E402
import sympy as sp                       # noqa: E402
import numpy as np                       # noqa: E402

ARGS = sys.argv[1:]
SUPP = os.environ.get('MOC4_SUPP')
if ARGS and os.path.isdir(ARGS[0]):
    SUPP = ARGS.pop(0)
BASE = os.path.join(SUPP, 'proofs', 'manin_transfer_20261021') if SUPP else None
SECTIONS = ARGS or ['ROLES', 'SRC', 'GRAM44', 'GRAM46', 'GRAM47', 'ALG', 'NUM']

NPASS = 0
NFAIL = 0
T0 = time.time()


def check(name, ok, info=''):
    global NPASS, NFAIL
    if ok:
        NPASS += 1
    else:
        NFAIL += 1
    print(('PASS ' if ok else 'FAIL ') + name + ((' :: ' + str(info)) if info != '' else ''), flush=True)
    return ok


def head(s):
    print('\n=== ' + s + '  [t=%.1fs]' % (time.time() - T0), flush=True)


# =====================================================================================
# Group action on the 82 generators; orbit membership with roles
# =====================================================================================
VAL2IDX = {tuple(v): j for j, v in enumerate(L.F_VALUES)}
GMAPS = []
for g in L.G_ELEMENTS:
    t = L.g_index_map(g)
    row = []
    for j in range(82):
        out = [None] * 24
        for s in range(24):
            out[t[s]] = L.F_VALUES[j][s]
        row.append(VAL2IDX[tuple(out)])
    GMAPS.append(row)


def find_role_map(src, dst):
    """Some g (index) with g.f_src[i] = f_dst[i] for all i, else None."""
    for gi, m in enumerate(GMAPS):
        if all(m[a] == b for a, b in zip(src, dst)):
            return gi
    return None


# Atom representatives (paper f_j indices) with roles, from so4-atoms.md (proved families).
REP = {
    'SOC': ((1, 9), 22),          # (i, j; k), pair unordered
    'H': ((12, 2), 70),           # (aux, main; k), ordered
    'chain': ((2, 0, 12), 70),    # (centre; leaf, leaf; k), leaves unordered
    'four': ((0, 6, 10), 54),     # (a, b; c; d), a, b unordered
}


def in_family(fam, roles, k):
    (r, kr) = REP[fam]
    tgt = list(roles) + [k]
    cands = [list(r) + [kr]]
    if fam == 'SOC':
        cands.append([r[1], r[0], kr])
    if fam == 'chain':
        cands.append([r[0], r[2], r[1], kr])
    if fam == 'four':
        cands.append([r[1], r[0], r[2], kr])
    for c in cands:
        if find_role_map(c, tgt) is not None:
            return True
    return False


# =====================================================================================
# Exact polynomial Gram engine
# =====================================================================================
VARS = ('p0', 'p1', 'p2', 'p3', 'q0', 'q1', 'q2', 'q3')
CTX = flint.fmpq_mpoly_ctx.get(VARS, 'deglex')
Y = CTX.gens()
PQ, QQ = Y[:4], Y[4:]
NP = sum(x * x for x in PQ)
NQ = sum(x * x for x in QQ)
NREL = NP - NQ


def red(f):
    return divmod(f, NREL)[1]


def qmul(a, b):
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    return (a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3,
            a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
            a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1,
            a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0)


def O_pxq(p, q, zero, one):
    E = [tuple(one if k == j else zero for k in range(4)) for j in range(4)]
    cols = [qmul(qmul(p, E[j]), q) for j in range(4)]
    return [[cols[j][i] for j in range(4)] for i in range(4)]


def leibniz(vals, O):
    out = 0
    for v, p, s in zip(vals, L.PERMS, L.SGN):
        if v:
            out = out + (s * v) * O[0][p[0]] * O[1][p[1]] * O[2][p[2]] * O[3][p[3]]
    return out


def parse_poly(s):
    """Sum of monomials such as '-p0**2*q3 + 3*p1*q0' (data; parsed, never evaluated)."""
    s = s.replace(' ', '')
    if s[0] not in '+-':
        s = '+' + s
    out = CTX.from_dict({})
    consumed = 0
    for m in re.finditer(r'([+-])([^+-]+)', s):
        coef = Fraction(-1 if m.group(1) == '-' else 1)
        expo = [0] * 8
        for t in m.group(2).replace('**', '^').split('*'):
            v, e = (t.split('^') + ['1'])[:2]
            if v in VARS:
                expo[VARS.index(v)] += int(e)
            else:
                assert re.fullmatch(r'\d+', v), t
                coef *= Fraction(int(v)) ** int(e)
        out += CTX.from_dict({tuple(expo): flint.fmpq(coef.numerator, coef.denominator)})
        consumed += len(m.group(0))
    assert consumed == len(s), s
    return out


def apply_map(f, smap):
    new = [s * Y[i] for (i, s) in smap]        # y'_k = s * y_i for row k = [i, s]
    return f.compose(*new)


class Face:
    def __init__(self, label, cubic_file, gens, letters):
        d = json.load(open(os.path.join(BASE, cubic_file)))
        self.label, self.data, self.gens, self.letters = label, d, list(gens), letters
        self.feat = [[parse_poly(s) for s in d['features'][par]] for par in range(2)]
        n = 8
        assert all(len(f) == n for f in self.feat)
        TF = [[apply_map(f, sm) for f in self.feat[par]] for par in range(2) for sm in d['left_four']]
        self.TF = TF
        self.pairs = [(a, b) for a in range(n) for b in range(a, n)]
        cols = []
        for (a, b) in self.pairs:
            s = CTX.from_dict({})
            for vec in TF:
                s += vec[a] * vec[b]
            cols.append(red((NP + NQ) * (s if a == b else 2 * s)))
        zero = CTX.from_dict({})
        O = O_pxq(PQ, QQ, zero, zero + 1)
        self.O = O
        costs = [red(4 * leibniz(L.F_VALUES[f], O)) for f in self.gens]
        mons = sorted(set().union(*[set(c.monoms()) for c in cols + costs]))
        idx = {m: i for i, m in enumerate(mons)}
        nc = len(cols) + len(costs)
        A = flint.fmpq_mat(len(mons), nc)
        for j, c in enumerate(cols):
            for m, v in zip(c.monoms(), c.coeffs()):
                A[idx[m], j] = v
        for j, c in enumerate(costs):
            for m, v in zip(c.monoms(), c.coeffs()):
                A[idx[m], len(cols) + j] = -v
        R, rk = A.rref()
        self.rank = rk
        self.rel = [[sp.Rational(int(R[i, j].p), int(R[i, j].q)) for j in range(nc)] for i in range(rk)]

        def sub_rank(js):
            B = flint.fmpq_mat(len(mons), len(js))
            for i in range(len(mons)):
                for jj, j in enumerate(js):
                    B[i, jj] = A[i, j]
            return B.rref()[1]
        self.gram_rank = sub_rank(list(range(len(cols))))
        self.cost_rank = sub_rank(list(range(len(cols), nc)))
        self.nmons = len(mons)
        self.W = None

    def set_W(self, W, gorder):
        self.W = W
        self.gorder = gorder          # letter of each g coordinate

    def residuals(self, M, coef):
        """M: 8x8 sympy Gram in g coordinates (rows in self.gorder); coef: dict letter->value."""
        MF = self.W.T * M * self.W
        vec = [MF[a, b] for (a, b) in self.pairs] + [coef.get(l, 0) for l in self.letters]
        out = []
        for row in self.rel:
            out.append(sp.cancel(sp.together(sum(c * z for c, z in zip(row, vec) if c != 0))))
        return [z for z in out if z != 0]

    def gram_from_vectors(self, X, Gb):
        """X: dict letter -> coefficient row (list) on a base with Gram Gb."""
        rows = [X.get(l, [0] * Gb.shape[0]) for l in self.gorder]
        Xm = sp.Matrix(rows)
        return sp.simplify(Xm * Gb * Xm.T)

    def exact_point_check(self, M, coef, subsd, seed):
        """Exact evaluation at rational unit quaternions: SOS(M_F) == 2 Phi_c(O)."""
        rng = random.Random(seed)

        def runit():
            a = [Fraction(rng.randint(-9, 9), rng.randint(1, 9)) for _ in range(3)]
            n2 = sum(x * x for x in a)
            return [(1 - n2) / (1 + n2)] + [2 * x / (1 + n2) for x in a]
        p, q = runit(), runit()
        y = p + q
        MF = (self.W.T * M * self.W).subs(subsd)
        sos = Fraction(0)
        for vec in self.TF:
            fv = [ev(f, y) for f in vec]
            for a in range(8):
                for b in range(8):
                    mab = MF[a, b]
                    if mab != 0:
                        mq = sp.Rational(sp.cancel(mab))
                        sos += Fraction(int(mq.p), int(mq.q)) * fv[a] * fv[b]
        O = O_pxq(p, q, Fraction(0), Fraction(1))
        cv = [0] * 24
        for l, f in zip(self.letters, self.gens):
            clv = sp.Rational(sp.cancel(sp.sympify(coef.get(l, 0)).subs(subsd)))
            cl = Fraction(int(clv.p), int(clv.q))
            for s in range(24):
                cv[s] += cl * L.F_VALUES[f][s]
        phi = leibniz(cv, O)
        return sos == 2 * phi, (float(sos), float(2 * phi))


def ev(f, y):
    tot = Fraction(0)
    for m, c in zip(f.monoms(), f.coeffs()):
        t = Fraction(int(c.p), int(c.q))
        for e, v in zip(m, y):
            if e:
                t *= v ** int(e)
        tot += t
    return tot


def load_T(fname, key='T'):
    d = json.load(open(os.path.join(BASE, fname)))
    return sp.Matrix([[sp.Rational(x) for x in r] for r in d[key]]), d


def probe_diagonal(face):
    """For each g coordinate i, the cost (as letters) represented by the unit Gram E_ii, or None."""
    cs = sp.symbols('c0:8')
    out = []
    for i in range(8):
        E = sp.zeros(8)
        E[i, i] = 1
        MF = face.W.T * E * face.W
        vec = [MF[a, b] for (a, b) in face.pairs] + list(cs)
        eqs = [sum(c * z for c, z in zip(row, vec)) for row in face.rel]
        sol = sp.solve(eqs, cs, dict=True)
        if not sol:
            out.append(None)
        else:
            v = [sol[0].get(c, c) for c in cs]
            out.append('+'.join(('%s*' % x if x != 1 else '') + face.letters[j] for j, x in enumerate(v) if x != 0))
    return out


def face_header(face):
    check(face.label + ': costs lie in the span of the Gram columns; ranks (gram, cost, total) = (18, 8, 18)',
          (face.gram_rank, face.cost_rank, face.rank) == (18, 8, 18),
          '(%d, %d, %d), %d monomials, solution space dim %d' % (face.gram_rank, face.cost_rank, face.rank,
                                                                face.nmons, 44 - face.rank))


def gram_claim(face, name, M, coef, subsd=None, seed=0):
    res = face.residuals(M, coef)
    ok = check(face.label + ' ' + name + ': Leibniz sum == quadratic form (all relations, symbolic)', not res,
               'nonzero residuals: %d' % len(res))
    if subsd is not None:
        ok2, vals = face.exact_point_check(M, coef, subsd, seed)
        check(face.label + ' ' + name + ': exact rational-point check SOS = 2*Phi', ok2, vals)
    return ok


# =====================================================================================
# Section ROLES: orbit memberships and transports
# =====================================================================================
def sec_roles():
    head('ROLES: atom memberships with the roles used (G acting on f_j)')
    check('old-H representatives (12,2;70) and (3,15;26) in one orbit', find_role_map([12, 2, 70], [3, 15, 26]) is not None)
    atoms = [
        # G1 T44
        ('T44 SOC(1,53;18)', 'SOC', (1, 53), 18), ('T44 H(15,53;18)', 'H', (15, 53), 18),
        ('T44 H(8,53;18)', 'H', (8, 53), 18),
        # G2 T45
        ('T45 (2) SOC(2,59;70)', 'SOC', (2, 59), 70), ('T45 (3) chain(2;0,12;70)', 'chain', (2, 0, 12), 70),
        ('T45 (4) chain(59;0,1;70)', 'chain', (59, 0, 1), 70), ('T45 (5) chain(59;67,72;70)', 'chain', (59, 67, 72), 70),
        ('T45 (6) H(1,59;70)', 'H', (1, 59), 70), ('T45 (6) H(72,59;70)', 'H', (72, 59), 70),
        # G3 T46
        ('T46 four(4,6;2;28)', 'four', (4, 6, 2), 28),
        # G4 T47
        ('T47 H(9,1;35)', 'H', (9, 1), 35), ('T47 SOC(45,1;35)', 'SOC', (45, 1), 35),
        ('T47 chain(1;7,9;35)', 'chain', (1, 7, 9), 35), ('T47 four(0,7;5;35)', 'four', (0, 7, 5), 35),
        # G5
        ('T28 (3) H(45,1;34)', 'H', (45, 1), 34), ('T28 (3) chain(1;0,45;34)', 'chain', (1, 0, 45), 34),
        ('T28 (3) four(6,0;5;34)', 'four', (6, 0, 5), 34),
        ('T38 H(2,5;40)', 'H', (2, 5), 40), ('T38 four(0,7;12;40)', 'four', (0, 7, 12), 40),
        ('T38 four(0,67;2;40)', 'four', (0, 67, 2), 40), ('T38 SOC(2,5;70)', 'SOC', (2, 5), 70),
        ('T27 four(67,41;39;40)', 'four', (67, 41, 39), 40), ('T27 four(0,67;2;40)', 'four', (0, 67, 2), 40),
        ('T27 H(2,5;40)', 'H', (2, 5), 40),
        ('T42 SOC(5,2;54)', 'SOC', (5, 2), 54),
        ('T36 (6) SOC(45,1;35)', 'SOC', (45, 1), 35), ('T36 (6) H(5,45;35)', 'H', (5, 45), 35),
        ('T36 (6) four(7,0;5;35)', 'four', (7, 0, 5), 35), ('T36 SOC(5,45;40)', 'SOC', (5, 45), 40),
        ('T48 H(2,5;40)', 'H', (2, 5), 40), ('T48 H(29,5;40)', 'H', (29, 5), 40),
        ('T48 chain(5;0,2;40)', 'chain', (5, 0, 2), 40), ('T48 chain(5;7,29;40)', 'chain', (5, 7, 29), 40),
        ('T48 chain(5;41,29;40)', 'chain', (5, 41, 29), 40),
    ]
    # H/chain allocation lemma role tuples (b, c, p, q, a, k): H(b,a;k), H(c,a;k), chain(a;p,b;k), chain(a;c,q;k)
    HC = {'T28 Lemma 2': (33, 5, 4, 6, 2, 46), 'T38 gen.70': (12, 39, 0, 67, 2, 70),
          'T42': (3, 23, 0, 6, 5, 54), 'T36 (10)': (2, 29, 0, 7, 5, 40)}
    for nm, (b, c, p, q, a, k) in HC.items():
        atoms += [(nm + ' H(b=%d,a=%d;%d)' % (b, a, k), 'H', (b, a), k), (nm + ' H(c=%d,a=%d;%d)' % (c, a, k), 'H', (c, a), k),
                  (nm + ' chain(a;p,b)', 'chain', (a, p, b), k), (nm + ' chain(a;c,q)', 'chain', (a, c, q), k)]
    for nm, fam, roles, k in atoms:
        check('membership ' + nm, in_family(fam, roles, k))
    # reversed old-H roles must fail (role order meaningful)
    check('control: reversed old-H roles (53,15;18) NOT in orbit', not in_family('H', (53, 15), 18))
    check('control: reversed old-H roles (5,29;40) NOT in orbit', not in_family('H', (5, 29), 40))

    head('ROLES: transports')
    tr = [
        ('T44 J*: (1,53,12,7,18) -> (53,1,8,7,18)', [1, 53, 12, 7, 18], [53, 1, 8, 7, 18]),
        ('T46 core -> T47: [1,59,6,27,28] -> [45,5,7,0,35]', [1, 59, 6, 27, 28], [45, 5, 7, 0, 35]),
        ('T24 core -> T46: [59,53,55,67,54] -> [1,59,6,27,28]', [59, 53, 55, 67, 54], [1, 59, 6, 27, 28]),
        ('six-ray T46 -> T47: [27,1,59,5,28,6] -> [0,45,5,9,35,7]', [27, 1, 59, 5, 28, 6], [0, 45, 5, 9, 35, 7]),
        ('T24 core -> T38: [59,53,55,67,54] -> [5,2,67,0,40]', [59, 53, 55, 67, 54], [5, 2, 67, 0, 40]),
        ('T24 core -> T27: [59,53,55,67,54] -> [5,2,0,67,40]', [59, 53, 55, 67, 54], [5, 2, 0, 67, 40]),
        ('T24 core -> T36: (67,55,53,54,59,2) -> (7,0,5,35,45,1)', [67, 55, 53, 54, 59, 2], [7, 0, 5, 35, 45, 1]),
        ('T24 core -> T47: (59,53,55,67,54,2) -> (45,5,7,0,35,1)', [59, 53, 55, 67, 54, 2], [45, 5, 7, 0, 35, 1]),
        ('T28 face -> T36: (4,2,33,5,6,46,45) -> (0,5,2,29,7,40,39)', [4, 2, 33, 5, 6, 46, 45], [0, 5, 2, 29, 7, 40, 39]),
        ('T36 stabiliser (7,0,2,29,45,1,5,35) -> (0,7,29,2,45,1,5,35)', [7, 0, 2, 29, 45, 1, 5, 35], [0, 7, 29, 2, 45, 1, 5, 35]),
    ]
    for nm, s, d in tr:
        check('transport ' + nm, find_role_map(s, d) is not None)
    # T42: some g preserving the contact set exchanges 22 and 54
    T42 = frozenset([1, 2, 4, 6, 7, 8])
    ok = any(GMAPS[gi][22] == 54 and GMAPS[gi][54] == 22 and L.g_act_set(L.G_ELEMENTS[gi], T42) == T42 for gi in range(1152))
    check('T42: a contact-set stabiliser exchanges f22 and f54', ok)
    T27 = frozenset([2, 3, 5, 6, 9, 12])
    ok = any(GMAPS[gi][40] == 68 and GMAPS[gi][68] == 40 and L.g_act_set(L.G_ELEMENTS[gi], T27) == T27 for gi in range(1152))
    check('T27: a contact-set stabiliser exchanges f40 and f68', ok)
    T36 = frozenset([0, 5, 6, 7, 9, 12])
    ok = any(all(GMAPS[gi][a] == b for a, b in zip([7, 0, 2, 29, 45, 1, 5, 35], [0, 7, 29, 2, 45, 1, 5, 35]))
             and L.g_act_set(L.G_ELEMENTS[gi], T36) == T36 for gi in range(1152))
    check('T36: the (a,b,q,p) exchange is a contact-set stabiliser', ok)
    # T45 zero faces
    T45 = [3, 4, 5, 9, 10, 12]
    for nm, src, img, extra in [('a=0: T33', [2, 3, 5, 9, 10, 12], [3, 5, 9, 10, 12, 13], 13),
                                ('c=0: T6', [0, 1, 2, 4, 7, 8], [2, 3, 4, 5, 9, 10], 2),
                                ('r=0: T14', [1, 2, 3, 4, 6, 7], [3, 4, 5, 9, 10, 16], 16)]:
        g = L.find_g(src, img)
        check('T45 zero face ' + nm + ' maps onto %s inside T45+{%d}' % (img, extra),
              g is not None and set(img) <= set(T45) | {extra})


# =====================================================================================
# Section SRC: source identities, generator lists, contact nullspaces
# =====================================================================================
def vec_of(combo):
    v = [0] * 24
    for c, f in combo:
        for s in range(24):
            v[s] += c * L.F_VALUES[f][s]
    return v


def sec_src():
    head('SRC: face generator lists and response identities (values at 24 permutations; P is injective)')
    for nm, T, gens in [('T28', [0, 4, 5, 6, 7, 8], [0, 1, 2, 3, 4, 5, 6, 33, 34, 45, 46]),
                        ('T42', [1, 2, 4, 6, 7, 8], [0, 1, 2, 3, 4, 5, 6, 21, 22, 23, 53, 54]),
                        ('T38', [2, 3, 5, 7, 9, 12], [0, 1, 2, 5, 7, 12, 39, 40, 67, 70]),
                        ('T27', [2, 3, 5, 6, 9, 12], [0, 1, 2, 5, 7, 9, 39, 40, 41, 67, 68]),
                        ('T44', [1, 2, 3, 7, 11, 12], None), ('T46', [1, 3, 5, 6, 7, 8], None),
                        ('T47', [0, 4, 5, 6, 9, 12], None), ('T48', [0, 3, 5, 6, 9, 12], None),
                        ('T45', [3, 4, 5, 9, 10, 12], None), ('T36', [0, 5, 6, 7, 9, 12], None)]:
        G = L.generators(T)
        nonminor = [j for j in G if not L.F_IS_MINOR[j]]
        if gens is not None:
            check(nm + ' generators of contact set %s' % T, G == gens, G)
        else:
            print('INFO %s generators %s nonminor %s' % (nm, G, nonminor))
    faces = {'T44': [2, 1, 53, 12, 8, 18, 7, 15], 'T46': [27, 1, 59, 5, 2, 28, 6, 4], 'T47': [2, 9, 1, 45, 5, 35, 0, 7],
             'T48': [7, 5, 29, 39, 2, 40, 41, 0], 'T45': [0, 1, 2, 12, 59, 67, 70, 72],
             'T28 g34': [3, 6, 1, 33, 45, 34, 0, 5], 'T28 g46': [4, 3, 2, 33, 5, 6, 46, 45],
             'T38 g40': [39, 12, 5, 7, 2, 40, 67, 0], 'T38 g70': [0, 1, 2, 12, 39, 67, 70, 5],
             'T42': [0, 53, 5, 3, 23, 6, 54, 2], 'T27': [7, 39, 5, 41, 2, 40, 0, 67],
             'T36 (2)': [2, 29, 1, 45, 5, 35, 0, 7], 'T36 (3)': [0, 39, 5, 2, 29, 7, 40, 45]}
    contact = {'T44': [1, 2, 3, 7, 11, 12], 'T46': [1, 3, 5, 6, 7, 8], 'T47': [0, 4, 5, 6, 9, 12], 'T48': [0, 3, 5, 6, 9, 12],
               'T45': [3, 4, 5, 9, 10, 12], 'T28 g34': [0, 4, 5, 6, 7, 8], 'T28 g46': [0, 4, 5, 6, 7, 8],
               'T38 g40': [2, 3, 5, 7, 9, 12], 'T38 g70': [2, 3, 5, 7, 9, 12], 'T42': [1, 2, 4, 6, 7, 8],
               'T27': [2, 3, 5, 6, 9, 12], 'T36 (2)': [0, 5, 6, 7, 9, 12], 'T36 (3)': [0, 5, 6, 7, 9, 12]}
    for nm, gs in faces.items():
        M = [list(L.F_VECS[j]) for j in gs]
        rk = L.rank(M)
        vanish = all(L.F_VALUES[j][s] == 0 for j in gs for s in contact[nm])
        nmin = [j for j in gs if not L.F_IS_MINOR[j]]
        check('%s face: 8 generators, rank 8, vanish on contacts, one non-minor %s' % (nm, nmin),
              len(gs) == 8 and rk == 8 and vanish and len(nmin) == 1)
    ids = [
        ('T28', [(1, 34), (1, 46)], [(1, 3), (2, 6), (1, 33), (1, 45)]),
        ('T28', [(1, 2), (1, 34)], [(1, 0), (1, 6), (1, 33)]),
        ('T28', [(1, 4), (1, 34)], [(1, 1), (1, 3), (1, 6), (1, 33)]),
        ('T28', [(1, 0), (1, 46)], [(1, 2), (1, 3), (1, 6), (1, 45)]),
        ('T28', [(1, 1), (1, 46)], [(1, 4), (1, 6), (1, 45)]),
        ('T42', [(1, 22), (1, 54)], [(1, 0), (1, 4), (1, 5), (1, 23), (1, 53)]),
        ('T42', [(1, 54), (1, 1)], [(1, 0), (1, 6), (1, 53)]),
        ('T42', [(1, 54), (1, 4)], [(1, 2), (1, 3), (1, 6), (1, 53)]),
        ('T42', [(1, 54), (1, 21)], [(1, 0), (1, 2), (1, 23), (1, 53)]),
        ('T38', [(1, 40), (1, 70)], [(1, 0), (1, 5), (1, 12), (1, 39), (1, 67)]),
        ('T38', [(1, 1), (1, 40)], [(1, 0), (1, 7), (1, 39)]),
        ('T38', [(1, 7), (1, 70)], [(1, 1), (1, 5), (1, 12), (1, 67)]),
        ('T36', [(1, 35), (1, 40)], [(1, 0), (1, 7), (1, 2), (1, 29), (1, 45)]),
        ('T36', [(1, 35), (1, 39)], [(1, 1), (1, 2), (1, 29), (1, 45)]),
        ('T36', [(1, 40), (1, 1)], [(1, 0), (1, 7), (1, 39)]),
        ('T27', [(1, 40), (1, 68)], [(1, 0), (1, 7), (1, 41), (1, 67)]),
        ('T27', [(1, 40), (1, 1)], [(1, 0), (1, 7), (1, 39)]),
        ('T27', [(1, 40), (1, 9)], [(1, 2), (1, 7), (1, 41)]),
    ]
    for nm, lhs, rhs in ids:
        check('%s identity %s = %s' % (nm, '+'.join('f%d' % f for _, f in lhs), '+'.join(('%d' % c if c != 1 else '') + 'f%d' % f for c, f in rhs)),
              vec_of(lhs) == vec_of(rhs))
    # T36 contact nullspace and decompositions (1)-(3)
    T36 = [0, 5, 6, 7, 9, 12]
    Nb = L.nullspace([L.P[s] for s in T36])
    check('T36 contact nullspace dimension 8', len(Nb) == 8, len(Nb))
    ok = [True] * 4
    for c in Nb:
        g = L.values(c)
        d35, d40 = g[3] - g[1], g[14] - g[16]
        ok[0] &= (d35 + d40 == -g[17])
        e1 = vec_of([(g[18], 0), (-d40, 1), (g[13], 2), (g[11], 5), (g[14], 7), (g[2], 29), (-d35, 39), (g[3], 45)])
        e2 = vec_of([(g[15], 2), (g[4], 29), (g[17], 1), (g[1], 45), (g[11], 5), (d35, 35), (g[18], 0), (g[14], 7)])
        e3 = vec_of([(g[23], 0), (g[17], 39), (g[11], 5), (g[13], 2), (g[2], 29), (g[16], 7), (d40, 40), (g[3], 45)])
        ok[1] &= (e1 == g)
        ok[2] &= (e2 == g)
        ok[3] &= (e3 == g)
    check('T36: delta35 + delta40 = -g17 on the nullspace', ok[0])
    check('T36: decomposition (1) holds on the nullspace', ok[1])
    check('T36: decomposition (2) holds on the nullspace', ok[2])
    check('T36: decomposition (3) holds on the nullspace', ok[3])
    # T45 (1)
    T45 = [3, 4, 5, 9, 10, 12]
    Nb = L.nullspace([L.P[s] for s in T45])
    gens = [0, 1, 2, 12, 59, 67, 70, 72]
    ev_at = [18, 16, 13, 8, 2, 14, 11, 7]
    mat = [[L.F_VALUES[f][s] for f in gens] for s in ev_at]
    check('T45 (1): nullspace dim 8; generators vanish on contacts; values at (18,16,13,8,2,14,11,7) = identity',
          len(Nb) == 8 and all(L.F_VALUES[f][s] == 0 for f in gens for s in T45) and mat == [[int(i == j) for j in range(8)] for i in range(8)])
    check('T45 (1): f70 is the only non-minor', [f for f in gens if not L.F_IS_MINOR[f]] == [70])


# =====================================================================================
# Gram sections
# =====================================================================================
def M_from(face, X, Gb):
    rows = [X.get(l, [0] * Gb.shape[0]) for l in face.gorder]
    Xm = sp.Matrix(rows)
    return Xm * Gb * Xm.T


def sec_gram44():
    head('GRAM44: historical T44 face q f2 + p f1 + j f53 + c f12 + e f8 + d f18 + b f7 + a f15')
    F = Face('T44', 'blind_moc4_T44_cubic_reduction_20261022.json', [2, 1, 53, 12, 8, 18, 7, 15], 'qpjcedba')
    face_header(F)
    T, _ = load_T('blind_moc4_T44_bad18_budget_PSD8_20261022.json')
    F.set_W(T.inv(), 'jceabqpd')
    pr = probe_diagonal(F)
    print('INFO T44 unit Gram E_ii -> cost (g order J,C,E,A,B,Q,P,D):', pr)
    r, k, v, w = sp.symbols('r k v w')
    x = (v * (2 * r - k) - 1) / (r - k)
    y = k * (1 - r * v) / (w * (r - k))
    X = {'j': [r * v, r, 0], 'c': [k * v, k, 0], 'e': [0, v, 0], 'a': [0, v, 0], 'b': [0, 0, w],
         'q': [0, 0, 0], 'p': [0, x, y], 'd': [1, -v, 0]}
    M = M_from(F, X, sp.eye(3))
    coef = {'q': 0, 'p': x ** 2 + y ** 2, 'j': r ** 2 * (1 + v ** 2), 'c': k ** 2 * (1 + v ** 2), 'e': v ** 2,
            'd': 1, 'b': w ** 2, 'a': v ** 2}
    gram_claim(F, 'Prop. 1 (vectors J,C,E,A,B,Q,P,D in R^3)', M, coef,
               {r: sp.Rational(3, 2), k: sp.Rational(-2, 3), v: sp.Rational(1, 3), w: sp.Rational(5, 4)}, 1)
    # control: a wrong coefficient (d-diagonal read as d instead of d + a) must fail
    bad = dict(coef)
    bad['d'] = 1 + v ** 2
    check('T44 control: Prop. 1 Gram with d read off as |D|^2 = 1+v^2 is rejected', bool(F.residuals(M, bad)))
    # b = 0 direct Gram
    J = sp.symbols('J', positive=True)
    X0 = {'j': [0, J, 0], 'p': [0, -1 / J, 0], 'd': [1, 0, 0]}
    M0 = M_from(F, X0, sp.eye(3))
    gram_claim(F, 'Prop. 2 b=0 direct Gram (p=1/j)', M0, {'j': J ** 2, 'p': 1 / J ** 2, 'd': 1}, {J: sp.Rational(7, 5)}, 2)
    return F


def sec_gram46():
    head('GRAM46: historical T46 face q f27 + p f1 + j f59 + c f5 + e f2 + d f28 + b f6 + a f4')
    F = Face('T46', 'final_referee_T46_cubic_reduction_20261022.json', [27, 1, 59, 5, 2, 28, 6, 4], 'qpjcedba')
    face_header(F)
    T, _ = load_T('final_referee_T46_budget_PSD8_20261022.json')
    F.set_W(T.inv(), 'qpjcedba')
    print('INFO T46 unit Gram E_ii -> cost (g order q,p,j,c,e,d,b,a):', probe_diagonal(F))
    q, b, be, t, j, v = sp.symbols('q b beta t j v')
    # core Thm 1: base (Q, B, J0, D)
    G0 = sp.Matrix([[q, be, -t, 0], [be, b, -t, 0], [-t, -t, j - t ** 2, 0], [0, 0, 0, 1]])
    s = q + b + 2 * be
    h = 1 + be
    Dl = s * j - (s + 4) * t ** 2
    lam = (-h * (j - t ** 2) + 2 * t * (t - 1)) / Dl
    mu = (-2 * t * h + s * (t - 1)) / Dl
    X = {'q': [1, 0, 0, 0], 'b': [0, 1, 0, 0], 'j': [0, 0, 1, t], 'd': [0, 0, 0, 1], 'p': [lam, lam, mu, 0]}
    M = M_from(F, X, G0)
    p0 = ((j - t ** 2) * h ** 2 - 4 * t * h * (t - 1) + s * (t - 1) ** 2) / Dl
    check('T46 core: |P|^2 equals threshold (1)', sp.cancel(M[1, 1] - p0) == 0)
    sub = {q: sp.Rational(3, 2), b: sp.Rational(5, 7), be: sp.Rational(-1, 3), t: sp.Rational(1, 5), j: sp.Rational(2, 1)}
    gram_claim(F, 'core coherent Thm 1', M, {'q': q, 'p': p0, 'j': j, 'd': 1, 'b': b}, sub, 3)
    # PD criterion: leading minors
    l = q * b - be ** 2
    m3 = sp.cancel(G0[:3, :3].det())
    check('T46 core: det of leaf+J0 block = l*(j - (1+(q+b-2beta)/l) t^2)', sp.cancel(m3 - l * (j - (1 + (q + b - 2 * be) / l) * t ** 2)) == 0)
    Gm = sp.Matrix([[s, -2 * t], [-2 * t, j - t ** 2]])
    check('T46 core: Gram of (Q+B, J0) has det Delta', sp.cancel(Gm.det() - Dl) == 0)
    # six-ray Thm 1: base (Q, B, J0, Z)
    G6 = sp.Matrix([[q, be, t, 0], [be, b, v, 0], [t, v, j - t ** 2, 0], [0, 0, 0, 1]])
    D6 = s * (j - t ** 2) - (t + v) ** 2
    lam6 = (-h * (j - t ** 2) + (t + v) * (1 + v)) / D6
    mu6 = (h * (t + v) - s * (1 + v)) / D6
    kc = (t - v) / (j - t ** 2 - v ** 2 / b)
    X6 = {'q': [1, 0, 0, 0], 'b': [0, 1, 0, 0], 'j': [0, 0, 1, t], 'p': [lam6, lam6, mu6, 0],
          'c': [0, -kc * v / b, kc, 0], 'd': [0, -kc * v / b, kc, -1]}
    M6 = M_from(F, X6, G6)
    ps = ((j - t ** 2) * h ** 2 - 2 * h * (t + v) * (1 + v) + s * (1 + v) ** 2) / D6
    cs = (t - v) ** 2 / (j - t ** 2 - v ** 2 / b)
    check('T46 six-ray: |P|^2 = p* and |C|^2 = c*', sp.cancel(M6[1, 1] - ps) == 0 and sp.cancel(M6[3, 3] - cs) == 0)
    sub6 = {q: sp.Rational(3, 2), b: sp.Rational(5, 7), be: sp.Rational(-1, 3), t: sp.Rational(1, 5), v: sp.Rational(-2, 7), j: 2}
    gram_claim(F, 'six-ray Thm 1', M6, {'q': q, 'p': ps, 'j': j, 'c': cs, 'd': 1, 'b': b}, sub6, 4)
    m3 = sp.cancel(G6[:3, :3].det())
    check('T46 six-ray: det of (Q,B,J0) block = l*(j - t^2) - (b t^2 - 2 beta t v + q v^2)',
          sp.cancel(m3 - (l * (j - t ** 2) - (b * t ** 2 - 2 * be * t * v + q * v ** 2))) == 0)
    # singular leaf: q = u^2, b = w^2, beta = -u w, v = -r t/q (r = u w)
    u, w = sp.symbols('u w', positive=True)
    rr = u * w
    subsl = {q: u ** 2, b: w ** 2, be: -rr, v: -rr * t / u ** 2}
    Msl = M6.subs(subsl)
    kap = 1 + 1 / u ** 2
    rho = 1 + rr / u ** 2
    eta = (1 - w ** 2) / (u ** 2 - rr)
    kk = u ** 2 * (1 - rr) ** 2 / (u ** 2 - rr) ** 2
    psl = kk + (1 - eta * t) ** 2 / (j - kap * t ** 2)
    csl = rho ** 2 * t ** 2 / (j - kap * t ** 2)
    check('T46 singular leaf: p* = k + (1-eta t)^2/(j - kappa t^2), c* = rho^2 t^2/(j - kappa t^2)',
          sp.cancel(Msl[1, 1] - psl) == 0 and sp.cancel(Msl[3, 3] - csl) == 0)
    check('T46 singular leaf: Delta = (q+b-2r)(j - kappa t^2) after substitution (so Delta > 0 when q != b)',
          sp.cancel(D6.subs(subsl) - (u - w) ** 2 * (j - kap * t ** 2)) == 0)
    check('T46 singular leaf: j - t^2 - v^2/b = j - kappa t^2',
          sp.cancel((j - t ** 2 - v ** 2 / b).subs(subsl) - (j - kap * t ** 2)) == 0)
    Gsl = G6.subs(subsl)
    # PSD: leaf block rank one with (t,v) in its range; generalized Schur complement j - kappa t^2
    check('T46 singular leaf: leaf block singular, (t,v) = (t/u)(u,-w) in its range',
          sp.simplify(Gsl[:2, :2].det()) == 0 and sp.simplify(Gsl[0, 2] * (-w) - Gsl[1, 2] * u) == 0)
    gram_claim(F, 'singular-leaf six-ray (beta=-r, v=-rt/q)', Msl, {'q': u ** 2, 'p': psl, 'j': j, 'c': csl, 'd': 1, 'b': w ** 2},
               {u: sp.Rational(3, 2), w: sp.Rational(1, 2), t: sp.Rational(1, 5), j: 2}, 5)
    return F


def sec_gram47():
    head('GRAM47: historical T47 face q f2 + p f9 + j f1 + c f45 + e f5 + d f35 + b f0 + a f7')
    F = Face('T47', 'independent_moc4_T47_cubic_reduction_20261022.json', [2, 9, 1, 45, 5, 35, 0, 7], 'qpjcedba')
    face_header(F)
    T, _ = load_T('independent_moc4_T47_budget_PSD8_20261022.json')
    nj = json.load(open(os.path.join(BASE, 'independent_moc4_T47_native_PSD8_20261022.json')))
    V = sp.Matrix(nj['V'])
    check('T47: det V = +-1 and det T != 0 (W = V^-1 T^-1 invertible)', abs(V.det()) == 1 and T.det() != 0, (V.det(), T.det()))
    for nm, W in [('T^-1', T.inv()), ('V^-1 T^-1', V.inv() * T.inv())]:
        F.set_W(W, 'qpjcedba')
        print('INFO T47 W = %s: unit Gram E_ii -> cost:' % nm, probe_diagonal(F))
    F.set_W(V.inv() * T.inv(), 'qpjcedba')
    x, y, z, w, Lp, j, c, e = sp.symbols('x y z w L j c e')
    idx = {l: i for i, l in enumerate('qpjcedba')}
    M = sp.zeros(8)

    def put(a, b, val):
        M[idx[a], idx[b]] = val
        M[idx[b], idx[a]] = val
    for l, val in zip('qpjcedba', [0, Lp, j, c, e, 1, Lp, 0]):
        put(l, l, val)
    put('p', 'b', Lp)
    put('p', 'j', y); put('b', 'j', y); put('j', 'd', -y)
    put('p', 'c', -x); put('b', 'c', -x); put('c', 'e', -x)
    put('j', 'c', 1 - x + y)
    put('p', 'e', z); put('b', 'e', z); put('e', 'd', -z)
    put('j', 'e', w)
    sub = {x: sp.Rational(1, 3), y: sp.Rational(-1, 4), z: sp.Rational(2, 5), w: sp.Rational(1, 7), Lp: sp.Rational(3, 2),
           j: 2, c: sp.Rational(5, 3), e: sp.Rational(7, 4)}
    gram_claim(F, 'PCBE Thm 1 matrix (identically in x,y,z,w,L,j,c,e)', M, {'p': Lp, 'j': j, 'c': c, 'e': e, 'd': 1, 'b': Lp}, sub, 6)
    # same matrix in the budget coordinates W = T^-1 must NOT satisfy the relations (the coordinate system matters)
    F.set_W(T.inv(), 'qpjcedba')
    nbad = len(F.residuals(M, {'p': Lp, 'j': j, 'c': c, 'e': e, 'd': 1, 'b': Lp}))
    print('INFO PCBE matrix read in budget coordinates T^-1: %d of 18 relations violated' % nbad)
    F.set_W(V.inv() * T.inv(), 'qpjcedba')
    # Schur complement over {B, D} (P removed as a duplicate of B)
    keep = [idx[l] for l in 'jec']
    elim = [idx['b'], idx['d']]
    Mr = M.extract(keep + elim, keep + elim)
    S = sp.simplify(Mr[:3, :3] - Mr[:3, 3:] * Mr[3:, 3:].inv() * Mr[3:, :3])
    kap = 1 + 1 / Lp
    Sclaim = sp.Matrix([[j - kap * y ** 2, w - kap * y * z, 1 - x + (1 + x / Lp) * y],
                        [w - kap * y * z, e - kap * z ** 2, -x + x * z / Lp],
                        [1 - x + (1 + x / Lp) * y, -x + x * z / Lp, c - x ** 2 / Lp]])
    check('PCBE: vectors P and B coincide (equal rows, M_PP = M_BB = M_PB = L)',
          all(M[idx['p'], k] == M[idx['b'], k] for k in range(8) if k not in (idx['p'], idx['b']))
          and M[idx['p'], idx['p']] == M[idx['b'], idx['b']] == M[idx['p'], idx['b']] == Lp)
    check('PCBE: Schur complement over {B,D} equals the displayed 3x3 matrix', sp.simplify(S - Sclaim) == sp.zeros(3))
    # completion: w = kappa y z + A_J A_E / C0 makes the Schur complement of C0 diagonal
    C0 = c - x ** 2 / Lp
    AJ = 1 - x + (1 + x / Lp) * y
    AE = -x + x * z / Lp
    Sc = Sclaim.subs(w, kap * y * z + AJ * AE / C0)
    S2 = sp.simplify(Sc[:2, :2] - Sc[:2, 2:] * Sc[2:, 2:].inv() * Sc[2:, :2])
    check('PCBE: chosen JE entry makes the Schur complement of C0 diagonal', sp.simplify(S2[0, 1]) == 0)
    # large_pc Lemma 1: base (P, C, A, D)
    a, p, d = sp.symbols('a p d', positive=True)
    Gb = sp.Matrix([[p, -a, 0, 0], [-a, c, 0, 0], [0, 0, a, 0], [0, 0, 0, d]])
    X = {'p': [1, 0, 0, 0], 'c': [0, 1, 0, 0], 'a': [0, 0, 1, 0], 'd': [0, 0, 0, 1], 'j': [0, d / (a + c), -d / (a + c), 0]}
    Mlp = M_from(F, X, Gb)
    check('large_pc: |J|^2 = d^2/(a+c); <J,P>=<J,A>=-da/(a+c); <J,C>=dc/(a+c)',
          sp.simplify(Mlp[idx['j'], idx['j']] - d ** 2 / (a + c)) == 0 and sp.simplify(Mlp[idx['j'], idx['p']] + d * a / (a + c)) == 0
          and sp.simplify(Mlp[idx['j'], idx['a']] + d * a / (a + c)) == 0 and sp.simplify(Mlp[idx['j'], idx['c']] - d * c / (a + c)) == 0)
    gram_claim(F, 'IIIb_plus_large_pc Lemma 1', Mlp, {'p': p, 'j': d ** 2 / (a + c), 'c': c, 'd': d, 'a': a},
               {a: sp.Rational(1, 2), p: sp.Rational(2, 3), c: sp.Rational(3, 4), d: sp.Rational(5, 4)}, 7)
    check('large_pc: base Gram of (P,C,A) PSD iff pc >= a^2 (det = a(pc - a^2))', sp.factor(Gb[:3, :3].det() - a * (p * c - a ** 2)) == 0)
    # merged-C Lemma 1: six Grams
    Jr, Er, Ll = sp.symbols('Jr Er Lleaf', positive=True)   # j = Jr^2, e = Er^2, leaf coefficient Ll
    for leaf, sg in [('a', 1), ('b', -1)]:
        # crossing: base (J, leaf, D), E = -(Er/Jr) J
        xc = d * Er / (Jr + Er)
        Gb = sp.diag(Jr ** 2, Ll, d)
        X = {'j': [1, 0, 0], 'e': [-Er / Jr, 0, 0], leaf: [0, 1, 0], 'd': [0, 0, 1],
             'c': [(d - xc) / Jr ** 2, sg * xc / Ll, 0]}
        Mc = M_from(F, X, Gb)
        ccross = Mc[idx['c'], idx['c']]
        check('merged-C crossing leaf %s: <C,J>=d-x, <C,E>=-x, <C,leaf>=%sx with x = d sqrt(e)/(sqrt j+sqrt e)' % (leaf, '+' if sg > 0 else '-'),
              sp.simplify(Mc[idx['c'], idx['j']] - (d - xc)) == 0 and sp.simplify(Mc[idx['c'], idx['e']] + xc) == 0
              and sp.simplify(Mc[idx['c'], idx[leaf]] - sg * xc) == 0)
        check('merged-C crossing leaf %s: |C|^2 = d^2(L+e)/(L(sqrt j+sqrt e)^2)' % leaf,
              sp.simplify(ccross - d ** 2 * (Ll + Er ** 2) / (Ll * (Jr + Er) ** 2)) == 0)
        gram_claim(F, 'merged-C crossing leaf ' + leaf, Mc, {'j': Jr ** 2, 'e': Er ** 2, leaf: Ll, 'd': d, 'c': ccross},
                   {Jr: sp.Rational(3, 2), Er: sp.Rational(2, 3), Ll: sp.Rational(5, 4), d: sp.Rational(7, 5)}, 8)
        # stationary: base (J, E, leaf, D) with <J,E> = -L
        jj, ee = sp.symbols('jj ee', positive=True)
        Gb = sp.Matrix([[jj, -Ll, 0, 0], [-Ll, ee, 0, 0], [0, 0, Ll, 0], [0, 0, 0, d]])
        X = {'j': [1, 0, 0, 0], 'e': [0, 1, 0, 0], leaf: [0, 0, 1, 0], 'd': [0, 0, 0, 1],
             'c': [d / (jj + Ll), 0, sg * d / (jj + Ll), 0]}
        Ms = M_from(F, X, Gb)
        check('merged-C stationary leaf %s: |C|^2 = d^2/(j+L)' % leaf, sp.simplify(Ms[idx['c'], idx['c']] - d ** 2 / (jj + Ll)) == 0)
        gram_claim(F, 'merged-C stationary leaf ' + leaf, Ms, {'j': jj, 'e': ee, leaf: Ll, 'd': d, 'c': d ** 2 / (jj + Ll)},
                   {jj: sp.Rational(3, 2), ee: sp.Rational(2, 1), Ll: sp.Rational(5, 4), d: sp.Rational(7, 5)}, 9)
        # half-energy: <J,E> = -j
        Gb = sp.Matrix([[jj, -jj, 0, 0], [-jj, ee, 0, 0], [0, 0, Ll, 0], [0, 0, 0, d]])
        X = {'j': [1, 0, 0, 0], 'e': [0, 1, 0, 0], leaf: [0, 0, 1, 0], 'd': [0, 0, 0, 1],
             'c': [d / (2 * jj), 0, sg * d / (2 * Ll), 0]}
        Mh = M_from(F, X, Gb)
        ch = d ** 2 / (4 * jj) + d ** 2 / (4 * Ll)
        check('merged-C half-energy leaf %s: |C|^2 = d^2/(4j) + d^2/(4L)' % leaf, sp.simplify(Mh[idx['c'], idx['c']] - ch) == 0)
        gram_claim(F, 'merged-C half-energy leaf ' + leaf, Mh, {'j': jj, 'e': ee, leaf: Ll, 'd': d, 'c': ch},
                   {jj: sp.Rational(3, 2), ee: sp.Rational(2, 1), Ll: sp.Rational(5, 4), d: sp.Rational(7, 5)}, 10)
    return F


# =====================================================================================
# Section ALG: closed-form algebra of the statements (thresholds, specialisations, allocations)
# =====================================================================================
def Z(expr):
    return sp.simplify(sp.cancel(sp.together(expr))) == 0


def sec_alg():
    head('ALG G1: T44 Prop. 2, scalar test, SOC/Hsum')
    r, k, v, w = sp.symbols('r k v w')
    x = (v * (2 * r - k) - 1) / (r - k)
    y = k * (1 - r * v) / (w * (r - k))
    J, zz, bb, cc, jj = sp.symbols('J z b c j', positive=True)
    pz = (x ** 2 + y ** 2).subs({v: 0, r: J, k: -zz, w: sp.sqrt(bb)})
    check('Prop.2: Prop.1 at v=0, r=sqrt j, k=-z, w=sqrt b gives (1+z^2/b)/(sqrt j+z)^2', Z(pz - (1 + zz ** 2 / bb) / (J + zz) ** 2))
    dz = sp.diff((1 + zz ** 2 / bb) / (J + zz) ** 2, zz)
    check('Prop.2: derivative = 2(z sqrt j - b)/(b (sqrt j+z)^3)', Z(dz - 2 * (zz * J - bb) / (bb * (J + zz) ** 3)))
    C = sp.symbols('C', positive=True)
    check('Prop.2: at z = sqrt c value (b+c)/(b(sqrt j+sqrt c)^2)', Z(((1 + zz ** 2 / bb) / (J + zz) ** 2).subs(zz, C) - (bb + C ** 2) / (bb * (J + C) ** 2)))
    check('Prop.2: at z = b/sqrt j value 1/(j+b)', Z(((1 + zz ** 2 / bb) / (J + zz) ** 2).subs(zz, bb / J) - 1 / (J ** 2 + bb)))
    check('Prop.2: branches agree at jc = b^2', Z(((bb + C ** 2) / (bb * (J + C) ** 2)).subs(bb, J * C) - 1 / (J ** 2 + J * C)))
    check('Prop.2: j=0 case (r=0,k=-sqrt c): threshold 1/c+1/b', Z((x ** 2 + y ** 2).subs({v: 0, r: 0, k: -C, w: sp.sqrt(bb)}) - (1 / C ** 2 + 1 / bb)))
    check('Prop.2 example r=1,k=-1,v=0,w=1: p = 1/2, j=c=b=1', (x ** 2 + y ** 2).subs({r: 1, k: -1, v: 0, w: 1}) == sp.Rational(1, 2))
    t, Jq, Cq, m = sp.symbols('t Jq Cq m', positive=True)
    sub = {v: t / sp.sqrt(1 - t ** 2), r: Jq * sp.sqrt(1 - t ** 2), k: -Cq * sp.sqrt(1 - t ** 2), w: sp.sqrt(bb)}
    rhs = (((2 * Jq + Cq) * t - 1) ** 2 / (1 - t ** 2) + Cq ** 2 / bb * (1 - Jq * t) ** 2) / (Jq + Cq) ** 2
    check('scalar test (1): Prop.1 substitution gives the displayed p-threshold', Z((x ** 2 + y ** 2).subs(sub) - rhs))
    check('scalar test (1): j and c thresholds become exactly j, c', Z((r ** 2 * (1 + v ** 2)).subs(sub) - Jq ** 2) and Z((k ** 2 * (1 + v ** 2)).subs(sub) - Cq ** 2))
    check('scalar test (1): v^2 <= m iff t^2 <= m/(1+m)', Z((v ** 2).subs(sub) - t ** 2 / (1 - t ** 2)) and Z((m / (1 + m)) / (1 - m / (1 + m)) - m))
    a, e, d = sp.symbols('a e d', positive=True)
    d1, d2 = d * a / (a + e), d * e / (a + e)
    check('Hsum: old-H atoms (15,53;18),(8,53;18) at threshold with d split a:e sum to j = d + d^2/(a+e)',
          Z((d1 + d1 ** 2 / a) + (d2 + d2 ** 2 / e) - (d + d ** 2 / (a + e))))

    head('ALG G3: T46 core (2), t*, Cor. 3; six-ray Cor. 2-4 and singular-leaf section')
    q, b, be, t, j = sp.symbols('q b beta t j')
    s = q + b + 2 * be
    h = 1 + be
    kq = 2 - q - b
    l = q * b - be ** 2
    u = (q + b - 2 * be) / l
    p1 = ((j - t ** 2) * h ** 2 - 4 * t * h * (t - 1) + s * (t - 1) ** 2) / (s * j - (s + 4) * t ** 2)
    p2 = h ** 2 / s + (s + kq * t) ** 2 / (s * (s * j - (s + 4) * t ** 2))
    check('T46 core (2): threshold (1) = h^2/s + (s+kt)^2/(s(sj-(s+4)t^2))', Z(p1 - p2))
    check('T46 core: derivative of (2) = 2(s+kt)(kj+(s+4)t)/(sj-(s+4)t^2)^2',
          Z(sp.diff(p2, t) - 2 * (s + kq * t) * (kq * j + (s + 4) * t) / (s * j - (s + 4) * t ** 2) ** 2))
    check('T46 core: u - 4/s = (q-b)^2/(s l)', Z(u - 4 / s - (q - b) ** 2 / (s * l)))
    check('T46 core: t* = -kj/(s+4) gives 1/j + (s+2k)/4 + k^2/(4(s+4))',
          Z(p2.subs(t, -kq * j / (s + 4)) - (1 / j + (s + 2 * kq) / 4 + kq ** 2 / (4 * (s + 4)))))
    check('T46 core: at t=0 threshold 1/j + (1+beta)^2/(q+b+2beta); derivative at 0 = 2(2-q-b)/((q+b+2beta) j)',
          Z(p1.subs(t, 0) - (1 / j + h ** 2 / s)) and Z(sp.diff(p1, t).subs(t, 0) - 2 * kq / (s * j)))
    check('T46 Cor.3: beta = 1-S, t = 0 gives 1/j + 2 - S', Z(p1.subs({t: 0, be: 1 - q - b}) - (1 / j + 2 - q - b)))
    check('T46 example q=1/8,b=4,j=4,beta=-1/2,t=1/4: threshold 399/1543',
          p1.subs({q: sp.Rational(1, 8), b: 4, j: 4, be: sp.Rational(-1, 2), t: sp.Rational(1, 4)}) == sp.Rational(399, 1543))
    kk_ = 2 - q - b
    check('T46 core: beta=1-q-b, t*: 1/j + k(k+3)/(k+4) with s=k',
          Z((1 / j + (s + 2 * kq) / 4 + kq ** 2 / (4 * (s + 4))).subs(be, 1 - q - b) - (1 / j + kk_ * (kk_ + 3) / (kk_ + 4))))
    # Cor.3 boundary: every point of (S-1)^2 = qb with q,b >= 0 has nonzero gradient of qb-(S-1)^2
    qq, bq = sp.symbols('qq bq', real=True)
    fbd = qq * bq - (qq + bq - 1) ** 2
    sol = sp.solve([sp.diff(fbd, qq), sp.diff(fbd, bq)], [qq, bq], dict=True)
    check('T46 Cor.3 closure: only critical point of qb-(q+b-1)^2 is (2/3,2/3), value 1/3 != 0, so every boundary point is a limit of strict points',
          sol == [{qq: sp.Rational(2, 3), bq: sp.Rational(2, 3)}] and fbd.subs(sol[0]) == sp.Rational(1, 3))
    # six-ray Cor 2-4 in T47 roles (b>1, 0<a<1/b): use the singular-leaf formulas with (q,b) -> (b47, a47)
    B, A, P, E, T = sp.symbols('B A P E T', positive=True)   # b47 = B^2, a47 = A^2
    r = A * B
    kap = 1 + 1 / B ** 2
    rho = 1 + r / B ** 2
    eta = (1 - A ** 2) / (B ** 2 - r)
    kk = B ** 2 * (1 - r) ** 2 / (B ** 2 - r) ** 2
    cst = lambda tt: kk + (1 - eta * tt) ** 2 / (E - kap * tt ** 2)
    pst = lambda tt: rho ** 2 * tt ** 2 / (E - kap * tt ** 2)
    D0 = kap + rho ** 2 / P
    Aa = eta ** 2 * E
    N = D0 ** 2 + Aa * (2 * D0 - kap + Aa)
    t3 = eta * E / (D0 + Aa)
    check('six-ray Cor.3 (T47 roles): t = eta e/M0 gives consumed c = k + D0^2/(eN)', Z(cst(t3) - (kk + D0 ** 2 / (E * N))))
    check('six-ray Cor.3: e - kappa t^2 = eN/M0^2 and N - rho^2 A/p = D0^2 + D0 A + A^2',
          Z(E - kap * t3 ** 2 - E * N / (D0 + Aa) ** 2) and Z(N - rho ** 2 * Aa / P - (D0 ** 2 + D0 * Aa + Aa ** 2)))
    check('six-ray Cor.3: consumed p = rho^2 A/N (< p by the previous identity)', Z(pst(t3) - rho ** 2 * Aa / N))
    t4 = eta * E / D0
    check('six-ray Cor.4 branch A<=D0: consumed c = k + (D0-A)^2/(e(D0^2-kappa A)), unused p = p D0(D0-A)/(D0^2-kappa A)',
          Z(cst(t4) - (kk + (D0 - Aa) ** 2 / (E * (D0 ** 2 - kap * Aa)))) and Z(P - pst(t4) - P * D0 * (D0 - Aa) / (D0 ** 2 - kap * Aa)))
    t5 = 1 / eta
    check('six-ray Cor.4 branch A>=D0: t=1/eta gives consumed c = k and consumed p = rho^2/(A-kappa)',
          Z(cst(t5) - kk) and Z(pst(t5) - rho ** 2 / (Aa - kap)))
    check('six-ray Cor.4: residual difference vs Cor.3 = A^2(2D0(D0-kappa)+kappa A - A^2)/(e(D0^2-kappa A)((D0+A)^2-kappa A))',
          Z((kk + D0 ** 2 / (E * N)) - (kk + (D0 - Aa) ** 2 / (E * (D0 ** 2 - kap * Aa)))
            - Aa ** 2 * (2 * D0 * (D0 - kap) + kap * Aa - Aa ** 2) / (E * (D0 ** 2 - kap * Aa) * ((D0 + Aa) ** 2 - kap * Aa))))
    t2 = eta * E / kap
    check('six-ray Cor.2: t = eta e/kappa gives c* = 1/e + k - eta^2/kappa and p* = rho^2 eta^2 e/(kappa(kappa - eta^2 e))',
          Z(cst(t2) - (1 / E + kk - eta ** 2 / kap)) and Z(pst(t2) - rho ** 2 * eta ** 2 * E / (kap * (kap - eta ** 2 * E))))
    check('six-ray Cor.2: p* <= p  <=>  e <= p kappa^2/(eta^2 (p kappa + rho^2))  (difference identity)',
          Z((P - pst(t2)) * kap * (kap - eta ** 2 * E) - (P * kap ** 2 - eta ** 2 * E * (P * kap + rho ** 2))))

    # T47 "old stated core" c >= 1/e + k: T46 core Thm 1 at t = 0, beta = -r (singular leaf), T47 roles
    Gs = sp.Matrix([[B ** 2, -r, 0, 0], [-r, A ** 2, 0, 0], [0, 0, E, 0], [0, 0, 0, 1]])
    pcore = (1 / E + (1 - r) ** 2 / (B ** 2 + A ** 2 - 2 * r))
    check('T47 old core c >= 1/e + k: core threshold at t=0, beta=-sqrt(ab) equals 1/e + b(1-r)^2/(b-r)^2',
          Z(p1.subs({t: 0, be: -r, q: B ** 2, b: A ** 2, j: E}) - (1 / E + kk)) and Z(pcore - (1 / E + kk)))
    check('T47 old core: Gram at beta=-r, t=0 is PSD with rank-one leaf block (det 0), and Delta = s e > 0 for a != b',
          sp.simplify(Gs[:2, :2].det()) == 0 and Z((s * j - (s + 4) * t ** 2).subs({t: 0, be: -r, q: B ** 2, b: A ** 2, j: E}) - (B - A) ** 2 * E))
    # PCBE_low_e chart (T47 PCBE_low_e_q_monotone p.1): failure of (2) with L = b, e <= b(b+1)
    bb_, ee_, sv_, cv_ = sp.symbols('bb ee sv cv', positive=True)
    Bl = 1 / ee_ + 1 / (bb_ + 1)
    check('PCBE_low_e chart: B = 1/e + 1/(b+1) >= 1/b  iff  e <= b(b+1)', Z((Bl - 1 / bb_) - (bb_ * (bb_ + 1) - ee_) / (ee_ * bb_ * (bb_ + 1))))
    cfail = sv_ ** 2 * Bl
    Kc = (cfail + 1) / (1 - sv_) ** 2      # sqrt((c+1)/K) = 1 - s  defines the boundary K
    Jst = 1 / (Kc - 1 / (bb_ + 1))
    check('PCBE_low_e chart: with c = s^2 B, failure of (2) is j < J* = (b+1)(1-s)^2/(b + 2s + (b+1)s^2/e)',
          Z(Jst - (bb_ + 1) * (1 - sv_) ** 2 / (bb_ + 2 * sv_ + (bb_ + 1) * sv_ ** 2 / ee_)))

    head('ALG G4: PCBE scalar identity (3), simplification, x*; merged-C allocation')
    Aq, Bq, Jv, kp, tt = sp.symbols('A B J kappa t')
    check('PCBE (3): (A+Bt)^2/(J-kt^2) - (A^2/J - B^2/k) = (kAt+BJ)^2/(kJ(J-kt^2))',
          Z((Aq + Bq * tt) ** 2 / (Jv - kp * tt ** 2) - (Aq ** 2 / Jv - Bq ** 2 / kp) - (kp * Aq * tt + Bq * Jv) ** 2 / (kp * Jv * (Jv - kp * tt ** 2))))
    x, Lp, j, e = sp.symbols('x L j e', positive=True)
    kap = 1 + 1 / Lp
    K = 1 / j + 1 / (Lp + 1)
    check('PCBE: x^2/L + (1-x)^2/j - (1+x/L)^2/kappa = K(1-x)^2 - 1', Z(x ** 2 / Lp + (1 - x) ** 2 / j - (1 + x / Lp) ** 2 / kap - (K * (1 - x) ** 2 - 1)))
    check('PCBE: x^2/L + x^2(1/e - 1/(L(L+1))) = x^2(1/e + 1/(L+1))', Z(x ** 2 / Lp + x ** 2 * (1 / e - 1 / (Lp * (Lp + 1))) - x ** 2 * (1 / e + 1 / (Lp + 1))))
    Kk, Bk = sp.symbols('K Bk', positive=True)
    xs = (Kk - 1) / (Kk + sp.sqrt(Kk + Bk * (Kk - 1)))
    check('PCBE: x* solves K(1-x)^2 - 1 = B x^2', Z(sp.expand(Kk * (1 - xs) ** 2 - 1 - Bk * xs ** 2)))
    d, a_, b_, jv, ev_ = sp.symbols('d a b jv ev', positive=True)
    A0 = a_ + b_
    da, db = a_ / A0, b_ / A0
    cross = lambda dd, jj, ee, LL: dd ** 2 * (LL + ee) / (LL * (sp.sqrt(jj) + sp.sqrt(ee)) ** 2)
    stat = lambda dd, jj, LL: dd ** 2 / (jj + LL)
    half = lambda dd, jj, LL: dd ** 2 / (4 * jj) + dd ** 2 / (4 * LL)
    check('merged-C: proportional merge of crossing bounds gives the A0 formula',
          Z(cross(da * d, da * jv, da * ev_, a_) + cross(db * d, db * jv, db * ev_, b_) - cross(d, jv, ev_, A0)))
    check('merged-C: proportional merge of stationary bounds gives d^2/(j+A0); j_i e_i - L_i^2 = delta_i^2 (je - A0^2)',
          Z(stat(da * d, da * jv, a_) + stat(db * d, db * jv, b_) - stat(d, jv, A0)) and Z((da * jv) * (da * ev_) - a_ ** 2 - da ** 2 * (jv * ev_ - A0 ** 2)))
    check('merged-C: proportional merge of half-energy bounds gives d^2/(4j)+d^2/(4A0)',
          Z(half(da * d, da * jv, a_) + half(db * d, db * jv, b_) - half(d, jv, A0)))
    xL = sp.symbols('xL', positive=True)
    LL = sp.symbols('LL', positive=True)
    Jr, Er = sp.symbols('Jr Er', positive=True)
    CL = lambda xx: xx ** 2 / LL + sp.Max(xx ** 2 / Er ** 2, (d - xx) ** 2 / Jr ** 2)
    xc = d * Er / (Jr + Er)
    check('merged-C: at crossing x the two max-branches agree', Z((xc ** 2 / Er ** 2) - (d - xc) ** 2 / Jr ** 2))
    check('merged-C: stationary choice x = dL/(j+L) gives x^2/L + (d-x)^2/j = d^2/(j+L)',
          Z(((d * LL / (Jr ** 2 + LL)) ** 2 / LL + (d - d * LL / (Jr ** 2 + LL)) ** 2 / Jr ** 2) - d ** 2 / (Jr ** 2 + LL)))

    head('ALG G2: T45 (6) split')
    rr, ss, dd = sp.symbols('r s d', positive=True)
    d1, d2 = dd * rr / (rr + ss), dd * ss / (rr + ss)
    check('T45 (6): H(1,59;70)+H(72,59;70) at threshold with d split r:s sum to c = d + d^2/(r+s)',
          Z(d1 + d1 ** 2 / rr + d2 + d2 ** 2 / ss - (dd + dd ** 2 / (rr + ss))))

    head('ALG G5: H/chain allocation lemma (6), T48 allocations, eighth chamber')
    b1, b2, de1, de2, bb, cc, pp, qq = sp.symbols('beta1 beta2 delta1 delta2 b c p q', positive=True)
    eta1, eta2 = b1 - de1, b2 - de2
    bH, bC = bb * eta1 / b1, bb * de1 / b1
    cH, cC = cc * eta2 / b2, cc * de2 / b2
    # atom requirements on the central coefficient a (each atom exactly at threshold)
    H1 = eta1 + eta1 ** 2 / bH          # old-H aux b_H, nonminor eta1
    H2 = eta2 + eta2 ** 2 / cH
    C1 = de1 ** 2 / pp + de1 ** 2 / bC  # chain leaves p, b_C, nonminor delta1
    C2 = de2 ** 2 / cC + de2 ** 2 / qq
    target = 1 - de1 - de2 + b1 ** 2 / bb + b2 ** 2 / cc + de1 ** 2 / pp + de2 ** 2 / qq
    check('H/chain (6): sum of the four atoms at threshold uses b, c exactly and needs a = RHS of (6) (beta1+beta2=1)',
          Z((H1 + H2 + C1 + C2 - target).subs(b2, 1 - b1)) and Z(bH + bC - bb) and Z(cH + cC - cc)
          and Z((eta1 + eta2 + de1 + de2 - 1).subs(b2, 1 - b1)))
    # T48: e-side current s, j-side 1-s; e-side split h_e (old-H aux e) + alpha (chain a,e); j-side h_j + chains (q,j),(b,j)
    sv, al, ga, aa, ee, jj, Qs, qv, bv = sp.symbols('s alpha gamma a e j Q qv bv', positive=True)
    he, hj = sv - al, (1 - sv) - ga
    # e split proportional to currents: e_H = e*he/s, e_C = e*alpha/s; chains (q,j),(b,j) currents gamma*q/Q, gamma*b/Q
    eside = he + he ** 2 / (ee * he / sv) + al ** 2 / aa + al ** 2 / (ee * al / sv)
    jside = hj + hj ** 2 / (jj * hj / (1 - sv)) + sum((ga * x_ / Qs) ** 2 / x_ + (ga * x_ / Qs) ** 2 / (jj * (ga * x_ / Qs) / (1 - sv)) for x_ in (qv, bv))
    check('T48: e-side/j-side atom sums = s - alpha + alpha^2/a + s^2/e, and (1-s) - gamma + gamma^2/Q + (1-s)^2/j (Q=q+b)',
          Z(eside - (sv - al + al ** 2 / aa + sv ** 2 / ee)) and Z((jside - ((1 - sv) - ga + ga ** 2 / Qs + (1 - sv) ** 2 / jj)).subs(Qs, qv + bv)))
    phi_lo = lambda tt, AA: tt - tt ** 2 / AA
    check('T48 Cor.2: optimal alpha = min(s, a/2): value s - phi_a(s) with phi_a(t) = t - t^2/a (t<=a/2), a/4 (t>=a/2)',
          Z((sv - sv + sv ** 2 / aa) - (sv - phi_lo(sv, aa))) and Z((sv - aa / 2 + (aa / 2) ** 2 / aa) - (sv - aa / 4)))
    # allocation I: u = lam Ea, v = lam JQ
    Ea, JQ = aa * ee / (aa + ee), jj * Qs / (jj + Qs)
    lam = 1 / (Ea + JQ)
    uu, vv = lam * Ea, lam * JQ
    check('T48 allocation I: u+v = 1 and u^2(1/a+1/e) + v^2(1/j+1/Q) = lambda', Z(uu + vv - 1) and Z(uu ** 2 * (1 / aa + 1 / ee) + vv ** 2 * (1 / jj + 1 / Qs) - lam))
    u2 = (1 + 2 / jj) / (2 * (1 / aa + 1 / ee + 1 / jj))
    v2 = 1 - u2
    pst = u2 ** 2 * (1 / aa + 1 / ee) + v2 - Qs / 4 + v2 ** 2 / jj
    check('T48 allocation II: p* = 1 - u/2 + v/j - Q/4', Z(pst - (1 - u2 / 2 + v2 / jj - Qs / 4)))
    u3 = (1 + 2 / ee) / (2 * (1 / Qs + 1 / jj + 1 / ee))
    v3 = 1 - u3
    pst3 = u3 ** 2 * (1 / Qs + 1 / jj) + v3 - aa / 4 + v3 ** 2 / ee
    check('T48 allocation III: p* = 1 - u/2 + v/e - a/4', Z(pst3 - (1 - u3 / 2 + v3 / ee - aa / 4)))
    S = ee + jj
    hE, hJ, alp, bet = ee / S - aa / 2, jj / S - Qs / 2, aa / 2, Qs / 2
    p4 = hE + hE / S + alp ** 2 / aa + alp / S + hJ + hJ / S + bet ** 2 / Qs + bet / S
    check('T48 allocation IV: currents sum to 1 and p* = 1 + 1/S - (a+Q)/4', Z(hE + hJ + alp + bet - 1) and Z(p4 - (1 + 1 / S - (aa + Qs) / 4)))
    # eighth chamber (T27 (10)-(11) / T38): four atom at threshold and coherent scaling
    de, A_, c_, b_, e_ = sp.symbols('delta a c b e', positive=True)
    p10 = de * (A_ + 1) * (c_ + de) / (A_ * c_ - de)
    check('eighth chamber: four atom (a2=delta a, c, p, d2=delta) at p(delta) is exactly on ab(c-d) = d^2(a+b+c+d)',
          Z((de * A_) * c_ * (p10 - de) - de ** 2 * (de * A_ + c_ + p10 + de)))
    r0 = 1 - de
    Er0 = r0 * (A_ + 1) * (b_ + r0) / (A_ * b_ - r0)
    Eab = lambda aa_, bb_: (aa_ + 1) * (bb_ + 1) / (aa_ * bb_ - 1)
    # T24 coherent: j >= 1/e - 1/E_ab at d=1; homogeneous rescaling to (a1, b, e, j, d1) = (r0 a, b, e, j, r0)
    check('eighth chamber: rescaled coherent bound r0*(r0/e - 1/E(a, b/r0)) = r0^2/e - r0^2/E_r0',
          Z(r0 * (r0 / e_ - 1 / Eab(A_, b_ / r0)) - (r0 ** 2 / e_ - r0 ** 2 / Er0)))
    check('eighth chamber: e >= E_r0 is the rescaled four-coefficient condition', Z(r0 * Eab(A_, b_ / r0) - Er0))
    check('eighth chamber: p(delta) increasing on [0,1) for ac >= 1 (derivative numerator positive)',
          sp.simplify(sp.numer(sp.together(sp.diff(p10, de))) - (A_ + 1) * (A_ * c_ ** 2 + 2 * A_ * c_ * de - de ** 2)) == 0)


# =====================================================================================
# Section NUM: attacks over SO(4)
# =====================================================================================
PERM_ARR = np.array(L.PERMS)
SGN_ARR = np.array(L.SGN, dtype=float)
FV = np.array(L.F_VALUES, dtype=float)


def O_batch(Zp):
    p = Zp[:, :4] / np.linalg.norm(Zp[:, :4], axis=1, keepdims=True)
    q = Zp[:, 4:] / np.linalg.norm(Zp[:, 4:], axis=1, keepdims=True)

    def qm(a, b):
        a0, a1, a2, a3 = a.T
        b0, b1, b2, b3 = b.T
        return np.stack([a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3, a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
                         a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1, a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0], axis=1)
    cols = []
    for jx in range(4):
        e = np.zeros_like(p)
        e[:, jx] = 1
        cols.append(qm(qm(p, e), q))
    return np.stack(cols, axis=2)          # O[n, i, j]


def phi_batch(cvals, Zp):
    O = O_batch(Zp)
    W = np.ones((Zp.shape[0], 24))
    for i in range(4):
        W *= O[:, i, PERM_ARR[:, i]]
    return W @ (SGN_ARR * cvals)


def upd(d, extra):
    out = dict(d)
    out.update(extra)
    return out


def coef_vals(cdict):
    v = np.zeros(24)
    for f, c in cdict.items():
        v += c * FV[f]
    return v


def minimize_SO4(cvals, rng, nstart=300, iters=160):
    Zp = rng.normal(size=(nstart, 8))
    lr = 0.05
    m = np.zeros_like(Zp)
    s2 = np.zeros_like(Zp)
    hh = 1e-6
    for it in range(iters):
        g = np.zeros_like(Zp)
        for kx in range(8):
            dZ = np.zeros_like(Zp)
            dZ[:, kx] = hh
            g[:, kx] = (phi_batch(cvals, Zp + dZ) - phi_batch(cvals, Zp - dZ)) / (2 * hh)
        m = 0.9 * m + 0.1 * g
        s2 = 0.999 * s2 + 0.001 * g * g
        Zp = Zp - lr * (m / (1 - 0.9 ** (it + 1))) / (np.sqrt(s2 / (1 - 0.999 ** (it + 1))) + 1e-12)
        Zp[:, :4] /= np.linalg.norm(Zp[:, :4], axis=1, keepdims=True)
        Zp[:, 4:] /= np.linalg.norm(Zp[:, 4:], axis=1, keepdims=True)
        if it == iters // 2:
            lr = 0.01
    vals = phi_batch(cvals, Zp)
    i = int(np.argmin(vals))
    # polish best 5 with scipy
    from scipy.optimize import minimize
    best = (vals[i], Zp[i])
    for i2 in np.argsort(vals)[:5]:
        res = minimize(lambda z: phi_batch(cvals, z[None, :])[0], Zp[i2], method='BFGS', options={'gtol': 1e-13, 'maxiter': 400})
        if res.fun < best[0]:
            best = (res.fun, res.x)
    return best


def mp_phi(cdict, z):
    import mpmath as mp
    mp.mp.dps = 40
    zz = [mp.mpf(float(x)) for x in z]
    np_ = mp.sqrt(sum(x * x for x in zz[:4]))
    nq_ = mp.sqrt(sum(x * x for x in zz[4:]))
    p = [x / np_ for x in zz[:4]]
    q = [x / nq_ for x in zz[4:]]
    O = O_pxq(p, q, mp.mpf(0), mp.mpf(1))
    cv = [0] * 24
    for f, c in cdict.items():
        for s_ in range(24):
            cv[s_] += mp.mpf(c) * L.F_VALUES[f][s_]
    return leibniz(cv, O)


def attack(name, sampler, nsamp, seed, results, control=None):
    rng = np.random.default_rng(seed)
    worst = (np.inf, None, None)
    for i in range(nsamp):
        cd = sampler(rng)
        cv = coef_vals(cd)
        scale = max(abs(c) for c in cd.values())
        val, z = minimize_SO4(cv, rng)
        val /= scale
        if val < worst[0]:
            worst = (val, z, cd)
    mpv = float(mp_phi(worst[2], worst[1])) / max(abs(c) for c in worst[2].values())
    ok = check('NUM %s: %d boundary samples x 300 starts, min Phi/max|coef| = %.3e (mpmath %.3e)' % (name, nsamp, worst[0], mpv),
               worst[0] > -1e-9 and mpv > -1e-9)
    results[name] = (worst[0], mpv)
    if control is not None:
        cname, csampler = control
        rng = np.random.default_rng(seed + 7)
        cw = np.inf
        for i in range(3):
            cd = csampler(rng)
            val, z = minimize_SO4(coef_vals(cd), rng)
            cw = min(cw, val / max(abs(c) for c in cd.values()))
        check('NUM control %s (%s) is detected: min = %.3e' % (name, cname, cw), cw < -1e-8)
        results[name + ' control'] = cw


def sec_num():
    head('NUM: attacks on the stated regions, coefficients exactly on the thresholds')
    R = {}
    U = lambda rng, lo=0.1, hi=3.0: float(np.exp(rng.uniform(np.log(lo), np.log(hi))))

    # ---- G1 T44: q f2, p f1, j f53, c f12, e f8, d f18, b f7, a f15
    def T44(qv, p, j, c, e, b, a):
        return {2: qv, 1: p, 53: j, 12: c, 8: e, 18: 1.0, 7: b, 15: a}

    def s_P1(rng, frac=1.0):
        while True:
            r_, k_, v_ = rng.normal(size=3)
            w_ = U(rng)
            if abs(r_ - k_) > 0.2:
                break
        x_ = (v_ * (2 * r_ - k_) - 1) / (r_ - k_)
        y_ = k_ * (1 - r_ * v_) / (w_ * (r_ - k_))
        return T44(0, frac * (x_ ** 2 + y_ ** 2), r_ ** 2 * (1 + v_ ** 2), k_ ** 2 * (1 + v_ ** 2), v_ ** 2, w_ ** 2, v_ ** 2)
    attack('G1 T44 Prop.1', s_P1, 8, 101, R, ('p and j set to 0', lambda rng: upd(s_P1(rng), {1: 0.0, 53: 0.0})))

    def p2thr(b, j, c):
        if b == 0:
            return 1 / j
        if j * c <= b * b:
            return (b + c) / (b * (math.sqrt(j) + math.sqrt(c)) ** 2)
        return 1 / (j + b)

    def s_P2(rng):
        b, j, c = U(rng), U(rng), U(rng)
        mode = rng.integers(0, 4)
        if mode == 1:
            j = 0.0
        if mode == 2:
            b = 0.0
        return T44(0, p2thr(b, j, c), j, c, 0, b, 0)
    attack('G1 T44 Prop.2 (incl. j=0, b=0)', s_P2, 8, 102, R, ('p set to 0', lambda rng: upd(s_P2(rng), {1: 0.0})))

    def s_J(rng):
        b, p, e = U(rng), U(rng), U(rng)
        Js = (b + e) / (b * (math.sqrt(p) + math.sqrt(e)) ** 2) if p * e <= b * b else 1 / (p + b)
        return {53: Js, 1: p, 8: e, 7: b, 18: 1.0}
    attack('G1 T44 J* (transported Prop.2)', s_J, 6, 103, R)

    def s_SOC(rng):
        p = U(rng)
        return {1: p, 53: 1 / p, 18: 1.0}
    attack('G1 T44 SOC pj=1', s_SOC, 4, 104, R, ('pj=0.9', lambda rng: (lambda p: {1: p, 53: 0.9 / p, 18: 1.0})(U(rng))))

    def s_Hs(rng):
        a, e = U(rng), U(rng)
        return {53: 1 + 1 / (a + e), 15: a, 8: e, 18: 1.0}
    attack('G1 T44 Hsum j=1+1/(a+e)', s_Hs, 4, 105, R)

    # ---- G3 T46: q f27, p f1, j f59, c f5, e f2, d f28, b f6, a f4
    def s_core(rng):
        q_, b_ = U(rng), U(rng)
        be_ = rng.uniform(-0.98, 0.98) * math.sqrt(q_ * b_)
        l_ = q_ * b_ - be_ ** 2
        j_ = U(rng)
        tmax = math.sqrt(j_ / (1 + (q_ + b_ - 2 * be_) / l_))
        t_ = rng.uniform(-0.95, 0.95) * tmax
        s_ = q_ + b_ + 2 * be_
        h_ = 1 + be_
        Dl = s_ * j_ - (s_ + 4) * t_ ** 2
        p_ = ((j_ - t_ ** 2) * h_ ** 2 - 4 * t_ * h_ * (t_ - 1) + s_ * (t_ - 1) ** 2) / Dl
        return {27: q_, 1: p_, 59: j_, 6: b_, 28: 1.0}
    attack('G3 T46 core Thm 1', s_core, 8, 106, R, ('p set to 0', lambda rng: upd(s_core(rng), {1: 0.0})))

    def s_cor3(rng):
        q_ = rng.uniform(0.0, 4 / 3)
        # boundary (q+b-1)^2 = q b : b^2 + b(q-2) ... solve b^2 + (2q-2-q) b + (q-1)^2 = 0
        disc = (q_ - 2) ** 2 - 4 * (q_ - 1) ** 2
        b_ = ((2 - q_) + (1 if rng.random() < 0.5 else -1) * math.sqrt(max(disc, 0))) / 2
        j_ = U(rng)
        return {27: q_, 1: 1 / j_ + 2 - q_ - b_, 59: j_, 6: max(b_, 0.0), 28: 1.0}
    attack('G3 T46 Cor.3 on the boundary (S-1)^2 = qb', s_cor3, 6, 107, R)

    def s_six(rng):
        while True:
            q_, b_ = U(rng), U(rng)
            be_ = rng.uniform(-0.95, 0.95) * math.sqrt(q_ * b_)
            l_ = q_ * b_ - be_ ** 2
            t_, v_ = rng.normal(size=2) * 0.7
            j_ = t_ ** 2 + (b_ * t_ ** 2 - 2 * be_ * t_ * v_ + q_ * v_ ** 2) / l_ + U(rng, 0.05, 2)
            break
        s_ = q_ + b_ + 2 * be_
        h_ = 1 + be_
        D6 = s_ * (j_ - t_ ** 2) - (t_ + v_) ** 2
        ps = ((j_ - t_ ** 2) * h_ ** 2 - 2 * h_ * (t_ + v_) * (1 + v_) + s_ * (1 + v_) ** 2) / D6
        cs = (t_ - v_) ** 2 / (j_ - t_ ** 2 - v_ ** 2 / b_)
        return {27: q_, 1: ps, 59: j_, 5: cs, 28: 1.0, 6: b_}
    attack('G3 T46 six-ray Thm 1', s_six, 8, 108, R, ('p set to 0', lambda rng: upd(s_six(rng), {1: 0.0})))

    def s_sing(rng):
        while True:
            q_, b_ = U(rng), U(rng)
            if abs(q_ - b_) > 0.1:
                break
        c_, j_ = U(rng), U(rng)
        r_ = math.sqrt(q_ * b_)
        kap, rho, eta = 1 + 1 / q_, 1 + r_ / q_, (1 - b_) / (q_ - r_)
        kk = q_ * (1 - r_) ** 2 / (q_ - r_) ** 2
        D0 = kap + rho ** 2 / c_
        N = D0 ** 2 + eta ** 2 * j_ * (2 * D0 - kap + eta ** 2 * j_)
        return {27: q_, 1: kk + D0 ** 2 / (j_ * N), 59: j_, 5: c_, 28: 1.0, 6: b_}
    attack('G3 T46 singular-leaf section p = P*', s_sing, 8, 109, R)

    # ---- G4 T47: q f2, p f9, j f1, c f45, e f5, d f35, b f0, a f7
    def s_cor4(rng):
        b_ = 1 + U(rng, 0.05, 3)
        a_ = rng.uniform(0.02, 0.98) / b_
        p_, e_ = U(rng), U(rng)
        r_ = math.sqrt(a_ * b_)
        kap, rho, eta = 1 + 1 / b_, 1 + r_ / b_, (1 - a_) / (b_ - r_)
        kk = b_ * (1 - r_) ** 2 / (b_ - r_) ** 2
        A_ = eta ** 2 * e_
        D0 = kap + rho ** 2 / p_
        c_ = kk + ((D0 - A_) ** 2 / (e_ * (D0 ** 2 - kap * A_)) if A_ <= D0 else 0.0)
        return {0: b_, 7: a_, 9: p_, 5: e_, 45: c_, 35: 1.0}
    attack('G3 six-ray Cor.4 (T47 roles)', s_cor4, 6, 110, R)

    def pcbe_c(p_, b_, j_, e_):
        L_ = min(p_, b_)
        K = 1 / j_ + 1 / (L_ + 1)
        B = max(1 / e_ + 1 / (L_ + 1), 1 / L_)
        if K <= 1:
            return 0.0
        xs = (K - 1) / (K + math.sqrt(K + B * (K - 1)))
        return B * xs * xs

    def s_pcbe(rng):
        p_, b_, j_, e_ = U(rng, 0.05, 2), U(rng, 0.05, 2), U(rng, 0.1, 2), U(rng, 0.1, 3)
        return {9: p_, 0: b_, 1: j_, 5: e_, 45: pcbe_c(p_, b_, j_, e_), 35: 1.0}
    attack('G4 T47 PCBE Thm 1', s_pcbe, 8, 111, R, ('p and b set to 0', lambda rng: upd(s_pcbe(rng), {9: 0.0, 0: 0.0})))

    def s_lpc(rng):
        a_, p_ = U(rng), U(rng)
        c_ = a_ * a_ / p_
        return {7: a_, 9: p_, 45: c_, 1: 1 / (a_ + c_), 35: 1.0}
    attack('G4 T47 large_pc Lemma 1 (pc=a^2, j=1/(a+c))', s_lpc, 6, 112, R, ('j set to 0', lambda rng: upd(s_lpc(rng), {1: 0.0})))

    def s_mc(rng):
        a_, b_, j_, e_ = U(rng), U(rng), U(rng), U(rng)
        A0 = a_ + b_
        mode = rng.integers(0, 3)
        if mode == 2:
            e_ = j_ + U(rng, 0.01, 2)
            c_ = 1 / (4 * j_) + 1 / (4 * A0)
        elif j_ * e_ <= A0 * A0:
            c_ = (A0 + e_) / (A0 * (math.sqrt(j_) + math.sqrt(e_)) ** 2)
        else:
            c_ = 1 / (j_ + A0)
        return {7: a_, 0: b_, 1: j_, 5: e_, 45: c_, 35: 1.0}
    attack('G4 T47 merged-C Lemma 1 (crossing/stationary/half-energy)', s_mc, 9, 113, R)

    # ---- G2 T45: p f0, r f1, a f2, b f12, c f59, q f67, d f70, s f72
    def s_45(rng):
        r_, s_ = U(rng), U(rng)
        return {1: r_, 72: s_, 59: 1 + 1 / (r_ + s_), 70: 1.0}
    attack('G2 T45 (6) c = 1 + 1/(r+s)', s_45, 4, 114, R)

    def s_45ch(rng):
        p_, b_, r_, q_, s_ = U(rng), U(rng), U(rng), U(rng), U(rng)
        k_ = rng.integers(0, 4)
        if k_ == 0:
            a_ = U(rng)
            return {2: a_, 59: 1 / a_, 70: 1.0}
        if k_ == 1:
            return {0: p_, 12: b_, 2: 1 / p_ + 1 / b_, 70: 1.0}
        if k_ == 2:
            return {0: p_, 1: r_, 59: 1 / p_ + 1 / r_, 70: 1.0}
        return {67: q_, 72: s_, 59: 1 / q_ + 1 / s_, 70: 1.0}
    attack('G2 T45 (2)-(5)', s_45ch, 8, 115, R)

    # ---- G5 H/chain allocation (6) with a at the minimised energy J
    from scipy.optimize import minimize_scalar

    def Jmin(p_, b_, c_, q_):
        def en(b1):
            b2 = 1 - b1
            d1, d2 = min(p_ / 2, b1), min(q_ / 2, b2)
            return 1 - d1 - d2 + b1 ** 2 / b_ + b2 ** 2 / c_ + d1 ** 2 / p_ + d2 ** 2 / q_
        res = minimize_scalar(en, bounds=(0, 1), method='bounded', options={'xatol': 1e-12})
        return res.fun
    HC = {'T28': (33, 5, 4, 6, 2, 46), 'T38': (12, 39, 0, 67, 2, 70), 'T42': (3, 23, 0, 6, 5, 54), 'T36': (2, 29, 0, 7, 5, 40)}
    for nm, (fb, fc, fp, fq, fa, fk) in HC.items():
        def s_hc(rng, fb=fb, fc=fc, fp=fp, fq=fq, fa=fa, fk=fk):
            p_, b_, c_, q_ = U(rng), U(rng), U(rng), U(rng)
            return {fb: b_, fc: c_, fp: p_, fq: q_, fa: Jmin(p_, b_, c_, q_), fk: 1.0}
        ctrl = None
        if nm == 'T42':
            ctrl = ('a set to 0', lambda rng, s=s_hc, fa=fa: upd(s(rng), {fa: 0.0}))
        attack('G5 H/chain (6) %s roles' % nm, s_hc, 4, 116 + len(R), R, ctrl)

    # ---- G5 T48: q f7, p f5, j f29, c f39, e f2, d f40, b f41, a f0 ; p = min_s T(s)
    def s_48(rng):
        a_, e_, j_, q_, b_ = U(rng), U(rng), U(rng), U(rng), U(rng)
        Q_ = q_ + b_
        ph = lambda t, A: t - t * t / A if t <= A / 2 else A / 4
        res = minimize_scalar(lambda s: 1 + s * s / e_ + (1 - s) ** 2 / j_ - ph(s, a_) - ph(1 - s, Q_), bounds=(0, 1),
                              method='bounded', options={'xatol': 1e-12})
        return {7: q_, 5: res.fun, 29: j_, 2: e_, 40: 1.0, 41: b_, 0: a_}
    attack('G5 T48 allocations I-IV (p = min_s T(s))', s_48, 6, 140, R)

    # ---- sharp control: old-H lowered by 10% (T44 roles (15,53;18)): j = 0.9(1 + 1/a)
    def s_Hlow(rng):
        a_ = U(rng)
        return {15: a_, 53: 0.9 * (1 + 1 / a_), 18: 1.0}
    rng = np.random.default_rng(150)
    cw = np.inf
    for i in range(4):
        cd = s_Hlow(rng)
        v_, z_ = minimize_SO4(coef_vals(cd), rng)
        cw = min(cw, v_ / max(cd.values()))
    check('NUM sharp control: old-H (15,53;18) with main coefficient lowered 10%% is detected: min = %.3e' % cw, cw < -1e-7)
    return R


# =====================================================================================
if __name__ == '__main__':
    print('atoms2b-verify: own re-verification of so4-atoms-part2b; sections', SECTIONS)
    if 'ROLES' in SECTIONS:
        sec_roles()
    if 'SRC' in SECTIONS:
        sec_src()
    if 'GRAM44' in SECTIONS:
        sec_gram44()
    if 'GRAM46' in SECTIONS:
        sec_gram46()
    if 'GRAM47' in SECTIONS:
        sec_gram47()
    if 'ALG' in SECTIONS:
        sec_alg()
    if 'NUM' in SECTIONS:
        sec_num()
    print('\nTOTAL: %d PASS, %d FAIL, %.1f s' % (NPASS, NFAIL, time.time() - T0))
    sys.exit(1 if NFAIL else 0)
