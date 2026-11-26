"""Part 2a of the SO(4) atom audit: rank-two family (R), T24 coherent / parallel / double-SOC
Grams (F1), T25 partial-four Gram (F2), T41 bad17 transported bounds (F3), inventory56 capacities (F4).

Own code. Uses ../lib/moc4lib.py. Nothing from the supplement is executed or imported. The only
supplement input is JSON DATA (cubic features, sign maps, congruence T) of three faces, read with
json.load from $MOC4_SUPP/proofs/manin_transfer_20261021/ (or argv[1] = supplement root):
  inventory56 face {0,1,5,7,9,12,17}: final_referee_inventory56_boundary17_{cubic_reduction,budget_PSD7}_20261023.json
  T24 face {3,4,7,8,10,12}:          root_moc4_T24_{cubic_reduction,budget_PSD8}_20261022.json
  T25 face {2,3,5,7,8,12}:           final_referee_T25_{cubic_reduction,budget_PSD8}_20261022.json
From that data I use only: the cubic polynomials ('features'), the four signed coordinate maps
('left_four') and the matrix T. Everything else (Gram columns, relation matrices R, budget
diagonals) is re-derived here from the definitions, and the supplement's R is used only for a
cross-comparison.

Model (plain convention, as in so4_atoms.py A5): O(p,q): x -> p x q, p, q quaternions, entries
bilinear; on |p| = |q| = 1 this is onto SO(4). N = |p|^2 - |q|^2. For a face with sources s_0..s_{m-1}
and cubics g = T^{-1} F (per parity), M_ik = sum_{parity, j} g_i(S_j y) g_k(S_j y). The Gram map is
   Mlin(G) = |p|^2 * sum_{i,k} G_ik M_ik      (mod N),
and the claim to verify is Mlin(G) = 2 Phi_{cost(G)} (mod N), with cost(G) = Binv diag(G).
If G is PSD then the right side is |p|^2 * (sum of squares), so Phi_cost >= 0 on SO(4).

Sections: A exact model; B face data and own relation spaces; C constructions (exact identities,
exact LDL^T PSD certificates); D orbit / transport checks; E numeric attack.
Run: MOC4_SUPP=<supplement root> python -I -B so4_atoms_part2a.py
"""
import os, sys, json, time, math
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lib'))
import moc4lib as L
import flint
import sympy as sp
import numpy as np
from scipy.optimize import minimize
import mpmath

T0 = time.time()
SUPP = os.environ.get('MOC4_SUPP') or (sys.argv[1] if len(sys.argv) > 1 else None)
BASE = os.path.join(SUPP, 'proofs', 'manin_transfer_20261021') if SUPP else None
RES = {'exact': [], 'constructions': {}, 'numeric': {}, 'orbits': {}}
NFAIL = 0


def check(name, ok, info=''):
    global NFAIL
    RES['exact'].append({'check': name, 'pass': bool(ok), 'info': info})
    if not ok:
        NFAIL += 1
    print(('PASS ' if ok else 'FAIL ') + name + ((' :: ' + info) if info else ''), flush=True)


# =============================================================================================
# A. quaternion model, plain convention (exact)
# =============================================================================================
NAMES = ('p0', 'p1', 'p2', 'p3', 'q0', 'q1', 'q2', 'q3')
QC = flint.fmpq_mpoly_ctx.get(NAMES, 'deglex')
Y = QC.gens()
ZERO = QC.from_dict({})
ONE = QC.from_dict({(0,) * 8: 1})
pv, qv = Y[:4], Y[4:]
NP = sum(x * x for x in pv)
NQ = sum(x * x for x in qv)


def qmul(a, b):
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    return (a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3, a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
            a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1, a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0)


cols = []
for c in range(4):
    e = [ZERO] * 4
    e[c] = ONE
    cols.append(qmul(qmul(pv, e), qv))
OM = [[cols[c][r] for c in range(4)] for r in range(4)]
W = []
for s, sg in zip(L.PERMS, L.SGN):
    t = OM[0][s[0]] * OM[1][s[1]] * OM[2][s[2]] * OM[3][s[3]]
    W.append(t if sg > 0 else -t)
PHI = [sum((v * w for v, w in zip(L.F_VALUES[j], W) if v), ZERO) for j in range(82)]
S2 = NP * NQ
okO = all(sum((OM[k][i] * OM[k][j] for k in range(4)), ZERO) == (S2 if i == j else ZERO)
          for i in range(4) for j in range(4)) and sum(W, ZERO) == S2 * S2
check('A0 plain model: O^T O = |p|^2|q|^2 I, det O = (|p|^2|q|^2)^2', okO)
SR = NQ - pv[1] ** 2 - pv[2] ** 2 - pv[3] ** 2


def redN(f):
    """Normal form modulo N = |p|^2 - |q|^2 (N is monic of degree 2 in p0): p0^2 -> SR."""
    out = ZERO
    for m, c in f.to_dict().items():
        k = m[0]
        out += QC.from_dict({(0,) + tuple(m[1:]): c}) * (pv[0] ** (k % 2)) * (SR ** (k // 2))
    return out


SY = sp.symbols(NAMES)


def parse(s):
    e = sp.Poly(sp.sympify(s, locals=dict(zip(NAMES, SY))), *SY)
    return QC.from_dict({m: flint.fmpq(int(c.p), int(c.q)) for m, c in zip(e.monoms(), e.coeffs())})


def fq(x):
    x = sp.Rational(x)
    return flint.fmpq(int(x.p), int(x.q))


# =============================================================================================
# B. faces: own Gram columns, own cost columns, own relation spaces
# =============================================================================================
class Face:
    def __init__(self, label, cubfile, budfile, sources, Binv, contacts):
        self.label = label
        C = json.load(open(os.path.join(BASE, cubfile)))
        D = json.load(open(os.path.join(BASE, budfile)))
        self.src = list(sources)
        m = self.m = len(sources)
        check(f'B[{label}] data source order equals the stated face order {sources}',
              list(D['source_order']) == self.src, str(D['source_order']))
        gens = L.generators(frozenset(contacts))
        check(f'B[{label}] contact set {sorted(contacts)}: generators {gens}; the {m} Gram sources are among them',
              set(self.src) <= set(gens) and list(C['contacts']) == sorted(contacts))
        self.gens = gens
        F = [[parse(s) for s in par] for par in C['features']]
        Tm = sp.Matrix([[sp.Rational(v) for v in row] for row in D['T']])
        Ti = Tm.inv()
        g = [[sum((fq(Ti[i, k]) * F[par][k] for k in range(m)), ZERO) for i in range(m)] for par in range(2)]
        gs = [[[gi.compose(*[s * Y[ix] for ix, s in sm]) for gi in g[par]] for sm in C['left_four']] for par in range(2)]
        self.entries = [(i, k) for i in range(m) for k in range(i, m)]
        self.Lcol = {}
        for (i, k) in self.entries:
            Mik = sum((gs[par][j][i] * gs[par][j][k] for par in range(2) for j in range(4)), ZERO)
            self.Lcol[(i, k)] = redN(NP * Mik) * (1 if i == k else 2)
        self.Ccol = {f: redN(2 * PHI[f]) for f in self.src}
        self.Binv = Binv
        # Delta(G) = Mlin(G) - 2 Phi_{Binv diag G}: matrix over monomials x entries
        cols_ = []
        for (i, k) in self.entries:
            col = self.Lcol[(i, k)]
            if i == k:
                for fi, f in enumerate(self.src):
                    if Binv[fi][i]:
                        col = col - fq(Binv[fi][i]) * self.Ccol[f]
            cols_.append(col.to_dict())
        mons = sorted(set(mm for c in cols_ for mm in c))
        rows = [[Fr(int(c.get(mm, 0).p), int(c.get(mm, 0).q)) if mm in c else Fr(0) for c in cols_] for mm in mons]
        Rr, piv = L.rref(rows)
        self.R = [r for r in Rr[:len(piv)]]
        self.nmons = len(mons)
        # ranks of the Gram map alone and of the costs
        gm = [self.Lcol[e_].to_dict() for e_ in self.entries]
        cm = [self.Ccol[f].to_dict() for f in self.src]
        mons2 = sorted(set(mm for c in gm + cm for mm in c))
        GR = L.rank([[Fr(int(c[mm].p), int(c[mm].q)) if mm in c else Fr(0) for c in gm] for mm in mons2])
        CR = L.rank([[Fr(int(c[mm].p), int(c[mm].q)) if mm in c else Fr(0) for c in cm] for mm in mons2])
        GCR = L.rank([[Fr(int(c[mm].p), int(c[mm].q)) if mm in c else Fr(0) for c in gm + cm] for mm in mons2])
        self.ranks = (GR, CR, GCR, len(self.R))
        check(f'B[{label}] own Gram map: {len(mons2)} monomial rows, rank {GR}; costs rank {CR}; '
              f'costs inside the Gram image: {GCR == GR}; relation space (Delta) rank {len(self.R)}',
              GCR == GR and CR == m)
        if 'R' in D:
            Rs = [[Fr(sp.Rational(v).p, sp.Rational(v).q) for v in row] for row in D['R']]
            same = L.rank(self.R + Rs) == len(self.R) == L.rank(Rs)
            check(f'B[{label}] own relation space equals the row space of the supplement R matrix (cross-check of data)', same)

    def holds(self, G, cost):
        """Exact symbolic: R vec(G) = 0 and Binv diag(G) = cost (sympy, simplify)."""
        vec = [G[i, k] for (i, k) in self.entries]
        bad = [sp.simplify(sum((sp.Rational(r[t].numerator, r[t].denominator) * vec[t] for t in range(len(vec)) if r[t]), 0))
               for r in self.R]
        bad = [b for b in bad if b != 0]
        dg = [G[i, i] for i in range(self.m)]
        cst = [sp.simplify(sum(sp.Rational(self.Binv[fi][i]) * dg[i] for i in range(self.m)) - cost[fi]) for fi in range(self.m)]
        return (not bad) and all(c == 0 for c in cst), bad, cst

    def poly_check(self, Gnum, costnum):
        """Full polynomial identity at a rational point: Mlin(G) == 2 Phi_cost (mod N)."""
        lhs = ZERO
        for (i, k) in self.entries:
            if Gnum[i][k] != 0:
                lhs += fq(Gnum[i][k]) * self.Lcol[(i, k)]
        rhs = ZERO
        for fi, f in enumerate(self.src):
            if costnum[fi] != 0:
                rhs += fq(costnum[fi]) * self.Ccol[f]
        return lhs == rhs


ident = lambda m: [[1 if i == k else 0 for k in range(m)] for i in range(m)]

if BASE is None:
    print('MOC4_SUPP not set; cannot load face data'); sys.exit(2)

# ---------------- (R) inventory56 face: (q,j,c,e,b,a,d) = (f0,f2,f5,f7,f12,f29,f35)
F56 = Face('inv56', 'final_referee_inventory56_boundary17_cubic_reduction_20261023.json',
           'final_referee_inventory56_boundary17_budget_PSD7_20261023.json',
           [0, 2, 5, 7, 12, 29, 35], ident(7), {0, 1, 5, 7, 9, 12, 17})
# printed ten relations (rank-two PDF p.2), Gram order (Q,J,C,E,B,A,D) = 0..6
pos = {(i, k): t for t, (i, k) in enumerate(F56.entries)}


def rel(terms):
    r = [Fr(0)] * len(F56.entries)
    for (i, k), c in terms:
        r[pos[(min(i, k), max(i, k))]] += c
    return r


Q_, J_, C_, E_, B_, A_, D_ = range(7)
PRINTED = [rel([((Q_, J_), 1), ((Q_, B_), 1), ((J_, B_), 1), ((B_, D_), -1)]),
           rel([((Q_, C_), 1), ((Q_, B_), -1), ((C_, A_), 1), ((B_, A_), -1), ((B_, D_), 1), ((D_, D_), -1)]),
           rel([((Q_, E_), 1), ((Q_, B_), -1), ((E_, A_), 1), ((B_, A_), -1), ((B_, D_), 1), ((D_, D_), -1)]),
           rel([((Q_, D_), 1)]), rel([((J_, D_), 1)]), rel([((E_, D_), 1)]), rel([((A_, D_), 1)]),
           rel([((J_, C_), 1), ((J_, B_), -1), ((C_, A_), -1), ((C_, D_), -1), ((B_, A_), 1), ((D_, D_), 1)]),
           rel([((C_, E_), 1), ((C_, A_), 1), ((C_, D_), 1), ((E_, A_), -1)]),
           rel([((E_, B_), 1), ((E_, A_), -1), ((B_, A_), 1), ((B_, D_), 1)])]
check('B[inv56] the ten relations printed on rank-two PDF p.2 span exactly my relation space',
      L.rank(PRINTED) == 10 and L.rank(PRINTED + F56.R) == 10 == len(F56.R))

# ---------------- (F1) T24 face: (q,p,j,c,e,d,b,a) = (f0,f2,f59,f6,f53,f54,f55,f67); diag = (q,p,j,c,e,q+d,b,a)
BinvT24 = ident(8)
BinvT24[5][0] = -1           # d = diag5 - diag0
F24 = Face('T24', 'root_moc4_T24_cubic_reduction_20261022.json', 'root_moc4_T24_budget_PSD8_20261022.json',
           [0, 2, 59, 6, 53, 54, 55, 67], BinvT24, {3, 4, 7, 8, 10, 12})
# ---------------- (F2) T25 face: (a,e,p,j,b,c,q,d) = (f0,f1,f2,f5,f6,f39,f67,f70); diag = coefficients
F25 = Face('T25', 'final_referee_T25_cubic_reduction_20261022.json', 'final_referee_T25_budget_PSD8_20261022.json',
           [0, 1, 2, 5, 6, 39, 67, 70], ident(8), {2, 3, 5, 7, 8, 12})
print('face setup done %.1f s' % (time.time() - T0), flush=True)


# =============================================================================================
# C. constructions: exact identity + exact LDL^T certificate of PSD
# =============================================================================================
def ldl(G, order):
    """Exact LDL^T of the principal submatrix on `order`. A zero pivot is allowed only if the
    whole remaining row is zero. Returns (pivots, ok) with ok = (P G P^T == L D L^T)."""
    n = len(order)
    A = sp.Matrix(n, n, lambda i, k: G[order[i], order[k]])
    Lm = sp.eye(n)
    Dv = []
    S = A.copy()
    for c in range(n):
        piv = sp.simplify(S[c, c])
        Dv.append(piv)
        if piv == 0:
            if any(sp.simplify(S[c, k]) != 0 for k in range(c + 1, n)):
                return Dv, False
            continue
        for i in range(c + 1, n):
            Lm[i, c] = sp.simplify(S[i, c] / piv)
        for i in range(c + 1, n):
            for k in range(c + 1, n):
                S[i, k] = sp.simplify(S[i, k] - Lm[i, c] * piv * Lm[k, c])
    Dm = sp.diag(*Dv)
    ok = all(sp.simplify(x) == 0 for x in (A - Lm * Dm * Lm.T))
    return Dv, ok


def record(name, face, G, cost, order, region, pivots_expected=None, numpt=None):
    ok, bad, cst = face.holds(G, cost)
    check(f'C {name}: Gram G satisfies my relation space and Binv diag(G) = stated cost', ok,
          '' if ok else f'nonzero relations {bad[:3]} cost diff {cst}')
    piv, okl = ldl(G, order)
    pivs = [sp.factor(x) for x in piv]
    check(f'C {name}: exact LDL^T on rows {order}; pivots {pivs}', okl)
    if pivots_expected is not None:
        okp = all(sp.simplify(a - b) == 0 for a, b in zip(piv, pivots_expected))
        check(f'C {name}: pivots equal the stated slack expressions', okp)
    zero_rows = [i for i in range(face.m) if i not in order]
    check(f'C {name}: rows outside the LDL order are identically zero', all(G[i, k] == 0 for i in zero_rows for k in range(face.m)))
    if numpt is not None:
        Gn = [[sp.nsimplify(G[i, k].subs(numpt)) for k in range(face.m)] for i in range(face.m)]
        cn = [sp.nsimplify(sp.sympify(c).subs(numpt)) for c in cost]
        check(f'C {name}: full polynomial identity Mlin(G) = 2 Phi_cost (mod N) at rational point {numpt}',
              face.poly_check(Gn, cn))
    RES['constructions'][name] = {'face': face.label, 'sources': face.src, 'cost': [str(c) for c in cost],
                                  'pivots': [str(x) for x in pivs], 'region': region}


# ---- (R) rank-two: G4 = v v^T + d x w w^T, order (Q,J,C,E,B,A,D)
d, x = sp.symbols('d x', positive=True)
v = sp.Matrix([0, 2 * x, 0, 0, d, -2 * x, 0])
w = sp.Matrix([0, 0, 0, 0, 1, 0, 2])
G4 = v * v.T + d * x * w * w.T
record('R rank-two G4 = vv^T + dx ww^T', F56, G4, [0, 4 * x ** 2, 0, 0, d * x + d ** 2, 4 * x ** 2, 4 * d * x],
       [B_, J_, A_, D_], 'd, x > 0', numpt={d: sp.Rational(3, 2), x: sp.Rational(2, 5)})
# PSD here is also evident: two rank-one terms with nonnegative weights.

# ---- (F1) T24 (5): order of Gram positions (Q,P,J,C,E,D,B,A) = 0..7
a, b, e, j, xx = sp.symbols('a b e j x', real=True)
Q4, P4, J4, C4, E4, D4, B4, A4 = range(8)


def sym8(entries):
    G = sp.zeros(8)
    for (i, k), val in entries.items():
        G[i, k] = val
        G[k, i] = val
    return G


G5 = sym8({(J4, J4): j, (E4, E4): e, (D4, D4): 1, (B4, B4): b, (A4, A4): a, (B4, A4): 1,
           (E4, D4): xx, (E4, B4): xx, (E4, A4): -xx, (E4, J4): xx - 1})
Lc = (a + 1) * (b + 1) / (a * b - 1)
emin5 = (xx - 1) ** 2 / j + Lc * xx ** 2
record('F1 T24 (5) coherent [general x]', F24, G5, [0, 0, j, 0, e, 1, b, a], [J4, D4, B4, A4, E4],
       'ab > 1, j > 0; PSD iff e >= (1-x)^2/j + L x^2', pivots_expected=[j, 1, b, a - 1 / b, e - emin5],
       numpt={a: 3, b: 2, j: sp.Rational(1, 3), e: 7, xx: sp.Rational(2, 7)})
xs = 1 / (1 + j * Lc)
check('C F1 T24 (5): at x = 1/(1+jL), (1-x)^2/j + L x^2 = 1/(j + 1/L) (so e >= 1/(j+1/L) <=> j >= 1/e - 1/L)',
      sp.simplify(emin5.subs(xx, xs) - 1 / (j + 1 / Lc)) == 0)
check('C F1 T24 (5): homogeneous form j >= d^2/e - d(ab-d^2)/((a+d)(b+d)) at d = 1 equals j >= 1/e - 1/L',
      sp.simplify((1 / e - (a * b - 1) / ((a + 1) * (b + 1))) - (1 / e - 1 / Lc)) == 0)
G5j0 = sym8({(E4, E4): e, (D4, D4): 1, (B4, B4): b, (A4, A4): a, (B4, A4): 1, (E4, D4): 1, (E4, B4): 1, (E4, A4): -1})
record('F1 T24 (5) at j = 0 (x = 1)', F24, G5j0, [0, 0, 0, 0, e, 1, b, a], [D4, B4, A4, E4],
       'ab > 1; PSD iff e >= L (= K >= 0)', pivots_expected=[1, b, a - 1 / b, e - Lc],
       numpt={a: 3, b: 2, e: 5})
# double-SOC endpoint: (J,E) and (B,A) blocks; valid whenever ej >= 1 and ab >= 1 (x = 0 member of (5))
GdS = sym8({(J4, J4): j, (E4, E4): e, (E4, J4): -1, (D4, D4): 1, (B4, B4): b, (A4, A4): a, (B4, A4): 1})
record('F1 T24 double-SOC', F24, GdS, [0, 0, j, 0, e, 1, b, a], [J4, E4, D4, B4, A4],
       'ej >= 1, ab >= 1 (stated for ab = 1)', pivots_expected=[j, e - 1 / j, 1, b, a - 1 / b],
       numpt={a: 1, b: 1, j: 2, e: sp.Rational(1, 2)})
# ---- (F1) T24 (6): b = r^2, a = s^2, r != s, J.v = tau = (1-rs)/(r-s)
r, s = sp.symbols('r s', positive=True)
tau = (1 - r * s) / (r - s)
G6 = sym8({(J4, J4): j, (E4, E4): e, (E4, J4): -1, (D4, D4): 1, (B4, B4): r ** 2, (A4, A4): s ** 2,
           (B4, A4): r * s, (J4, B4): r * tau, (J4, A4): s * tau})
record('F1 T24 (6) parallel, order with B first (b > a or a = 0)', F24, G6, [0, 0, j, 0, e, 1, r ** 2, s ** 2],
       [E4, D4, B4, J4, A4], 'b = r^2 > 0, a = s^2, r != s; PSD iff j >= 1/e + tau^2',
       pivots_expected=[e, 1, r ** 2, j - 1 / e - tau ** 2, 0],
       numpt={r: 2, s: sp.Rational(1, 3), e: 3, j: 5})
piv, okl = ldl(G6, [E4, D4, A4, J4, B4])
check('C F1 T24 (6): with a > b the same G is PSD (LDL order with A first), pivots e, 1, s^2, j - 1/e - tau^2, 0',
      okl and all(sp.simplify(p_ - q_) == 0 for p_, q_ in zip(piv, [e, 1, s ** 2, j - 1 / e - tau ** 2, 0])))
check('C F1 T24 (6): tau^2 = b(1-w)^2/(b-w)^2 = a(1-w)^2/(a-w)^2 with w = rs = sqrt(ab) (formula symmetric in a, b)',
      sp.simplify(tau ** 2 - r ** 2 * (1 - r * s) ** 2 / (r ** 2 - r * s) ** 2) == 0 and
      sp.simplify(tau ** 2 - s ** 2 * (1 - r * s) ** 2 / (s ** 2 - r * s) ** 2) == 0)
# homogeneous form of (6): j >= d^2/e + b (d - w)^2/(b - w)^2 ; scaling check
dd = sp.symbols('dd', positive=True)
lhs6 = dd * (1 / (e / dd) + (r ** 2 / dd) * (1 - r * s / dd) ** 2 / ((r ** 2 - r * s) / dd) ** 2)
check('C F1 T24 (6): d * [threshold at normalised coefficients] = d^2/e + b(d-w)^2/(b-w)^2 (printed homogeneous form)',
      sp.simplify(lhs6 - (dd ** 2 / e + r ** 2 * (dd - r * s) ** 2 / (r ** 2 - r * s) ** 2)) == 0)

# ---- (F2) T25 partial-four: positions (A,E,P,J,B,C,Q,D) = 0..7, sources (f0,f1,f2,f5,f6,f39,f67,f70)
A5_, E5_, P5_, J5_, B5_, C5_, Q5_, D5_ = range(8)
a, b, c, q, t, z, kap = sp.symbols('a b c q t z kappa', real=True)
rr = 1 - t
hb = b - t ** 2 / a
Gbase = sp.diag(sp.Matrix([[a, t, 0], [t, b, z], [0, z, q]]), sp.Matrix([[c]]), sp.Matrix([[1]]))   # basis (A,B,Q,C,D)


def T25G(kappa_expr, P_C):
    K = sp.zeros(8, 5)
    K[A5_, 0] = 1; K[B5_, 1] = 1; K[Q5_, 2] = 1; K[C5_, 3] = 1; K[D5_, 4] = 1
    den = a * b - t ** 2
    K[E5_, 0] = t * (b + t) / den; K[E5_, 1] = -t * (a + t) / den; K[E5_, 4] = t
    K[P5_, 0] = -kappa_expr * t / a; K[P5_, 1] = kappa_expr; K[P5_, 2] = -kappa_expr; K[P5_, 3] = P_C
    return (K * Gbase * K.T).applyfunc(sp.simplify)


et = t ** 2 + t ** 2 * (a + b + 2 * t) / (a * b - t ** 2)
Ggen = T25G(kap, rr / c)
check('C F2 T25: E.A = t, E.B = -t, E.D = t, |E|^2 = e(t), |V|^2 = hb + q - 2z, |P|^2 = r^2/c + kappa^2(hb+q-2z)',
      all(sp.simplify(u_) == 0 for u_ in [Ggen[E5_, A5_] - t, Ggen[E5_, B5_] + t, Ggen[E5_, D5_] - t, Ggen[E5_, E5_] - et,
                                           Ggen[P5_, P5_] - (rr ** 2 / c + kap ** 2 * (hb + q - 2 * z))]))
okfree, badfree, _ = F25.holds(Ggen, [a, et, Ggen[P5_, P5_], 0, b, c, q, 1])
check('C F2 T25 with kappa free: exactly one relation fails, and it is kappa (hb + q - 2z) = r - z (forces (12))',
      (not okfree) and len(badfree) == 1 and
      sp.simplify(sp.solve(badfree[0], kap)[0] - (rr - z) / (hb + q - 2 * z)) == 0, str(badfree))
k12 = (rr - z) / (hb + q - 2 * z)
G12 = Ggen.subs(kap, k12).applyfunc(sp.simplify)
record('F2 T25 Lemma 1 general (free z, kappa from (12))', F25, G12,
       [a, et, rr ** 2 / c + (rr - z) ** 2 / (hb + q - 2 * z), 0, b, c, q, 1], [A5_, B5_, Q5_, C5_, D5_, E5_, P5_],
       'a, c > 0, ab > t^2, hb + q - 2z != 0; PSD iff q hb >= z^2 (G = K Gbase K^T)',
       numpt={a: 2, b: 3, c: 5, q: 7, t: sp.Rational(1, 2), z: sp.Rational(1, 3)})
pb, okb = ldl(Gbase, [0, 1, 2, 3, 4])
check('C F2 T25: Gbase pivots a, hb, q - z^2/hb, c, 1 (exact LDL^T)',
      okb and all(sp.simplify(p_ - q_) == 0 for p_, q_ in zip(pb, [a, hb, q - z ** 2 / hb, c, 1])))
# (9): z = r, P = (r/c) C ; (10): z = wr, q = w^2 r^2/hb, kappa = (r-z)/(hb+q-2z)
G9 = T25G(0, rr / c).subs(z, rr)
record('F2 T25 Lemma 1 (9)', F25, G9, [a, et, rr ** 2 / c, 0, b, c, q, 1], [A5_, B5_, Q5_, C5_, D5_, E5_, P5_],
       'q hb >= r^2', numpt={a: 2, b: 3, c: 5, q: 7, t: sp.Rational(1, 2)})
ww = sp.symbols('w', real=True)
q10 = ww ** 2 * rr ** 2 / hb
k10 = (rr - ww * rr) / (hb + q10 - 2 * ww * rr)
check('C F2 T25 (10): hb + q - 2z = (hb - wr)^2/hb and kappa^2 (hb+q-2z) = hb r^2 (1-w)^2/(hb-wr)^2',
      sp.simplify(hb + q10 - 2 * ww * rr - (hb - ww * rr) ** 2 / hb) == 0 and
      sp.simplify(k10 ** 2 * (hb + q10 - 2 * ww * rr) - hb * rr ** 2 * (1 - ww) ** 2 / (hb - ww * rr) ** 2) == 0)
check('C F2 T25 (10): at z = wr, q = w^2 r^2/hb the third Gbase pivot q - z^2/hb is 0 (PSD boundary)',
      sp.simplify(q10 - (ww * rr) ** 2 / hb) == 0)
G10 = Ggen.subs(kap, k10).subs({z: ww * rr}).subs(q, q10).applyfunc(sp.simplify)
record('F2 T25 Lemma 1 (10)', F25, G10,
       [a, et, rr ** 2 / c + hb * rr ** 2 * (1 - ww) ** 2 / (hb - ww * rr) ** 2, 0, b, c, q10, 1],
       [A5_, B5_, Q5_, C5_, D5_, E5_, P5_], 'q = w^2 r^2/hb, hb > wr',
       numpt={a: 2, b: 3, c: 5, t: sp.Rational(1, 2), ww: sp.Rational(1, 3)})

# =============================================================================================
# D. orbits and transports (exact, lib group action)
# =============================================================================================
VAL = {tuple(vv): jj for jj, vv in enumerate(L.F_VALUES)}
GPERM = sorted(set(tuple(VAL[tuple(L.g_act_values(g_, list(L.F_VALUES[jj])))] for jj in range(82)) for g_ in L.G_ELEMENTS))
check('D G permutes the 82 f_j (1152 distinct permutations)', len(GPERM) == 1152)
REPS = {'SOC': (1, 9, 22), 'H': (3, 15, 26), 'chain': (2, 0, 12, 70), 'four': (0, 6, 10, 54), 'rank2': (2, 29, 12, 35)}


def keyf(t_, fam):
    if fam == 'SOC':
        return (frozenset(t_[:2]), t_[2])
    if fam == 'H':
        return tuple(t_)
    if fam == 'chain':
        return (t_[0], frozenset(t_[1:3]), t_[3])
    return (frozenset(t_[:2]), t_[2], t_[3])    # four (a,b;c;d) and rank2 (leaf,leaf;centre;bad)


ORB = {fam: set(keyf(tuple(m_[x_] for x_ in rep), fam) for m_ in GPERM) for fam, rep in REPS.items()}
for fam in ORB:
    print(f'orbit {fam}: {len(ORB[fam])}')
USES = [('SOC', (2, 59, 54), 'T24 (4) pj >= d^2'), ('H', (53, 59, 54), 'T24 (4) e(j-d) >= d^2'),
        ('four', (67, 55, 53, 54), 'T24 (4) K >= 0'),
        ('SOC', (2, 5, 70), 'T25 (4) pj >= 1'), ('H', (39, 2, 70), 'T25 (4) c(p-1) >= 1'),
        ('four', (0, 6, 1, 70), 'T25 (5) first'), ('four', (6, 67, 39, 70), 'T25 (5) second'),
        ('SOC', (2, 5, 54), 'T41 bad54 (1) pj >= 1'), ('H', (15, 2, 54), 'T41 bad54 (1) a(p-1) >= 1'),
        ('H', (3, 5, 54), 'T41 bad54 (1) b(j-1) >= 1'), ('chain', (2, 0, 15, 54), 'T41 bad54 (1) pae >= a+e'),
        ('chain', (5, 0, 3, 54), 'T41 bad54 (1) jbe >= b+e'),
        ('chain', (12, 0, 2, 35), 'inv56 boundary17 (1) chain leaves (0,2)'),
        ('chain', (12, 7, 29, 35), 'inv56 boundary17 (1) chain leaves (7,29)'),
        ('four', (7, 0, 5, 35), 'inv56 boundary17 (3) / large_product (2): qe(c-d) >= d^2(e+q+c+d)'),
        ('H', (29, 12, 35), 'rank-two p.3 old-H'), ('H', (2, 12, 35), 'rank-two p.3 old-H'),
        ('rank2', (2, 29, 12, 35), 'rank-two representative')]
for fam, t_, where in USES:
    ok = keyf(t_, fam) in ORB[fam]
    RES['orbits'][f'{fam}{t_}'] = {'where': where, 'in_orbit_with_roles': ok}
    check(f'D {fam} tuple {t_} [{where}] lies in the proved orbit of {REPS[fam]} with the printed roles', ok)
# the a(p-1) >= 1 instance of T41 bad54 (aux f15 coefficient a, main f2 coefficient p):
# T24 core transports. Roles: J (f59), E (f53), {B, A} (f55, f67, symmetric), D (f54), SOC partner P (f2)
CORE = (59, 53, 55, 67, 54, 2)


def transports(img, use_p=True, setAB=True):
    """exists g in G with g(59)=img0, g(53)=img1, {g(55),g(67)} = {img2,img3} (or ordered), g(54)=img4 [, g(2)=img5]"""
    for m_ in GPERM:
        if m_[59] != img[0] or m_[53] != img[1] or m_[54] != img[4]:
            continue
        if setAB:
            if {m_[55], m_[67]} != {img[2], img[3]}:
                continue
        elif (m_[55], m_[67]) != (img[2], img[3]):
            continue
        if use_p and len(img) > 5 and m_[2] != img[5]:
            continue
        return True
    return False


TR = [((1, 53, 4, 16, 17), 'T41 bad17 p.2 / large_product p.1 [59,53,55,67,54] -> [1,53,4,16,17]'),
      ((5, 2, 0, 67, 40), 'T27 (paper T_28) p.2 (f59,f53,f55,f67,f54) -> (f5,f2,f0,f67,f40)'),
      ((5, 2, 67, 0, 40), 'T38 (paper T_39) p.1 (f59,f53,f55,f67,f54) -> (f5,f2,f67,f0,f40)'),
      ((45, 5, 0, 7, 35, 1), 'T36 (paper T_37) p.2 (67,55,53,54,59,2) -> (7,0,5,35,45,1)'),
      ((59, 5, 6, 0, 34, 1), 'T43 (paper T_44) p.2 (59,53,55,67,54,2) -> (59,5,6,0,34,1)'),
      ((5, 3, 67, 0, 54, 2), 'T23 g17 (paper T_24) p.2 (59,53,55,67,54,2) -> (5,3,67,0,54,2)'),
      ((2, 39, 6, 67, 70),'T25 (6),(7): j->p (f2), e->c (f39), {b,a}->{b,q} (f6,f67), f54->f70')]
for img, where in TR:
    ok_ordered = transports(img, setAB=False)
    ok_set = transports(img, setAB=True)
    RES['orbits'][f'T24core->{img}'] = {'where': where, 'ordered': ok_ordered, 'unordered_AB': ok_set}
    check(f'D T24 core transport {where}: some g in G realises it (roles J, E, D[, P] exact, {{B,A}} as a set); '
          f'ordered (B,A) also realised: {ok_ordered}', ok_set)
# gram faces vs contact sets
check('D T24 face generators = {0,1,2,6,53,54,55,59,67,70} (doc p.1)', F24.gens == [0, 1, 2, 6, 53, 54, 55, 59, 67, 70], str(F24.gens))
check('D T25 face generators = {0,1,2,5,6,39,67,70} (doc p.1)', F25.gens == [0, 1, 2, 5, 6, 39, 67, 70], str(F25.gens))
check('D inv56 seven-contact face generators = {0,2,5,7,12,29,35} (doc p.1)', F56.gens == [0, 2, 5, 7, 12, 29, 35], str(F56.gens))
check('D only f54, f70 (T24), f70 (T25), f35 (inv56) are non-minor among the face generators',
      [g_ for g_ in F24.gens if not L.F_IS_MINOR[g_]] == [54, 70] and [g_ for g_ in F25.gens if not L.F_IS_MINOR[g_]] == [70]
      and [g_ for g_ in F56.gens if not L.F_IS_MINOR[g_]] == [35])
# inventory56 chain allocation (boundary17 p.1, large_product p.1)
qq, jj_, ee, aa, D_s = sp.symbols('q j e a d', positive=True)
al, be = qq * jj_ / (qq + jj_), ee * aa / (ee + aa)
d1, d2 = D_s * al / (al + be), D_s * be / (al + be)
alloc = d1 ** 2 * (1 / qq + 1 / jj_) + d2 ** 2 * (1 / ee + 1 / aa)
Bch = D_s ** 2 * (qq + jj_) * (ee + aa) / (qq * jj_ * (ee + aa) + ee * aa * (qq + jj_))
check('D inv56 (1): split d = d1 + d2 in ratio qj/(q+j) : ea/(e+a); centre allocations d1^2(1/q+1/j) + d2^2(1/e+1/a) = Bch',
      sp.simplify(alloc - Bch) == 0 and sp.simplify(d1 + d2 - D_s) == 0)
check('D inv56 large_product (1): b(qj/(q+j) + ea/(e+a)) >= d^2 is b >= Bch',
      sp.simplify(D_s ** 2 / (al + be) - Bch) == 0)
# rank-two pointwise form from the coefficient form: inf over t of t X + Y/(4t) = sqrt(X Y)
X_, Y_, t_s = sp.symbols('X Y t', positive=True)
expr = t_s * X_ + Y_ / (4 * t_s)
tm = sp.sqrt(Y_ / (4 * X_))
check('D rank-two (2): min_t>0 [t(f2+f29) + f12/(4t)] = sqrt(f12 (f2+f29)), attained at t = sqrt(f12/(4(f2+f29)))',
      sp.simplify(expr.subs(t_s, tm) - sp.sqrt(X_ * Y_)) == 0 and sp.simplify(sp.diff(expr, t_s).subs(t_s, tm)) == 0)
print('exact part done, %.1f s, failures %d' % (time.time() - T0, NFAIL), flush=True)
if os.environ.get('PART2A_EXACT_ONLY'):
    sys.exit(0 if NFAIL == 0 else 1)

# =============================================================================================
# E. numeric attack (NUMERIC evidence only)
# =============================================================================================
FV = np.array(L.F_VALUES, dtype=float)
SGNA = np.array(L.SGN, dtype=float)
PIDX = np.array(L.PERMS)
E4m = np.eye(4)


def qmul_np(a_, b_):
    a0, a1, a2, a3 = a_
    b0, b1, b2, b3 = b_
    return np.array([a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3, a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
                     a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1, a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0])


def O_np(zv):
    p_ = zv[:4] / np.linalg.norm(zv[:4])
    q_ = zv[4:] / np.linalg.norm(zv[4:])
    return np.array([qmul_np(qmul_np(p_, E4m[k]), q_) for k in range(4)]).T


def wts(O):
    return SGNA * O[0, PIDX[:, 0]] * O[1, PIDX[:, 1]] * O[2, PIDX[:, 2]] * O[3, PIDX[:, 3]]


def multistart(fun, starts, seed, maxiter=2500):
    rng = np.random.default_rng(seed)
    best = (np.inf, None)
    for _ in range(starts):
        rres = minimize(fun, rng.standard_normal(8), method='Nelder-Mead',
                        options={'xatol': 1e-11, 'fatol': 1e-15, 'maxiter': maxiter, 'maxfev': maxiter})
        if rres.fun < best[0]:
            best = (float(rres.fun), rres.x.copy())
    return best


def cost_fun(terms):
    idx = np.array([t_[0] for t_ in terms]); co = np.array([t_[1] for t_ in terms], dtype=float)
    M = FV[idx]
    return lambda zv: co @ (M @ wts(O_np(zv)))


mpmath.mp.dps = 50


def mp_cost(terms, zv):
    zz = [mpmath.mpf(float(u)) for u in zv]
    p_, q_ = zz[:4], zz[4:]
    n1 = mpmath.sqrt(sum(u * u for u in p_)); n2 = mpmath.sqrt(sum(u * u for u in q_))
    p_ = [u / n1 for u in p_]; q_ = [u / n2 for u in q_]

    def qm(a_, b_):
        a0, a1, a2, a3 = a_; b0, b1, b2, b3 = b_
        return [a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3, a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
                a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1, a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0]
    cl = [qm(qm(p_, [1 if i == k else 0 for i in range(4)]), q_) for k in range(4)]
    O = [[cl[k][i] for k in range(4)] for i in range(4)]
    wv = [sg * O[0][s_[0]] * O[1][s_[1]] * O[2][s_[2]] * O[3][s_[3]] for s_, sg in zip(L.PERMS, L.SGN)]
    return sum(mpmath.mpf(float(co)) * sum(vv * w_ for vv, w_ in zip(L.F_VALUES[f], wv)) for f, co in terms)


def attack(name, sampler, nsamp, starts, seed, scale=1.0):
    rng = np.random.default_rng(seed)
    vals = []
    for i in range(nsamp):
        terms = sampler(rng, i)
        if terms is None:
            continue
        terms = [(f, co * (scale if tag == 'bnd' else 1.0)) for f, co, tag in terms]
        vmin, zbest = multistart(cost_fun(terms), starts, seed=seed * 1000 + i)
        vals.append((vmin, terms, zbest))
    worst = min(vals, key=lambda u: u[0])
    mpv = mp_cost(worst[1], worst[2])
    out = {'samples': len(vals), 'starts_each': starts, 'seed': seed, 'scale_on_threshold_coefficient': scale,
           'min': worst[0], 'mp50_at_worst': mpmath.nstr(mpv, 8),
           'worst_terms': [(int(f), float(co)) for f, co in worst[1]], 'max_of_minima': max(u[0] for u in vals)}
    RES['numeric'][name] = out
    print(f'E {name}: {len(vals)} coefficient points x {starts} starts, min = {worst[0]:.3e} '
          f'(mpmath50 {mpmath.nstr(mpv, 6)}) at {[(f, round(co, 4)) for f, co in worst[1]]}', flush=True)
    return out


def pick(rng, lo, hi):   # log-uniform, with occasional extreme values
    u = rng.uniform(math.log(lo), math.log(hi))
    return math.exp(u)


# (R) boundary b = d/4 + d^2/(4 min(a,j)), d = 1; threshold coefficient is b (tag 'bnd')
def sR(rng, i):
    a_, j_ = pick(rng, 0.02, 50), pick(rng, 0.02, 50)
    if i % 4 == 0:
        j_ = a_
    return [(12, 0.25 + 1 / (4 * min(a_, j_)), 'bnd'), (2, j_, ''), (29, a_, ''), (35, 1.0, '')]


attack('R rank-two boundary (1)', sR, 24, 4, 11)


# pointwise (2): Phi35 + Phi12/4 + sqrt(Phi12 (Phi2 + Phi29))
def pw2(zv):
    wv = wts(O_np(zv))
    p12, p2, p29, p35 = FV[12] @ wv, FV[2] @ wv, FV[29] @ wv, FV[35] @ wv
    return p35 + p12 / 4 + math.sqrt(max(0.0, p12) * max(0.0, p2 + p29))


vpw, zpw = multistart(pw2, 32, seed=7)
RES['numeric']['R rank-two pointwise (2)'] = {'starts': 32, 'seed': 7, 'min': vpw}
print(f'E R rank-two pointwise (2): 32 starts, min = {vpw:.3e}', flush=True)


# (F1) T24 (5): j = 1/e - 1/L (when > 0) ; tag j as threshold
def s5(rng, i):
    a_ = pick(rng, 0.05, 30); b_ = pick(rng, 0.05, 30)
    if a_ * b_ <= 1.0:
        b_ = (1.0 + pick(rng, 0.001, 5)) / a_
    Lc_ = (a_ + 1) * (b_ + 1) / (a_ * b_ - 1)
    e_ = Lc_ * rng.uniform(0.02, 0.999) if i % 5 else Lc_ * 0.999999
    j_ = 1 / e_ - 1 / Lc_
    return [(59, j_, 'bnd'), (53, e_, ''), (54, 1.0, ''), (55, b_, ''), (67, a_, '')]


attack('F1 T24 (5) boundary', s5, 24, 4, 21)


def s6(rng, i):
    wv_ = rng.uniform(0.02, 0.98)
    rat = pick(rng, 1.001, 200)              # b/a
    b_ = wv_ * math.sqrt(rat); a_ = wv_ / math.sqrt(rat)
    if i % 2:
        a_, b_ = b_, a_                      # symmetric use (a > b), as in T27/T38/T36/T43/T25
    e_ = pick(rng, 0.02, 50)
    j_ = 1 / e_ + max(a_, b_) * (1 - wv_) ** 2 / (max(a_, b_) - wv_) ** 2
    return [(59, j_, 'bnd'), (53, e_, ''), (54, 1.0, ''), (55, b_, ''), (67, a_, '')]


attack('F1 T24 (6) boundary (both orientations)', s6, 24, 4, 22)


def sdS(rng, i):
    a_ = pick(rng, 0.05, 20); e_ = pick(rng, 0.05, 20)
    return [(59, 1 / e_, 'bnd'), (53, e_, ''), (54, 1.0, ''), (55, 1 / a_, ''), (67, a_, '')]


attack('F1 T24 double-SOC ab = 1, ej = 1', sdS, 16, 4, 23)


# (F2) T25 (9) and (10); threshold coefficient p (f2)
def s9(rng, i):
    a_, b_, c_ = pick(rng, 0.05, 20), pick(rng, 0.05, 20), pick(rng, 0.05, 20)
    t_ = rng.uniform(0, 0.98) if i % 4 else 0.0
    if a_ * b_ <= t_ * t_:
        a_ = (t_ * t_ + 0.1) / b_
    h_ = b_ - t_ * t_ / a_; r_ = 1 - t_
    q_ = r_ * r_ / h_
    e_ = t_ * t_ + t_ * t_ * (a_ + b_ + 2 * t_) / (a_ * b_ - t_ * t_)
    return [(0, a_, ''), (1, e_, ''), (2, r_ * r_ / c_, 'bnd'), (6, b_, ''), (39, c_, ''), (67, q_, ''), (70, 1.0, '')]


attack('F2 T25 Lemma 1 (9) boundary q hb = r^2', s9, 24, 4, 31)


def s10(rng, i):
    a_, b_, c_ = pick(rng, 0.05, 20), pick(rng, 0.05, 20), pick(rng, 0.05, 20)
    t_ = rng.uniform(0, 0.98)
    if a_ * b_ <= t_ * t_:
        a_ = (t_ * t_ + 0.1) / b_
    h_ = b_ - t_ * t_ / a_; r_ = 1 - t_
    w_ = rng.uniform(0, 1) if i % 4 else 0.0
    if h_ <= w_ * r_:
        return None
    q_ = w_ * w_ * r_ * r_ / h_
    e_ = t_ * t_ + t_ * t_ * (a_ + b_ + 2 * t_) / (a_ * b_ - t_ * t_)
    p_ = r_ * r_ / c_ + h_ * r_ * r_ * (1 - w_) ** 2 / (h_ - w_ * r_) ** 2
    return [(0, a_, ''), (1, e_, ''), (2, p_, 'bnd'), (6, b_, ''), (39, c_, ''), (67, q_, ''), (70, 1.0, '')]


attack('F2 T25 Lemma 1 (10) boundary', s10, 24, 4, 32)


# (F3) T41 bad17 transported: cost p f1 + j f53 + a f4 + b f16 + f17 (q, c, e = 0)
def s3(rng, i):
    a_ = pick(rng, 0.05, 30); b_ = pick(rng, 0.05, 30)
    if a_ * b_ < 1:
        b_ = (1 + pick(rng, 1e-4, 5)) / a_
    Rr_ = (a_ + 1) * (b_ + 1); A_ = a_ * b_
    j_ = pick(rng, 0.05, 30)
    p_ = max(0.0, 1 / j_ - (A_ - 1) / Rr_)
    if p_ == 0.0:
        j_ = Rr_ / (A_ - 1); p_ = 0.0
    return [(1, p_, 'bnd'), (53, j_, ''), (4, a_, ''), (16, b_, ''), (17, 1.0, '')]


attack('F3 T41 bad17 (3) transported boundary', s3, 24, 4, 41)


def s4(rng, i):
    z_ = rng.uniform(0.02, 0.98); rho = pick(rng, 1.001, 300)
    a_, b_ = z_ * math.sqrt(rho), z_ / math.sqrt(rho)   # a > b
    j_ = pick(rng, 0.05, 30)
    k0 = (1 - z_) ** 2 / (a_ + b_ - 2 * z_)
    return [(1, 1 / j_ + k0, 'bnd'), (53, j_, ''), (4, a_, ''), (16, b_, ''), (17, 1.0, '')]


attack('F3 T41 bad17 (4) transported boundary', s4, 24, 4, 42)


# (F4) inventory56 (1) chain-sum boundary and (3) four boundary
def sch(rng, i):
    q_, j_, e_, a_ = (pick(rng, 0.05, 30) for _ in range(4))
    Bc = (q_ + j_) * (e_ + a_) / (q_ * j_ * (e_ + a_) + e_ * a_ * (q_ + j_))
    return [(0, q_, ''), (2, j_, ''), (7, e_, ''), (29, a_, ''), (12, Bc, 'bnd'), (35, 1.0, '')]


attack('F4 inv56 (1) two-chain boundary b = Bch', sch, 20, 4, 51)


def sfour(rng, i):
    q_, e_ = pick(rng, 0.05, 30), pick(rng, 0.05, 30)
    if q_ * e_ <= 1:
        e_ = (1 + pick(rng, 1e-3, 5)) / q_
    c_ = (q_ + 1) * (e_ + 1) / (q_ * e_ - 1)          # boundary qe(c-1) = q+e+c+1
    return [(0, q_, ''), (7, e_, ''), (5, c_, 'bnd'), (35, 1.0, '')]


attack('F4 inv56 (3) four-coefficient boundary', sfour, 16, 4, 52)

# sensitivity: shrink the threshold coefficient by 10 percent (negative = the test detects violations)
SENS = {}
for nm, smp, sd in [('R rank-two', sR, 11), ('F1 T24 (5)', s5, 21), ('F1 T24 (6)', s6, 22), ('F2 T25 (9)', s9, 31),
                    ('F2 T25 (10)', s10, 32), ('F3 T41 (4)', s4, 42), ('F4 inv56 (1)', sch, 51), ('F4 inv56 (3)', sfour, 52)]:
    for sc in (0.9, 0.5, 0.0):
        o = attack(nm + f' threshold x{sc} (outside region)', smp, 6, 3, sd + 100 + int(10 * sc), scale=sc)
        SENS[f'{nm} x{sc}'] = o['min']
RES['numeric']['sensitivity_summary'] = SENS

# =============================================================================================
# output
# =============================================================================================
PR = 'PROVED HERE (EXACT + PROOF-READ)'
INVENTORY = [
    {'id': 'R', 'name': 'rank-two capacity', 'generators': {'b': 12, 'j': 2, 'a': 29, 'd (non-minor)': 35},
     'region': 'a, j, d > 0, b >= d/4 + d^2/(4 min(a, j)); pointwise: Phi35 + Phi12/4 + |x12| sqrt(Phi2 + Phi29) >= 0',
     'status': PR, 'constructions': ['R rank-two G4 = vv^T + dx ww^T'],
     'numeric': ['R rank-two boundary (1)', 'R rank-two pointwise (2)'],
     'docs': ['final_referee_inventory56_rank_two_capacity_20261023 p.1-2 (statement, proof)',
              'final_referee_inventory56_boundary17_uniform_20261023 p.1-2 (condition (2), B_r2)']},
    {'id': 'F1-5', 'name': 'T24 coherent capacity (5)', 'generators': {'j': 59, 'e': 53, 'd (non-minor)': 54, 'b': 55, 'a': 67},
     'region': 'ab > d^2, e > 0, j >= d^2/e - d(ab - d^2)/((a+d)(b+d)), j >= 0; minor coefficients q (f0), p (f2), c (f6) >= 0 arbitrary',
     'status': PR, 'constructions': ['F1 T24 (5) coherent [general x]', 'F1 T24 (5) at j = 0 (x = 1)'],
     'numeric': ['F1 T24 (5) boundary'], 'docs': ['root_moc4_T24_uniform_20261022 p.2 Lemma 1 (5)', 'uses: see transports']},
    {'id': 'F1-6', 'name': 'T24 parallel capacity (6)', 'generators': {'j': 59, 'e': 53, 'd (non-minor)': 54, 'b': 55, 'a': 67},
     'region': 'ab < d^2, a != b (printed as b > sqrt(ab); the threshold is symmetric in a, b), e > 0, j >= d^2/e + b(d - sqrt(ab))^2/(b - sqrt(ab))^2',
     'status': PR, 'constructions': ['F1 T24 (6) parallel, order with B first (b > a or a = 0)'],
     'numeric': ['F1 T24 (6) boundary (both orientations)'], 'docs': ['root_moc4_T24_uniform_20261022 p.2 Lemma 1 (6)']},
    {'id': 'F1-dS', 'name': 'T24 double-SOC', 'generators': {'j': 59, 'e': 53, 'd (non-minor)': 54, 'b': 55, 'a': 67},
     'region': 'ab = d^2 (proved for ab >= d^2), ej >= d^2', 'status': PR, 'constructions': ['F1 T24 double-SOC'],
     'numeric': ['F1 T24 double-SOC ab = 1, ej = 1'], 'docs': ['root_moc4_T24_uniform_20261022 p.2']},
    {'id': 'F2-9', 'name': 'T25 partial-four (9)',
     'generators': {'a': 0, 'e': 1, 'p': 2, 'j': 5, 'b': 6, 'c': 39, 'q': 67, 'd (non-minor)': 70},
     'region': 'd = 1 normalised; a, b, c > 0, 0 <= t < 1, ab > t^2, q hb >= r^2, e >= e(t), p >= r^2/c, j >= 0',
     'status': PR, 'constructions': ['F2 T25 Lemma 1 (9)'], 'numeric': ['F2 T25 Lemma 1 (9) boundary q hb = r^2'],
     'docs': ['final_referee_T25_uniform_20261022 p.2 Lemma 1; used p.5 (18), (20)']},
    {'id': 'F2-10', 'name': 'T25 partial-four (10)',
     'generators': {'a': 0, 'e': 1, 'p': 2, 'j': 5, 'b': 6, 'c': 39, 'q': 67, 'd (non-minor)': 70},
     'region': 'd = 1; a, b, c > 0, 0 <= t < 1, ab > t^2, q = w^2 r^2/hb, 0 <= w <= 1, hb > wr, e >= e(t), p >= r^2/c + hb r^2 (1-w)^2/(hb - wr)^2',
     'status': PR, 'constructions': ['F2 T25 Lemma 1 (10)', 'F2 T25 Lemma 1 general (free z, kappa from (12))'],
     'numeric': ['F2 T25 Lemma 1 (10) boundary'], 'docs': ['final_referee_T25_uniform_20261022 p.2 Lemma 1; used p.5 (21)']},
    {'id': 'F2-67', 'name': 'T25 (6), (7): transported T24 (5), (6), double-SOC',
     'generators': {'p': 2, 'c': 39, '{b, q}': [6, 67], 'd (non-minor)': 70},
     'region': 'bq > 1: p >= 1/c - 1/E2; bq < 1, b != q: p >= 1/c + b(1-sqrt(bq))^2/(b-sqrt(bq))^2; bq = 1: pc >= 1',
     'status': 'REDUCES TO PROVED ATOMS (F1-5, F1-6, F1-dS via an exact group element)', 'docs': ['final_referee_T25_uniform_20261022 p.2']},
    {'id': 'F3', 'name': 'T41 bad17 (3), (4): transported T24 (5), (6), double-SOC',
     'generators': {'p': 1, 'j': 53, 'a': 4, 'b': 16, 'd (non-minor)': 17},
     'region': '(3) ab >= 1: 1/j <= p + (ab-1)/((a+1)(b+1)); (4) ab = z^2 < 1, a > z: 1/j <= p - a(1-z)^2/(a-z)^2',
     'status': 'REDUCES TO PROVED ATOMS (F1-5, F1-6, F1-dS via an exact group element)',
     'numeric': ['F3 T41 bad17 (3) transported boundary', 'F3 T41 bad17 (4) transported boundary'],
     'docs': ['final_referee_T41_bad17_uniform_20261022 p.2', 'final_referee_T41_bad17_large_product_uniform_20261022 p.1',
              'final_referee_T41_bad54_uniform_20261022 p.1 (SOC, old-H, chain only)',
              'final_referee_T41_negative_Heron_region_20261022 (no universal SO(4) inequality used)']},
    {'id': 'F4-1', 'name': 'inventory56 (1): two ordinary chains centred at f12',
     'generators': {'q': 0, 'j': 2, 'e': 7, 'a': 29, 'b (centre)': 12, 'd (non-minor)': 35},
     'region': 'b >= Bch = d^2 (q+j)(e+a)/(qj(e+a) + ea(q+j)), q, j, e, a > 0',
     'status': 'REDUCES TO PROVED ATOMS (chain (12;0,2;35) + chain (12;7,29;35), exact allocation)',
     'numeric': ['F4 inv56 (1) two-chain boundary b = Bch'],
     'docs': ['final_referee_inventory56_boundary17_uniform_20261023 p.1', 'final_referee_inventory56_boundary17_large_product_20261023 p.1 (1)']},
    {'id': 'F4-2', 'name': 'inventory56 (2): rank-two', 'status': 'see R', 'docs': ['final_referee_inventory56_boundary17_uniform_20261023 p.2']},
    {'id': 'F4-3', 'name': 'inventory56 (3): four-coefficient (7, 0; 5; 35)',
     'generators': {'e': 7, 'q': 0, 'c': 5, 'd (non-minor)': 35}, 'region': 'qe(c - d) >= d^2(e + q + c + d)',
     'status': 'REDUCES TO PROVED ATOMS (four-coefficient family, first pass)',
     'numeric': ['F4 inv56 (3) four-coefficient boundary'],
     'docs': ['final_referee_inventory56_boundary17_uniform_20261023 p.1', 'final_referee_inventory56_boundary17_large_product_20261023 p.1 (2)']},
    {'id': 'F4-nh', 'name': 'inventory56 native harmonic', 'status': 'no capacity used; only minor nonnegativity L_O >= d Phi_35(O)',
     'docs': ['final_referee_inventory56_native_harmonic_20261023 p.5']},
]
for it in INVENTORY:
    it['numeric_min'] = {k: RES['numeric'][k]['min'] for k in it.get('numeric', []) if k in RES['numeric']}
RES['inventory'] = INVENTORY
RES['faces'] = {f.label: {'sources': f.src, 'generators': f.gens, 'gram_rank, cost_rank, combined, relations': f.ranks,
                          'monomial_rows': f.nmons} for f in (F56, F24, F25)}
RES['n_exact_checks'] = len(RES['exact'])
RES['n_exact_failures'] = NFAIL
RES['seconds'] = round(time.time() - T0, 1)
OUTJ = os.environ.get('SO4_PART2A_RAW_JSON')
if OUTJ:
    with open(OUTJ, 'w') as fh:
        json.dump(RES, fh, indent=1, default=str)
print('TOTAL exact checks %d, failures %d, %.1f s' % (len(RES['exact']), NFAIL, time.time() - T0))
