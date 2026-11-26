"""Gap check for blind_moc4_T44_coherent_chain_capacity_20261022, Proposition 1 Gram step.

Own code. Builds the Leibniz sums Phi_f(O) of the eight sources from the paper's
definitions (via moc4lib), parametrizes SO(4) by unit quaternions, and checks EXACTLY
that a concrete set of eight vector-valued cubic node functions m_K(p,q) realizes the
Gram identity of Proposition 1 (and the b = 0 Gram of Proposition 2).

The candidate node functions are built from the supplement's JSON data
(blind_moc4_T44_cubic_reduction_20261022.json: 'features', 'left_four';
blind_moc4_T44_bad18_budget_PSD8_20261022.json: 'T'), loaded with json.load and
parsed by a small parser below (no eval, nothing of theirs is executed or imported).
They are used only as a candidate certificate; every identity is re-derived here.

Usage: MOC4_SUPP=<supplement root> python -I -B this.py   (or pass the root as argv[1])
"""
import os, sys, re, json
from fractions import Fraction
import flint
import sympy as S

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lib'))
import moc4lib as L

SUPP = sys.argv[1] if len(sys.argv) > 1 else os.environ['MOC4_SUPP']
DD = os.path.join(SUPP, 'proofs', 'manin_transfer_20261021')
FAIL = 0
def check(name, ok):
    global FAIL
    print(('PASS ' if ok else 'FAIL ') + name)
    if not ok: FAIL += 1

# ---- polynomial ring: lex with q0 first so that q0^2 leads (q0^2+q1^2+q2^2+q3^2-np) ----
NAMES = ('q0', 'q1', 'q2', 'q3', 'p0', 'p1', 'p2', 'p3')
ctx = flint.fmpq_mpoly_ctx.get(NAMES, 'lex')
q0, q1, q2, q3, p0, p1, p2, p3 = ctx.gens()
P = [p0, p1, p2, p3]; Q = [q0, q1, q2, q3]
VAR = {'p0': p0, 'p1': p1, 'p2': p2, 'p3': p3, 'q0': q0, 'q1': q1, 'q2': q2, 'q3': q3}
ORDER8 = [p0, p1, p2, p3, q0, q1, q2, q3]    # the JSON's variable order (p0..p3, q0..q3)
NP = sum(x * x for x in P); NQ = sum(x * x for x in Q)
REL = NQ - NP                                   # leading term q0^2 in lex

def red(f):
    return f % REL

def parse(s):
    """Integer polynomial string in p0..p3,q0..q3 -> flint poly (tokenizer, no eval)."""
    s = s.replace(' ', '')
    if s[0] not in '+-': s = '+' + s
    out = ctx.from_dict({}) if hasattr(ctx, 'from_dict') else 0 * p0
    for t in re.findall(r'[+-][^+-]+', s):
        sign = -1 if t[0] == '-' else 1
        body = t[1:]
        toks = re.findall(r'(\d+|[pq]\d)(\*\*\d+)?', body)
        assert '*'.join(a + b for a, b in toks) == body, body
        term = flint.fmpq(sign)
        for a, b in toks:
            e = int(b[2:]) if b else 1
            term = term * ((VAR[a] if a[0] in 'pq' else flint.fmpq(int(a))) ** e)
        out = out + term
    return out

def qmul(a, b):
    a0, a1, a2, a3 = a; b0, b1, b2, b3 = b
    return [a0*b0 - a1*b1 - a2*b2 - a3*b3, a0*b1 + a1*b0 + a2*b3 - a3*b2,
            a0*b2 - a1*b3 + a2*b0 + a3*b1, a0*b3 + a1*b2 - a2*b1 + a3*b0]

# O(p,q): x -> p x q, column j = p e_j q; orthogonal with det +1 times |p|^2|q|^2-scaling
one, zero = flint.fmpq(1) + 0 * p0, 0 * p0
E = [[one if i == j else zero for i in range(4)] for j in range(4)]
cols = [qmul(qmul(P, E[j]), Q) for j in range(4)]
O = [[cols[j][i] for j in range(4)] for i in range(4)]

# sanity: O^T O = np*nq*I
ok = True
for i in range(4):
    for j in range(4):
        s = sum(O[k][i] * O[k][j] for k in range(4))
        ok &= (s - (NP * NQ if i == j else 0)).is_zero()
check('O(p,q)^T O(p,q) = |p|^2|q|^2 I', ok)
check('det O(p,q) = (|p|^2|q|^2)^2 (so O/|p||q| is in SO(4))',
      (L.det_leibniz4(O) - (NP * NQ) ** 2).is_zero())

W = L.leibniz_weights(O)                      # sgn(sigma) prod_i O[i][sigma(i)]
def Phi_f(j):
    return sum((W[s] * L.F_VALUES[j][s] for s in range(24) if L.F_VALUES[j][s] != 0), zero)

SRC = {'J': 53, 'C': 12, 'E': 8, 'A': 15, 'B': 7, 'Q': 2, 'P': 1, 'D': 18}
for K, j in SRC.items():
    print('  node', K, '-> f%d' % j, L.f_name(j), 'minor' if L.F_IS_MINOR[j] else 'NON-minor')
check('all sources except f18 are minor vectors (Phi >= 0 by (4.2))',
      all(L.F_IS_MINOR[j] for K, j in SRC.items() if K != 'D') and not L.F_IS_MINOR[18])
F = {K: Phi_f(j) for K, j in SRC.items()}

# ---- candidate node functions from the JSON data ----
Dj = json.load(open(os.path.join(DD, 'blind_moc4_T44_cubic_reduction_20261022.json')))
Mj = json.load(open(os.path.join(DD, 'blind_moc4_T44_bad18_budget_PSD8_20261022.json')))
check('JSON source order is q,p,j,c,e,d,b,a = f2,f1,f53,f12,f8,f18,f7,f15',
      Mj['source_order'] == [2, 1, 53, 12, 8, 18, 7, 15])
feats = [[parse(s) for s in blk] for blk in Dj['features']]
Tm = S.Matrix([[S.Rational(x) for x in r] for r in Mj['T']]); Ti = Tm.inv()
copies = []
for blk in feats:
    for g in Dj['left_four']:
        sub = [g[i][1] * ORDER8[g[i][0]] for i in range(8)]   # new var i := sign * old var
        # compose takes images of ctx generators in ctx order (q0..q3,p0..p3)
        args = sub[4:8] + sub[0:4]        # images in ctx generator order (q0..q3, p0..p3)
        copies.append([f.compose(*args) for f in blk])
NODES = 'J C E A B Q P D'.split()
m = {}
for K_i, K in enumerate(NODES):
    m[K] = [sum((flint.fmpq(int(Ti[K_i, l].p), int(Ti[K_i, l].q)) * c[l] for l in range(8) if Ti[K_i, l] != 0), zero)
            for c in copies]
def Mx(K, Lb):
    return sum((a * b for a, b in zip(m[K], m[Lb])), zero)

# Relation "X <-> Y" means  np * X - 2 * Y  == 0  mod (nq - np), i.e.  X = 2 Y / np on the
# variety np = nq, i.e. X(p,q) = 2 np^3 * Y(O/np) when np = nq.
def rel(name, Xpoly, Ypoly):
    check(name, red(NP * Xpoly - 2 * Ypoly).is_zero())

Mm = {(a, b): Mx(a, b) for a in NODES for b in NODES if NODES.index(a) <= NODES.index(b)}
def MM(a, b):
    return Mm[(a, b)] if (a, b) in Mm else Mm[(b, a)]
rel('R1  M_BB <-> Phi_B', MM('B', 'B'), F['B'])
rel('R2  M_CC <-> Phi_C', MM('C', 'C'), F['C'])
rel('R3  M_JJ <-> Phi_J', MM('J', 'J'), F['J'])
rel('R4  M_PP <-> Phi_P', MM('P', 'P'), F['P'])
rel('R5  M_JC <-> 0', MM('J', 'C'), zero)
rel('R6  2M_JP - M_DD <-> -Phi_18', 2 * MM('J', 'P') - MM('D', 'D'), -F['D'])
rel('R7  M_CP - M_BP + M_JP <-> 0', MM('C', 'P') - MM('B', 'P') + MM('J', 'P'), zero)
rel('R8  M_JE + M_JA + M_DD <-> Phi_18', MM('J', 'E') + MM('J', 'A') + MM('D', 'D'), F['D'])
rel('R9  2(M_CE + M_CA + M_BP) - M_DD <-> -Phi_18',
    2 * (MM('C', 'E') + MM('C', 'A') + MM('B', 'P')) - MM('D', 'D'), -F['D'])
rel('R10 M_EP + M_AP - M_PD <-> 0', MM('E', 'P') + MM('A', 'P') - MM('P', 'D'), zero)
rel('R11 |m_E + m_A - m_D|^2 <-> Phi_A + Phi_E',
    MM('E', 'E') + MM('A', 'A') + MM('D', 'D') + 2 * MM('E', 'A') - 2 * MM('E', 'D') - 2 * MM('A', 'D'),
    F['A'] + F['E'])

# ---- symbolic: the relations R1..R11 imply the Proposition 1 identity for all r,k,v,w ----
r, k, v, w = S.symbols('r k v w')
x = (v * (2 * r - k) - 1) / (r - k); y = k * (1 - r * v) / (w * (r - k))
V = S.Matrix([[r * v, r, 0], [k * v, k, 0], [0, v, 0], [0, v, 0], [0, 0, w], [0, 0, 0], [0, x, y], [1, -v, 0]])
G = V * V.T
Msym = {}
for i, a in enumerate(NODES):
    for j2, b in enumerate(NODES):
        if i <= j2: Msym[(a, b)] = S.Symbol('M_' + a + b)
def MS(a, b): return Msym[(a, b)] if (a, b) in Msym else Msym[(b, a)]
Ph = {K: S.Symbol('Phi_' + K) for K in NODES}
form = sum(G[i, j] * MS(a, b) for i, a in enumerate(NODES) for j, b in enumerate(NODES))
target = (r**2 * (1 + v**2) * Ph['J'] + k**2 * (1 + v**2) * Ph['C'] + v**2 * Ph['E'] + v**2 * Ph['A']
          + w**2 * Ph['B'] + (x**2 + y**2) * Ph['P'] + Ph['D'])
# impose R1..R11 (solved for M_BB, M_CC, M_JJ, M_PP, M_JC, M_DD, M_CP, M_JE, M_CE, M_EP, M_EE)
subs = {MS('B', 'B'): Ph['B'], MS('C', 'C'): Ph['C'], MS('J', 'J'): Ph['J'], MS('P', 'P'): Ph['P'],
        MS('J', 'C'): 0}
subs[MS('D', 'D')] = 2 * MS('J', 'P') + Ph['D']
subs[MS('C', 'P')] = MS('B', 'P') - MS('J', 'P')
subs[MS('J', 'E')] = Ph['D'] - subs[MS('D', 'D')] - MS('J', 'A')
subs[MS('C', 'E')] = (subs[MS('D', 'D')] - Ph['D']) / 2 - MS('C', 'A') - MS('B', 'P')
subs[MS('E', 'P')] = MS('P', 'D') - MS('A', 'P')
subs[MS('E', 'E')] = (Ph['A'] + Ph['E'] - MS('A', 'A') - subs[MS('D', 'D')] - 2 * MS('E', 'A')
                      + 2 * MS('E', 'D') + 2 * MS('A', 'D'))
diff = S.together((form - target).subs(subs))
check('symbolic: R1..R11 imply  sum_KL G_KL M_KL = cost(G)  for all r,k,v (r!=k), w!=0',
      S.expand(S.numer(diff)) == 0)
# the diagonal stated in the document
diag = [S.simplify(G[i, i]) for i in range(8)]
check('Gram diagonal = (r^2(1+v^2), k^2(1+v^2), v^2, v^2, w^2, 0, x^2+y^2, 1+v^2)',
      all(S.simplify(a - b) == 0 for a, b in zip(diag, [r**2*(1+v**2), k**2*(1+v**2), v**2, v**2, w**2, 0, x**2+y**2, 1+v**2])))
# b = 0 Gram of Proposition 2: J=(0,s,0), P=(0,-1/s,0), D=(1,0,0), s = sqrt(j)
s_ = S.Symbol('s')
V0 = S.zeros(8, 3); V0[0, 1] = s_; V0[6, 1] = -1 / s_; V0[7, 0] = 1
G0 = V0 * V0.T
form0 = sum(G0[i, j] * MS(a, b) for i, a in enumerate(NODES) for j, b in enumerate(NODES))
t0 = s_**2 * Ph['J'] + Ph['P'] / s_**2 + Ph['D']
check('symbolic: b=0 Gram gives  j f53 + (1/j) f1 + f18  under R3, R4, R6',
      S.expand(S.numer(S.together((form0 - t0).subs(subs)))) == 0)

# ---- direct exact spot checks of the full identity at rational parameters (no use of R*) ----
def direct(rv, kv, vv, wv):
    xv = (vv * (2 * rv - kv) - 1) / (rv - kv); yv = kv * (1 - rv * vv) / (wv * (rv - kv))
    Vn = [[rv*vv, rv, 0], [kv*vv, kv, 0], [0, vv, 0], [0, vv, 0], [0, 0, wv], [0, 0, 0], [0, xv, yv], [1, -vv, 0]]
    fq = lambda z: flint.fmpq(z.numerator, z.denominator)
    sos = zero
    for c in range(len(copies)):
        for col in range(3):
            lin = sum((fq(Fraction(Vn[i][col])) * m[K][c] for i, K in enumerate(NODES) if Vn[i][col] != 0), zero)
            sos = sos + lin * lin
    coeff = {'J': rv**2*(1+vv**2), 'C': kv**2*(1+vv**2), 'E': vv**2, 'A': vv**2, 'B': wv**2,
             'P': xv**2 + yv**2, 'D': Fraction(1)}
    cost = sum((fq(Fraction(cf)) * F[K] for K, cf in coeff.items()), zero)
    return red(NP * sos - 2 * cost).is_zero()
for prm in [(Fraction(1), Fraction(-1), Fraction(0), Fraction(1)),
            (Fraction(3, 2), Fraction(-2, 5), Fraction(1, 3), Fraction(7, 4)),
            (Fraction(-2), Fraction(5, 3), Fraction(-3, 7), Fraction(1, 2))]:
    check('direct: np*SOS == 2*Phi_cost mod (nq-np) at (r,k,v,w)=%s' % (tuple(str(z) for z in prm),), direct(*prm))

print('RESULT', 'ALL EXACT CHECKS PASS' if FAIL == 0 else '%d FAILURES' % FAIL)
sys.exit(1 if FAIL else 0)
