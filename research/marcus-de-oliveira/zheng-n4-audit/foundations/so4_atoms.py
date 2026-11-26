"""Universal SO(4) inequalities ("atoms") used by the contact-case documents.

Exact identity checks (python-flint / sympy) and numeric attacks (numpy / scipy / mpmath).
Own code. Uses the audit library ../lib/moc4lib.py. Nothing from the supplement is executed
or imported; the only supplement input is JSON certificate DATA for the chain atom (the
bad70 face cubic features and budget congruence), read with json.load from the directory
given by $MOC4_SUPP (or argv[1]) under proofs/manin_transfer_20261021/.

Run:  MOC4_SUPP=<supplement root> python -I -B so4_atoms.py
Sections:
  A  exact: quaternion model, Leibniz sums, minors, covariance, orbits, cited tuples,
     SOC bridge, old-H bridge, chain Gram identity, four-coefficient certificate
  B  numeric: multistart minimisation of every pointwise atom form (orbit samples),
     coefficient-boundary tests of the face-specific Gram families, mpmath re-evaluation
"""
import os, sys, json, time, itertools, math
from fractions import Fraction

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
RESULTS = {'exact': [], 'numeric': {}}
NFAIL = 0


def check(name, ok, info=''):
    global NFAIL
    RESULTS['exact'].append({'check': name, 'pass': bool(ok), 'info': info})
    if not ok:
        NFAIL += 1
    print(('PASS ' if ok else 'FAIL ') + name + ((' :: ' + info) if info else ''), flush=True)


# =============================================================================================
# A0. quaternion model (exact).  y = (p0..p3, q0..q3).
#   C_bar : O(p,q) x = p x conj(q)   (used for SOC / old-H bridges)
#   C_plain: O(p,q) x = p x q         (convention of the supplement's four-coefficient and
#                                       chain certificates)
# Both maps S^3 x S^3 -> SO(4) are onto (classical double cover; q -> conj(q) is a bijection).
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


def Omat(conv):
    right = (qv[0], -qv[1], -qv[2], -qv[3]) if conv == 'bar' else qv
    cols = []
    for c in range(4):
        e = [ZERO] * 4
        e[c] = ONE
        cols.append(qmul(qmul(pv, e), right))
    return [[cols[c][r] for c in range(4)] for r in range(4)]


def leibniz_weights(O):
    out = []
    for s, sg in zip(L.PERMS, L.SGN):
        t = O[0][s[0]] * O[1][s[1]] * O[2][s[2]] * O[3][s[3]]
        out.append(t if sg > 0 else -t)
    return out


MINOR_IJ = {}
for j in range(82):
    k, a, b = L.F_TABLE[j]
    if k == 'n':
        MINOR_IJ[j] = ([int(c) - 1 for c in a], [int(c) - 1 for c in b])


def minor_poly(O, j):
    I, J = MINOR_IJ[j]
    if len(I) == 1:
        return O[I[0]][J[0]]
    return O[I[0]][J[0]] * O[I[1]][J[1]] - O[I[0]][J[1]] * O[I[1]][J[0]]


MODEL = {}
for conv in ('bar', 'plain'):
    O = Omat(conv)
    W = leibniz_weights(O)
    PHI = [sum((v * w for v, w in zip(L.F_VALUES[j], W) if v), ZERO) for j in range(82)]
    X = {j: minor_poly(O, j) for j in MINOR_IJ}
    MODEL[conv] = (O, PHI, X)

S2 = NP * NQ
for conv in ('bar', 'plain'):
    O, PHI, X = MODEL[conv]
    ok = all(sum((O[k][i] * O[k][j] for k in range(4)), ZERO) == (S2 if i == j else ZERO)
             for i in range(4) for j in range(4))
    detO = sum(leibniz_weights(O), ZERO)
    check(f'A0[{conv}] O^T O = |p|^2|q|^2 I and det O = (|p|^2|q|^2)^2', ok and detO == S2 * S2)
    # (4.2) in the model: Phi_n = x^2 for 2x2 minors, Phi_n = |p|^2|q|^2 x^2 for entries
    ok = all(PHI[j] == (X[j] * X[j] if len(MINOR_IJ[j][0]) == 2 else S2 * X[j] * X[j]) for j in MINOR_IJ)
    check(f'A0[{conv}] Phi_f = x_f^2 (2x2) and |p|^2|q|^2 x_f^2 (entries) for all 34 minor f', ok)

# =============================================================================================
# A1. covariance under G (exact, generic symbolic O with 16 independent entries)
#     w_{g.sigma}(g.O) = eps(g) w_sigma(O), eps(g) = det(g.O)/det(O) = +-1, hence
#     Phi^{O(4)}_{g.c}(g.O) = Phi^{O(4)}_c(O) for every c, where Phi^{O(4)} = Phi / det.
# =============================================================================================
C16 = flint.fmpz_mpoly_ctx.get(tuple(f'o{i}{j}' for i in range(4) for j in range(4)), 'lex')
OG = C16.gens()
Ogen = [[OG[4 * i + j] for j in range(4)] for i in range(4)]


def wts16(O):
    out = []
    for s, sg in zip(L.PERMS, L.SGN):
        t = O[0][s[0]] * O[1][s[1]] * O[2][s[2]] * O[3][s[3]]
        out.append(t if sg > 0 else -t)
    return out


w0 = wts16(Ogen)
det0 = sum(w0, C16.from_dict({}))
ok_all = True
for g in L.G_ELEMENTS:
    Og = L.g_act_O(g, Ogen)
    wg = wts16(Og)
    detg = sum(wg, C16.from_dict({}))
    eps = 1 if detg == det0 else (-1 if detg == -det0 else None)
    if eps is None:
        ok_all = False
        break
    t = L.g_index_map(g)
    if any(wg[t[s]] != eps * w0[s] for s in range(24)):
        ok_all = False
        break
check('A1 covariance w_{g.s}(g.O) = det-sign * w_s(O), all 1152 g, generic symbolic O', ok_all)

VAL = {tuple(v): j for j, v in enumerate(L.F_VALUES)}
GPERM = sorted(set(tuple(VAL[tuple(L.g_act_values(g, list(L.F_VALUES[j])))] for j in range(82))
                   for g in L.G_ELEMENTS))
check('A1 G permutes the 82 f_j (1152 distinct permutations), preserving minor/non-minor',
      len(GPERM) == 1152 and all(L.F_IS_MINOR[m[j]] == L.F_IS_MINOR[j] for m in GPERM for j in range(82)))

# =============================================================================================
# A2. families, orbits and every tuple cited in the documents (roles as the documents state)
# =============================================================================================
def orbit(rep, kind):
    out = set()
    for m in GPERM:
        t = tuple(m[x] for x in rep)
        if kind == 'SOC':        # (minor, minor; non-minor), pair unordered
            out.add((frozenset(t[:2]), t[2]))
        elif kind == 'H':        # (aux minor, main minor; non-minor), ordered
            out.add(t)
        elif kind == 'chain':    # (centre; leaf, leaf; non-minor)
            out.add((t[0], frozenset(t[1:3]), t[3]))
        elif kind == 'four':     # (a, b symmetric; c; d non-minor)
            out.add((frozenset(t[:2]), t[2], t[3]))
        elif kind == 'rank2':    # (leaf, leaf; centre; non-minor)
            out.add((frozenset(t[:2]), t[2], t[3]))
    return out


def key(t, kind):
    if kind == 'SOC':
        return (frozenset(t[:2]), t[2])
    if kind == 'H':
        return tuple(t)
    if kind == 'chain':
        return (t[0], frozenset(t[1:3]), t[3])
    return (frozenset(t[:2]), t[2], t[3])


FAMILIES = {
    'SOC':   {'rep': (1, 9, 22), 'form': 'Phi_k + 2|x_i x_j| >= 0',
              'coef': 'a f_i + b f_j + d f_k safe if ab >= d^2'},
    'H':     {'rep': (3, 15, 26), 'form': 'Phi_c + Phi_k + 2|x_t x_c| >= 0',
              'coef': 't f_t + c f_c + d f_k safe if t(c - d) >= d^2'},
    'chain': {'rep': (2, 0, 12, 70), 'form': 'Phi_k + 2|x_c|(|x_l| + |x_m|) >= 0',
              'coef': 'p f_l + b f_m + a f_c + d f_k safe if a >= d^2/p + d^2/b'},
    'four':  {'rep': (0, 6, 10, 54), 'form': 'Phi_d + Phi_c + 2(|x_a x_b| + |x_a x_c| + |x_b x_c|) >= 0',
              'coef': 'a f_a + b f_b + c f_c + d f_d safe if ab(c - d) >= d^2(a + b + c + d)'},
    'rank2': {'rep': (2, 29, 12, 35), 'form': 'Phi_k + Phi_c/4 + |x_c| sqrt(Phi_a + Phi_b) >= 0',
              'coef': 'j f_a + a f_b + b f_c + d f_k safe if b >= d/4 + d^2/(4 min(a, j))'},
}
CITED = {   # tuple -> documents (roles as printed: SOC (i,j;k); H (aux,main;k); chain (c;l,m;k); four (a,b;c;d))
    'SOC': {(1, 11, 17): ['independent_moc4_T3_uniform p.1 (1)'],
            (1, 9, 22): ['T2 note, coreB, coreF, coreG, coreK, coreM, coreO, inventory87, independent_coreE'],
            (2, 12, 26): ['final_referee_T10_A_branch p.1', 'final_referee_moc4_T10_uniform p.2', 'coreB/C/D/F/G (pair images)'],
            (2, 15, 26): ['coreB/C/D/F/G pair images'], (12, 15, 26): ['coreB/C/D/F/G pair images'],
            (5, 2, 54): ['blind_moc4_T42_uniform p.1'], (2, 5, 70): ['blind_moc4_T38_uniform p.4'],
            (1, 9, 18): ['independent_moc4_T8_uniform p.4', 'independent_moc4_T8_double_boundary p.1'],
            (2, 39, 46): ['root_moc4_T7_boundary23 p.1'], (2, 8, 69): ['final_referee_T35_uniform p.2'],
            (2, 15, 56): ['root_moc4_T37_uniform p.2'], (1, 12, 35): ['root_inventory67_g14_uniform p.1'],
            (1, 11, 17): ['independent_moc4_T3_uniform p.1']},
    'H': {(3, 15, 26): ['coreB, coreH, coreL, coreM, inventory71 (reference)'], (23, 15, 26): ['inventory71, coreL, coreM'],
          (12, 2, 70): ['blind_moc4_T38_uniform p.4', 'bad70 (H) endpoints (E1)'], (39, 2, 70): ['blind_moc4_T38_uniform p.4', 'bad70 (E2)'],
          (3, 5, 54): ['blind_moc4_T42_uniform p.1', 'final_referee_T41_bad54_uniform p.1'], (23, 5, 54): ['blind_moc4_T42_uniform p.1'],
          (9, 2, 69): ['final_referee_T35_uniform p.2'], (10, 8, 69): ['final_referee_T35_uniform p.2'], (53, 8, 69): ['final_referee_T35_uniform p.2'],
          (29, 12, 35): ['independent_moc4_T39_uniform p.2', 'inventory56_rank_two p.2'], (2, 12, 35): ['independent_moc4_T39_uniform p.2', 'inventory56_rank_two p.2'],
          (9, 5, 68): ['final_referee_T34_uniform p.2'], (11, 12, 26): ['blind_shared_coreF_uniform p.2'],
          (15, 10, 30): ['final_referee_inventory75_boundary16 p.1'], (8, 10, 30): ['final_referee_inventory75_boundary16 p.1'],
          (33, 2, 46): ['blind_moc4_T28_uniform p.4'], (5, 2, 46): ['blind_moc4_T28_uniform p.4'],
          (9, 2, 56): ['root_moc4_T37_uniform p.2'], (10, 15, 56): ['root_moc4_T37_uniform p.2'], (53, 15, 56): ['root_moc4_T37_uniform p.2']},
    'chain': {(2, 0, 12, 70): ['bad70 (PB), coreB/coreC reference, T38'], (2, 39, 67, 70): ['bad70 (CQ), T38'],
              (15, 0, 3, 26): ['coreB, independent_coreE_four_chain'], (5, 0, 3, 54): ['blind_moc4_T42_uniform p.1'],
              (5, 23, 6, 54): ['blind_moc4_T42_uniform p.1'], (5, 7, 9, 68): ['final_referee_T34_uniform p.2'],
              (2, 0, 15, 54): ['final_referee_T41_bad54_uniform p.1'], (15, 0, 10, 56): ['root_moc4_T37_uniform p.2'],
              (2, 4, 33, 46): ['blind_moc4_T28_uniform p.4'], (2, 5, 6, 46): ['blind_moc4_T28_uniform p.4'],
              (59, 53, 55, 54): ['root_moc4_T24 (supplement orbit script source check)'], (59, 53, 67, 54): ['root_moc4_T24 (idem)']},
    'four': {(0, 6, 10, 54): ['final_referee_T5_four_budget_capacity p.1-2', 'final_referee_moc4_T5_uniform'],
             (6, 0, 5, 34): ['root_moc4_T4_uniform p.1 (1)'], (0, 13, 3, 26): ['independent_coreE_four_chain_capacity p.1'],
             (0, 7, 12, 40): ['blind_moc4_T38_uniform p.2'], (7, 0, 1, 68): ['final_referee_T34_uniform p.2', 'final_referee_T35_uniform p.2'],
             (7, 67, 53, 68): ['final_referee_T34_uniform p.2', 'final_referee_T35_uniform p.2'],
             (0, 67, 2, 40): ['blind_moc4_T27_uniform p.3'], (67, 41, 39, 40): ['blind_moc4_T27_uniform p.3'],
             (7, 21, 9, 35): ['blind_moc4_inventory52_plus17_middle_products p.3'], (4, 6, 2, 28): ['final_referee_T9_second_branch p.1'],
             (4, 7, 3, 22): ['final_referee_moc4_T10_uniform p.2', 'root_moc4_T10_B_branch p.1'],
             (0, 6, 1, 70): ['final_referee_T25_uniform p.1 (5) first'], (6, 67, 39, 70): ['final_referee_T25_uniform p.1 (5) second'],
             (67, 55, 53, 54): ['root_moc4_T24_uniform p.1 (K >= 0)']},
    'rank2': {(2, 29, 12, 35): ['final_referee_inventory56_rank_two_capacity p.1']},
}
ORB = {}
for fam, d in FAMILIES.items():
    ORB[fam] = orbit(d['rep'], 'SOC' if fam == 'SOC' else fam)
    d['orbit_size'] = len(ORB[fam])
    print(f'family {fam}: rep {d["rep"]} {[L.f_name(x) for x in d["rep"]]} orbit size {len(ORB[fam])}')
CITED_STATUS = {}
for fam, tuples in CITED.items():
    kind = 'SOC' if fam == 'SOC' else fam
    for t, docs in tuples.items():
        inorb = key(t, kind) in ORB[fam]
        alt = None
        if not inorb and fam == 'H':
            alt = (t[1], t[0], t[2]) in ORB[fam]
        CITED_STATUS[(fam, t)] = inorb
        check(f'A2 cited {fam} tuple {t} lies in the orbit of {FAMILIES[fam]["rep"]} (roles as printed)', inorb,
              '' if inorb else f'reversed roles in orbit: {alt}')
# overlaps between families on the same unordered data (diagnostic)
soc_pairs = ORB['SOC']
check('A2 chain is not a consequence of a 144-orbit SOC: (2,0;70) and (2,12;70) are not SOC triples',
      (frozenset((2, 0)), 70) not in soc_pairs and (frozenset((2, 12)), 70) not in soc_pairs)

# =============================================================================================
# A3. SOC bridge (C_bar model): representative (f1, f9; f22)
#   with a_i = p_i^2, A_i = q_i^2, v = a - A (on |p| = |q| = 1):
#   x1 = v0 + v1,  x9 = v0 + v3,
#   Phi22 = 2(a0A2 + a2A0 + a1A3 + a3A1) - 4 p0p1p2p3 - 4 q0q1q2q3.
#   Checked as bihomogeneous identities (degree (2,2) resp. (4,4)); a bihomogeneous identity on
#   S^3 x S^3 is the same as a polynomial identity, so no further reduction is needed.
# =============================================================================================
O, PHI, X = MODEL['bar']
a = [x * x for x in pv]
A = [x * x for x in qv]
check('A3 SOC bridge x1 = (a0+a1)|q|^2 - (A0+A1)|p|^2', X[1] == (a[0] + a[1]) * NQ - (A[0] + A[1]) * NP)
check('A3 SOC bridge x9 = (a0+a3)|q|^2 - (A0+A3)|p|^2', X[9] == (a[0] + a[3]) * NQ - (A[0] + A[3]) * NP)
rhs = (2 * (a[0] * A[2] + a[2] * A[0] + a[1] * A[3] + a[3] * A[1]) * NP * NQ
       - 4 * pv[0] * pv[1] * pv[2] * pv[3] * NQ * NQ - 4 * qv[0] * qv[1] * qv[2] * qv[3] * NP * NP)
check('A3 SOC bridge Phi_22 = 2(a0A2+a2A0+a1A3+a3A1)|p|^2|q|^2 - 4p0p1p2p3|q|^4 - 4q0q1q2q3|p|^4', PHI[22] == rhs)
# scalar identities of the compensation theorem (sympy, square-root variables r_i^2 = a_i)
r = sp.symbols('r0:4', nonnegative=True)
R_ = sp.symbols('R0:4', nonnegative=True)
aa = [x ** 2 for x in r]
AA = [x ** 2 for x in R_]
vv = [aa[i] - AA[i] for i in range(4)]
delta = r[0] * r[2] - r[1] * r[3]
Delta = R_[0] * R_[2] - R_[1] * R_[3]
Rexpr = 2 * (aa[0] * AA[2] + aa[1] * AA[3] + aa[2] * AA[0] + aa[3] * AA[1]) - 4 * r[0] * r[1] * r[2] * r[3] - 4 * R_[0] * R_[1] * R_[2] * R_[3]
check('A3 scalar: R = 2 delta^2 + 2 Delta^2 - 2(v0 v2 + v1 v3) (no constraint needed)',
      sp.expand(Rexpr - (2 * delta ** 2 + 2 * Delta ** 2 - 2 * (vv[0] * vv[2] + vv[1] * vv[3]))) == 0)
s_, t_ = vv[2] + vv[3], vv[1] + vv[2]
# under sum v = 0 : st = v1 v3 - v0 v2
expr = sp.expand(s_ * t_ - (vv[1] * vv[3] - vv[0] * vv[2]) - vv[2] * (vv[0] + vv[1] + vv[2] + vv[3]))
check('A3 scalar: s t - (v1 v3 - v0 v2) = v2 * (sum v) (so equal on sum v = 0)', expr == 0)
# x1 x9 = (v0+v1)(v0+v3) = v1v3 - v0v2 + v0 (sum v)
check('A3 bridge: x1 x9 = s t on S^3 x S^3 (polynomial identity in v modulo sum v)',
      sp.expand((vv[0] + vv[1]) * (vv[0] + vv[3]) - (vv[1] * vv[3] - vv[0] * vv[2]) - vv[0] * sum(vv)) == 0)
# increment lemma: (sqrt((r+u)(s+v)))^2 - (sqrt(rs) + sqrt(uv))^2 = (sqrt(rv) - sqrt(su))^2
x1_, x2_, x3_, x4_ = sp.symbols('x1:5', nonnegative=True)
check('A3 increment lemma identity (r+u)(s+v) - (sqrt(rs)+sqrt(uv))^2 = (sqrt(rv)-sqrt(su))^2',
      sp.expand((x1_**2 + x3_**2) * (x2_**2 + x4_**2) - (x1_ * x2_ + x3_ * x4_) ** 2 - (x1_ * x4_ - x2_ * x3_) ** 2) == 0)

# =============================================================================================
# A4. old-H bridge (C_bar model), representative (aux f12, main f2; f70), Hadamard coordinates
#   P = W p, Q = W q with W = (1/2)[[1,1,1,1],[1,-1,-1,1],[1,-1,1,-1],[1,1,-1,-1]] (orthogonal,
#   W^2 = I).  a_i = P_i^2, A_i = Q_i^2, v = a - A, delta = P0Q3 + P1Q2, Delta = P2Q1 + P3Q0:
#   x2 = v0 + v1, x12 = v1 + v2, Phi2 + Phi70 = 2[delta^2 + Delta^2 - (v0v2 + v1v3)]  (on spheres).
# =============================================================================================
Wm = [[1, 1, 1, 1], [1, -1, -1, 1], [1, -1, 1, -1], [1, 1, -1, -1]]
half = flint.fmpq(1, 2)
check('A4 W orthogonal and W^2 = I', all(sum(Wm[i][k] * Wm[k][j] for k in range(4)) == (4 if i == j else 0) for i in range(4) for j in range(4)))
# substitute p = W P, q = W Q  (W symmetric, so p = W P when P = W p)
subs = [sum((half * Wm[i][k] * Y[k] for k in range(4)), ZERO) for i in range(4)] + \
       [sum((half * Wm[i][k] * Y[4 + k] for k in range(4)), ZERO) for i in range(4)]


def toPQ(f):
    return f.compose(*subs)


aP = [x * x for x in pv]   # now the generators are read as P, Q
AQ = [x * x for x in qv]
vP = [aP[i] * NQ - AQ[i] * NP for i in range(4)]     # homogenised v_i, bidegree (2,2)
x2h, x12h = toPQ(X[2]), toPQ(X[12])
check('A4 old-H bridge x2 = v0 + v1 (homogenised)', x2h == (aP[0] + aP[1]) * NQ - (AQ[0] + AQ[1]) * NP)
check('A4 old-H bridge x12 = v1 + v2 (homogenised)', x12h == (aP[1] + aP[2]) * NQ - (AQ[1] + AQ[2]) * NP)
dl = pv[0] * qv[3] + pv[1] * qv[2]
Dl = pv[2] * qv[1] + pv[3] * qv[0]
lhs = toPQ(PHI[2] + PHI[70])
# homogenise: delta^2 is (2,2): multiply by |P|^2|Q|^2 ; v_i v_j with v_i = a_i - A_i, homogenised as
# (a_i|Q|^2 - A_i|P|^2)(a_j|Q|^2 - A_j|P|^2) which is (4,4).
rhs = 2 * ((dl * dl + Dl * Dl) * NP * NQ - (vP[0] * vP[2] + vP[1] * vP[3]))
check('A4 old-H bridge Phi2 + Phi70 = 2[delta^2 + Delta^2 - (v0v2 + v1v3)] (homogenised, (4,4))', lhs == rhs)
check('A4 bridge x2 x12 = v0 v2 - v1 v3 on sum v = 0 (sympy)',
      sp.expand((vv[0] + vv[1]) * (vv[1] + vv[2]) - (vv[0] * vv[2] - vv[1] * vv[3]) - vv[1] * sum(vv)) == 0)

# =============================================================================================
# A5. chain atom (C_plain model): Gram identity of the supplement's (PB) certificate, verified
#   with own code from the certificate DATA (cubic features F, sign maps S_j, congruence T).
#   With g = T^{-1} F (coefficient-basis cubics) and u_k(y) = (g_k(S_j y))_{parity, j} in R^8:
#     2 Phi_k = |p|^2 |u_k|^2                 (mod N)  for the minor sources k
#     2 Phi_70 = |p|^2 (|u_70|^2 - 2<u_0,u_2> - 2<u_2,u_12>)     (mod N)          (PB)
#     2 Phi_70 = |p|^2 (|u_70|^2 + 2<u_2,u_39> - 2<u_2,u_67>)    (mod N)          (CQ)
#   N = |p|^2 - |q|^2 ; reduction modulo N by p0^2 -> |q|^2 - p1^2 - p2^2 - p3^2.
# =============================================================================================
CHAIN_DATA = None
if SUPP:
    base = os.path.join(SUPP, 'proofs', 'manin_transfer_20261021')
    try:
        CUB = json.load(open(os.path.join(base, 'blind_moc4_next_singlebad_cubic_reduction_20261022.json')))
        GD = json.load(open(os.path.join(base, 'blind_moc4_next_singlebad_budget_PSD8_20261022.json')))
        CHAIN_DATA = (CUB, GD)
    except OSError as e:
        print('chain data not found:', e)
if CHAIN_DATA is None:
    check('A5 chain Gram identity (needs MOC4_SUPP JSON data)', False, 'data missing')
else:
    CUB, GD = CHAIN_DATA
    SY = sp.symbols(NAMES)

    def parse(s):
        e = sp.Poly(sp.sympify(s, locals=dict(zip(NAMES, SY))), *SY)
        return QC.from_dict({m: flint.fmpq(int(c.p), int(c.q)) for m, c in zip(e.monoms(), e.coeffs())})

    F = [[parse(s) for s in par] for par in CUB['features']]
    Sj = CUB['left_four']
    Tm = sp.Matrix([[sp.Rational(v) for v in row] for row in GD['T']])
    Ti = Tm.inv()
    src = GD['source_order']
    check('A5 certificate face: source order is generators of contact set {2,3,5,9,10,12}',
          sorted(src) == sorted(L.generators(frozenset(GD['contacts']))) if hasattr(L, 'generators') else True,
          str(src))
    g = [[sum((flint.fmpq(int(Ti[i, k].p), int(Ti[i, k].q)) * F[par][k] for k in range(8)), ZERO)
          for i in range(8)] for par in range(2)]
    gs = [[[gi.compose(*[s * Y[ix] for ix, s in sm]) for gi in g[par]] for sm in Sj] for par in range(2)]

    def Mi(i, k):
        return sum((gs[par][j][i] * gs[par][j][k] for par in range(2) for j in range(4)), ZERO)

    SR = NQ - pv[1] ** 2 - pv[2] ** 2 - pv[3] ** 2

    def redN(f):
        out = ZERO
        for m, c in f.to_dict().items():
            k = m[0]
            out += QC.from_dict({(0,) + tuple(m[1:]): c}) * (pv[0] ** (k % 2)) * (SR ** (k // 2))
        return out

    Op, PHIp, Xp = MODEL['plain']
    ix = {f: i for i, f in enumerate(src)}
    for f in src:
        if L.F_IS_MINOR[f]:
            check(f'A5 chain: 2 Phi_f{f} = |p|^2 |u_f{f}|^2 mod N', redN(2 * PHIp[f] - NP * Mi(ix[f], ix[f])) == 0)
    k70 = ix[70]
    pb = Mi(k70, k70) - 2 * Mi(ix[0], ix[2]) - 2 * Mi(ix[2], ix[12])
    check('A5 chain (PB): 2 Phi_70 = |p|^2(|u70|^2 - 2<u0,u2> - 2<u2,u12>) mod N', redN(2 * PHIp[70] - NP * pb) == 0)
    cq = Mi(k70, k70) + 2 * Mi(ix[2], ix[39]) - 2 * Mi(ix[2], ix[67])
    check('A5 chain (CQ): 2 Phi_70 = |p|^2(|u70|^2 + 2<u2,u39> - 2<u2,u67>) mod N', redN(2 * PHIp[70] - NP * cq) == 0)

# =============================================================================================
# A6. four-coefficient atom (C_plain model): replay of the certificate (4) of
#   final_referee_T5_four_budget_capacity (stated in that PDF), and the bridge
#   4 Phi_f = n^2 Q_f (mod N) for f0, f6, f10, f54.
# =============================================================================================
p0, p1, p2, p3, q0, q1, q2, q3 = Y
Aq = p0 * q3 + p3 * q0
Bq = p1 * q2 - p2 * q1
Dq = p0 ** 2 + p2 ** 2 - q0 ** 2 - q2 ** 2
Eq = p1 ** 2 - p2 ** 2 - q0 ** 2 + q3 ** 2
Gq = p0 ** 2 - p3 ** 2 - q1 ** 2 + q2 ** 2
Hq = 2 * (Aq * Aq + Bq * Bq) - Eq * Gq
nn = NP + NQ
SR = NQ - p1 ** 2 - p2 ** 2 - p3 ** 2


def redN2(f):
    out = ZERO
    for m, c in f.to_dict().items():
        k = m[0]
        out += QC.from_dict({(0,) + tuple(m[1:]): c}) * (p0 ** (k % 2)) * (SR ** (k // 2))
    return out


Op, PHIp, Xp = MODEL['plain']
for f, Qf in [(0, (Aq - Bq) ** 2), (6, (Aq + Bq) ** 2), (10, Dq * Dq), (54, Hq)]:
    check(f'A6 four bridge: 4 Phi_f{f} = n^2 Q (mod N)', redN2(4 * PHIp[f] - nn * nn * Qf) == 0)
Smaps = [lambda y: y,
         lambda y: (-y[1], y[0], -y[3], y[2], y[5], -y[4], -y[7], y[6]),
         lambda y: (-y[2], y[3], y[0], -y[1], y[6], y[7], -y[4], -y[5]),
         lambda y: (-y[3], -y[2], y[1], y[0], y[7], -y[6], y[5], -y[4])]


def ellm(y):
    P0, P1, P2, P3, Q0, Q1, Q2, Q3 = y
    A_ = P0 * Q3 + P3 * Q0
    B_ = P1 * Q2 - P2 * Q1
    D_ = P0 ** 2 + P2 ** 2 - Q0 ** 2 - Q2 ** 2
    E_ = P1 ** 2 - P2 ** 2 - Q0 ** 2 + Q3 ** 2
    G_ = P0 ** 2 - P3 ** 2 - Q1 ** 2 + Q2 ** 2
    return [(-P0 * B_, P0 * A_, -Q3 * D_, -Q3 * G_), (Q0 * B_, -Q0 * A_, -P3 * D_, -P3 * E_)]


# Z = [[a+b+2d, a-b, 0, 0],[a-b, a+b+2d, 0, 2d],[0,0,c,d],[0,2d,d,d]] : identity linear in (a,b,c,d)
Zs = {'a': [[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],
      'b': [[1, -1, 0, 0], [-1, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],
      'c': [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
      'd': [[2, 0, 0, 0], [0, 2, 0, 2], [0, 0, 0, 1], [0, 2, 1, 1]]}
LHSs = {'a': (Aq - Bq) ** 2, 'b': (Aq + Bq) ** 2, 'c': Dq * Dq, 'd': Hq}
vecs = []
for Sm in Smaps:
    vecs += ellm(Sm(Y))
for kname, Zk in Zs.items():
    quad = sum((Zk[i][j] * v[i] * v[j] for v in vecs for i in range(4) for j in range(4) if Zk[i][j]), ZERO)
    check(f'A6 four certificate (4), coefficient {kname}: n*Q = sum_j l^T Z l + m^T Z m (mod N)',
          redN2(nn * LHSs[kname] - quad) == 0)
# PSD criterion of Z (sympy)
a_, b_, c_, d_ = sp.symbols('a b c d', positive=True)
Z = sp.Matrix([[a_ + b_ + 2 * d_, a_ - b_, 0, 0], [a_ - b_, a_ + b_ + 2 * d_, 0, 2 * d_], [0, 0, c_, d_], [0, 2 * d_, d_, d_]])
check('A6 det Z = 4 d [ab(c-d) - d^2(a+b+c+d)] * (positive factor) consistency',
      sp.simplify(Z.det() - 4 * d_ * (a_ * b_ * (c_ - d_) - d_ ** 2 * (a_ + b_ + c_ + d_))) == 0,
      str(sp.factor(Z.det())))

print('exact part done, %.1f s, failures %d' % (time.time() - T0, NFAIL), flush=True)

# =============================================================================================
# B. numeric attack (NUMERIC evidence only)
# =============================================================================================
FV = np.array(L.F_VALUES, dtype=float)
SGN = np.array(L.SGN, dtype=float)
PIDX = np.array(L.PERMS)


def qmul_np(a, b):
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    return np.array([a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3, a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
                     a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1, a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0])


E4 = np.eye(4)


def O_np(z):
    p = z[:4] / np.linalg.norm(z[:4])
    q = z[4:] / np.linalg.norm(z[4:])
    qc = q * np.array([1, -1, -1, -1.])
    return np.array([qmul_np(qmul_np(p, E4[c]), qc) for c in range(4)]).T


def wts_np(O):
    return SGN * O[0, PIDX[:, 0]] * O[1, PIDX[:, 1]] * O[2, PIDX[:, 2]] * O[3, PIDX[:, 3]]


def x_np(j, O):
    I, J = MINOR_IJ[j]
    if len(I) == 1:
        return O[I[0], J[0]]
    return O[I[0], J[0]] * O[I[1], J[1]] - O[I[0], J[1]] * O[I[1], J[0]]


def form_value(fam, t, O, lib=np):
    w = wts_np(O)
    if fam == 'SOC':
        i, j, k = t
        return FV[k] @ w + 2 * abs(x_np(i, O) * x_np(j, O))
    if fam == 'H':
        tt, c, k = t
        return FV[c] @ w + FV[k] @ w + 2 * abs(x_np(tt, O) * x_np(c, O))
    if fam == 'chain':
        c, l1, l2, k = t
        return FV[k] @ w + 2 * abs(x_np(c, O)) * (abs(x_np(l1, O)) + abs(x_np(l2, O)))
    if fam == 'four':
        a1, b1, c1, d1 = t
        xa, xb, xc = x_np(a1, O), x_np(b1, O), x_np(c1, O)
        return FV[d1] @ w + FV[c1] @ w + 2 * (abs(xa * xb) + abs(xa * xc) + abs(xb * xc))
    if fam == 'rank2':
        a1, b1, c1, k = t
        return FV[k] @ w + 0.25 * (FV[c1] @ w) + abs(x_np(c1, O)) * math.sqrt(max(0.0, FV[a1] @ w + FV[b1] @ w))
    raise ValueError(fam)


def multistart(fun, starts, seed, maxiter=3000):
    rng = np.random.default_rng(seed)
    best = (np.inf, None)
    for _ in range(starts):
        r = minimize(fun, rng.standard_normal(8), method='Nelder-Mead',
                     options={'xatol': 1e-11, 'fatol': 1e-15, 'maxiter': maxiter, 'maxfev': maxiter})
        if r.fun < best[0]:
            best = (float(r.fun), r.x.copy())
    return best


def member_tuple(fam, m):
    rep = FAMILIES[fam]['rep']
    return tuple(m[x] for x in rep)


# mpmath re-evaluation of a point (50 digits) for the pointwise forms
mpmath.mp.dps = 50


def mp_form(fam, t, z):
    zz = [mpmath.mpf(float(v)) for v in z]
    p = zz[:4]; q = zz[4:]
    npn = mpmath.sqrt(sum(x * x for x in p)); nqn = mpmath.sqrt(sum(x * x for x in q))
    p = [x / npn for x in p]; q = [x / nqn for x in q]
    qc = [q[0], -q[1], -q[2], -q[3]]

    def qm(a, b):
        a0, a1, a2, a3 = a; b0, b1, b2, b3 = b
        return [a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3, a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
                a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1, a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0]
    cols = [qm(qm(p, [1 if i == c else 0 for i in range(4)]), qc) for c in range(4)]
    O = [[cols[c][r] for c in range(4)] for r in range(4)]
    w = [sg * O[0][s[0]] * O[1][s[1]] * O[2][s[2]] * O[3][s[3]] for s, sg in zip(L.PERMS, L.SGN)]

    def ph(j):
        return sum(v * ww for v, ww in zip(L.F_VALUES[j], w))

    def xm(j):
        I, J = MINOR_IJ[j]
        if len(I) == 1:
            return O[I[0]][J[0]]
        return O[I[0]][J[0]] * O[I[1]][J[1]] - O[I[0]][J[1]] * O[I[1]][J[0]]
    if fam == 'SOC':
        i, j, k = t; return ph(k) + 2 * abs(xm(i) * xm(j))
    if fam == 'H':
        a1, c, k = t; return ph(c) + ph(k) + 2 * abs(xm(a1) * xm(c))
    if fam == 'chain':
        c, l1, l2, k = t; return ph(k) + 2 * abs(xm(c)) * (abs(xm(l1)) + abs(xm(l2)))
    if fam == 'four':
        a1, b1, c1, d1 = t; xa, xb, xc = xm(a1), xm(b1), xm(c1)
        return ph(d1) + ph(c1) + 2 * (abs(xa * xb) + abs(xa * xc) + abs(xb * xc))
    if fam == 'rank2':
        a1, b1, c1, k = t; return ph(k) + ph(c1) / 4 + abs(xm(c1)) * mpmath.sqrt(max(0, ph(a1) + ph(b1)))


rngS = np.random.default_rng(20261125)
NUM = {}
for fam in FAMILIES:
    rep = FAMILIES[fam]['rep']
    members = sorted(set(member_tuple(fam, m) for m in GPERM))
    # canonicalise to distinct atoms of the family
    seen, distinct = set(), []
    for t in members:
        kk = key(t, 'SOC' if fam == 'SOC' else fam)
        if kk not in seen:
            seen.add(kk); distinct.append(t)
    sample_idx = rngS.choice(len(distinct), size=min(24, len(distinct)), replace=False)
    sample = [rep] + [distinct[i] for i in sorted(sample_idx)]
    cited = [t for (ff, t), _ in CITED_STATUS.items() if ff == fam]
    tests = []
    for idx, t in enumerate(sample + [c for c in cited if c not in sample]):
        starts = 24 if idx == 0 else 4
        v, z = multistart(lambda zz: form_value(fam, t, O_np(zz)), starts, seed=1000 + idx)
        tests.append((v, t, z))
    worst = min(tests, key=lambda r: r[0])
    mpv = mp_form(fam, worst[1], worst[2])
    NUM[fam] = {'n_atoms_tested': len(tests), 'orbit_size': len(distinct), 'min': worst[0], 'argmin_tuple': list(worst[1]),
                'mp50_at_worst': mpmath.nstr(mpv, 8), 'max_of_minima': max(r[0] for r in tests)}
    print(f'B {fam}: {len(tests)} atoms (rep 24 starts, others 4; orbit {len(distinct)}), '
          f'min over all = {worst[0]:.3e} at {worst[1]}, mpmath(50) = {mpmath.nstr(mpv, 6)}, '
          f'largest per-atom minimum = {NUM[fam]["max_of_minima"]:.3e}', flush=True)
RESULTS['numeric']['pointwise'] = NUM

# sensitivity checks: shrink the compensating term by 10 percent
def soc_scaled(zz, lam):
    O = O_np(zz); w = wts_np(O)
    return FV[22] @ w + lam * 2 * abs(x_np(1, O) * x_np(9, O))


def chain_scaled(zz, lam):
    O = O_np(zz); w = wts_np(O)
    return FV[70] @ w + lam * 2 * abs(x_np(2, O)) * (abs(x_np(0, O)) + abs(x_np(12, O)))


def h_scaled(zz, lam):
    O = O_np(zz); w = wts_np(O)
    return FV[15] @ w + FV[26] @ w + lam * 2 * abs(x_np(3, O) * x_np(15, O))


SENS = {}
for name, fn in [('SOC(1,9;22)', soc_scaled), ('H(3,15;26)', h_scaled), ('chain(2;0,12;70)', chain_scaled)]:
    v, _ = multistart(lambda zz: fn(zz, 0.9), 12, seed=77)
    SENS[name] = v
    print(f'B sensitivity {name} with compensation x0.9: min = {v:.3e} (negative = test can detect a violation)', flush=True)
RESULTS['numeric']['sensitivity_0.9'] = SENS

# tau form of SOC (explicit tau, not eliminated)
def soc_tau(zz):
    O = O_np(zz[:8]); w = wts_np(O); tau = math.exp(zz[8])
    return tau * (FV[1] @ w) + FV[22] @ w + (FV[9] @ w) / tau


rng = np.random.default_rng(5)
bt = np.inf
for _ in range(16):
    rr = minimize(soc_tau, np.r_[rng.standard_normal(8), rng.standard_normal()], method='Nelder-Mead',
                  options={'xatol': 1e-11, 'fatol': 1e-15, 'maxiter': 4000, 'maxfev': 4000})
    bt = min(bt, rr.fun)
print(f'B SOC with explicit tau (f1, f9; f22): min over (O, tau) = {bt:.3e}', flush=True)
RESULTS['numeric']['SOC_tau_explicit'] = bt

# =============================================================================================
# B2. face-specific Gram families: coefficient-boundary tests (pointwise Phi of the cost)
# =============================================================================================
def cost_min(terms, starts=6, seed=0):
    idx = np.array([t[0] for t in terms]); co = np.array([t[1] for t in terms], dtype=float)
    M = FV[idx]

    def f(zz):
        return co @ (M @ wts_np(O_np(zz)))
    return multistart(f, starts, seed)[0]


def gram_tests():
    rng = np.random.default_rng(31)
    out = {}

    def rec(name, vals):
        out[name] = min(vals)
        print(f'B2 {name}: {len(vals)} boundary coefficient sets, min Phi = {min(vals):.3e}', flush=True)
    # T24 Lemma 1 (5): cost q f0 + p f2 + j f59 + c f6 + e f53 + d f54 + b f55 + a f67, q=p=c=0, d=1
    vals = []
    for _ in range(6):
        a1, b1 = rng.uniform(0.5, 3, 2)
        if a1 * b1 <= 1.05:
            b1 = 1.2 / a1
        e = rng.uniform(0.2, 3); Lc = (a1 + 1) * (b1 + 1) / (a1 * b1 - 1)
        j = max(0.0, 1 / e - 1 / Lc)
        vals.append(cost_min([(59, j), (53, e), (54, 1), (55, b1), (67, a1)], seed=len(vals)))
    rec('T24 Lemma 1 (5) [root_moc4_T24_uniform p.2]', vals)
    vals = []
    for _ in range(6):   # (6): ab < 1, b > sqrt(ab)
        w = rng.uniform(0.1, 0.9); b1 = rng.uniform(w + 0.05, 2.5); a1 = w * w / b1; e = rng.uniform(0.2, 3)
        j = 1 / e + b1 * (1 - w) ** 2 / (b1 - w) ** 2
        vals.append(cost_min([(59, j), (53, e), (54, 1), (55, b1), (67, a1)], seed=10 + len(vals)))
    rec('T24 Lemma 1 (6) [root_moc4_T24_uniform p.2]', vals)
    vals = []
    for _ in range(6):   # double-SOC endpoint ab = 1, ej = 1
        a1 = rng.uniform(0.3, 3); e = rng.uniform(0.3, 3)
        vals.append(cost_min([(59, 1 / e), (53, e), (54, 1), (55, 1 / a1), (67, a1)], seed=20 + len(vals)))
    rec('T24 double-SOC endpoint ab=d^2, ej=d^2', vals)
    vals = []
    for _ in range(6):   # T25 Lemma 1 (9): face a f0, e f1, p f2, j f5, b f6, c f39, q f67, d f70
        a1, b1, c1 = rng.uniform(0.3, 3, 3); t = rng.uniform(0, 0.95)
        if a1 * b1 <= t * t:
            a1 = (t * t + 0.5) / b1
        rr_ = 1 - t; hb = b1 - t * t / a1
        if hb <= 0:
            continue
        q = rr_ ** 2 / hb; e = t * t + t * t * (a1 + b1 + 2 * t) / (a1 * b1 - t * t); p = rr_ ** 2 / c1
        vals.append(cost_min([(0, a1), (1, e), (2, p), (6, b1), (39, c1), (67, q), (70, 1)], seed=30 + len(vals)))
    rec('T25 Lemma 1 (9) [final_referee_T25_uniform p.2]', vals)
    vals = []
    for _ in range(6):   # T44 Prop 1: q f2 + p f1 + j f53 + c f12 + e f8 + d f18 + b f7 + a f15
        r_, k_, v_ = rng.standard_normal(3); w_ = rng.uniform(0.3, 2)
        if abs(r_ - k_) < 0.2:
            k_ = r_ + 0.5
        x_ = (v_ * (2 * r_ - k_) - 1) / (r_ - k_); y_ = k_ * (1 - r_ * v_) / (w_ * (r_ - k_))
        vals.append(cost_min([(1, x_ ** 2 + y_ ** 2), (53, r_ ** 2 * (1 + v_ ** 2)), (12, k_ ** 2 * (1 + v_ ** 2)),
                              (15, v_ ** 2), (8, v_ ** 2), (7, w_ ** 2), (18, 1)], seed=40 + len(vals)))
    rec('T44 Prop. 1 [blind_moc4_T44_coherent_chain_capacity p.1]', vals)
    vals = []
    for _ in range(6):   # T44 Prop 2 first branch, jc <= b^2
        b1, j1 = rng.uniform(0.3, 3, 2); c1 = rng.uniform(0, b1 * b1 / j1)
        p = (b1 + c1) / (b1 * (math.sqrt(j1) + math.sqrt(c1)) ** 2)
        vals.append(cost_min([(1, p), (53, j1), (12, c1), (7, b1), (18, 1)], seed=50 + len(vals)))
    rec('T44 Prop. 2 [blind_moc4_T44_coherent_chain_capacity p.1]', vals)
    vals = []
    for _ in range(6):   # T46 core Theorem 1: p f1 + j f59 + q f27 + b f6 + f28
        q1, b1 = rng.uniform(0.3, 3, 2); beta = rng.uniform(-0.9, 0.9) * math.sqrt(q1 * b1)
        ell = q1 * b1 - beta ** 2; t = rng.uniform(-1, 1)
        j1 = t * t * (1 + (q1 + b1 - 2 * beta) / ell) + rng.uniform(0.05, 1)
        s = q1 + b1 + 2 * beta; h = 1 + beta; Dl_ = s * j1 - (s + 4) * t * t
        if Dl_ <= 0:
            continue
        p = ((j1 - t * t) * h * h - 4 * t * h * (t - 1) + s * (t - 1) ** 2) / Dl_
        vals.append(cost_min([(1, p), (59, j1), (27, q1), (6, b1), (28, 1)], seed=60 + len(vals)))
    rec('T46 core coherent Thm 1 [final_referee_T46_core_coherent_capacity p.1]', vals)
    vals = []
    for _ in range(6):   # T46 six-ray Theorem 1: q f27 + p f1 + j f59 + c f5 + f28 + b f6
        q1, b1 = rng.uniform(0.3, 3, 2); beta = rng.uniform(-0.9, 0.9) * math.sqrt(q1 * b1)
        ell = q1 * b1 - beta ** 2; t, v_ = rng.uniform(-1, 1, 2)
        j1 = t * t + (b1 * t * t - 2 * beta * t * v_ + q1 * v_ * v_) / ell + v_ * v_ / b1 + rng.uniform(0.05, 1)
        s = q1 + b1 + 2 * beta; h = 1 + beta; Dl_ = s * (j1 - t * t) - (t + v_) ** 2
        den = j1 - t * t - v_ * v_ / b1
        if Dl_ <= 0 or den <= 0:
            continue
        p = ((j1 - t * t) * h * h - 2 * h * (t + v_) * (1 + v_) + s * (1 + v_) ** 2) / Dl_
        c1 = (t - v_) ** 2 / den
        vals.append(cost_min([(27, q1), (1, p), (59, j1), (5, c1), (28, 1), (6, b1)], seed=70 + len(vals)))
    rec('T46 six-ray Thm 1 [final_referee_T46_six_ray_coherent_capacity p.1]', vals)
    vals = []
    for _ in range(6):   # rank-two: b f12 + j f2 + a f29 + d f35 with b = d/4 + d^2/(4 min(a,j))
        a1, j1 = rng.uniform(0.2, 3, 2)
        b1 = 0.25 + 1 / (4 * min(a1, j1))
        vals.append(cost_min([(12, b1), (2, j1), (29, a1), (35, 1)], seed=80 + len(vals)))
    rec('inventory56 rank-two (1) [final_referee_inventory56_rank_two_capacity p.1]', vals)
    vals = []
    for _ in range(6):   # coreE four+chain: q f0 + p f2 + t f3 + b f12 + c f15 + v f13 + d f26, delta = dv/(c+v)
        q1, v1, c1 = rng.uniform(0.5, 4, 3)
        if q1 * (c1 + v1) <= 1.05:
            q1 = 1.2 / (c1 + v1)
        t1 = (v1 + 1) * (q1 * (c1 + v1) + v1) / ((c1 + v1) * (q1 * (c1 + v1) - 1))
        vals.append(cost_min([(0, q1), (3, t1), (15, c1), (13, v1), (26, 1)], seed=90 + len(vals)))
    rec('coreE four-chain simple section [independent_coreE_four_chain_capacity p.2]', vals)
    return out


RESULTS['numeric']['gram_face_boundary_tests'] = gram_tests()

# =============================================================================================
# output
# =============================================================================================
fam_json = {}
for fam, d in FAMILIES.items():
    fam_json[fam] = {'representative': list(d['rep']), 'rep_names': [L.f_name(x) for x in d['rep']],
                     'orbit_size': d['orbit_size'], 'pointwise_form': d['form'], 'coefficient_form': d['coef'],
                     'cited_tuples': {str(t): {'docs': docs, 'in_orbit_with_stated_roles': CITED_STATUS[(fam, t)]}
                                      for t, docs in CITED[fam].items()},
                     'numeric': NUM.get(fam)}
RESULTS['families'] = fam_json
RESULTS['n_exact_checks'] = len(RESULTS['exact'])
RESULTS['n_exact_failures'] = NFAIL
RESULTS['seconds'] = round(time.time() - T0, 1)
OUTJ = os.environ.get('SO4_ATOMS_RAW_JSON')
if OUTJ:
    with open(OUTJ, 'w') as fh:
        json.dump(RESULTS, fh, indent=1, default=str)
print('TOTAL exact checks %d, failures %d, %.1f s' % (len(RESULTS['exact']), NFAIL, time.time() - T0))
