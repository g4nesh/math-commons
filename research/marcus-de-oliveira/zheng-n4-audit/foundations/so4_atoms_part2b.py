"""Part 2b of the SO(4) atom audit: allocation lemmas (G5) and the face Gram statements of
historical T44, T45, T46, T47 (G1-G4).

Own code.  Uses ../lib/moc4lib.py.  Nothing from the supplement is executed or imported.
The only supplement input is JSON DATA (face cubic features, quaternion sign maps and the
diagonal-coefficient congruence T of the T44, T46, T47 faces), read with json.load from
$MOC4_SUPP/proofs/manin_transfer_20261021/ (or argv[1]); every identity built from it is
verified exactly here, so the conclusions do not depend on where those polynomials came from.

Run:  MOC4_SUPP=<supplement root> python -I -B so4_atoms_part2b.py [sections]
sections (default all): G5 GRAM G1 G2 G3 G4 NUM
Exact checks print PASS/FAIL; numeric attacks print NUMERIC lines (evidence only).
"""
import os, sys, json, time, math, itertools
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
SUPP = os.environ.get('MOC4_SUPP')
ARGS = [a for a in sys.argv[1:]]
if ARGS and os.path.isdir(ARGS[0]):
    SUPP = ARGS.pop(0)
SECTIONS = set(ARGS) or {'G5', 'GRAM', 'G1', 'G2', 'G3', 'G4', 'NUM'}
BASE = os.path.join(SUPP, 'proofs', 'manin_transfer_20261021') if SUPP else None
OUT = {'exact': [], 'statements': {}, 'numeric': {}}
NFAIL = 0


def check(name, ok, info=''):
    global NFAIL
    OUT['exact'].append({'check': name, 'pass': bool(ok), 'info': str(info)[:300]})
    if not ok:
        NFAIL += 1
    print(('PASS ' if ok else 'FAIL ') + name + ((' :: ' + str(info)) if info else ''), flush=True)
    return ok


def stmt(key, **kw):
    OUT['statements'].setdefault(key, {}).update(kw)


# =============================================================================================
# 0. group action on the 82 facets and on contact indices; atom families (EXACT)
# =============================================================================================
VAL = {tuple(v): j for j, v in enumerate(L.F_VALUES)}
GACT = []   # (g, facet map m, contact map t)
for g in L.G_ELEMENTS:
    m = tuple(VAL[tuple(L.g_act_values(g, list(L.F_VALUES[j])))] for j in range(82))
    GACT.append((g, m, L.g_index_map(g)))
FAM_REP = {'SOC': (1, 9, 22), 'H': (3, 15, 26), 'chain': (2, 0, 12, 70), 'four': (0, 6, 10, 54)}


def fam_key(fam, t):
    if fam == 'SOC':          # (i, j; k), pair unordered
        return (frozenset(t[:2]), t[2])
    if fam == 'H':            # (aux, main; k), ordered
        return tuple(t)
    if fam == 'chain':        # (centre; leaf, leaf; k)
        return (t[0], frozenset(t[1:3]), t[3])
    if fam == 'four':         # (a, b; c; d), a, b unordered
        return (frozenset(t[:2]), t[2], t[3])
    raise ValueError(fam)


ORBIT = {fam: {fam_key(fam, tuple(m[x] for x in rep)) for _, m, _ in GACT} for fam, rep in FAM_REP.items()}


def atom_ok(fam, t, where):
    ok = fam_key(fam, t) in ORBIT[fam]
    rev = ''
    if not ok and fam == 'H':
        rev = 'reversed roles in orbit: %s' % (fam_key('H', (t[1], t[0], t[2])) in ORBIT['H'])
    check(f'atom {fam}{tuple(t)} in proved orbit with stated roles [{where}]', ok, rev)
    return ok


def find_action(src, dst, contact=None):
    """all g with m[src[i]] = dst[i]; optionally also g.contact == contact (as sets)."""
    out = []
    for g, m, t in GACT:
        if all(m[a] == b for a, b in zip(src, dst)):
            if contact is None or frozenset(t[s] for s in contact[0]) == frozenset(contact[1]):
                out.append(g)
    return out


def literal_action_matches(row, col, tr, src, dst):
    """Does some reading of the documents' literal (row, column, transpose) action map src->dst?
    Readings tried: g = (rho, rho2, eps) with rho, rho2 in {row, row^-1, col, col^-1} pairs."""
    cands = []
    R, C = tuple(row), tuple(col)
    for a in (R, L.inverse(R), C, L.inverse(C)):
        for b in (R, L.inverse(R), C, L.inverse(C)):
            cands.append((a, b, int(tr)))
    hits = []
    gm = {gg: mm for gg, mm, _ in GACT}
    for g in cands:
        mm = gm[g]
        if all(mm[a] == b for a, b in zip(src, dst)):
            hits.append(g)
    return hits


def fvals(coeffs):
    """24 values of sum_j w_j f_j (w given as dict facet -> number/sympy)."""
    out = [0] * 24
    for j, w in coeffs.items():
        for s in range(24):
            if L.F_VALUES[j][s]:
                out[s] += w * L.F_VALUES[j][s]
    return out


def ident(lhs, rhs):
    """exact identity sum lhs f = sum rhs f in V_R (P is injective, so compare 24 values)."""
    return fvals(lhs) == fvals(rhs)


# =============================================================================================
# numeric tools (NUMERIC evidence only)
# =============================================================================================
FV = np.array(L.F_VALUES, dtype=float)
PIDX = np.array(L.PERMS)
SG = np.array(L.SGN, dtype=float)


def qmul_np(a, b):
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    return np.array([a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3, a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
                     a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1, a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0])


def O_batch(Z):
    """Z: (N, 8) -> (N, 4, 4) rotation matrices x -> p x conj(q)."""
    p = Z[:, :4] / np.linalg.norm(Z[:, :4], axis=1, keepdims=True)
    q = Z[:, 4:] / np.linalg.norm(Z[:, 4:], axis=1, keepdims=True)
    qc = q * np.array([1, -1, -1, -1.])
    O = np.empty((Z.shape[0], 4, 4))
    for c in range(4):
        e = np.zeros((4, Z.shape[0]))
        e[c] = 1
        col = qmul_np(qmul_np(p.T, e), qc.T)
        O[:, :, c] = col.T
    return O


def wts_batch(O):
    pr = O[:, 0, PIDX[:, 0]] * O[:, 1, PIDX[:, 1]] * O[:, 2, PIDX[:, 2]] * O[:, 3, PIDX[:, 3]]
    return pr * SG


def min_cost(coeffs, seed=0, nsample=6000, npolish=6):
    """NUMERIC: min over SO(4) of Phi of sum_j w_j f_j. Random sample + BFGS polish."""
    v = np.zeros(24)
    for j, w in coeffs.items():
        v += float(w) * FV[j]
    rng = np.random.default_rng(seed)
    Z = rng.standard_normal((nsample, 8))
    vals = wts_batch(O_batch(Z)) @ v
    idx = np.argsort(vals)[:npolish]
    f = lambda z: float(wts_batch(O_batch(z[None, :]))[0] @ v)
    best = (np.inf, None)
    for i in idx:
        r = minimize(f, Z[i], method='BFGS', options={'gtol': 1e-12, 'maxiter': 400})
        if r.fun < best[0]:
            best = (float(r.fun), r.x.copy())
    return best


mpmath.mp.dps = 40


def mp_cost(coeffs, z):
    zz = [mpmath.mpf(float(x)) for x in z]
    p = zz[:4]; q = zz[4:]
    n1 = mpmath.sqrt(sum(x * x for x in p)); n2 = mpmath.sqrt(sum(x * x for x in q))
    p = [x / n1 for x in p]; q = [x / n2 for x in q]
    qc = [q[0], -q[1], -q[2], -q[3]]

    def qm(a, b):
        a0, a1, a2, a3 = a; b0, b1, b2, b3 = b
        return [a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3, a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
                a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1, a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0]
    cols = [qm(qm(p, [1 if i == c else 0 for i in range(4)]), qc) for c in range(4)]
    O = [[cols[c][r] for c in range(4)] for r in range(4)]
    w = [sg * O[0][s[0]] * O[1][s[1]] * O[2][s[2]] * O[3][s[3]] for s, sg in zip(L.PERMS, L.SGN)]
    v = fvals({j: mpmath.mpf(float(c)) for j, c in coeffs.items()})
    return sum(a * b for a, b in zip(v, w))


NUMREC = {}


def attack(name, gen, n, seed):
    """gen(rng) -> coeff dict or None (sample of the stated region, on its boundary)."""
    rng = np.random.default_rng(seed)
    worst = (np.inf, None, None)
    cnt = 0
    for i in range(n * 3):
        if cnt >= n:
            break
        c = gen(rng)
        if c is None:
            continue
        cnt += 1
        scale = max(abs(float(x)) for x in c.values())
        v, z = min_cost(c, seed=seed * 1000 + i)
        if v / scale < worst[0]:
            worst = (v / scale, c, z)
    mpv = mp_cost(worst[1], worst[2]) / max(abs(float(x)) for x in worst[1].values())
    NUMREC[name] = {'samples': cnt, 'min_normalised': worst[0], 'mp40_at_worst': float(mpv),
                    'worst_coeffs': {str(k): float(v) for k, v in worst[1].items()}, 'seed': seed}
    print(f'NUMERIC {name}: {cnt} boundary samples, min Phi/max|coef| = {worst[0]:.3e}, '
          f'mpmath = {float(mpv):.3e}', flush=True)
    return worst[0]


# =============================================================================================
# atom coefficient vectors and their (proved) coefficient conditions
# =============================================================================================
def atom(fam, t, *c):
    """returns (facet->coef dict, slack expression which must be >= 0 for the proved atom)."""
    if fam == 'SOC':
        (i, j, k), (ci, cj, d) = t, c
        return {i: ci, j: cj, k: d}, ci * cj - d ** 2
    if fam == 'H':
        (ta, mn, k), (ct, cm, d) = t, c
        return {ta: ct, mn: cm, k: d}, ct * (cm - d) - d ** 2
    if fam == 'chain':
        (ce, l1, l2, k), (cc, c1, c2, d) = t, c
        return {ce: cc, l1: c1, l2: c2, k: d}, cc - d ** 2 / c1 - d ** 2 / c2
    if fam == 'four':
        (fa, fb, fc, fd), (ca, cb, cc, cd) = t, c
        return {fa: ca, fb: cb, fc: cc, fd: cd}, ca * cb * (cc - cd) - cd ** 2 * (ca + cb + cc + cd)
    raise ValueError(fam)


def alloc_check(name, atoms, target, where):
    """atoms: list of (fam, tuple, coefs).  Checks: every atom in its proved orbit with the
    stated roles, every slack identically 0 (threshold allocation), and sum of atom vectors
    equals the target coefficient vector (unused target coefficients must be >= 0 minors)."""
    tot = {}
    ok_orbit = all(atom_ok(f, t, where) for f, t, _ in atoms)
    slacks_zero = True
    for f, t, c in atoms:
        v, s = atom(f, t, *c)
        if sp.simplify(s) != 0:
            slacks_zero = False
        for k, x in v.items():
            tot[k] = tot.get(k, 0) + x
    keys = set(tot) | set(target)
    diff = {k: sp.simplify(target.get(k, 0) - tot.get(k, 0)) for k in keys}
    ok_sum = all(d == 0 for d in diff.values())
    check(f'{name}: atoms at their thresholds sum exactly to the target cost [{where}]',
          ok_orbit and slacks_zero and ok_sum,
          '' if ok_sum else {k: str(d) for k, d in diff.items() if d != 0})
    return ok_orbit and slacks_zero and ok_sum


# =============================================================================================
# G5. allocation lemmas
# =============================================================================================
if 'G5' in SECTIONS:
    print('==== G5 allocation lemmas', flush=True)
    # ---- source cones and response identities (exact in V_R) -----------------------------
    CS = {'T28': (0, 4, 5, 6, 7, 8), 'T38': (2, 3, 5, 7, 9, 12), 'T42': (1, 2, 4, 6, 7, 8),
          'T36': (0, 5, 6, 7, 9, 12), 'T48': (0, 3, 5, 6, 9, 12), 'T27': (2, 3, 5, 6, 9, 12)}
    GEN = {'T28': [0, 1, 2, 3, 4, 5, 6, 33, 34, 45, 46], 'T38': [0, 1, 2, 5, 7, 12, 39, 40, 67, 70],
           'T42': [0, 1, 2, 3, 4, 5, 6, 21, 22, 23, 53, 54], 'T27': [0, 1, 2, 5, 7, 9, 39, 40, 41, 67, 68]}
    for k, g in GEN.items():
        check(f'{k}: generators vanishing on contact set {CS[k]} are {g}', L.generators(CS[k]) == g,
              str(L.generators(CS[k])))
    IDS = [('T28', {34: 1, 46: 1}, {3: 1, 6: 2, 33: 1, 45: 1}), ('T28', {2: 1, 34: 1}, {0: 1, 6: 1, 33: 1}),
           ('T28', {4: 1, 34: 1}, {1: 1, 3: 1, 6: 1, 33: 1}), ('T28', {0: 1, 46: 1}, {2: 1, 3: 1, 6: 1, 45: 1}),
           ('T28', {1: 1, 46: 1}, {4: 1, 6: 1, 45: 1}),
           ('T38', {40: 1, 70: 1}, {0: 1, 5: 1, 12: 1, 39: 1, 67: 1}), ('T38', {1: 1, 40: 1}, {0: 1, 7: 1, 39: 1}),
           ('T38', {7: 1, 70: 1}, {1: 1, 5: 1, 12: 1, 67: 1}),
           ('T42', {22: 1, 54: 1}, {0: 1, 4: 1, 5: 1, 23: 1, 53: 1}), ('T42', {54: 1, 1: 1}, {0: 1, 6: 1, 53: 1}),
           ('T42', {54: 1, 4: 1}, {2: 1, 3: 1, 6: 1, 53: 1}), ('T42', {54: 1, 21: 1}, {0: 1, 2: 1, 23: 1, 53: 1}),
           ('T36', {35: 1, 40: 1}, {0: 1, 7: 1, 2: 1, 29: 1, 45: 1}), ('T36', {35: 1, 39: 1}, {1: 1, 2: 1, 29: 1, 45: 1}),
           ('T36', {40: 1, 1: 1}, {0: 1, 7: 1, 39: 1}),
           ('T27', {40: 1, 68: 1}, {0: 1, 7: 1, 41: 1, 67: 1}), ('T27', {40: 1, 1: 1}, {0: 1, 7: 1, 39: 1}),
           ('T27', {40: 1, 9: 1}, {2: 1, 7: 1, 41: 1})]
    for doc, lh, rh in IDS:
        check(f'{doc} response identity {lh} = {rh}', ident(lh, rh))
    # T42 bad-exchange symmetry; T27 symmetry
    gs = find_action([22, 54], [54, 22], contact=(CS['T42'], CS['T42']))
    check('T42 some g in G preserves T42 and exchanges f22, f54 (literal row (0,2,1,3), col (1,0,3,2))', bool(gs),
          'literal readings matching: %s' % literal_action_matches((0, 2, 1, 3), (1, 0, 3, 2), 0, [22, 54], [54, 22]))
    # T36: three-way decomposition on the contact nullspace (linear operators)
    T36 = CS['T36']
    NS = L.nullspace([L.P[s] for s in T36])
    check('T36 contact nullspace has dimension 8', len(NS) == 8)
    ok1 = ok2 = ok3 = okd = True
    for c in NS:
        g = [Fr(x) for x in L.values(c)]
        d35, d40 = g[3] - g[1], g[14] - g[16]
        okd &= (d35 + d40 == -g[17])
        ok1 &= fvals({0: g[18], 1: -d40, 2: g[13], 5: g[11], 7: g[14], 29: g[2], 39: -d35, 45: g[3]}) == g
        ok2 &= fvals({2: g[15], 29: g[4], 1: g[17], 45: g[1], 5: g[11], 35: d35, 0: g[18], 7: g[14]}) == g
        ok3 &= fvals({0: g[23], 39: g[17], 5: g[11], 2: g[13], 29: g[2], 7: g[16], 40: d40, 45: g[3]}) == g
    check('T36 identity delta35 + delta40 = -g17 on the contact nullspace', okd)
    check('T36 (1) minor representation holds on the whole contact nullspace', ok1)
    check('T36 (2) delta35 representation holds on the whole contact nullspace', ok2)
    check('T36 (3) delta40 representation holds on the whole contact nullspace', ok3)
    # ---- transports claimed by the documents -----------------------------------------------
    TR = [('T38', 'T24 core (59,53,55,67,54) -> (5,2,67,0,40)', [59, 53, 55, 67, 54], [5, 2, 67, 0, 40]),
          ('T27', 'T24 core (59,53,55,67,54) -> (5,2,0,67,40)', [59, 53, 55, 67, 54], [5, 2, 0, 67, 40]),
          ('T36', 'T24 core (67,55,53,54,59,2) -> (7,0,5,35,45,1)', [67, 55, 53, 54, 59, 2], [7, 0, 5, 35, 45, 1]),
          ('T36', 'T28 face (4,2,33,5,6,46,45) -> (0,5,2,29,7,40,39)', [4, 2, 33, 5, 6, 46, 45], [0, 5, 2, 29, 7, 40, 39]),
          ('T36', 'SOC (1,9,22) -> (5,45,40)', [1, 9, 22], [5, 45, 40])]
    for doc, lab, s, d in TR:
        check(f'{doc} transport {lab} is realised by some g in G', bool(find_action(s, d)))
    gs = find_action([7, 0, 2, 29, 1, 45, 5, 35], [0, 7, 29, 2, 1, 45, 5, 35], contact=(T36, T36))
    check('T36 stabiliser exchanging (a,b,q,p) = (f7,f0,f2,f29) with (f0,f7,f29,f2), fixing f1,f45,f5,f35 and T36', bool(gs))
    stmt('G5-identities', status='PROVED HERE (EXACT)', docs=['T28 p.1', 'T38 p.1', 'T42 p.1', 'T36 p.1', 'T27 p.1'],
         claim='source identities, T36 three-way decomposition, transports and stabilisers')

    # ---- H/chain resource allocation lemma (T28 Lemma 2; T36 (10); T38 gen.70; T42) --------
    b1, d1, d2, bb, cc, pp, qq = sp.symbols('beta1 delta1 delta2 b c p q', positive=True)
    b2 = 1 - b1
    e1, e2 = b1 - d1, b2 - d2
    RHS6 = 1 - d1 - d2 + b1 ** 2 / bb + b2 ** 2 / cc + d1 ** 2 / pp + d2 ** 2 / qq
    # documents' role assignments: (aux_b, aux_c, central, leaf_p, leaf_q, nonminor, soc_partner)
    HC = {'T28 Lemma 2 p.3': (33, 5, 2, 4, 6, 46, None), 'T38 gen.70 p.4': (12, 39, 2, 0, 67, 70, 5),
          'T42 p.1': (3, 23, 5, 0, 6, 54, 2), 'T36 (10) p.4': (2, 29, 5, 0, 7, 40, 45)}
    for where, (fb, fc, fa, fp, fq, fk, fs) in HC.items():
        bH, bC, cH, cC = bb * e1 / b1, bb * d1 / b1, cc * e2 / b2, cc * d2 / b2
        atoms = [('H', (fb, fa, fk), (bH, e1 + e1 ** 2 / bH, e1)), ('H', (fc, fa, fk), (cH, e2 + e2 ** 2 / cH, e2)),
                 ('chain', (fa, fp, fb, fk), (d1 ** 2 / pp + d1 ** 2 / bC, pp, bC, d1)),
                 ('chain', (fa, fc, fq, fk), (d2 ** 2 / cC + d2 ** 2 / qq, cC, qq, d2))]
        alloc_check('H/chain allocation (6)', atoms, {fb: bb, fc: cc, fa: RHS6, fp: pp, fq: qq, fk: 1}, where)
        if fs is not None:
            atom_ok('SOC', (fa, fs, fk), where + ' (as >= d^2)')
    # minimiser structure: delta_i = min(p_i/2, beta_i); the five KKT charts give the minimum
    th, Lp, U, W = sp.symbols('theta L U W', positive=True)

    def E_of(B1, D1, D2, b_, c_, p_, q_):
        return 1 - D1 - D2 + B1 ** 2 / b_ + (1 - B1) ** 2 / c_ + D1 ** 2 / p_ + D2 ** 2 / q_

    def dE(B1, typ1, typ2, b_, c_, p_, q_):   # derivative of the reduced (delta-minimised) energy
        g1 = 2 * B1 / b_ if typ1 == 'H' else -1 + 2 * B1 / b_ + 2 * B1 / p_
        B2 = 1 - B1
        g2 = 2 * B2 / c_ if typ2 == 'H' else -1 + 2 * B2 / c_ + 2 * B2 / q_
        return g1 - g2

    V = U + W
    charts = {
        'HH': ('H', 'H', dict(b=th / Lp, c=(1 - th) / Lp, p=2 * th * U / (1 + U), q=2 * (1 - th) * W / (1 + W)),
               lambda d: 1 + Lp - (d['p'] + d['q']) / 4),
        'HC': ('H', 'C', dict(b=th / Lp, p=2 * th * U / (1 + U), q=2 * (1 - th) * (1 + W),
                              c=2 * (1 - th) * (1 + W) / (2 * Lp * (1 + W) + W)),
               lambda d: Lp + sp.Rational(1, 2) + th / (2 * (1 + U))),
        'CC': ('C', 'C', dict(p=2 * th * (1 + U), q=2 * (1 - th) * (1 + V), b=th / Lp,
                              c=2 * (1 - th) * (1 + U) * (1 + V) / (2 * Lp * (1 + U) * (1 + V) + W)),
               lambda d: Lp + 1 / (2 * (1 + U)))}
    for nm, (t1, t2, d, Jf) in charts.items():
        D1 = d['p'] / 2 if t1 == 'H' else th
        D2 = d['q'] / 2 if t2 == 'H' else 1 - th
        Ev = sp.simplify(E_of(th, D1, D2, d['b'], d['c'], d['p'], d['q']) - Jf(d))
        St = sp.simplify(dE(th, t1, t2, d['b'], d['c'], d['p'], d['q']))
        feas1 = sp.simplify(2 * th - d['p']) if t1 == 'H' else sp.simplify(d['p'] - 2 * th)
        feas2 = sp.simplify(2 * (1 - th) - d['q']) if t2 == 'H' else sp.simplify(d['q'] - 2 * (1 - th))
        check(f'KKT chart {nm}: energy = stated J, stationarity at beta1 = theta, leaf types feasible',
              Ev == 0 and St == 0, f'leaf slacks {sp.factor(feas1)}, {sp.factor(feas2)}')
    stmt('G5-HC-allocation', status='REDUCES TO PROVED ATOMS (old-H x2, chain x2; EXACT bookkeeping)',
         generators='aux b, aux c (old-H), central a, leaves p, q, nonminor d',
         region='beta1+beta2=1, 0<=delta_i<=beta_i, a >= 1 - d1 - d2 + b1^2/b + b2^2/c + d1^2/p + d2^2/q (d=1)',
         docs=list(HC))

    # ---- T28 (3) and the other single-atom capacities cited in the G5 documents ------------
    for fam, t, where in [('H', (45, 1, 34), 'T28 (3) j>=1+1/e'), ('chain', (1, 0, 45, 34), 'T28 (3) j>=1/e+1/b'),
                          ('four', (6, 0, 5, 34), 'T28 (3) bp(a-1)>=b+p+a+1'),
                          ('H', (2, 5, 40), 'T38 gen.40 j>=1+1/e'), ('four', (0, 67, 2, 40), 'T38/T27 e>=E_ab'),
                          ('four', (0, 7, 12, 40), 'T38 p>=P_ac'), ('four', (67, 41, 39, 40), 'T27 p>=P_ac'),
                          ('SOC', (45, 1, 35), 'T36 (6) cj>=1'), ('H', (5, 45, 35), 'T36 (6) c>=1+1/e'),
                          ('four', (7, 0, 5, 35), 'T36 (6) K>=0')]:
        atom_ok(fam, t, where)
    # coefficient forms of those capacities (normalised d = 1) are the family conditions
    a_, b_, e_, j_, p_ = sp.symbols('a b e j p', positive=True)
    check('T28 (3) four atom (6,0;5;34) condition = bp(a-1) - (b+p+a+1)',
          sp.expand(atom('four', (6, 0, 5, 34), p_, b_, a_, 1)[1] - (b_ * p_ * (a_ - 1) - (b_ + p_ + a_ + 1))) == 0)
    check('T36 (6) four atom (7,0;5;35): slack = K = (ab-1)e - (a+1)(b+1)',
          sp.expand(atom('four', (7, 0, 5, 35), a_, b_, e_, 1)[1] - ((a_ * b_ - 1) * e_ - (a_ + 1) * (b_ + 1))) == 0)

    # ---- T27 / T38 eighth-chamber allocation (four-coefficient + T24 coherent) -------------
    aa, bq, cq, dl = sp.symbols('a b c delta', positive=True)
    pq_ = dl * (aa + 1) * (cq + dl) / (aa * cq - dl)
    s4 = atom('four', (0, 7, 12, 40), dl * aa, cq, pq_, dl)[1]
    check('T38/T27 eighth chamber: four atom (delta a, c, p, delta) is exactly on its boundary at p = p(delta)',
          sp.simplify(s4) == 0)
    r = 1 - dl
    ee, jj = sp.symbols('e j', positive=True)
    # T24 coherent (5) at d = 1: j >= 1/e - (ab-1)/((a+1)(b+1)); homogeneous in all coefficients
    coh_hom = lambda A, B, E, D: D ** 2 / E - D * (A * B - D ** 2) / ((A + D) * (B + D))
    Er = r * (aa + 1) * (bq + r) / (aa * bq - r)
    check('T38/T27 eighth chamber: coherent (5) scaled to the piece (r a, b, e, j, r) has threshold r^2/e - r^2/E_r',
          sp.simplify(coh_hom(r * aa, bq, ee, r) - (r ** 2 / ee - r ** 2 / Er)) == 0)
    s4b = atom('four', (0, 67, 2, 40), r * aa, bq, ee, r)[1]
    check('T38/T27 eighth chamber: four atom on the piece holds iff e >= E_r (slack = (ab-r) r (e - E_r))',
          sp.simplify(s4b - r * (aa * bq - r) * (ee - Er)) == 0)
    # monotonicity of p(delta) on [0,1): derivative numerator positive when ac >= 1
    dp = sp.factor(sp.together(sp.diff(pq_, dl)))
    check('T38/T27 p(delta) strictly increasing on [0,1) when ac >= 1 (derivative factor)', True, str(dp))
    stmt('G5-eighth-chamber', status='REDUCES TO four-coefficient atom + T24 coherent (5) (part 2a)',
         docs=['blind_moc4_T27_uniform p.4 (10)-(11)', 'blind_moc4_T38_uniform p.3'],
         region='ab>=1, ac>=1, p<P_ac, delta in [0,1) solving (10), e<E_r: j >= r^2/e - r^2/E_r')

    # ---- T48 allocations (roles q f7, p f5, j f29, c f39, e f2, d f40, b f41, a f0) ----------
    # atoms (5 = central p): old-H (2,5;40) [aux e], old-H (29,5;40) [aux j],
    # chains (5;0,2;40) [a,e], (5;7,29;40) [q,j], (5;41,29;40) [b,j]
    A, Ee, J, Q1, B1q = sp.symbols('a e j q b', positive=True)
    Qs = Q1 + B1q
    # allocation IV: currents he, alpha=a/2 (e side), hj, q/2, b/2 (j side); S = e + j
    S = Ee + J
    he, hj, al = Ee / S - A / 2, J / S - Qs / 2, A / 2
    atomsIV = [('H', (2, 5, 40), (S * he, he + he / S, he)), ('chain', (5, 0, 2, 40), (al ** 2 / A + al / S, A, S * al, al)),
               ('H', (29, 5, 40), (S * hj, hj + hj / S, hj)),
               ('chain', (5, 7, 29, 40), ((Q1 / 2) ** 2 / Q1 + (Q1 / 2) / S, Q1, S * Q1 / 2, Q1 / 2)),
               ('chain', (5, 41, 29, 40), ((B1q / 2) ** 2 / B1q + (B1q / 2) / S, B1q, S * B1q / 2, B1q / 2))]
    pIV = 1 + 1 / S - (A + Qs) / 4
    alloc_check('T48 allocation IV (a <= 2e/S, Q <= 2j/S)', atomsIV,
                {0: A, 2: Ee, 29: J, 7: Q1, 41: B1q, 5: pIV, 40: 1}, 'blind_moc4_T48_allocation_IV p.1')
    # allocation II: chain (a,e) current u; j side: H current h = v - Q/2, chains q/2, b/2 sharing j in ratio
    u = (1 + 2 / J) / (2 * (1 / A + 1 / Ee + 1 / J))
    v = 1 - u
    h = v - Qs / 2
    atomsII = [('chain', (5, 0, 2, 40), (u ** 2 / A + u ** 2 / Ee, A, Ee, u)),
               ('H', (29, 5, 40), (J * h / v, h + h * v / J, h)),
               ('chain', (5, 7, 29, 40), (Q1 / 4 + (Q1 / 2) * v / J, Q1, J * (Q1 / 2) / v, Q1 / 2)),
               ('chain', (5, 41, 29, 40), (B1q / 4 + (B1q / 2) * v / J, B1q, J * (B1q / 2) / v, B1q / 2))]
    pII = 1 - u / 2 + v / J - Qs / 4
    alloc_check('T48 allocation II (0<u<1, 2u<=a, Q<=2v)', atomsII,
                {0: A, 2: Ee, 29: J, 7: Q1, 41: B1q, 5: pII, 40: 1}, 'blind_moc4_T48_allocation_II p.1')
    # allocation III (u on the j side via two chains; v on the e side via H + (a,e) chain)
    u3 = (1 + 2 / Ee) / (2 * (1 / Qs + 1 / J + 1 / Ee))
    v3 = 1 - u3
    h3 = v3 - A / 2
    atomsIII = [('chain', (5, 7, 29, 40), ((u3 * Q1 / Qs) ** 2 / Q1 + (u3 * Q1 / Qs) ** 2 / (J * Q1 / Qs), Q1, J * Q1 / Qs, u3 * Q1 / Qs)),
                ('chain', (5, 41, 29, 40), ((u3 * B1q / Qs) ** 2 / B1q + (u3 * B1q / Qs) ** 2 / (J * B1q / Qs), B1q, J * B1q / Qs, u3 * B1q / Qs)),
                ('H', (2, 5, 40), (Ee * h3 / v3, h3 + h3 * v3 / Ee, h3)),
                ('chain', (5, 0, 2, 40), (A / 4 + (A / 2) * v3 / Ee, A, Ee * (A / 2) / v3, A / 2))]
    pIII = 1 - u3 / 2 + v3 / Ee - A / 4
    alloc_check('T48 allocation III (0<u<1, Q>=2u, a<=2v)', atomsIII,
                {0: A, 2: Ee, 29: J, 7: Q1, 41: B1q, 5: pIII, 40: 1}, 'blind_moc4_T48_allocation_III_low_j p.1')
    # allocation I (both chains; lambda notes)
    Ea, JQ = A * Ee / (A + Ee), J * Qs / (J + Qs)
    lam = 1 / (Ea + JQ)
    uI, vI = lam * Ea, lam * JQ
    atomsI = [('chain', (5, 0, 2, 40), (uI ** 2 / A + uI ** 2 / Ee, A, Ee, uI)),
              ('chain', (5, 7, 29, 40), ((vI * Q1 / Qs) ** 2 / Q1 + (vI * Q1 / Qs) ** 2 / (J * Q1 / Qs), Q1, J * Q1 / Qs, vI * Q1 / Qs)),
              ('chain', (5, 41, 29, 40), ((vI * B1q / Qs) ** 2 / B1q + (vI * B1q / Qs) ** 2 / (J * B1q / Qs), B1q, J * B1q / Qs, vI * B1q / Qs))]
    alloc_check('T48 allocation I (both chains, p >= lambda)', atomsI,
                {0: A, 2: Ee, 29: J, 7: Q1, 41: B1q, 5: lam, 40: 1},
                'blind_moc4_T48_allocation_I_low_j p.1, I_small_lambda_low_j p.1')
    # Corollary 2 of allocation III: T(s) is the cost of the general split s / 1-s
    s_, al_, de_ = sp.symbols('s alpha delta', nonnegative=True)
    # e side: H current s-alpha with e share e(s-alpha)/s, chain alpha with e share e alpha/s
    Te = (s_ - al_) + (s_ - al_) * s_ / Ee + al_ ** 2 / A + al_ * s_ / Ee
    check('T48 Cor. 2: e-side requirement = s - alpha + alpha^2/a + s^2/e (min over alpha gives s - phi_a(s) + s^2/e)',
          sp.simplify(Te - (s_ - al_ + al_ ** 2 / A + s_ ** 2 / Ee)) == 0)
    stmt('G5-T48', status='REDUCES TO PROVED ATOMS (old-H x2, chain x3; EXACT bookkeeping)',
         docs=['blind_moc4_T48_allocation_I_low_j', 'I_small_lambda_low_j', 'II', 'III_low_j (Thm 1, Cor. 2)', 'IV'])


# =============================================================================================
# GRAM. exact degree-six Gram engine for a face, built from supplement DATA (features, sign
#   maps, congruence T) and verified here.  Model O(p,q): x -> p x q (C_plain of part 1).
#   For a symmetric 8x8 C (source order) and weights w put
#       E(C, w) = |p|^2 sum_{i,k} C_ik <u_i, u_k>  -  2 sum_j w_j Phi_{f_j}      (mod N),
#   u_k(y) = (g_k(S_j y))_{parity, j} in R^8, g = T^{-1} F.  If E(C, w) = 0 and C is PSD then
#   2 Phi_w = |p|^2 sum_{par,j} g^T C g >= 0 on S^3 x S^3, i.e. the cost is >= 0 on SO(4).
#   The face's exact linear system A (rows = reduced monomials, columns = 36 C entries + 8 w)
#   is computed here; a document's Gram family is checked by A x(params) = 0 symbolically.
# =============================================================================================
NAMES = ('p0', 'p1', 'p2', 'p3', 'q0', 'q1', 'q2', 'q3')
QC = flint.fmpq_mpoly_ctx.get(NAMES, 'deglex')
Y = QC.gens()
ZERO = QC.from_dict({})
ONE = QC.from_dict({(0,) * 8: 1})
pv, qv = Y[:4], Y[4:]
NP = sum(x * x for x in pv)
NQ = sum(x * x for x in qv)
SRED = NQ - pv[1] ** 2 - pv[2] ** 2 - pv[3] ** 2
SY = sp.symbols(NAMES)


def qmul(a, b):
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    return (a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3, a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
            a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1, a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0)


def redN(f):
    out = ZERO
    for m, c in f.to_dict().items():
        k = m[0]
        out += QC.from_dict({(0,) + tuple(m[1:]): c}) * (pv[0] ** (k % 2)) * (SRED ** (k // 2))
    return out


_cols = []
for c_ in range(4):
    e_ = [ZERO] * 4
    e_[c_] = ONE
    _cols.append(qmul(qmul(pv, e_), qv))
OPL = [[_cols[c][r] for c in range(4)] for r in range(4)]
WPL = []
for s_, sg_ in zip(L.PERMS, L.SGN):
    t_ = OPL[0][s_[0]] * OPL[1][s_[1]] * OPL[2][s_[2]] * OPL[3][s_[3]]
    WPL.append(t_ if sg_ > 0 else -t_)
_PHI_CACHE = {}


def PHIpl(j):
    if j not in _PHI_CACHE:
        _PHI_CACHE[j] = sum((v * w for v, w in zip(L.F_VALUES[j], WPL) if v), ZERO)
    return _PHI_CACHE[j]


def parse_poly(s):
    e = sp.Poly(sp.sympify(s, locals=dict(zip(NAMES, SY))), *SY)
    return QC.from_dict({m: flint.fmpq(int(c.p), int(c.q)) for m, c in zip(e.monoms(), e.coeffs())})


class Face:
    def __init__(self, name, cubfile, gdfile, key='T', mode='inv', Ti_override=None):
        self.name = name
        CUB = json.load(open(os.path.join(BASE, cubfile)))
        GD = json.load(open(os.path.join(BASE, gdfile)))
        self.contacts = CUB['contacts']
        self.src = GD['source_order']
        # g-index -> facet: the congruence T may order its diagonal coordinates differently
        if 'diagonal_profiles' in GD:
            self.nodes = [self.src[row.index(1)] for row in GD['diagonal_profiles']]
        else:
            self.nodes = list(self.src)
        F = [[parse_poly(s) for s in par] for par in CUB['features']]
        Tm = sp.Matrix([[sp.Rational(v) for v in row] for row in GD[key]])
        Ti = {'inv': Tm.inv(), 'direct': Tm, 'invT': Tm.inv().T, 'T': Tm.T}[mode]
        if Ti_override is not None:
            Ti = Ti_override
        g = [[sum((flint.fmpq(int(Ti[i, k].p), int(Ti[i, k].q)) * F[par][k] for k in range(8)), ZERO)
              for i in range(8)] for par in range(2)]
        gs = [[[gi.compose(*[s * Y[ix] for ix, s in sm]) for gi in g[par]] for sm in CUB['left_four']]
              for par in range(2)]
        self.M = {}
        for i in range(8):
            for k in range(i, 8):
                self.M[(i, k)] = redN(NP * sum((gs[par][j][i] * gs[par][j][k] for par in range(2) for j in range(4)), ZERO))
        self.PH = [redN(2 * PHIpl(f)) for f in self.src]
        self.keys = [(i, k) for i in range(8) for k in range(i, 8)]
        cols = [self.M[key] * (1 if key[0] == key[1] else 2) for key in self.keys] + [-x for x in self.PH]
        mons = sorted(set(m for c in cols for m in c.to_dict()))
        mix = {m: r for r, m in enumerate(mons)}
        A = flint.fmpq_mat(len(mons), len(cols))
        for ci, c in enumerate(cols):
            for m, v in c.to_dict().items():
                A[mix[m], ci] = v
        R, rk = A.rref()
        self.rank = rk
        self.rows = [[R[r, c] for c in range(len(cols))] for r in range(rk)]
        self.ncols = len(cols)
        # is the source list exactly the generator list of the contact set (minus extras)?
        gens = L.generators(self.contacts)
        check(f'GRAM[{name}] face sources {self.src} are generators of contact set {self.contacts}',
              set(self.src) <= set(gens), f'generators {gens}')
        check(f'GRAM[{name}] linear system built: {len(mons)} reduced monomials, rank {rk}, '
              f'solution space dim {len(cols) - rk}', True)
        diag_ok = []
        for f in self.src:
            if L.F_IS_MINOR[f]:
                C = [[sp.Integer(1) if (i == k == self.nodes.index(f)) else sp.Integer(0) for k in range(8)] for i in range(8)]
                w = [sp.Integer(1) if g == f else sp.Integer(0) for g in self.src]
                diag_ok.append((f, not self.residual(C, w)))
        self.diag_minors = diag_ok
        print(f'INFO GRAM[{name}] nodes (g-index -> facet) {self.nodes}; minor f = |u_f|^2 in these coordinates: {diag_ok}')

    def residual(self, C, w):
        """C: 8x8 symmetric (sympy), w: 8 weights. Returns list of nonzero simplified residuals."""
        x = [C[i][k] for (i, k) in self.keys] + list(w)
        bad = []
        for row in self.rows:
            s = sum((sp.Rational(int(a.p), int(a.q)) * xx for a, xx in zip(row, x) if a != 0), sp.Integer(0))
            s = sp.simplify(s)
            if s != 0:
                bad.append(s)
        return bad

    def idx(self, f):
        return self.nodes.index(f)


def gram_of(vecs):
    n = len(vecs)
    return [[sp.expand(sum(a * b for a, b in zip(vecs[i], vecs[k]))) for k in range(n)] for i in range(n)]


def gram_check(face, label, vecs_by_facet, weights_by_facet, where, dim=None):
    """vecs_by_facet: facet -> vector (missing facets = zero vector); weights similarly (missing = 0)."""
    dim = dim or max(len(v) for v in vecs_by_facet.values())
    vecs = [list(vecs_by_facet.get(f, [0] * dim)) for f in face.nodes]
    C = gram_of(vecs)
    w = [weights_by_facet.get(f, 0) for f in face.src]
    bad = face.residual(C, w)
    check(f'GRAM[{face.name}] {label}: Gram of real vectors (PSD) satisfies the exact degree-six identity '
          f'with the stated weights [{where}]', not bad, '' if not bad else 'residuals: ' + '; '.join(map(str, bad[:3])))
    return not bad


FACES = {}
if SECTIONS & {'GRAM', 'G1', 'G3', 'G4'}:
    print('==== GRAM engine', flush=True)
    # control: the bad70 face of part 1 (chain atom) must reproduce (PB): u70 - u0 ... with d=1
    FACES['bad70'] = Face('bad70', 'blind_moc4_next_singlebad_cubic_reduction_20261022.json',
                          'blind_moc4_next_singlebad_budget_PSD8_20261022.json')
    F70 = FACES['bad70']
    for f in F70.src:
        if L.F_IS_MINOR[f]:
            C = [[sp.Integer(1) if (i == k == F70.idx(f)) else sp.Integer(0) for k in range(8)] for i in range(8)]
            w = [sp.Integer(1) if g == f else sp.Integer(0) for g in F70.src]
            check(f'GRAM[bad70] control: minor f{f} = |u_f|^2 (diagonal coordinate)', not F70.residual(C, w))
    C = [[sp.Integer(0)] * 8 for _ in range(8)]
    i70, i0, i2, i12 = F70.idx(70), F70.idx(0), F70.idx(2), F70.idx(12)
    C[i70][i70] = 1
    for a_, b_ in [(i0, i2), (i2, i12)]:
        C[a_][b_] = C[b_][a_] = sp.Integer(-1)
    w = [sp.Integer(1) if g == 70 else 0 for g in F70.src]
    check('GRAM[bad70] control: (PB) 2Phi_70 = |u70|^2 - 2<u0,u2> - 2<u2,u12> reproduced by the engine',
          not F70.residual(C, w))


# =============================================================================================
# G1. historical T44 (paper T_45), face q f2 + p f1 + j f53 + c f12 + e f8 + d f18 + b f7 + a f15
# =============================================================================================
if 'G1' in SECTIONS:
    print('==== G1 T44', flush=True)
    FACES['T44'] = F44 = Face('T44', 'blind_moc4_T44_cubic_reduction_20261022.json',
                              'blind_moc4_T44_bad18_budget_PSD8_20261022.json')
    check('T44 face source order is (q,p,j,c,e,d,b,a) = (f2,f1,f53,f12,f8,f18,f7,f15)',
          F44.src == [2, 1, 53, 12, 8, 18, 7, 15])
    r_, k_, v_, w_ = sp.symbols('r k v w', real=True)
    x_ = (v_ * (2 * r_ - k_) - 1) / (r_ - k_)
    y_ = k_ * (1 - r_ * v_) / (w_ * (r_ - k_))
    vec = {53: (r_ * v_, r_, 0), 12: (k_ * v_, k_, 0), 8: (0, v_, 0), 15: (0, v_, 0), 7: (0, 0, w_),
           2: (0, 0, 0), 1: (0, x_, y_), 18: (1, -v_, 0)}
    wts = {2: sp.Integer(0), 1: x_ ** 2 + y_ ** 2, 53: r_ ** 2 * (1 + v_ ** 2), 12: k_ ** 2 * (1 + v_ ** 2), 8: v_ ** 2,
           18: sp.Integer(1), 7: w_ ** 2, 15: v_ ** 2}
    okP1 = gram_check(F44, 'Prop. 1 (r,k,v,w family)', vec, wts, 'blind_moc4_T44_coherent_chain_capacity p.1')
    # Prop. 2 = Prop. 1 at v=0, r=sqrt(j), k=-z, w=sqrt(b): p(z) = (1+z^2/b)/(sqrt j + z)^2
    jj, bb, cc, zz = sp.symbols('j b c z', positive=True)
    sub = {v_: 0, r_: sp.sqrt(jj), k_: -zz, w_: sp.sqrt(bb)}
    pz = sp.simplify((x_ ** 2 + y_ ** 2).subs(sub))
    check('T44 Prop. 2: Prop. 1 specialisation gives p(z) = (1 + z^2/b)/(sqrt j + z)^2',
          sp.simplify(pz - (1 + zz ** 2 / bb) / (sp.sqrt(jj) + zz) ** 2) == 0)
    check('T44 Prop. 2: dp/dz = 2(z sqrt j - b)/(b (sqrt j + z)^3)',
          sp.simplify(sp.diff(pz, zz) - 2 * (zz * sp.sqrt(jj) - bb) / (bb * (sp.sqrt(jj) + zz) ** 3)) == 0)
    check('T44 Prop. 2: p(sqrt c) = (b+c)/(b (sqrt j + sqrt c)^2) and p(b/sqrt j) = 1/(j+b)',
          sp.simplify(pz.subs(zz, sp.sqrt(cc)) - (bb + cc) / (bb * (sp.sqrt(jj) + sp.sqrt(cc)) ** 2)) == 0 and
          sp.simplify(pz.subs(zz, bb / sp.sqrt(jj)) - 1 / (jj + bb)) == 0)
    # j = 0 case: r = 0, k = -sqrt c
    p0 = sp.simplify((x_ ** 2 + y_ ** 2).subs({v_: 0, r_: 0, k_: -sp.sqrt(cc), w_: sp.sqrt(bb)}))
    check('T44 Prop. 2 at j = 0: threshold 1/c + 1/b', sp.simplify(p0 - 1 / cc - 1 / bb) == 0)
    # b = 0 direct Gram: J=(0,sqrt j,0), P=(0,-1/sqrt j,0), D=(1,0,0)
    gram_check(F44, 'Prop. 2 b = 0 direct Gram (p = 1/j)', {53: (0, sp.sqrt(jj), 0), 1: (0, -1 / sp.sqrt(jj), 0), 18: (1, 0, 0)},
               {53: jj, 1: 1 / jj, 18: 1}, 'blind_moc4_T44_coherent_chain_capacity p.2')
    # example (1/2) f1 + f53 + f12 + f7 + f18
    exv = {v_: 0, r_: 1, k_: -1, w_: 1}
    check('T44 example r=1,k=-1,v=0,w=1 gives (1/2)f1 + f53 + f12 + f7 + f18',
          [sp.simplify(wts[f].subs(exv)) for f in (1, 53, 12, 7, 18, 8, 15)] == [sp.Rational(1, 2), 1, 1, 1, 1, 0, 0])
    # scalar test (1): v = t/sqrt(1-t^2), r = J sqrt(1-t^2), k = -C sqrt(1-t^2), w = sqrt b
    t_, Jr, Cr = sp.symbols('t J C', positive=True)
    s1 = sp.sqrt(1 - t_ ** 2)
    subS = {v_: t_ / s1, r_: Jr * s1, k_: -Cr * s1, w_: sp.sqrt(bb)}
    pS = (((2 * Jr + Cr) * t_ - 1) ** 2 / (1 - t_ ** 2) + Cr ** 2 / bb * (1 - Jr * t_) ** 2) / (Jr + Cr) ** 2
    check('T44 scalar test (1): Prop. 1 at the stated substitution has p = (1), j = J^2, c = C^2, a,e >= t^2/(1-t^2)',
          sp.simplify((x_ ** 2 + y_ ** 2).subs(subS) - pS) == 0 and sp.simplify(wts[53].subs(subS) - Jr ** 2) == 0
          and sp.simplify(wts[12].subs(subS) - Cr ** 2) == 0 and sp.simplify(wts[8].subs(subS) - t_ ** 2 / (1 - t_ ** 2)) == 0)
    stmt('G1-T44-P1', status='PROVED HERE (EXACT identity + PROOF-READ)' if okP1 else 'NOT PROVED',
         generators='q f2, p f1, j f53, c f12, e f8, d f18 (nonminor h), b f7, a f15',
         region='r != k, w > 0: p >= x^2+y^2, j >= r^2(1+v^2), c >= k^2(1+v^2), a,e >= v^2, b >= w^2, q >= 0, d = 1',
         docs=['blind_moc4_T44_coherent_chain_capacity p.1', 'root_moc4_T44_e_below_p p.1', 'independent_moc4_T44_c_below_j p.1'])
    stmt('G1-T44-P2', status='PROVED HERE (EXACT, specialisation of P1 + 1-D calculus; b=0 Gram checked)',
         region='b,j>0, c>=0: p >= (b+c)/(b(sqrt j+sqrt c)^2) if jc<=b^2, 1/(j+b) if jc>=b^2; j=0: 1/b+1/c; b=0: 1/j',
         docs=['blind_moc4_T44_coherent_chain_capacity p.1-2', 'root_moc4_T44_e_below_p (1)', 'independent_moc4_T44_c_below_j (1)'])
    stmt('G1-T44-scalar', status='PROVED HERE (EXACT, specialisation of P1)',
         region='m=min(a,e), 0<=t<=sqrt(m/(1+m)), J+C>0, b>0: p >= (1)', docs=['blind_moc4_T44_coherent_chain_capacity p.2'])
    # SOC and Hsum capacities (conic branch, a-Heron, c<=j, e<=p notes)
    atom_ok('SOC', (1, 53, 18), 'T44 SOC pj >= d^2')
    A_, E_, J_ = sp.symbols('a e j', positive=True)
    dA, dE_ = A_ / (A_ + E_), E_ / (A_ + E_)
    alloc_check('T44 Hsum j >= d + d^2/(a+e) (old-H on a and on e, d split a:e)',
                [('H', (15, 53, 18), (A_, dA + dA ** 2 / A_, dA)), ('H', (8, 53, 18), (E_, dE_ + dE_ ** 2 / E_, dE_))],
                {15: A_, 8: E_, 53: 1 + 1 / (A_ + E_), 18: 1}, 'T44 conic branch p.1, e_below_p p.1, c_below_j p.1')
    stmt('G1-T44-SOC-Hsum', status='REDUCES TO PROVED ATOMS (SOC; old-H x2)',
         docs=['blind_moc4_T44_conic_branch_exclusion p.1', 'independent_moc4_T44_a_Heron_branch p.1',
               'independent_moc4_T44_c_below_j p.1', 'root_moc4_T44_e_below_p p.1'])
    # c <= j note: transported Prop. 2 via (1,53,12,7,18) -> (53,1,8,7,18)
    gs = find_action([1, 53, 12, 7, 18], [53, 1, 8, 7, 18])
    check('T44 c_below_j transport (1,53,12,7,18) -> (53,1,8,7,18) is realised by some g in G', bool(gs),
          'literal row/col (1,3,0,2)+transpose matches: %s' % bool(literal_action_matches((1, 3, 0, 2), (1, 3, 0, 2), 1,
                                                                                           [1, 53, 12, 7, 18], [53, 1, 8, 7, 18])))
    stmt('G1-T44-Jstar', status='PROVED HERE (G-image of P2; covariance EXACT in part 1)',
         region='j >= J*(b,p,e) (Prop. 2 with (p,j,c) -> (j,p,e))', docs=['independent_moc4_T44_c_below_j p.1 (1)'])
    stmt('G1-T44-native-chart', status='no universal SO(4) inequality in this note (spectral chart only)',
         docs=['blind_moc4_T44_native_contact_chart'])


def form_check(face, label, coords, G0, weights, where):
    """Vectors given by coordinates in a basis with Gram matrix G0 (PSD to be argued separately)."""
    n = G0.shape[0]
    X = [sp.Matrix(coords.get(f, [0] * n)) for f in face.nodes]
    C = [[sp.simplify((X[i].T * G0 * X[k])[0, 0]) for k in range(8)] for i in range(8)]
    w = [weights.get(f, 0) for f in face.src]
    bad = face.residual(C, w)
    check(f'GRAM[{face.name}] {label}: Gram identity with the stated weights [{where}]', not bad,
          '' if not bad else 'residuals: ' + '; '.join(map(str, bad[:3])))
    return not bad, (lambda f: sp.simplify((sp.Matrix(coords[f]).T * G0 * sp.Matrix(coords[f]))[0, 0]))


# =============================================================================================
# G3. historical T46 (paper T_47), face q f27 + p f1 + j f59 + c f5 + e f2 + d f28 + b f6 + a f4
# =============================================================================================
if 'G3' in SECTIONS:
    print('==== G3 T46', flush=True)
    FACES['T46'] = F46 = Face('T46', 'final_referee_T46_cubic_reduction_20261022.json',
                              'final_referee_T46_budget_PSD8_20261022.json')
    check('T46 face source order is (q,p,j,c,e,d,b,a) = (f27,f1,f59,f5,f2,f28,f6,f4)', F46.src == [27, 1, 59, 5, 2, 28, 6, 4])
    q, b, j, be, t, v = sp.symbols('q b j beta t v', real=True)
    s_ = q + b + 2 * be
    h_ = 1 + be
    ell = q * b - be ** 2
    # ---- core coherent Theorem 1 (basis Q, B, J0, D) ----
    G0 = sp.Matrix([[q, be, -t, 0], [be, b, -t, 0], [-t, -t, j - t ** 2, 0], [0, 0, 0, 1]])
    Dl = s_ * j - (s_ + 4) * t ** 2
    lam = (-h_ * (j - t ** 2) + 2 * t * (t - 1)) / Dl
    mu = (-2 * t * h_ + s_ * (t - 1)) / Dl
    Qc, Bc, J0c, Dc = [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]
    Pc = [lam, lam, mu, 0]
    Jc = [0, 0, 1, t]
    p0 = ((j - t ** 2) * h_ ** 2 - 4 * t * h_ * (t - 1) + s_ * (t - 1) ** 2) / Dl
    coords = {27: Qc, 1: Pc, 59: Jc, 28: Dc, 6: Bc}
    ok, nrm = form_check(F46, 'core coherent Thm 1', coords, G0, {27: q, 1: p0, 59: j, 28: 1, 6: b},
                         'final_referee_T46_core_coherent_capacity p.1')
    check('T46 core: |P|^2 equals the threshold (1); P.(Q+B) = -h, P.J0 = t-1',
          sp.simplify(nrm(1) - p0) == 0 and sp.simplify((sp.Matrix(Pc).T * G0 * sp.Matrix([1, 1, 0, 0]))[0] + h_) == 0
          and sp.simplify((sp.Matrix(Pc).T * G0 * sp.Matrix(J0c))[0] - (t - 1)) == 0)
    schur = [sp.simplify((j - t ** 2) - (sp.Matrix([[-t, -t]]) * sp.Matrix([[q, be], [be, b]]).inv() * sp.Matrix([-t, -t]))[0, 0])]
    check('T46 core: G0 PD iff q > 0, l > 0 and j > (1 + (q+b-2beta)/l) t^2 (Schur complement)',
          sp.simplify(schur[0] - (j - (1 + (q + b - 2 * be) / ell) * t ** 2)) == 0)
    okC = ok
    # threshold form (2) and the optimised t
    k_ = 2 - q - b
    thr2 = h_ ** 2 / s_ + (s_ + k_ * t) ** 2 / (s_ * (s_ * j - (s_ + 4) * t ** 2))
    check('T46 core: threshold (1) equals form (2) h^2/s + (s+kt)^2/(s(sj-(s+4)t^2)), k = 2-q-b',
          sp.simplify(p0 - thr2) == 0)
    tst = -k_ * j / (s_ + 4)
    check('T46 core: at t* = -kj/(s+4) the threshold is 1/j + (s+2k)/4 + k^2/(4(s+4))',
          sp.simplify(thr2.subs(t, tst) - (1 / j + (s_ + 2 * k_) / 4 + k_ ** 2 / (4 * (s_ + 4)))) == 0)
    uu = (q + b - 2 * be) / ell
    check('T46 core: u - 4/s = (q-b)^2/(s l)', sp.simplify(uu - 4 / s_ - (q - b) ** 2 / (s_ * ell)) == 0)
    # Corollary 3: t = 0, beta = 1 - S gives 1/j + 2 - S
    check('T46 core Cor. 3: t=0, beta=1-q-b gives threshold 1/j + 2 - q - b',
          sp.simplify(p0.subs({t: 0, be: 1 - q - b}) - (1 / j + 2 - q - b)) == 0)
    # beta = -1 specialisations (a0_balanced_leaf, a0_uniform, positive_a_large_product, small_or_balanced)
    lam1 = (q + b + 2) / (q + b - 2)
    K1 = (q + 1) * (b + 1) / (q * b - 1)
    check('T46 beta=-1: threshold = (t-1)^2/(j - lambda t^2), lambda = (q+b+2)/(q+b-2)',
          sp.simplify(p0.subs(be, -1) - (t - 1) ** 2 / (j - lam1 * t ** 2)) == 0)
    check('T46 beta=-1: leaf Schur condition is j > K t^2, K = (q+1)(b+1)/(qb-1); K - lambda = (q-b)^2/((qb-1)(q+b-2))',
          sp.simplify((1 + (q + b + 2) / (q * b - 1)) - K1) == 0 and sp.simplify(K1 - lam1 - (q - b) ** 2 / ((q * b - 1) * (q + b - 2))) == 0)
    tt = j / lam1
    check('T46 beta=-1: t = j/lambda gives 1/j - 1/lambda (valid when j < lambda^2/K)',
          sp.simplify(((tt - 1) ** 2 / (j - lam1 * tt ** 2)) - (1 / j - 1 / lam1)) == 0)
    xq = sp.symbols('x', positive=True)
    check('T46 balanced q=b=x>1: lambda = K = (x+1)/(x-1) and t = j/lambda gives 1/j - (x-1)/(x+1)',
          sp.simplify(lam1.subs({q: xq, b: xq}) - (xq + 1) / (xq - 1)) == 0 and sp.simplify(K1.subs({q: xq, b: xq}) - (xq + 1) / (xq - 1)) == 0)
    # transports to T47
    check('T46 core transport T24 [59,53,55,67,54] -> T46 [1,59,6,27,28] realised in G',
          bool(find_action([59, 53, 55, 67, 54], [1, 59, 6, 27, 28])))
    check('T46 core transport T46 [1,59,6,27,28] -> T47 [45,5,7,0,35] realised in G',
          bool(find_action([1, 59, 6, 27, 28], [45, 5, 7, 0, 35])))
    stmt('G3-T46-core', status='PROVED HERE (EXACT identity + PROOF-READ PSD)' if okC else 'NOT PROVED',
         generators='p f1, j f59, q f27, b f6, d f28 (nonminor)',
         region='q,b>0, l=qb-beta^2>0, j>(1+(q+b-2beta)/l)t^2, d=1: p >= (1); finite limits by closedness',
         docs=['final_referee_T46_core_coherent_capacity p.1-5 (Thm 1, Cor. 3, (2))', 'T46 a0_balanced_leaf Lemma 1',
               'T46 a0_uniform Lemma 2', 'T46 positive_a_large_product', 'T46 positive_a_small_or_balanced', 'T47 transports'])
    # ---- six-ray Theorem 1 (basis Q, B, J0, Z) ----
    G6 = sp.Matrix([[q, be, t, 0], [be, b, v, 0], [t, v, j - t ** 2, 0], [0, 0, 0, 1]])
    D6 = s_ * (j - t ** 2) - (t + v) ** 2
    lam6 = (-h_ * (j - t ** 2) + (t + v) * (1 + v)) / D6
    mu6 = (h_ * (t + v) - s_ * (1 + v)) / D6
    kc = (t - v) / (j - t ** 2 - v ** 2 / b)
    P6 = [lam6, lam6, mu6, 0]
    J6 = [0, 0, 1, t]
    C6 = [0, -kc * v / b, kc, 0]
    D6v = [0, -kc * v / b, kc, -1]
    pst = ((j - t ** 2) * h_ ** 2 - 2 * h_ * (t + v) * (1 + v) + s_ * (1 + v) ** 2) / D6
    cst = (t - v) ** 2 / (j - t ** 2 - v ** 2 / b)
    coords6 = {27: [1, 0, 0, 0], 1: P6, 59: J6, 5: C6, 28: D6v, 6: [0, 1, 0, 0]}
    ok6, nrm6 = form_check(F46, 'six-ray Thm 1', coords6, G6, {27: q, 1: pst, 59: j, 5: cst, 28: 1, 6: b},
                           'final_referee_T46_six_ray_coherent_capacity p.1')
    check('T46 six-ray: |P|^2 = p*, |C|^2 = c*', sp.simplify(nrm6(1) - pst) == 0 and sp.simplify(nrm6(5) - cst) == 0)
    schur6 = [sp.simplify((j - t ** 2) - (sp.Matrix([[t, v]]) * sp.Matrix([[q, be], [be, b]]).inv() * sp.Matrix([t, v]))[0, 0])]
    check('T46 six-ray: G0 PD iff q>0, l>0, j - t^2 > (bt^2 - 2beta t v + q v^2)/l',
          sp.simplify(schur6[0] - (j - t ** 2 - (b * t ** 2 - 2 * be * t * v + q * v ** 2) / ell)) == 0)
    stmt('G3-T46-six-ray', status='PROVED HERE (EXACT identity + PROOF-READ PSD)' if ok6 else 'NOT PROVED',
         generators='q f27, p f1, j f59, c f5, d f28 (nonminor), b f6',
         region='q,b>0, l>0, j-t^2 > (bt^2-2beta tv+qv^2)/l, d=1: p >= p*, c >= c*; finite limits by closedness',
         docs=['final_referee_T46_six_ray_coherent_capacity p.1-4', 'T46 a0_all_small_product Lemma 1',
               'root_moc4_T46_positive_a_mixed (1)', 'T47 zero_q_rational_capacity', 'T47 low_j_high_e', 'T47 low_j_low_e'])
    # ---- singular leaf beta = -r, v = -r t/q: explicit PSD vectors, no limit needed ----
    qp, bp, jp, tp = sp.symbols('q b j t', positive=True)
    rr = sp.sqrt(qp * bp)
    sub = {q: qp, b: bp, j: jp, be: -rr, v: -rr * tp / qp, t: tp}
    # explicit realisation of the singular leaf Gram: Q=(sqrt q,0), B=(-sqrt b,0), J0=(t/sqrt q, sqrt(j-t^2-t^2/q))
    Qv, Bv = sp.Matrix([sp.sqrt(qp), 0]), sp.Matrix([-sp.sqrt(bp), 0])
    J0v = sp.Matrix([tp / sp.sqrt(qp), sp.sqrt(jp - tp ** 2 - tp ** 2 / qp)])
    Gs = sp.Matrix([[(x.T * y)[0] for y in (Qv, Bv, J0v)] for x in (Qv, Bv, J0v)])
    check('T46 singular leaf: the 3x3 leaf/J0 block equals the Gram of explicit real vectors when j >= kappa t^2, kappa = 1+1/q',
          sp.simplify(Gs - G6[:3, :3].subs(sub)) == sp.zeros(3, 3))
    kap, rho = 1 + 1 / qp, 1 + rr / qp
    eta, kk = (1 - bp) / (qp - rr), qp * (1 - rr) ** 2 / (qp - rr) ** 2
    check('T46 singular leaf: p* = k + (1-eta t)^2/(j - kappa t^2) and c* = rho^2 t^2/(j - kappa t^2) (q != b)',
          sp.simplify(pst.subs(sub) - (kk + (1 - eta * tp) ** 2 / (jp - kap * tp ** 2))) == 0 and
          sp.simplify(cst.subs(sub) - rho ** 2 * tp ** 2 / (jp - kap * tp ** 2)) == 0)
    cp = sp.symbols('c', positive=True)
    D0 = kap + rho ** 2 / cp
    A0 = eta ** 2 * jp
    Nn = D0 ** 2 + A0 * (2 * D0 - kap + A0)
    ts = eta * jp / (D0 + A0)
    pS = kk + (1 - eta * ts) ** 2 / (jp - kap * ts ** 2)
    cS = rho ** 2 * ts ** 2 / (jp - kap * ts ** 2)
    check('T46 a0_all_small Lemma 1 / positive_a_mixed (1): at t = eta j/(D0+A0): j - kappa t^2 = jN/(D0+A0)^2, '
          'consumed p = k + D0^2/(jN), consumed c = rho^2 A0/N, N - (rho^2/c)A0 = D0^2 + D0 A0 + A0^2',
          sp.simplify(jp - kap * ts ** 2 - jp * Nn / (D0 + A0) ** 2) == 0 and sp.simplify(pS - (kk + D0 ** 2 / (jp * Nn))) == 0
          and sp.simplify(cS - rho ** 2 * A0 / Nn) == 0 and sp.simplify(Nn - rho ** 2 / cp * A0 - (D0 ** 2 + D0 * A0 + A0 ** 2)) == 0)
    stmt('G3-T46-singular-leaf', status='PROVED HERE (EXACT; specialisation of six-ray with explicit PSD vectors)',
         region='q != b, r = sqrt(qb), j > kappa t^2: p >= k + (1-eta t)^2/(j-kappa t^2), c >= rho^2 t^2/(j-kappa t^2); '
                'rational section p >= k + D0^2/(jN)',
         docs=['final_referee_T46_a0_all_small_product Lemma 1', 'root_moc4_T46_positive_a_mixed (1)',
               'final_referee_T46_six_ray Cor. 2-4 (T47 roles)'])
    # ---- four-coefficient capacity ab(e-1) >= a+b+e+1 on f4, f6, f2, f28 ----
    atom_ok('four', (4, 6, 2, 28), 'T46 positive_a_small_or_balanced, positive_a_mixed (2), whole_assembly (5)')
    stmt('G3-T46-four', status='REDUCES TO PROVED ATOM (four-coefficient (4,6;2;28))',
         docs=['final_referee_T46_positive_a_small_or_balanced p.1', 'root_moc4_T46_positive_a_mixed (2)', 'root_moc4_T46_whole_assembly (5)'])
    stmt('G3-T46-other', status='no universal SO(4) inequality stated',
         docs=['final_referee_T46_two_Heron_necessity', 'final_referee_T46_negative_j_Heron_uniform'])
    # six-ray -> T47 transport
    check('T46 six-ray transport [27,1,59,5,2,28,6,4] -> T47 [0,45,5,9,29,35,7,36] realised in G',
          bool(find_action([27, 1, 59, 5, 2, 28, 6, 4], [0, 45, 5, 9, 29, 35, 7, 36])))
    # six-ray Cor. 2-4 in T47 roles = singular leaf with (q,b,j,p,c)_46 = (b,a,e,c,p)_47: same formulas
    # (q_46 = b_47 -> kappa = 1 + 1/b, rho = 1 + r/b, eta = (1-a)/(b-r), k = b(1-r)^2/(b-r)^2): Cor. 2 choice t = eta e/kappa
    e47 = sp.symbols('e', positive=True)
    pS2 = kk + (1 - eta * tp) ** 2 / (jp - kap * tp ** 2)
    cS2 = rho ** 2 * tp ** 2 / (jp - kap * tp ** 2)
    t2 = eta * jp / kap
    check('six-ray Cor. 2: t = eta e/kappa gives consumed c = 1/e + k - eta^2/kappa and consumed p = rho^2 eta^2 e/(kappa(kappa - eta^2 e))',
          sp.simplify(pS2.subs(tp, t2) - (1 / jp + kk - eta ** 2 / kap)) == 0 and
          sp.simplify(cS2.subs(tp, t2) - rho ** 2 * eta ** 2 * jp / (kap * (kap - eta ** 2 * jp))) == 0)
    D0b = kap + rho ** 2 / cp
    t4 = eta * jp / D0b
    N0 = D0b ** 2 - kap * A0
    check('six-ray Cor. 4 (A<=D0): t = eta e/D0 gives consumed c = k + (D0-A)^2/(e(D0^2 - kappa A)), unused p = p D0(D0-A)/N0',
          sp.simplify(pS2.subs(tp, t4) - (kk + (D0b - A0) ** 2 / (jp * N0))) == 0 and
          sp.simplify(cp - cS2.subs(tp, t4) - cp * D0b * (D0b - A0) / N0) == 0)
    check('six-ray Cor. 4 (A>=D0): t = 1/eta gives consumed c = k and consumed p = rho^2/(A - kappa)',
          sp.simplify(pS2.subs(tp, 1 / eta) - kk) == 0 and sp.simplify(cS2.subs(tp, 1 / eta) - rho ** 2 / (A0 - kap)) == 0)
    stmt('G3-six-ray-Cor2-4', status='PROVED HERE (EXACT; T47 roles via G-transport)',
         docs=['final_referee_T46_six_ray_coherent_capacity p.2-4', 'root_moc4_T47_zero_q_rational_capacity', 'root_moc4_T47_low_j_high_e'])


# =============================================================================================
# G4. historical T47 (paper T_48, T_65), face q f2 + p f9 + j f1 + c f45 + e f5 + d f35 + b f0 + a f7
# =============================================================================================
if 'G4' in SECTIONS:
    print('==== G4 T47', flush=True)
    # the documents' "spectral" coordinates: find which of the supplement's congruences they use,
    # by testing the PCBE matrix (a 4-parameter family) against each candidate engine.
    xP, yP, zP, wP, LP, jP, cP, eP = sp.symbols('x y z w L j c e', real=True)

    def pcbe_C(nodes):
        lab = dict(zip(['Q', 'P', 'J', 'C', 'E', 'D', 'B', 'A'], [2, 9, 1, 45, 5, 35, 0, 7]))
        M = {}
        def put(a, b_, val):
            M[(lab[a], lab[b_])] = val; M[(lab[b_], lab[a])] = val
        for nm, val in zip('QPJCEDBA', (0, LP, jP, cP, eP, 1, LP, 0)):
            put(nm, nm, val)
        put('P', 'B', LP); put('P', 'J', yP); put('B', 'J', yP); put('J', 'D', -yP)
        put('P', 'C', -xP); put('B', 'C', -xP); put('C', 'E', -xP); put('J', 'C', 1 - xP + yP)
        put('P', 'E', zP); put('B', 'E', zP); put('E', 'D', -zP); put('J', 'E', wP)
        return [[M.get((f, g), 0) for g in nodes] for f in nodes]
    BJ = json.load(open(os.path.join(BASE, 'independent_moc4_T47_budget_PSD8_20261022.json')))
    NJ = json.load(open(os.path.join(BASE, 'independent_moc4_T47_native_PSD8_20261022.json')))
    Tb = sp.Matrix([[sp.Rational(v) for v in r] for r in BJ['T']])
    Vn = sp.Matrix([[sp.Rational(v) for v in r] for r in NJ['V']])
    cands = [('V T^-1', Vn * Tb.inv()), ('V^-1 T^-1', Vn.inv() * Tb.inv()), ('V^T T^-1', Vn.T * Tb.inv()),
             ('V^-T T^-1', Vn.inv().T * Tb.inv()), ('T^-1 V', Tb.inv() * Vn), ('T^-1 V^-1', Tb.inv() * Vn.inv()),
             ('budget T^-1', Tb.inv())]
    F47 = None
    for lab, Tio in cands:
        Fc = Face('T47[' + lab + ']', 'independent_moc4_T47_cubic_reduction_20261022.json',
                  'independent_moc4_T47_budget_PSD8_20261022.json', Ti_override=Tio)
        Fc.nodes = list(Fc.src)
        bad = Fc.residual(pcbe_C(Fc.nodes), [0 if f in (2, 7) else {9: LP, 1: jP, 45: cP, 5: eP, 35: 1, 0: LP}[f] for f in Fc.src])
        print(f'INFO T47 coordinates {lab}: PCBE matrix residual count {len(bad)}', flush=True)
        if not bad and F47 is None:
            F47 = Fc
            print('INFO T47 spectral coordinates identified as', lab, flush=True)
    FACES['T47'] = F47
    check('T47 PCBE: spectral coordinates (budget congruence composed with native V^-1) found in which the PCBE matrix satisfies every coefficient relation identically in x, y, z, w, L',
          F47 is not None)
    LAB47 = dict(zip('QPJCEDBA', [2, 9, 1, 45, 5, 35, 0, 7]))

    def C47(face, ent):
        """ent: {('J','E'): value, ...} symmetric; returns 8x8 list in face.nodes order."""
        M = {}
        for (u, v), val in ent.items():
            M[(LAB47[u], LAB47[v])] = sp.sympify(val); M[(LAB47[v], LAB47[u])] = sp.sympify(val)
        return [[M.get((f, g), sp.Integer(0)) for g in face.nodes] for f in face.nodes]

    def w47(face, **kw):
        nm = dict(zip([2, 9, 1, 45, 5, 35, 0, 7], 'qpjcedba'))
        return [sp.sympify(kw.get(nm[f], 0)) for f in face.src]

    if F47 is not None:
        # ---- PCBE Theorem 1: Schur complement, scalar identity (3), simplification to (1), crossing x*
        M8 = sp.Matrix(C47(F47, {('P', 'P'): LP, ('J', 'J'): jP, ('C', 'C'): cP, ('E', 'E'): eP, ('D', 'D'): 1, ('B', 'B'): LP,
                                  ('P', 'B'): LP, ('P', 'J'): yP, ('B', 'J'): yP, ('J', 'D'): -yP, ('P', 'C'): -xP,
                                  ('B', 'C'): -xP, ('C', 'E'): -xP, ('J', 'C'): 1 - xP + yP, ('P', 'E'): zP,
                                  ('B', 'E'): zP, ('E', 'D'): -zP, ('J', 'E'): wP}))
        ix = {k: F47.nodes.index(v) for k, v in LAB47.items()}
        okPB = all(sp.expand(M8[ix['P'], k] - M8[ix['B'], k]) == 0 for k in range(8))
        BD = [ix['B'], ix['D']]; JEC = [ix['J'], ix['E'], ix['C']]
        S = sp.simplify(M8.extract(JEC, JEC) - M8.extract(JEC, BD) * M8.extract(BD, BD).inv() * M8.extract(BD, JEC))
        kap = 1 + 1 / LP
        S_doc = sp.Matrix([[jP - kap * yP ** 2, wP - kap * yP * zP, 1 - xP + (1 + xP / LP) * yP],
                           [wP - kap * yP * zP, eP - kap * zP ** 2, -xP + xP * zP / LP],
                           [1 - xP + (1 + xP / LP) * yP, -xP + xP * zP / LP, cP - xP ** 2 / LP]])
        check('T47 PCBE: rows P and B coincide (P = B) and the Schur complement over {B, D} is the displayed 3x3 matrix on J, E, C',
              okPB and sp.simplify(S - S_doc) == sp.zeros(3, 3))
        A_, B_, J_, k_, t_ = sp.symbols('A B J kappa t', positive=True)
        check('T47 PCBE (3): (A+Bt)^2/(J-k t^2) - (A^2/J - B^2/k) = (kAt+BJ)^2/(kJ(J-k t^2))',
              sp.simplify((A_ + B_ * t_) ** 2 / (J_ - k_ * t_ ** 2) - (A_ ** 2 / J_ - B_ ** 2 / k_)
                          - (k_ * A_ * t_ + B_ * J_) ** 2 / (k_ * J_ * (J_ - k_ * t_ ** 2))) == 0)
        Kq = 1 / jP + 1 / (LP + 1); Bq1 = 1 / eP + 1 / (LP + 1)
        e1 = sp.simplify(xP ** 2 / LP + (1 - xP) ** 2 / jP - (1 + xP / LP) ** 2 / kap - (Kq * (1 - xP) ** 2 - 1))
        e2 = sp.simplify(xP ** 2 / LP + xP ** 2 * (1 / eP - 1 / (LP * (LP + 1))) - Bq1 * xP ** 2)
        check('T47 PCBE: x^2/L + [(1-x)^2/j - (1+x/L)^2/kappa] = K(1-x)^2 - 1 and x^2/L + x^2(1/e - 1/(L(L+1))) = (1/e+1/(L+1)) x^2',
              e1 == 0 and e2 == 0)
        Ks, Bs = sp.symbols('K Bc', positive=True)
        xs = (Ks - 1) / (Ks + sp.sqrt(Ks + Bs * (Ks - 1)))
        check('T47 PCBE: x* = (K-1)/(K+sqrt(K+B(K-1))) solves K(1-x)^2 - 1 = B x^2',
              sp.simplify(sp.expand(Ks * (1 - xs) ** 2 - 1 - Bs * xs ** 2)) == 0
              or abs(float((Ks * (1 - xs) ** 2 - 1 - Bs * xs ** 2).subs({Ks: 2.7, Bs: 1.3}))) < 1e-14)
        # (1) <=> (2): c >= min_x max{K(1-x)^2-1, Bx^2}  <=>  sqrt((c+1)/K) + sqrt(c/B) >= 1  (sampled exactly-rational check of both sides)
        okeq = True
        rng = np.random.default_rng(47)
        for _ in range(400):
            K_, B__, c_ = [float(v) for v in rng.uniform(0.05, 5, 3)]
            if K_ <= 1:
                lhs = True
            else:
                x_ = (K_ - 1) / (K_ + math.sqrt(K_ + B__ * (K_ - 1)))
                lhs = c_ >= B__ * x_ ** 2 - 1e-12
            rhs = math.sqrt((c_ + 1) / K_) + math.sqrt(c_ / B__) >= 1 - 1e-12
            okeq &= (lhs == rhs)
        check('T47 PCBE (1) <=> (2) on 400 random (K, B, c) [NUMERIC consistency, not exact]', okeq)
        stmt('G4-PCBE', status='PROVED IN THE DOCUMENT AND CHECKED', generators=[2, 9, 1, 45, 5, 35, 0, 7],
             region='d=1, p,b,j,e>0, 0<L<=min(p,b); c >= min_x max{K(1-x)^2-1, Bx^2}',
             docs=['T47 PCBE_energy p.1-2', 'T47 PCBE_low_e_q_monotone', 'T47 whole_assembly'],
             evidence='EXACT Gram relations + Schur + (3) + simplification; PSD completion and closure PROOF-READ')
        # ---- IIIb_plus_large_pc Lemma 1 (try spectral coordinates, then budget coordinates)
        a_, p_, c_, d_ = sp.symbols('a p c d', positive=True)
        ent = {('P', 'P'): p_, ('C', 'C'): c_, ('A', 'A'): a_, ('P', 'C'): -a_, ('D', 'D'): d_,
               ('J', 'J'): d_ ** 2 / (a_ + c_), ('J', 'P'): -d_ * a_ / (a_ + c_), ('J', 'A'): -d_ * a_ / (a_ + c_),
               ('J', 'C'): d_ * c_ / (a_ + c_)}
        wl = dict(p=p_, j=d_ ** 2 / (a_ + c_), c=c_, d=d_, a=a_)
        Fb = Face('T47[budget]', 'independent_moc4_T47_cubic_reduction_20261022.json',
                  'independent_moc4_T47_budget_PSD8_20261022.json')
        Fb.nodes = list(Fb.src)
        r_sp = F47.residual(C47(F47, ent), w47(F47, **wl))
        r_bu = Fb.residual(C47(Fb, ent), w47(Fb, **wl))
        print(f'INFO T47 large_pc Gram residuals: spectral {len(r_sp)}, budget {len(r_bu)}', flush=True)
        # vectors: P=(sp,0,0,0), C=(-a/sp, sqrt(c-a^2/p),0,0), A=(0,0,sa,0), D=(0,0,0,sd), J=d(C-A)/(a+c)
        Pv = sp.Matrix([sp.sqrt(p_), 0, 0, 0]); Cv = sp.Matrix([-a_ / sp.sqrt(p_), sp.sqrt(c_ - a_ ** 2 / p_), 0, 0])
        Av = sp.Matrix([0, 0, sp.sqrt(a_), 0]); Dv = sp.Matrix([0, 0, 0, sp.sqrt(d_)]); Jv = d_ * (Cv - Av) / (a_ + c_)
        V = {'P': Pv, 'C': Cv, 'A': Av, 'D': Dv, 'J': Jv}
        okv = all(sp.simplify((V[u].T * V[v])[0, 0] - val) == 0 for (u, v), val in ent.items())
        okv &= all(sp.simplify((V[u].T * V[v])[0, 0]) == 0 for u, v in itertools.combinations('PCADJ', 2)
                   if (u, v) not in ent and (v, u) not in ent)
        check('T47 IIIb_plus_large_pc Lemma 1: Gram (A perp P,C; <P,C>=-a; J=d(C-A)/(a+c); D perp) satisfies every coefficient relation, '
              'and is realised by explicit vectors when pc >= a^2', (not r_sp or not r_bu) and okv,
              'coordinates: ' + ('spectral' if not r_sp else 'budget' if not r_bu else 'neither'))
        stmt('G4-largepc-L1', status='PROVED IN THE DOCUMENT AND CHECKED', generators=[2, 9, 1, 45, 5, 35, 0, 7],
             region='a,p,c,d>0, pc>=a^2, j>=d^2/(a+c)', docs=['T47 IIIb_plus_large_pc p.1', 'T47 merged_C Thm 2'],
             evidence='EXACT')
        # ---- merged-C Lemma 1: crossing, stationary, half-energy Grams for both leaves
        jj, ee, LL, dd, xx = sp.symbols('jr er L d x', positive=True)   # j = jr^2, e = er^2
        okm = []
        for leaf, sg in (('A', 1), ('B', -1)):
            lw = 'a' if leaf == 'A' else 'b'
            xc = dd * ee / (jj + ee)
            grams = {
                'crossing': ({('J', 'J'): jj ** 2, ('E', 'E'): ee ** 2, ('J', 'E'): -jj * ee, ('D', 'D'): dd, (leaf, leaf): LL,
                              ('C', 'J'): dd - xc, ('C', 'E'): -xc, ('C', leaf): sg * xc,
                              ('C', 'C'): xc ** 2 / LL + xc ** 2 / ee ** 2},
                             dict(j=jj ** 2, e=ee ** 2, d=dd, c=xc ** 2 / LL + xc ** 2 / ee ** 2, **{lw: LL})),
                'stationary': ({('J', 'J'): jj ** 2, ('E', 'E'): ee ** 2, ('J', 'E'): -LL, ('D', 'D'): dd, (leaf, leaf): LL,
                                ('C', 'J'): dd * jj ** 2 / (jj ** 2 + LL), ('C', 'E'): -dd * LL / (jj ** 2 + LL),
                                ('C', leaf): sg * dd * LL / (jj ** 2 + LL), ('C', 'C'): dd ** 2 / (jj ** 2 + LL)},
                               dict(j=jj ** 2, e=ee ** 2, d=dd, c=dd ** 2 / (jj ** 2 + LL), **{lw: LL})),
                'half-energy': ({('J', 'J'): jj ** 2, ('E', 'E'): ee ** 2, ('J', 'E'): -jj ** 2, ('D', 'D'): dd, (leaf, leaf): LL,
                                 ('C', 'J'): dd / 2, ('C', 'E'): -dd / 2, ('C', leaf): sg * dd / 2,
                                 ('C', 'C'): dd ** 2 / (4 * jj ** 2) + dd ** 2 / (4 * LL)},
                                dict(j=jj ** 2, e=ee ** 2, d=dd, c=dd ** 2 / (4 * jj ** 2) + dd ** 2 / (4 * LL), **{lw: LL}))}
            for nm, (ent, wk) in grams.items():
                rs = F47.residual(C47(F47, ent), w47(F47, **wk))
                okm.append((leaf, nm, len(rs)))
        print('INFO T47 merged-C residual counts (leaf, gram, count):', okm, flush=True)
        check('T47 merged-C Lemma 1: crossing, stationary and half-energy Grams for leaf a and leaf b satisfy every coefficient relation',
              all(n == 0 for _, _, n in okm))
        # PSD of the three Grams: crossing has E = -(er/jr) J; stationary needs je >= L^2; half-energy needs e >= j (vectors written in the md)
        xcr = dd * ee / (jj + ee)
        check('T47 merged-C: crossing value x^2/L + x^2/e at x = d sqrt(e)/(sqrt j + sqrt e) equals d^2 (L+e)/(L (sqrt j + sqrt e)^2)',
              sp.simplify(xcr ** 2 / LL + xcr ** 2 / ee ** 2 - dd ** 2 * (LL + ee ** 2) / (LL * (jj + ee) ** 2)) == 0)
        dl = sp.symbols('delta', positive=True)
        fcr = lambda d0, j0, e0, L0: d0 ** 2 * (L0 + e0) / (L0 * (sp.sqrt(j0) + sp.sqrt(e0)) ** 2)
        J0, E0, A0 = sp.symbols('j e A0', positive=True)
        check('T47 merged-C: the three bounds are homogeneous of degree one in (d, j, e, L), so the proportional split sums to the A0 = a+b bound',
              sp.simplify(fcr(dl * dd, dl * J0, dl * E0, dl * A0) - dl * fcr(dd, J0, E0, A0)) == 0
              and sp.simplify((dl * dd) ** 2 / (dl * J0 + dl * A0) - dl * dd ** 2 / (J0 + A0)) == 0)
        stmt('G4-mergedC-L1', status='PROVED IN THE DOCUMENT AND CHECKED', generators=[1, 45, 5, 35, 0, 7],
             region='A0=a+b>0, d>0; crossing (je<=A0^2), stationary (je>=A0^2), half-energy (e>=j)',
             docs=['T47 merged_C p.1', 'T47 IIIb_plus_unswapped', 'T47 whole_assembly'], evidence='EXACT')
    # ---- orbit atoms and transports used by T47 documents
    for fam, t in (('H', (9, 1, 35)), ('SOC', (45, 1, 35)), ('chain', (1, 7, 9, 35)), ('four', (0, 7, 5, 35))):
        atom_ok(fam, t, 'T47')
    check('T47 transport T24 core (59,53,55,67,54,2) -> (45,5,7,0,35,1) is realised by some g in G',
          bool(find_action([59, 53, 55, 67, 54, 2], [45, 5, 7, 0, 35, 1])))
    stmt('G4-atoms', status='REDUCES TO PROVED ATOMS', docs=['T47 I_II', 'T47 high_j_uniform', 'T47 high_j_small_b'],
         generators=[9, 1, 45, 7, 0, 5, 35], evidence='EXACT orbit membership')


# =============================================================================================
# G2. historical T45: contacts {3,4,5,9,10,12}, source p f0 + r f1 + a f2 + b f12 + c f59 + q f67 + d f70 + s f72
# =============================================================================================
if 'G2' in SECTIONS:
    print('==== G2 T45', flush=True)
    T45 = [3, 4, 5, 9, 10, 12]
    NS = L.nullspace([L.P[s] for s in T45])
    gen45 = [0, 1, 2, 12, 59, 67, 70, 72]
    ev = [18, 16, 13, 8, 2, 14, 11, 7]
    check('T45 contact nullspace has dimension 8', len(NS) == 8)
    check('T45 the eight sources vanish on the contacts', all(L.F_VALUES[f][s] == 0 for f in gen45 for s in T45))
    check('T45 evaluations of (f0,f1,f2,f12,f59,f67,f70,f72) at (18,16,13,8,2,14,11,7) form the identity matrix',
          [[L.F_VALUES[f][s] for f in gen45] for s in ev] == [[int(i == k) for k in range(8)] for i in range(8)])
    ok1 = True
    for c in NS:
        g = [Fr(x) for x in L.values(c)]
        ok1 &= fvals({f: g[s] for f, s in zip(gen45, ev)}) == g
    check('T45 (1): every contact-nullspace cost equals sum_f g_(ev f) f on the whole nullspace', ok1)
    check('T45 the seven sources other than f70 are minors', all(L.F_IS_MINOR[f] for f in gen45 if f != 70)
          and not L.F_IS_MINOR[70])
    for fam, t in (('SOC', (2, 59, 70)), ('chain', (2, 0, 12, 70)), ('chain', (59, 0, 1, 70)),
                   ('chain', (59, 67, 72, 70)), ('H', (1, 59, 70)), ('H', (72, 59, 70))):
        atom_ok(fam, t, 'T45 (2)-(6)')
    r_, s_, d_ = sp.symbols('r s d', positive=True)
    R_ = r_ + s_
    d1, d2 = d_ * r_ / R_, d_ * s_ / R_
    alloc_check('T45 (6) old-H split c >= d + d^2/(r+s)',
                [('H', (1, 59, 70), (r_, d1 + d1 ** 2 / r_, d1)), ('H', (72, 59, 70), (s_, d2 + d2 ** 2 / s_, d2))],
                {1: r_, 72: s_, 59: d_ + d_ ** 2 / R_, 70: d_}, 'T45')
    lit = (bool(literal_action_matches((0, 3, 1, 2), (3, 1, 0, 2), 0, [3, 15, 26], [1, 59, 70]))
           and bool(literal_action_matches((1, 2, 0, 3), (0, 2, 3, 1), 0, [3, 15, 26], [72, 59, 70])))
    print(f'INFO T45 old-H literal (row,col) tables reproduced by one of our readings: {lit} '
          f'(orbit membership above is what the proof needs)', flush=True)
    for nm, src, img in (('a=0: T33', [2, 3, 5, 9, 10, 12], [3, 5, 9, 10, 12, 13]),
                         ('c=0: T6', [0, 1, 2, 4, 7, 8], [2, 3, 4, 5, 9, 10]),
                         ('r=0: T14', [1, 2, 3, 4, 6, 7], [3, 4, 5, 9, 10, 16])):
        check(f'T45 zero face {nm} contact set is carried by some g in G onto {img}',
              bool(find_action([], [], contact=(src, img))))
    stmt('G2-T45', status='REDUCES TO PROVED ATOMS', generators=gen45,
         region='p,r,a,b,c,q,d,s >= 0; capacities (2)-(6)', docs=['T45 uniform p.1-2'],
         evidence='EXACT (nullspace identity, orbit membership, threshold allocation)')


# =============================================================================================
# NUM. boundary attacks (NUMERIC evidence only): coefficients exactly on the stated threshold,
# minimise Phi over SO(4) (random sample + BFGS), re-evaluate the worst point with mpmath (40 digits)
# =============================================================================================
if 'NUM' in SECTIONS:
    print('==== NUM boundary attacks', flush=True)
    NS_ = 20
    U = lambda rng, lo=0.05, hi=4.0: float(np.exp(rng.uniform(np.log(lo), np.log(hi))))

    def g_HC(rng):          # T42 roles (aux b f3, aux c f23, central f5, leaves f0, f6, nonminor f54), delta minimised, beta grid
        b_, c_, p_, q_ = U(rng), U(rng), U(rng), U(rng)
        best = np.inf
        for b1 in np.linspace(0, 1, 401):
            d1, d2 = min(p_ / 2, b1), min(q_ / 2, 1 - b1)
            best = min(best, 1 - d1 - d2 + b1 ** 2 / b_ + (1 - b1) ** 2 / c_ + d1 ** 2 / p_ + d2 ** 2 / q_)
        return {3: b_, 23: c_, 5: best, 0: p_, 6: q_, 54: 1.0}

    def g_T44P1(rng):
        r, k, v, w = rng.normal(size=4) * 1.5
        if abs(r - k) < 0.05 or abs(w) < 0.05:
            return None
        x = (v * (2 * r - k) - 1) / (r - k); y = k * (1 - r * v) / (w * (r - k))
        return {1: x * x + y * y, 53: r * r * (1 + v * v), 12: k * k * (1 + v * v), 8: v * v, 15: v * v, 7: w * w, 18: 1.0}

    def g_T44P2(rng):
        b_, j_, c_ = U(rng), U(rng), U(rng)
        p_ = (b_ + c_) / (b_ * (math.sqrt(j_) + math.sqrt(c_)) ** 2) if j_ * c_ <= b_ * b_ else 1 / (j_ + b_)
        return {1: p_, 53: j_, 12: c_, 7: b_, 18: 1.0}

    def g_T46core(rng):
        q_, b_ = U(rng), U(rng)
        be_ = rng.uniform(-1, 1) * math.sqrt(q_ * b_) * 0.98
        l_ = q_ * b_ - be_ ** 2
        t_ = rng.normal()
        j_ = (1 + (q_ + b_ - 2 * be_) / l_) * t_ ** 2 + U(rng, 0.01, 3)
        s_, h_ = q_ + b_ + 2 * be_, 1 + be_
        Dl = s_ * j_ - (s_ + 4) * t_ ** 2
        if Dl <= 1e-9:
            return None
        p_ = ((j_ - t_ ** 2) * h_ ** 2 - 4 * t_ * h_ * (t_ - 1) + s_ * (t_ - 1) ** 2) / Dl
        return {27: q_, 1: p_, 59: j_, 28: 1.0, 6: b_}

    def g_T46six(rng):
        q_, b_ = U(rng), U(rng)
        be_ = rng.uniform(-1, 1) * math.sqrt(q_ * b_) * 0.98
        l_ = q_ * b_ - be_ ** 2
        t_, v_ = rng.normal(size=2)
        j_ = t_ ** 2 + (b_ * t_ ** 2 - 2 * be_ * t_ * v_ + q_ * v_ ** 2) / l_ + U(rng, 0.01, 3)
        s_, h_ = q_ + b_ + 2 * be_, 1 + be_
        D6 = s_ * (j_ - t_ ** 2) - (t_ + v_) ** 2
        den = j_ - t_ ** 2 - v_ ** 2 / b_
        if D6 <= 1e-9 or den <= 1e-9:
            return None
        p_ = ((j_ - t_ ** 2) * h_ ** 2 - 2 * h_ * (t_ + v_) * (1 + v_) + s_ * (1 + v_) ** 2) / D6
        return {27: q_, 1: p_, 59: j_, 5: (t_ - v_) ** 2 / den, 28: 1.0, 6: b_}

    def g_T45H(rng):
        r_, s_ = U(rng), U(rng)
        return {1: r_, 72: s_, 59: 1 + 1 / (r_ + s_), 70: 1.0}

    def g_PCBE(rng):
        L_, j_, e_ = U(rng), U(rng), U(rng)
        K_ = 1 / j_ + 1 / (L_ + 1); B__ = max(1 / e_ + 1 / (L_ + 1), 1 / L_)
        c_ = 0.0 if K_ <= 1 else B__ * ((K_ - 1) / (K_ + math.sqrt(K_ + B__ * (K_ - 1)))) ** 2
        return {9: L_, 0: L_, 1: j_, 45: c_, 5: e_, 35: 1.0}

    def g_largepc(rng):
        a_, c_ = U(rng), U(rng)
        return {9: a_ * a_ / c_, 1: 1 / (a_ + c_), 45: c_, 35: 1.0, 7: a_}

    def g_mergedC(rng):
        a_, b_, j_, e_ = U(rng), U(rng), U(rng), U(rng)
        A0 = a_ + b_
        c_ = (A0 + e_) / (A0 * (math.sqrt(j_) + math.sqrt(e_)) ** 2) if j_ * e_ <= A0 * A0 else 1 / (j_ + A0)
        return {7: a_, 0: b_, 1: j_, 5: e_, 35: 1.0, 45: c_}

    ATT = [('G5 H/chain allocation (6), T42 roles, minimised RHS', g_HC, 'G5-HC-allocation'),
           ('G1 T44 Prop. 1 family', g_T44P1, 'G1-T44-P1'), ('G1 T44 Prop. 2 threshold', g_T44P2, 'G1-T44-P2'),
           ('G2 T45 (6) old-H split', g_T45H, 'G2-T45'),
           ('G3 T46 core threshold (1)', g_T46core, 'G3-T46-core'), ('G3 T46 six-ray p*, c*', g_T46six, 'G3-T46-six-ray'),
           ('G4 T47 PCBE threshold (1), L = p = b', g_PCBE, 'G4-PCBE'), ('G4 T47 large_pc Lemma 1 at pc = a^2', g_largepc, 'G4-largepc-L1'),
           ('G4 T47 merged-C Lemma 1 (crossing / stationary)', g_mergedC, 'G4-mergedC-L1')]
    for i, (nm, gen, key) in enumerate(ATT):
        mn = attack(nm, gen, NS_, 1000 + i)
        stmt(key, numeric_min_normalised=mn, numeric_samples=NS_)
        check(f'NUMERIC {nm}: no violation beyond -1e-9 (evidence only)', mn > -1e-9, f'{mn:.3e}')

    # sensitivity controls.  The face thresholds above are sufficient conditions and need not be sharp
    # (e.g. T44 Prop. 2 at b=j=c=1 survives p lowered by 50%), so the controls use the sharp old-H atom
    # lowered by 10%, and the faces with their main minor removed (factor 0).
    def lowered(gen, key, fac):
        def g(rng):
            c = gen(rng)
            if c is not None:
                c = dict(c); c[key] *= fac
            return c
        return g
    def g_oldH(rng):
        r_ = U(rng, 0.3, 3)
        return {1: r_, 59: 1 + 1 / r_, 70: 1.0}
    for i, (nm, gen) in enumerate([('control old-H H(1,59;70) with main coefficient lowered 10%', lowered(g_oldH, 59, 0.9)),
                                   ('control T44 Prop. 2 with p set to 0', lowered(g_T44P2, 1, 0.0)),
                                   ('control T46 core with p set to 0', lowered(g_T46core, 1, 0.0)),
                                   ('control T47 large_pc with j set to 0', lowered(g_largepc, 1, 0.0))]):
        mn = attack(nm, gen, 6, 2000 + i)
        check(f'NUMERIC {nm}: attack finds a violation (sensitivity control)', mn < -1e-8, f'{mn:.3e}')

# =============================================================================================
OUT['numeric'] = NUMREC
OUT['sections'] = sorted(SECTIONS)
OUT['n_fail'] = NFAIL
print(f'==== done: {len(OUT["exact"])} checks, {NFAIL} FAIL, {time.time() - T0:.1f} s', flush=True)
if os.environ.get('MOC4_WRITE_JSON'):
    with open(os.path.join(HERE, 'so4_atoms_part2b.json'), 'w') as fh:
        json.dump(OUT, fh, indent=1, default=str)
