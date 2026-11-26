"""Second-pass exact replay for final_referee_T46_core_coherent_capacity_20261023.

Usage: python -I -B this.py [SUPP]   (or set MOC4_SUPP to the supplement root)

What is our own: the SO(4) parametrization O(p,q): v -> p v q (unit quaternions), every
Leibniz sum Phi_f (from the lib's Appendix D values), the Gram matrix of Theorem 1 (built
from the PDF text), the reduction modulo |p|^2 - |q|^2, and all linear algebra.
What is taken from the supplement as DATA only (json.load, parsed by our own parser, never
executed): the 16 cubic polynomials ("features"), the 4 signed coordinate permutations
("left_four") and the 8x8 change of basis T.  These are only candidate objects; the
identities below decide whether they work.  The certificate is valid because of the
identities, whatever their origin.

Prints one line per check; exits 1 if an EXACT check fails.
"""
import itertools
import json
import os
import re
import sys
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "lib"))
import moc4lib as L  # noqa: E402
import flint  # noqa: E402
import sympy as sp  # noqa: E402

SUPP = Path(sys.argv[1] if len(sys.argv) > 1 else os.environ["MOC4_SUPP"])
DATA = SUPP / "proofs" / "manin_transfer_20261021"
FAIL = []


def report(name, ok, level="EXACT", info=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {level} {name} {info}".rstrip())
    if not ok and level == "EXACT":
        FAIL.append(name)


# ---------------------------------------------------------------------------------
# 0. safe polynomial parser for the feature strings (no eval)
VARS = ['p0', 'p1', 'p2', 'p3', 'q0', 'q1', 'q2', 'q3']


def parse_poly(s):
    s = s.replace(' ', '')
    if s[0] not in '+-':
        s = '+' + s
    out = {}
    for t in re.findall(r'[+-][^+-]+', s):
        body = t[1:]
        assert re.sub(r'[a-z]\d(\*\*\d+)?|\d+(/\d+)?', '', body).replace('*', '') == '', t
        coef = Fraction(-1 if t[0] == '-' else 1)
        ex = [0] * 8
        for v, e, num in re.findall(r'([a-z]\d)(?:\*\*(\d+))?|(\d+(?:/\d+)?)', body):
            if v:
                ex[VARS.index(v)] += int(e) if e else 1
            else:
                coef *= Fraction(num)
        out[tuple(ex)] = out.get(tuple(ex), 0) + coef
    return {k: v for k, v in out.items() if v}


ctx = flint.fmpq_mpoly_ctx.get(tuple(VARS), 'degrevlex')
X = ctx.gens()
ZERO = ctx.from_dict({})


def from_dict(d):
    return ctx.from_dict({k: flint.fmpq(v.numerator, v.denominator) for k, v in d.items()})


# ---------------------------------------------------------------------------------
# 1. SO(4) parametrization O(p,q) v = p v q
def qmul(a, b):
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    return [a0*b0 - a1*b1 - a2*b2 - a3*b3, a0*b1 + a1*b0 + a2*b3 - a3*b2,
            a0*b2 - a1*b3 + a2*b0 + a3*b1, a0*b3 + a1*b2 - a2*b1 + a3*b0]


P4, Q4 = X[:4], X[4:]
E = [[ctx.from_dict({(0,)*8: 1}) if i == j else ZERO for i in range(4)] for j in range(4)]
cols = [qmul(qmul(P4, E[j]), Q4) for j in range(4)]
O = [[cols[j][i] for j in range(4)] for i in range(4)]       # O[i][j] = (p e_j q)_i
NP = sum(x*x for x in P4)
NQ = sum(x*x for x in Q4)
ok = all((sum(O[i][k]*O[j][k] for k in range(4)) - (NP*NQ if i == j else ZERO)).is_zero()
         for i in range(4) for j in range(4))
report("O(p,q) O(p,q)^T = |p|^2|q|^2 I (so O in O(4) on the unit torus)", ok)
# det O = (|p|^2|q|^2)^2 (so det = +1 on unit locus)
detO = sum(L.SGN[k]*O[0][s[0]]*O[1][s[1]]*O[2][s[2]]*O[3][s[3]] for k, s in enumerate(L.PERMS))
report("det O(p,q) = (|p|^2|q|^2)^2", (detO - (NP*NQ)**2).is_zero())

LW = [L.SGN[k]*O[0][s[0]]*O[1][s[1]]*O[2][s[2]]*O[3][s[3]] for k, s in enumerate(L.PERMS)]


def Phi_hom(vals):
    out = ZERO
    for v, w in zip(vals, LW):
        if v:
            out += int(v)*w
    return out


# ---------------------------------------------------------------------------------
# 2. normal form modulo (|p|^2 - |q|^2): replace p0^2 by |q|^2 - p1^2 - p2^2 - p3^2
NN = NQ - X[1]**2 - X[2]**2 - X[3]**2
NNpow = [ctx.from_dict({(0,)*8: 1})]
for _ in range(8):
    NNpow.append(NNpow[-1]*NN)


def NF(f):
    parts = {}
    for mon, c in f.to_dict().items():
        k = mon[0]
        parts.setdefault(k, {})[(0,) + tuple(mon[1:])] = c
    out = ZERO
    for k, d in parts.items():
        g = ctx.from_dict(d)
        out += (X[0] if k % 2 else 1) * NNpow[k // 2] * g
    return out


# ---------------------------------------------------------------------------------
# 3. candidate cubic features (data) and the averaged Gram functions
cub = json.load(open(DATA / "final_referee_T46_cubic_reduction_20261022.json"))
bud = json.load(open(DATA / "final_referee_T46_budget_PSD8_20261022.json"))
feats = [[from_dict(parse_poly(f)) for f in blk] for blk in cub['features']]
group = cub['left_four']
NATIVE = [27, 1, 59, 5, 2, 28, 6, 4]           # roles Q,P,J,C,E,D,B,A
assert bud['source_order'] == NATIVE
T = sp.Matrix([[sp.Rational(x) for x in row] for row in bud['T']])
Ti = T.inv()
u_all = []
for blk in feats:
    for g in group:
        subs = [s*X[i] for i, s in g]
        u = [f.compose(*subs) for f in blk]
        # native cubic vector psi = Ti u
        psi = []
        for k in range(8):
            acc = ZERO
            for m in range(8):
                if Ti[k, m] != 0:
                    acc += flint.fmpq(int(Ti[k, m].p), int(Ti[k, m].q)) * u[m]
            psi.append(acc)
        u_all.append(psi)
ENT = [(k, l) for k in range(8) for l in range(k, 8)]
W = []
for (k, l) in ENT:
    acc = ZERO
    for psi in u_all:
        acc += psi[k]*psi[l]
    W.append(acc if k == l else 2*acc)
# identity sought on |p|=|q|=1:  sum_e G_e W_e = 2 Phi_c ; homogeneous form: |p|^2 W vs 2 Phi_hom
cols_G = [NF(NP*w) for w in W]
cols_c = [NF(2*Phi_hom(L.F_VALUES[r])) for r in NATIVE]
allcols = cols_G + [-c for c in cols_c]
mons = sorted({m for c in allcols for m in c.to_dict()})
idx = {m: i for i, m in enumerate(mons)}
Mrows = [[flint.fmpq(0)]*len(allcols) for _ in mons]
for j, c in enumerate(allcols):
    for m, v in c.to_dict().items():
        Mrows[idx[m]][j] = v
M = flint.fmpq_mat(Mrows)
R, rk = M.rref()
print(f"[INFO] coefficient rows {len(mons)}, unknowns {len(allcols)} (36 Gram + 8 cost), rank {rk}")
Rr = [[sp.Rational(int(R[i, j].p), int(R[i, j].q)) for j in range(len(allcols))] for i in range(rk)]

# ---------------------------------------------------------------------------------
# 4. derived completion relations and cost map, in native coordinates
gsym = sp.symbols('g0:36')
csym = sp.symbols('c0:8')
vec = list(gsym) + list(csym)
Rm = sp.Matrix(Rr)
# eliminate cost unknowns: rref has the cost columns; relations on G = rows whose cost part is zero
A = Rm[:, :36]
B = Rm[:, 36:]
Nl = B.T.nullspace()
rel = sp.Matrix.vstack(*[(n.T*A) for n in Nl]) if Nl else sp.zeros(0, 36)
relr = rel.rref()[0]
relrank = rel.rank()
print(f"[INFO] number of independent completion relations on G: {relrank}")
# cost map: solve B c = -A g for c (B has full column rank 8?)
report("cost columns independent (cost of a Gram is unique)", B.rank() == 8)
Gm = sp.zeros(8)
for e, (k, l) in enumerate(ENT):
    Gm[k, l] = Gm[l, k] = gsym[e]
# compare with the doc's budget diagonal on the allowed set:
# c27=G00, c1=G11, c59=G22, c5=G33, c2=G44, c28=G55-G33, c6=G66, c4=G77
budget = [Gm[0, 0], Gm[1, 1], Gm[2, 2], Gm[3, 3], Gm[4, 4], Gm[5, 5] - Gm[3, 3], Gm[6, 6], Gm[7, 7]]
# (G, budget(G)) must satisfy all rref equations whenever G satisfies the relations:
# i.e. each row of R applied to (G, budget(G)) lies in the row span of rel.
okb = True
for i in range(rk):
    row = sp.expand(sum(Rm[i, j]*gsym[j] for j in range(36)) + sum(Rm[i, 36+m]*budget[m] for m in range(8)))
    coeffs = sp.Matrix([[row.coeff(gsym[j]) for j in range(36)]])
    if relrank:
        if sp.Matrix.vstack(rel, coeffs).rank() != relrank:
            okb = False
    elif any(coeffs):
        okb = False
report("cost of an allowed Gram = budget diagonal (q,p,j,c,e,c+d,b,a)", okb)

# doc's nontrivial relations for the Theorem 1 support and Theorem 2 core
def ent(k, l):
    return gsym[ENT.index((min(k, l), max(k, l)))]
Qi, Pi, Ji, Ci, Ei, Di, Bi, Ai = range(8)
core = [Qi, Pi, Ji, Di, Bi]
zero_rows = [Ci, Ei, Ai]
sub0 = {ent(z, m): 0 for z in zero_rows for m in range(8)}
rel_core = sp.Matrix([[sp.expand((rel[i, :]*sp.Matrix(gsym))[0].subs(sub0)).coeff(gsym[j]) for j in range(36)]
                      for i in range(rel.rows)])
seven = [ent(Qi, Pi) + ent(Qi, Bi) + ent(Pi, Bi) + ent(Di, Di), ent(Qi, Ji) + ent(Ji, Di), ent(Qi, Di),
         ent(Pi, Ji) + ent(Ji, Bi) + ent(Di, Di), ent(Pi, Di), ent(Ji, Di) + ent(Ji, Bi), ent(Di, Bi)]
S7 = sp.Matrix([[sp.expand(e).coeff(gsym[j]) for j in range(36)] for e in seven])
r1, r2, r12 = rel_core.rank(), S7.rank(), sp.Matrix.vstack(rel_core, S7).rank()
report("Thm 2: restricted relations = the seven displayed equations (ranks)", r1 == r2 == r12 == 7,
       info=f"ranks {r1},{r2},{r12}")

# ---------------------------------------------------------------------------------
# 5. Theorem 1 Gram (built from the PDF) satisfies every equation exactly
q, b, j, be, t = sp.symbols('q b j beta t', real=True)
s_ = q + b + 2*be
h = 1 + be
ell = q*b - be**2
u_ = (q + b - 2*be)/ell
Dl = s_*j - (s_ + 4)*t**2
metric = sp.Matrix([[q, be, -t, 0], [be, b, -t, 0], [-t, -t, j - t**2, 0], [0, 0, 0, 1]])
lam = (-h*(j - t**2) + 2*t*(t - 1))/Dl
mu = (-2*t*h + s_*(t - 1))/Dl
# vectors in basis (Q,B,J0,D): Q, P=lam(Q+B)+mu J0, J=J0+tD, 0,0, D, B, 0
VV = sp.Matrix([[1, 0, 0, 0], [lam, lam, mu, 0], [0, 0, 1, t], [0]*4, [0]*4, [0, 0, 0, 1], [0, 1, 0, 0], [0]*4])
Gn = (VV*metric*VV.T).applyfunc(sp.cancel)
p0 = ((j - t**2)*h**2 - 4*t*h*(t - 1) + s_*(t - 1)**2)/Dl
report("Thm1: P.(Q+B) = -h", sp.cancel((sp.Matrix([[lam, lam, mu, 0]])*metric*sp.Matrix([1, 1, 0, 0]))[0] + h) == 0)
report("Thm1: P.J0 = t-1", sp.cancel((sp.Matrix([[lam, lam, mu, 0]])*metric*sp.Matrix([0, 0, 1, 0]))[0] - (t - 1)) == 0)
report("Thm1: diag = (q,p0,j,0,0,1,b,0)",
       all(sp.cancel(Gn[i, i] - v) == 0 for i, v in enumerate([q, p0, j, 0, 0, 1, b, 0])))
report("Thm1: Gram of (Q+B,J0) = [[s,-2t],[-2t,j-t^2]], det = Delta",
       sp.cancel(s_*(j - t**2) - 4*t**2 - Dl) == 0)
# PD of metric under hypotheses: leading minors q, ell, ell*(j-(1+u)t^2), same*1
m3 = sp.cancel(metric[:3, :3].det() - ell*(j - (1 + u_)*t**2))
report("Thm1: det of 3x3 block = ell*(j-(1+u)t^2) (Schur complement)", m3 == 0)
report("Thm1: u - 4/s = (q-b)^2/(s ell)  (so Delta>0 follows)", sp.cancel(u_ - 4/s_ - (q - b)**2/(s_*ell)) == 0)
vals = {gsym[e]: sp.cancel(Dl**2*Gn[k, l]) for e, (k, l) in enumerate(ENT)}
costv = [Dl**2*x for x in [q, p0, j, 0, 0, 1, b, 0]]
okT1 = True
for i in range(rk):
    expr = sum(Rm[i, e]*vals[gsym[e]] for e in range(36)) + sum(Rm[i, 36+m]*costv[m] for m in range(8))
    if sp.cancel(expr) != 0:
        okT1 = False
report("Thm1 CERTIFICATE: |p|^2 * sum G psi psi == 2 Delta^2 Phi_c mod (|p|^2-|q|^2), "
       "c = p0 f1 + j f59 + q f27 + b f6 + f28, all 5 parameters symbolic", okT1)

# ---------------------------------------------------------------------------------
# 6. elementary algebra of pages 2-5
report("t=0 threshold", sp.cancel(p0.subs(t, 0) - (1/j + h**2/s_)) == 0)
report("d/dt at 0", sp.cancel(sp.diff(p0, t).subs(t, 0) - 2*(2 - q - b)/(s_*j)) == 0)
wit = {q: sp.Rational(1, 8), b: 4, j: 4, be: sp.Rational(-1, 2), t: sp.Rational(1, 4)}
report("witness hypotheses and 399/1543 < 13/50",
       ell.subs(wit) > 0 and (j - (1 + u_)*t**2).subs(wit) > 0 and sp.cancel(p0.subs(wit)) == sp.Rational(399, 1543)
       and sp.Rational(399, 1543) < sp.Rational(13, 50))
old = sp.Rational(1, 4) + (1 - 1/sp.sqrt(2))**2/(sp.Rational(33, 8) - sp.sqrt(2))
xr = sp.expand(sp.radsimp(old - sp.Rational(13, 50)))
a0, cs = xr.as_independent(sp.sqrt(2))
a1 = sp.simplify(cs/sp.sqrt(2))
sgn_ok = a0.is_Rational and a1.is_Rational and ((a0 > 0 and a1 >= 0) or (a0 > 0 and a1 < 0 and a0**2 > 2*a1**2)
                                                or (a0 <= 0 and a1 > 0 and 2*a1**2 > a0**2))
report("old parallel threshold > 13/50 (exact, a + c sqrt2 sign test)", bool(sgn_ok), info=f"a={a0}, c={a1}")
kk = 2 - q - b
AA, BB, CC, DD, EE = s_ - h**2 - 4*h, 4*h - 2*s_, j*h**2 + s_, s_*j, s_ + 4
report("threshold = (At^2+Bt+C)/(D-Et^2)", sp.cancel(p0 - (AA*t**2 + BB*t + CC)/(DD - EE*t**2)) == 0)
report("stationary quadratic", sp.cancel(sp.diff(p0, t)*(DD - EE*t**2)**2 - (BB*EE*t**2 + 2*(AA*DD + EE*CC)*t + BB*DD)) == 0)
report("(2) decomposition", sp.cancel(p0 - h**2/s_ - (s_ + kk*t)**2/(s_*Dl)) == 0)
report("(2) derivative", sp.cancel(sp.diff(p0, t) - 2*(s_ + kk*t)*(kk*j + (s_ + 4)*t)/Dl**2) == 0)
report("product of the two critical |t| = sj/(s+4)", sp.cancel((s_/kk)*(kk*j/(s_ + 4)) - s_*j/(s_ + 4)) == 0)
tst = -kk*j/(s_ + 4)
report("threshold at t*", sp.cancel(p0.subs(t, tst) - (1/j + (s_ + 2*kk)/4 + kk**2/(4*(s_ + 4)))) == 0)
sb = {be: 1 - q - b}
report("beta=1-q-b: s=h=k; t* threshold = 1/j + k(k+3)/(k+4)",
       sp.cancel((s_ - kk).subs(sb)) == 0 and sp.cancel((h - kk).subs(sb)) == 0 and
       sp.cancel(p0.subs(sb).subs(t, -kk*j/(kk + 4)) - (1/j + kk*(kk + 3)/(kk + 4))) == 0)
report("t* allowed <=> j < (k+4)^2/(k^2(1+u)) (at beta=1-q-b)",
       sp.cancel(((kk*j/(kk + 4))**2 - j/(1 + u_)).subs(sb) - (j/(1 + u_.subs(sb)))*((kk**2*(1 + u_.subs(sb))*j)/(kk + 4)**2 - 1)) == 0)
report("Cor3: t=0, beta=1-S threshold = 1/j + 2 - S", sp.cancel(p0.subs({t: 0, be: 1 - q - b}) - (1/j + 2 - q - b)) == 0)
report("Cor3: d/dbeta (1+beta)^2/(S+2beta)", sp.cancel(sp.diff(h**2/s_, be) - 2*(1 + be)*(q + b + be - 1)/s_**2) == 0)
report("Cor3: discriminant b(4-3b)", sp.expand(sp.discriminant((q + b - 1)**2 - q*b, q) - b*(4 - 3*b)) == 0)
eps = sp.symbols('eps', positive=True)
gq = (q*b - (q + b - 1)**2).subs({q: 1 - eps, b: 1 - eps})
report("Cor3: q=b=1-eps strict for small eps (= 2eps-3eps^2)", sp.expand(gq - (2*eps - 3*eps**2)) == 0)
tmpl = sp.Matrix([[1, -sp.Rational(1, 2), 0, 0, 0], [-sp.Rational(1, 2), 2, -1, 0, -sp.Rational(1, 2)],
                  [0, -1, 1, 0, 0], [0, 0, 0, 1, 0], [0, -sp.Rational(1, 2), 0, 0, 1]])   # order Q,P,J,D,B
tpos = all(tmpl[:k, :k].det() > 0 for k in range(1, 6))
tv = {ent(core[a_], core[b_]): tmpl[a_, b_] for a_ in range(5) for b_ in range(5)}
teq = all(sp.expand(e.subs(tv)) == 0 for e in seven)
report("Thm2 template positive definite and satisfies the seven equations", tpos and teq,
       info=str([tmpl[:k, :k].det() for k in range(1, 6)]))

# ---------------------------------------------------------------------------------
# 7. role transports through the symmetry group G (lib convention)
FV = {tuple(v): jx for jx, v in enumerate(L.F_VALUES)}


def img(g, r):
    return FV.get(tuple(int(x) for x in L.g_act_values(g, L.F_VALUES[r])))


def find(src, dst):
    return [g for g in L.G_ELEMENTS if [img(g, r) for r in src] == dst]


g1 = find([59, 53, 55, 67, 54], [1, 59, 6, 27, 28])
g2 = find([1, 59, 6, 27, 28], [45, 5, 7, 0, 35])
report("T24 roles [59,53,55,67,54] -> [1,59,6,27,28] by some g in G", len(g1) > 0, info=f"{len(g1)} elements")
report("T46 roles [1,59,6,27,28] -> [45,5,7,0,35] by some g in G", len(g2) > 0, info=f"{len(g2)} elements")
lab = lambda r1, r2, e: (L.perm_from_oneline(r1), L.perm_from_oneline(r2), e)
print("[INFO] group elements (rho, rho', eps) in lib convention, 1-based one-line:",
      [(L.perm_to_oneline(a), L.perm_to_oneline(c), e) for a, c, e in g1], [(L.perm_to_oneline(a), L.perm_to_oneline(c), e) for a, c, e in g2])
report("T47 roles lie in the generators of paper T_65",
       set([45, 5, 7, 0, 35]) <= set(L.generators(L.CONTACT_SETS[65])), info=str(L.generators(L.CONTACT_SETS[65])))
report("T46 roles lie in the generators of paper T_47",
       set([1, 59, 6, 27, 28]) <= set(L.generators(L.CONTACT_SETS[47])))

# ---------------------------------------------------------------------------------
# 8. NUMERIC falsification
import numpy as np  # noqa: E402
rng = np.random.default_rng(20261125)
worst = np.inf
nsamp = 0
f = sp.lambdify((q, b, j, be, t), p0)
uf = sp.lambdify((q, b, be), u_)
while nsamp < 40:
    qq, bb = np.exp(rng.uniform(-2.5, 2.5, 2))
    rr = np.sqrt(qq*bb)
    bet = rng.uniform(-0.98, 0.98)*rr
    jj = np.exp(rng.uniform(-2.5, 2.5))
    T_ = np.sqrt(jj/(1 + uf(qq, bb, bet)))
    tt = rng.uniform(-0.98, 0.98)*T_
    pp = f(qq, bb, jj, bet, tt)
    c = [pp*float(x) for x in L.F_VECS[1]]
    for coef, r in [(jj, 59), (qq, 27), (bb, 6), (1.0, 28)]:
        c = [ci + coef*float(x) for ci, x in zip(c, L.F_VECS[r])]
    val, _ = L.min_phi_over_SO4(c, restarts=12, seed=nsamp)
    worst = min(worst, val/(1 + pp + jj + qq + bb))
    nsamp += 1
report("Thm1 at p = p0: min over SO(4) of Phi_c / (sum of coefficients)", worst > -1e-9, level="NUMERIC",
       info=f"{nsamp} samples, 12 BFGS restarts each, seed 20261125, worst {worst:.3e}")
# pointwise identity check of the certificate at random unit (p,q) and random parameters
def evalp(poly, x):
    return float(sum(float(cf)*np.prod([x[i]**int(e) for i, e in enumerate(m)]) for m, cf in poly.to_dict().items()))
mx = 0.0
for trial in range(5):
    pv = rng.normal(size=4); pv /= np.linalg.norm(pv); qv = rng.normal(size=4); qv /= np.linalg.norm(qv)
    x = np.concatenate([pv, qv])
    qq, bb = 0.7, 2.3; bet = -0.4; jj = 1.9; tt = 0.3
    Gnum = np.array(Gn.subs({q: qq, b: bb, be: bet, j: jj, t: tt}).evalf(), dtype=float)
    S = sum(float(Gnum[k, l]*(1 if k == l else 2))*evalp(sum([ps[k]*ps[l] for ps in u_all], ZERO), x)
            for (k, l) in ENT)
    Om = np.array([[evalp(O[i][jx], x) for jx in range(4)] for i in range(4)])
    c = [float(f(qq, bb, jj, bet, tt))*float(v) for v in L.F_VECS[1]]
    for coef, r in [(jj, 59), (qq, 27), (bb, 6), (1.0, 28)]:
        c = [ci + coef*float(v) for ci, v in zip(c, L.F_VECS[r])]
    mx = max(mx, abs(S - 2*L.Phi_np(c, Om)))
report("pointwise S_G(p,q) = 2 Phi_c(O(p,q)) at 5 random unit points", mx < 1e-9, level="NUMERIC", info=f"max err {mx:.2e}")

# Corollary 3 on its boundary (S-1)^2 = qb (d = 1), including q = b = 1 and q = 0, b = 1
worst3 = np.inf
pts3 = [(1.0, 1.0), (0.0, 1.0), (1.0, 0.0)]
for _ in range(17):
    bb = rng.uniform(0.0, 4/3)
    # roots in q of q^2 + (b-2)q + (b-1)^2 = 0
    disc = bb*(4 - 3*bb)
    for sgn in (1, -1):
        qq = ((2 - bb) + sgn*np.sqrt(disc))/2
        if qq >= 0:
            pts3.append((qq, bb))
for n3, (qq, bb) in enumerate(pts3):
    jj = np.exp(rng.uniform(-2, 2))
    pp = 1/jj + 2 - qq - bb
    c = [pp*float(x) for x in L.F_VECS[1]]
    for coef, r in [(jj, 59), (qq, 27), (bb, 6), (1.0, 28)]:
        c = [ci + coef*float(x) for ci, x in zip(c, L.F_VECS[r])]
    val, _ = L.min_phi_over_SO4(c, restarts=12, seed=100 + n3)
    worst3 = min(worst3, val/(1 + abs(pp) + jj + qq + bb))
report("Cor 3 on the boundary (S-1)^2 = qb: min Phi_c / (sum of |coefficients|)", worst3 > -1e-9, level="NUMERIC",
       info=f"{len(pts3)} points, 12 restarts each, worst {worst3:.3e}")

print("SUMMARY:", "ALL EXACT CHECKS PASS" if not FAIL else f"FAILED: {FAIL}")
sys.exit(1 if FAIL else 0)
