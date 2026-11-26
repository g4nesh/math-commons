"""Gap check for blind_moc4_T44_conic_branch_exclusion_20261022: identity (2) (discriminant
and the two square-class identities mod G) and the X = 0 octic identity of Lemma 2.

Own code. From the paper's definitions (moc4lib: f_j vectors, cubic C of (A.2), kappa):
  g = q f2 + p f1 + j f53 + c f12 + e f8 + d f18 + b f7 + a f15,  w = g + h kappa,
  C(w) = C + hQ,  H(w) = Hessian of the cubic at w (14x14),
  det H(w)|_{h=-C/Q} = -4 G^2 / Q^2  (eq. (1)),  det H[I,I] = (q+d)^2 P,  I = {2..13}.
Phi is computed from the document's displayed formula; Psi := (B^2 - 4 A C0)/(e^2 Phi).
The existential polynomials S, T, Omega (not printed in the document) are taken from the
supplement's JSON data file blind_moc4_T44_principal_squareclass_20261022.json, loaded
with json.load as a candidate certificate only; U := A*P. The radial polynomial G is
recomputed here; the JSON copies of G, P, Phi, Psi are only compared against.
Nothing of the supplement is executed or imported.

Usage: MOC4_SUPP=<supplement root> python -I -B this.py   (or pass the root as argv[1])
"""
import os, sys, json, time
import flint

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lib'))
import moc4lib as L

SUPP = sys.argv[1] if len(sys.argv) > 1 else os.environ['MOC4_SUPP']
DD = os.path.join(SUPP, 'proofs', 'manin_transfer_20261021')
FAIL = 0
T0 = time.time()
def check(name, ok):
    global FAIL
    print(('PASS ' if ok else 'FAIL ') + name + '   [%.1fs]' % (time.time() - T0), flush=True)
    if not ok: FAIL += 1

ctx = flint.fmpq_mpoly_ctx.get(('q', 'p', 'j', 'c', 'e', 'd', 'b', 'a', 'h'), 'lex')
q, p, j, c, e, d, b, a, h = ctx.gens(); Z = 0 * q; ONE = Z + 1
GENS = ctx.gens()

def fq(x): return flint.fmpq(x.numerator, x.denominator)
def rd(ts):
    out = Z
    for ex, co in ts:
        num, den = (co.split('/') + ['1'])[:2]
        t = flint.fmpq(int(num), int(den))
        for i, ee in enumerate(ex):
            if ee: t = t * GENS[i] ** ee
        out = out + t
    return out
def exact_div(x, y):
    qq, rr = divmod(x, y); assert rr.is_zero(); return qq
def coeff_in(f, idx, k):
    return ctx.from_dict({tuple(0 if z == idx else x for z, x in enumerate(ex)): co
                          for ex, co in f.to_dict().items() if ex[idx] == k})
def deg(f, idx): return int(f.degrees()[idx])
def subst_rat(f, idx, num, den):
    """den^deg_idx(f) * f(var_idx = num/den)."""
    n = deg(f, idx)
    return sum((coeff_in(f, idx, k) * num ** k * den ** (n - k) for k in range(n + 1)), Z)

# ---- cost vector, cubic restriction, Hessian ----
SRC = {2: q, 1: p, 53: j, 12: c, 8: e, 18: d, 7: b, 15: a}
g = [sum((fq(L.F_VECS[f][i]) * x for f, x in SRC.items()), Z) for i in range(14)]
check('kappa = e_13', list(L.KAPPA) == [0] * 13 + [1])
w = [g[i] + (h if i == 13 else 0) for i in range(14)]
Cw = L.cubic_C(w); Cg = L.cubic_C(g)
Qg = exact_div(Cw - Cg, h)
check('C(g + h kappa) = C + h Q with C, Q independent of h', deg(Qg, 8) == 0 and deg(Cg, 8) == 0)
rad = json.load(open(os.path.join(DD, 'blind_moc4_T44_full_radial_20261022.json')))
check('JSON variable order q,p,j,c,e,d,b,a,(z)', rad['variables'][:8] == ['q', 'p', 'j', 'c', 'e', 'd', 'b', 'a'])
check('own C, Q equal the JSON C_terms, Q_terms (comparison only)', Cg == rd(rad['C_terms']) and Qg == rd(rad['Q_terms']))

mc = flint.fmpq_mpoly_ctx.get(tuple('m%d' % i for i in range(14)), 'lex'); ms = mc.gens()
Cm = L.cubic_C(list(ms))
H = [[None] * 14 for _ in range(14)]
for i in range(14):
    for k in range(14):
        lin = Cm.derivative(i).derivative(k)
        val = Z
        for ex, co in lin.to_dict().items():
            t = flint.fmpq(co)
            for idx, ee in enumerate(ex):
                if ee: t = t * w[idx] ** ee
            val = val + t
        H[i][k] = val

def bareiss(M):
    M = [row[:] for row in M]; n = len(M); prev = ONE; sign = 1
    for kk in range(n - 1):
        if M[kk][kk].is_zero():
            sw = next((r for r in range(kk + 1, n) if not M[r][kk].is_zero()), None)
            if sw is None: return Z
            M[kk], M[sw] = M[sw], M[kk]; sign = -sign
        for i in range(kk + 1, n):
            for jj in range(kk + 1, n):
                M[i][jj] = exact_div(M[i][jj] * M[kk][kk] - M[i][kk] * M[kk][jj], prev)
        prev = M[kk][kk]
    return sign * M[n - 1][n - 1]

I = list(range(2, 14))
D12 = bareiss([[H[i][k] for k in I] for i in I])
check('principal minor on I={2..13} is independent of h', deg(D12, 8) == 0)
P = exact_div(D12, (q + d) ** 2)
hp = json.load(open(os.path.join(DD, 'blind_moc4_T44_height_pencil_20261022.json')))
check('P := det H[I,I]/(q+d)^2 is a polynomial; degree %d, %d terms; equals JSON P10' %
      (P.total_degree(), len(P.to_dict())), P == rd(hp['polynomials']['P10']['terms']))

D14 = bareiss(H)
check('full det H(g + h kappa) has h-degree <= 4 (got %d)' % deg(D14, 8), deg(D14, 8) <= 4)
# Q^4 * det at h = -C/Q
forced = subst_rat(D14, 8, -Cg, Qg) * Qg ** (4 - deg(D14, 8))      # = Q^4 det(-C/Q)
G2 = exact_div(forced, -4 * Qg ** 2)                                 # should be G^2
Gm = G2.sqrt()
G_json = -rd(next(f['terms'] for f in rad['radial_factors'] if f['degree'] == 9))
check('eq.(1): Q^4 det H(w(g,-C/Q)) = -4 Q^2 G^2 with G a polynomial (%d terms, degree %d); '
      'G = JSON radial factor up to sign' % (len(Gm.to_dict()), Gm.total_degree()),
      Gm * Gm == G2 and (Gm == G_json or Gm == -G_json))
G = G_json            # the supplement's sign convention; every identity below is sign-symmetric
check('G has q-degree 2', deg(G, 0) == 2)
A = coeff_in(G, 0, 2); B = coeff_in(G, 0, 1); C0 = coeff_in(G, 0, 0)
check('G = A q^2 + B q + C0', A * q ** 2 + B * q + C0 == G)
check('gcd(A, B, C0) = 1 (q-content of G is one)', A.gcd(B).gcd(C0) == 1)

# ---- (2), first line ----
W = c + d + b
R = (d + a) * (j - c) * p + j * c * d + j * (d + b) * (d - a) + c * b * (d + a)
Phi = R ** 2 - 4 * d * j * (b - a) * W * (c * (d + a) - j * a)
sq = json.load(open(os.path.join(DD, 'blind_moc4_T44_principal_squareclass_20261022.json')))
PP = sq['polynomials']
check('Phi (document formula) equals JSON F6 (comparison only)', Phi == rd(PP['F6']['terms']))
disc = B ** 2 - 4 * A * C0
Psi = exact_div(disc, e ** 2 * Phi)
check('(2a) B^2 - 4 A C0 = e^2 Phi Psi with Psi a polynomial: degree %d, %d terms, q-free' %
      (Psi.total_degree(), len(Psi.to_dict())), Psi * e ** 2 * Phi == disc and deg(Psi, 0) == 0)
check('Psi equals JSON F8 (comparison only)', Psi == rd(PP['F8']['terms']))

# ---- (2), square classes; S, T, Omega from the JSON certificate; U := A P ----
S = rd(PP['H8']['terms']); T = rd(PP['K8']['terms']); Om = rd(PP['W17']['terms'])
check('S, T octics, Omega degree 17, all q-free (and h-free)',
      S.total_degree() == 8 and T.total_degree() == 8 and Om.total_degree() == 17
      and deg(S, 0) == deg(T, 0) == deg(Om, 0) == 0 and deg(S, 8) == deg(T, 8) == deg(Om, 8) == 0)
check('S, T, Omega nonzero', not S.is_zero() and not T.is_zero() and not Om.is_zero())
U = A * P
lhs1 = P * (e ** 2 * Phi * S ** 2) - (U + Om) ** 2
q1, r1 = divmod(lhs1, G)
check('(2b) P (e^2 Phi S^2) - (U + Omega)^2 is divisible by G', r1.is_zero() and q1 * G == lhs1)
lhs2 = P * (Psi * T ** 2) - (U - Om) ** 2
q2, r2 = divmod(lhs2, G)
check('(2c) P (Psi T^2) - (U - Omega)^2 is divisible by G', r2.is_zero() and q2 * G == lhs2)

# ---- X = 0 octic identity of Lemma 2 ----
Zs = p * j - p * c - d * W
E1 = (p - b) * Zs - j * d * (p + W); E0 = p * (d + b) * Zs
check('Psi has a-degree 2 (so W^2 Psi(a=pc/W) is the full clearing)', deg(Psi, 7) == 2)
check('W^2 Psi|_{a=pc/W} = c^2 (E1 e + E0)^2', subst_rat(Psi, 7, p * c, W) == c ** 2 * (E1 * e + E0) ** 2)

# ---- the other displayed identities of Lemma 2 (not in dispute; cheap) ----
Dd = c * b * (d + a) - a * j * (d + b); X = p * c - a * W; Y = -p * Dd + a * j * d * W
Mterm = d * j * (b - a) * W * (c * (d + a) - j * a)
check('Lemma 2: D R = D^2 + dj(b-a)W{c(d+a)-ja} - (d+a)(j-c) Y', Dd * R == Dd ** 2 + Mterm - (d + a) * (j - c) * Y)
check('Lemma 2: D^2 Phi|_{p=ajdW/D} = (D^2 - dj(b-a)W{...})^2',
      subst_rat(Phi, 1, a * j * d * W, Dd) == (Dd ** 2 - Mterm) ** 2)
M0 = p * c + d * W; L1 = M0 ** 2 + d * j * W ** 2 + (d + b) * W * Zs
L0 = d * L1 + p * c * (d + b) * Zs; T3 = p * Zs - d * (p * c + (j + d) * W)
GX = subst_rat(G, 7, p * c, W)
check('Lemma 2: G has a-degree %d; W^2 E1^2 G|_{a=pc/W, e=-E0/E1} = d(d+b) Z T3 (q L1 + L0)^2' % deg(G, 7),
      deg(G, 7) == 2 and deg(GX, 4) == 2 and subst_rat(GX, 4, -E0, E1) == d * (d + b) * Zs * T3 * (q * L1 + L0) ** 2)
check('Lemma 2: L1 + W E1 = p{pjW - (d+b)(pc+dW)}', L1 + W * E1 == p * (p * j * W - (d + b) * M0))
# conic (3)
A2 = coeff_in(G, 4, 2); B1 = exact_div(coeff_in(G, 4, 1), q + d); K = exact_div(coeff_in(G, 4, 0), (q + d) ** 2)
check('G = e^2 A2 + e (q+d) B1 + (q+d)^2 K, A2 quadratic and B1 linear in q',
      G == e ** 2 * A2 + e * (q + d) * B1 + (q + d) ** 2 * K and deg(A2, 0) == 2 and deg(B1, 0) == 1)
AA = [coeff_in(A2, 0, k) for k in range(3)]; BB = [coeff_in(B1, 0, k) for k in range(2)]
Hc = [[2 * K, BB[0], BB[1]], [BB[0], 2 * AA[0], AA[1]], [BB[1], AA[1], 2 * AA[2]]]
det3 = (Hc[0][0] * Hc[1][1] * Hc[2][2] + 2 * Hc[0][1] * Hc[0][2] * Hc[1][2]
        - Hc[0][0] * Hc[1][2] ** 2 - Hc[1][1] * Hc[0][2] ** 2 - Hc[2][2] * Hc[0][1] ** 2)
check('(3) det H_conic = 2 j X Y Phi^2', det3 == 2 * j * X * Y * Phi ** 2)

print('RESULT', 'ALL EXACT CHECKS PASS' if FAIL == 0 else '%d FAILURES' % FAIL)
sys.exit(1 if FAIL else 0)
