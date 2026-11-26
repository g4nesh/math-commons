"""Build an exact rational certificate that  beta - Phi(O) >= 0  for every real orthogonal 5x5 matrix O  (default beta = 13/10).

Shape of the certificate:   beta - Phi  =  sum over the 225 character blocks of  v_b^T G_b v_b   modulo the orthogonality
relations, where v_b is the vector of cubic monomials of block b and every G_b is a positive definite rational matrix.
The Gram entries are constant on orbits of the symmetry group (890 unknowns theta).

Steps.  (1) Exact linear algebra in the space of GAMMA-invariant polynomials of degree <= 6 with trivial sign character:
the subspace J spanned by the averages of  (monomial of degree <= 4) x (orthogonality relation)  and the quotient by J.
(2) A semidefinite programme in floating point for theta, maximising the smallest eigenvalue.  (3) Rounding to rationals and an
exact correction so that the identity holds exactly modulo J.  (4) Multipliers: an exact expression of the remainder as a
combination of the generators of J.  verify_certificate.py re-checks the result without this script.
usage: build_certificate.py [beta numerator] [beta denominator] [rounding denominator]"""
import itertools, json, sys, time
from fractions import Fraction as Fr
import numpy as np, flint
from common import *
t0 = time.time(); BETA = Fr(int(sys.argv[1]) if len(sys.argv) > 1 else 13, int(sys.argv[2]) if len(sys.argv) > 2 else 10); DEN = int(sys.argv[3]) if len(sys.argv) > 3 else 10**7
blocks = cubic_blocks(); key = pair_orbits(blocks); orbit_keys = sorted(set(key.values())); OI = {k: i for i, k in enumerate(orbit_keys)}; NO = len(orbit_keys)
print("cubic monomials %d, blocks %d, pair orbits %d  [%ds]" % (len(CUBICS), len(blocks), NO, time.time() - t0), flush=True)
# (1) orbit-feature polynomials F_o = sum over pairs in the orbit of w m m'  (w = 1 on the diagonal, 2 off it), in orbit coordinates
F = [dict() for _ in range(NO)]
for ms in blocks.values():
    for a in range(len(ms)):
        for b in range(a, len(ms)):
            o = OI[key[(ms[a], ms[b])]]; m = tuple(sorted(ms[a] + ms[b])); rep, size = canon(m); F[o][rep] = F[o].get(rep, 0) + Fr(1 if a == b else 2, size)
phi_c = average(phi_poly()); one_c = {(): Fr(1)}
# generators of J
rels = relations(); mons_by_char = {}
for d in range(0, 5):
    for m in itertools.combinations_with_replacement(range(25), d): mons_by_char.setdefault(char(m), []).append(m)
gens = []; gen_id = []
for ri, r in enumerate(rels):
    ch = char(next(m for m in r if m != ()))
    for m in mons_by_char.get(ch, []):
        g = average(pmul({m: Fr(1)}, r))
        if g: gens.append(g); gen_id.append((ri, m))
print("generators of J: %d  [%ds]" % (len(gens), time.time() - t0), flush=True)
coords = sorted(set().union(*[set(g) for g in gens], *[set(f) for f in F], set(phi_c), {()})); CI = {m: i for i, m in enumerate(coords)}; D = len(coords)
print("invariant trivial-character monomial orbits of degree <= 6 in play: %d" % D, flush=True)
# exact reduced row echelon form of J (rows = generators), keeping only independent generators
def vec(p):
    v = [Fr(0)] * D
    for m, c in p.items(): v[CI[m]] = c
    return v
# choose independent generators greedily by elimination (exact, fraction arithmetic through flint)
A = flint.fmpq_mat(len(gens), D, [flint.fmpq(c.numerator, c.denominator) for g in gens for c in vec(g)])
R, rank = A.rref(); print("rank of J: %d ; quotient dimension %d  [%ds]" % (rank, D - rank, time.time() - t0), flush=True)
piv = []
for i in range(rank):
    j = 0
    while R[i, j] == 0: j += 1
    piv.append(j)
free = [j for j in range(D) if j not in set(piv)]
def reduce(v):
    """coordinates of v modulo J on the free columns"""
    v = v[:]
    for i, j in enumerate(piv):
        if v[j]:
            c = v[j]
            for jj in range(j, D):
                rij = R[i, jj]
                if rij != 0: v[jj] -= c * Fr(int(rij.p), int(rij.q))
    return [v[j] for j in free]
QF = [reduce(vec(f)) for f in F]; q1 = reduce(vec(one_c)); qphi = reduce(vec(phi_c)); Q = len(free)
print("reduced features  [%ds]" % (time.time() - t0), flush=True)
# (2) floating-point SDP: sum_o theta_o QF[o] = beta q1 - qphi ; G_b >= t I ; maximise t
import cvxpy as cp
Aeq = np.array([[float(QF[o][k]) for o in range(NO)] for k in range(Q)]); beq = np.array([float(BETA * q1[k] - qphi[k]) for k in range(Q)])
Ua, sv, Vt = np.linalg.svd(Aeq, full_matrices=False); rk = int((sv > 1e-9 * sv[0]).sum()); Ar = Vt[:rk]; br = (Ua[:, :rk].T @ beq) / sv[:rk]
print("identity conditions: %d independent of %d ; rhs consistent to %.1e" % (rk, Q, np.abs(Aeq @ (Ar.T @ br) - beq).max()), flush=True)
reps = {}
for ch, ms in blocks.items():
    rep_ch = min(char(act(g, ms[0])) for g in GAMMA)
    if rep_ch == ch: reps[ch] = ms
theta = cp.Variable(NO); tt = cp.Variable(); cons = [Ar @ theta == br]
for ch, ms in reps.items():
    s = len(ms); M = np.zeros((s * s, NO))
    for a in range(s):
        for b in range(s): M[a * s + b, OI[key[(ms[min(a, b)], ms[max(a, b)])]]] = 1.0
    cons.append(cp.reshape(M @ theta, (s, s), order='C') >> tt * np.eye(s))
# the smallest beta with a positive semidefinite certificate of this shape (floating point, for the record)
bvar = cp.Variable(); b1 = (Ua[:, :rk].T @ np.array([float(q1[k]) for k in range(Q)])) / sv[:rk]; bp = (Ua[:, :rk].T @ np.array([float(qphi[k]) for k in range(Q)])) / sv[:rk]
pb = cp.Problem(cp.Minimize(bvar), [Ar @ theta == bvar * b1 - bp] + [c_ for c_ in cons[1:]] + [tt == 0]); pb.solve(solver=cp.CLARABEL)
print("smallest beta with a degree-6 certificate of this shape (floating point): %.6f (%s)" % (bvar.value, pb.status), flush=True)
prob = cp.Problem(cp.Maximize(tt), cons + [tt <= 10]); prob.solve(solver=cp.CLARABEL)
print("SDP status %s ; smallest eigenvalue margin t = %.6f  [%ds]" % (prob.status, tt.value, time.time() - t0), flush=True)
if tt.value is None or tt.value <= 1e-6: raise SystemExit("no positive definite certificate at this beta")
# (3) rationalise and correct exactly: theta = round(theta) + delta with A delta = residual, delta of least norm in exact arithmetic
th = [Fr(int(round(x * DEN)), DEN) for x in theta.value]
Aex = flint.fmpq_mat(Q, NO, [flint.fmpq(QF[o][k].numerator, QF[o][k].denominator) for k in range(Q) for o in range(NO)])
res = [BETA * q1[k] - qphi[k] - sum(QF[o][k] * th[o] for o in range(NO) if QF[o][k]) for k in range(Q)]
# solve A delta = res exactly with delta supported on pivot columns of rref([A | res])
Aug = flint.fmpq_mat(Q, NO + 1, [x for k in range(Q) for x in ([Aex[k, o] for o in range(NO)] + [flint.fmpq(res[k].numerator, res[k].denominator)])])
Rr, rr = Aug.rref(); delta = [Fr(0)] * NO
for i in range(rr):
    j = 0
    while Rr[i, j] == 0: j += 1
    if j == NO: raise SystemExit("the identity is inconsistent in exact arithmetic")
    delta[j] = Fr(int(Rr[i, NO].p), int(Rr[i, NO].q))
th = [a + b for a, b in zip(th, delta)]; print("largest exact correction: %.2e  [%ds]" % (max(abs(float(x)) for x in delta), time.time() - t0), flush=True)
assert all(BETA * q1[k] - qphi[k] == sum(QF[o][k] * th[o] for o in range(NO) if QF[o][k]) for k in range(Q)), "identity not exact"
ok = all(ldl_positive([[th[OI[key[(ms[min(a, b)], ms[max(a, b)])]]] for b in range(len(ms))] for a in range(len(ms))]) for ms in blocks.values())
print("all 225 Gram blocks positive definite (exact):", ok, " [%ds]" % (time.time() - t0), flush=True)
if not ok: raise SystemExit("rounded certificate is not positive definite; increase the rounding denominator")
# (4) multipliers: h = beta - Phi - sum theta_o F_o lies in J; express avg(h) in the independent generators
h = dict(one_c); h[()] = BETA; padd(h, phi_c, -1)
for o in range(NO):
    if th[o]: padd(h, F[o], -th[o])
hv = vec(h)
# pick `rank` independent generators and solve  sum mu_k gens_k = h  exactly
At = flint.fmpq_mat(D, len(gens), [A[i, j] for j in range(D) for i in range(len(gens))])        # transpose: columns = generators
AugT = flint.fmpq_mat(D, len(gens) + 1, [x for j in range(D) for x in ([At[j, i] for i in range(len(gens))] + [flint.fmpq(hv[j].numerator, hv[j].denominator)])])
Rt, rt = AugT.rref(); mult = {}
for i in range(rt):
    j = 0
    while Rt[i, j] == 0: j += 1
    if j == len(gens): raise SystemExit("remainder is not in J")
    c = Fr(int(Rt[i, len(gens)].p), int(Rt[i, len(gens)].q))
    if c: mult[j] = c
print("multipliers: %d generators used  [%ds]" % (len(mult), time.time() - t0), flush=True)
cert = {"beta": [BETA.numerator, BETA.denominator],
        "orbit_keys": [[list(k[0]), list(k[1])] for k in orbit_keys],
        "theta": [[x.numerator, x.denominator] for x in th],
        "multipliers": [[gen_id[j][0], list(gen_id[j][1]), c.numerator, c.denominator] for j, c in sorted(mult.items())]}
json.dump(cert, open("certificate.json", "w"), separators=(",", ":")); print("wrote certificate.json  [%ds]" % (time.time() - t0))
