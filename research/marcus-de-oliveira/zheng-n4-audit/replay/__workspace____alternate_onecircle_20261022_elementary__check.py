"""Exact replay and numeric falsification for the note
"An elementary compensation inequality for two rotations" (alternate_onecircle_20261022_elementary).

Own code only; nothing from the supplement is read or executed.  The three copies of this note
in the supplement (__workspace__, principal_minor_result_20261022, moc_principal_minors_20261022)
have identical mathematical content, so the sibling check scripts call run_all() here.

Every printed line is one check.  Exit status is nonzero iff an EXACT check fails.
"""
import sys
import random
import math
import sympy as sp

FAIL = []


def check(name, ok, level="EXACT"):
    print(f"[{level}] {'PASS' if ok else 'FAIL'}  {name}")
    if not ok and level == "EXACT":
        FAIL.append(name)


def zero(e):
    return sp.expand(e) == 0


def exact_checks():
    # ---------------- Lemma 1 --------------------------------------------------------
    # u,x,v,y >= 0 written as squares U^2.. so that all square roots are polynomial.
    U, X, V, Y = sp.symbols("U X V Y", nonnegative=True)
    lhs = (U**2 + X**2) * (V**2 + Y**2) - (U * V + X * Y) ** 2
    check("Lemma 1: (u+x)(v+y) - (sqrt(uv)+sqrt(xy))^2 = (sqrt(uy)-sqrt(vx))^2",
          zero(lhs - (U * Y - V * X) ** 2))

    # ---------------- Theorem 1 ------------------------------------------------------
    # square-root variables: a = ra^2 etc.
    ra, rb, rc, rd, RA, RB, RC, RD = sp.symbols("ra rb rc rd RA RB RC RD", nonnegative=True)
    a, b, c, d, A, B, C, D = [v**2 for v in (ra, rb, rc, rd, RA, RB, RC, RD)]
    Rf = a*b + c*d + A*B + C*D - 2*rb*rc*RA*RD - 2*ra*rd*RB*RC
    s = b + c - B - C
    t = a + b - A - B
    tot = a + b + c + d - A - B - C - D

    def subst_perm(expr, perm):
        lo = [ra, rb, rc, rd]
        up = [RA, RB, RC, RD]
        m = {}
        for i, j in enumerate(perm):
            m[lo[i]] = lo[j]
            m[up[i]] = up[j]
        return expr.xreplace(m)

    # (a,b,c,d) -> (b,a,d,c): variable a is replaced by b etc.
    p1 = (1, 0, 3, 2)
    p2 = (3, 2, 1, 0)
    check("Thm 1 symmetry (b,a,d,c): R invariant", zero(subst_perm(Rf, p1) - Rf))
    # s -> -s modulo the total-sum relation: s' + s = (a+d-A-D)+(b+c-B-C) = tot
    check("Thm 1 symmetry (b,a,d,c): s' = -s + (total difference)",
          zero(subst_perm(s, p1) + s - tot))
    check("Thm 1 symmetry (b,a,d,c): t' = t", zero(subst_perm(t, p1) - t))
    check("Thm 1 symmetry (d,c,b,a): R invariant", zero(subst_perm(Rf, p2) - Rf))
    check("Thm 1 symmetry (d,c,b,a): s' = s", zero(subst_perm(s, p2) - s))
    check("Thm 1 symmetry (d,c,b,a): t' = -t + (total difference)",
          zero(subst_perm(t, p2) + t - tot))
    # the two symmetries generate sign changes of s and t independently: check on a grid
    # of the four combinations that each achieves the claimed sign pattern (done above).

    h1, h2 = rb * RA, ra * RB          # h = (sqrt(bA), sqrt(aB))
    k1, k2 = rc * RD, rd * RC          # k = (sqrt(cD), sqrt(dC))
    m_ = b - B
    n_ = D - d
    F = Rf + s**2 + s*t
    nhk = (h1 - k1)**2 + (h2 - k2)**2
    nhJk = (h1 - k2)**2 + (h2 - k1)**2
    # (3) holds modulo the total-sum relation: F - (|h-k|^2 + 2mn) is a multiple of tot.
    diff3 = sp.expand(F - nhk - 2*m_*n_)
    # eliminate D (i.e. RD^2) via D = a+b+c+d-A-B-C; RD appears only squared in diff3?
    q, r = sp.div(sp.Poly(diff3, RD), sp.Poly(RD**2 - (a + b + c + d - A - B - C), RD))
    check("Thm 1 eq (3): F = |h-k|^2 + 2mn modulo a+b+c+d = A+B+C+D", r.as_expr() == 0)
    check("Thm 1 eq (4): |h-k|^2 + 2mn = |h-Jk|^2 + 2(mn - (h1-h2)(k1-k2))",
          zero(nhk + 2*m_*n_ - (nhJk + 2*(m_*n_ - (h1 - h2)*(k1 - k2)))))

    # case analysis: parametrize each case by free variables and check the claimed relations
    A_, B_, C_, D_, M, N, al, be, s_, t_ = sp.symbols("A B C D M N alpha beta s t")
    # Case b<B, d<D: free A,B,C, M>0, s,t>=0; b = B - M, a = A + (M+t), c = C + (M+s),
    # d determined by total equality.
    bb = B_ - M
    aa = A_ + M + t_
    cc = C_ + M + s_
    DD = sp.Symbol("Dfree")
    dd = A_ + B_ + C_ + DD - aa - bb - cc
    check("Case b<B,d<D: s = b+c-B-C", zero((bb + cc - B_ - C_) - s_))
    check("Case b<B,d<D: t = a+b-A-B", zero((aa + bb - A_ - B_) - t_))
    check("Case b<B,d<D: N := D-d = M+s+t", zero((DD - dd) - (M + s_ + t_)))
    Nc = M + s_ + t_
    check("Case b<B,d<D: alpha*beta = MN + st", zero((M + t_)*(M + s_) - (M*Nc + s_*t_)))
    check("Case b<B,d<D: mn = -MN", zero((bb - B_)*(DD - dd) + M*Nc))
    # Case b>B, d>D: b = B + M, a = A - (M - t), c = C - (M - s)
    bb = B_ + M
    aa = A_ - (M - t_)
    cc = C_ - (M - s_)
    dd = A_ + B_ + C_ + DD - aa - bb - cc
    check("Case b>B,d>D: s = b+c-B-C", zero((bb + cc - B_ - C_) - s_))
    check("Case b>B,d>D: t = a+b-A-B", zero((aa + bb - A_ - B_) - t_))
    Nc = dd - DD
    check("Case b>B,d>D: M = N+s+t", zero(M - (Nc + s_ + t_)))
    check("Case b>B,d>D: alpha = M-t = N+s and beta = M-s = N+t",
          zero((M - t_) - (Nc + s_)) and zero((M - s_) - (Nc + t_)))
    check("Case b>B,d>D: alpha*beta = MN + st", zero((M - t_)*(M - s_) - (M*Nc + s_*t_)))
    check("Case b>B,d>D: mn = -MN", zero((bb - B_)*(DD - dd) + M*Nc))
    # Lemma-1 instances: in case 1, h2 = sqrt(aB) = sqrt((A+alpha)(b+M)), h1 = sqrt(Ab);
    # k1 = sqrt(cD) = sqrt((C+beta)(d+N)), k2 = sqrt(Cd).  These are substitutions (a = A+alpha,
    # B = b+M, c = C+beta, D = d+N), which hold by the definitions of alpha, M, beta, N.
    check("Case analysis is exhaustive given mn<0: (b-B)(D-d)<0 <=> b-B, d-D same strict sign",
          True, level="PROOF-READ")

    # ---------------- Theorem 2 ------------------------------------------------------
    w, x, y, z, W, Xq, Yq, Z = sp.symbols("w x y z W X Y Z", real=True)

    def rotmat(w, x, y, z):
        # rotation v -> q v q^* for q = w + x i + y j + z k (computed, not copied)
        qv = sp.Matrix([w, x, y, z])

        def qmul(p, q):
            p0, p1, p2, p3 = p
            q0, q1, q2, q3 = q
            return [p0*q0 - p1*q1 - p2*q2 - p3*q3,
                    p0*q1 + p1*q0 + p2*q3 - p3*q2,
                    p0*q2 - p1*q3 + p2*q0 + p3*q1,
                    p0*q3 + p1*q2 - p2*q1 + p3*q0]
        cols = []
        for e in ([0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]):
            r = qmul(qmul(list(qv), e), [w, -x, -y, -z])
            cols.append([sp.expand(r[1]), sp.expand(r[2]), sp.expand(r[3])])
        return sp.Matrix(3, 3, lambda i, j: cols[j][i])

    Pm = rotmat(w, x, y, z)
    Qm = rotmat(W, Xq, Yq, Z)
    nrm = w**2 + x**2 + y**2 + z**2
    # on the unit sphere, nrm = 1; the doc's formulas use 1 - 2(...) i.e. nrm^2 = 1 already used.
    # Check doc's convention against q v q^* (homogeneous version: replace 1 by nrm).
    conv = {
        (0, 0): nrm - 2*(y**2 + z**2), (2, 2): nrm - 2*(x**2 + y**2),
        (0, 1): 2*(x*y - w*z), (1, 0): 2*(x*y + w*z),
        (1, 2): 2*(y*z - w*x), (2, 1): 2*(y*z + w*x)}
    check("Thm 2: doc's quaternion convention equals v -> q v q^* (homogeneous form)",
          all(zero(Pm[i, j] - v) for (i, j), v in conv.items()))
    check("Thm 2: q v q^* matrix is orthogonal with det 1 on |q|=1 (homogeneous: R R^t = nrm^2 I, det = nrm^3)",
          all(zero(e) for e in (Pm * Pm.T - nrm**2 * sp.eye(3))) and zero(Pm.det() - nrm**3))
    # Use the doc's (unit-norm) formulas
    def doc_entries(w, x, y, z):
        return {"11": 1 - 2*(y**2 + z**2), "33": 1 - 2*(x**2 + y**2),
                "12": 2*(x*y - w*z), "21": 2*(x*y + w*z),
                "23": 2*(y*z - w*x), "32": 2*(y*z + w*x)}
    Pd = doc_entries(w, x, y, z)
    Qd = doc_entries(W, Xq, Yq, Z)
    d0 = Pd["11"] - Qd["11"]
    e0 = Pd["33"] - Qd["33"]
    G = 2*d0**2 + Pd["12"]**2 + Qd["12"]**2 + Pd["21"]**2 + Qd["21"]**2 + 2*(Pd["32"]*Qd["32"] - Pd["23"]*Qd["23"])
    aq, bq, cq, dq = x**2, y**2, z**2, w**2
    Aq, Bq, Cq, Dq = Xq**2, Yq**2, Z**2, W**2
    sq = bq + cq - Bq - Cq
    tq = aq + bq - Aq - Bq
    check("Thm 2: d0 = -2s (identity, no norm relation needed)", zero(d0 + 2*sq))
    check("Thm 2: e0 = -2t (identity)", zero(e0 + 2*tq))
    check("Thm 2: G/8 = s^2 + ab + cd + AB + CD + 2(yzWX + wxYZ) (identity)",
          zero(G/8 - (sq**2 + aq*bq + cq*dq + Aq*Bq + Cq*Dq + 2*(y*z*W*Xq + w*x*Yq*Z))))
    check("Thm 2: |yzWX|^2 = bcAD and |wxYZ|^2 = adBC",
          zero((y*z*W*Xq)**2 - bq*cq*Aq*Dq) and zero((w*x*Yq*Z)**2 - aq*dq*Bq*Cq))
    check("Thm 2: 2|d0 e0|/8 = |st| (since d0 e0 = 4 st)", zero(d0*e0 - 4*sq*tq))


def numeric_checks(seed=20261125):
    import numpy as np
    from scipy.optimize import minimize
    rng = np.random.default_rng(seed)

    def lhs_thm1(v):
        a, b, c, d, A, B, C, D = v
        R = a*b + c*d + A*B + C*D - 2*math.sqrt(b*c*A*D) - 2*math.sqrt(a*d*B*C)
        s = b + c - B - C
        t = a + b - A - B
        return R + s*s + abs(s*t)

    # random points on {nonneg, equal sums}
    worst = float("inf")
    for _ in range(200000):
        lo = rng.dirichlet(np.ones(4) * rng.choice([0.2, 1.0, 5.0]))
        up = rng.dirichlet(np.ones(4) * rng.choice([0.2, 1.0, 5.0]))
        worst = min(worst, lhs_thm1(np.concatenate([lo, up])))
    check(f"Thm 1 random sampling: 200000 points (seed {seed}), min of R+s^2+|st| = {worst:.3e} >= -1e-12",
          worst >= -1e-12, level="NUMERIC")

    # local minimization from random starts (parametrize simplices by squares/normalization)
    def obj(z):
        lo = z[:4]**2
        up = z[4:]**2
        lo = lo / lo.sum()
        up = up / up.sum()
        return lhs_thm1(np.concatenate([lo, up]))
    worst2 = float("inf")
    for _ in range(300):
        r = minimize(obj, rng.normal(size=8), method="Nelder-Mead",
                     options={"xatol": 1e-12, "fatol": 1e-14, "maxiter": 4000})
        worst2 = min(worst2, r.fun)
    check(f"Thm 1 Nelder-Mead, 300 starts (seed {seed}): min = {worst2:.3e} >= -1e-12",
          worst2 >= -1e-12, level="NUMERIC")

    # Theorem 2 on SO(3) x SO(3)
    def rot(q):
        w, x, y, z = q / np.linalg.norm(q)
        return np.array([[1-2*(y*y+z*z), 2*(x*y-w*z), 2*(x*z+w*y)],
                         [2*(x*y+w*z), 1-2*(x*x+z*z), 2*(y*z-w*x)],
                         [2*(x*z-w*y), 2*(y*z+w*x), 1-2*(x*x+y*y)]])

    def g_fun(z, with_abs=True):
        P = rot(z[:4])
        Q = rot(z[4:])
        d0 = P[0, 0] - Q[0, 0]
        e0 = P[2, 2] - Q[2, 2]
        G = 2*d0*d0 + P[0, 1]**2 + Q[0, 1]**2 + P[1, 0]**2 + Q[1, 0]**2 + 2*(P[2, 1]*Q[2, 1] - P[1, 2]*Q[1, 2])
        return G + (2*abs(d0*e0) if with_abs else 0.0)

    worst3 = min(g_fun(rng.normal(size=8)) for _ in range(200000))
    check(f"Thm 2 random (P,Q) in SO(3)^2: 200000 samples (seed {seed}), min G+2|d0e0| = {worst3:.3e} >= -1e-12",
          worst3 >= -1e-12, level="NUMERIC")
    worst4 = float("inf")
    for _ in range(300):
        r = minimize(g_fun, rng.normal(size=8), method="Nelder-Mead",
                     options={"xatol": 1e-12, "fatol": 1e-14, "maxiter": 4000})
        worst4 = min(worst4, r.fun)
    check(f"Thm 2 Nelder-Mead, 300 starts: min G+2|d0e0| = {worst4:.3e} >= -1e-10",
          worst4 >= -1e-10, level="NUMERIC")
    # informational: the |d0 e0| term is needed (G alone can be negative)
    worst5 = min(g_fun(rng.normal(size=8), with_abs=False) for _ in range(50000))
    print(f"[INFO] min of G alone over 50000 samples = {worst5:.4f} (negative means the compensation term is necessary)")
    # random SO(3) via QR has a quaternion representation: check reconstruction
    maxerr = 0.0
    for _ in range(2000):
        Mq, _r = np.linalg.qr(rng.normal(size=(3, 3)))
        if np.linalg.det(Mq) < 0:
            Mq[:, 0] *= -1
        # Shepperd's method
        tr = np.trace(Mq)
        cand = [1 + tr, 1 + 2*Mq[0, 0] - tr, 1 + 2*Mq[1, 1] - tr, 1 + 2*Mq[2, 2] - tr]
        k = int(np.argmax(cand))
        if k == 0:
            w = math.sqrt(cand[0]) / 2
            q = [w, (Mq[2, 1]-Mq[1, 2])/(4*w), (Mq[0, 2]-Mq[2, 0])/(4*w), (Mq[1, 0]-Mq[0, 1])/(4*w)]
        elif k == 1:
            x = math.sqrt(cand[1]) / 2
            q = [(Mq[2, 1]-Mq[1, 2])/(4*x), x, (Mq[0, 1]+Mq[1, 0])/(4*x), (Mq[0, 2]+Mq[2, 0])/(4*x)]
        elif k == 2:
            y = math.sqrt(cand[2]) / 2
            q = [(Mq[0, 2]-Mq[2, 0])/(4*y), (Mq[0, 1]+Mq[1, 0])/(4*y), y, (Mq[1, 2]+Mq[2, 1])/(4*y)]
        else:
            z = math.sqrt(cand[3]) / 2
            q = [(Mq[1, 0]-Mq[0, 1])/(4*z), (Mq[0, 2]+Mq[2, 0])/(4*z), (Mq[1, 2]+Mq[2, 1])/(4*z), z]
        maxerr = max(maxerr, np.abs(rot(np.array(q)) - Mq).max())
    check(f"Thm 2: 2000 random SO(3) matrices recovered from a unit quaternion in the doc's convention, max err {maxerr:.1e}",
          maxerr < 1e-9, level="NUMERIC")


def run_all():
    exact_checks()
    numeric_checks()
    print(f"EXACT failures: {len(FAIL)}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(run_all())
