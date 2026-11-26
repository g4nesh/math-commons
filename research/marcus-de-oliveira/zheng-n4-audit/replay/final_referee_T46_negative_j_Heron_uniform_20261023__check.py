"""Exact replay for final_referee_T46_negative_j_Heron_uniform_20261023 (historical T46 = paper T_47).

Own code only. Nothing from the supplement is read or executed. Everything is rebuilt from the
paper's definitions through the shared library moc4lib (cubic C of (A.2), generators f_j of App. D,
kappa = e_13). Every printed line is one check; exit status is nonzero iff an EXACT check fails.

G9 is computed here from scratch: u = sum coef*f_j + h*kappa for the eight canonical rays
(q,p,j,c,e,d,b,a) = (f27,f1,f59,f5,f2,f28,f6,f4); det Hess C(u) is computed by fraction-free
elimination over Z[q,...,a,h]; h is eliminated through C(u) = C0 + hQ.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lib'))
import moc4lib as L
import flint

FAIL = []
def check(name, ok, level='EXACT'):
    print(f"[{level}] {'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok and level == 'EXACT': FAIL.append(name)

NAMES = ('q','p','j','c','e','d','b','a','h')
CTX = flint.fmpz_mpoly_ctx.get(NAMES, 'degrevlex')
ROLES = dict(q=27,p=1,j=59,c=5,e=2,d=28,b=6,a=4)
Z = CTX.from_dict({})

def bareiss(A):
    n = len(A); A = [r[:] for r in A]; sign = 1; prev = None
    for k in range(n-1):
        best = None
        for i in range(k, n):
            for jj in range(k, n):
                x = A[i][jj]
                if not x.is_zero() and (best is None or len(x) < best[0]): best = (len(x), i, jj)
        if best is None: return Z
        _, pi, pj = best
        if pi != k: A[k], A[pi] = A[pi], A[k]; sign = -sign
        if pj != k:
            for row in A: row[k], row[pj] = row[pj], row[k]
            sign = -sign
        piv = A[k][k]
        for i in range(k+1, n):
            for jj in range(k+1, n):
                num = piv*A[i][jj] - A[i][k]*A[k][jj]
                A[i][jj] = num if prev is None else num / prev
        prev = piv
    return sign*A[n-1][n-1]

def build():
    g = dict(zip(NAMES, CTX.gens()))
    u = []
    for i in range(14):
        s = Z
        for k, v in ROLES.items():
            fr = L.F_VECS[v][i]; assert fr.denominator == 1
            s += int(fr.numerator)*g[k]
        if i == 13: s += g['h']
        u.append(s)
    mctx = flint.fmpz_mpoly_ctx.get(tuple('m%d' % i for i in range(14)), 'degrevlex')
    Cm = L.cubic_C(list(mctx.gens()))
    def tolin(pm):
        s = Z
        for mon, cf in pm.to_dict().items():
            l = [i for i, x in enumerate(mon) if x]; s += int(cf)*u[l[0]]
        return s
    H = [[tolin(Cm.derivative(i).derivative(k)) for k in range(14)] for i in range(14)]
    return u, H

def coeffs_in(P, var):
    out = {}
    for mon, cf in P.to_dict().items():
        k = mon[var]; m2 = list(mon); m2[var] = 0
        out.setdefault(k, {})[tuple(m2)] = cf
    return {k: CTX.from_dict(v) for k, v in out.items()}

def g9():
    q,p,j,c,e,d,b,a,h = CTX.gens()
    u, H = build()
    D_ = a+d; M = (c+d)*(a+d+e); Lr = M + b*D_
    C0 = -(q+d)*Lr - j*M; Q = j*p - j*(c+e+d) - (c+d)*(e+d)
    ok_height = (L.cubic_C(u) == C0 + h*Q)
    det = bareiss(H)
    ck = coeffs_in(det, 8); hd = max(ck)
    R = sum((ck[k]*(-C0)**k*Q**(hd-k) for k in ck), Z)       # = Q^hd det Hess C(u) at h = -C0/Q
    qq, rem = divmod(R, Q**(hd-2))
    G = (-qq/4).sqrt() if rem.is_zero() else None
    ok_sq = G is not None and -4*G*G*Q**(hd-2) == R
    return dict(G=G, C0=C0, Q=Q, H=H, ok_height=ok_height, ok_sq=ok_sq, hd=hd)


def numeric_stratum(seeds, margin):
    """NUMERIC: points with contact set containing T (lib least-squares sampler), their canonical coordinates,
    relative size of G9, sign of H6 when d > 0, and the region tests of the two positive-a notes."""
    import numpy as np
    R = g9(); GD = {tuple(int(x) for x in m): int(c) for m, c in R['G'].to_dict().items()}
    eight = [27,1,59,5,2,28,6,4]
    Fm = np.array([[float(v) for v in L.F_VECS[k]] for k in eight]).T; Pn = L.P_numpy()
    stats = dict(points=0, dpos=0, worstG=0.0, minH6_dpos=None, region_viol=0)
    for seed in seeds:
        r = L.random_point_on_contact_stratum([1,3,5,6,7,8], seed=seed, tries=8, margin=margin)
        if r is None: continue
        c = np.linalg.lstsq(Pn, r['gaps'], rcond=None)[0]; x = np.linalg.lstsq(Fm, c, rcond=None)[0]
        nrm = sum(abs(cf)*np.prod([abs(x[k])**m[k] for k in range(8)]) for m, cf in GD.items())
        gv = sum(cf*np.prod([x[k]**m[k] for k in range(8)]) for m, cf in GD.items())
        stats['points'] += 1; stats['worstG'] = max(stats['worstG'], abs(gv)/max(nrm, 1e-300))
        if x[5] > 1e-6:
            stats['dpos'] += 1
            q,p,j,cc,e,d,b,a = x/x[5]
            D_ = a+1; M = (cc+1)*(a+1+e); Lr = M+b*D_; U = cc+e
            H6 = (Lr*(p+q) - b*D_*U)**2 - 4*p*q*M*Lr + 4*cc*e*b*D_*M
            stats['minH6_dpos'] = H6 if stats['minH6_dpos'] is None else min(stats['minH6_dpos'], H6)
            four = a*b*(e-1) >= a+b+e+1 - 1e-9
            core = None
            if q*b > 1 and abs(q-b) > 1e-9:
                lam = (q+b+2)/(q+b-2); K = (q+1)*(b+1)/(q*b-1)
                core = j >= K or (j < lam**2/K and p >= 1/j - 1/lam - 1e-9) or (j >= lam**2/K and p >= (1-np.sqrt(j/K))**2/(j*(1-lam/K)) - 1e-9)
            if (q <= 1 and b <= 1 and not four) or (core is False and not four): stats['region_viol'] += 1
    return stats

if __name__ == '__main__':
    if 'numeric' in sys.argv:
        for lo, hi, mg in ((0, 40, 0.0), (100, 140, 0.02)):
            st = numeric_stratum(range(lo, hi), mg)
            check('stratum sampler seeds %d-%d margin %g: %d points (%d with d>0), max |G|/sum|terms| = %.1e, min H6 at d>0 points = %s, region violations of the positive-a notes = %d' % (lo, hi-1, mg, st['points'], st['dpos'], st['worstG'], st['minH6_dpos'], st['region_viol']), st['region_viol'] == 0 and (st['minH6_dpos'] is None or st['minH6_dpos'] > -1e-9), level='NUMERIC')
        sys.exit(0)
    q,p,j,c,e,d,b,a,h = CTX.gens()
    T46 = frozenset([1,3,5,6,7,8])
    check('contact set {1,3,5,6,7,8} is paper T_47 (identity group element)', L.find_g(L.CONTACT_SETS[47], T46) == ((0,1,2,3),(0,1,2,3),0))
    check('generators(T) = {0,1,2,3,4,5,6,27,28,59}; the eight canonical rays are among them; only f28 non-minor',
          sorted(L.generators(T46)) == [0,1,2,3,4,5,6,27,28,59] and [k for k, v in ROLES.items() if not L.F_IS_MINOR[v]] == ['d'])
    R = g9(); G = R['G']
    check('C(sum coef f + h kappa) == C0 + h Q with C0 = -(q+d)L - jM, Q = jp - j(c+e+d) - (c+d)(e+d)', R['ok_height'])
    check('det Hess C at h = -C0/Q equals -4 G^2 / Q^2 for a polynomial G (deg_h det = %d)' % R['hd'], R['ok_sq'])
    check('G has 1556 terms, total degree 9, degree 2 in j', len(G) == 1556 and G.total_degree() == 9 and G.degrees()[2] == 2)
    cj = coeffs_in(G, 2); A, B, C = cj[2], cj[1], cj[0]
    D_ = a+d; M = (c+d)*(a+d+e); Lr = M + b*D_; U = c+e
    H6 = (Lr*(p+q) - b*D_*U)**2 - 4*p*q*M*Lr + 4*c*e*b*D_*M
    check('H6 (printed formula) has 150 terms and is irreducible over Q', len(H6) == 150 and len(H6.factor()[1]) == 1 and H6.factor()[1][0][1] == 1)
    qq, rem = divmod(B*B - 4*A*C, (q+d)**2*H6); H8 = qq
    check('(1) B^2 - 4AC = (q+d)^2 H6 H8 with H8 a polynomial of 279 terms', rem.is_zero() and len(H8) == 279)
    W0 = (q+d)*Lr*p + b*((q-d)*Lr - a*e*(c+d)); Xi = q*(c+e+d+b) - c*e
    check('(2) H8 = W0^2 + 4 d D b^2 L Xi', H8 == W0**2 + 4*d*D_*b**2*Lr*Xi)
    check('(3) dH8/dp = 2 W0 (q+d) L when the formula (2) is differentiated (L, D, Xi free of p)', H8.derivative(1) == 2*W0*(q+d)*Lr)
    check('(3) at Xi = 0, W0 = 0: dH8/dq = 4 d D b^2 L (c+e+d+b) (identity dH8/dq = 2 W0 W0_q + 4dDb^2 L Xi_q)',
          H8.derivative(0) == 2*W0*W0.derivative(0) + 4*d*D_*b**2*Lr*(c+e+d+b))
    check('(4) gcd(A,B,C) = 1', A.gcd(B).gcd(C).is_one())
    Gj = G.derivative(2)
    check('(4) gcd(G, G_j) = 1', G.gcd(Gj).is_one())
    check('(5) G_j^2 - 4 A G = (q+d)^2 H6 H8', Gj*Gj - 4*A*G == (q+d)**2*H6*H8)
    Hm = R['H']; idx = [i for i in range(14) if i not in (0, 1)]
    Pm = bareiss([[Hm[i][k] for k in idx] for i in idx])
    P10, rem = divmod(-Pm, j*j)
    check('principal 12-minor deleting coordinates 0,1 is height-free and equals -j^2 P10, P10 quadratic in j',
          Pm.degrees()[8] == 0 and rem.is_zero() and P10.degrees()[2] == 2)
    cp = coeffs_in(P10, 2); al, be, ga = cp[2], cp[1], cp[0]
    u1 = A*be - al*B; u0 = A*ga - al*C
    W = (u0*u0 - u1*(B*ga - be*C)).sqrt()
    check('u0^2 - u1(B gamma - beta C) = W^2 for a polynomial W independent of j', W is not None and W.degrees()[2] == 0)
    Sp = 2*A*u0 - B*u1 + 2*A*W; Sm = 2*A*u0 - B*u1 - 2*A*W
    if not divmod(-Sp, H8)[1].is_zero(): W = -W; Sp, Sm = Sm, Sp
    R8 = divmod(-Sp, H8); T8 = divmod(-Sm, (q+d)**2*H6)
    R8r = R8[0].sqrt() if R8[1].is_zero() else None; T8r = T8[0].sqrt() if T8[1].is_zero() else None
    check('S+ = 2Au0 - Bu1 + 2AW = -H8 R8^2 (sign of W fixed by this), R8 independent of j', R8r is not None and R8r.degrees()[2] == 0)
    check('S- = 2Au0 - Bu1 - 2AW = -(q+d)^2 H6 T8^2, T8 independent of j', T8r is not None and T8r.degrees()[2] == 0)
    check('P10 S+ == (u1 j + u0 + W)^2 mod G (exact division)', divmod(P10*Sp - (u1*j + u0 + W)**2, G)[1].is_zero())
    check('P10 S- == (u1 j + u0 - W)^2 mod G (exact division)', divmod(P10*Sm - (u1*j + u0 - W)**2, G)[1].is_zero())
    print('summary: %d EXACT failures' % len(FAIL))
    sys.exit(1 if FAIL else 0)
