"""Exact replay for final_referee_T46_positive_a_small_or_balanced_20261023 (historical T46 = paper T_47).

Own code only; nothing from the supplement is read or executed. G9 is rebuilt from the paper's
definitions by the sibling replay final_referee_T46_negative_j_Heron_uniform_20261023__check.py
(det Hess C at the restored height = -4 G^2/Q^2). The chart substitutions of the document are
applied with d = 1; each variable's denominator is cleared to the power deg_x G and positive-
coefficient denominators are then divided out while they divide. Every printed line is one check;
exit status is nonzero iff an EXACT check fails.
"""
import os, sys, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import final_referee_T46_negative_j_Heron_uniform_20261023__check as base
import flint, sympy as sp
FAIL = []
def check(name, ok, level='EXACT'):
    print(f"[{level}] {'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok and level == 'EXACT': FAIL.append(name)

def substitute(GD, DEG, subs, newvars):
    nctx = flint.fmpz_mpoly_ctx.get(tuple(str(v) for v in newvars), 'degrevlex')
    def tof(expr):
        P = sp.Poly(sp.expand(expr), *newvars)
        return nctx.from_dict({m: int(cf) for m, cf in P.terms()})
    ND = []; dens_pos = True
    for nm in base.NAMES[:8]:
        num, den = sp.fraction(sp.together(sp.sympify(subs[nm])))
        Nf, Df = tof(num), tof(den)
        if any(int(x) < 0 for x in Df.coeffs()): Nf, Df = -Nf, -Df
        dens_pos &= all(int(x) > 0 for x in Df.coeffs())
        ND.append((Nf, Df))
    cache = {}
    def pw(k, w, e):
        if (k, w, e) not in cache: cache[(k, w, e)] = ND[k][w]**e
        return cache[(k, w, e)]
    tot = nctx.from_dict({})
    for mon, cf in GD.items():
        t = nctx.from_dict({(0,)*len(newvars): int(cf)})
        for k in range(8): t *= pw(k, 0, mon[k])*pw(k, 1, DEG[k]-mon[k])
        tot += t
    for k in range(8):
        Df = ND[k][1]
        if Df.is_constant(): continue
        while True:
            qq, r = divmod(tot, Df)
            if r.is_zero(): tot = qq
            else: break
    ct = tot.content()
    if ct != 0: tot = tot/ct
    return tot, dens_pos, nctx

def sign_report(tot, names, optional):
    d = tot.to_dict(); opt = [names.index(o) for o in optional]
    npos = sum(1 for c in d.values() if int(c) > 0); nneg = len(d) - npos
    free = [int(c) for m, c in d.items() if all(m[i] == 0 for i in opt)]
    return len(d), npos, nneg, sum(1 for c in free if c > 0), sum(1 for c in free if c < 0)

if __name__ == '__main__':
    R = base.g9(); G = R['G']
    check('G9 rebuilt (1556 terms; det Hess = -4G^2/Q^2)', R['ok_sq'] and len(G) == 1556)
    GD = {tuple(int(x) for x in m): int(c) for m, c in G.to_dict().items()}; DEG = [int(x) for x in G.degrees()]
    q,p,j,c,e,d,b,a,h = base.CTX.gens()
    D_ = a+d; M = (c+d)*(a+d+e); Lr = M + b*D_
    C0 = -(q+d)*Lr - j*M
    check('C0 = -(q+d)L - jM has only negative coefficients (so C0 < 0 when d, b > 0 and the rest >= 0)',
          all(int(x) < 0 for x in C0.coeffs()))
    V,W,P,J,C,U,A = sp.symbols('V W P J C U A', positive=True)
    def echarts(bb):
        e0 = 1 + 1/bb; a4 = lambda ee: (bb + ee + 1)/(bb*(ee - 1) - 1)
        return [('e<e0', e0*U/(1+U), A), ('e=e0', e0, A), ('e>e0', e0 + U, a4(e0 + U)*A/(1+A))]
    def run(label, qq, bb, pp, jj, vl, expect=None):
        for i, (nm, ee, aa) in enumerate(echarts(bb)):
            tot, dp, _ = substitute(GD, DEG, dict(q=qq, p=pp, j=jj, c=C, e=ee, d=1, b=bb, a=aa), vl)
            n, npos, nneg, fp, fn = sign_report(tot, [str(v) for v in vl], ('A', 'C'))
            tag = '' if expect is None else ' (document: %d terms)' % expect[i]
            check('%s, %s: every coefficient negative, denominators positive, strict negative monomial free of A, C; %d terms%s'
                  % (label, nm, n, tag), npos == 0 and nneg == n and dp and fn > 0 and (expect is None or n == expect[i]))
    run('q,b<1', V/(1+V), W/(1+W), P, J, (V,W,P,J,C,U,A), expect=(4916, 1815, 10809))
    run('q=1,b<1 wall', 1, W/(1+W), P, J, (W,P,J,C,U,A))
    run('q<1,b=1 wall', V/(1+V), 1, P, J, (V,P,J,C,U,A))
    run('q=b=1 wall', 1, 1, P, J, (P,J,C,U,A))
    x = 1 + V; lam = (V+2)/V
    run('q=b>1 outside core', x, x, V*P/((V+2)*J*(1+P)), lam*J/(1+J), (V,P,J,C,U,A), expect=(4095, 1377, 8398))
    # chart coverage identities (scalar)
    vv, ww, uu, AA, JJ, PP = sp.symbols('v w u A J P', positive=True)
    bb = ww/(1+ww); e0 = 1 + 1/bb; ee = e0 + uu
    a4 = (bb + ee + 1)/(bb*(ee - 1) - 1)
    check('e > e0 iff b(e-1) - 1 > 0 (so a4 > 0 exactly on the e > e0 chart): b(e0+u-1)-1 = b u',
          sp.simplify(bb*(ee - 1) - 1 - bb*uu) == 0)
    xx = 1 + vv; lm = (xx + xx + 2)/(xx + xx - 2); Kq = (xx + 1)**2/(xx*xx - 1)
    check('balanced: lambda = K = (V+2)/V, core threshold max(1/j - (x-1)/(x+1), 0)',
          sp.simplify(lm - Kq) == 0 and sp.simplify(lm - (vv + 2)/vv) == 0)
    jj = lm*JJ/(1+JJ); cap = 1/jj - 1/lm
    check('balanced chart: p = V P/((V+2) J (1+P)) equals cap * P/(1+P) with cap = 1/j - 1/lambda',
          sp.simplify(cap*PP/(1+PP) - vv*PP/((vv+2)*JJ*(1+PP))) == 0)
    print('summary: %d EXACT failures' % len(FAIL))
    sys.exit(1 if FAIL else 0)
