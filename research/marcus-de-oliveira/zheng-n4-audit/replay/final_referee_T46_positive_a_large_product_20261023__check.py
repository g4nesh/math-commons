"""Exact replay for final_referee_T46_positive_a_large_product_20261023 (historical T46 = paper T_47).

Own code only; nothing from the supplement is read or executed (the document's 82-square list was
NOT used: for the one mixed-sign table an independent AM-GM certificate is built and verified here).
G9 is rebuilt from the paper's definitions by the sibling replay of the negative_j_Heron note. Our
clearing of denominators is not the document's minimal one, so our term counts differ from the
document's table; sign conclusions are unaffected because every removed or retained factor has
positive coefficients. Usage: python -I -B this.py [scalar] [A] [B] [C]  (A: q>b tables, B: b>q
tables except the last, C: last table + certificate). Exit nonzero iff an EXACT check fails.
"""
import os, sys, itertools
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import final_referee_T46_negative_j_Heron_uniform_20261023__check as base
import final_referee_T46_positive_a_small_or_balanced_20261023__check as ch
import sympy as sp
FAIL = []
def check(name, ok, level='EXACT'):
    print(f"[{level}] {'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok and level == 'EXACT': FAIL.append(name)

def amgm_certificate(Pd, names, optional):
    """Pd: dict monomial -> int coefficient of a polynomial claimed negative. Build -P = R + sum w (X^N1 - 2X^M + X^N2),
    M = (N1+N2)/2, w > 0, R >= 0 coefficientwise; verify exactly. AM-GM makes each bracket >= 0 for X >= 0."""
    negP = {m: -int(c) for m, c in Pd.items()}
    bad = {m: -c for m, c in negP.items() if c < 0}
    avail = {m: Fr(c) for m, c in negP.items() if c > 0}
    nv = len(next(iter(Pd)))
    offs = [o for o in itertools.product(range(-3, 4), repeat=nv) if any(o) and o > tuple(-x for x in o)]
    cert = []; covered = True
    for M, beta in sorted(bad.items(), key=lambda kv: -kv[1]):
        rem = Fr(beta, 2); cands = []
        for o in offs:
            N1 = tuple(m+x for m, x in zip(M, o)); N2 = tuple(m-x for m, x in zip(M, o))
            if min(N1) >= 0 and min(N2) >= 0 and N1 in avail and N2 in avail: cands.append((min(avail[N1], avail[N2]), N1, N2))
        for cap, N1, N2 in sorted(cands, reverse=True):
            if rem <= 0: break
            w = min(rem, avail[N1], avail[N2])
            if w <= 0: continue
            avail[N1] -= w; avail[N2] -= w; rem -= w; cert.append((M, N1, N2, w))
        covered &= rem <= 0
    R = {m: Fr(c) for m, c in negP.items()}
    for M, N1, N2, w in cert:
        assert all(2*x == y + z for x, y, z in zip(M, N1, N2))
        R[N1] -= w; R[N2] -= w; R[M] += 2*w
    opt = [names.index(o) for o in optional]; jx = names.index('J')
    strict = [m for m, v in R.items() if v > 0 and all(m[i] == 0 for i in opt)]
    return len(bad), len(cert), covered, all(v >= 0 for v in R.values()), len(strict), sum(1 for m in strict if m[jx] == 0)

def tables(which):
    R = base.g9(); G = R['G']
    check('G9 rebuilt (1556 terms; det Hess = -4G^2/Q^2)', R['ok_sq'] and len(G) == 1556)
    GD = {tuple(int(x) for x in m): int(c) for m, c in G.to_dict().items()}; DEG = [int(x) for x in G.degrees()]
    V, Ap, J, S, C, E, H = sp.symbols('V Ap J S C E H', positive=True)
    vl = (V, Ap, J, S, C, E, H); names = [str(v) for v in vl]
    r = 1 + V; big = r + Ap; small = r**2/(r + Ap)
    doc = {('q>b','low'): (196914, 65863, 422404), ('q>b','middle'): (524880, 176835, 1066536),
           ('b>q','low'): (178890, 59765, 344636), ('b>q','middle'): (462252, 155158, 869439)}
    for orient, (q, b) in (('q>b', (big, small)), ('b>q', (small, big))):
        lam = (q + b + 2)/(q + b - 2); K = (q + 1)*(b + 1)/(q*b - 1)
        for rng in ('low', 'middle'):
            if rng == 'low': j = lam**2*J/(K*(1+J)); cap = 1/j - 1/lam
            else: j = (lam + K*J)**2/(K*(1+J)**2); cap = (K - lam)/(lam + K*J)**2
            p = cap*S/(1+S); e0 = 1 + 1/b
            a4 = lambda ee: (b + ee + 1)/(b*(ee - 1) - 1)
            for i, (nm, ee, aa) in enumerate((('e<e0', e0*E/(1+E), H), ('e=e0', e0, H), ('e>e0', e0 + E, a4(e0 + E)*H/(1+H)))):
                last = (orient, rng, nm) == ('b>q', 'middle', 'e>e0')
                part = 'A' if orient == 'q>b' else ('C' if last else 'B')
                if part not in which: continue
                tot, dp, _ = ch.substitute(GD, DEG, dict(q=q, p=p, j=j, c=C, e=ee, d=1, b=b, a=aa), vl)
                n, npos, nneg, fp, fn = ch.sign_report(tot, names, ('C', 'H'))
                d = tot.to_dict()
                jfree = [int(c) for m, c in d.items() if m[2] == 0]
                label = '%s %s %s (our %d terms; document %d)' % (orient, rng, nm, n, doc[(orient, rng)][i])
                if not last:
                    check(label + ': all coefficients negative, denominators positive, strict negative monomial free of C,H'
                          + ('' if rng == 'low' else ', and one with J-exponent 0'),
                          npos == 0 and dp and fn > 0 and (rng == 'low' or sum(1 for c in jfree if c < 0) > 0))
                else:
                    check(label + ': %d positive coefficients (document: 80 in its clearing)' % npos, dp and npos > 0, level='INFO')
                    nb, nc, cov, Rnn, nstrict, nstrictJ = amgm_certificate({tuple(int(x) for x in m): int(c) for m, c in d.items()}, names, ('C', 'H'))
                    check(label + ': auditor AM-GM certificate, %d deficits covered by %d pairs, remainder >= 0, %d strict monomials free of C,H (%d with J^0)'
                          % (nb, nc, nstrict, nstrictJ), cov and Rnn and nstrict > 0 and nstrictJ > 0)

def scalar():
    q, b, j, t, JJ = sp.symbols('q b j t J', positive=True)
    lam = (q + b + 2)/(q + b - 2); K = (q + 1)*(b + 1)/(q*b - 1)
    check('K - lambda = (q-b)^2/((qb-1)(q+b-2))', sp.simplify(K - lam - (q - b)**2/((q*b - 1)*(q + b - 2))) == 0)
    s = q + b - 2; h = 0
    thr = ((j - t**2)*h**2 - 4*t*h*(t - 1) + s*(t - 1)**2)/(s*j - (s + 4)*t**2)
    check('core Thm 1 at beta = -1: threshold = (t-1)^2/(j - lambda t^2)', sp.simplify(thr - (t - 1)**2/(j - lam*t**2)) == 0)
    check('core Thm 1 domain at beta = -1: (1 + (q+b-2beta)/l) = K', sp.simplify(1 + (q + b + 2)/(q*b - 1) - K) == 0)
    f = (t - 1)**2/(j - lam*t**2)
    check('t = j/lambda is stationary and gives 1/j - 1/lambda', sp.simplify(sp.diff(f, t).subs(t, j/lam)) == 0 and sp.simplify(f.subs(t, j/lam) - (1/j - 1/lam)) == 0)
    check('t = sqrt(j/K) gives (1 - sqrt(j/K))^2 / (j(1 - lambda/K))', sp.simplify(f.subs(t, sp.sqrt(j/K)) - (1 - sp.sqrt(j/K))**2/(j*(1 - lam/K))) == 0)
    jm = (lam + K*JJ)**2/(K*(1 + JJ)**2)
    capm = (1 - sp.sqrt(jm/K))**2/(jm*(1 - lam/K))
    # sqrt(jm/K) = (lam + K J)/(K (1+J)) for positive quantities
    capm2 = capm.subs(sp.sqrt(jm/K), (lam + K*JJ)/(K*(1 + JJ)))
    check('middle chart cap = (K - lambda)/(lambda + K J)^2 (with sqrt(j/K) = (lambda+KJ)/(K(1+J)))',
          sp.simplify(((1 - (lam + K*JJ)/(K*(1 + JJ)))**2/(jm*(1 - lam/K))) - (K - lam)/(lam + K*JJ)**2) == 0)
    check('low chart j = lambda^2 J/(K(1+J)) ranges over (0, lambda^2/K); middle j(J=0) = lambda^2/K, j(J->oo) = K',
          sp.simplify(jm.subs(JJ, 0) - lam**2/K) == 0 and sp.simplify(sp.limit(jm.subs({q: 3, b: 2}), JJ, sp.oo) - K.subs({q: 3, b: 2})) == 0)

if __name__ == '__main__':
    which = sys.argv[1:] or ['scalar', 'A', 'B', 'C']
    if 'scalar' in which: scalar()
    tables([w for w in which if w in 'ABC'])
    print('summary: %d EXACT failures' % len(FAIL))
    sys.exit(1 if FAIL else 0)
