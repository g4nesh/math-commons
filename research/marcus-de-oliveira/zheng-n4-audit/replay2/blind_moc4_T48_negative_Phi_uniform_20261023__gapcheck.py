"""Gap check (second pass) for blind_moc4_T48_negative_Phi_uniform_20261023, gap 1 (p.2, rank charts and
signed quadratic classes).

Own code only. The source vector m = q f7 + p f5 + j f29 + c f39 + e f2 + d f40 + b f41 + a f0 is built from the
paper's f_j (shared lib ../lib/moc4lib.py, paper numbering), the cubic is the paper's (A.2), w = m + h*kappa.
The only supplement input is the polynomial G (and C, Q for comparison), read AS DATA from
blind_moc4_T48_native_full_radial_20261023.json.
Usage:  MOC4_SUPP=<supplement root> python -I -B this.py PART     with PART in {A, B, C}
  A: F(w)=C+hQ, P (12-principal deleting 0,5), eq. (3) and H6, all gcds, irreducibility of G,
     square-class tests at 60 random rational specializations, numeric sign tests at real roots
  B: eq. (2) det H(w) = -4G^2/Q^2 at h=-C/Q, and the three 13-minors divisible by G      (about 250 s)
  C: complete symbolic proof that P*H6 is a square in the function field of G            (about 230 s)
Exit status nonzero iff an EXACT check fails.
"""
import os, sys, json, time, random, math
from fractions import Fraction
import flint, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lib'))
import moc4lib as L
SUPP = os.environ['MOC4_SUPP']
JS = os.path.join(SUPP, 'proofs', 'manin_transfer_20261021', 'blind_moc4_T48_native_full_radial_20261023.json')
PART = sys.argv[1]
fails = []
def check(name, cond):
    print(('PASS ' if cond else 'FAIL ') + name, flush=True); (None if cond else fails.append(name))
VARS = ('q', 'p', 'j', 'c', 'e', 'd', 'b', 'a', 'h')
CTX = flint.fmpq_mpoly_ctx.get(VARS, 'degrevlex')
q, p, j, c, e, d, b, a, h = CTX.gens()
ZERO = 0*q
ORDER = [7, 5, 29, 39, 2, 40, 41, 0]
FV = L.F_VALUES
add = lambda *ks: tuple(sum(FV[k][s] for k in ks) for s in range(24))
check('profile identities f40+f1=f0+f7+f39, f40+f9=f2+f7+f41', add(40, 1) == add(0, 7, 39) and add(40, 9) == add(2, 7, 41))
check('generators of T={0,3,5,6,9,12} = {0,1,2,5,7,9,29,39,40,41}',
      sorted(L.generators(frozenset([0, 3, 5, 6, 9, 12]))) == [0, 1, 2, 5, 7, 9, 29, 39, 40, 41])
m = [ZERO]*14
for k, cf in zip(ORDER, [q, p, j, c, e, d, b, a]):
    for i, x in enumerate(L.F_VECS[k]):
        x = Fraction(x)
        if x: m[i] = m[i]+cf*flint.fmpq(x.numerator, x.denominator)
w = list(m); w[13] = w[13]+h
ms, Cs = L.cubic_sympy(); Hs = sp.hessian(Cs, ms)
def lin(expr):
    out = ZERO
    for mon, cf in sp.Poly(expr, *ms).terms():
        t = flint.fmpq(int(sp.numer(cf)), int(sp.denom(cf)))
        for i, e_ in enumerate(mon):
            if e_: t = t*w[i]**e_
        out = out+t
    return out
H = [[lin(Hs[i, k]) for k in range(14)] for i in range(14)]
Fw = lin(Cs).to_dict()
C0 = CTX.from_dict({k: v for k, v in Fw.items() if k[8] == 0})
Qh = CTX.from_dict({k[:8]+(0,): v for k, v in Fw.items() if k[8] == 1})
check('F(w) has degree 1 in h, Q = (c+d)(d+e+j) - p(e+j)', max(k[8] for k in Fw) == 1 and Qh == (c+d)*(d+e+j)-p*(e+j))
def load(terms):
    dd = {}
    for ex, cf in terms:
        assert ex[8] == 0
        fr = Fraction(cf); dd[tuple(ex[:8])+(0,)] = flint.fmpq(fr.numerator, fr.denominator)
    return CTX.from_dict(dd)
dat = json.load(open(JS))
check('supplement C and Q equal the rebuilt ones', load(dat['C_terms']) == C0 and load(dat['Q_terms']) == Qh)
G = [load(f['terms']) for f in dat['radial_factors'] if f['degree'] == 9][0]
check('G: 1772 terms, total degree 9', len(G.to_dict()) == 1772 and G.total_degree() == 9)
def bareiss(M):
    M = [row[:] for row in M]; n = len(M); sign = 1; prev = None
    for k in range(n-1):
        piv = next((i for i in range(k, n) if not M[i][k].is_zero()), None)
        if piv is None: return ZERO
        if piv != k: M[k], M[piv] = M[piv], M[k]; sign = -sign
        for i in range(k+1, n):
            for jj in range(k+1, n):
                t = M[i][jj]*M[k][k]-M[i][k]*M[k][jj]
                M[i][jj] = t if prev is None else t/prev
            M[i][k] = ZERO
        prev = M[k][k]
    return sign*M[n-1][n-1]
def coeffs(F, vi):
    out = {}
    for k, v in F.to_dict().items():
        kk = list(k); e_ = kk[vi]; kk[vi] = 0; out.setdefault(e_, {})[tuple(kk)] = v
    return {e_: CTX.from_dict(dd) for e_, dd in out.items()}
def forced(F):
    cs = coeffs(F, 8); n = max(cs)
    return sum((cs.get(k, ZERO))*(-C0)**k*Qh**(n-k) for k in range(n+1)), n
idx = [i for i in range(14) if i not in (0, 5)]
P = bareiss([[H[i][k] for k in idx] for i in idx])
check('P (12-principal deleting 0,5) is height-free with 2784 terms', P.degrees()[8] == 0 and len(P.to_dict()) == 2784)
gc = coeffs(G, 2); A, B, Cg = gc[2], gc.get(1, ZERO), gc[0]
pc = coeffs(P, 2); P2, P1, P0 = pc.get(2, ZERO), pc.get(1, ZERO), pc[0]
check('deg_j G = 2 and deg_j P = 2', max(gc) == 2 and max(pc) == 2)
Y = b+c+d; V = a+d+e; T = q*V-(a+d)*(d+e)
Phi = ((q+d)*V*b+p*T)**2+4*d*p**2*(a+d)*(q-e)*V
D = B**2-4*A*Cg
H6 = D/(Y**2*(q+d)**2*Phi)
check('(3) B^2-4AC0 = (b+c+d)^2 (q+d)^2 Phi H6 with H6 a j-free polynomial of 159 terms',
      Y**2*(q+d)**2*Phi*H6 == D and H6.degrees()[2] == 0 and len(H6.to_dict()) == 159)
if PART == 'A':
    check('content gcd(A,B,C0) = 1 (G has no j-free factor)', A.gcd(B).gcd(Cg) == 1)
    check('gcd(G,P) = gcd(G,G_j) = gcd(G,H6) = gcd(G,Phi) = 1',
          G.gcd(P) == 1 and G.gcd(G.derivative('j')) == 1 and G.gcd(H6) == 1 and G.gcd(Phi) == 1)
    fl = G.factor()[1]
    check('G is irreducible over Q', len(fl) == 1 and fl[0][1] == 1)
    r1 = A*P1-B*P2; r0 = A*P0-Cg*P2
    def ev(F, pt):
        v = F.subs({k: flint.fmpq(x.numerator, x.denominator) for k, x in pt.items()}).to_dict()
        val = v.get((0,)*9, 0); return Fraction(int(val.p), int(val.q)) if hasattr(val, 'p') else Fraction(str(val))
    def issq(x):
        x = Fraction(x)
        return x >= 0 and math.isqrt(x.numerator)**2 == x.numerator and math.isqrt(x.denominator)**2 == x.denominator
    random.seed(1); names = ['q', 'p', 'c', 'e', 'd', 'b', 'a']; okc = badc = 0
    for trial in range(60):
        pt = {n_: Fraction(random.randint(1, 40), random.randint(1, 12)) for n_ in names}
        Av, Bv, Cv = ev(A, pt), ev(B, pt), ev(Cg, pt); Dv = Bv*Bv-4*Av*Cv
        if Av == 0 or Dv == 0 or issq(Dv): continue
        for fac in (ev(H6, pt), ev(Phi, pt)):
            U = 2*fac*(2*Av*ev(r0, pt)-Bv*ev(r1, pt)); Vv = 2*fac*ev(r1, pt); n2 = U*U-Dv*Vv*Vv
            good = False
            if issq(n2):
                nn = Fraction(math.isqrt(n2.numerator), math.isqrt(n2.denominator))
                good = issq((U+nn)/2) or issq((U-nn)/2)
            okc += good; badc += (not good)
    check('specializations: P*H6 and P*Phi are squares in Q(sqrt D) at %d random rational points' % okc, badc == 0 and okc > 0)
    import mpmath as mp
    mp.mp.dps = 40; cnt = neg = viol = 0; negok = True
    for trial in range(400):
        pt = {n_: Fraction(random.randint(1, 60), random.randint(1, 15)) for n_ in names}
        if trial % 2 == 0: pt['e'] = pt['q']+Fraction(random.randint(1, 40), random.randint(1, 10))
        Av, Bv, Cv = ev(A, pt), ev(B, pt), ev(Cg, pt); Dv = Bv*Bv-4*Av*Cv
        if Dv <= 0 or Av == 0: continue
        H6v, Phv = ev(H6, pt), ev(Phi, pt); co = [ev(P2, pt), ev(P1, pt), ev(P0, pt)]
        f = lambda x: mp.mpf(x.numerator)/x.denominator
        for sg in (1, -1):
            jr = (-f(Bv)+sg*mp.sqrt(f(Dv)))/(2*f(Av)); Pv = f(co[0])*jr**2+f(co[1])*jr+f(co[2]); cnt += 1
            if Pv*H6v < -mp.mpf(10)**-25*(abs(Pv*H6v)+1): viol += 1
            if Phv < 0:
                neg += 1; negok = negok and Pv < 0
    print('NUMERIC real roots of G in j (seed 1, 40 digits): %d points, %d with Phi<0, %d violations of P*H6>=0, '
          'P<0 at every Phi<0 point: %s' % (cnt, neg, viol, negok))
elif PART == 'B':
    t0 = time.time()
    num, n = forced(bareiss(H))
    check('(2) det H(w) at h=-C/Q equals -4 G^2/Q^2', num == -4*G**2*Qh**(n-2))
    for name, (dr, dc), cnt in [('(0,0)', (0, 0), 872), ('(5,5)', (5, 5), 421), ('(0,5)', (0, 5), 613)]:
        rows = [i for i in range(14) if i != dr]; cols = [k for k in range(14) if k != dc]
        nm, _ = forced(bareiss([[H[i][k] for k in cols] for i in rows]))
        try: quo = nm/G; ok = True
        except Exception: ok = False
        check('13-minor deleting row/col %s: forced numerator divisible by G, quotient %d terms' % (name, cnt),
              ok and len(quo.to_dict()) == cnt)
    print('part B time %.0f s' % (time.time()-t0))
elif PART == 'C':
    t0 = time.time()
    r1 = A*P1-B*P2; r0 = A*P0-Cg*P2; u0 = 2*A*r0-B*r1
    Res = (A*r0**2-B*r0*r1+Cg*r1**2)/A          # = Res_j(G,P)
    def sq(F):
        try: return F.sqrt()
        except Exception: return None
    rho = sq(Res)
    check('Res_j(G,P) is a perfect square', rho is not None)
    hits = []
    if rho is not None:
        for e_ in (1, -1):
            X = H6*(u0+2*e_*A*rho)
            hits.append((e_, sq(X) is not None, sq(D*X) is not None))
    print('  (sign, H6(u0+-2A rho) square, D*H6(u0+-2A rho) square):', hits)
    check('P*H6 is a square in the function field Q(q,p,c,e,d,b,a)[j]/(G)', any(x[1] or x[2] for x in hits))
    print('part C time %.0f s' % (time.time()-t0))
print('FAILED:', fails if fails else 'none'); sys.exit(1 if fails else 0)
