"""Gap check (second pass) for independent_true_T48_c0_native_Heron_20261023, gap 1 (p.2 eqs (3)-(5)).

Rebuilds from the definitions in the companion note independent_true_T48_native_harmonic_conditional
(chart B=(t,0,1,inf), A=(y,z,w,x), eq. (4); J, W, M, ZD) the 24 normalized gaps H_sigma, the last-contact
curve, and the short amplitudes; loads the supplement's affine c amplitude AS DATA from
independent_true_T48_native_harmonic_conditional_exact_20261023.json and checks that it is a valid
representative (differs from H_17 by a function of (r,s,u,h) on the contact surface); then builds the three
real rows (3), and replays (4), (5), (6), (7), (9), (10) exactly with sympy.
Usage: python -I -B this.py   (supplement root from env MOC4_SUPP).  Nothing from the supplement is executed.
"""
import os, sys, json, itertools
import sympy as sp
SUPP = os.environ.get('MOC4_SUPP') or sys.argv[1]
JS = os.path.join(SUPP, 'proofs', 'manin_transfer_20261021',
                  'independent_true_T48_native_harmonic_conditional_exact_20261023.json')
fails = []
def check(name, cond):
    print(('PASS ' if cond else 'FAIL ') + name); (None if cond else fails.append(name))
PERMS = sorted(itertools.permutations(range(4)))
z, w, r, s, u, h, J, W = sp.symbols('z w r s u h J W')
# --- chart (harmonic note eq. (4))
V = r*(z-w)+s*(w-1); Dt = -V-u*z*(w-1)+w*z; Dy = V+(w-z)*(w-u*(w-1))
x = s*z*(w-1)/V; t = w*z/Dt; y = w*(V+w-z)/Dy
A = [y, z, w, x]; B = [t, 0, 1, None]
def F(sig):
    pr = 1
    for i in range(4):
        if B[sig[i]] is not None: pr *= (A[i]-B[sig[i]])
    return pr
Fs = [F(sg) for sg in PERMS]
E = Fs[0]-Fs[6]
H = [sp.cancel((f-Fs[0])/E) for f in Fs]
check('contact values H0=0,H3=s-r,H5=-r,H6=-1,H12=-u',
      [sp.simplify(H[k]-v) for k, v in [(0, 0), (3, s-r), (5, -r), (6, -1), (12, -u)]] == [0]*5)
d_ = r-s; M = d_*W+s; ZD = J*(W-1)-r
HJ = [sp.cancel(hh.subs({z: -M/ZD, w: W})) for hh in H]
R_ = d_+h
num, _ = sp.fraction(sp.together(HJ[9]-h)); Cx = sp.expand(num)
lead = sp.expand(-R_*J*((d_-1)*J+d_*(u-1)))
Cx = sp.expand(Cx*sp.cancel(lead/sp.expand(Cx.coeff(W, 2))))
check('[W^2]C, [J^2]C, C(J,0) equal the harmonic note eq. (6)',
      sp.expand(Cx.coeff(W, 2)-lead) == 0 and sp.expand(Cx.coeff(J, 2)+((d_-1)*W+s)*(R_*W-R_-1)) == 0
      and sp.expand(Cx.subs(W, 0)-s*(J+u)*((R_+1)*J+d_)) == 0)
# exclusive gaps: j=g2, b=g7, c=g17; c+d = g1-g7 (profiles of f7,f5,f29,f39,f2,f40,f41,f0)
Zs = J*(1-W)
check('j_amp = H_2 = J', sp.simplify(HJ[2]-J) == 0)
check('(c+d)_amp = H_1 - H_7 = (J-s)/Z', sp.simplify(HJ[1]-HJ[7]-(J-s)/Zs) == 0)
check('b_amp = H_7 = ((r-1)J+s)/Z-(r-s)', sp.simplify(HJ[7]-(((r-1)*J+s)/Zs-(r-s))) == 0)
H9 = HJ[9]
def onS(e): return sp.cancel(e.subs(h, H9))
check('(a+d)_amp = H19-H13 = -(u-1)(J-s)/J + (u-1)', sp.cancel(onS(HJ[19]-HJ[13]+(u-1)*(J-s)/J)-(u-1)) == 0)
check('(q+d)_amp = H14 = -(RZ+J+r-s)/r + (1-u) on the surface', sp.cancel(onS(HJ[14]+(R_*Zs+J+r-s)/r)-(1-u)) == 0)
# p: difference from the note's formula is a function of (r,s,u,h) on the surface (Jacobian criterion)
pd = onS(HJ[11]-(R_*(Zs+r)/(J+r-s)-(r+h-1)))
check('p_amp = H11 equals the note formula up to a real function of (r,s,u,h) (Jacobian test)',
      sp.cancel(sp.diff(pd, J)*sp.diff(H9, W)-sp.diff(pd, W)*sp.diff(H9, J)) == 0)
# --- supplied c representative (data)
dat = json.load(open(JS)); vs = [J, W, r, s, u, h]
def poly(terms): return sp.Add(*[sp.Rational(c)*sp.Mul(*[v**e for v, e in zip(vs, ex)]) for ex, c in terms])
cj = poly(dat['affine_amplitudes']['c']['num'])/poly(dat['affine_amplitudes']['c']['den'])
cd = onS(cj-HJ[17])
check('supplied c amplitude = H_17 + real function of (r,s,u,h) on the surface (Jacobian test)',
      sp.cancel(sp.diff(cd, J)*sp.diff(H9, W)-sp.diff(cd, W)*sp.diff(H9, J)) == 0)
check('supplied curve array equals the rebuilt curve C', sp.cancel(poly(dat['curve'])/Cx) == 1)
# --- rows (3)
Zc = sp.Symbol('Zc')
X, yy, zr, v, RR, KK, rr, ss = sp.symbols('x y zr v R K rr ss', real=True)
P = sp.cancel(sp.cancel(Cx.subs(W, 1-Zc/J)*J**2)/J)
c = sp.cancel(cj.subs(W, 1-Zc/J)); CD = Zc*(r*J-(r-s)*Zc); CN = sp.cancel(c*CD)
check('P(J,Z) and CN are polynomials', sp.fraction(sp.together(P))[1] == 1 and sp.fraction(sp.together(CN))[1] == 1)
sub = {J: X+sp.I*yy, Zc: zr+sp.I*yy*v, u: 1-KK/ss, h: RR-rr+ss, r: rr, s: ss}
def ex(e): return sp.expand(sp.expand(e).subs(sub, simultaneous=True))
Pe = ex(P); CNe = ex(CN); CDe = ex(CD)
rows = [sp.expand(sp.re(Pe)), sp.expand(sp.cancel(sp.im(Pe)/yy)),
        sp.expand(sp.cancel(sp.im(sp.expand(CNe*sp.conjugate(CDe)))/yy))]
Y = sp.Symbol('Y')
Mrows = []; brow = []; lin = True
for e in rows:
    nm, dn = sp.fraction(sp.together(e)); nm = sp.expand(nm.subs(yy, sp.sqrt(Y)))
    if nm.has(sp.sqrt(Y)): lin = False
    pe = sp.Poly(nm, RR, KK)
    if not set(pe.monoms()) <= {(1, 1), (1, 0), (0, 1), (0, 0)}: lin = False
    Mrows.append([pe.coeff_monomial(RR*KK), pe.coeff_monomial(RR), pe.coeff_monomial(KK)]); brow.append(-pe.coeff_monomial(1))
check('rows are polynomials in Y=y^2 and linear in (RK,R,K)', lin)
M3 = sp.Matrix(Mrows); b3 = sp.Matrix(brow)
x_, z_, r_, s_ = X, zr, rr, ss
N = z_**2+Y*v**2; Aa = x_**2+Y; t_ = z_-v*x_; m = r_+v-1
Fp = (r_-s_)*v-r_-t_; TN = (r_*x_-(r_-s_)*z_)**2+Y*(r_-(r_-s_)*v)**2
Dn = t_+s_*v; Hh = (r_-1)*t_-v*(r_*(r_-s_)+s_)
det = sp.expand(M3.det(method='berkowitz'))
check('(4) det M3 = -s N^2 m Fp TN (row normalization as above)', sp.expand(det+s_*N**2*m*Fp*TN) == 0)
def cram(k):
    Mk = M3.copy(); Mk[:, k] = b3; return sp.expand(Mk.det(method='berkowitz'))
NRK, NR, NK = cram(0), cram(1), cram(2)
R0 = sp.cancel(NR/(-s_*N*Fp*TN))
check('N_R = -s N Fp TN R0 with R0 a polynomial', sp.fraction(sp.together(R0))[1] == 1)
check('(5) m N + v R0 = t H', sp.expand(m*N+v*R0-t_*Hh) == 0)
check('(6) H+(r-1)Fp = -rm and (r-1)(Dn+r^2)-H = r^2 m',
      sp.expand(Hh+(r_-1)*Fp+r_*m) == 0 and sp.expand((r_-1)*(Dn+r_**2)-Hh-r_**2*m) == 0)
Lc_RK = -t_*v*N; Lc_K = -t_*N+Dn*N; Lc_R = Aa*t_*v*N; Lc_0 = Aa*t_*N-Aa*Dn**2
check('(9) m*(L numerator at the Cramer solution)*det = -A t Dn (m-H) det',
      sp.expand(m*(Lc_RK*NRK+Lc_R*NR+Lc_K*NK+Lc_0*det)+Aa*t_*Dn*(m-Hh)*det) == 0)
Ly3 = -t_*Dn*(m-Hh)/(m*N**2); By2 = -t_**2*Hh/(m*N**2)
check('(10) F6/y^6 = L^2/y^6 - 4 (d/y)^2 B/y^2 = t^2 Dn^2 (m+H)^2/(m^2 N^4)',
      sp.cancel(Ly3**2-4*(Dn/N)**2*By2-t_**2*Dn**2*(m+Hh)**2/(m**2*N**4)) == 0)
print('FAILED:', fails if fails else 'none'); sys.exit(1 if fails else 0)
