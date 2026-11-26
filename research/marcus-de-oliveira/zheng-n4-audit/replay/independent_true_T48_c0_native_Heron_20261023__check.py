"""Replay for the note 'A Heron identity and orientation condition on the T48 c = 0 face'
(independent_true_T48_c0_native_Heron_20261023, supplement).

Own code only; nothing from the supplement is read or executed.  Argument 1: path of moc4lib.
Exact: the short coordinate functions, the c = 0 consequences (1), identities (6), the
derivation of (8) from (1) and of (10) from (7) and (9).  The three-row system (3)-(5) is
representation-specific and is tested NUMERICALLY through its consequences (7) and (9) on
points of the seven-contact stratum constructed from the explicit spectral inverse.
Exit status nonzero iff an EXACT check fails.
"""
import os, sys, time, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "lib"))
import moc4lib as L
import sympy as sp
import numpy as np
from scipy.optimize import fsolve

FAIL = []
def check(name, ok, level="EXACT"):
    print(f"[{level}] {'PASS' if ok else 'FAIL'}  {name}")
    if not ok and level == "EXACT":
        FAIL.append(name)
t0 = time.time()
T7 = frozenset([0, 3, 5, 6, 9, 12, 17])
check("contacts {0,3,5,6,9,12,17} = Table 5 set #23 (paper T_49 plus 17)", T7 == L.AUX7[23] and (T7 - {17}) == L.CONTACT_SETS[49])
G8 = [7, 5, 29, 39, 2, 40, 41, 0]
names = ['q', 'p', 'j', 'c', 'e', 'd', 'b', 'a']
check("canonical rays (7,5,29,39,2,40,41,0); the only one nonzero at 17 is f39 (c), so contact 17 <=> c = 0",
      [names[k] for k, jj in enumerate(G8) if L.F_VALUES[jj][17] != 0] == ['c'] and L.F_VALUES[39][17] == 1)

# ---------------- explicit inverse chart (as in the companion harmonic note, eq. (9)) ----------
r, s, u, h, J, W, R, Z, K = sp.symbols('r s u h J W R Z K')
dl = r - s
M = dl*W + s; ZD = J*(W - 1) - r; Kc = J*W - s
BD = W - (J + u)*(W - 1)
YN = W*(Kc - J*M); YD = Kc*(W - u*(W - 1)) - J*M
Ah = [(YN, YD), (-M, ZD), (W, 1), (s, J)]
Bh = [(W, BD), (0, 1), (1, 1), (1, 0)]
br = lambda p_, q_: p_[0]*q_[1] - p_[1]*q_[0]
Pi = [sp.expand(sp.prod([br(Ah[i], Bh[sg[i]]) for i in range(4)])) for sg in L.PERMS]
E = sp.expand(Pi[0] - Pi[6])
Hs = [sp.cancel((Pi[k] - Pi[0])/E) for k in range(24)]
check("contacts: H0=0, H3=s-r, H5=-r, H6=-1, H12=-u identically",
      [sp.cancel(Hs[k] - v) for k, v in ((0, 0), (3, s - r), (5, -r), (6, -1), (12, -u))] == [0]*5)
num9, den9 = sp.fraction(sp.cancel(Hs[9]))
C = sp.expand(num9 - h*den9)
# amplitudes (dual basis of the 8 canonical generators, completed) as in the companion check
from fractions import Fraction
basis = [list(L.F_VECS[jj]) for jj in G8]
for i in range(14):
    ev = [Fraction(0)]*14; ev[i] = Fraction(1)
    if L.rank(basis + [ev]) > len(basis):
        basis.append(ev)
Mf = L.matmul(L.inv(L.transpose(basis)), L.inv([L.P[rw] for rw in L.P_PIVOT_ROWS]))
Qr = lambda f: sp.Rational(f.numerator, f.denominator)
amp = {names[k]: sp.cancel(sum(Qr(Mf[k][i])*Hs[L.P_PIVOT_ROWS[i]] for i in range(14))) for k in range(8)}
print("  time", round(time.time() - t0, 1))
# doc-2 short coordinate functions, with W = 1 - Z/J, u = 1 - K/s, h = R - r + s
H9 = Hs[9]
doc2 = {'j': J, 'c+d': (J - s)/Z, 'a+d': -(u - 1)*(J - s)/J, 'q+d': -(R*Z + J + r - s)/r,
        'b': ((r - 1)*J + s)/Z - (r - s), 'p': R*(Z + r)/(J + r - s) - (r + h - 1)}
for nm, f in doc2.items():
    fW = f.subs(Z, J*(1 - W)).subs(R, dl + h).subs(h, H9)
    D = sp.cancel(sum(amp[x] for x in nm.split('+')) - fW)
    c9 = sp.cancel(sp.diff(D, J)/sp.diff(H9, J))
    rest = sp.cancel(D - c9*H9)
    check(f"short coordinate function {nm}_amp equals the printed one up to a real shift ({c9}*h + {rest})",
          c9.is_number and not (rest.has(J) or rest.has(W)))
# (1) and (8), (10): algebra in real variables
x, y, z, v = sp.symbols('x y z v', real=True)
rr, ss, RR, KK = sp.symbols('r s R K', real=True)
Jc = x + sp.I*y; Zc = z + sp.I*y*v
uu = 1 - KK/ss; hh = RR - rr + ss
imag = lambda e: sp.simplify(sp.im(sp.expand_complex(e)))
sub2 = {J: Jc, Z: Zc, r: rr, s: ss, R: RR, u: uu, h: hh}
cd_ = imag(doc2['c+d'].subs(sub2)); ad_ = imag(doc2['a+d'].subs(sub2)); qd_ = imag(doc2['q+d'].subs(sub2))
b_ = imag(doc2['b'].subs(sub2)); p_ = imag(doc2['p'].subs(sub2))
t_ = z - v*x; m_ = rr + v - 1; N_ = z**2 + y**2*v**2; A_ = x**2 + y**2
Dn = t_ + ss*v; Fp = (rr - ss)*v - rr - t_; Hh = (rr - 1)*t_ - v*(rr*(rr - ss) + ss)
d_ = cd_          # c = 0
check("(1) d/y = Dn/N on c = 0", sp.simplify(d_/y - Dn/N_) == 0)
check("(1) (b+d)/y = r t/N", sp.simplify((b_ + d_)/y - rr*t_/N_) == 0)
check("(1) (q+d)/y = -(Rv+1)/r", sp.simplify(qd_/y + (RR*v + 1)/rr) == 0)
check("(1) (a+d)/y = K/A", sp.simplify(ad_/y - KK/A_) == 0)
check("(1) p/y = R Fp/|J + r - s|^2", sp.simplify(p_/y - RR*Fp/((x + rr - ss)**2 + y**2)) == 0)
check("(6) H + (r-1)Fp = -r m", sp.expand(Hh + (rr - 1)*Fp + rr*m_) == 0)
check("(6) (r-1)(Dn + r^2) - H = r^2 m", sp.expand((rr - 1)*(Dn + rr**2) - Hh - rr**2*m_) == 0)
# (8) from (1): L = a(B + dj) + (d - j)B, j = y
a_ = ad_ - d_; B_ = (b_ + d_)*qd_
Lexpr = a_*(B_ + d_*y) + (d_ - y)*B_
Bn = -t_*(RR*v + 1)
Lnum = (KK - A_)*Bn*N_ + KK*Dn*N_ - A_*Dn**2
check("(8) L/y^3 = [(K-A) Bn N + K Dn N - A Dn^2]/(A N^2) with Bn = -t(Rv+1)", sp.simplify(Lexpr/y**3 - Lnum/(A_*N_**2)) == 0)
# (10) from (7) and (9): substitute Rv + 1 = tH/(mN) and L/y^3 = -t Dn (m-H)/(m N^2)
Bq = sp.simplify((B_/y**2).subs(RR, (t_*Hh/(m_*N_) - 1)/v))
check("(10) B/y^2 = -t^2 H/(m N^2) given (7)", sp.simplify(Bq + t_**2*Hh/(m_*N_**2)) == 0)
L3 = -t_*Dn*(m_ - Hh)/(m_*N_**2)
F6 = L3**2 - 4*(Dn/N_)**2*(-t_**2*Hh/(m_*N_**2))
check("(10) F6/y^6 = t^2 Dn^2 (m+H)^2/(m^2 N^4) given (7),(9)", sp.simplify(F6 - t_**2*Dn**2*(m_ + Hh)**2/(m_**2*N_**4)) == 0)
print("  time", round(time.time() - t0, 1))

# ---------------- NUMERIC: points of the seven-contact stratum --------------------------
Pi_f = [sp.lambdify((J, W, r, s, u), Pi[k], 'numpy') for k in range(24)]
GV = np.array([L.F_VALUES[jj] for jj in G8], float).T     # 24 x 8
def Hvec(Jv, Wv, rv, sv, uv):
    P_ = np.array([complex(f(Jv, Wv, rv, sv, uv)) for f in Pi_f])
    Ev = P_[0] - P_[6]
    return (P_ - P_[0])/Ev
rng = np.random.default_rng(31)
pts = 0; pos = 0; bad7 = bad9 = bad1 = 0; worst_F6 = np.inf; sign_bad = 0; r_pos_pos = 0
maxres = 0.0
for trial in range(4000):
    rv, sv, uv = rng.uniform(-3, 3, 3)
    Jv = complex(rng.uniform(-3, 3), rng.uniform(0.05, 3))
    def eqs(X):
        Wv = complex(X[0], X[1])
        Hv = Hvec(Jv, Wv, rv, sv, uv)
        return [Hv[9].imag, Hv[17].imag]
    X0 = rng.uniform(-3, 3, 2)
    try:
        sol, info, ier, msg = fsolve(eqs, X0, full_output=True, xtol=1e-13)
    except Exception:
        continue
    if ier != 1: continue
    Wv = complex(sol[0], sol[1])
    Hv = Hvec(Jv, Wv, rv, sv, uv)
    if max(abs(Hv[9].imag), abs(Hv[17].imag)) > 1e-10 or abs(Wv.imag) < 1e-6: continue
    hv = Hv[9].real
    g = Hv.imag            # eta = -1: coefficients = imaginary parts of amplitudes; j = Im J > 0
    coef, *_ = np.linalg.lstsq(GV, g, rcond=None)
    resid = np.max(np.abs(GV @ coef - g)); maxres = max(maxres, resid)
    q_, p_, j_, c_, e_, d_, b_, a_ = coef
    pts += 1
    xv, yv = Jv.real, Jv.imag
    Zv = Jv*(1 - Wv); zv = Zv.real; vv = Zv.imag/yv
    Kv = sv*(1 - uv); Rv = rv - sv + hv
    tv = zv - vv*xv; mv = rv + vv - 1; Nv = abs(Zv)**2; Av = abs(Jv)**2
    Dnv = tv + sv*vv; Hhv = (rv - 1)*tv - vv*(rv*(rv - sv) + sv)
    sc = 1 + abs(mv*Nv*(Rv*vv + 1)) + abs(tv*Hhv)
    if abs(mv*Nv*(Rv*vv + 1) - tv*Hhv) > 1e-6*sc: bad7 += 1
    if abs(d_/yv - Dnv/Nv) > 1e-6*(1 + abs(Dnv/Nv)) or abs(j_ - yv) > 1e-8 or abs(c_) > 1e-8: bad1 += 1
    Bv = (b_ + d_)*(q_ + d_)
    Lv = a_*(Bv + d_*j_) + (d_ - j_)*Bv
    L9 = -tv*Dnv*(mv - Hhv)/(mv*Nv**2)*yv**3
    if abs(Lv - L9) > 1e-6*(1 + abs(Lv)): bad9 += 1
    if min(p_, j_, d_) > 1e-9 and min(q_, b_) > -1e-12:
        pos += 1
        F6v = Lv**2 - 4*d_**2*j_**2*Bv
        worst_F6 = min(worst_F6, F6v/(1 + Lv**2))
        if np.sign(Lv) != -np.sign(rv): sign_bad += 1
        if rv > 0 and j_ > d_: r_pos_pos += 1
print(f"  numeric (seed 31): {pts} points of the seven-contact stratum (max decomposition residual {maxres:.1e}); "
      f"(1) failures {bad1}, (7) failures {bad7}, (9) failures {bad9}")
print(f"  numeric: {pos} of them satisfy p,j,d > 0, q,b >= 0; worst normalized F6 = {worst_F6:.3e}; "
      f"sign L != -sign r in {sign_bad}; points in the open branch r > 0, j > d: {r_pos_pos}")
check("NUMERIC: identities (1), (7), (9) hold at all sampled stratum points", pts > 0 and bad1 == bad7 == bad9 == 0, "NUMERIC")
check("NUMERIC: F6 >= 0 and sign L = -sign r at all sampled supported points", sign_bad == 0 and (pos == 0 or worst_F6 > -1e-9), "NUMERIC")
print("  total time", round(time.time() - t0, 1))
print("FAILED:", FAIL if FAIL else "none")
sys.exit(1 if FAIL else 0)
