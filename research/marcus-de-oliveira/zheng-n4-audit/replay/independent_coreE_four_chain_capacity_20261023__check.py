"""Exact replay for independent_coreE_four_chain_capacity_20261023 (a sufficient coefficient bound
on the seven-contact set E = {0,1,2,9,10,16,17}).

All checks are own code (sympy, moc4lib).  Nothing of the supplement is executed or read.
Run: python -I -B <this file> [--numeric].  Exits 1 if an EXACT check fails.
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L
import sympy as sp

FAIL = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + ("  " + str(info) if info else ""), flush=True)
    if not ok:
        FAIL.append(name)

E = [0, 1, 2, 9, 10, 16, 17]
check("generators(E) = {f0,f2,f3,f12,f13,f15,f26}; the displayed source uses exactly these",
      L.generators(E) == [0, 2, 3, 12, 13, 15, 26])

q, v, c, d, t, dl = sp.symbols('q v c d t delta', positive=True)
lam = d - dl
A = dl + dl**2 / v; B = lam**2 / c; K = A + B; W = q + dl - K
Pd = c * v * t * (q + dl) - (t + q + dl) * (c * v * dl + c * dl**2 + v * (d - dl)**2)
check("W = q - delta^2/v - B = q + delta - K", sp.cancel(W - (q - dl**2 / v - B)) == 0)
check("(1) P(delta) = c v (t W - (q+delta) K)  [so P >= 0 iff t W >= (q+delta) K]",
      sp.expand(sp.cancel(Pd - c * v * (t * W - (q + dl) * K))) == 0)
X = W * A / K
qF = dl**2 / v + X; qC = q - qF; tF = A + A**2 / X; tC = B + B**2 / (W - X)
check("q_F v (t_F - delta) = delta^2 (q_F + v + t_F + delta)", sp.cancel(qF * v * (tF - dl) - dl**2 * (qF + v + tF + dl)) == 0)
check("lambda^2/q_C + lambda^2/t_C = c", sp.cancel(lam**2 / qC + lam**2 / tC - c) == 0)
check("t_F + t_C = K + K^2/W = (q+delta)K/(q+delta-K)", sp.cancel(tF + tC - K - K**2 / W) == 0 and sp.cancel(tF + tC - (q + dl) * K / W) == 0)
check("q_F v - delta^2 = X v and q_C - B = W - X = W B/K", sp.cancel(qF * v - dl**2 - X * v) == 0 and sp.cancel(qC - B - W * B / K) == 0)
check("delta = 0: P(0) = v{c q t - d^2 (q + t)}", sp.expand(Pd.subs(dl, 0) - v * (c * q * t - d**2 * (q + t))) == 0)
check("delta = d: P(d) = c{q v (t - d) - d^2 (q + v + t + d)}", sp.expand(Pd.subs(dl, d) - c * (q * v * (t - d) - d**2 * (q + v + t + d))) == 0)
M = sp.Matrix([[q, dl, dl, d - dl], [dl, t - dl, -dl, -(d - dl)], [dl, -dl, v, 0], [d - dl, -(d - dl), 0, c]])
Sch = M[:2, :2] - M[:2, 2:] * M[2:, 2:].inv() * M[2:, :2]
check("capacity matrix: Schur complement in the last two diagonals is [[W, K], [K, t - K]]",
      all(sp.cancel(x) == 0 for x in (Sch - sp.Matrix([[W, K], [K, t - K]]))))
check("capacity matrix: det = P(delta)", sp.expand(M.det() - Pd) == 0)
sub = {dl: d * v / (c + v)}
check("delta = dv/(c+v): W = q - d^2/(c+v)", sp.cancel(W.subs(sub) - (q - d**2 / (c + v))) == 0)
R = d * (v + d) * (q * (c + v) + d * v) / ((c + v) * (q * (c + v) - d**2))
check("delta = dv/(c+v): (q+delta)K/W equals the displayed simple threshold", sp.cancel(((q + dl) * K / W).subs(sub) - R) == 0)
check("delta = dv/(c+v) lies in (0, d)", sp.simplify(d - d * v / (c + v) - d * c / (c + v)) == 0)
# minimization claim
Xs, As, Bs, Ws = sp.symbols('X A B W', positive=True)
check("A^2/X + B^2/(W-X) - (A+B)^2/W = (A W - (A+B) X)^2 / (X (W-X) W) >= 0, equality at X = W A/(A+B)",
      sp.cancel(As**2 / Xs + Bs**2 / (Ws - Xs) - (As + Bs)**2 / Ws - (As * Ws - (As + Bs) * Xs)**2 / (Xs * (Ws - Xs) * Ws)) == 0)
# example
dd, qq, vv, cc = Fr(1), Fr(100), Fr(1, 50), Fr(1, 50)
dF, dC = Fr(1, 2), Fr(1, 2); qF_, tF_, qC_, tC_ = Fr(40), Fr(20), Fr(60), Fr(300, 19)
check("example: four inequality strict: qF v (tF - 1/2) > (1/4)(qF + v + tF + 1/2)",
      qF_ * vv * (tF_ - dF) > dF**2 * (qF_ + vv + tF_ + dF), (qF_ * vv * (tF_ - dF), dF**2 * (qF_ + vv + tF_ + dF)))
check("example: (1/2)^2/60 + (1/2)^2/(300/19) = 1/50", dC**2 / qC_ + dC**2 / tC_ == cc)
check("example: allocations exhaust q = 100 and total t = 680/19", qF_ + qC_ == qq and tF_ + tC_ == Fr(680, 19))
TFour = dd * (qq + dd) * (vv + dd) / (qq * vv - dd**2)
TChain = dd**2 / (cc - dd**2 / qq)
Rex = dd * (vv + dd) * (qq * (cc + vv) + dd * vv) / ((cc + vv) * (qq * (cc + vv) - dd**2))
print(f"INFO example: four alone needs t >= {TFour} = {float(TFour):.3f}; chain alone needs t >= {TChain}; "
      f"combined simple allocation 680/19 = {680/19:.3f}; rational section delta=dv/(c+v) gives {Rex} = {float(Rex):.3f}")
check("example: 680/19 is below both separate thresholds (claim 'beyond separately tested capacities')",
      Fr(680, 19) < TFour and Fr(680, 19) < TChain)

if "--numeric" in sys.argv:
    import numpy as np
    from scipy.optimize import minimize
    Pn = L.P_numpy(); FV = {j: np.array([float(x) for x in L.F_VECS[j]]) for j in (0, 2, 3, 12, 13, 15, 26)}
    rng = np.random.default_rng(11)
    def minphi(co, restarts=8):
        vals = Pn @ sum(w * FV[j] for j, w in co.items())
        f = lambda th: float(vals @ L.leibniz_weights_np(L.so4_from_params(th)))
        return min(minimize(f, rng.uniform(-np.pi, np.pi, 6), method='BFGS').fun for _ in range(restarts)) / np.abs(vals).max()
    worst = (np.inf, None); N = 60
    for k in range(N):
        qx, vx, cx = np.exp(rng.uniform(-2.5, 2.5, 3)); dx = rng.uniform(0.02, 0.98)
        A_ = dx + dx**2 / vx; B_ = (1 - dx)**2 / cx; K_ = A_ + B_; W_ = qx + dx - K_
        if W_ <= 0: continue
        tx = (qx + dx) * K_ / W_          # P(delta) = 0: boundary of the hypothesis
        m = minphi({0: qx, 13: vx, 15: cx, 3: tx, 26: 1.0})
        if m < worst[0]: worst = (m, (qx, vx, cx, dx, tx))
    print(f"NUMERIC P(delta)=0 boundary, p=b=0, d=1: {N} draws (seed 11, 8 BFGS restarts each): worst min Phi/max|Pc| = {worst[0]:.2e} at (q,v,c,delta,t) = {np.round(worst[1], 4)}")
    m = minphi({0: 100.0, 13: 0.02, 15: 0.02, 3: 680 / 19, 26: 1.0}, 30)
    print(f"NUMERIC example q=100, v=c=1/50, t=680/19, d=1: min Phi/max|Pc| = {m:.2e} (30 restarts)")

print("SUMMARY", "ALL EXACT CHECKS PASS" if not FAIL else f"{len(FAIL)} FAILED: {FAIL}")
sys.exit(1 if FAIL else 0)
