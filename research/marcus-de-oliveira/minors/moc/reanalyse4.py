"""Re-analysis of the archived n = 4 targeted search (moc/runs/bad4_s11..19.json), prompted by the independent review
in ../reviews/mdo-minors.md (finding 1).  Each file holds 60 records: the first 30 unitaries were 'deepened' by gradient
ascent on the violation, the last 30 are raw Haar violators.  For every record U is FIXED and only (a, b, theta) vary.

Reported here:
 (a) the symmetry-invariant depth  min over P_alpha U P_beta and transposes of the exotic functional g, per record;
 (b) correlations of J with the stored (non-invariant) violation and with the invariant depth, overall and within groups;
 (c) the three best stored J, re-optimised at fixed U by restarted Nelder-Mead (seed 0)."""
import sys, json, itertools, numpy as np
from scipy.optimize import minimize
sys.path.insert(0, 'moc')
from bsearch import Problem, normalise
PERMS = list(itertools.permutations(range(4)))
def g(U):
    D = lambda S, T: abs(U[S[0], T[0]] * U[S[1], T[1]] - U[S[0], T[1]] * U[S[1], T[0]])**2
    return abs(U[0, 0])**2 + abs(U[2, 2])**2 + D((0, 3), (2, 3)) - D((0, 1), (0, 3))
def depth(U):
    return min(g(V) for al in PERMS for be in PERMS for V in (U[np.ix_(al, be)], U[np.ix_(al, be)].T))
recs = []
for s in range(11, 20):
    for k, r in enumerate(json.load(open('moc/runs/bad4_s%d.json' % s))):
        U = np.array(r['U_re']) + 1j * np.array(r['U_im']); r['U'] = U; r['deep'] = k < 30; r['depth'] = depth(U); r['file'] = s; recs.append(r)
J = np.array([r['J'] for r in recs]); viol = np.array([r['viol'] for r in recs]); dep = np.array([r['depth'] for r in recs]); deep = np.array([r['deep'] for r in recs])
print("records: %d (deepened %d, raw %d)" % (len(recs), deep.sum(), (~deep).sum()))
print("(a) invariant depth of the deepened unitaries: min %.6f, max %.6f  (global minimum of g over U(4): -0.067442)" % (dep[deep].min(), dep[deep].max()))
prof = np.array([np.sort((np.abs(r['U'])**2).ravel()) for r in recs if r['deep']])
print("    sorted |u_ij|^2 profiles of the deepened unitaries: largest standard deviation over the 16 positions %.1e" % prof.std(0).max())
print("    invariant depth of the raw unitaries: from %.6f to %.6f" % (dep[~deep].min(), dep[~deep].max()))
print("    stored (non-invariant) violation of the deepened unitaries: from %.4f to %.4f" % (viol[deep].min(), viol[deep].max()))
c = lambda x, y: float(np.corrcoef(x, y)[0, 1])
print("(b) corr(J, stored violation): all %.2f ; raw only %.2f ; deepened only %.2f" % (c(J, viol), c(J[~deep], viol[~deep]), c(J[deep], viol[deep])))
print("    corr(J, -invariant depth): all %.2f ; raw only %.2f" % (c(J, -dep), c(J[~deep], -dep[~deep])))
print("    J over the deepened copies of the extremal point: from %.4f to %.4f ; J over raw points: from %.2e to %.2e" % (J[deep].min(), J[deep].max(), J[~deep].min(), J[~deep].max()))
pr = Problem(4); rng = np.random.default_rng(0)
def Jval(x, U):
    a = (x[0:4] + 1j * x[4:8])[None]; b = (x[8:12] + 1j * x[12:16])[None]; a, b = normalise(a, b)
    return float(pr.evaluate(U[None], a, b, np.array([x[16]]), None, grad=False)[0][0])
print("(c) best stored J re-optimised at the same fixed U (Nelder-Mead, 12 restarts each):")
for i in np.argsort(-J)[:3]:
    r = recs[i]; x0 = np.concatenate([r['a'][0], r['a'][1], r['b'][0], r['b'][1], [r['th']]]); best = Jval(x0, r['U']); xb = x0
    for rs in range(12):
        x = xb + (0 if rs == 0 else 1e-3 * rng.standard_normal(17))
        res = minimize(lambda z: -Jval(z, r['U']), x, method='Nelder-Mead', options={'maxiter': 20000, 'maxfev': 20000, 'xatol': 1e-12, 'fatol': 1e-15})
        if -res.fun > best: best = -res.fun; xb = res.x
    a = xb[0:4] + 1j * xb[4:8]; b = xb[8:12] + 1j * xb[12:16]; gap = lambda v: min(abs(v[p] - v[q]) for p in range(4) for q in range(p + 1, 4)) / np.abs(v - v.mean()).max()
    print("    file s%d, %s point, invariant depth %.4f: stored J = %.3e -> re-optimised J = %.3e (still negative: %s) ; smallest relative gap among a: %.1e, among b: %.1e"
          % (r['file'], 'deepened' if r['deep'] else 'raw', r['depth'], r['J'], best, best < 0, gap(a), gap(b)))
