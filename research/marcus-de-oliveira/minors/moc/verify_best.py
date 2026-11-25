"""Re-evaluate the best configurations in 50-digit arithmetic and diagnose them."""
import sys, json, itertools, glob, numpy as np, mpmath as mp
mp.mp.dps = 50
def analyse(f):
    d = json.load(open(f)); n = d['n']
    U = np.array(d['U_re']) + 1j * np.array(d['U_im'])
    a = np.array(d['a'][0]) + 1j * np.array(d['a'][1]); b = np.array(d['b'][0]) + 1j * np.array(d['b'][1])
    th = d['th']
    # re-unitarise U in high precision (polar decomposition)
    Um = mp.matrix([[mp.mpc(U[i, j].real, U[i, j].imag) for j in range(n)] for i in range(n)])
    for _ in range(3):   # Newton polar iteration U <- (U + U^{-H})/2
        Um = (Um + (Um**-1).transpose_conj()) / 2
    am = [mp.mpc(x.real, x.imag) for x in a]; bm = [mp.mpc(x.real, x.imag) for x in b]
    M = mp.diag(am) + Um * mp.diag(bm) * Um.transpose_conj()
    dd = mp.det(M)
    e = mp.exp(-1j * th)
    perms = list(itertools.permutations(range(n)))
    z = [mp.fprod([am[i] + bm[s[i]] for i in range(n)]) for s in perms]
    v = [mp.re(e * x) for x in z]
    vmax = max(v); W = vmax - min(v)
    J = (mp.re(e * dd) - vmax) / W
    # diagnostics: closeness of |U|^2 to a permutation matrix; top two vertices in direction
    M1 = np.abs(U)**2
    perm_dist = min(np.abs(M1 - np.eye(n)[list(s)]).sum() for s in perms) if n <= 6 else np.nan
    order = np.argsort([-float(x) for x in v])
    s1, s2 = perms[order[0]], perms[order[1]]
    diff = sum(1 for i in range(n) if s1[i] != s2[i])
    gap12 = float((v[order[0]] - v[order[1]]) / W)
    return n, float(J), perm_dist, diff, gap12
for f in sorted(glob.glob('moc/runs/b_n*_s*.json')):
    if '_s99' in f: continue
    n, J, pd, diff, g = analyse(f)
    print("%s n=%d J(50dig)=%+.3e  L1 dist |U|^2 to perm=%.2e  top-2 vertices differ in %d positions, their gap/W=%.2e" % (f, n, J, pd, diff, g))
