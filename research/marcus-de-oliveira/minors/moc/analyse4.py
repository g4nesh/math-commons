import sys, json, glob, numpy as np, itertools
sys.path.insert(0, 'moc')
from bad4 import viol
R = []
for f in glob.glob('moc/runs/bad4_s1*.json'): R += json.load(open(f))
J = np.array([r['J'] for r in R]); V = np.array([r['viol'] for r in R])
print("total bad U:", len(R), " max J %.3e" % J.max())
print("corr(J, violation depth) = %.3f" % np.corrcoef(J, V)[0, 1])
for q in [0.0, 0.01, 0.02, 0.03]:
    m = V > q
    if m.any(): print("  U with violation > %.2f: %d, their max J %.3e" % (q, m.sum(), J[m].max()))
P = np.array(list(itertools.permutations(range(4))))
for r in sorted(R, key=lambda r: -r['J'])[:4]:
    a = np.array(r['a'][0]) + 1j * np.array(r['a'][1]); b = np.array(r['b'][0]) + 1j * np.array(r['b'][1])
    z = np.prod(a[None, :] + b[P], axis=1)
    # collinearity measures: spread of a and b around best-fit lines
    def lin(v):
        X = np.column_stack([v.real - v.real.mean(), v.imag - v.imag.mean()]); s = np.linalg.svd(X, compute_uv=False); return s[1] / s[0]
    zx = np.column_stack([z.real, z.imag]); sz = np.linalg.svd(zx - zx.mean(0), compute_uv=False)
    print("J=%.3e viol=%.4f | a-collinearity ratio %.3f, b-collinearity %.3f, hull-points flatness %.3f, min|a_i - a_j| %.3f, min|b_i-b_j| %.3f" %
          (r['J'], r['viol'], lin(a), lin(b), sz[1] / sz[0], min(abs(a[i]-a[j]) for i in range(4) for j in range(i)), min(abs(b[i]-b[j]) for i in range(4) for j in range(i))))
