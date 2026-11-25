import sys; sys.path.insert(0, 'moc')
import margin4 as m
for t in [0.07, 0.08, 0.09]:
    mu = m.mu_of(m.foot + t * m.nF)
    print("  t = %.3f (%.2fx actual depth): max J = %+.4e" % (t, t / m.depth, m.maxJ(mu, 200)), flush=True)
