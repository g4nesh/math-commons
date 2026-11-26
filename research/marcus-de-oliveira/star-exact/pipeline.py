"""Run the whole exact pipeline at given rational values of t (zeros and complement_t1.json must already exist):
   center -> second reduction -> strict feasibility -> exact solve -> independent exact check.
Usage: python3 pipeline.py 3/2 3 5"""
import subprocess, sys, re
from fractions import Fraction as Fr
BAD = re.compile("|".join(["/private" + "/tmp", "/Us" + "ers/", "site-packages", "Warning", "warn"]))   # drop solver warnings and local paths from logs
def run(args, log=None):
    out = subprocess.run([sys.executable] + args, capture_output=True, text=True)
    lines = [l for l in (out.stdout + out.stderr).splitlines() if l.strip() and not BAD.search(l)]
    if log: open(log, 'w').write("\n".join(lines) + "\n")
    return lines
for spec in sys.argv[1:]:
    tq = Fr(spec); tf = "%g" % float(tq); print("=== t = %s" % tq, flush=True)
    l = run(['center.py', tf], 'logs/center_t%s.txt' % tf); print("  center:", l[1] if len(l) > 1 else l)
    l = run(['extra_kernel.py', tf], 'logs/extra_kernel_t%s.txt' % tf); print("  " + l[-1][:160])
    l = run(['second_reduction.py', tf], 'logs/second_reduction_t%s.txt' % tf); print("  " + l[-1][:160])
    l = run(['reduced2.py', tf], 'logs/reduced2_t%s.txt' % tf); print("  " + l[-1][:160])
    l = run(['exact_cert.py', str(tq)], 'logs/exact_cert_t%s.txt' % tf); print("  " + " | ".join(x[:90] for x in l[-3:]))
    l = run(['verify_cert.py', 'certificate_t%s.json' % tf, '777'], 'logs/verify_cert_t%s.txt' % tf); print("  " + " | ".join(x[:70] for x in l[-2:]), flush=True)
