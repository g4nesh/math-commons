"""Whole degree-4 pipeline at one rational t (given as p/q):  central solution, second reduction, exact solve, and the
independent check with two fresh seeds.  Logs go to logs/*_t<p>-<q>.txt ; the certificate to certificate4_t<p>-<q>.json."""
import os, subprocess, sys
from fractions import Fraction as Fr
tq = Fr(sys.argv[1]); t = float(tq); tag = "%d-%d" % (tq.numerator, tq.denominator) if tq.denominator != 1 else "%d" % tq.numerator
def run(args, log):
    out = subprocess.run([sys.executable] + args, capture_output=True, text=True)
    text = "\n".join(l for l in (out.stdout + out.stderr).splitlines() if 'Warning' not in l and 'warn' not in l)
    open('logs/%s_t%s.txt' % (log, tag), 'w').write(text + "\n"); print(text.splitlines()[-1] if text.strip() else "(no output)", flush=True)
    if out.returncode != 0: sys.exit("step %s failed" % log)
os.environ['OMP_NUM_THREADS'] = '1'
run(['center4.py', '-', repr(t)], 'center4'); run(['second4.py', str(tq)], 'second4'); run(['solve4.py', str(tq)], 'solve4')
cert = 'certificate4_t%s.json' % tag; os.replace('certificate4_t%g.json' % t, cert)
run(['verify4.py', cert, '777'], 'verify4'); run(['verify4.py', cert, '20261125'], 'verify4_seed2')
