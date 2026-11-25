"""Numerical search, exact solve and independent check (two seeds) for one interval.
Usage: python3 run_interval.py <deg> <a> <b> <d>      Logs: logs/build_interval<tag>.txt, logs/verify_interval<tag>*.txt"""
import subprocess, sys
deg, a, b, d = sys.argv[1:5]; tag = '%s_%s_%s' % (deg, a.replace('/', '-'), b.replace('/', '-'))
def run(args):
    out = subprocess.run([sys.executable] + args, capture_output=True, text=True)
    text = "\n".join(l for l in (out.stdout + out.stderr).splitlines() if 'Warning' not in l and 'warn' not in l)
    print(text.splitlines()[-1][:230] if text.strip() else "(no output)", flush=True)
    if out.returncode != 0 or 'A CHECK FAILED' in text: open('logs/failed_interval%s.txt' % tag, 'w').write(text + "\n"); sys.exit("failed: " + " ".join(args))
    return text
build = run(['interval_sdp.py', deg, a, b, d]) + "\n" + run(['interval_exact.py', deg, a, b, d]); open('logs/build_interval%s.txt' % tag, 'w').write(build + "\n")
open('logs/verify_interval%s.txt' % tag, 'w').write(run(['verify_interval.py', 'interval%s.json' % tag, '777']) + "\n")
open('logs/verify_interval%s_seed2.txt' % tag, 'w').write(run(['verify_interval.py', 'interval%s.json' % tag, '20261125']) + "\n")
