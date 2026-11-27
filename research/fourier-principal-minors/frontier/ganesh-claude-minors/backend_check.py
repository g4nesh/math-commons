#!/usr/bin/env python3
"""Compare the two builds of pmscan.c (tables and polynomial arithmetic) on the
same charts: the lists of vanishing principal minors must be identical.

Usage: python3 backend_check.py BIN_DIR M:p [M:p ...]
"""
import os
import subprocess
import sys

sys.path.insert(0, ".")
import frontier as fr


def zeros(exe, M, p, f, h):
    cmd = [exe, str(M), str(p), str(f)] + [str(c) for c in h[:-1]] + ["list"]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    return sorted(int(l.split()[1], 16) for l in out.splitlines() if l.startswith("Z "))


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(2)
    bin_dir = sys.argv[1]
    for arg in sys.argv[2:]:
        M, p = map(int, arg.split(":"))
        f, facs = fr.cyclotomic_factors(M, p)
        for h in facs:
            a = zeros(os.path.join(bin_dir, "pmscan"), M, p, f, h)
            b = zeros(os.path.join(bin_dir, "pmscan_poly"), M, p, f, h)
            assert a == b, (M, p, h)
            print("(M,p)=(%d,%d) h=%s: %d vanishing minors in both builds" % (M, p, h, len(a)))
    print("both builds agree")
