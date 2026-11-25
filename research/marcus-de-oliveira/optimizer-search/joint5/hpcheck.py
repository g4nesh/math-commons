"""Re-evaluate the best configurations of a run file at high precision.
Usage: python hpcheck.py RUN.json [K] [DPS]"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
from moc_hp import hp_margin  # noqa: E402
import mpmath as mp  # noqa: E402

def load(rec):
    a = np.array([complex(*v) for v in rec["a"]])
    b = np.array([complex(*v) for v in rec["b"]])
    U = np.array(rec["U_re"]) + 1j * np.array(rec["U_im"])
    return a, b, U

if __name__ == "__main__":
    f = sys.argv[1]
    K = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    dps = int(sys.argv[3]) if len(sys.argv) > 3 else 60
    R = json.load(open(f))
    for rec in R["best"][:K]:
        a, b, U = load(rec)
        h = hp_margin(a, b, U, dps)
        print("restart %s start %s float_m %.3e hp_m %s udef %s nv %d dp %.3g db %.3g sep %s vd %.2e label %s" % (
            rec["restart"], rec["start"], rec["m"], mp.nstr(h["m"], 8), mp.nstr(h["unitarity_defect"], 3),
            h["n_vertices"], rec["dist_genperm"], rec["dist_blockdiag"],
            [round(s, 4) for s in rec["sep"]], rec["vertex_dist"], rec["label"]))
