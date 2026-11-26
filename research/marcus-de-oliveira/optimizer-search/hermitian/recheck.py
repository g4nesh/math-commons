"""Recompute every saved record of every run file with robust_margin; print the per-run max."""
import glob, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from robust import robust_margin
rows = []
for f in sorted(glob.glob(os.path.join(HERE, "runs", "n*.json"))):
    j = json.load(open(f))
    best = None
    flips = 0
    for r in j["best"]:
        a = np.array([complex(*v) for v in r["a"]]); b = np.array([complex(*v) for v in r["b"]])
        U = np.array(r["U_re"]) + 1j * np.array(r["U_im"])
        rm = robust_margin(a, b, U)["m"]
        r["m_robust"] = rm
        if (rm > 0) != (r["m"] > 0) and abs(r["m"]) > 1e-10:
            flips += 1
        if best is None or rm > best[0]:
            best = (rm, r)
    s = j["summary"]
    print("%-18s done=%4d  run max_m=%+.3e  robust max over saved=%+.3e (r=%d dgp=%.3f dbd=%.3f dz=%.1e sep=%.3f,%.3f) sign-flips=%d" % (
        os.path.basename(f)[:-5], s["restarts_done"], s.get("max_m", s.get("max_m_robust")), best[0], best[1]["restart"],
        best[1]["dist_genperm"], best[1]["dist_blockdiag"], best[1]["d_to_nearest_z"], best[1]["sep"][0], best[1]["sep"][1], flips))
    json.dump(j, open(f, "w"), indent=1)
