"""Tau continuation on saved constrained records: does the constrained optimum approach m = 0
(a contact away from the trivial sets) as the log-sum-exp temperature goes to 0?
Usage: python cont.py runs/<tag>.json top out.json"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from hsearch import Cls, s_smooth, record
f, top, outp = sys.argv[1], int(sys.argv[2]), sys.argv[3]
j = json.load(open(f)); s = j["summary"]
C = Cls(s["cls"], s["n"], k=s["k"], mult=tuple(s["mult"]) if s["mult"] else None, R=s["R"])
opts = {k: v for k, v in s["opts"].items() if k != "maxiter"}
res = []
for r in j["best"][:top]:
    y = np.array(r["y"]); U = np.array(r["U_re"]) + 1j * np.array(r["U_im"])
    row = [r["m_robust"]]
    for taus in ((3e-5,), (1e-5,), (3e-6,), (1e-6,)):
        y, U, th, _ = s_smooth(C, y, U, taus=taus, maxiter=600, **opts)
        rec = record(C, y, U, dict(seed=s["seed"], restart=r["restart"], m_smooth=None, m_polish=None, taus_done=taus))
        row.append(rec["m_robust"])
    a = rec
    print("%s r=%d m(tau=1e-4,3e-5,1e-5,3e-6,1e-6) = %s  dgp=%.3f dbd=%.3f dz=%.2e nact=%d sep=%s" % (
        os.path.basename(f), r["restart"], " ".join("%.2e" % v for v in row), a["dist_genperm"], a["dist_blockdiag"],
        a["d_to_nearest_z"], a["n_active"], ["%.3f" % v for v in a["sep"]]))
    rec["continuation_m"] = row
    res.append(rec)
json.dump(dict(source=os.path.basename(f), summary=s, records=res), open(outp, "w"), indent=1)
