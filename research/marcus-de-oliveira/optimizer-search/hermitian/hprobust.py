"""Usage: python hprobust.py runs/<tag>.json [top] [dps]: 60-digit check of top saved records,
both moc_hp.hp_margin and the hull-independent lower bound from robust.robust_margin_hp."""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from robust import robust_margin, robust_margin_hp
from moc_hp import hp_margin
import mpmath as mp
f = sys.argv[1]; top = int(sys.argv[2]) if len(sys.argv) > 2 else 3; dps = int(sys.argv[3]) if len(sys.argv) > 3 else 60
j = json.load(open(f))
recs = sorted(j["best"], key=lambda r: -r.get("m_robust", r["m"]))[:top]
for r in recs:
    a = np.array([complex(*v) for v in r["a"]]); b = np.array([complex(*v) for v in r["b"]])
    U = np.array(r["U_re"]) + 1j * np.array(r["U_im"])
    rf = robust_margin(a, b, U)
    h = hp_margin(a, b, U, dps=dps)
    lb = robust_margin_hp(a, b, U, [rf["theta"]], dps=dps)
    r["hp_m"] = h["m_str"]; r["hp_lower_bound_at_theta"] = mp.nstr(lb, 20); r["hp_dps"] = dps
    print("%s r=%d robust=%+.4e hp_margin=%s hp f(theta)/diam=%s dgp=%.2e sep=%s" % (
        os.path.basename(f), r["restart"], rf["m"], mp.nstr(h["m"], 8), mp.nstr(lb, 8), r["dist_genperm"],
        ["%.2e" % v for v in r["sep"]]))
json.dump(j, open(f, "w"), indent=1)
