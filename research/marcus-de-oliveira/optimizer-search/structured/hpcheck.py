"""Re-evaluate at 50 digits (lib/moc_hp.hp_margin, all n! points, any n) every saved
record whose float margin exceeds 1e-7, plus the best non-degenerate guarded record per file.
Usage: python hpcheck.py OUT.json RUN.json [RUN.json ...]"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fixedu  # noqa: E402,F401  (puts lib on the path)
from moc_hp import hp_margin  # noqa: E402

out = []
t0 = time.time()
for fn in sys.argv[2:]:
    D = json.load(open(fn))
    units = {(u["family"], u["n"], u["instance"]): np.array(u["U_re"]) + 1j * np.array(u["U_im"]) for u in D["units"]}
    nd = [r for r in D["best"] if r["mode"] == "guard" and r["hull_width"] >= 1e-3]
    todo = [r for r in D["best"] if r["m"] > 1e-7] + nd[:1]
    seen = set()
    for r in todo:
        key = (r["family"], r["n"], r["instance"], r["restart"], r["mode"])
        if key in seen:
            continue
        seen.add(key)
        U = units[(r["family"], r["n"], r["instance"])]
        a = np.array([complex(*v) for v in r["a"]])
        b = np.array([complex(*v) for v in r["b"]])
        h = hp_margin(a, b, U, dps=50)
        rec = dict(file=os.path.basename(fn), family=r["family"], n=r["n"], instance=r["instance"],
                   restart=r["restart"], mode=r["mode"], m_float=r["m"], m_hp50=h["m_str"],
                   hull_width=r["hull_width"], max_abs_a=float(np.abs(a).max()), max_abs_b=float(np.abs(b).max()),
                   unitarity_defect_hp=str(h["unitarity_defect"]))
        out.append(rec)
        print(rec, "t=%.0f" % (time.time() - t0), flush=True)
json.dump(out, open(sys.argv[1], "w"), indent=1)
