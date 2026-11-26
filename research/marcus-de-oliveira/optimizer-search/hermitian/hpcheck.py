"""60-digit re-evaluation of the top records of a run file.
Usage: python hpcheck.py runs/<tag>.json [top] [dps]"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
from moc_hp import hp_margin
from moc_core import margin
f = sys.argv[1]
top = int(sys.argv[2]) if len(sys.argv) > 2 else 3
dps = int(sys.argv[3]) if len(sys.argv) > 3 else 60
j = json.load(open(f))
out = []
for r in j["best"][:top]:
    a = np.array([complex(*v) for v in r["a"]])
    b = np.array([complex(*v) for v in r["b"]])
    U = np.array(r["U_re"]) + 1j * np.array(r["U_im"])
    h = hp_margin(a, b, U, dps=dps)
    fm = margin(a, b, U)["m"]
    print("%s r=%d float m=%.6e  hp(%d) m=%s udef=%s" % (os.path.basename(f), r["restart"], fm, dps, h["m_str"],
          str(h["unitarity_defect"])[:10]))
    out.append(dict(restart=r["restart"], float_m=fm, hp_m=h["m_str"], dps=dps))
