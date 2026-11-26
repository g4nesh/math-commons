"""High-precision re-evaluation and classification of saved joint67 configurations.

Usage (from this folder): python hpcheck.py OUT.json RUN1.json [RUN2.json ...]
For each run file: the TOP best configurations (by float m) and every saved
nontrivial near-contact are recomputed with moc_hp.hp_margin at 60 digits and
classified:
  collinear  - all z_sigma on one line (hull width / diam < 1e-8); this happens
               when a and -b lie on one circle or line (Moebius image of the
               Hermitian case), a known zero-margin family
  vertex     - one active permutation: d(U) equals a hull vertex z_sigma
  edge       - two active permutations
  other      - anything else
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
from moc_core import z_all, margin, dist_genperm, dist_blockdiag, min_separation  # noqa
from moc_hp import hp_margin  # noqa

TOP = 5


def load(rec):
    a = np.array([complex(*v) for v in rec["a"]])
    b = np.array([complex(*v) for v in rec["b"]])
    U = np.array(rec["U_re"]) + 1j * np.array(rec["U_im"])
    return a, b, U


def circle_residual(pts):
    """Relative residual of the best algebraic circle/line fit
    A|p|^2 + D x + E y + F = 0 (smallest singular value / largest)."""
    x, y = pts.real, pts.imag
    M = np.column_stack([x * x + y * y, x, y, np.ones_like(x)])
    s = np.linalg.svd(M, compute_uv=False)
    return float(s[-1] / s[0])


def classify(a, b, U):
    r = margin(a, b, U)
    z = z_all(a, b)
    zc = np.column_stack([z.real - z.real.mean(), z.imag - z.imag.mean()])
    sv = np.linalg.svd(zc, compute_uv=False)
    width = float(sv[-1] / sv[0])
    conc = circle_residual(np.concatenate([a, -b]))
    if width < 1e-8:
        kind = "collinear"
    elif len(r["active"]) == 1:
        kind = "vertex"
    elif len(r["active"]) == 2:
        kind = "edge"
    else:
        kind = "other"
    return kind, dict(z_width_ratio=width, concyclic_residual_a_minus_b=conc, n_active=len(r["active"]),
                      active=r["active"][:20], float_m=r["m"], dist_genperm=dist_genperm(U)[0],
                      dist_blockdiag=dist_blockdiag(U)[0], min_sep=min(min_separation(a, b)))


def main():
    out, files = sys.argv[1], sys.argv[2:]
    res = []
    for f in files:
        j = json.load(open(f))
        recs = j["best"][:TOP] + [c for c in j["near_contacts"]
                                  if c["restart"] not in [x["restart"] for x in j["best"][:TOP]]]
        for rec in recs:
            a, b, U = load(rec)
            kind, info = classify(a, b, U)
            h = hp_margin(a, b, U, dps=60)
            item = dict(run=os.path.basename(f), n=len(a), seed=rec["seed"], restart=rec["restart"],
                        start=rec["start"], kind=kind, hp_m_60=h["m_str"],
                        hp_unitarity_defect=str(h["unitarity_defect"])[:12], **info,
                        theta=rec["theta"], a=rec["a"], b=rec["b"], U_re=rec["U_re"], U_im=rec["U_im"])
            res.append(item)
            print(item["run"], item["restart"], kind, "float %.2e" % info["float_m"], "hp", h["m_str"][:14],
                  "dperm %.3f dblock %.3f sep %.3f conc %.1e" % (info["dist_genperm"], info["dist_blockdiag"],
                                                               info["min_sep"], info["concyclic_residual_a_minus_b"]),
                  flush=True)
    json.dump(res, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
