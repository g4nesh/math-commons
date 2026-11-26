"""Second-pass classification of the saved m ~ 0 configurations that joint5.py
labelled NONTRIVIAL (not a generalized permutation, not block-diagonal, simple
spectra, d not at any z_sigma). Writes contacts/classification.json."""
import glob, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
from moc_core import z_all, hull_indices  # noqa: E402

out, counts = [], {}
for f in sorted(glob.glob(os.path.join(HERE, "runs", "*.json"))):
    R = json.load(open(f))
    for i, rec in enumerate(R["best"]):
        if rec["label"] != "NONTRIVIAL":
            continue
        a = np.array([complex(*v) for v in rec["a"]]); b = np.array([complex(*v) for v in rec["b"]])
        z = z_all(a, b); V = z[hull_indices(z)]
        width = min(max(abs((np.conj(-1j * (V[(k + 1) % len(V)] - V[k]) / abs(V[(k + 1) % len(V)] - V[k])) * (z - V[k])).real))
                    for k in range(len(V))) / rec["diam"]
        S = np.abs(a[:, None] + b[None, :])
        pairs = [[int(i0), int(j0), float(S[i0, j0])] for i0, j0 in zip(*np.nonzero(S < 1e-6))]
        if width < 1e-8:
            cls = "sliver-hull (z nearly collinear: degenerate hull, trivial)"
        elif pairs:
            cls = "zero-factor edge contact (a_i + b_j ~ 0, d on a hull edge through 0)"
        else:
            cls = "unexplained"
        counts[cls] = counts.get(cls, 0) + 1
        out.append(dict(file=os.path.relpath(f, HERE), best_index=i, seed=rec["seed"], restart=rec["restart"],
                        m_float=rec["m"], cls=cls, hull_width_over_diam=width, zero_pairs=pairs,
                        dist_genperm=rec["dist_genperm"], dist_blockdiag=rec["dist_blockdiag"],
                        sep=rec["sep"], vertex_dist=rec["vertex_dist"]))
os.makedirs(os.path.join(HERE, "contacts"), exist_ok=True)
json.dump(dict(counts=counts, records=out), open(os.path.join(HERE, "contacts", "classification.json"), "w"), indent=1)
print(counts)
