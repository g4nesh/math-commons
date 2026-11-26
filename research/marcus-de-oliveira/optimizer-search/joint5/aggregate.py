"""Aggregate runs/*.json into summary.json; for every saved m ~ 0 configuration
labelled NONTRIVIAL, count the z_sigma lying on the line of the hull edge nearest d."""
import glob, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
from moc_core import z_all, hull_indices, perms  # noqa: E402

def edge_line_count(rec, tol=1e-10):
    a = np.array([complex(*v) for v in rec["a"]]); b = np.array([complex(*v) for v in rec["b"]])
    d = complex(*rec["d"]); z = z_all(a, b); V = z[hull_indices(z)]
    best = None
    for k in range(len(V)):
        u, w = V[k], V[(k + 1) % len(V)]
        t = ((d - u) * np.conj(w - u)).real / abs(w - u) ** 2
        q = u + min(1, max(0, t)) * (w - u)
        if best is None or abs(d - q) < best[0]:
            best = (abs(d - q), u, w, t)
    _, u, w, t = best
    nrm = -1j * (w - u) / abs(w - u)
    dist = np.abs((np.conj(nrm) * (z - u)).real) / rec["diam"]
    on = np.nonzero(dist < tol)[0]
    return int(len(on)), float(t), [list(map(int, perms(5)[k])) for k in on]

runs, tot_r, tot_e, labels, best = [], 0, 0, {}, None
for f in sorted(glob.glob(os.path.join(HERE, "runs", "*.json"))):
    R = json.load(open(f)); s = R["summary"]
    tot_r += s["restarts_done"]; tot_e += s["objective_evals"]
    for k, v in s["labels_all"].items():
        labels[k] = labels.get(k, 0) + v
    nontriv = []
    for i, rec in enumerate(R["best"]):
        if rec["label"] == "NONTRIVIAL":
            c, t, on = edge_line_count(rec)
            nontriv.append(dict(index=i, restart=rec["restart"], m=rec["m"], n_z_on_edge_line=c, t=t,
                                dp=rec["dist_genperm"], db=rec["dist_blockdiag"], sep=min(rec["sep"]), vd=rec["vertex_dist"]))
    runs.append(dict(file=os.path.relpath(f, HERE), family=s["family"], seed=s["seed"], opts=s["opts"],
                     restarts=s["restarts_done"], evals=s["objective_evals"], max_m=s["max_m"],
                     n_m_gt_1e7=s["n_m_gt_1e7"], labels=s["labels_all"], seconds=round(s["seconds"]),
                     nontrivial_saved=nontriv))
    if best is None or s["max_m"] > best[0]:
        best = (s["max_m"], os.path.relpath(f, HERE))
hist = {}
for r in runs:
    for x in r["nontrivial_saved"]:
        hist[x["n_z_on_edge_line"]] = hist.get(x["n_z_on_edge_line"], 0) + 1
out = dict(total_restarts=tot_r, total_objective_evals=tot_e, labels=labels, best_float_m=best,
           nontrivial_edge_line_point_counts=hist, runs=runs)
json.dump(out, open(os.path.join(HERE, "summary.json"), "w"), indent=1)
print(json.dumps({k: out[k] for k in ("total_restarts", "total_objective_evals", "labels", "best_float_m", "nontrivial_edge_line_point_counts")}))
for r in runs:
    print(r["file"], r["opts"], r["restarts"], "max_m %.2e" % r["max_m"], r["labels"])
