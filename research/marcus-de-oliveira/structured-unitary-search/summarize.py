"""Print the per-family Markdown table of README.md from results/n*.json.  Usage: python -I -B summarize.py"""
import json, glob, os
HERE = os.path.dirname(os.path.abspath(__file__))


def cell(t):
    if not t:
        return "none feasible", ""
    sep = min(t["sep_a"], t["sep_b"], t["sep_ab"])
    which = min(("a-a", t["sep_a"]), ("b-b", t["sep_b"]), ("a+b", t["sep_ab"]), key=lambda u: u[1])[0]
    hp = (", 50d " + "%.3g" % float(t["hp50_m"])) if "hp50_m" in t else ""
    return "%.2e%s" % (t["m"], hp), "%d active; min sep %.1e (%s)" % (t["n_active"], sep, which)


print("| family | n | params | restarts | best m, free spectra | contact (free) | best m, separated spectra | contact (separated) |")
print("|---|---|---|---|---|---|---|---|")
for path in sorted(glob.glob(os.path.join(HERE, "results", "n*.json"))):
    for name, f in json.load(open(path))["families"].items():
        mf, cf = cell(f["best_free"][0] if f["best_free"] else None)
        mp, cp = cell(f["best_penalized_feasible"][0] if f["best_penalized_feasible"] else None)
        print("| %s | %d | %d | %d | %s | %s | %s | %s |" % (name, f["n"], f["params"], f["restarts"], mf, cf, mp, cp))
