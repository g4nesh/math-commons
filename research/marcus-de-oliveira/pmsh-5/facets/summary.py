"""Markdown table of the 34 orbits of hard facets (README), from orbits.json, decompose.json and the logs of
unitary_valid.py, family_min.py and inherit_orbits.py.

usage: summary.py
"""
import glob
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
LOGS = os.path.join(HERE, "logs")


def parse(pattern, regex):
    out = {}
    for fn in glob.glob(os.path.join(LOGS, pattern)):
        for line in open(fn):
            m = re.match(regex, line)
            if m:
                out[m.group(1)] = float(m.group(2))
    return out


def main():
    orbits = json.load(open(os.path.join(HERE, "orbits.json")))["orbits"]
    dec = {r["orbit"]: r for r in json.load(open(os.path.join(HERE, "decompose.json")))["decompositions"]}
    umin = parse("unitary_valid*.txt", r"(\w+): min over U\(5\) ([-+.\de]+)")
    fam = parse("family_min_*.txt", r"orbit (\d+): .* max over U\(5\) - 1 = ([-+.\de]+)")
    grids = {}  # orbit -> {grid: max over R_inh}
    for name, fn in (("82", "inherit_orbits_sym82.txt"), ("32", "inherit_n3sym.txt"), ("50", "inherit_dense.txt")):
        path = os.path.join(LOGS, fn)
        if os.path.exists(path):
            for line in open(path):
                k, _, v = line.split()
                grids.setdefault(int(k), {})[name] = float(v)
    valid = {"A"} | {lab for lab, v in umin.items() if v > -1e-8}  # A: Lemma A (proved)
    print("| orbit | size | found | vertices | overshoot | follows linearly from PMSH(4) "
          "| non-minor Q5 facets (weight) | squared minors | status | family max - 1 |")
    print("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for o in orbits:
        k = o["orbit"]
        terms = dec[k]["terms"]
        kinds = [t["kind"] for t in terms if t["kind"] != "minor"]
        nm = ", ".join(f"{t['kind']} ({t['weight']})" for t in terms if t["kind"] != "minor")
        nminor = sum(t["kind"] == "minor" for t in terms)
        if set(kinds) == {"A"}:
            status = "**proved** (Lemma A)"
        elif set(kinds) <= valid:
            status = "follows from " + " and ".join(sorted(set(kinds) - {"A"})) + " >= 0 (numerical)"
        else:
            bad = sorted(set(kinds) - valid)
            status = "obstacle: " + ", ".join(bad) + (" fails" if len(bad) == 1 else " fail") + " on U(5)"
        res = grids.get(k, {})
        cut = [g for g, v in res.items() if v <= 1 + 1e-6]
        if cut:
            inh = "yes (" + ", ".join(f"{g} shifts" for g in sorted(cut, key=int)) + ")"
        elif "inherited" in o["status"]:
            inh = "yes (a member, first N3 run)"
        else:
            inh = "not cut: " + ", ".join(f"{v:.3f} ({g})" for g, v in sorted(res.items(), key=lambda x: int(x[0])))
        fm = f"{fam[str(k)]:+.1e}" if str(k) in fam else "not run"
        print(f"| {k} | {o['size']} | {o['harvested']} | {o['vertices']} | {o['relaxation'][0] - 1:.3f} | {inh} "
              f"| {nm} | {nminor} | {status} | {fm} |")
    print()
    print("Minimum over U(5) of each non-minor facet of Q5 (unitary_valid.py, 40 starts):",
          ", ".join(f"{lab} {umin[lab]:+.4f}" for lab in sorted(umin, key=lambda x: (len(x), x))))


if __name__ == "__main__":
    main()
