"""Standard-library certificate replay; exact enumeration uses the pinned lrslib run."""
from pathlib import Path
import importlib.util
import json
import re

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("split_ray_algebra", HERE / "algebra.py")
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)


def expected_input(kernel):
    rows = [[0] + list(k) for k in kernel]
    rows += [[-1] + [1] * 33]
    rows += [[0] + [int(i == j) for i in range(33)] for j in range(33)]
    return ("split_dual_normalized\nH-representation\nlinearity 13 "
            + " ".join(map(str, range(1, 14))) + "\nbegin\n46 34 rational\n"
            + "\n".join(" ".join(map(str, r)) for r in rows)
            + "\nend\nmaxoutput 10000\n")


def validate_output(path):
    text = path.read_text()
    a.require("lrsgmp.h,mini-gmp" in text, "arbitrary precision tool transcript")
    totals = re.findall(r"\*Totals: vertices=(\d+) rays=(\d+) bases=(\d+)", text)
    a.require(totals == [("935", "0", "199641")], "completed enumeration totals")
    a.require("\nend\n" in text, "terminated V representation")
    a.require("*elapsed time:" in text, "normal tool footer")
    lines = text.splitlines()
    start = lines.index("begin")
    end = lines.index("end", start)
    a.require(lines[start + 1] == "***** 34 rational", "33 dimensional rational V format")
    data = lines[start + 2:end]
    a.require(len(data) == 935, "935 output vertices, below 10000 guard")
    a.require(all(len(s.split()) == 34 and s.split()[0] == "1" for s in data),
              "only normalized finite vertices")
    rays = a.read_lrs(path)
    a.require(len(rays) == 935, "no duplicate normalized vertices")
    return rays


def check():
    tables, keys, matrix, kernel = a.structure()
    expected = dict(schema="mdo.split-dual-system.v1", row_colors=list(a.ROW),
                    column_colors=list(a.COL), tables=tables, monomials=keys,
                    matrix=matrix, kernel=kernel)
    # Normalize tuples to JSON arrays before comparing archived independent data.
    a.require(json.loads(json.dumps(expected)) == json.loads((HERE / "system.json").read_text()),
              "independent coefficient system reconstruction")
    a.require((HERE / "bounded.ine").read_text() == expected_input(kernel),
              "exact normalized H input")
    rays = validate_output(HERE / "lrs-output.txt")
    stored = json.loads((HERE / "rays.json").read_text())
    a.require(stored == [list(r) for r in sorted(rays)], "complete sorted ray archive")
    group = a.actions(tables)
    rows = {tuple(r) for r in matrix}
    a.require(all(a.act(r, g) in rows for r in rows for g in group),
              "physical group preserves the coefficient rows")
    known = {a.primitive(r) for r in matrix}
    partition = {}
    for ray in rays:
        representative = min(a.act(ray, g) for g in group)
        partition.setdefault(representative, set()).add(ray)
    records = []
    for index, (ray, members) in enumerate(sorted(partition.items())):
        a.verify_ray(ray, kernel)
        a.require({a.act(ray, g) for g in group} == members, "full physical orbit")
        inv = a.invariant(ray, tables, group)
        a.require(all(a.invariant(r, tables, group) == inv for r in members),
                  "claimed physical invariants")
        records.append(dict(id=index, ray=list(ray), invariant=inv,
                            grouped_coefficient=bool(members & known)))
    result = dict(rays=len(rays), physical_orbits=len(records),
                  grouped_coefficient_orbits=sum(r["grouped_coefficient"] for r in records),
                  records=records)
    a.require(result == json.loads((HERE / "orbits.json").read_text()), "orbit certificate")
    a.require(len(records) == 62 and result["grouped_coefficient_orbits"] == 10,
              "orbit classification totals")
    return {"tables": 33, "coefficient_rows": 192, "rank": 21, "kernel_rank": 12,
            "rays": 935, "physical_group_order": 24, "physical_orbits": 62,
            "coefficient_orbits": 10, "other_orbits": 52,
            "extremality": "exact support-nullity, all representatives and group images",
            "completeness": "pinned exact lrslib reverse-search completion; rerun command in README",
            "physical_universal_signs": "not established here"}


if __name__ == "__main__":
    print(json.dumps(check(), indent=2))
