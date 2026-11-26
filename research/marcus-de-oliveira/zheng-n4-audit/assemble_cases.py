#!/usr/bin/env python3
"""Assemble the status of the ninety contact cases from the audit's own files.

Reads (all inside this directory, nothing from the supplement):
  manifest.json                  document ids and PDF paths
  deps/dependency-map.json       per case: contact set, group, stage, primary and branch documents, prerequisites
  docs/<id>.json                 first-pass report of each proof document
  docs2/<id>__gapcheck.md        a flagged step settled by the second pass
  docs2/<id>__certificates.md    large certificates replayed after the first pass
Writes:
  cases/status.json              machine-readable status per case and per auxiliary statement
  cases/TABLE.md                 the tables
and, if PROOF.md exists, replaces the text between the marker pairs
  <!-- BEGIN GENERATED STATUS SUMMARY --> ... <!-- END GENERATED STATUS SUMMARY --> and
  <!-- BEGIN GENERATED CASE STATUS --> ... <!-- END GENERATED CASE STATUS -->.

Run from this directory:  python3 -I -B assemble_cases.py
A case is "first-passed" when every document of its closure in the supplement's catalogue has a
first-pass report with verdict confirmed or confirmed-modulo-citations, counting a document whose
flagged steps were settled, or whose certificates were replayed, by the second pass. This is a
bookkeeping status, not a proof: see METHOD.md for what a verdict means.
"""
import glob
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))


def load(path):
    with open(os.path.join(HERE, path), encoding="utf-8") as fh:
        return json.load(fh)


def main():
    manifest = load("manifest.json")
    deps = load("deps/dependency-map.json")
    pdf_to_id = {}
    for batch in manifest["batches"]:
        for doc in batch["docs"]:
            pdf_to_id[doc["pdf"]] = doc["id"]
    all_ids = sorted(pdf_to_id.values())

    report = {}
    for path in glob.glob(os.path.join(HERE, "docs", "*.json")):
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
        report[os.path.basename(path)[:-5]] = data
    settled = {os.path.basename(p)[:-len("__gapcheck.md")] for p in glob.glob(os.path.join(HERE, "docs2", "*__gapcheck.md"))}
    replayed = {os.path.basename(p)[:-len("__certificates.md")] for p in glob.glob(os.path.join(HERE, "docs2", "*__certificates.md"))}

    def doc_state(doc_id):
        """Return (state, note): state in ok / flagged / incomplete / missing."""
        data = report.get(doc_id)
        if data is None:
            return "missing", "no first-pass report"
        verdict = data.get("verdict")
        real_gaps = [g for g in (data.get("gaps") or []) if isinstance(g, dict) and g.get("severity") in ("substantive", "fatal")]
        if verdict in ("confirmed", "confirmed-modulo-citations") and not real_gaps:
            return "ok", verdict
        if verdict == "could-not-check":
            if doc_id in replayed:
                return "ok", "certificates replayed in the second pass"
            return "incomplete", "certificates not replayed"
        if doc_id in settled:
            return "ok", "flagged step settled in the second pass"
        if verdict == "gap-found" and not real_gaps:
            return "ok", "presentation remarks only"
        return "flagged", "flagged step not yet settled"

    states = {i: doc_state(i) for i in all_ids}
    counts = {}
    for st, _ in states.values():
        counts[st] = counts.get(st, 0) + 1

    def short(node):
        return node.replace("T^(7)_", "S").replace("B_", "B:")

    cases = []
    for c in sorted(deps["cases"], key=lambda c: int(c["case"].split("_")[1])):
        k = int(c["case"].split("_")[1])
        own = [pdf_to_id[d["document"]] for d in c.get("own_documents", []) if d["document"] in pdf_to_id]
        prim = pdf_to_id.get(c["primary"]["document"])
        if prim and prim not in own:
            own.insert(0, prim)
        closure = list(own)
        for d in c.get("catalog_branches", []):
            i = pdf_to_id.get(d["document"])
            if i and i not in closure:
                closure.append(i)
        missing = [i for i in closure if states[i][0] == "missing"]
        flagged = [i for i in closure if states[i][0] == "flagged"]
        incomplete = [i for i in closure if states[i][0] == "incomplete"]
        notes = [i for i in closure if states[i][0] == "ok" and states[i][1] not in ("confirmed", "confirmed-modulo-citations")]
        if missing:
            status = "awaiting first pass"
        elif flagged or incomplete:
            status = "first pass has an open item"
        else:
            status = "first pass complete"
        cases.append({
            "case": k, "contact_set": c["contact_set"], "group": c.get("paper_group") or c.get("group"),
            "stage": c.get("stage_paper"), "primary": prim, "primary_page": c["primary"].get("page"),
            "own_documents": own, "closure": closure,
            "prerequisites": [short(x) for x in c.get("graph_direct_prerequisites", [])],
            "missing": missing, "flagged": flagged, "incomplete": incomplete, "settled_or_replayed": notes,
            "status": status,
        })

    aux = []
    for a in deps.get("auxiliary_and_branch_nodes", []):
        docs = []
        for p in a.get("proving_documents", []):
            path = p.get("document") if isinstance(p, dict) else p
            if path in pdf_to_id and pdf_to_id[path] not in docs:
                docs.append(pdf_to_id[path])
        missing = [i for i in docs if states[i][0] == "missing"]
        bad = [i for i in docs if states[i][0] in ("flagged", "incomplete")]
        aux.append({"node": short(a["node"]), "contacts": a.get("contacts"), "stage": a.get("stage_paper"),
                    "documents": docs, "missing": missing, "open": bad,
                    "status": "awaiting first pass" if missing else ("first pass has an open item" if bad else ("first pass complete" if docs else "no document identified"))})

    out = {"documents": {"total": len(all_ids), **counts},
           "cases": cases, "auxiliary": aux,
           "documents_without_report": sorted(i for i in all_ids if states[i][0] == "missing"),
           "documents_open": sorted(i for i in all_ids if states[i][0] in ("flagged", "incomplete"))}
    os.makedirs(os.path.join(HERE, "cases"), exist_ok=True)
    with open(os.path.join(HERE, "cases", "status.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)

    done = [c for c in cases if c["status"] == "first pass complete"]
    waiting = [c for c in cases if c["status"] == "awaiting first pass"]
    openitem = [c for c in cases if c["status"] == "first pass has an open item"]
    own_missing = [c for c in cases if any(i in c["missing"] for i in c["own_documents"])]

    def tl(cs):
        return ", ".join(f"$T_{{{c['case']}}}$" for c in cs) if cs else "none"

    L = []
    L.append(f"**Documents.** {len(all_ids)} proof documents: {counts.get('ok', 0)} with a first-pass report and no open item, "
             f"{counts.get('flagged', 0) + counts.get('incomplete', 0)} with an open item, {counts.get('missing', 0)} without a first-pass report.")
    L.append("")
    L.append(f"**Cases.** {len(done)} of 90 have a complete first pass over every document of their closure. "
             f"{len(openitem)} have an open item. {len(waiting)} wait for at least one document.")
    L.append("")
    L.append(f"- Cases whose own documents lack a first-pass report: {tl(own_missing)}.")
    L.append(f"- Cases that only wait for a document of another case or a shared document: {tl([c for c in waiting if c not in own_missing])}.")
    L.append(f"- Cases with an open item: {tl(openitem)}.")
    L.append("")
    summary = "\n".join(L) + "\n"
    L = []
    L.append("| Case | Contact set | Group | Stage | Primary proof document | Cited cases and auxiliary statements | Documents (reported / closure) | Status |")
    L.append("| --- | --- | --- | --- | --- | --- | --- | --- |")
    for c in cases:
        n = len(c["closure"])
        rep = n - len(c["missing"])
        extra = ""
        if c["settled_or_replayed"]:
            extra = f"; {len(c['settled_or_replayed'])} settled or replayed in the second pass"
        if c["missing"]:
            own_m = [i for i in c["missing"] if i in c["own_documents"]]
            extra += f"; own documents missing: {len(own_m)}"
        prim = f"[`{c['primary']}`](docs/{c['primary']}.md)" if c["primary"] in report else f"`{c['primary']}`"
        L.append(f"| $T_{{{c['case']}}}$ | {{{','.join(map(str, c['contact_set']))}}} | {c['group']} | {c['stage']} | {prim}, p. {c['primary_page']} | "
                 f"{', '.join(c['prerequisites']) or '—'} | {rep} / {n} | {c['status']}{extra} |")
    L.append("")
    L.append("In the column of cited statements, `T_k` is a contact case, `S` followed by a number or name is one of the 27 seven-element contact sets of the paper's Table 5, and `B:` marks one of the four branch statements.")
    L.append("")
    L.append("| Auxiliary statement | Contact set | Stage | Proving documents | Status |")
    L.append("| --- | --- | --- | --- | --- |")
    for a in aux:
        docs = ", ".join(f"[`{i}`](docs/{i}.md)" if i in report else f"`{i}`" for i in a["documents"]) or "—"
        cs = "{" + ",".join(map(str, a["contacts"])) + "}" if a.get("contacts") else "—"
        L.append(f"| {a['node']} | {cs} | {a['stage']} | {docs} | {a['status']} |")
    table = "\n".join(L) + "\n"
    with open(os.path.join(HERE, "cases", "TABLE.md"), "w", encoding="utf-8") as fh:
        fh.write("# Status of the ninety contact cases\n\nGenerated by `assemble_cases.py`. Do not edit by hand.\n\n" + summary + "\n" + table)

    proof = os.path.join(HERE, "PROOF.md")
    if os.path.exists(proof):
        with open(proof, encoding="utf-8") as fh:
            text = fh.read()
        for name, block in (("STATUS SUMMARY", summary), ("CASE STATUS", table)):
            a, b = f"<!-- BEGIN GENERATED {name} -->", f"<!-- END GENERATED {name} -->"
            if a in text and b in text:
                i, j = text.index(a) + len(a), text.index(b)
                text = text[:i] + "\n" + block + text[j:]
        with open(proof, "w", encoding="utf-8") as fh:
            fh.write(text)
    print(f"documents: {counts}; cases complete {len(done)}, open {len(openitem)}, waiting {len(waiting)}")


if __name__ == "__main__":
    main()
