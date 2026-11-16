"""Convert the structured agent outputs stored in a workflow journal (journal.jsonl) into Markdown review records.
Usage:
  python journal_to_md.py review <workflow_dir> <cycle> <out.md>      # 10-reviewer cycle (labels cN:R*, cN:consolidate)
  python journal_to_md.py audit  <workflow_dir> <claims_out.md> <citations_out.md>
Agent outputs are copied verbatim (no editing); the raw JSON is saved next to the Markdown file for traceability."""
import json, os, sys


def load(wdir):
    labels, results = {}, {}
    for line in open(os.path.join(wdir, "journal.jsonl")):
        j = json.loads(line)
        if j.get("type") == "started":
            labels[j["agentId"]] = j.get("label")
        elif j.get("type") == "result":
            r = j.get("result")
            if isinstance(r, str):
                try:
                    r = json.loads(r)
                except Exception:
                    pass
            results[labels.get(j["agentId"], j["agentId"])] = r
    return results


def esc(s):
    return str(s).replace("|", "\\|").replace("\n", " ")


SEV = {"fatal": 0, "critical": 1, "major": 2, "moderate": 3, "minor": 4}


def review_md(res, cycle):
    revs = {k: v for k, v in res.items() if k and k.startswith(f"c{cycle}:R") and isinstance(v, dict)}
    cons = res.get(f"c{cycle}:consolidate")
    L = [f"# Review cycle {cycle:02d}", "",
         f"Ten independent reviewer agents (R1–R10, distinct lenses) reviewed the manuscript, Supplementary Material, code and result files; "
         f"a consolidation agent merged and prioritised their issues. Agent outputs are reproduced verbatim below; the response and the "
         f"actions taken are recorded at the end of this file.", "",
         "## Recommendations", "", "| Reviewer | Recommendation | Issues (fatal/critical/major/moderate/minor) |", "|---|---|---|"]
    for k in sorted(revs):
        v = revs[k]; iss = v.get("issues", [])
        cnt = "/".join(str(sum(1 for i in iss if i.get("severity") == s)) for s in SEV)
        L.append(f"| {esc(v.get('reviewer', k))} | {esc(v.get('recommendation', ''))} | {cnt} |")
    missing = [f"R{i}" for i in range(1, 11) if not any(k.startswith(f"c{cycle}:R{i}_") for k in revs)]
    if missing:
        L += ["", f"Reviews not returned in this cycle: {', '.join(missing)}."]
    if cons:
        L += ["", "## Consolidated action list", "", f"**Overall:** {cons.get('overall', '')}", "",
              "| ID | Priority | Problem | Raised by | Fix | Experiment | Feasible now |", "|---|---|---|---|---|---|---|"]
        for a in sorted(cons.get("actions", []), key=lambda a: SEV.get(a.get("priority"), 9)):
            L.append(f"| {esc(a.get('action_id'))} | {esc(a.get('priority'))} | {esc(a.get('problem'))} | {esc(', '.join(a.get('raised_by', [])))} | "
                     f"{esc(a.get('fix'))} | {esc(a.get('experiment_needed'))} | {a.get('feasible_now')} |")
    ver = {k: v for k, v in res.items() if k and k.startswith(f"c{cycle}:verify:") and isinstance(v, dict)}
    if ver:
        L += ["", "## Adversarial verification of fatal/critical/major actions", "", "| Action | Verdict | Evidence | Minimal correct fix |", "|---|---|---|---|"]
        for k in sorted(ver):
            v = ver[k]
            L.append(f"| {esc(v.get('action_id', k))} | {esc(v.get('verdict'))} | {esc(v.get('evidence'))} | {esc(v.get('minimal_correct_fix'))} |")
    L += ["", "## Individual reviews"]
    for k in sorted(revs):
        v = revs[k]
        L += ["", f"### {v.get('reviewer', k)} — {v.get('recommendation', '')}", "", v.get("summary", ""), ""]
        for i in sorted(v.get("issues", []), key=lambda i: SEV.get(i.get("severity"), 9)):
            L += [f"* **[{i.get('severity')}] {i.get('id')}: {i.get('title')}** ({i.get('location')}; confidence: {i.get('confidence')})",
                  f"  * Evidence: {i.get('evidence')}", f"  * Corrective action: {i.get('corrective_action')}"]
            if i.get("required_experiment") and str(i.get("required_experiment")).strip().lower() not in ("none", "n/a", "-", ""):
                L.append(f"  * Required experiment: {i.get('required_experiment')}")
        if v.get("verified_previous_fixes"):
            L += ["", "Verified previous fixes: " + "; ".join(v["verified_previous_fixes"])]
    L += ["", "## Responses and actions taken", "", "(to be completed after implementation)", ""]
    return "\n".join(L)


def audit_md(res):
    cl = res.get("audit:claims") or {}; ci = res.get("audit:citations") or {}
    A = ["# Claim–evidence audit (manuscript v1)", "", f"Claims checked: {cl.get('n_claims_checked', '?')}. {cl.get('summary', '')}", "",
         "| ID | Section | Claim (quote) | Evidence | Status | Problem | Fix |", "|---|---|---|---|---|---|---|"]
    order = {"CONTRADICTED": 0, "UNSUPPORTED": 1, "OVERSTATED": 2, "PARTIALLY_SUPPORTED": 3, "SUPPORTED": 4}
    for c in sorted(cl.get("claims", []), key=lambda c: order.get(c.get("status"), 9)):
        A.append("| " + " | ".join(esc(c.get(k, "")) for k in ["id", "section", "quote", "evidence", "status", "problem", "fix"]) + " |")
    A += ["", "## Resolution", "", "(to be completed after implementation)", ""]
    B = ["# Citation audit (manuscript v1)", "", ci.get("summary", ""), "", "## Reference integrity", "",
         "| Key | DOI | Resolves | Metadata OK | Correction/retraction | Problem | Fix |", "|---|---|---|---|---|---|---|"]
    for r in ci.get("references", []):
        B.append("| " + " | ".join(esc(r.get(k, "")) for k in ["key", "doi", "doi_resolves", "metadata_ok", "retraction_or_correction", "problem", "fix"]) + " |")
    B += ["", "## Use of citations", "", "| Key | Statement | Status | Basis | Fix |", "|---|---|---|---|---|"]
    order = {"UNSUPPORTED": 0, "PARTIALLY_SUPPORTED": 1, "UNVERIFIABLE": 2, "SUPPORTED": 3}
    for u in sorted(ci.get("citation_uses", []), key=lambda u: order.get(u.get("status"), 9)):
        B.append("| " + " | ".join(esc(u.get(k, "")) for k in ["key", "statement", "status", "basis", "fix"]) + " |")
    B += ["", "## Missing citations / competing studies", ""] + [f"* {m}" for m in ci.get("missing_citations", [])]
    B += ["", "## Resolution", "", "(to be completed after implementation)", ""]
    return "\n".join(A), "\n".join(B)


if __name__ == "__main__":
    mode, wdir = sys.argv[1], sys.argv[2]
    res = load(wdir)
    if mode == "review":
        cycle, out = int(sys.argv[3]), sys.argv[4]
        open(out, "w").write(review_md(res, cycle))
        json.dump(res, open(out.replace(".md", "_raw.json"), "w"), indent=1, ensure_ascii=False)
    else:
        a, b = audit_md(res)
        open(sys.argv[3], "w").write(a); open(sys.argv[4], "w").write(b)
        json.dump(res, open(os.path.join(os.path.dirname(sys.argv[3]), "audit_raw.json"), "w"), indent=1, ensure_ascii=False)
    print("written")
