#!/usr/bin/env python3
"""ACA journal-quality benchmark (Analytica Chimica Acta, ISSN 0003-2670 / 1873-4324).

Retrieves, ranks and tabulates ACA research articles (publication years 2023-2026) in the data-treatment
themes relevant to the planned paper (distribution-free / conformal retention-time prediction intervals
for LC-HRMS non-target screening on a new chromatographic system with few calibrants).

Data sources (no ScienceDirect scraping; it returns HTTP 403):
  * NCBI PubMed E-utilities (esearch/efetch, retmode=xml): titles, abstracts, keywords, DOIs, PMC ids,
    publication types and article-history dates (received / revised / accepted).
  * Crossref REST API /works/<doi>: bibliographic verification and `is-referenced-by-count`.
  * PMC full text (efetch db=pmc) only for selected papers that have a PMC id.

Stages (each caches its raw downloads under ../data/aca_benchmark_raw/, so re-runs are offline unless --force):
  search    PubMed esearch per theme (ACA, 2023:2026[dp], research articles only) -> esearch_<theme>.json
  fetch     PubMed efetch XML for the union of PMIDs                         -> pubmed_efetch_<k>.xml
  crossref  Crossref /works/<doi> for every candidate                         -> crossref_works.jsonl
  rank      parse + merge all candidates                                      -> candidates_all.csv,
                                                                                 top_cited_by_theme.json
  coverage  completeness checks: PubMed vs Crossref yearly totals and Crossref's 400 most-cited ACA articles
            flagged by membership in the candidate pool                      -> coverage_check_*.json,
                                                                                 crossref_top400_flagged.csv
  context   PubMed records of highly cited on-theme tutorials (excluded as non-research articles)
                                                                              -> pubmed_efetch_context_tutorials.xml
  pmc       PMC full text for selected papers that have a PMC id              -> pmc_<PMCID>.xml
  build     join selected papers with the manual coding sheet
            (../data/aca_benchmark_manual_coding.json; coded by reading each abstract / PMC full text)
            -> ../aca_benchmark_papers.csv and aca_benchmark_raw/summary_stats.json

Usage:  python aca_benchmark.py [search|fetch|crossref|rank|pmc|build|all] [--force]

Notes
  * Requests carry a generic User-Agent only (no names / e-mail addresses).
  * NCBI calls sleep 0.4 s; Crossref calls sleep 0.25 s (public-pool etiquette).
  * `received_to_accepted_days` = PubMed History accepted - received (calendar days).
  * Crossref citation counts are a snapshot (date stored in crossref_works.jsonl) and are age-biased:
    older papers have had more time to be cited, so citations/year since online publication are also given.
"""
from __future__ import annotations

import csv
import datetime as dt
import json
import re
import statistics
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import requests

BASE = Path(__file__).resolve().parents[1]            # .../journal_selection
RAW = BASE / "data" / "aca_benchmark_raw"
CODING = BASE / "data" / "aca_benchmark_manual_coding.json"
OUT_CSV = BASE / "aca_benchmark_papers.csv"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
CROSSREF = "https://api.crossref.org/works/"
UA = {"User-Agent": "aca-benchmark-script/1.0 (python-requests; literature metadata study)"}
NCBI_SLEEP, CROSSREF_SLEEP = 0.4, 0.25

JOURNAL = '"Anal Chim Acta"[ta]'
YEARS = "2023:2026[dp]"
EXCL = (' AND journal article[pt] NOT (review[pt] OR "systematic review"[pt] OR editorial[pt] OR comment[pt] '
        'OR letter[pt] OR erratum[pt] OR "published erratum"[pt] OR "retraction of publication"[pt])')

# Theme queries (title/abstract). Deliberately broad; relevance is decided by manual triage of the ranked list.
THEMES = {
    "rt_qsrr": (
        '("retention time"[tiab] OR "retention times"[tiab] OR "retention index"[tiab] OR "retention indices"[tiab] '
        'OR "retention factor"[tiab] OR "retention factors"[tiab] OR "retention behavior"[tiab] '
        'OR "retention behaviour"[tiab] OR QSRR[tiab] OR "structure-retention"[tiab] '
        'OR "quantitative structure retention"[tiab] OR "collision cross section"[tiab] '
        'OR "collision cross-section"[tiab]) AND (predict*[tiab] OR "machine learning"[tiab] OR QSRR[tiab] '
        'OR "in silico"[tiab] OR model*[tiab])'),
    "nts_annotation": (
        '("non-target"[tiab] OR "nontarget"[tiab] OR "non-targeted"[tiab] OR "nontargeted"[tiab] '
        'OR "suspect screening"[tiab] OR "suspect analysis"[tiab] OR "untargeted"[tiab] OR "annotation"[tiab] '
        'OR "structure elucidation"[tiab] OR "structural elucidation"[tiab] OR "molecular networking"[tiab] '
        'OR "identification confidence"[tiab] OR "compound identification"[tiab]) '
        'AND ("high-resolution mass spectrometry"[tiab] OR "high resolution mass spectrometry"[tiab] OR HRMS[tiab] '
        'OR "LC-HRMS"[tiab] OR "LC-MS"[tiab] OR "LC-MS/MS"[tiab] OR Orbitrap[tiab] OR QTOF[tiab] OR "Q-TOF"[tiab] '
        'OR "mass spectrometry"[tiab])'),
    "ml_dl": (
        '("machine learning"[tiab] OR "deep learning"[tiab] OR "neural network"[tiab] OR "neural networks"[tiab] '
        'OR "random forest"[tiab] OR "support vector"[tiab] OR "artificial intelligence"[tiab] '
        'OR "gradient boosting"[tiab] OR XGBoost[tiab] OR convolutional[tiab] OR transformer[tiab] '
        'OR "graph neural"[tiab] OR "Gaussian process"[tiab])'),
    "chemometrics": (
        '(chemometric*[tiab] OR "multivariate data analysis"[tiab] OR "multivariate analysis"[tiab] '
        'OR "multivariate calibration"[tiab] OR "multivariate curve resolution"[tiab] OR "MCR-ALS"[tiab] '
        'OR PARAFAC[tiab] OR "partial least squares"[tiab] OR ASCA[tiab] OR "variable selection"[tiab] '
        'OR "data fusion"[tiab])'),
    "uncertainty_validation": (
        '("prediction interval"[tiab] OR "prediction intervals"[tiab] OR "conformal prediction"[tiab] '
        'OR conformal[tiab] OR "uncertainty quantification"[tiab] OR "measurement uncertainty"[tiab] '
        'OR "uncertainty estimation"[tiab] OR "predictive uncertainty"[tiab] OR "applicability domain"[tiab] '
        'OR "false discovery rate"[tiab] OR "confidence interval"[tiab] OR "confidence intervals"[tiab] '
        'OR "Monte Carlo"[tiab] OR Bayesian[tiab] OR bootstrap*[tiab])'),
    "transfer_domain_shift": (
        '("calibration transfer"[tiab] OR "model transfer"[tiab] OR "transfer learning"[tiab] '
        'OR "domain adaptation"[tiab] OR "domain shift"[tiab] OR "batch effect"[tiab] OR "batch effects"[tiab] '
        'OR "batch correction"[tiab] OR "signal drift"[tiab] OR "instrumental drift"[tiab] '
        'OR "inter-laboratory"[tiab] OR interlaboratory[tiab] OR "between-laboratory"[tiab] '
        'OR transferability[tiab] OR "retention time alignment"[tiab] OR "retention time correction"[tiab] '
        'OR "retention time mapping"[tiab] OR "retention time projection"[tiab] OR "retention time transfer"[tiab] '
        'OR "retention time shift"[tiab] OR "retention time shifts"[tiab])'),
}


# ----------------------------------------------------------------------------------------------- helpers
def http_get(url, params=None, retries=5, sleep=0.4, ok404=False):
    for k in range(retries):
        try:
            r = requests.get(url, params=params, headers=UA, timeout=90)
            if r.status_code == 200:
                time.sleep(sleep)
                return r
            if r.status_code == 404 and ok404:
                time.sleep(sleep)
                return None
            if r.status_code in (429, 500, 502, 503, 504):
                time.sleep(3 * (k + 1))
                continue
            print(f"  HTTP {r.status_code} for {url}", file=sys.stderr)
        except requests.RequestException as e:
            print(f"  request error {e!r}", file=sys.stderr)
        time.sleep(2 * (k + 1))
    raise RuntimeError(f"GET failed: {url} {params}")


def http_post(url, data, retries=5, sleep=0.4):
    for k in range(retries):
        try:
            r = requests.post(url, data=data, headers=UA, timeout=180)
            if r.status_code == 200:
                time.sleep(sleep)
                return r
            print(f"  HTTP {r.status_code} for {url}", file=sys.stderr)
        except requests.RequestException as e:
            print(f"  request error {e!r}", file=sys.stderr)
        time.sleep(3 * (k + 1))
    raise RuntimeError(f"POST failed: {url}")


def text_of(node):
    return re.sub(r"\s+", " ", "".join(node.itertext())).strip() if node is not None else ""


def hist_date(node):
    try:
        return dt.date(int(node.findtext("Year")), int(node.findtext("Month")), int(node.findtext("Day")))
    except Exception:
        return None


def cr_date(msg, key):
    parts = (msg.get(key) or {}).get("date-parts") or [[None]]
    p = parts[0] + [1, 1]
    try:
        return dt.date(int(p[0]), int(p[1] or 1), int(p[2] or 1))
    except Exception:
        return None


# ----------------------------------------------------------------------------------------------- stages
def stage_search(force=False):
    RAW.mkdir(parents=True, exist_ok=True)
    for theme, q in THEMES.items():
        f = RAW / f"esearch_{theme}.json"
        if f.exists() and not force:
            print(f"[search] cached {f.name}")
            continue
        term = f"{JOURNAL} AND {YEARS} AND {q}{EXCL}"
        r = http_get(EUTILS + "esearch.fcgi", {"db": "pubmed", "term": term, "retmax": 5000, "retmode": "json",
                                                "tool": "aca_benchmark_script"}, sleep=NCBI_SLEEP)
        js = r.json()
        js["_query_term"] = term
        js["_retrieved_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        f.write_text(json.dumps(js, indent=1))
        print(f"[search] {theme}: count={js['esearchresult']['count']}")


def union_pmids():
    themes_of = {}
    for theme in THEMES:
        js = json.loads((RAW / f"esearch_{theme}.json").read_text())
        for pmid in js["esearchresult"]["idlist"]:
            themes_of.setdefault(pmid, []).append(theme)
    return themes_of


def stage_fetch(force=False):
    pmids = sorted(union_pmids(), key=int)
    existing = sorted(RAW.glob("pubmed_efetch_*.xml"))
    if existing and not force:
        print(f"[fetch] cached {len(existing)} batch files")
        return
    for old in existing:
        old.unlink()
    for k, i in enumerate(range(0, len(pmids), 200)):
        r = http_post(EUTILS + "efetch.fcgi", {"db": "pubmed", "id": ",".join(pmids[i:i + 200]), "retmode": "xml",
                                                "tool": "aca_benchmark_script"}, sleep=NCBI_SLEEP)
        (RAW / f"pubmed_efetch_{k:02d}.xml").write_bytes(r.content)
        print(f"[fetch] batch {k}: {len(pmids[i:i + 200])} PMIDs")


def parse_pubmed():
    recs = {}
    for f in sorted(RAW.glob("pubmed_efetch_*.xml")):
        root = ET.fromstring(f.read_bytes())
        for art in root.findall(".//PubmedArticle"):
            mc = art.find("MedlineCitation")
            pmid = mc.findtext("PMID")
            a = mc.find("Article")
            abst = []
            for t in a.findall("Abstract/AbstractText"):
                lab = t.get("Label")
                abst.append((f"{lab}: " if lab else "") + text_of(t))
            ids = {x.get("IdType"): (x.text or "").strip() for x in art.findall("PubmedData/ArticleIdList/ArticleId")}
            doi = ids.get("doi") or next((text_of(e) for e in a.findall("ELocationID") if e.get("EIdType") == "doi"), "")
            hist = {}
            for p in art.findall("PubmedData/History/PubMedPubDate"):
                hist.setdefault(p.get("PubStatus"), hist_date(p))
            pd = a.find("Journal/JournalIssue/PubDate")
            year = pd.findtext("Year") or (pd.findtext("MedlineDate") or "")[:4]
            ra = (hist["accepted"] - hist["received"]).days if hist.get("received") and hist.get("accepted") else None
            recs[pmid] = {
                "pmid": pmid, "doi": doi.lower(), "pmc": ids.get("pmc", ""), "pii": ids.get("pii", ""),
                "title": text_of(a.find("ArticleTitle")), "abstract": "\n".join(abst),
                "structured_abstract": any(t.get("Label") for t in a.findall("Abstract/AbstractText")),
                "keywords": "; ".join(text_of(k) for k in mc.findall("KeywordList/Keyword")),
                "pub_types": "; ".join(text_of(p) for p in a.findall("PublicationTypeList/PublicationType")),
                "journal_issue_year": year, "volume": a.findtext("Journal/JournalIssue/Volume") or "",
                "received": str(hist.get("received") or ""), "revised": str(hist.get("revised") or ""),
                "accepted": str(hist.get("accepted") or ""), "received_to_accepted_days": ra,
                "n_authors": len(a.findall("AuthorList/Author")),
                "countries_affil_text": "",  # not used (avoid personal data); kept for schema stability
            }
    return recs


def stage_crossref(force=False):
    recs = parse_pubmed()
    f = RAW / "crossref_works.jsonl"
    done = {}
    if f.exists() and not force:
        for line in f.read_text().splitlines():
            if line.strip():
                o = json.loads(line)
                done[o["doi"]] = o
    todo = [r["doi"] for r in recs.values() if r["doi"] and r["doi"] not in done]
    print(f"[crossref] {len(done)} cached, {len(todo)} to fetch")
    with f.open("a") as fh:
        for n, doi in enumerate(todo, 1):
            r = http_get(CROSSREF + doi, sleep=CROSSREF_SLEEP, ok404=True)
            msg = r.json()["message"] if r is not None else None
            o = {"doi": doi, "retrieved_utc": dt.datetime.now(dt.timezone.utc).isoformat(), "message": msg}
            fh.write(json.dumps(o) + "\n")
            if n % 50 == 0:
                print(f"[crossref] {n}/{len(todo)}", flush=True)


def load_crossref():
    out = {}
    f = RAW / "crossref_works.jsonl"
    if f.exists():
        for line in f.read_text().splitlines():
            if line.strip():
                o = json.loads(line)
                out[o["doi"]] = o
    return out


def merged_candidates():
    recs, cr, themes_of = parse_pubmed(), load_crossref(), union_pmids()
    rows = []
    for pmid, r in recs.items():
        o = cr.get(r["doi"], {})
        m = o.get("message") or {}
        online = cr_date(m, "published-online") or cr_date(m, "published-print") or cr_date(m, "issued")
        snap = (o.get("retrieved_utc") or "")[:10]
        cites = m.get("is-referenced-by-count")
        years_online = None
        if online and snap:
            years_online = max((dt.date.fromisoformat(snap) - online).days / 365.25, 0.25)
        issued = cr_date(m, "published-print") or cr_date(m, "issued")
        rows.append({
            **r, "themes": ";".join(themes_of.get(pmid, [])),
            "crossref_title": (m.get("title") or [""])[0], "crossref_container": (m.get("container-title") or [""])[0],
            "crossref_issn": ";".join(m.get("ISSN") or []), "crossref_type": m.get("type", ""),
            "crossref_volume": m.get("volume", ""), "crossref_article_number": m.get("article-number", ""),
            "crossref_published_online": str(online or ""), "crossref_issue_year": issued.year if issued else "",
            "crossref_citations": cites, "crossref_snapshot_date": snap,
            "citations_per_year": round(cites / years_online, 2) if (cites is not None and years_online) else None,
            "crossref_reference_count": m.get("reference-count"),
            "crossref_license": ";".join(x.get("URL", "") for x in (m.get("license") or [])),
        })
    rows.sort(key=lambda x: (-(x["crossref_citations"] or -1), x["pmid"]))
    return rows


def stage_rank(force=False):
    rows = merged_candidates()
    keys = ["pmid", "doi", "pmc", "crossref_issue_year", "journal_issue_year", "title", "themes", "crossref_citations",
            "citations_per_year", "crossref_published_online", "crossref_snapshot_date", "received", "revised",
            "accepted", "received_to_accepted_days", "pub_types", "structured_abstract", "keywords",
            "crossref_container", "crossref_issn", "crossref_type", "crossref_volume", "crossref_article_number",
            "crossref_reference_count", "crossref_license"]
    with (RAW / "candidates_all.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    top = {}
    for theme in THEMES:
        sel = [r for r in rows if theme in r["themes"].split(";")]
        top[theme] = {"n_candidates": len(sel),
                      "top40_by_crossref_citations": [{k: r[k] for k in ("pmid", "doi", "crossref_issue_year", "title",
                                                                         "crossref_citations", "citations_per_year")}
                                                      for r in sel[:40]]}
    (RAW / "top_cited_by_theme.json").write_text(json.dumps(top, indent=1))
    ra = [r["received_to_accepted_days"] for r in rows if r["received_to_accepted_days"] is not None]
    print(f"[rank] {len(rows)} candidates; with history dates: {len(ra)}; "
          f"median received->accepted {statistics.median(ra) if ra else 'NA'} d")


def stage_pmc(force=False):
    if not CODING.exists():
        print("[pmc] no coding sheet yet; skipping")
        return
    sel = {p["pmid"] for p in json.loads(CODING.read_text())["papers"]}
    recs = parse_pubmed()
    for pmid in sorted(sel):
        pmc = recs.get(pmid, {}).get("pmc")
        if not pmc:
            continue
        f = RAW / f"pmc_{pmc}.xml"
        if f.exists() and not force:
            continue
        r = http_get(EUTILS + "efetch.fcgi", {"db": "pmc", "id": pmc.replace("PMC", ""), "retmode": "xml",
                                               "tool": "aca_benchmark_script"}, sleep=NCBI_SLEEP)
        f.write_bytes(r.content)
        print(f"[pmc] {pmid} -> {pmc} ({len(r.content)} bytes)")


def stage_coverage(force=False):
    """Completeness checks: (1) PubMed vs Crossref yearly ACA totals; (2) Crossref's 400 most-cited ACA journal
    articles (print years 2023-2026), flagged by whether they are in the PubMed candidate pool."""
    f1, f2 = RAW / "coverage_check_pubmed_vs_crossref.json", RAW / "crossref_aca_top400_cited_2023-2026.json"
    if not f1.exists() or force:
        out = {}
        for y in (2023, 2024, 2025, 2026):
            r = http_get(EUTILS + "esearch.fcgi", {"db": "pubmed", "term": f"{JOURNAL} AND {y}[dp]", "retmode": "json",
                                                    "retmax": 0, "tool": "aca_benchmark_script"}, sleep=NCBI_SLEEP)
            r2 = http_get("https://api.crossref.org/journals/0003-2670/works",
                          {"filter": f"from-print-pub-date:{y}-01-01,until-print-pub-date:{y}-12-31,type:journal-article",
                           "rows": 0}, sleep=0.5)
            out[y] = {"pubmed_all_records": int(r.json()["esearchresult"]["count"]),
                      "crossref_journal_articles_print_year": r2.json()["message"]["total-results"]}
        f1.write_text(json.dumps(out, indent=1))
    if not f2.exists() or force:
        p = {"filter": "from-print-pub-date:2023-01-01,until-print-pub-date:2026-12-31,type:journal-article",
             "sort": "is-referenced-by-count", "order": "desc", "rows": 400,
             "select": "DOI,title,is-referenced-by-count,published-print,published-online"}
        js = http_get("https://api.crossref.org/journals/0003-2670/works", p, sleep=0.5).json()
        js["_params"] = p
        f2.write_text(json.dumps(js))
    pool = {r["doi"] for r in parse_pubmed().values()}
    items = json.loads(f2.read_text())["message"]["items"]
    rows = [{"rank_in_ACA_2023_2026": i + 1, "doi": it["DOI"].lower(), "crossref_citations": it["is-referenced-by-count"],
             "in_candidate_pool": it["DOI"].lower() in pool, "title": (it.get("title") or [""])[0]}
            for i, it in enumerate(items)]
    with (RAW / "crossref_top400_flagged.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"[coverage] yearly totals: {f1.read_text()}")
    print(f"[coverage] top-400 cited ACA articles in candidate pool: {sum(r['in_candidate_pool'] for r in rows)}")


# Highly cited on-theme ACA items that the research-article filter excludes (PubMed pt = Review/tutorial);
# identified by screening crossref_top400_flagged.csv. Fetched for context only (not coded as benchmark papers).
CONTEXT_DOIS = ["10.1016/j.aca.2023.341532", "10.1016/j.aca.2025.344838", "10.1016/j.aca.2023.341304",
                "10.1016/j.aca.2022.340711"]


def stage_context(force=False):
    f = RAW / "pubmed_efetch_context_tutorials.xml"
    if f.exists() and not force:
        print("[context] cached")
        return
    pmids = []
    for doi in CONTEXT_DOIS:
        r = http_get(EUTILS + "esearch.fcgi", {"db": "pubmed", "term": f"{doi}[doi]", "retmode": "json",
                                                "tool": "aca_benchmark_script"}, sleep=NCBI_SLEEP)
        pmids += r.json()["esearchresult"]["idlist"]
    r = http_post(EUTILS + "efetch.fcgi", {"db": "pubmed", "id": ",".join(pmids), "retmode": "xml",
                                            "tool": "aca_benchmark_script"}, sleep=NCBI_SLEEP)
    f.write_bytes(r.content)
    print(f"[context] fetched PMIDs {pmids}")


CSV_COLS = ["pmid", "doi", "year", "title", "theme", "research_problem", "research_gap", "novelty_claim",
            "method_contribution", "datasets", "experimental_design", "baselines_compared", "key_quantitative_results",
            "ablation", "statistical_validation", "robustness_eval", "generalization_eval",
            "computational_cost_reported", "limitations_stated", "code_or_data_availability",
            "received_to_accepted_days", "crossref_citations", "why_it_fits_ACA"]


def stage_build(force=False):
    coding = json.loads(CODING.read_text())
    cand = {r["pmid"]: r for r in merged_candidates()}
    out_rows, missing = [], []
    for p in coding["papers"]:
        c = cand.get(p["pmid"])
        if c is None:
            missing.append(p["pmid"])
            continue
        row = {k: p.get(k, "") for k in CSV_COLS}
        # bibliographic + metric fields are always taken from the downloaded metadata, never typed by hand
        row.update({"pmid": c["pmid"], "doi": c["doi"], "year": c["crossref_issue_year"] or c["journal_issue_year"],
                    "title": c["title"], "received_to_accepted_days": c["received_to_accepted_days"],
                    "crossref_citations": c["crossref_citations"]})
        assert c["crossref_container"] == "Analytica Chimica Acta", (c["pmid"], c["crossref_container"])
        out_rows.append(row)
    if missing:
        raise SystemExit(f"PMIDs in coding sheet but not among candidates: {missing}")
    out_rows.sort(key=lambda r: (r["theme"].split(" | ")[0], -(r["crossref_citations"] or 0)))
    with OUT_CSV.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_COLS)
        w.writeheader()
        w.writerows(out_rows)
    summarize(out_rows, cand, coding)
    print(f"[build] wrote {OUT_CSV} ({len(out_rows)} papers)")


def _share(rows, field, prefix):
    n = sum(1 for r in rows if str(r[field]).strip().lower().startswith(prefix))
    return {"n": n, "pct": round(100 * n / len(rows), 1)}


def _quart(v):
    v = [x for x in v if x is not None]
    if len(v) < 4:
        return {"n": len(v)}
    q = statistics.quantiles(v, n=4)
    return {"n": len(v), "median": statistics.median(v), "q1": q[0], "q3": q[2], "min": min(v), "max": max(v)}


def summarize(rows, cand, coding):
    """Tallies over the coded papers. Flags come from the coding sheet ('flags' per paper; abstract-based)."""
    flags = {p["pmid"]: p.get("flags", {}) for p in coding["papers"]}
    primary = lambda r: r["theme"].split(" | ")[0]
    n = len(rows)

    def tally(key, value=True):
        k = sum(1 for r in rows if flags[r["pmid"]].get(key) == value)
        return {"n": k, "pct": round(100 * k / n, 1)}

    s = {
        "n_papers": n,
        "by_year": {y: sum(1 for r in rows if str(r["year"]) == y) for y in sorted({str(r["year"]) for r in rows})},
        "by_primary_theme": {t: sum(1 for r in rows if primary(r) == t) for t in sorted({primary(r) for r in rows})},
        "received_to_accepted_days": _quart([r["received_to_accepted_days"] for r in rows]),
        "received_to_accepted_days_by_primary_theme": {
            t: _quart([r["received_to_accepted_days"] for r in rows if primary(r) == t]) for t in sorted({primary(r) for r in rows})},
        "crossref_citations": _quart([r["crossref_citations"] for r in rows]),
        "citations_per_year": _quart([cand[r["pmid"]]["citations_per_year"] for r in rows]),
        "structured_abstract": {"n": sum(1 for r in rows if cand[r["pmid"]]["structured_abstract"]),
                                "pct": round(100 * sum(1 for r in rows if cand[r["pmid"]]["structured_abstract"]) / n, 1)},
        "ablation": {"yes": _share(rows, "ablation", "yes"), "no": _share(rows, "ablation", "no"),
                     "unclear": _share(rows, "ablation", "unclear")},
        "statistical_validation_not_determinable": _share(rows, "statistical_validation", "not determinable"),
        "computational_cost_reported_yes": _share(rows, "computational_cost_reported", "yes"),
        "code_or_data_availability_not_determinable": _share(rows, "code_or_data_availability", "not determinable"),
        "flag_external_validation": tally("external_validation", "external"),
        "flag_internal_validation_only": tally("external_validation", "internal"),
        "flag_validation_none_or_unclear": tally("external_validation", "none_or_unclear"),
        "flag_multi_dataset": tally("multi_dataset"),
        "flag_real_sample_application": tally("real_sample_application"),
        "flag_formal_stats_in_abstract": tally("formal_stats"),
        "flag_named_baselines": tally("named_baselines"),
        "flag_code_or_data_stated_in_abstract": tally("code_data_stated"),
        "flag_public_dataset_used": tally("public_data"),
        "flag_uncertainty_output": tally("uncertainty_output"),
        "top10_by_crossref_citations": [
            {k: r[k] for k in ("pmid", "doi", "year", "crossref_citations", "title")}
            for r in sorted(rows, key=lambda r: -(r["crossref_citations"] or 0))[:10]],
        "flags_top8_vs_rest": {},
    }
    top8 = {r["pmid"] for r in sorted(rows, key=lambda r: -(r["crossref_citations"] or 0))[:8]}
    for key, val in [("external_validation", "external"), ("multi_dataset", True), ("formal_stats", True),
                     ("named_baselines", True), ("code_data_stated", True), ("real_sample_application", True)]:
        a = [flags[r["pmid"]].get(key) == val for r in rows if r["pmid"] in top8]
        b = [flags[r["pmid"]].get(key) == val for r in rows if r["pmid"] not in top8]
        s["flags_top8_vs_rest"][key] = {"top8_pct": round(100 * sum(a) / len(a), 1),
                                        "rest_pct": round(100 * sum(b) / len(b), 1)}
    (RAW / "summary_stats.json").write_text(json.dumps(s, indent=1))
    print(json.dumps({k: v for k, v in s.items() if k != "top10_by_crossref_citations"}, indent=1))


STAGES = {"search": stage_search, "fetch": stage_fetch, "crossref": stage_crossref, "rank": stage_rank,
          "coverage": stage_coverage, "context": stage_context, "pmc": stage_pmc, "build": stage_build}

if __name__ == "__main__":
    args = sys.argv[1:] or ["all"]
    force = "--force" in args
    names = [a for a in args if not a.startswith("--")]
    for name in (list(STAGES) if names == ["all"] else names):
        STAGES[name](force=force)
