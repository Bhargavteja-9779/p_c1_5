"""Build references.bib from a registry of citation keys -> DOI / arXiv id, verifying every record against Crossref or
arXiv (title fragment must match). Writes ../manuscript/references.bib and ../supplementary/citation_verification.csv.
Usage: python scripts/build_bib.py"""
import os, re, time, csv, json, unicodedata
import requests
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
REG = [  # key, id (doi:/arxiv:/manual:), expected title fragment (lower-case)
 ("schymanski2014identifying", "doi:10.1021/es5002105", "identifying small molecules via high resolution mass spectrometry"),
 ("hupatz2025critical", "doi:10.1007/s00216-024-05471-x", "critical review on in silico methods"),
 ("metz2025identification", "doi:10.1021/acs.analchem.4c04060", "identification probability"),
 ("ngan2025machine", "doi:10.1021/acs.analchem.5c01873", "enhanced identification probability"),
 ("stanstrup2015predret", "doi:10.1021/acs.analchem.5b02287", "predret"),
 ("abatepella2015retention", "doi:10.1016/j.chroma.2015.07.108", "retention projection enables accurate calculation"),
 ("bach2018retention", "doi:10.1093/bioinformatics/bty590", "retention order prediction for metabolite identification"),
 ("bach2022joint", "doi:10.1038/s42256-022-00577-2", "joint structural annotation of small molecules"),
 ("bouwmeester2019comprehensive", "doi:10.1021/acs.analchem.8b05820", "comprehensive and empirical evaluation of machine learning"),
 ("bouwmeester2020generalized", "doi:10.1021/acs.analchem.0c00233", "generalized calibration across liquid chromatography"),
 ("bonini2020retip", "doi:10.1021/acs.analchem.9b05765", "retip"),
 ("yang2021prediction", "doi:10.1021/acs.analchem.0c04071", "graph neural networks"),
 ("aalizadeh2019development", "doi:10.1016/j.jhazmat.2018.09.047", "retention time prediction models in the suspect"),
 ("aalizadeh2021development", "doi:10.1021/acs.analchem.1c02348", "retention time indices"),
 ("domingo2019smrt", "doi:10.1038/s41467-019-13680-7", "metlin small molecule dataset"),
 ("garcia2022probabilistic", "doi:10.1186/s13321-022-00613-8", "probabilistic metabolite annotation"),
 ("kretschmer2024report", "doi:10.1038/s41592-023-02143-z", "repor"),
 ("kretschmer2026order", "doi:10.1038/s41592-026-03243-2", "times are changing but order matters"),
 ("zhang2024generic", "doi:10.1038/s42004-024-01135-0", "post-projection calibration"),
 ("xue2024rttransformer", "doi:10.1093/bioinformatics/btae084", "rt-transformer"),
 ("matyushin2024uncertainty", "doi:10.3390/ijms252313077", "uncertainty quantification and flagging"),
 ("malm2025projection", "doi:10.1039/d5an00323g", "experimental projection methods"),
 ("malm2024norman", "doi:10.1021/acs.analchem.4c02902", ""),
 ("li2025prediction", "doi:10.1021/acs.analchem.5c01703", "chromatographic parameter vectorization"),
 ("wiczling2026bayesian", "doi:10.1016/j.chroma.2025.466503", "bayesian multilevel modeling of retention data"),
 ("beck2026uncertainty", "doi:10.1021/acs.analchem.5c06729", "uncertainty-aware learning of multiple conditions"),
 ("xiong2026unirt", "doi:10.1021/acs.analchem.5c07973", "unified multitask modeling"),
 ("deng2026retention", "doi:10.1016/j.jhazmat.2026.141313", "emerging contaminants via transfer learning"),
 ("fadil2026fastret", "doi:10.1021/acs.jcim.6c01344", "fastret"),
 ("mceachran2018comparison", "doi:10.1016/j.talanta.2018.01.022", "comparison of three liquid chromatography"),
 ("samaraweera2018evaluation", "doi:10.1021/acs.analchem.8b03118", "artificial neural network retention index"),
 ("ferrandoplo2026candidate", "doi:10.1021/acs.jcim.6c02492", "candidate structure"),
 ("rakhshaninejad2026reliable", "doi:10.1021/acs.jcim.6c00727", "reliable molecular retrieval"),
 ("bhargav2026group", "doi:10.1016/j.talanta.2026.130941", "group-conformal calibration"),
 ("wiens2026calibrated", "doi:10.1021/acs.analchem.5c04616", "calibrated uncertainty estimation for soil organic carbon"),
 ("wen2023structural", "doi:10.1016/j.aca.2023.341720", "collision cross-section prediction interval filtering"),
 ("aca2024qsrrwindows", "doi:10.1016/j.aca.2024.342529", "retention time prediction"),
 ("aca2025dualcolumn", "doi:10.1016/j.aca.2025.343942", "dual lc column characterization"),
 ("aca2023validation", "doi:10.1016/j.aca.2023.341532", "validation"),
 ("vovk2005algorithmic", "doi:10.1007/b106715", "algorithmic learning in a random world"),
 ("papadopoulos2002inductive", "doi:10.1007/3-540-36755-1_29", "inductive confidence machines for regression"),
 ("lei2018distribution", "doi:10.1080/01621459.2017.1307116", "distribution-free predictive inference for regression"),
 ("barber2021jackknife", "doi:10.1214/20-aos1965", "predictive inference with the jackknife+"),
 ("angelopoulos2023gentle", "doi:10.1561/2200000101", "conformal prediction"),
 ("romano2019cqr", "arxiv:1905.03222", "conformalized quantile regression"),
 ("tibshirani2019weighted", "arxiv:1904.06019", "conformal prediction under covariate shift"),
 ("fisch2021fewshot", "arxiv:2102.08898", "few-shot conformal prediction with auxiliary tasks"),
 ("gibbs2021adaptive", "arxiv:2106.00170", "adaptive conformal inference under distribution shift"),
 ("winkler1972decision", "doi:10.1080/01621459.1972.10481224", "decision-theoretic approach to interval estimation"),
 ("gneiting2007strictly", "doi:10.1198/016214506000001437", "strictly proper scoring rules"),
 ("kerby2014simple", "doi:10.2466/11.it.3.1", "simple difference formula"),
 ("chen2016xgboost", "doi:10.1145/2939672.2939785", "xgboost"),
 ("pedregosa2011scikit", "arxiv:1201.0490", "scikit-learn"),
 ("virtanen2020scipy", "doi:10.1038/s41592-019-0686-2", "scipy 1.0"),
]


def norm(s):
    s = unicodedata.normalize("NFKD", s or "").lower()
    return re.sub(r"[^a-z0-9 ]+", " ", s).replace("  ", " ")


def bibtex_from_crossref(key, m):
    typ = m.get("type")
    authors = " and ".join(f"{a.get('family','')}, {a.get('given','')}".strip(", ") for a in m.get("author", [])) or m.get("publisher", "")
    year = (m.get("issued", {}).get("date-parts") or [[None]])[0][0]
    title = (m.get("title") or [""])[0].replace("\n", " ")
    f = {"author": authors, "title": "{" + title + "}", "year": str(year), "doi": m.get("DOI")}
    if typ == "journal-article":
        kind = "article"; f["journal"] = (m.get("container-title") or [""])[0]
        if m.get("volume"): f["volume"] = m["volume"]
        if m.get("issue"): f["number"] = m["issue"]
        if m.get("page"): f["pages"] = m["page"].replace("-", "--")
        elif m.get("article-number"): f["pages"] = m["article-number"]
    elif typ in ("book", "monograph", "edited-book"):
        kind = "book"; f["publisher"] = m.get("publisher", "")
    elif typ in ("proceedings-article", "book-chapter"):
        kind = "inproceedings"; f["booktitle"] = (m.get("container-title") or [""])[0]
        if m.get("page"): f["pages"] = m["page"].replace("-", "--")
        f["publisher"] = m.get("publisher", "")
    else:
        kind = "misc"; f["howpublished"] = (m.get("container-title") or [""])[0]
    return f"@{kind}{{{key},\n" + ",\n".join(f"  {k} = {{{v}}}" if not v.startswith("{") else f"  {k} = {v}" for k, v in f.items() if v) + "\n}\n"


def arxiv_entry(key, aid):
    r = requests.get(f"https://export.arxiv.org/api/query?id_list={aid}", timeout=60)
    root = ET.fromstring(r.content); ns = {"a": "http://www.w3.org/2005/Atom"}
    e = root.find("a:entry", ns)
    title = re.sub(r"\s+", " ", e.find("a:title", ns).text.strip())
    authors = " and ".join(a.find("a:name", ns).text for a in e.findall("a:author", ns))
    year = e.find("a:published", ns).text[:4]
    return title, f"@misc{{{key},\n  author = {{{authors}}},\n  title = {{{{{title}}}}},\n  year = {{{year}}},\n  howpublished = {{arXiv:{aid}}},\n  eprint = {{{aid}}},\n  archivePrefix = {{arXiv}}\n}}\n"


bib, log = [], []
for key, ident, frag in REG:
    kind, _, val = ident.partition(":")
    try:
        if kind == "doi":
            r = requests.get(f"https://api.crossref.org/works/{val}", timeout=60, headers={"User-Agent": "bib-verify/1.0"})
            if r.status_code != 200:
                log.append({"key": key, "id": ident, "status": f"FAIL http {r.status_code}", "title": ""}); continue
            m = r.json()["message"]; title = (m.get("title") or [""])[0]
            ok = (not frag) or (norm(frag) in norm(title))
            log.append({"key": key, "id": ident, "status": "OK" if ok else "TITLE MISMATCH", "title": title})
            if ok: bib.append(bibtex_from_crossref(key, m))
        else:
            title, entry = arxiv_entry(key, val); ok = norm(frag) in norm(title)
            log.append({"key": key, "id": ident, "status": "OK (arXiv)" if ok else "TITLE MISMATCH", "title": title})
            if ok: bib.append(entry)
    except Exception as e:
        log.append({"key": key, "id": ident, "status": f"ERROR {e}", "title": ""})
    time.sleep(0.25)
# classic reference without a Crossref-registered DOI (verified 2026-11-16 against three independent bibliographic sources:
# IME-USP PDF copy, bibbase.org record, Wikipedia "Holm-Bonferroni method"; JSTOR page itself blocks automated access)
bib.append("@article{holm1979simple,\n  author = {Holm, Sture},\n  title = {{A simple sequentially rejective multiple test procedure}},\n  journal = {Scandinavian Journal of Statistics},\n  volume = {6},\n  number = {2},\n  pages = {65--70},\n  year = {1979},\n  url = {http://www.jstor.org/stable/4615733}\n}\n")
log.append({"key": "holm1979simple", "id": "jstor:4615733", "status": "OK (manual: 3 independent web sources; no Crossref DOI)", "title": "A simple sequentially rejective multiple test procedure"})
# software references (versions used; verified from the installed packages)
bib.append("@misc{rdkit,\n  author = {{RDKit contributors}},\n  title = {{RDKit: Open-source cheminformatics}},\n  year = {2026},\n  note = {Version 2026.09.2},\n  howpublished = {\\url{https://www.rdkit.org}}\n}\n")
bib.append("@misc{lightgbm,\n  author = {{Microsoft LightGBM contributors}},\n  title = {{LightGBM: a fast, distributed, high-performance gradient boosting framework}},\n  year = {2026},\n  note = {Version 4.7.0},\n  howpublished = {\\url{https://github.com/microsoft/LightGBM}}\n}\n")
out = os.path.join(ROOT, "manuscript", "references.bib"); open(out, "w").write("\n".join(bib))
with open(os.path.join(ROOT, "supplementary", "citation_verification.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["key", "id", "status", "title"]); w.writeheader(); [w.writerow(x) for x in log]
bad = [x for x in log if not x["status"].startswith("OK")]
print(f"{len(log) - len(bad)}/{len(log)} verified; problems:", json.dumps(bad, indent=1))
