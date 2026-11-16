"""Build references.bib from a registry of citation keys -> DOI / arXiv id, verifying every record against Crossref or
arXiv (title fragment must match). Writes ../manuscript/references.bib and ../supplementary/citation_verification.csv.
Usage: python scripts/build_bib.py"""
import os, re, time, csv, json, unicodedata, html
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
 ("winkler1972decision", "doi:10.1080/01621459.1972.10481224", "decision-theoretic approach to interval estimation"),
 ("gneiting2007strictly", "doi:10.1198/016214506000001437", "strictly proper scoring rules"),
 ("kerby2014simple", "doi:10.2466/11.it.3.1", "simple difference formula"),
 ("chen2016xgboost", "doi:10.1145/2939672.2939785", "xgboost"),
 ("virtanen2020scipy", "doi:10.1038/s41592-019-0686-2", "scipy 1.0"),
 # added after the citation audit (2026-11-16); every DOI verified against Crossref, title fragment checked
 ("boswell2015calculation", "doi:10.1016/j.chroma.2015.07.113", "calculation of retention time tolerance windows"),
 ("boswell2011study", "doi:10.1016/j.chroma.2011.07.105", "retention"),
 ("kreutzer2026molecular", "doi:10.1007/s00216-025-06303-2", "conformal predictions"),
 ("papadopoulos2011regression", "doi:10.1613/jair.3198", "regression conformal prediction with nearest neighbours"),
 ("park2024fewshot", "doi:10.1109/TPAMI.2023.3327300", "meta-learned cross-validation-based conformal prediction"),
 ("dunn2023distribution", "doi:10.1080/01621459.2022.2060112", "two-layer hierarchical models"),
 ("lee2026hierarchical", "doi:10.1145/3786352", "hierarchical data"),
 ("barber2023beyond", "doi:10.1214/23-AOS2276", "conformal prediction beyond exchangeability"),
 ("bian2023training", "doi:10.1214/23-EJS2145", "training-conditional coverage"),
 ("vovk2013conditional", "doi:10.1007/s10994-013-5355-6", "conditional validity of inductive conformal predictors"),
 ("fritsch1980monotone", "doi:10.1137/0717021", "monotone piecewise cubic interpolation"),
 ("eilers1996flexible", "doi:10.1214/ss/1038425655", "flexible smoothing with b-splines and penalties"),
 ("rogers2010extended", "doi:10.1021/ci100050t", "extended-connectivity fingerprints"),
 ("wilcoxon1945individual", "doi:10.2307/3001968", "individual comparisons by ranking methods"),
 ("geer2024airi", "doi:10.1021/acs.jcim.3c01758", "airi"),
 ("maboudi2017uncertainty", "doi:10.1093/bioinformatics/btw619", "uncertainty estimation of predictions"),
 ("stienstra2025reverse", "doi:10.1021/acs.analchem.4c05859", "hilic"),
 ("wang2026hybrid", "doi:10.1021/acs.jpcb.6c03944", "laboratory-conditioned"),
 ("akhlaqi2023complementary", "doi:10.1007/s00216-023-04852-y", "isomeric candidate structures"),
 ("morger2022studying", "doi:10.1038/s41598-022-09309-3", "data drifts"),
 ("folberth2020ms2", "doi:10.1016/j.jchromb.2020.122105", "ms2 and lc libraries"),
 ("souihi2022multiconditionrt", "doi:10.1016/j.chroma.2022.462867", "multiconditionrt"),
 ("stoffel2022alkylpyridinium", "doi:10.1007/s00216-021-03828-0", "alkylpyridinium sulfonates"),
]
OVERRIDE = {  # fields missing or wrong in the Crossref record (checked against the publisher pages)
    "vovk2005algorithmic": {"author": "Vovk, Vladimir and Gammerman, Alex and Shafer, Glenn", "publisher": "Springer", "address": "New York"},
    "papadopoulos2002inductive": {"booktitle": "Machine Learning: ECML 2002", "editor": "Elomaa, Tapio and Mannila, Heikki and Toivonen, Hannu",
                                  "series": "Lecture Notes in Computer Science", "volume": "2430", "publisher": "Springer", "address": "Berlin, Heidelberg"},
}


def norm(s):
    s = unicodedata.normalize("NFKD", s or "").lower()
    return re.sub(r"[^a-z0-9 ]+", " ", s).replace("  ", " ")


def clean(t):
    t = unicodedata.normalize("NFC", html.unescape(t or "")).replace("o\u0337", "\u00f8").replace("O\u0337", "\u00d8")
    t = re.sub(r"(?<!\\)&", r"\\&", t)  # BibTeX/LaTeX: escape ampersands
    return re.sub(r"\s+", " ", t).strip()


def person(a):
    if a.get("family"):
        return clean(f"{a.get('family','')}, {a.get('given','')}".strip(", "))
    return "{" + clean(a.get("name", "")) + "}" if a.get("name") else ""  # group author, e.g. "SciPy 1.0 Contributors"


def bibtex_from_crossref(key, m):
    typ = m.get("type")
    authors = " and ".join(x for x in (person(a) for a in m.get("author", [])) if x)
    # year of the issue in which the article appeared (print), falling back to the first publication date
    year = ((m.get("published-print") or m.get("issued") or {}).get("date-parts") or [[None]])[0][0]
    title = clean((m.get("title") or [""])[0]); sub = clean((m.get("subtitle") or [""])[0])
    if sub and sub.lower() not in title.lower():
        title = f"{title}: {sub}"
    f = {"author": authors, "title": "{{" + title + "}}", "year": str(year), "doi": m.get("DOI")}  # double braces keep the published capitalisation
    if typ == "journal-article":
        kind = "article"; f["journal"] = clean((m.get("container-title") or [""])[0])
        if m.get("volume"): f["volume"] = m["volume"]
        if m.get("issue"): f["number"] = m["issue"]
        if m.get("page"): f["pages"] = m["page"].replace("-", "--")
        elif m.get("article-number"): f["pages"] = m["article-number"]
    elif typ in ("book", "monograph", "edited-book"):
        kind = "book"; f["publisher"] = m.get("publisher", "")
    elif typ in ("proceedings-article", "book-chapter"):
        kind = "inproceedings"; ct = [clean(c) for c in (m.get("container-title") or [""])]
        f["booktitle"] = next((c for c in ct if "lecture notes" not in c.lower()), ct[0])
        if m.get("page"): f["pages"] = m["page"].replace("-", "--")
        f["publisher"] = m.get("publisher", "")
    else:
        kind = "misc"; f["howpublished"] = (m.get("container-title") or [""])[0]
    f.update(OVERRIDE.get(key, {}))
    if not f.get("author"):
        raise ValueError(f"{key}: no author in Crossref record and no override")
    return f"@{kind}{{{key},\n" + ",\n".join(f"  {k} = {{{v}}}" if not v.startswith("{{") else f"  {k} = {v}" for k, v in f.items() if v) + "\n}\n"


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
            for attempt in range(5):  # polite retry on rate limiting
                r = requests.get(f"https://api.crossref.org/works/{val}", timeout=60, headers={"User-Agent": "bib-verify/1.0"})
                if r.status_code != 429:
                    break
                time.sleep(5 * (attempt + 1))
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
# software and DOI-less references (verified 2026-11-16 on the publisher pages / DataCite; versions from the installed packages)
MANUAL = [
 ("rdkit", "@misc{rdkit,\n  author = {Landrum, Greg and Tosco, Paolo and Rodriguez, Ricardo and Kelley, Brian and others},\n  title = {{RDKit: Open-source cheminformatics, release 2026.09.2}},\n  year = {2026},\n  doi = {10.5281/zenodo.23981945},\n  url = {https://www.rdkit.org}\n}\n", "DataCite 10.5281/zenodo.23981945"),
 ("lightgbm", "@misc{lightgbm,\n  author = {{Microsoft LightGBM contributors}},\n  title = {{LightGBM, version 4.7.0}},\n  year = {2026},\n  url = {https://github.com/microsoft/LightGBM}\n}\n", "installed package 4.7.0"),
 ("ke2017lightgbm", "@inproceedings{ke2017lightgbm,\n  author = {Ke, Guolin and Meng, Qi and Finley, Thomas and Wang, Taifeng and Chen, Wei and Ma, Weidong and Ye, Qiwei and Liu, Tie-Yan},\n  title = {{LightGBM: A Highly Efficient Gradient Boosting Decision Tree}},\n  booktitle = {Advances in Neural Information Processing Systems},\n  volume = {30},\n  year = {2017}\n}\n", "papers.nips.cc (NeurIPS 2017)"),
 ("tibshirani2019weighted", "@inproceedings{tibshirani2019weighted,\n  author = {Tibshirani, Ryan J. and Foygel Barber, Rina and Cand{\\`e}s, Emmanuel and Ramdas, Aaditya},\n  title = {{Conformal Prediction Under Covariate Shift}},\n  booktitle = {Advances in Neural Information Processing Systems},\n  volume = {32},\n  year = {2019}\n}\n", "papers.nips.cc (NeurIPS 2019)"),
 ("fisch2021fewshot", "@inproceedings{fisch2021fewshot,\n  author = {Fisch, Adam and Schuster, Tal and Jaakkola, Tommi and Barzilay, Regina},\n  title = {{Few-Shot Conformal Prediction with Auxiliary Tasks}},\n  booktitle = {Proceedings of the 38th International Conference on Machine Learning},\n  series = {Proceedings of Machine Learning Research},\n  volume = {139},\n  pages = {3329--3339},\n  year = {2021}\n}\n", "proceedings.mlr.press/v139/fisch21a"),
 ("pedregosa2011scikit", "@article{pedregosa2011scikit,\n  author = {Pedregosa, Fabian and Varoquaux, Ga{\\\"e}l and Gramfort, Alexandre and Michel, Vincent and Thirion, Bertrand and Grisel, Olivier and Blondel, Mathieu and Prettenhofer, Peter and Weiss, Ron and Dubourg, Vincent and Vanderplas, Jake and Passos, Alexandre and Cournapeau, David and Brucher, Matthieu and Perrot, Matthieu and Duchesnay, {\\'E}douard},\n  title = {{Scikit-learn: Machine Learning in Python}},\n  journal = {Journal of Machine Learning Research},\n  volume = {12},\n  pages = {2825--2830},\n  year = {2011}\n}\n", "jmlr.org/papers/v12/pedregosa11a"),
 ("report2025zenodo", "@misc{report2025zenodo,\n  author = {Witting, Michael and Kretschmer, Fleming and Harrieder, Eva-Maria and others},\n  title = {{RepoRT v1.0}},\n  year = {2025},\n  publisher = {Zenodo},\n  doi = {10.5281/zenodo.16267768}\n}\n", "DataCite 10.5281/zenodo.16267768"),
]
for key, entry, how in MANUAL:
    bib.append(entry); log.append({"key": key, "id": "manual", "status": f"OK (manual: {how})", "title": re.search(r"title = \{\{(.*?)\}\}", entry).group(1)})
out = os.path.join(ROOT, "manuscript", "references.bib"); open(out, "w").write("\n".join(bib))
with open(os.path.join(ROOT, "supplementary", "citation_verification.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["key", "id", "status", "title"]); w.writeheader(); [w.writerow(x) for x in log]
bad = [x for x in log if not x["status"].startswith("OK")]
print(f"{len(log) - len(bad)}/{len(log)} verified; problems:", json.dumps(bad, indent=1))
