"""Compute submission(received)->acceptance intervals from PubMed article-history metadata.

For each journal, retrieves PubMed records of type 'Journal Article' (excluding reviews,
editorials, comments, letters) with publication date in a window, fetches XML, and parses
<History><PubMedPubDate PubStatus="received"|"accepted">. Reports N, median, IQR, mean,
and the share of articles accepted within 30/60/100 days.

Usage: python pubmed_timelines.py out.csv  (journal list defined below)
"""
import sys, time, json, statistics, datetime as dt
import xml.etree.ElementTree as ET
import requests

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
JOURNALS = [
    # (label, pubmed query fragment, publisher)
    ("Briefings in Bioinformatics", '"Brief Bioinform"[ta]', "OUP"),
    ("Bioinformatics", '"1367-4811"[is]', "OUP"),
    ("JAMIA", '"J Am Med Inform Assoc"[ta]', "OUP"),
    ("GigaScience", '"Gigascience"[ta]', "OUP"),
    ("NAR Genomics and Bioinformatics", '"NAR Genom Bioinform"[ta]', "OUP"),
    ("Bioinformatics Advances", '"Bioinform Adv"[ta]', "OUP"),
    ("European Heart Journal - Digital Health", '"Eur Heart J Digit Health"[ta]', "OUP"),
    ("PNAS Nexus", '"PNAS Nexus"[ta]', "OUP"),
    ("Computers in Biology and Medicine", '"Comput Biol Med"[ta]', "Elsevier"),
    ("Artificial Intelligence in Medicine", '"Artif Intell Med"[ta]', "Elsevier"),
    ("Computer Methods and Programs in Biomedicine", '"Comput Methods Programs Biomed"[ta]', "Elsevier"),
    ("Journal of Biomedical Informatics", '"J Biomed Inform"[ta]', "Elsevier"),
    ("Neural Networks", '"Neural Netw"[ta]', "Elsevier"),
    ("Medical Image Analysis", '"Med Image Anal"[ta]', "Elsevier"),
    ("Computational and Structural Biotechnology Journal", '"Comput Struct Biotechnol J"[ta]', "Elsevier"),
    ("International Journal of Medical Informatics", '"Int J Med Inform"[ta]', "Elsevier"),
    ("Patterns", '"Patterns (N Y)"[ta]', "Cell Press/Elsevier"),
    ("Journal of Cheminformatics", '"J Cheminform"[ta]', "Springer Nature"),
    ("Journal of Chemical Information and Modeling", '"J Chem Inf Model"[ta]', "ACS"),
    ("IEEE Journal of Biomedical and Health Informatics", '"IEEE J Biomed Health Inform"[ta]', "IEEE"),
    ("IEEE Transactions on Neural Networks and Learning Systems", '"IEEE Trans Neural Netw Learn Syst"[ta]', "IEEE"),
    ("IEEE Transactions on Biomedical Engineering", '"IEEE Trans Biomed Eng"[ta]', "IEEE"),
    ("npj Digital Medicine", '"NPJ Digit Med"[ta]', "Springer Nature"),
    ("Scientific Reports", '"Sci Rep"[ta]', "Springer Nature"),
    ("Journal of Neural Engineering", '"J Neural Eng"[ta]', "IOP"),
    ("Physiological Measurement", '"Physiol Meas"[ta]', "IOP"),
    ("Biomedical Signal Processing and Control", '"Biomed Signal Process Control"[ta]', "Elsevier"),
]
JOURNALS2 = [
    ("Chemical Science", '"Chem Sci"[ta]', "RSC"),
    ("Digital Discovery", '"Digit Discov"[ta]', "RSC"),
    ("Journal of Medicinal Chemistry", '"J Med Chem"[ta]', "ACS"),
    ("Analytical Chemistry", '"Anal Chem"[ta]', "ACS"),
    ("Environmental Science & Technology", '"Environ Sci Technol"[ta]', "ACS"),
    ("Journal of Chemical Theory and Computation", '"J Chem Theory Comput"[ta]', "ACS"),
    ("Molecular Pharmaceutics", '"Mol Pharm"[ta]', "ACS"),
    ("Chemical Research in Toxicology", '"Chem Res Toxicol"[ta]', "ACS"),
    ("Journal of Proteome Research", '"J Proteome Res"[ta]', "ACS"),
    ("ACS Central Science", '"ACS Cent Sci"[ta]', "ACS"),
    ("Journal of Physical Chemistry Letters", '"J Phys Chem Lett"[ta]', "ACS"),
    ("Journal of Hazardous Materials", '"J Hazard Mater"[ta]', "Elsevier"),
    ("Water Research", '"Water Res"[ta]', "Elsevier"),
    ("Science of the Total Environment", '"Sci Total Environ"[ta]', "Elsevier"),
    ("Environment International", '"Environ Int"[ta]', "Elsevier"),
    ("Food Chemistry", '"Food Chem"[ta]', "Elsevier"),
    ("Analytica Chimica Acta", '"Anal Chim Acta"[ta]', "Elsevier"),
    ("Talanta", '"Talanta"[ta]', "Elsevier"),
    ("Spectrochimica Acta Part A", '"Spectrochim Acta A Mol Biomol Spectrosc"[ta]', "Elsevier"),
    ("European Journal of Medicinal Chemistry", '"Eur J Med Chem"[ta]', "Elsevier"),
    ("Environmental Research", '"Environ Res"[ta]', "Elsevier"),
    ("Journal of Environmental Management", '"J Environ Manage"[ta]', "Elsevier"),
    ("Advanced Science", '"Adv Sci (Weinh)"[ta]', "Wiley"),
    ("Molecular Informatics", '"Mol Inform"[ta]', "Wiley"),
    ("Journal of Computational Chemistry", '"J Comput Chem"[ta]', "Wiley"),
    ("Archives of Toxicology", '"Arch Toxicol"[ta]', "Springer Nature"),
    ("Green Chemistry", '"Green Chem"[ta]', "RSC"),
    ("Computational Biology and Chemistry", '"Comput Biol Chem"[ta]', "Elsevier"),
    ("Methods", '"Methods"[ta]', "Elsevier"),
    ("Journal of Molecular Biology", '"J Mol Biol"[ta]', "Elsevier"),
    ("Protein Science", '"Protein Sci"[ta]', "Wiley"),
    ("Journal of Translational Medicine", '"J Transl Med"[ta]', "Springer Nature"),
    ("Molecular Systems Biology", '"Mol Syst Biol"[ta]', "EMBO/Springer Nature"),
    ("Cell Reports Methods", '"Cell Rep Methods"[ta]', "Cell Press/Elsevier"),
]
EXCL = ' NOT (review[pt] OR editorial[pt] OR comment[pt] OR letter[pt] OR erratum[pt] OR "published erratum"[pt] OR "retraction of publication"[pt])'


def get(url, params, retries=4):
    for k in range(retries):
        try:
            r = requests.get(url, params=params, timeout=60)
            if r.status_code == 200:
                return r
        except requests.RequestException:
            pass
        time.sleep(2 * (k + 1))
    raise RuntimeError(f"failed {url} {params}")


def parse_date(node):
    try:
        y = int(node.findtext("Year")); m = int(node.findtext("Month")); d = int(node.findtext("Day"))
        return dt.date(y, m, d)
    except Exception:
        return None


def journal_stats(query, mindate, maxdate, cap=400):
    term = f"{query} AND journal article[pt]{EXCL}"
    r = get(EUTILS + "esearch.fcgi", {"db": "pubmed", "term": term, "retmax": cap, "retmode": "json",
                                       "datetype": "pdat", "mindate": mindate, "maxdate": maxdate, "sort": "pub_date"})
    js = r.json()["esearchresult"]
    ids = js.get("idlist", []); total = int(js.get("count", 0))
    days, rev = [], 0
    for i in range(0, len(ids), 200):
        time.sleep(0.4)
        x = get(EUTILS + "efetch.fcgi", {"db": "pubmed", "id": ",".join(ids[i:i + 200]), "retmode": "xml"})
        root = ET.fromstring(x.content)
        for art in root.findall(".//PubmedArticle"):
            hist = {}
            for p in art.findall(".//History/PubMedPubDate"):
                hist.setdefault(p.get("PubStatus"), parse_date(p))
            if hist.get("received") and hist.get("accepted"):
                dd = (hist["accepted"] - hist["received"]).days
                if 0 <= dd < 2000:
                    days.append(dd)
                    rev += 1 if hist.get("revised") else 0
    out = {"total_in_window": total, "sampled": len(ids), "n_with_dates": len(days)}
    if days:
        q = statistics.quantiles(days, n=4)
        out.update({"median": statistics.median(days), "q1": q[0], "q3": q[2], "mean": round(statistics.mean(days), 1),
                    "pct_le_30": round(100 * sum(d <= 30 for d in days) / len(days), 1),
                    "pct_le_60": round(100 * sum(d <= 60 for d in days) / len(days), 1),
                    "pct_le_100": round(100 * sum(d <= 100 for d in days) / len(days), 1)})
    return out


if __name__ == "__main__":
    out_csv = sys.argv[1]
    mindate, maxdate = (sys.argv[2], sys.argv[3]) if len(sys.argv) > 3 else ("2025/11/01", "2026/10/31")
    rows = []
    for label, q, pub in (JOURNALS2 if len(sys.argv) > 4 and sys.argv[4] == "batch2" else JOURNALS):
        time.sleep(0.5)
        try:
            s = journal_stats(q, mindate, maxdate)
        except Exception as e:
            s = {"error": str(e)}
        s.update({"journal": label, "publisher": pub, "query": q, "window": f"{mindate}-{maxdate}"})
        rows.append(s); print(json.dumps(s), flush=True)
    import csv
    keys = ["journal", "publisher", "window", "total_in_window", "sampled", "n_with_dates", "median", "q1", "q3",
            "mean", "pct_le_30", "pct_le_60", "pct_le_100", "query", "error"]
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore"); w.writeheader(); [w.writerow(r) for r in rows]
