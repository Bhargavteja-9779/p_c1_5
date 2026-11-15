"""Received->accepted intervals restricted to topic-relevant articles (ML/chemometrics/data analysis)
in candidate journals, published 2024-01-01..2026-10-31, from PubMed history dates."""
import sys, json, time, statistics, csv
sys.path.insert(0, "scripts")
from pubmed_timelines import get, parse_date, EUTILS, EXCL
import xml.etree.ElementTree as ET
J0 = {"Analytical Chemistry (ACS)": '"Anal Chem"[ta]', "Analytica Chimica Acta (Elsevier)": '"Anal Chim Acta"[ta]',
     "Talanta (Elsevier)": '"Talanta"[ta]', "J Chem Inf Model (ACS)": '"J Chem Inf Model"[ta]',
     "Chemical Science (RSC)": '"Chem Sci"[ta]', "Water Research (Elsevier)": '"Water Res"[ta]',
     "J Hazardous Materials (Elsevier)": '"J Hazard Mater"[ta]'}
J = {k: v for k, v in J0.items() if k.split(" (")[0] in ("Analytical Chemistry", "Analytica Chimica Acta", "Talanta", "J Chem Inf Model")}
TOPIC = ('("machine learning"[tiab] OR "deep learning"[tiab] OR "neural network"[tiab] OR chemometric*[tiab] '
         'OR "partial least squares"[tiab] OR "retention time"[tiab] OR "uncertainty"[tiab] OR "prediction interval*"[tiab])')
MIN, MAX = (sys.argv[2], sys.argv[3]) if len(sys.argv) > 3 else ("2024/01/01", "2026/10/31")
rows = []
for name, jq in J.items():
    term = f"{jq} AND {TOPIC} AND journal article[pt]{EXCL}"
    r = get(EUTILS + "esearch.fcgi", {"db": "pubmed", "term": term, "retmax": 600, "retmode": "json", "datetype": "pdat",
                                      "mindate": MIN, "maxdate": MAX, "sort": "pub_date"})
    ids = r.json()["esearchresult"]["idlist"]; total = int(r.json()["esearchresult"]["count"])
    days = []
    for i in range(0, len(ids), 200):
        time.sleep(0.4)
        root = ET.fromstring(get(EUTILS + "efetch.fcgi", {"db": "pubmed", "id": ",".join(ids[i:i+200]), "retmode": "xml"}).content)
        for art in root.findall(".//PubmedArticle"):
            h = {}
            for p in art.findall(".//History/PubMedPubDate"):
                h.setdefault(p.get("PubStatus"), parse_date(p))
            if h.get("received") and h.get("accepted"):
                d = (h["accepted"] - h["received"]).days
                if 0 <= d < 2000: days.append(d)
    q = statistics.quantiles(days, n=4)
    row = {"journal": name, "topic_total": total, "sampled": len(ids), "n_with_dates": len(days), "median": statistics.median(days),
           "q1": q[0], "q3": q[2], "p90": statistics.quantiles(days, n=10)[8], "pct_le_60": round(100*sum(d<=60 for d in days)/len(days),1),
           "pct_le_100": round(100*sum(d<=100 for d in days)/len(days),1)}
    rows.append(row); print(json.dumps(row), flush=True)
with open(sys.argv[1], "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); [w.writerow(r) for r in rows]
