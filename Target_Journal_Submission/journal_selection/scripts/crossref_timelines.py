"""Compute received->accepted intervals from Crossref 'assertion' metadata (deposited by Springer
Nature, Taylor & Francis, IOP and some others) for recent journal articles.
Usage: python crossref_timelines.py out.csv [from-pub-date] [until-pub-date]
"""
import sys, time, json, csv, statistics, datetime as dt, requests

JOURNALS = [
    ("Journal of Big Data", "2196-1115", "Springer Nature"),
    ("Complex & Intelligent Systems", "2198-6053", "Springer Nature"),
    ("Artificial Intelligence Review", "1573-7462", "Springer Nature"),
    ("Machine Learning", "1573-0565", "Springer Nature"),
    ("Data Mining and Knowledge Discovery", "1573-756X", "Springer Nature"),
    ("Neural Computing and Applications", "1433-3058", "Springer Nature"),
    ("Journal of Intelligent Manufacturing", "1572-8145", "Springer Nature"),
    ("Cognitive Computation", "1866-9964", "Springer Nature"),
    ("Applied Intelligence", "1573-7497", "Springer Nature"),
    ("Int. J. Machine Learning and Cybernetics", "1868-808X", "Springer Nature"),
    ("Machine Intelligence Research", "2731-5398", "Springer Nature"),
    ("J. King Saud Univ. Comput. Inf. Sci.", "2213-1248", "Springer Nature"),
    ("EPJ Data Science", "2193-1127", "Springer Nature"),
    ("Journal of Cheminformatics", "1758-2946", "Springer Nature"),
    ("Engineering with Computers", "1435-5663", "Springer Nature"),
    ("Cluster Computing", "1573-7543", "Springer Nature"),
    ("Statistics and Computing", "1573-1375", "Springer Nature"),
    ("Journal of Cloud Computing", "2192-113X", "Springer Nature"),
    ("Cybersecurity", "2523-3246", "Springer Nature"),
    ("Energy Informatics", "2520-8942", "Springer Nature"),
    ("Journal of Materials Science & Technology", "1005-0302", "Elsevier"),
    ("Machine Learning: Science and Technology", "2632-2153", "IOP"),
    ("Int. J. Production Research", "1366-588X", "Taylor & Francis"),
    ("Int. J. Geographical Information Science", "1362-3087", "Taylor & Francis"),
    ("GIScience & Remote Sensing", "1943-7226", "Taylor & Francis"),
    ("Int. J. Digital Earth", "1753-8955", "Taylor & Francis"),
    ("Enterprise Information Systems", "1751-7583", "Taylor & Francis"),
    ("Connection Science", "1360-0494", "Taylor & Francis"),
    ("Technometrics", "1537-2723", "Taylor & Francis"),
    ("J. Computational and Graphical Statistics", "1537-2715", "Taylor & Francis"),
    ("Virtual and Physical Prototyping", "1745-2767", "Taylor & Francis"),
    ("Journal of Computational Design and Engineering", "2288-5048", "OUP"),
    ("Data-Centric Engineering", "2632-6736", "Cambridge UP"),
    ("Structural Health Monitoring", "1741-3168", "SAGE"),
    ("IET Renewable Power Generation", "1752-1424", "IET/Wiley"),
    ("Advanced Intelligent Systems", "2640-4567", "Wiley"),
    ("International Journal of Intelligent Systems", "1098-111X", "Wiley"),
]


def parse(v):
    for fmt in ("%d %B %Y", "%d %b %Y", "%Y-%m-%d", "%B %d, %Y", "%d/%m/%Y"):
        try:
            return dt.datetime.strptime(v.strip(), fmt).date()
        except Exception:
            pass
    return None


def stats(issn, frm, until, rows=500):
    url = f"https://api.crossref.org/journals/{issn}/works"
    params = {"filter": f"from-pub-date:{frm},until-pub-date:{until},type:journal-article", "rows": rows,
              "select": "DOI,assertion,title,type"}
    for k in range(4):
        try:
            r = requests.get(url, params=params, timeout=90, headers={"User-Agent": "journal-timeline-study/0.1"})
            if r.status_code == 200:
                break
        except requests.RequestException:
            pass
        time.sleep(3 * (k + 1))
    else:
        return {"error": "request failed"}
    msg = r.json()["message"]; items = msg["items"]
    days, have_assert = [], 0
    for it in items:
        a = {x.get("name", "").lower(): x.get("value") for x in it.get("assertion", []) if x.get("value")}
        if a: have_assert += 1
        rec = a.get("received"); acc = a.get("accepted")
        if rec and acc:
            d1, d2 = parse(rec), parse(acc)
            if d1 and d2 and 0 <= (d2 - d1).days < 2000:
                days.append((d2 - d1).days)
    out = {"total_in_window": msg.get("total-results"), "sampled": len(items), "with_assertions": have_assert,
           "n_with_dates": len(days)}
    if len(days) >= 5:
        q = statistics.quantiles(days, n=4)
        out.update({"median": statistics.median(days), "q1": q[0], "q3": q[2], "mean": round(statistics.mean(days), 1),
                    "pct_le_30": round(100 * sum(d <= 30 for d in days) / len(days), 1),
                    "pct_le_60": round(100 * sum(d <= 60 for d in days) / len(days), 1),
                    "pct_le_100": round(100 * sum(d <= 100 for d in days) / len(days), 1)})
    return out


if __name__ == "__main__":
    out_csv = sys.argv[1]; frm = sys.argv[2] if len(sys.argv) > 2 else "2025-11-01"; until = sys.argv[3] if len(sys.argv) > 3 else "2026-10-31"
    res = []
    for label, issn, pub in JOURNALS:
        s = stats(issn, frm, until); s.update({"journal": label, "issn": issn, "publisher": pub, "window": f"{frm}..{until}"})
        res.append(s); print(json.dumps(s), flush=True); time.sleep(1)
    keys = ["journal", "publisher", "issn", "window", "total_in_window", "sampled", "with_assertions", "n_with_dates", "median", "q1", "q3", "mean", "pct_le_30", "pct_le_60", "pct_le_100", "error"]
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore"); w.writeheader(); [w.writerow(r) for r in res]
