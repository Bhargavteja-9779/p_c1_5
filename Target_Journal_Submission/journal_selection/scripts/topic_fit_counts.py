"""Count recent (2024-2026) PubMed records per candidate journal matching topic queries (scope-fit evidence)."""
import requests, time, json, sys
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
J = {"Analytical Chemistry": '"Anal Chem"[ta]', "Analytica Chimica Acta": '"Anal Chim Acta"[ta]', "Talanta": '"Talanta"[ta]',
     "J Chem Inf Model": '"J Chem Inf Model"[ta]', "Chemical Science": '"Chem Sci"[ta]', "Water Research": '"Water Res"[ta]',
     "J Hazardous Materials": '"J Hazard Mater"[ta]', "Environmental Research": '"Environ Res"[ta]',
     "J Chem Theory Comput": '"J Chem Theory Comput"[ta]', "Spectrochim Acta A": '"Spectrochim Acta A Mol Biomol Spectrosc"[ta]'}
T = {"machine learning (tiab)": '("machine learning"[tiab] OR "deep learning"[tiab] OR "neural network"[tiab])',
     "chemometric* (tiab)": 'chemometric*[tiab]',
     "retention time prediction": '("retention time"[tiab] AND predict*[tiab])',
     "non-target screening": '("non-target"[tiab] OR "nontarget"[tiab] OR "suspect screening"[tiab])',
     "calibration transfer": '("calibration transfer"[tiab] OR "model transfer"[tiab] OR "domain adaptation"[tiab])',
     "uncertainty/conformal": '("conformal prediction"[tiab] OR "prediction interval*"[tiab] OR "uncertainty quantification"[tiab])'}
out = {}
for jn, jq in J.items():
    out[jn] = {}
    for tn, tq in T.items():
        q = f"{jq} AND {tq} AND 2024:2026[dp]"
        for k in range(4):
            try:
                r = requests.get(E, params={"db": "pubmed", "term": q, "retmode": "json", "retmax": 0}, timeout=60); break
            except Exception: time.sleep(2)
        out[jn][tn] = int(r.json()["esearchresult"]["count"]); time.sleep(0.35)
    tot = requests.get(E, params={"db": "pubmed", "term": f"{jq} AND 2024:2026[dp]", "retmode": "json", "retmax": 0}, timeout=60).json()["esearchresult"]["count"]
    out[jn]["ALL 2024-2026"] = int(tot); time.sleep(0.35)
    print(jn, json.dumps(out[jn]), flush=True)
json.dump(out, open(sys.argv[1], "w"), indent=1)
