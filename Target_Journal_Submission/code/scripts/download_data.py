"""Download the two public data sets and verify their SHA-256 checksums:
  (1) RepoRT v1.0 (Zenodo record 16267768; CC BY-SA 4.0) -> data/raw/repoRT/
  (2) NORMAN interlaboratory NTS retention-time data: Electronic Supplementary Information of Malm & Kruve, Analyst 150 (2025)
      (doi 10.1039/D5AN00323G; open access, PMC12268317), via the Europe PMC supplementary-files service -> data/raw/norman_ilc/
Usage: python scripts/download_data.py --dest data/raw [--no-norman]"""
import argparse, hashlib, os, sys, zipfile, urllib.request
URL = "https://zenodo.org/api/records/16267768/files/f-kretschmer/RepoRT-v1.0.zip/content"
SHA256 = "4eeeb46a65ec4950eca67f652caed5cffbdf5bcb5dffc45e496edae0552c6091"
NORMAN_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12268317/supplementaryFiles"
# Europe PMC assembles the supplementary zip on request, so the archive bytes (and their hash) change between downloads while
# its contents do not; the checksums are therefore verified for the extracted data files that the analysis reads.
NORMAN_FILES = {"AN-150-D5AN00323G-s002.xlsx": "f5e1e021bff5dec7ae8ad1139c5ef6847b9f94b558ec089f3bb8b31f16936afa",
                "SI/data/NORMAN_interlab_RT_data.csv": "7b146dc202892973b2708e42db8f4073f829617482d991cf1c37d121a8ab420e",
                "SI/data/gradient_lenght_codes.csv": "ce4475ce9009bd9a0bee58c96537ef96a108a487edf6b0e03f4d009925e77bc1"}
ap = argparse.ArgumentParser(); ap.add_argument("--dest", default="data/raw"); ap.add_argument("--no-norman", action="store_true")
a = ap.parse_args()


def fetch(url, path):
    print("downloading", url)
    req = urllib.request.Request(url, headers={"User-Agent": "conformrt-reproduction/1.0"})
    with urllib.request.urlopen(req) as r, open(path, "wb") as f:
        while True:
            b = r.read(1 << 20)
            if not b:
                break
            f.write(b)
os.makedirs(os.path.join(a.dest, "repoRT"), exist_ok=True)
zp = os.path.join(a.dest, "repoRT", "RepoRT-v1.0.zip")
if not os.path.exists(zp):
    fetch(URL, zp)
h = hashlib.sha256(open(zp, "rb").read()).hexdigest()
if h != SHA256:
    sys.exit(f"checksum mismatch: {h} != {SHA256}")
print("checksum OK")
root = os.path.join(a.dest, "repoRT")
with zipfile.ZipFile(zp) as z:
    members = [m for m in z.namelist() if "/raw_data/" not in m and not m.endswith(".pdf")]
    z.extractall(root, members)
print("extracted to", root, "->", [d for d in os.listdir(root) if d.startswith("f-kretschmer")])

if not a.no_norman:
    nd = os.path.join(a.dest, "norman_ilc"); os.makedirs(os.path.join(nd, "extracted"), exist_ok=True)
    nz = os.path.join(nd, "malm_kruve_2025_supp.zip")
    if not os.path.exists(nz):
        fetch(NORMAN_URL, nz)
    ex = os.path.join(nd, "extracted")
    with zipfile.ZipFile(nz) as z:
        z.extract("AN-150-D5AN00323G-s002.xlsx", ex); z.extract("AN-150-D5AN00323G-s001.zip", ex)
    with zipfile.ZipFile(os.path.join(ex, "AN-150-D5AN00323G-s001.zip")) as z:
        z.extractall(ex, [m for m in z.namelist() if m.startswith("SI/data/") and m.endswith(".csv")])
    for rel, ref in NORMAN_FILES.items():
        h = hashlib.sha256(open(os.path.join(ex, rel), "rb").read()).hexdigest()
        if h != ref:
            sys.exit(f"NORMAN checksum mismatch for {rel}: {h} != {ref}")
    print("NORMAN checksums OK; data extracted to", ex)
