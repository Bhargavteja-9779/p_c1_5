"""Download the two public data sets and verify their SHA-256 checksums:
  (1) RepoRT v1.0 (Zenodo record 16267768; CC BY-SA 4.0) -> data/raw/repoRT/
  (2) NORMAN interlaboratory NTS retention-time data: Electronic Supplementary Information of Malm & Kruve, Analyst 150 (2025)
      (doi 10.1039/D5AN00323G; open access, PMC12268317), via the Europe PMC supplementary-files service -> data/raw/norman_ilc/
Usage: python scripts/download_data.py --dest data/raw [--no-norman]"""
import argparse, hashlib, os, sys, zipfile, urllib.request
URL = "https://zenodo.org/api/records/16267768/files/f-kretschmer/RepoRT-v1.0.zip/content"
SHA256 = "4eeeb46a65ec4950eca67f652caed5cffbdf5bcb5dffc45e496edae0552c6091"
NORMAN_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12268317/supplementaryFiles"
NORMAN_SHA256 = "f96093af95a496979fcecc35009d6e8b93cf3fcf3d77f546a05c228cc9da4c15"
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
    h = hashlib.sha256(open(nz, "rb").read()).hexdigest()
    if h != NORMAN_SHA256:
        sys.exit(f"NORMAN checksum mismatch: {h} != {NORMAN_SHA256}")
    print("NORMAN checksum OK")
    ex = os.path.join(nd, "extracted")
    with zipfile.ZipFile(nz) as z:
        z.extract("AN-150-D5AN00323G-s002.xlsx", ex); z.extract("AN-150-D5AN00323G-s001.zip", ex)
    with zipfile.ZipFile(os.path.join(ex, "AN-150-D5AN00323G-s001.zip")) as z:
        z.extractall(ex, [m for m in z.namelist() if m.startswith("SI/data/") and m.endswith(".csv")])
    print("NORMAN data extracted to", ex)
