"""Download RepoRT v1.0 from Zenodo and verify its SHA-256 checksum.
Usage: python scripts/download_data.py --dest data/raw"""
import argparse, hashlib, os, sys, zipfile, urllib.request
URL = "https://zenodo.org/api/records/16267768/files/f-kretschmer/RepoRT-v1.0.zip/content"
SHA256 = "4eeeb46a65ec4950eca67f652caed5cffbdf5bcb5dffc45e496edae0552c6091"
ap = argparse.ArgumentParser(); ap.add_argument("--dest", default="data/raw"); a = ap.parse_args()
os.makedirs(os.path.join(a.dest, "repoRT"), exist_ok=True)
zp = os.path.join(a.dest, "repoRT", "RepoRT-v1.0.zip")
if not os.path.exists(zp):
    print("downloading", URL)
    req = urllib.request.Request(URL, headers={"User-Agent": "conformrt-reproduction/1.0"})
    with urllib.request.urlopen(req) as r, open(zp, "wb") as f:
        while True:
            b = r.read(1 << 20)
            if not b: break
            f.write(b)
h = hashlib.sha256(open(zp, "rb").read()).hexdigest()
if h != SHA256:
    sys.exit(f"checksum mismatch: {h} != {SHA256}")
print("checksum OK")
root = os.path.join(a.dest, "repoRT")
with zipfile.ZipFile(zp) as z:
    members = [m for m in z.namelist() if "/raw_data/" not in m and not m.endswith(".pdf")]
    z.extractall(root, members)
print("extracted to", root, "->", [d for d in os.listdir(root) if d.startswith("f-kretschmer")])
