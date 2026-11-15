"""Build the processed RepoRT tables and the compound feature cache.
Usage: python scripts/build_dataset.py --repo <RepoRT root> --out <dir>"""
import argparse, os, sys, time, json, hashlib
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from conformrt.data.repo_rt import load_repo_rt
from conformrt.preprocessing.features import featurize

ap = argparse.ArgumentParser(); ap.add_argument("--repo", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--n_jobs", type=int, default=4); a = ap.parse_args()
os.makedirs(a.out, exist_ok=True); t = time.time()
rt, meta = load_repo_rt(a.repo, method_types=("RP", "HILIC"))
rt.to_parquet(f"{a.out}/rt_long.parquet"); meta.to_csv(f"{a.out}/systems.csv", index=False)
comp = rt.groupby("ik14").agg(smiles=("smiles", "first"), formula=("formula", "first"), superclass=("superclass", "first"),
                               cls=("cls", "first"), n_systems=("system", "nunique")).reset_index()
F = featurize(comp["smiles"], n_jobs=a.n_jobs)
comp["feat_ok"] = F["ok"]
comp.to_parquet(f"{a.out}/compounds.parquet")
np.savez_compressed(f"{a.out}/features.npz", desc=F["desc"], ion=F["ion"], morgan=F["morgan"], bits=F["bits"], elem=F["elem"],
                    desc_names=np.array(F["desc_names"]), ion_names=np.array(F["ion_names"]), elem_names=np.array(F["elem_names"]))
summ = {"n_rows": int(len(rt)), "n_systems": int(rt.system.nunique()), "n_compounds": int(len(comp)), "n_feat_fail": int((~F["ok"]).sum()),
        "seconds": round(time.time() - t, 1)}
json.dump(summ, open(f"{a.out}/build_summary.json", "w"), indent=1); print(summ)
