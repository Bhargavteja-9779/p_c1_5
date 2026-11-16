"""Command-line RT windows for candidate structures on a new LC system.
Usage: python scripts/predict_windows.py --calibrants cal.csv --candidates cand.csv [--system system.yaml] [--alpha 0.10]
       [--method jackknife+|normalised] [--bundle models/conformrt_rp_bundle.joblib] --out windows.csv
cal.csv: columns smiles, rt (min); cand.csv: columns id, smiles[, rt_obs]; system.yaml (normalised method only): ph_class,
modifier, usp, col_len, particle, temp, flow, grad_dur, t0 (missing values allowed). See src/conformrt/api.py."""
import argparse, os, sys
import pandas as pd, yaml
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from conformrt.api import load_bundle, predict_windows  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--calibrants", required=True); ap.add_argument("--candidates", required=True); ap.add_argument("--system")
ap.add_argument("--alpha", type=float, default=0.10); ap.add_argument("--method", default="jackknife+", choices=["jackknife+", "normalised"])
ap.add_argument("--bundle", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "models", "conformrt_rp_bundle.joblib"))
ap.add_argument("--out", required=True); a = ap.parse_args()
meta = yaml.safe_load(open(a.system)) if a.system else None
w = predict_windows(pd.read_csv(a.calibrants), pd.read_csv(a.candidates), load_bundle(a.bundle), alpha=a.alpha, method=a.method, system_meta=meta)
w.to_csv(a.out, index=False); print(f"{len(w)} windows written to {a.out}" + (f"; kept {int(w.keep.sum())}" if "keep" in w else ""))
