"""Train the deployable model bundle used by conformrt.api / scripts/predict_windows.py (reversed phase).
* base model: retention-order model trained on ALL RepoRT RP systems (+ 15,000-compound SMRT subsample), as in the article;
* difficulty model: trained on simulated calibrant protocols on all RP systems, each source group scored by the model trained
  without that group (inner cross-fitting, as in the article);
* novelty reference: InChIKey blocks and fingerprints of the base-model training compounds.
Output: models/conformrt_rp_bundle.joblib. Usage: python scripts/train_deployable_models.py [config/full.yaml]"""
import os, sys, time, json, platform
import numpy as np, pandas as pd, yaml, joblib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "src")); os.chdir(os.path.join(HERE, ".."))
from conformrt.evaluation.driver import Experiment  # noqa: E402
from conformrt.evaluation.protocol import simulate_calibration_rows  # noqa: E402
from conformrt.models.difficulty import DifficultyModel  # noqa: E402

cfg = yaml.safe_load(open(sys.argv[1] if len(sys.argv) > 1 else "config/full.yaml"))
t0 = time.time(); ex = Experiment(dict(cfg, name=cfg["name"] + "_deploy")); data = ex.data
base, train_ik = ex.fm.get(frozenset())  # nothing excluded
parts = []
for g in sorted(data.meta.fold.unique()):
    systems = [s for s in data.meta[data.meta.fold == g].index if s != data.smrt_id]
    nov = ex.scores_and_novelty(ex.fm.get(frozenset({g})), systems)
    scores, novelty = {}, {}
    for sid in systems:
        r = data.fs.rows(data.rt[data.rt.system == sid].ik14)
        scores[sid] = nov.loc[r, "score"].values; novelty[sid] = (nov.loc[r, "maxsim"].values, nov.loc[r, "seen"].values.astype(np.float32))
    parts.append(simulate_calibration_rows(data, systems, scores, novelty, n_list=tuple(cfg["sim_n_list"]), draws=cfg["sim_draws"], seed=cfg.get("seed", 0)))
    print(f"group {g}: {len(parts[-1])} simulated rows", flush=True)
rows = pd.concat(parts, ignore_index=True)
sig = DifficultyModel(n_jobs=cfg["n_jobs"], **cfg["difficulty_model"]).fit(rows, rows["target"].values); sig.model.set_params(n_jobs=1)
import xgboost, lightgbm, rdkit
bundle = {"base_model": base.booster, "train_ik": np.asarray(train_ik), "train_bits": data.fs.bits[data.fs.rows(train_ik)],
          "difficulty_model": sig, "separation_mode": "RP",
          "provenance": {"data": "RepoRT v1.0 (doi:10.5281/zenodo.16267768)", "config": cfg, "n_sigma_rows": int(len(rows)),
                         "versions": {"python": platform.python_version(), "xgboost": xgboost.__version__, "lightgbm": lightgbm.__version__,
                                      "rdkit": rdkit.__version__}, "build_seconds": None}}
bundle["provenance"]["build_seconds"] = round(time.time() - t0, 1)
os.makedirs("models", exist_ok=True); joblib.dump(bundle, "models/conformrt_rp_bundle.joblib", compress=3)
print("bundle written:", os.path.getsize("models/conformrt_rp_bundle.joblib") / 1e6, "MB;", bundle["provenance"]["build_seconds"], "s", flush=True)
