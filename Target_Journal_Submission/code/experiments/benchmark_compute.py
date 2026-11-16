"""E9: computational cost of the proposed workflow, measured single-threaded after all other runs (no contention).
Measures: base-model training (600 trees, all RP data), difficulty-model training, per-system calibration for n calibrants
(LOO monotone maps + LOO difficulty scales), per-candidate window computation, peak memory.
Usage: python experiments/benchmark_compute.py [config]"""
import sys, os, time, json, resource, yaml, platform
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
os.chdir(os.path.join(os.path.dirname(__file__), ".."))
from conformrt.evaluation.protocol import Data, simulate_calibration_rows
from conformrt.models.base_models import RetentionOrderModel
from conformrt.models.difficulty import DifficultyModel, encode_system
from conformrt.models.mapping import MonotoneMap
from conformrt.evaluation.system_eval import _sigma_matrix, _sigma_cal_loo
from conformrt.evaluation.intervals import jackknife_plus

cfg = yaml.safe_load(open(sys.argv[1] if len(sys.argv) > 1 else "config/full.yaml"))
data = Data(cfg["processed_dir"]); fs = data.fs; out = {"platform": platform.platform(), "python": platform.python_version(), "threads": 1}
tr = data.training_table(set(), cfg["smrt_n"])
t = time.perf_counter(); base = RetentionOrderModel(n_jobs=1, **cfg["base_model"]).fit(fs, tr); out["base_train_s_1thread"] = time.perf_counter() - t
out["base_train_rows"] = int(len(tr)); out["base_train_systems"] = int(tr.system.nunique())
# difficulty model trained on simulated rows of a subset of systems (timing only)
systems = [s for s in data.meta.index if s != data.smrt_id][:150]
rows = data.fs.rows(data.rt[data.rt.system.isin(systems)].ik14.unique()); sc = pd.Series(base.predict(fs, rows), index=rows)
scores = {s: sc.loc[fs.rows(data.rt[data.rt.system == s].ik14)].values for s in systems}
nov = {s: (np.ones(len(v), np.float32), np.ones(len(v), np.float32)) for s, v in scores.items()}
sim = simulate_calibration_rows(data, systems, scores, nov, n_list=(10, 20, 40), draws=3)
t = time.perf_counter(); dm = DifficultyModel(n_jobs=1, **cfg["difficulty_model"]).fit(sim, sim["target"].values); out["sigma_train_s_1thread"] = time.perf_counter() - t
out["sigma_train_rows"] = int(len(sim)); dm.model.set_params(n_jobs=1)
# per-system calibration and per-candidate windows on a large target system
sid = data.targets(60).sort_values("n_compounds").index[-1]
tab = data.rt[data.rt.system == sid]; r = fs.rows(tab.ik14); y = tab.rt.values; s = base.predict(fs, r)
C = fs.physchem(r); C["maxsim"] = 1.0; C["seen"] = 1.0; sysf = encode_system(data.meta.loc[sid].to_dict())
res = {}
for n in (10, 20, 40):
    rng = np.random.default_rng(0); cal = rng.choice(len(y), n, replace=False); te = np.setdiff1d(np.arange(len(y)), cal)[:1000]
    t = time.perf_counter()
    loo = [MonotoneMap().fit(np.delete(s[cal], i), np.delete(y[cal], i)) for i in range(n)]
    mu_c = np.array([loo[i].predict(s[cal][i:i + 1])[0] for i in range(n)])
    sl_cal = _sigma_cal_loo(dm, C.iloc[cal].reset_index(drop=True), s[cal], y[cal], sysf); R = np.abs(y[cal] - mu_c) / sl_cal
    t_cal = time.perf_counter() - t
    t = time.perf_counter()
    mu_t = np.vstack([m.predict(s[te]) for m in loo]); sl_te = _sigma_matrix(dm, C.iloc[te].reset_index(drop=True), s[te], s[cal], y[cal], sysf, loo=True)
    lo, hi = jackknife_plus(mu_t, R, sl_te, 0.1); t_win = time.perf_counter() - t
    res[n] = {"calibration_ms": 1000 * t_cal, "windows_for_candidates": int(len(te)), "windows_ms_total": 1000 * t_win, "windows_us_per_candidate": 1e6 * t_win / len(te)}
out["per_system"] = res; out["benchmark_system"] = sid
out["peak_rss_MB"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
os.makedirs(os.path.join(cfg["out_dir"], "experiment_logs"), exist_ok=True)
json.dump(out, open(os.path.join(cfg["out_dir"], "experiment_logs", "compute_benchmark.json"), "w"), indent=1); print(json.dumps(out, indent=1))
