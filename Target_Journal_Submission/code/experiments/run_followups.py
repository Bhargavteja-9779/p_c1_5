"""E4 robustness, E5 generalisation (SMRT base model; HILIC systems), E7 sensitivity (difficulty-model hyper-parameters).
Usage: python experiments/run_followups.py [config] [which ...]   which ∈ {robust, smrt, hilic, sens}"""
import sys, os, time, json, copy, yaml
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
os.chdir(os.path.join(os.path.dirname(__file__), ".."))
from conformrt.evaluation.driver import Experiment, save_results
cfg_path = sys.argv[1] if len(sys.argv) > 1 else "config/full.yaml"
which = sys.argv[2:] or ["robust", "smrt", "hilic", "sens"]
base = yaml.safe_load(open(cfg_path)); raw = os.path.join(base["out_dir"], "raw"); os.makedirs(raw, exist_ok=True)
quick = base["name"] == "quick"
meta = {}


def sub(**kw):
    c = copy.deepcopy(base); c.update(kw); return c


if "robust" in which:  # E4: n = 20, 10 draws, main sigma model only
    c = sub(name=base["name"] + "_robust", n_list=[base["primary_n"]], draws=(2 if quick else 10), decoy_draws=0, expensive_n=[base["primary_n"]])
    ex = Experiment(c)
    for strategy, perturb in [("span", None), ("common", None), ("random", "noise"), ("random", "outlier")]:
        t = time.time(); tag = f"robust_{strategy}_{perturb or 'none'}"
        res = ex.run(base_kind="order", strategy=strategy, perturb=perturb, sigma_which=["P"], systems_limit=c.get("systems_limit"), tag=tag, ckpt_dir=os.path.join(raw, 'checkpoints'))
        save_results(res, raw, tag); meta[tag] = round(time.time() - t, 1)

if "smrt" in which:  # E5a: second base model (SMRT regressor), n in {10, 20, 40}
    c = sub(name=base["name"] + "_smrt", n_list=[base["primary_n"]] if quick else [10, 20, 40], decoy_draws=0)
    ex = Experiment(c); t = time.time()
    res = ex.run(base_kind="smrt", sigma_which=["P"], methods=["B1_gauss_rmse", "B2_pspline_PI", "B3_gp_PI", "B4_split_conformal", "B5_jackknife_plus"],
                 systems_limit=c.get("systems_limit"), tag="gen_smrt", ckpt_dir=os.path.join(raw, "checkpoints"))
    save_results(res, raw, "gen_smrt"); meta["gen_smrt"] = round(time.time() - t, 1); meta["smrt_train_s"] = ex.timings.get("smrt_train_s")

if "hilic" in which:  # E5b: HILIC systems (different separation mechanism), base model trained on other-source HILIC systems
    c = sub(name=base["name"] + "_hilic", method_type="HILIC", n_list=[base["primary_n"]] if quick else [10, 20, 40], decoy_draws=0)
    ex = Experiment(c); t = time.time()
    res = ex.run(base_kind="order", sigma_which=["P"], systems_limit=c.get("systems_limit"), tag="gen_hilic", ckpt_dir=os.path.join(raw, "checkpoints"))
    save_results(res, raw, "gen_hilic"); meta["gen_hilic"] = round(time.time() - t, 1)

if "sens" in which:  # E7: difficulty-model hyper-parameters (n = 20, 10 draws)
    for nl, ne in [(15, 400), (63, 400), (31, 200), (31, 800)]:
        dm = dict(base["difficulty_model"]); dm.update(num_leaves=nl, n_estimators=ne)
        c = sub(name=base["name"] + f"_sens_l{nl}_e{ne}", n_list=[base["primary_n"]], draws=(2 if quick else 10), decoy_draws=0,
                difficulty_model=dm, methods=["B5_jackknife_plus"], expensive_n=[])
        ex = Experiment(c); t = time.time(); tag = f"sens_l{nl}_e{ne}"
        res = ex.run(base_kind="order", sigma_which=["P"], methods=["B5_jackknife_plus"], systems_limit=c.get("systems_limit"), tag=tag, ckpt_dir=os.path.join(raw, 'checkpoints'))
        save_results(res, raw, tag); meta[tag] = round(time.time() - t, 1)

os.makedirs(os.path.join(base["out_dir"], "experiment_logs"), exist_ok=True)
json.dump(meta, open(os.path.join(base["out_dir"], "experiment_logs", f"followups_{base['name']}_{'_'.join(which)}.json"), "w"), indent=1)
print("done", meta)
