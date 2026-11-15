import sys, os, time, yaml
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
os.chdir(os.path.join(os.path.dirname(__file__), ".."))
from conformrt.evaluation.driver import Experiment, save_results
cfg = yaml.safe_load(open("config/quick.yaml"))
t = time.time(); ex = Experiment(cfg)
res = ex.run(base_kind="order", systems_limit=cfg["systems_limit"], tag="smoke")
save_results(res, cfg["out_dir"], "smoke")
print({k: len(v) for k, v in res.items()}, f"{time.time()-t:.0f}s")
r = res["rows"]
print(r.groupby("method")[["coverage", "width_mean", "winkler"]].mean().round(3).to_string())
print(res["decoys"].groupby("method")[["decoy_elim_rate", "true_retained"]].mean().round(3).to_string() if len(res["decoys"]) else "no decoys")
print(res["strata"].groupby(["stratum"]).total.sum().to_dict())
