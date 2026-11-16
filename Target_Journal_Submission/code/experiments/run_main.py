"""E1/E2/E3/E6 main experiment (pre-specified full configuration). Usage: python experiments/run_main.py [config]"""
import sys, os, time, json, hashlib, subprocess, yaml
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
os.chdir(os.path.join(os.path.dirname(__file__), ".."))
from conformrt.evaluation.driver import Experiment, save_results
from conformrt.utils import detect_hardware
cfg_path = sys.argv[1] if len(sys.argv) > 1 else "config/full.yaml"
cfg = yaml.safe_load(open(cfg_path))
raw_dir = os.path.join(cfg["out_dir"], "raw"); os.makedirs(raw_dir, exist_ok=True)
manifest = {"config": cfg_path, "config_sha256": hashlib.sha256(open(cfg_path, "rb").read()).hexdigest(),
            "git_commit": subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip(),
            "hardware": detect_hardware(), "start": time.strftime("%Y-%m-%dT%H:%M:%S")}
t = time.time(); ex = Experiment(cfg)
res = ex.run(base_kind="order", strategy="random", systems_limit=cfg.get("systems_limit"), tag="main", ckpt_dir=os.path.join(raw_dir, "checkpoints"))
save_results(res, raw_dir, "main")
manifest.update({"end": time.strftime("%Y-%m-%dT%H:%M:%S"), "seconds": round(time.time() - t, 1),
                 "timings": {k: (v if not isinstance(v, list) else {"sum": float(sum(v)), "n": len(v)}) for k, v in ex.timings.items()},
                 "n_rows": {k: int(len(v)) for k, v in res.items()}})
os.makedirs(os.path.join(cfg["out_dir"], "experiment_logs"), exist_ok=True)
json.dump(manifest, open(os.path.join(cfg["out_dir"], "experiment_logs", f"manifest_main_{cfg['name']}.json"), "w"), indent=1)
print("done", manifest["seconds"])
