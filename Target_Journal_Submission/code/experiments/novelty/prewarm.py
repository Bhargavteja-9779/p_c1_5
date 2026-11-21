"""Re-create the cached fold base models and simulated difficulty-model rows (cache is not versioned)."""
import sys, os, yaml
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from conformrt.evaluation.driver import Experiment
cfg = yaml.safe_load(open(sys.argv[1]))
ex = Experiment(cfg)
for f in range(5):
    ex.fm.get(frozenset({f}))
    ex.difficulty_training_rows(f, "order")
    print("fold", f, "ready", flush=True)
print("PREWARM_DONE", flush=True)
