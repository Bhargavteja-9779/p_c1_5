"""Hardware detection, seeding, logging helpers."""
from __future__ import annotations
import os, json, logging, platform, random, time
import numpy as np


def detect_hardware() -> dict:
    info = {"cpu_count": os.cpu_count(), "platform": platform.platform(), "python": platform.python_version(),
            "cuda": False, "mps": False}
    try:  # optional: torch is not required by this project
        import torch  # noqa
        info["cuda"] = bool(torch.cuda.is_available())
        info["mps"] = bool(getattr(torch.backends, "mps", None) and torch.backends.mps.is_available())
    except Exception:
        pass
    # XGBoost can use CUDA if present; LightGBM and the conformal layer run on CPU. Apple MPS is not used by
    # tree libraries, so CPU is the fallback on Apple Silicon.
    info["xgb_device"] = "cuda" if info["cuda"] else "cpu"
    return info


def set_seed(seed: int):
    random.seed(seed); np.random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)


def get_logger(path: str | None = None, name="conformrt"):
    lg = logging.getLogger(name)
    if lg.handlers:
        return lg
    lg.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    h = logging.StreamHandler(); h.setFormatter(fmt); lg.addHandler(h)
    if path:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        fh = logging.FileHandler(path); fh.setFormatter(fmt); lg.addHandler(fh)
    return lg


class Timer:
    def __init__(self): self.t = time.perf_counter()
    def __call__(self): return time.perf_counter() - self.t
