"""One-command reproduction of the study.

  python run_all.py --mode quick   # pipeline check (tiny models, few systems; minutes) -> ../results_quick
  python run_all.py --mode full    # full pre-registered study (several CPU-hours on 4 cores) -> ../results, ../tables, ../figures

Stages: data (download + checksum + build), tests, experiments (main + follow-ups; resumable), analysis, manuscript.
"""
import argparse, os, subprocess, sys, glob, time

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable


def sh(cmd):
    print("+", " ".join(cmd), flush=True)
    r = subprocess.run(cmd, cwd=HERE)
    if r.returncode != 0:
        sys.exit(f"stage failed: {' '.join(cmd)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["quick", "full"], default="quick")
    ap.add_argument("--stages", default="data,tests,experiments,analysis,manuscript")
    a = ap.parse_args(); stages = a.stages.split(","); cfg = f"config/{a.mode}.yaml"
    t0 = time.time()
    raw = os.path.join(HERE, "data", "raw"); proc = os.path.join(HERE, "data", "processed")
    if "data" in stages:
        roots = glob.glob(os.path.join(raw, "repoRT", "f-kretschmer-RepoRT-*"))
        if not roots:
            sh([PY, "scripts/download_data.py", "--dest", "data/raw"])
            roots = glob.glob(os.path.join(raw, "repoRT", "f-kretschmer-RepoRT-*"))
        if not os.path.exists(os.path.join(proc, "features.npz")):
            sh([PY, "scripts/build_dataset.py", "--repo", roots[0], "--out", proc])
    if "tests" in stages:
        sh([PY, "-m", "pytest", "-q", "tests"])
    if "experiments" in stages:
        sh([PY, "experiments/run_main.py", cfg])
        sh([PY, "experiments/run_followups.py", cfg] + (["robust", "smrt", "hilic", "sens"]))
    out = "../results_quick" if a.mode == "quick" else ".."
    rawdir = "../results_quick/raw" if a.mode == "quick" else "../results/raw"
    if "analysis" in stages:
        sh([PY, "experiments/analyze.py", "--raw", rawdir, "--out", out] + (["--B", "500"] if a.mode == "quick" else []))
    if "manuscript" in stages and a.mode == "full" and os.path.exists(os.path.join(HERE, "scripts", "build_manuscript.py")):
        sh([PY, "scripts/build_manuscript.py"])
        sh([PY, "scripts/build_supplementary.py"])
    print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
