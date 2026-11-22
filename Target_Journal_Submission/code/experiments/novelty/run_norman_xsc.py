"""Addendum G3: external validation of cross-system calibration (XSC) on the 37 NORMAN interlaboratory RP methods.

Multipliers are calibrated on all RepoRT RP systems (each scored by the base model that excluded its own source group and
σ trained without that group), replaying (i) n random calibrants (n = 5, 10, 20) and (ii) the NORMAN calibration mix as
available on the system (fixed-mix rule). They are then applied without change to the NORMAN methods (base model and σ
trained on all RepoRT systems), with random calibrants (draws of the main NORMAN run) and with the laboratory mix
(calibration chemicals → suspects).
Usage: python experiments/novelty/run_norman_xsc.py [config]
"""
import sys, os, re, json, zlib, time
import numpy as np, pandas as pd, yaml
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_xsc as X  # noqa: E402
from joblib import Parallel, delayed
from scipy import stats
from rdkit import Chem
from conformrt.evaluation.system_eval import select_calibrants, _sigma_matrix, _sigma_cal_loo
from conformrt.evaluation.intervals import jackknife_plus, winkler
from conformrt.evaluation.crosssys import (loo_residuals, scale_pivot, pivot_summary, choose_avg, choose_lab, choose_lab_exact,
                                            C_GRID, EPS_LIST)
from conformrt.evaluation.protocol import simulate_calibration_rows
from conformrt.models.difficulty import DifficultyModel, encode_system
from conformrt.preprocessing.features import mol_from_smiles

cfg, data, fs, ex = X.cfg, X.data, X.fs, X.ex
quick = cfg["name"] == "quick"
OUT = os.path.join(cfg["out_dir"], "raw", "xsc")
# MASK=1: replay with the metadata the NORMAN methods lack (temperature, flow, dead time) blanked on the RepoRT systems, so that
# the multiplier is calibrated for the information the laboratory actually has (Addendum G6, after review).
MASK = os.environ.get("XSC_MASK_META", "0") == "1"; SUF = "_masked" if MASK else ""
N_RAND = [5, 10, 20]; FAMS = ("flat", "sigs"); DTR = 2 if quick else 8; DEV = 2 if quick else 20
# ---- reuse the NORMAN data preparation of experiments/run_norman.py (everything before the decoy section) ----
src = open("experiments/run_norman.py").read()
prefix = src.split("# ---------------- decoys")[0]
prefix = prefix.replace("ex = Experiment(dict(cfg, name=cfg[\"name\"] + \"_norman\")); data = ex.data; fs = data.fs", "")
g_ = {"__name__": "norman_prep", "__file__": os.path.abspath("experiments/run_norman.py"), "cfg": cfg, "ex": ex, "data": data, "fs": fs, "sys": sys}
exec(compile(prefix.replace('cfg = yaml.safe_load', 'cfg_unused = yaml.safe_load'), "run_norman_prefix", "exec"), g_)
d, uniq, pc, SYS, sig_all, model_all = g_["d"], g_["uniq"], g_["pc"], g_["SYS"], g_["sig"], g_["model"]
sim_all = g_["sim"]
NORMAN_CAL_IK = set(uniq.ik14[uniq.compound.isin(d[d.sample_type == "cal"].compound)])


def replay(job):
    y, sid = job["y"], job["sid"]; N = len(y); out = {}
    pool = np.flatnonzero(np.isin(job["ik"], list(NORMAN_CAL_IK)))
    plans = [("random", n) for n in N_RAND if N >= n + 20]
    if len(pool) >= 10 and N - len(pool) >= 10:
        plans.append(("norman", len(pool)))
    for rule, n in plans:
        acc = {}
        for r in range(1 if rule == "norman" else DTR):
            if rule == "norman":
                cal = pool
            else:
                rng = np.random.default_rng([int(sid), n, r, 0, 99]); cal = rng.choice(N, n, replace=False)
            te = np.setdiff1d(np.arange(N), cal)
            q = X.draw_quantities(job, cal, te)
            for fam in FAMS:
                acc.setdefault(fam, []).append(pivot_summary(q[fam][2]))
        for fam, lst in acc.items():
            out[(rule, n if rule == "random" else -1, fam)] = (np.mean([c for c, _ in lst], 0), {e: np.array([q_[e] for _, q_ in lst]) for e in EPS_LIST})
    return sid, out


if __name__ == "__main__":
    t0 = time.time()
    sim_path = os.path.join(OUT, f"norman_replay{SUF}.pkl")
    if os.path.exists(sim_path):
        outs = pd.read_pickle(sim_path)
    else:
        sim_all["sfold"] = data.meta.loc[sim_all.system.astype(str), "fold"].values
        jobs = []
        for g in sorted(data.meta.fold.unique()):
            sub = sim_all[sim_all.sfold != g]
            sig_g = DifficultyModel(n_jobs=cfg["n_jobs"], **cfg["difficulty_model"]).fit(sub, sub["target"].values)
            sig_g.model.set_params(n_jobs=1)
            systems = [s_ for s_ in data.meta[(data.meta.fold == g) & (data.meta.n_compounds >= 40)].index if s_ != data.smrt_id]
            if quick:
                systems = systems[:5]
            sc = ex.scores_and_novelty(ex.fm.get(frozenset({g})), systems)
            for sid in systems:
                j = X.make_job(sid, sc, sig_g, {g})
                if MASK:
                    j["sysf"] = dict(j["sysf"], temp=np.nan, flow=np.nan, t0=np.nan)
                j["ik"] = data.rt[data.rt.system == sid].ik14.values; jobs.append(j)
        outs = Parallel(n_jobs=cfg["n_jobs"], backend="loky")(delayed(replay)(j) for j in jobs)
        pd.to_pickle(outs, sim_path)
    keys = sorted(set(k for _, o in outs for k in o))
    ctab = []
    for k in keys:
        per = [o[k] for _, o in outs if k in o]
        row = {"rule": k[0], "n": k[1], "family": k[2], "S": len(per)}
        for a in (0.05, 0.10, 0.20):
            row[f"c_avg_{a}"] = choose_avg(np.array([p[0] for p in per]), a)
        for e in EPS_LIST:
            row[f"c_lab_{e}"] = choose_lab([p[1][e] for p in per], e, X.DELTA)
            row[f"c_labx_{e}"] = choose_lab_exact(np.array([p[1][e][0] for p in per]), X.DELTA)
        ctab.append(row)
    ctab = pd.DataFrame(ctab); ctab.to_csv(os.path.join(OUT, f"norman_multipliers{SUF}.csv"), index=False)
    cd = {(r["rule"], r["n"], r["family"]): r for r in ctab.to_dict("records")}
    print("replay done", time.time() - t0, flush=True)

    key = uniq.set_index("compound"); rows = []
    for cs, g in d.groupby("CS_code"):
        g = g[g.compound.isin(key.index)].groupby("compound", as_index=False).agg(RT=("RT", "mean"), sample_type=("sample_type", "first"))
        if len(g) < 20 or cs not in SYS:
            continue
        order = [uniq.index[uniq.compound == c][0] for c in g.compound]
        C = pc.iloc[order].reset_index(drop=True); y = g.RT.values.astype(float); s = uniq.score.values[order]
        pool = np.flatnonzero(g.sample_type.values == "cal"); sus = np.flatnonzero(g.sample_type.values == "sus")
        sysf = encode_system(dict(SYS[cs]))
        for c_ in ["col_len", "particle", "temp", "flow", "grad_dur", "t0"]:
            sysf[c_] = SYS[cs][c_]
        job = {"y": y, "s": s, "C": C, "sysf": sysf, "sig": sig_all}
        sid = zlib.crc32(cs.encode()) % 10 ** 6
        scen = [("fixed_mix", "norman", -1, [(pool, sus)])]
        for n in N_RAND:
            draws = []
            for r in range(DEV):
                rng = np.random.default_rng([int(sid), n, r, 0, cfg.get("seed", 0)])
                cal = rng.choice(len(y), n, replace=False); draws.append((cal, np.setdiff1d(np.arange(len(y)), cal)))
            scen.append(("random_all", "random", n, draws))
        for name, rule, n, draws in scen:
            for r, (cal, te) in enumerate(draws):
                q = X.draw_quantities(job, cal, te, want_P=len(cal) in (10, 20))
                if "mu_loo_te" not in q and len(cal) >= 10:  # fixed mix of variable size: jackknife+ comparators too
                    q["mu_loo_te"] = np.vstack([m.predict(s[te]) for m in q["loo"]])
                    q["P_sl_te"] = _sigma_matrix(sig_all, C.iloc[te].reset_index(drop=True), s[te], s[cal], y[cal], sysf, loo=True)
                    q["P_R"] = np.abs(q["loo_res"]) / _sigma_cal_loo(sig_all, C.iloc[cal].reset_index(drop=True), s[cal], y[cal], sysf)
                y_te = y[te]; z = stats.norm.ppf(0.95); rmse = np.sqrt(np.mean(q["loo_res"] ** 2)); mu = q["flat"][0]
                W = {"B1_gauss_rmse": (mu - z * rmse, mu + z * rmse)}
                if "mu_loo_te" in q:
                    W["B5_jackknife_plus"] = jackknife_plus(q["mu_loo_te"], np.abs(q["loo_res"]), None, 0.10)
                    W["P"] = jackknife_plus(q["mu_loo_te"], q["P_R"], q["P_sl_te"], 0.10)
                for fam in FAMS:
                    c = cd.get((rule, n, fam))
                    if c is None:
                        continue
                    m_, tau, _ = q[fam]; sc = tau * (1.0 if fam == "flat" else q["sig_te"])
                    for tgt in ("avg", "lab"):
                        k = c[f"c_{tgt}_0.1"]; W[f"XSC_{fam}_{tgt}"] = (m_ - k * sc, m_ + k * sc)
                for meth, (lo, hi) in W.items():
                    cov = (lo <= y_te) & (y_te <= hi); w = hi - lo; fin = bool(np.isfinite(w).all())
                    rows.append({"cs": cs, "scenario": name, "n": len(cal), "draw": r, "method": meth, "coverage": float(cov.mean()),
                                 "width_median": float(np.median(w)) if fin else np.inf,
                                 "winkler": float(winkler(lo, hi, y_te, 0.10).mean()) if fin else np.inf, "span": float(y.max() - y.min())})
    pd.DataFrame(rows).to_parquet(os.path.join(OUT, f"norman_eval{SUF}.parquet"))
    print("NORMAN_XSC_DONE", time.time() - t0, flush=True)
