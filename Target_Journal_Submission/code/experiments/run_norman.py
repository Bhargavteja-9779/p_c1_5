"""E5c (Addendum C): external validation on the NORMAN interlaboratory RT data (37 RP NTS methods).
Usage: python experiments/run_norman.py [config] [--norman_dir data/raw/norman_ilc/extracted]"""
import sys, os, re, json, yaml, time, zlib
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
os.chdir(os.path.join(os.path.dirname(__file__), ".."))
from joblib import Parallel, delayed
from rdkit import Chem, DataStructs
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from conformrt.evaluation.driver import Experiment, save_results
from conformrt.evaluation.protocol import simulate_calibration_rows
from conformrt.evaluation.system_eval import evaluate_system
from conformrt.models.difficulty import DifficultyModel, encode_system
from conformrt.preprocessing.features import featurize, mol_from_smiles
from conformrt.statistics.analysis import per_system

cfg = yaml.safe_load(open(sys.argv[1] if len(sys.argv) > 1 and sys.argv[1].endswith(".yaml") else "config/full.yaml"))
ND = "data/raw/norman_ilc/extracted"
quick = cfg["name"] == "quick"
ex = Experiment(dict(cfg, name=cfg["name"] + "_norman")); data = ex.data; fs = data.fs

# ---------------- NORMAN data and features ----------------
d = pd.read_csv(f"{ND}/SI/data/NORMAN_interlab_RT_data.csv")
d = d.dropna(subset=["SMILES", "RT"]); d = d[d.RT > 0]
uniq = d[["compound", "SMILES"]].drop_duplicates("compound").reset_index(drop=True)
mols = [mol_from_smiles(s) for s in uniq.SMILES]
uniq["ik14"] = [Chem.MolToInchiKey(m)[:14] if m is not None else None for m in mols]
uniq["formula"] = [CalcMolFormula(m) if m is not None else None for m in mols]
uniq["smiles_std"] = [Chem.MolToSmiles(m) if m is not None else None for m in mols]
uniq = uniq.dropna(subset=["ik14"]).reset_index(drop=True)
F = featurize(uniq["smiles_std"], n_jobs=cfg["n_jobs"])
X_base = np.hstack([F["desc"], F["morgan"].astype(np.float32)])
pc = pd.DataFrame(F["desc"][:, [F["desc_names"].index(k) for k in ["MolLogP", "TPSA", "MolWt", "NumHDonors", "NumHAcceptors", "NumRotatableBonds", "FractionCSP3", "NumAromaticRings", "RingCount", "HeavyAtomCount", "MolMR", "LabuteASA"]]],
                  columns=["MolLogP", "TPSA", "MolWt", "NumHDonors", "NumHAcceptors", "NumRotatableBonds", "FractionCSP3", "NumAromaticRings", "RingCount", "HeavyAtomCount", "MolMR", "LabuteASA"])
for j, n in enumerate(F["ion_names"]): pc["ion_" + n] = F["ion"][:, j]
for j, n in enumerate(F["elem_names"]): pc[n] = F["elem"][:, j]
pc["n_acidic_groups"] = pc[["ion_carboxylic_acid", "ion_sulfonic_acid", "ion_phosphate_acid", "ion_tetrazole_or_acidic_NH"]].sum(1)
pc["n_basic_groups"] = pc[["ion_aliphatic_amine_prim", "ion_aliphatic_amine_sec", "ion_aliphatic_amine_tert", "ion_quaternary_N", "ion_amidine_guanidine"]].sum(1)

# ---------------- base model on ALL RepoRT RP systems ----------------
model, train_ik = ex.fm.get(frozenset())
import xgboost as xgb
score = model.booster.predict(xgb.DMatrix(X_base))
ref = [data.bitvects()[r] for r in fs.rows(train_ik)]
qbv = []
for row in F["bits"]:
    v = DataStructs.ExplicitBitVect(1024); v.SetBitsFromList(np.flatnonzero(np.unpackbits(row)).tolist()); qbv.append(v)
seen = np.isin(uniq.ik14.values, train_ik)
maxsim = np.array([1.0 if seen[i] else max(DataStructs.BulkTanimotoSimilarity(q, ref)) for i, q in enumerate(qbv)], np.float32)
pc["maxsim"] = maxsim; pc["seen"] = seen.astype(np.float32)
uniq["score"] = score; uniq["maxsim"] = maxsim; uniq["seen"] = seen

# ---------------- difficulty model on all RepoRT systems (out-of-source scores) ----------------
rows = []
for f in sorted(data.meta.fold.unique()):
    outer = ex.fm.get(frozenset({f}))
    systems = [s for s in data.meta[data.meta.fold == f].index if s != data.smrt_id]
    nov = ex.scores_and_novelty(outer, systems)
    scores, novelty = {}, {}
    for sid in systems:
        r = fs.rows(data.rt[data.rt.system == sid].ik14)
        scores[sid] = nov.loc[r, "score"].values; novelty[sid] = (nov.loc[r, "maxsim"].values, nov.loc[r, "seen"].values.astype(np.float32))
    rows.append(simulate_calibration_rows(data, systems, scores, novelty, n_list=tuple(cfg["sim_n_list"]), draws=cfg["sim_draws"], seed=cfg.get("seed", 0)))
sim = pd.concat(rows, ignore_index=True)
sig = DifficultyModel(groups=["compound", "ion", "novelty", "geometry", "system"], n_jobs=cfg["n_jobs"], **cfg["difficulty_model"]).fit(sim, sim["target"].values)
sig.model.set_params(n_jobs=1)

# ---------------- system descriptors from ESI Table S1 ----------------
t1 = pd.read_excel(f"{ND}/AN-150-D5AN00323G-s002.xlsx", sheet_name="Table S1", header=None)
hdr = [str(h) for h in t1.iloc[1].tolist()]; t1 = t1.iloc[2:].reset_index(drop=True); t1.columns = hdr
glen = pd.read_csv(f"{ND}/SI/data/gradient_lenght_codes.csv", sep=r"\s+").drop_duplicates("CS_code").set_index("CS_code")["Lenghts"]


def sys_feats(row):
    mp = (str(row.get("Mobile phase A")) + " " + str(row.get("Mobile phase B"))).lower()
    acid = any(k in mp for k in ["formic acid", "acetic acid", "trifluoroacetic"])
    base = any(k in mp for k in ["ammonium hydroxide", "ammonia", "ammonium bicarbonate", "ammonium carbonate"])
    salt = any(k in mp for k in ["ammonium acetate", "ammonium formate", "ammonium fluoride"])
    ph = "acidic" if acid else ("basic" if base else ("buffered" if salt else "none"))
    meoh = "meoh" in mp or "methanol" in mp; acn = "mecn" in mp or "acetonitrile" in mp or "acn" in mp
    mod = "mixed" if (meoh and acn) else ("acn" if acn else ("meoh" if meoh else "other"))
    sp = str(row.get("Column stationary phase")).split(";")[0].strip().upper()
    usp = "L1" if "C18" in sp else ("L7" if "C8" in sp else ("L11" if "PHENYL" in sp else "other"))
    dims = str(row.get("Column dimensions")).split(";")[0]
    m = re.search(r"(\d+(?:\.\d+)?)\s*[×x]\s*(\d+(?:\.\d+)?)\s*mm,?\s*(\d+(?:\.\d+)?)\s*µm", dims)
    col_len, particle = (float(m.group(1)), float(m.group(3))) if m else (np.nan, np.nan)
    return {"ph_class": ph, "modifier": mod, "usp": usp, "col_len": col_len, "particle": particle, "temp": np.nan,
            "flow": np.nan, "grad_dur": float(glen.get(row["CS"], np.nan)), "t0": np.nan}


SYS = {r["CS"]: sys_feats(r) for _, r in t1.iterrows()}
json.dump({k: {kk: (None if (isinstance(vv, float) and np.isnan(vv)) else vv) for kk, vv in v.items()} for k, v in SYS.items()},
          open(os.path.join(cfg["out_dir"], "experiment_logs", "norman_system_descriptors.json"), "w"), indent=1)

# ---------------- decoys (RepoRT isomers of NORMAN compounds) ----------------
comp = fs.comp
fg = comp.groupby("formula").indices

# ---------------- jobs ----------------
key = uniq.set_index("compound")
alphas = cfg["alphas"]; draws = 2 if quick else 20
jobs = []
methods = ["B1_gauss_rmse", "B2_pspline_PI", "B3_gp_PI", "B4_split_conformal", "B5_jackknife_plus"]
for cs, g in d.groupby("CS_code"):
    g = g[g.compound.isin(key.index)].groupby("compound", as_index=False).agg(RT=("RT", "mean"), sample_type=("sample_type", "first"))
    if len(g) < 20 or cs not in SYS:
        continue
    k = key.loc[g.compound]; idx = uniq.index[uniq.compound.isin(g.compound)]
    order = [uniq.index[uniq.compound == c][0] for c in g.compound]
    C = pc.iloc[order].reset_index(drop=True); y = g.RT.values.astype(float); s = uniq.score.values[order]
    pool = np.flatnonzero(g.sample_type.values == "cal"); sus = np.flatnonzero(g.sample_type.values == "sus")
    sysf = encode_system(dict(SYS[cs]))
    for c_ in ["col_len", "particle", "temp", "flow", "grad_dur", "t0"]:
        sysf[c_] = SYS[cs][c_]
    strata = {"sample_type": g.sample_type.values.astype(object), "seen": np.where(uniq.seen.values[order], "seen", "unseen"),
              "ionisable": np.where((C.n_acidic_groups + C.n_basic_groups).values > 0, "ionisable", "neutral")}
    # decoys for suspects
    par, drows_ = [], []
    rng = np.random.default_rng(zlib.crc32(cs.encode()))
    for j in sus:
        frm = uniq.formula.values[order[j]]; ik = uniq.ik14.values[order[j]]
        if frm in fg:
            cand = fg[frm]; cand = cand[comp.ik14.values[cand] != ik]
            if len(cand) > cfg["max_decoys"]:
                cand = rng.choice(cand, cfg["max_decoys"], replace=False)
            par += [j] * len(cand); drows_ += list(cand)
    dec = None
    if drows_:
        drows_ = np.array(drows_); dsc = model.predict(fs, drows_); dnov = ex.scores_and_novelty((model, train_ik), [], extra_rows=drows_)
        Cd = fs.physchem(drows_); Cd["maxsim"] = dnov.loc[drows_, "maxsim"].values; Cd["seen"] = dnov.loc[drows_, "seen"].values.astype(np.float32)
        dec = {"parent": np.array(par), "s": dsc, "C": Cd}
    base_cfg = {"alphas": alphas, "min_test": 10, "seed": cfg.get("seed", 0), "methods": methods, "expensive_n": [], "strata": True}
    scen = [("fixed_mix", "fixed", [len(pool)], 1, pool, sus), ("pool_random", "pool_random", [10, 20], draws, pool, sus),
            ("random_all", "random", [10, 20], draws, None, None)]
    for name, strat, nl, dr, pl, tidx in scen:
        cfg_j = dict(base_cfg, n_list=nl, draws=dr, primary_n=nl[0] if strat == "fixed" else 20, strategy=strat,
                     strategy_code={"fixed": 3, "pool_random": 4, "random": 0}[strat], expensive_n=nl,
                     decoy_draws=(1 if strat == "fixed" else 0))
        job = {"sid": str(zlib.crc32(cs.encode()) % 10 ** 6), "cs": cs, "scenario": name, "y": y, "s": s, "C": C, "ik": uniq.ik14.values[order],
               "sysf": sysf, "sig_models": {"P": sig}, "ens": None, "strata": strata, "cfg": cfg_j, "pool": pl, "test_idx": tidx,
               "decoys": dec if strat == "fixed" else None}
        jobs.append(job)
print("NORMAN jobs:", len(jobs), "systems:", len({j['cs'] for j in jobs}), flush=True)
outs = Parallel(n_jobs=cfg["n_jobs"])(delayed(evaluate_system)(j) for j in jobs)
res = {k: [] for k in ["rows", "strata", "decoys", "errors", "timing"]}
for j, o in zip(jobs, outs):
    for k in res:
        for r in o[k]:
            r["cs"] = j["cs"]; r["scenario"] = j["scenario"]; r["fold"] = -1
        res[k] += o[k]
# optional (Addendum D): selection-aware system-level calibration factors learned on RepoRT with the NORMAN-mix rule
sasc_file = os.path.join(cfg["out_dir"], "raw", "sasc_norman_c.json")
if os.path.exists(sasc_file):
    from conformrt.evaluation.sasc import loo_quantities, windows
    from conformrt.evaluation.intervals import winkler as _wk
    from conformrt.evaluation.system_eval import _sigma_matrix as _sm, _sigma_cal_loo as _scl
    cN = json.load(open(sasc_file))
    for j in jobs:
        if j["scenario"] != "fixed_mix":
            continue
        y, s_, C_, sysf_, cal, te = j["y"], j["s"], j["C"], j["sysf"], np.asarray(j["pool"]), np.asarray(j["test_idx"])
        mu, R, sl = loo_quantities(s_[cal], y[cal], s_[te], lambda st, sc, yc: (_sm(sig, C_.iloc[te].reset_index(drop=True), st, sc, yc, sysf_, loo=True),
                                                                              _scl(sig, C_.iloc[cal].reset_index(drop=True), sc, yc, sysf_)))
        muJ, RJ, _ = loo_quantities(s_[cal], y[cal], s_[te], None)
        for a_ in alphas:
            for meth, (m_, R_, sl_, c_) in {"P+SASC": (mu, R, sl, cN["c_P"]), "J+SASC": (muJ, RJ, None, cN["c_J"])}.items():
                if a_ != 0.1:
                    continue  # factors were learned for alpha = 0.10
                lo, hi = windows(m_, R_, sl_, a_, c_); cov = (lo <= y[te]) & (y[te] <= hi)
                res["rows"].append({"system": j["sid"], "n": len(cal), "draw": 0, "method": meth, "alpha": a_, "m_test": len(te), "coverage": cov.mean(),
                                    "width_mean": (hi - lo).mean(), "width_median": np.median(hi - lo), "width_rel_mean": (hi - lo).mean() / (y.max() - y.min()),
                                    "winkler": _wk(lo, hi, y[te], a_).mean(), "cs": j["cs"], "scenario": "fixed_mix", "fold": -1})
res = {k: pd.DataFrame(v) for k, v in res.items()}
raw = os.path.join(cfg["out_dir"], "raw"); os.makedirs(raw, exist_ok=True)
save_results(res, raw, "norman")
meta = {"n_systems": int(res["rows"].cs.nunique()), "n_compounds": int(len(uniq)), "seen_share": float(seen.mean()),
        "sigma_rows": int(len(sim))}
json.dump(meta, open(os.path.join(cfg["out_dir"], "experiment_logs", "norman_meta.json"), "w"), indent=1)
print("done", meta)
