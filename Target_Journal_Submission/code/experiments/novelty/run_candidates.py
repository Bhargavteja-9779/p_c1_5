"""Addendum G2: candidate filtering against real candidate lists (PubChemLite for Exposomics, 2026-09-25).

For each test compound of a target system, the candidate list is every PubChemLite structure with the same molecular
formula (different InChIKey first block), as retrieved for an LC-HRMS feature by an accurate-mass/formula query.
At n = 20 random calibrants (draws 0-4 of the main experiment, α = ε = 0.10) we record, per method: whether the correct
structure is retained, the fraction of candidates eliminated, the number retained, and the rank of the correct structure
when candidates are ordered by the standardised RT deviation |t_obs − μ(c)| / scale(c).
Usage: python experiments/novelty/run_candidates.py [config]   (after run_xsc.py; uses its fold multipliers)
"""
import sys, os, time
import numpy as np, pandas as pd, yaml
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_xsc as X  # noqa: E402  (sets cwd to code/, builds the Experiment)
from joblib import Parallel, delayed
from scipy import stats
from rdkit import DataStructs
from conformrt.evaluation.system_eval import select_calibrants, _sigma_matrix, _sigma_cal_loo
from conformrt.evaluation.intervals import jackknife_plus
from conformrt.evaluation import sasc as S
from conformrt.evaluation.crosssys import loo_residuals, scale_pivot
from conformrt.preprocessing.features import featurize
from conformrt.evaluation.driver import _compound_frame
from conformrt.evaluation.protocol import max_similarity

cfg, data, fs, ex = X.cfg, X.data, X.fs, X.ex
PCL = "data/raw/pubchemlite/PubChemLite_exposomics_20260925.csv"
OUT = os.path.join(cfg["out_dir"], "raw", "xsc"); CACHE = os.path.join(X.ex.cache_root, "pubchemlite")
os.makedirs(CACHE, exist_ok=True)
MAX_CAND = 200       # per formula (random subset, fixed seed) to bound compute; recorded in the output
DRAWS = 5; N = 20; A = 0.10
quick = cfg["name"] == "quick"


def extend_feature_store():
    """Append featurised PubChemLite candidates (formulas of RP target compounds) to the FeatureStore."""
    path = os.path.join(CACHE, "cand_features.npz")
    targets = data.targets(cfg["min_target_compounds"])
    tik = data.rt[data.rt.system.isin(targets.index)].ik14.unique()
    forms = set(data.formula.reindex(tik).dropna())
    if not os.path.exists(path):
        p = pd.read_csv(PCL, usecols=["FirstBlock", "MolecularFormula", "SMILES"])
        p = p[p.MolecularFormula.isin(forms)].drop_duplicates("FirstBlock")
        p = p[~p.FirstBlock.isin(fs.index.index)]  # structures already in RepoRT keep their RepoRT features
        p = p.sample(frac=1.0, random_state=0)
        p = p[p.groupby("MolecularFormula").cumcount() < MAX_CAND]
        F = featurize(p.SMILES.reset_index(drop=True), n_jobs=cfg["n_jobs"])
        ok = F["ok"]
        np.savez_compressed(path, ik14=p.FirstBlock.values[ok], formula=p.MolecularFormula.values[ok], desc=F["desc"][ok],
                            morgan=F["morgan"][ok], bits=F["bits"][ok], ion=F["ion"][ok], elem=F["elem"][ok])
    z = np.load(path, allow_pickle=True)
    n0 = len(fs.comp)
    fs.desc = np.vstack([fs.desc, z["desc"].astype(np.float32)]); fs.morgan = np.vstack([fs.morgan, z["morgan"]])
    fs.bits = np.vstack([fs.bits, z["bits"]]); fs.ion = np.vstack([fs.ion, z["ion"].astype(np.float32)])
    fs.elem = np.vstack([fs.elem, z["elem"].astype(np.float32)])
    fs.index = pd.concat([fs.index, pd.Series(np.arange(n0, n0 + len(z["ik14"])), index=z["ik14"])])
    data._bitvects = None
    # candidate lists by formula: RepoRT compounds + PubChemLite additions present in PubChemLite
    p = pd.read_csv(PCL, usecols=["FirstBlock", "MolecularFormula"])
    p = p[p.MolecularFormula.isin(forms)].drop_duplicates("FirstBlock")
    p = p[p.FirstBlock.isin(fs.index.index)]
    rows = fs.index.loc[p.FirstBlock.values].values
    return {f_: r for f_, r in pd.Series(rows, index=p.MolecularFormula.values).groupby(level=0)}


def eval_system(job, ctab, cand):
    y, s, sid, ik = job["y"], job["s"], job["sid"], job["ik"]; Nn = len(y); out = []
    cd = {(r["rule"], r["n"], r["family"]): r for r in ctab.to_dict("records")}
    for r in range(DRAWS):
        rng = np.random.default_rng([int(sid), N, r, X.STRATEGY_CODE["random"], cfg.get("seed", 0)])
        cal = select_calibrants(rng, y, N, "random", job["nsys"], None, job["logp"]); te = np.setdiff1d(np.arange(Nn), cal)
        te = np.array([j for j in te if len(cand.get(j, ())) > 0], int)
        if len(te) == 0:
            continue
        par = np.concatenate([[j] * (1 + len(cand[j])) for j in te])           # true structure + its candidates
        is_true = np.concatenate([[True] + [False] * len(cand[j]) for j in te])
        cs = np.concatenate([[s[j]] + list(job["cand_s"][cand[j]]) for j in te])
        Cq = pd.concat([pd.concat([job["C"].iloc[[j]], job["cand_C"].iloc[cand[j]]]) for j in te]).reset_index(drop=True)
        yq = y[par]
        s_cal, y_cal = s[cal], y[cal]; Cc = job["C"].iloc[cal].reset_index(drop=True)
        loo_res, loo = loo_residuals(s_cal, y_cal)
        sig = job["sig"]; sysf = job["sysf"]
        sig_q = _sigma_matrix(sig, Cq, cs, s_cal, y_cal, sysf, loo=False)
        sig_cal_loo = _sigma_cal_loo(sig, Cc, s_cal, y_cal, sysf)
        mu, tau_f, _ = scale_pivot(s_cal, y_cal, cs, None, None, None, loo_res)
        _, tau_s, _ = scale_pivot(s_cal, y_cal, cs, None, sig_q, sig_cal_loo, loo_res, shrink=X.SHRINK)
        mu_loo = np.vstack([m.predict(cs) for m in loo])
        sl_q = _sigma_matrix(sig, Cq, cs, s_cal, y_cal, sysf, loo=True)
        z = stats.norm.ppf(1 - A / 2); rmse = np.sqrt(np.mean(loo_res ** 2))
        W = {"B1_gauss_rmse": (mu - z * rmse, mu + z * rmse, np.ones_like(mu)),
             "B5_jackknife_plus": jackknife_plus(mu_loo, np.abs(loo_res), None, A) + (np.ones_like(mu),),
             "P": jackknife_plus(mu_loo, np.abs(loo_res) / sig_cal_loo, sl_q, A) + (sig_q,)}
        for fam, tau, sc in (("flat", tau_f, np.ones_like(mu)), ("sigs", tau_s, sig_q)):
            c = cd[("random", N, fam)]
            for tgt in ("avg", "lab"):
                k = c[f"c_{tgt}_{A}"]
                W[f"XSC_{fam}_{tgt}"] = (mu - k * tau * sc, mu + k * tau * sc, sc)
        c = cd.get(("random", N, "jp"))
        if c is not None:
            W["JP_lab"] = S.windows(mu_loo, np.abs(loo_res), None, 0.10, c[f"c_lab_{A}"]) + (np.ones_like(mu),)
        df = pd.DataFrame({"par": par, "true": is_true})
        for meth, (lo, hi, sc) in W.items():
            keep = (lo <= yq) & (yq <= hi)
            dev = np.abs(yq - mu) / sc
            d = df.assign(keep=keep, dev=dev)
            g = d.groupby("par")
            tr = d[d.true].set_index("par")
            rank = g.apply(lambda q: int((q.dev[~q.true].values < q.dev[q.true].values[0]).sum()) + 1)
            n_c = g.true.size() - 1
            elim = g.apply(lambda q: float((~q.keep[~q.true]).mean()))
            kept = g.apply(lambda q: int(q.keep[~q.true].sum()))
            out.append(pd.DataFrame({"system": sid, "draw": r, "method": meth, "par": tr.index, "true_kept": tr.keep.values,
                                     "n_cand": n_c.loc[tr.index].values, "elim": elim.loc[tr.index].values,
                                     "n_kept": kept.loc[tr.index].values, "rank": rank.loc[tr.index].values}))
    return pd.concat(out) if out else pd.DataFrame()


if __name__ == "__main__":
    t0 = time.time()
    by_formula = extend_feature_store()
    X.log.info(f"[cand] candidate lists for {len(by_formula)} formulas ({time.time()-t0:.0f}s)")
    targets = data.targets(cfg["min_target_compounds"])
    for f in sorted(targets.fold.unique()):
        pth = os.path.join(OUT, f"fold{f}_candidates.parquet")
        if os.path.exists(pth):
            continue
        ctab = pd.read_csv(os.path.join(OUT, f"fold{f}_multipliers.csv"))
        sig_f = ex.train_sigma_models(f, "order", which=["P"])[0]["P"]
        model, train_ik = ex.fm.get(frozenset({f}))
        tsys = list(targets[targets.fold == f].index)[: (3 if quick else None)]
        sc_t = ex.scores_and_novelty((model, train_ik), tsys)
        # candidate scores and novelty for this fold's model
        allc = np.unique(np.concatenate([v for v in by_formula.values()]))
        cs = model.predict(fs, allc)
        seen = np.isin(np.asarray(fs.index.index)[allc], train_ik)
        ms = max_similarity(data, allc, fs.rows(train_ik), n_jobs=cfg["n_jobs"]); ms = np.where(seen, 1.0, ms)
        cand_s = pd.Series(cs, index=allc); cand_C = _compound_frame(data, allc, ms, seen); cand_C.index = allc
        jobs = []
        for sid in tsys:
            j = X.make_job(sid, sc_t, sig_f, {f})
            tab = data.rt[data.rt.system == sid].reset_index(drop=True)
            j["ik"] = tab.ik14.values
            forms = data.formula.reindex(tab.ik14.values).values
            own = set(fs.rows(tab.ik14.values))
            cand = {}
            for jj, fm in enumerate(forms):
                rows = by_formula.get(fm)
                if rows is None:
                    continue
                rows = rows[(np.asarray(fs.index.index)[rows] != tab.ik14.values[jj])]
                if len(rows):
                    cand[jj] = rows
            u = np.unique(np.concatenate(list(cand.values()))) if cand else np.array([], int)
            pos = {r_: k for k, r_ in enumerate(u)}
            j["cand"] = {jj: np.array([pos[r_] for r_ in v]) for jj, v in cand.items()}
            j["cand_s"] = cand_s.loc[u].values; j["cand_C"] = cand_C.loc[u].reset_index(drop=True)
            jobs.append(j)
        res = Parallel(n_jobs=cfg["n_jobs"], backend="loky")(delayed(eval_system)(j, ctab, j["cand"]) for j in jobs)
        df = pd.concat([r for r in res if len(r)]); df["fold"] = f; df["max_cand_per_formula"] = MAX_CAND
        df.to_parquet(pth)
        X.log.info(f"[cand] fold {f}: {len(tsys)} systems, {len(df)} rows ({time.time()-t0:.0f}s)")
    print("CANDIDATES_DONE", flush=True)
