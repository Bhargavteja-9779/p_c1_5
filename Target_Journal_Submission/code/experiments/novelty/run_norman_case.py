"""Addendum G5: real-world case study — suspect screening in the 37 NORMAN interlaboratory NTS methods.

Each laboratory calibrates with its own NORMAN calibration chemicals (as in the interlaboratory study). For every suspect
detected by the laboratory, the candidate list is every PubChemLite-for-Exposomics structure with the suspect's molecular
formula (cap 200 per formula, fixed seed). Windows: Gaussian ±z·RMSE, jackknife+, normalised jackknife+ (meta-learned σ),
XSC (σ ≡ 1 and shrunk σ; multipliers from the replayed 'NORMAN mix' rule on RepoRT, average and laboratory-level targets).
Records per laboratory and suspect: correct structure retained, candidates eliminated, number kept, rank of the correct
structure by standardised RT deviation. Usage: python experiments/novelty/run_norman_case.py [config]
(after run_norman_xsc.py, which writes norman_multipliers.csv)"""
import sys, os, zlib, time
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_norman_xsc as NX  # noqa: E402  (NORMAN preparation: uniq, pc, SYS, sig_all, model_all, d)
from scipy import stats
from rdkit import DataStructs
import xgboost as xgb
from conformrt.evaluation.system_eval import _sigma_matrix, _sigma_cal_loo
from conformrt.evaluation.intervals import jackknife_plus
from conformrt.evaluation.crosssys import loo_residuals, scale_pivot
from conformrt.models.difficulty import encode_system
from conformrt.preprocessing.features import featurize

X, cfg, data, fs = NX.X, NX.cfg, NX.data, NX.fs
uniq, pc, SYS, sig, model, d = NX.uniq, NX.pc, NX.SYS, NX.sig_all, NX.model_all, NX.d
OUT = NX.OUT; PCL = "data/raw/pubchemlite/PubChemLite_exposomics_20260925.csv"; MAX_CAND = 200; A = 0.10
PCKEYS = list(pc.columns)


def candidate_table():
    path = os.path.join(X.ex.cache_root, "pubchemlite", "norman_cand.parquet")
    if os.path.exists(path):
        return pd.read_parquet(path), np.load(path.replace(".parquet", "_feat.npz"))
    forms = set(uniq.formula.dropna())
    p = pd.read_csv(PCL, usecols=["FirstBlock", "MolecularFormula", "SMILES"])
    p = p[p.MolecularFormula.isin(forms)].drop_duplicates("FirstBlock").sample(frac=1.0, random_state=0)
    p = p[p.groupby("MolecularFormula").cumcount() < MAX_CAND].reset_index(drop=True)
    F = featurize(p.SMILES, n_jobs=cfg["n_jobs"]); ok = F["ok"]; p = p[ok].reset_index(drop=True)
    Xb = np.hstack([F["desc"][ok], F["morgan"][ok].astype(np.float32)])
    p["score"] = model.booster.predict(xgb.DMatrix(Xb))
    ref = [data.bitvects()[r] for r in fs.rows(NX.g_["train_ik"])]
    bits = F["bits"][ok]; ms = []
    for row in bits:
        v = DataStructs.ExplicitBitVect(1024); v.SetBitsFromList(np.flatnonzero(np.unpackbits(row)).tolist())
        ms.append(max(DataStructs.BulkTanimotoSimilarity(v, ref)))
    seen = p.FirstBlock.isin(NX.g_["train_ik"]).values
    p["maxsim"] = np.where(seen, 1.0, ms); p["seen"] = seen.astype(np.float32)
    keep = ["MolLogP", "TPSA", "MolWt", "NumHDonors", "NumHAcceptors", "NumRotatableBonds", "FractionCSP3", "NumAromaticRings", "RingCount", "HeavyAtomCount", "MolMR", "LabuteASA"]
    C = pd.DataFrame(F["desc"][ok][:, [F["desc_names"].index(k) for k in keep]], columns=keep)
    for j, n in enumerate(F["ion_names"]): C["ion_" + n] = F["ion"][ok][:, j]
    for j, n in enumerate(F["elem_names"]): C[n] = F["elem"][ok][:, j]
    C["n_acidic_groups"] = C[["ion_carboxylic_acid", "ion_sulfonic_acid", "ion_phosphate_acid", "ion_tetrazole_or_acidic_NH"]].sum(1)
    C["n_basic_groups"] = C[["ion_aliphatic_amine_prim", "ion_aliphatic_amine_sec", "ion_aliphatic_amine_tert", "ion_quaternary_N", "ion_amidine_guanidine"]].sum(1)
    C["maxsim"] = p.maxsim.values; C["seen"] = p.seen.values
    C = C[PCKEYS]
    p.to_parquet(path); np.savez_compressed(path.replace(".parquet", "_feat.npz"), C=C.to_numpy(np.float32))
    return p, np.load(path.replace(".parquet", "_feat.npz"))


if __name__ == "__main__":
    t0 = time.time()
    cand, cf = candidate_table(); CC = pd.DataFrame(cf["C"], columns=PCKEYS)
    by_formula = cand.groupby("MolecularFormula").indices
    cd = {(r["rule"], r["n"], r["family"]): r for r in pd.read_csv(os.path.join(OUT, f"norman_multipliers{NX.SUF}.csv")).to_dict("records")}
    key = uniq.set_index("compound"); out = []
    for cs, g in d.groupby("CS_code"):
        g = g[g.compound.isin(key.index)].groupby("compound", as_index=False).agg(RT=("RT", "mean"), sample_type=("sample_type", "first"))
        if len(g) < 20 or cs not in SYS:
            continue
        order = [uniq.index[uniq.compound == c][0] for c in g.compound]
        C = pc.iloc[order].reset_index(drop=True); y = g.RT.values.astype(float); s = uniq.score.values[order]
        cal = np.flatnonzero(g.sample_type.values == "cal"); sus = np.flatnonzero(g.sample_type.values == "sus")
        sysf = encode_system(dict(SYS[cs]))
        for c_ in ["col_len", "particle", "temp", "flow", "grad_dur", "t0"]:
            sysf[c_] = SYS[cs][c_]
        forms = uniq.formula.values[order]; iks = uniq.ik14.values[order]
        par, is_true, qs, qC = [], [], [], []
        for j in sus:
            rows = by_formula.get(forms[j], np.array([], int))
            rows = rows[cand.FirstBlock.values[rows] != iks[j]]
            if not len(rows):
                continue
            par += [j] * (len(rows) + 1); is_true += [True] + [False] * len(rows)
            qs.append(np.r_[s[j], cand.score.values[rows]]); qC.append(pd.concat([C.iloc[[j]], CC.iloc[rows]]))
        if not qs:
            continue
        par = np.array(par); is_true = np.array(is_true); qs = np.concatenate(qs); Cq = pd.concat(qC).reset_index(drop=True); yq = y[par]
        s_cal, y_cal = s[cal], y[cal]; Cc = C.iloc[cal].reset_index(drop=True)
        loo_res, loo = loo_residuals(s_cal, y_cal)
        sig_q = _sigma_matrix(sig, Cq, qs, s_cal, y_cal, sysf, loo=False); sig_cal_loo = _sigma_cal_loo(sig, Cc, s_cal, y_cal, sysf)
        mu, tau_f, _ = scale_pivot(s_cal, y_cal, qs, None, None, None, loo_res)
        _, tau_s, _ = scale_pivot(s_cal, y_cal, qs, None, sig_q, sig_cal_loo, loo_res, shrink=X.SHRINK)
        mu_loo = np.vstack([m.predict(qs) for m in loo]); sl_q = _sigma_matrix(sig, Cq, qs, s_cal, y_cal, sysf, loo=True)
        z = stats.norm.ppf(1 - A / 2); rmse = np.sqrt(np.mean(loo_res ** 2)); one = np.ones_like(mu)
        W = {"B1_gauss_rmse": (mu - z * rmse, mu + z * rmse, one),
             "RMSE_x2": (mu - 2 * rmse, mu + 2 * rmse, one),
             "B5_jackknife_plus": jackknife_plus(mu_loo, np.abs(loo_res), None, A) + (one,),
             "P": jackknife_plus(mu_loo, np.abs(loo_res) / sig_cal_loo, sl_q, A) + (sig_q,)}
        for fam, tau, sc in (("flat", tau_f, one), ("sigs", tau_s, sig_q)):
            c = cd.get(("norman", -1, fam))
            if c is None:
                continue
            for tgt in ("avg", "lab"):
                k = c[f"c_{tgt}_{A}"]; W[f"XSC_{fam}_{tgt}"] = (mu - k * tau * sc, mu + k * tau * sc, sc)
        df = pd.DataFrame({"par": par, "true": is_true})
        for meth, (lo, hi, sc) in W.items():
            dd = df.assign(keep=(lo <= yq) & (yq <= hi), dev=np.abs(yq - mu) / sc)
            for p_, q in dd.groupby("par"):
                tq = q[q.true].iloc[0]; fq = q[~q.true]
                out.append({"cs": cs, "method": meth, "suspect": int(p_), "n_cal": len(cal), "true_kept": bool(tq.keep),
                            "n_cand": len(fq), "elim": float((~fq.keep).mean()), "n_kept": int(fq.keep.sum()),
                            "rank": int((fq.dev.values < tq.dev).sum()) + 1})
    pd.DataFrame(out).to_parquet(os.path.join(OUT, f"norman_case{NX.SUF}.parquet"))
    print("NORMAN_CASE_DONE", time.time() - t0, flush=True)
