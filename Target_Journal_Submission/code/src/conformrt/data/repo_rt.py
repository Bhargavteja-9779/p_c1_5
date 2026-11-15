"""Loading and cleaning of RepoRT v1.0 (Kretschmer et al., Nat. Methods 2024; Zenodo 10.5281/zenodo.16267768).

Design decisions (documented in the manuscript, Section 2.1):
* Unit of a "chromatographic system" = one RepoRT dataset id (column + eluents + gradient + lab).
* Compound identity = first InChIKey block (connectivity; 14 chars). Stereoisomers are merged because
  achiral RP columns rarely separate enantiomers; if entries sharing a connectivity block differ in RT by
  more than `dup_tol_min` within a system (e.g. separated diastereomers or mis-assignments) the block is
  dropped from that system, otherwise RTs are averaged.
* Retention times are in minutes (RepoRT convention).
"""
from __future__ import annotations
import glob, os
import numpy as np
import pandas as pd

ACIDIC = ["formic", "acetic", "trifluoroacetic", "phosphor"]
BASIC = ["nh4oh", "nh4bicarb", "nh4carb", "trieth", "triprop", "tribut", "nndimethylhex"]
SALT = ["nh4ac", "nh4form", "nh4f"]


def _num(x):
    try:
        v = float(x)
        return v if np.isfinite(v) else 0.0
    except Exception:
        return 0.0


def classify_mobile_phase(md: pd.Series) -> dict:
    """Coarse, rule-based description of the mobile phase from RepoRT metadata.
    pH class: 'acidic' if any acid additive in A or B; 'basic' if any basic additive (and no acid);
    'buffered' if only ammonium salts; 'none' otherwise. Modifier: MeOH / ACN / mixed / other."""
    acid = any(_num(md.get(f"eluent.{e}.{a}")) > 0 for e in "AB" for a in ACIDIC)
    base = any(_num(md.get(f"eluent.{e}.{a}")) > 0 for e in "AB" for a in BASIC)
    salt = any(_num(md.get(f"eluent.{e}.{a}")) > 0 for e in "AB" for a in SALT)
    ph = "acidic" if acid else ("basic" if base else ("buffered" if salt else "none"))
    meoh = sum(_num(md.get(f"eluent.{e}.meoh")) for e in "ABCD")
    acn = sum(_num(md.get(f"eluent.{e}.acn")) for e in "ABCD")
    other = sum(_num(md.get(f"eluent.{e}.{s}")) for e in "ABCD" for s in ["iproh", "acetone", "hex", "chcl3", "ch2cl2", "hept"])
    if meoh > 0 and acn > 0:
        mod = "mixed"
    elif acn > 0:
        mod = "acn"
    elif meoh > 0:
        mod = "meoh"
    else:
        mod = "other"
    if other > 0 and mod != "other":
        mod = mod + "+other"
    return {"ph_class": ph, "modifier": mod}


def load_system(path: str, dup_tol_min: float = 0.2) -> tuple[pd.DataFrame, dict] | None:
    sid = os.path.basename(path)
    try:
        info = pd.read_csv(f"{path}/{sid}_info.tsv", sep="\t", dtype=str).iloc[0]
        md = pd.read_csv(f"{path}/{sid}_metadata.tsv", sep="\t").iloc[0]
        rt = pd.read_csv(f"{path}/{sid}_rtdata_canonical_success.tsv", sep="\t")
    except Exception:
        return None
    try:
        grad = pd.read_csv(f"{path}/{sid}_gradient.tsv", sep="\t")
        tcol = [c for c in grad.columns if c.lower().startswith("t")][0]
        grad_dur = float(pd.to_numeric(grad[tcol], errors="coerce").max())
    except Exception:
        grad_dur = np.nan
    rt = rt[["id", "name", "formula", "rt", "smiles.std", "inchikey.std", "classyfire.superclass", "classyfire.class"]].copy()
    rt["rt"] = pd.to_numeric(rt["rt"], errors="coerce")
    rt = rt.dropna(subset=["rt", "smiles.std", "inchikey.std"])
    rt = rt[rt["rt"] > 0]
    rt["ik14"] = rt["inchikey.std"].str[:14]
    g = rt.groupby("ik14")
    spread = g["rt"].agg(lambda s: s.max() - s.min())
    keep = spread[spread <= dup_tol_min].index
    agg = g.agg(rt=("rt", "mean"), smiles=("smiles.std", "first"), formula=("formula", "first"), name=("name", "first"),
                superclass=("classyfire.superclass", "first"), cls=("classyfire.class", "first"), n_entries=("rt", "size"))
    agg = agg.loc[agg.index.isin(keep)].reset_index()
    n_dropped_dups = int((spread > dup_tol_min).sum())
    mp = classify_mobile_phase(md)
    meta = {"system": sid, "name": info.get("name"), "method_type": info.get("method.type"), "source": info.get("source"),
            "column": md.get("column.name"), "usp": md.get("column.usp.code"), "col_len": _num(md.get("column.length")),
            "col_id": _num(md.get("column.id")), "particle": _num(md.get("column.particle.size")),
            "temp": _num(md.get("column.temperature")), "flow": _num(md.get("column.flowrate")), "t0": _num(md.get("column.t0")),
            "grad_dur": grad_dur, "n_raw": int(len(rt)), "n_compounds": int(len(agg)), "n_dropped_dup_blocks": n_dropped_dups,
            "rt_min": float(agg.rt.min()) if len(agg) else np.nan, "rt_max": float(agg.rt.max()) if len(agg) else np.nan, **mp}
    agg.insert(0, "system", sid)
    return agg, meta


def load_repo_rt(root: str, method_types=("RP",), min_compounds: int = 1, dup_tol_min: float = 0.2):
    """Return (long table of system x compound RTs, system metadata table)."""
    tabs, metas = [], []
    for p in sorted(glob.glob(os.path.join(root, "processed_data", "*"))):
        r = load_system(p, dup_tol_min)
        if r is None:
            continue
        tab, meta = r
        if method_types and meta["method_type"] not in method_types:
            continue
        if meta["n_compounds"] < min_compounds:
            continue
        tabs.append(tab); metas.append(meta)
    return pd.concat(tabs, ignore_index=True), pd.DataFrame(metas)
