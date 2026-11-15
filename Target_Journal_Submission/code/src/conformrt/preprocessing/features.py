"""Molecular featurisation (RDKit) for base retention models and the difficulty model.

* `rdkit_descriptors`: all RDKit 2D descriptors (Descriptors._descList, 217 in RDKit 2026.09).
* `morgan_counts`: folded Morgan count fingerprint, radius 2, 1024 bits (for base models) and
  bit-vector fingerprints for Tanimoto similarity (novelty / applicability-domain feature).
* `ionizable_counts`: SMARTS counts of acidic / basic functional groups, which govern pH-dependent
  selectivity changes between mobile phases (a known driver of retention-order changes).
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

IONIZABLE_SMARTS = {
    "carboxylic_acid": "[CX3](=O)[OX2H1,OX1-]",
    "sulfonic_acid": "[SX4](=O)(=O)[OX2H1,OX1-]",
    "phosphate_acid": "[PX4](=O)([OX2H1,OX1-])",
    "phenol": "c[OX2H1]",
    "tetrazole_or_acidic_NH": "[nH]1nnnc1",
    "aliphatic_amine_prim": "[NX3;H2;!$(NC=[O,S,N]);!$(N-a)][CX4]",
    "aliphatic_amine_sec": "[NX3;H1;!$(NC=[O,S,N]);!$(N-a)]([CX4])[CX4]",
    "aliphatic_amine_tert": "[NX3;H0;!$(NC=[O,S,N]);!$(N-a)]([CX4])([CX4])[CX4]",
    "quaternary_N": "[NX4+;!$(N~O)]",
    "amidine_guanidine": "[NX3][CX3]=[NX2]",
    "pyridine_like_N": "[nX2;r6]",
    "imidazole_like": "c1ncc[nH]1",
}
_PATTS = {k: Chem.MolFromSmarts(v) for k, v in IONIZABLE_SMARTS.items()}
_DESC_NAMES = [d[0] for d in Descriptors._descList]
_MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
_UNCHARGER = rdMolStandardize.Uncharger()


def mol_from_smiles(smi: str):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    try:
        m = _UNCHARGER.uncharge(m)
    except Exception:
        pass
    return m


def featurize_one(smi: str):
    m = mol_from_smiles(smi)
    if m is None:
        return None
    desc = []
    for n, f in Descriptors._descList:
        try:
            v = float(f(m))
        except Exception:
            v = np.nan
        desc.append(v)
    ion = [len(m.GetSubstructMatches(p)) for p in _PATTS.values()]
    counts = _MORGAN.GetCountFingerprintAsNumPy(m).astype(np.uint8)
    bits = np.packbits((counts > 0).astype(np.uint8))
    elems = {a.GetSymbol() for a in m.GetAtoms()}
    natoms = {s: sum(1 for a in m.GetAtoms() if a.GetSymbol() == s) for s in ["N", "O", "S", "P", "F", "Cl", "Br", "I"]}
    charge = Chem.GetFormalCharge(m)
    return np.array(desc, dtype=np.float32), np.array(ion, dtype=np.int16), counts, bits, [natoms[s] for s in natoms] + [charge, int("Si" in elems)]


def featurize(smiles: pd.Series, n_jobs: int = 4):
    """Return dict of aligned arrays for the given SMILES (index preserved)."""
    from joblib import Parallel, delayed
    res = Parallel(n_jobs=n_jobs, batch_size=256)(delayed(featurize_one)(s) for s in smiles.tolist())
    ok = np.array([r is not None for r in res])
    nd, ni, nm = len(_DESC_NAMES), len(_PATTS), 1024
    D = np.full((len(res), nd), np.nan, np.float32); I = np.zeros((len(res), ni), np.int16)
    M = np.zeros((len(res), nm), np.uint8); B = np.zeros((len(res), nm // 8), np.uint8); E = np.zeros((len(res), 10), np.int16)
    for k, r in enumerate(res):
        if r is None:
            continue
        D[k], I[k], M[k], B[k], E[k] = r[0], r[1], r[2], r[3], r[4]
    D[~np.isfinite(D)] = np.nan
    return {"ok": ok, "desc": D, "desc_names": _DESC_NAMES, "ion": I, "ion_names": list(_PATTS), "morgan": M, "bits": B,
            "elem": E, "elem_names": ["n_N", "n_O", "n_S", "n_P", "n_F", "n_Cl", "n_Br", "n_I", "formal_charge", "has_Si"]}


def tanimoto_max(query_bits: np.ndarray, ref_bits: np.ndarray, chunk: int = 2048) -> np.ndarray:
    """Max Tanimoto similarity of each packed-bit query fingerprint to a reference set (vectorised popcount)."""
    pc = np.unpackbits(np.arange(256, dtype=np.uint8)[:, None], axis=1).sum(1).astype(np.int32)
    ref_cnt = pc[ref_bits].sum(1)
    out = np.zeros(len(query_bits), np.float32)
    for s in range(0, len(query_bits), chunk):
        q = query_bits[s:s + chunk]
        q_cnt = pc[q].sum(1)
        inter = np.zeros((len(q), len(ref_bits)), np.int32)
        for b in range(q.shape[1]):
            inter += pc[np.bitwise_and(q[:, b][:, None], ref_bits[:, b][None, :])]
        union = q_cnt[:, None] + ref_cnt[None, :] - inter
        sim = np.where(union > 0, inter / np.maximum(union, 1), 0.0)
        out[s:s + chunk] = sim.max(1)
    return out
