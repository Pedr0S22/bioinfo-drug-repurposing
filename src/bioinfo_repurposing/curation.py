"""Molecular standardization, bioactivity curation, and leakage control.

Reuses and extends Assignment 5 (CAP5) Task A (`largest_fragment_smiles`, `curate_activity`).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Descriptors, inchi
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

_TO_NM = {"nM": 1.0, "uM": 1e3, "µM": 1e3, "M": 1e9}
_UNCHARGER = rdMolStandardize.Uncharger()
_ORGANIC_ATOMIC_NUMS = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 34, 35, 53}


def largest_fragment_smiles(smiles: str | None) -> str | None:
    """Canonical SMILES of the largest fragment (by heavy atoms) of `smiles`; None if invalid.

    Directly reused from Assignment 5 (CAP5).
    """
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return None
    largest = max(frags, key=lambda f: f.GetNumHeavyAtoms())
    return Chem.MolToSmiles(largest, canonical=True)


def standardize_smiles(smiles: str | None, neutralize: bool = True) -> str | None:
    """Strip salts/counter-ions (largest fragment), optionally neutralize charges, return canonical SMILES."""
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return None
    parent = max(frags, key=lambda f: f.GetNumHeavyAtoms())
    if neutralize:
        try:
            parent = _UNCHARGER.uncharge(parent)
            Chem.SanitizeMol(parent)
        except Exception:
            pass
    return Chem.MolToSmiles(parent, canonical=True)


def smiles_to_inchikey14(smiles: str | None) -> str | None:
    """First 14 characters of the InChIKey (connectivity layer) for robust cross-database matching."""
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        ik = inchi.MolToInchiKey(mol)
        return ik[:14] if ik else None
    except Exception:
        return None


def _prepare_numeric_ic50(
    raw: pd.DataFrame,
    standard_type: str = "IC50",
    allowed_relations: tuple[str, ...] = ("=",),
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Shared filtering funnel returning filtered records with `nM`, `pIC50_record`, `canonical_smiles`, `inchikey14`."""
    funnel: dict[str, int] = {"raw_records": int(len(raw))}
    df = raw.copy()

    if "data_validity_comment" in df.columns:
        valid_mask = df["data_validity_comment"].isna() | (df["data_validity_comment"].astype(str).str.strip() == "")
        df = df[valid_mask].copy()
        funnel["after_data_validity_check"] = int(len(df))

    if "assay_type" in df.columns:
        df = df[df["assay_type"] == "B"].copy()
        funnel["after_binding_assay_B"] = int(len(df))

    df = df[df["standard_type"] == standard_type].copy()
    funnel[f"after_type_{standard_type}"] = int(len(df))

    df = df[df["standard_relation"].isin(allowed_relations)].copy()
    funnel["after_relation_filter"] = int(len(df))

    df = df[df["standard_units"].isin(_TO_NM)].copy()
    funnel["after_molar_units"] = int(len(df))

    df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
    df = df[df["standard_value"].notna() & (df["standard_value"] > 0)].copy()
    df["nM"] = df["standard_value"] * df["standard_units"].map(_TO_NM)
    funnel["after_positive_values"] = int(len(df))

    df["canonical_smiles"] = [standardize_smiles(s) for s in df["smiles"]]
    df = df[df["canonical_smiles"].notna()].copy()
    funnel["after_valid_smiles_standardization"] = int(len(df))

    df["inchikey14"] = [smiles_to_inchikey14(s) for s in df["canonical_smiles"]]
    df = df[df["inchikey14"].notna()].copy()
    df["pIC50_record"] = 9.0 - np.log10(df["nM"])
    return df, funnel


def curate_activity(
    raw: pd.DataFrame,
    standard_type: str = "IC50",
    return_funnel: bool = False,
) -> pd.DataFrame | tuple[pd.DataFrame, dict[str, int]]:
    """Curate raw bioactivity records for exact-value pIC50 regression.

    Directly adapted from Assignment 5 Task A:
    - keep `standard_type` (default `IC50`) and `standard_relation == '='`
    - keep valid molar units (nM, uM, M) and convert to nM
    - drop non-positive values
    - standardize SMILES (largest fragment + neutralize) and drop unparseable structures
    - compute pIC50 = 9 - log10(IC50 in nM)
    - aggregate replicates per canonical SMILES using the median.
    """
    df, funnel = _prepare_numeric_ic50(raw, standard_type=standard_type, allowed_relations=("=",))
    grouped = (
        df.groupby("canonical_smiles")
        .agg(
            pIC50=("pIC50_record", "median"),
            n_measurements=("pIC50_record", "size"),
            pIC50_min=("pIC50_record", "min"),
            pIC50_max=("pIC50_record", "max"),
            inchikey14=("inchikey14", "first"),
        )
        .sort_index()
    )
    grouped["pIC50_range"] = grouped["pIC50_max"] - grouped["pIC50_min"]
    funnel["unique_curated_compounds"] = int(len(grouped))
    if return_funnel:
        return grouped, funnel
    return grouped


def replicate_disagreement_report(raw: pd.DataFrame, threshold: float = 2.0) -> pd.DataFrame:
    """Identify compounds whose exact IC50 replicates disagree by more than `threshold` log units (Assignment 5 Q1b)."""
    curated = curate_activity(raw)
    return curated[curated["pIC50_range"] > threshold].sort_values("pIC50_range", ascending=False)


def curate_activity_classification(
    raw: pd.DataFrame,
    active_threshold: float = 6.0,
    inactive_threshold: float = 5.0,
    rescue_censored_inactives: bool = True,
    standard_type: str = "IC50",
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Curate bioactivity table for binary classification (Active vs Inactive).

    - Exact '=' measurements: active if median pIC50 >= `active_threshold`, inactive if <= `inactive_threshold`,
      discarding the ambiguous grey zone between `inactive_threshold` and `active_threshold`.
    - Optionally rescues unambiguous censored inactives (`standard_relation` in {'>', '>='} with
      bound >= 10**(9 - inactive_threshold) nM, i.e., IC50 > 10 uM -> guaranteed pIC50 < 5.0) for compounds
      that do NOT have an exact '=' measurement.
    """
    exact_df, funnel = curate_activity(raw, standard_type=standard_type, return_funnel=True)
    exact_df = exact_df.copy()
    exact_df["source_label"] = "exact"

    exact_actives = exact_df[exact_df["pIC50"] >= active_threshold].copy()
    exact_actives["y"] = 1

    exact_inactives = exact_df[exact_df["pIC50"] <= inactive_threshold].copy()
    exact_inactives["y"] = 0

    grey_zone = exact_df[(exact_df["pIC50"] > inactive_threshold) & (exact_df["pIC50"] < active_threshold)]
    funnel["exact_actives"] = int(len(exact_actives))
    funnel["exact_inactives"] = int(len(exact_inactives))
    funnel["dropped_grey_zone"] = int(len(grey_zone))

    rescued_inactives = pd.DataFrame()
    if rescue_censored_inactives:
        cens_raw, _ = _prepare_numeric_ic50(raw, standard_type=standard_type, allowed_relations=(">", ">="))
        min_inactive_nM = 10.0 ** (9.0 - inactive_threshold)  # 10,000 nM for pIC50 <= 5.0
        cens_inact = cens_raw[cens_raw["nM"] >= min_inactive_nM].copy()
        if len(cens_inact) > 0:
            cens_grouped = (
                cens_inact.groupby("canonical_smiles")
                .agg(
                    pIC50=("pIC50_record", "min"),
                    n_measurements=("pIC50_record", "size"),
                    pIC50_min=("pIC50_record", "min"),
                    pIC50_max=("pIC50_record", "max"),
                    inchikey14=("inchikey14", "first"),
                )
            )
            cens_grouped["pIC50_range"] = cens_grouped["pIC50_max"] - cens_grouped["pIC50_min"]
            # Only keep censored compounds not already present in exact measurements (by canonical SMILES or InChIKey14)
            known_ik = set(exact_df["inchikey14"])
            cens_grouped = cens_grouped[
                (~cens_grouped.index.isin(exact_df.index)) & (~cens_grouped["inchikey14"].isin(known_ik))
            ].copy()
            cens_grouped["source_label"] = "censored_inactive"
            cens_grouped["y"] = 0
            rescued_inactives = cens_grouped

    funnel["rescued_censored_inactives"] = int(len(rescued_inactives))

    combined = pd.concat([exact_actives, exact_inactives, rescued_inactives], axis=0)
    # Final safety deduplication by InChIKey14 (keep exact measurement over censored if any collision)
    combined = combined.sort_values(["source_label", "n_measurements"], ascending=[False, False])
    combined = combined[~combined["inchikey14"].duplicated(keep="first")].sort_index()
    funnel["final_classification_compounds"] = int(len(combined))
    funnel["final_actives"] = int((combined["y"] == 1).sum())
    funnel["final_inactives"] = int((combined["y"] == 0).sum())
    return combined, funnel


def filter_approved_small_molecules(
    df_drugs: pd.DataFrame,
    min_mw: float = 100.0,
    max_mw: float = 900.0,
    min_heavy_atoms: int = 6,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Filter raw approved drug records to unique, valid, drug-like organic small molecules."""
    funnel: dict[str, int] = {"raw_approved_records": int(len(df_drugs))}
    df = df_drugs.copy()

    df["canonical_smiles"] = [standardize_smiles(s) for s in df["smiles"]]
    df = df[df["canonical_smiles"].notna()].copy()
    funnel["after_valid_smiles"] = int(len(df))

    mols = [Chem.MolFromSmiles(s) for s in df["canonical_smiles"]]
    mws = []
    heavy = []
    is_org = []
    for m in mols:
        if m is None:
            mws.append(np.nan)
            heavy.append(0)
            is_org.append(False)
            continue
        mws.append(Descriptors.MolWt(m))
        heavy.append(m.GetNumHeavyAtoms())
        atoms = {a.GetAtomicNum() for a in m.GetAtoms()}
        is_org.append(6 in atoms and atoms.issubset(_ORGANIC_ATOMIC_NUMS))

    df["MolWt"] = mws
    df["heavy_atoms"] = heavy
    df["is_organic"] = is_org

    df = df[df["is_organic"]].copy()
    funnel["after_organic_filter"] = int(len(df))

    df = df[(df["MolWt"] >= min_mw) & (df["MolWt"] <= max_mw) & (df["heavy_atoms"] >= min_heavy_atoms)].copy()
    funnel["after_mw_and_size_filter"] = int(len(df))

    df["inchikey14"] = [smiles_to_inchikey14(s) for s in df["canonical_smiles"]]
    df = df[df["inchikey14"].notna()].copy()

    # Deduplicate by connectivity layer (InChIKey14), preserving any positive-control flag if merged later
    df = df.drop_duplicates(subset=["inchikey14"], keep="first").reset_index(drop=True)
    funnel["unique_approved_small_molecules"] = int(len(df))
    return df, funnel


def remove_leakage(
    train_df: pd.DataFrame,
    approved_df: pd.DataFrame,
    train_fps: list,
    approved_fps: list,
    tanimoto_cutoff: float | None = None,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Remove approved drugs (exact InChIKey14 / canonical SMILES) and optionally close analogues from training data."""
    approved_ik = set(approved_df["inchikey14"].dropna())
    approved_smi = set(approved_df["canonical_smiles"].dropna())

    exact_mask = train_df["inchikey14"].isin(approved_ik) | train_df.index.isin(approved_smi)
    analogue_mask = np.zeros(len(train_df), dtype=bool)

    if tanimoto_cutoff is not None and len(approved_fps) > 0:
        for i, fp in enumerate(train_fps):
            sims = DataStructs.BulkTanimotoSimilarity(fp, approved_fps)
            if sims and max(sims) >= tanimoto_cutoff:
                analogue_mask[i] = True

    drop_mask = exact_mask.to_numpy() | analogue_mask
    clean_df = train_df.iloc[~drop_mask].copy()
    stats = {
        "n_before": int(len(train_df)),
        "n_exact_overlap": int(exact_mask.sum()),
        "n_analogue_overlap": int((analogue_mask & ~exact_mask.to_numpy()).sum()),
        "n_removed_total": int(drop_mask.sum()),
        "n_after": int(len(clean_df)),
    }
    return clean_df, stats
