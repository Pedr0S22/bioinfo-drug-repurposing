"""Download and cache raw ChEMBL bioactivities, ChEMBL approved small molecules, and positive controls."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from bioinfo_repurposing.data_fetch import (
    CANDIDATE_TARGETS,
    fetch_chembl_activities,
    fetch_chembl_approved_drugs,
    fetch_chembl_mechanisms_approved,
    fetch_opentargets_approved_drugs,
    get_chembl_status,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", required=True, help="Target symbol from CANDIDATE_TARGETS (e.g., DRD2, ABL1, AChE)")
    args = parser.parse_args()

    tmap = {t["symbol"].upper(): t for t in CANDIDATE_TARGETS}
    key = args.target.upper()
    if key not in tmap:
        raise ValueError(f"Unknown target {args.target}. Choose from: {sorted(tmap)}")
    target = tmap[key]
    sym = target["symbol"]

    raw_dir = Path("data/raw")
    raw_dir.mkdir(parents=True, exist_ok=True)

    status = get_chembl_status()
    print(f"ChEMBL release: {status.get('chembl_db_version')}")

    act_path = raw_dir / f"{sym.lower()}_chembl_activities_raw.csv"
    if act_path.exists():
        print(f"Cached bioactivities found: {act_path}")
        df_act = pd.read_csv(act_path)
    else:
        print(f"Downloading ChEMBL binding IC50 bioactivities for {sym} ({target['chembl_id']})...")
        df_act = fetch_chembl_activities(target["chembl_id"], standard_types=("IC50",), assay_type="B")
        df_act.to_csv(act_path, index=False)
        print(f"Saved {len(df_act):,} records to {act_path}")

    drugs_path = raw_dir / "chembl_approved_small_molecules_raw.csv"
    if drugs_path.exists():
        print(f"Cached approved drugs found: {drugs_path}")
        df_drugs = pd.read_csv(drugs_path)
    else:
        print("Downloading ChEMBL approved small-molecule drugs (max_phase == 4)...")
        df_drugs = fetch_chembl_approved_drugs()
        df_drugs.to_csv(drugs_path, index=False)
        print(f"Saved {len(df_drugs):,} approved drug records to {drugs_path}")

    print(f"Fetching Phase 4 positive controls for {sym} from ChEMBL mechanism + Open Targets...")
    mech = fetch_chembl_mechanisms_approved(target["chembl_id"])
    mech_p4 = mech[(mech["max_phase"] == 4) & (mech["molecule_type"] == "Small molecule")].copy() if len(mech) else pd.DataFrame()
    ot_p4 = fetch_opentargets_approved_drugs(target["ensembl"])
    ot_p4_sm = ot_p4[ot_p4["drug_type"] == "Small molecule"].copy() if len(ot_p4) else pd.DataFrame()

    pos_rows = {}
    if len(mech_p4):
        for _, r in mech_p4.iterrows():
            pos_rows[r["molecule_chembl_id"]] = {
                "molecule_chembl_id": r["molecule_chembl_id"],
                "pref_name": r["pref_name"],
                "mechanism_of_action": r["mechanism_of_action"],
                "in_chembl_mechanism": True,
                "in_opentargets": False,
            }
    if len(ot_p4_sm):
        for _, r in ot_p4_sm.iterrows():
            mid = r["molecule_chembl_id"]
            if mid in pos_rows:
                pos_rows[mid]["in_opentargets"] = True
            else:
                pos_rows[mid] = {
                    "molecule_chembl_id": mid,
                    "pref_name": r["pref_name"],
                    "mechanism_of_action": r["mechanism_of_action"],
                    "in_chembl_mechanism": False,
                    "in_opentargets": True,
                }

    df_pos = pd.DataFrame(list(pos_rows.values())).sort_values("pref_name", na_position="last")
    pos_path = raw_dir / f"{sym.lower()}_positive_controls_raw.csv"
    df_pos.to_csv(pos_path, index=False)
    print(f"Saved {len(df_pos)} Phase 4 small-molecule positive controls to {pos_path}")

    meta = {
        "target": target,
        "chembl_db_version": status.get("chembl_db_version"),
        "n_raw_activities": int(len(df_act)),
        "n_raw_approved_drugs": int(len(df_drugs)),
        "n_raw_positive_controls": int(len(df_pos)),
    }
    with open(raw_dir / f"{sym.lower()}_metadata.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)


if __name__ == "__main__":
    main()
