"""Screen all Slide 13 candidate targets + CDK2 & CDK9 via live ChEMBL, Open Targets, and RCSB PDB APIs."""

from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from bioinfo_repurposing.data_fetch import (
    CANDIDATE_TARGETS,
    fetch_chembl_mechanisms_approved,
    fetch_opentargets_approved_drugs,
    fetch_pdb_holo_count,
    get_chembl_activity_count,
    get_chembl_status,
)


def screen_one(t: dict[str, str]) -> tuple[dict, str, list[dict]]:
    sym = t["symbol"]
    cid = t["chembl_id"]
    ens = t["ensembl"]
    up = t["uniprot"]

    n_ic50 = get_chembl_activity_count(cid, "IC50")
    n_ki = get_chembl_activity_count(cid, "Ki")

    mech_df = fetch_chembl_mechanisms_approved(cid)
    mech_p4 = (
        mech_df[(mech_df["max_phase"] == 4) & (mech_df["molecule_type"] == "Small molecule")]
        if len(mech_df)
        else mech_df
    )
    ot_p4 = fetch_opentargets_approved_drugs(ens)
    ot_p4_sm = ot_p4[ot_p4["drug_type"] == "Small molecule"] if len(ot_p4) else ot_p4

    chembl_ids = set(mech_p4["molecule_chembl_id"]) if len(mech_p4) else set()
    ot_ids = set(ot_p4_sm["molecule_chembl_id"]) if len(ot_p4_sm) else set()
    union_ids = sorted(chembl_ids | ot_ids)

    names_map = {}
    for _, r in mech_p4.iterrows():
        names_map[r["molecule_chembl_id"]] = r["pref_name"]
    for _, r in ot_p4_sm.iterrows():
        names_map.setdefault(r["molecule_chembl_id"], r["pref_name"])

    n_pdb, top_pdbs = fetch_pdb_holo_count(up)

    pos_list = [{"molecule_chembl_id": mid, "pref_name": names_map.get(mid)} for mid in union_ids]
    row = {
        "symbol": sym,
        "source_list": t["source_list"],
        "chembl_id": cid,
        "uniprot": up,
        "ensembl": ens,
        "chembl_ic50_records": n_ic50,
        "chembl_ki_records": n_ki,
        "chembl_mech_phase4_sm": len(chembl_ids),
        "opentargets_phase4_sm": len(ot_ids),
        "approved_pos_controls_union": len(union_ids),
        "pdb_holo_structures": n_pdb,
        "top_pdb_ids": ", ".join(top_pdbs[:3]),
        "example_approved_drugs": ", ".join([str(names_map.get(m, m)) for m in union_ids[:6]]),
    }
    print(
        f"{sym:8s} ({t['source_list']:13s}) | IC50={n_ic50:6d} | Ki={n_ki:6d} | "
        f"ChEMBL_P4={len(chembl_ids):2d} | OT_P4={len(ot_ids):2d} | Union={len(union_ids):2d} | PDB={n_pdb:4d}",
        flush=True,
    )
    return row, sym, pos_list


def main() -> None:
    out_dir = Path("data/processed")
    out_dir.mkdir(parents=True, exist_ok=True)

    status = get_chembl_status()
    print(f"ChEMBL release: {status.get('chembl_db_version')} (status: {status.get('status')})", flush=True)

    rows = []
    pos_controls_all = {}
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = [ex.submit(screen_one, t) for t in CANDIDATE_TARGETS]
        for fut in as_completed(futs):
            row, sym, pos_list = fut.result()
            rows.append(row)
            pos_controls_all[sym] = pos_list

    df = pd.DataFrame(rows).sort_values(
        ["approved_pos_controls_union", "chembl_ic50_records"], ascending=[False, False]
    )
    df.to_csv(out_dir / "target_screening_summary.csv", index=False)
    with open(out_dir / "target_positive_controls.json", "w", encoding="utf-8") as f:
        json.dump(
            {"chembl_release": status.get("chembl_db_version"), "targets": pos_controls_all},
            f,
            indent=2,
        )
    print("\nSaved data/processed/target_screening_summary.csv and target_positive_controls.json", flush=True)


if __name__ == "__main__":
    main()
