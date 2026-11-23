"""Open-access API clients for ChEMBL, Open Targets Platform, UniProt, and RCSB PDB."""

from __future__ import annotations

import time
from typing import Any

import pandas as pd
import requests

CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
OPENTARGETS_API = "https://api.platform.opentargets.org/api/v4/graphql"
UNIPROT_API = "https://rest.uniprot.org/uniprotkb"
RCSB_SEARCH_API = "https://search.rcsb.org/rcsbsearch/v2/query"

# All Slide 13 targets from Bioinformatics.Projects.pdf + CDK2 & CDK9 (plus CDK4/CDK6 for CDK-family comparison)
CANDIDATE_TARGETS: list[dict[str, str]] = [
    {"symbol": "ABL1", "name": "Tyrosine-protein kinase ABL1", "chembl_id": "CHEMBL1862", "uniprot": "P00519", "ensembl": "ENSG00000097007", "source_list": "Slide 13"},
    {"symbol": "BRAF", "name": "Serine/threonine-protein kinase B-raf", "chembl_id": "CHEMBL5145", "uniprot": "P15056", "ensembl": "ENSG00000157764", "source_list": "Slide 13"},
    {"symbol": "PARP1", "name": "Poly [ADP-ribose] polymerase-1", "chembl_id": "CHEMBL3105", "uniprot": "P09874", "ensembl": "ENSG00000143799", "source_list": "Slide 13"},
    {"symbol": "AChE", "name": "Acetylcholinesterase", "chembl_id": "CHEMBL220", "uniprot": "P22303", "ensembl": "ENSG00000087085", "source_list": "Slide 13"},
    {"symbol": "DHFR", "name": "Dihydrofolate reductase", "chembl_id": "CHEMBL202", "uniprot": "P00374", "ensembl": "ENSG00000228716", "source_list": "Slide 13"},
    {"symbol": "ADORA2A", "name": "Adenosine A2a receptor", "chembl_id": "CHEMBL251", "uniprot": "P29274", "ensembl": "ENSG00000128271", "source_list": "Slide 13"},
    {"symbol": "DRD2", "name": "Dopamine D2 receptor", "chembl_id": "CHEMBL217", "uniprot": "P14416", "ensembl": "ENSG00000149295", "source_list": "Slide 13"},
    {"symbol": "SLC6A4", "name": "Serotonin transporter (SERT)", "chembl_id": "CHEMBL228", "uniprot": "P31645", "ensembl": "ENSG00000108576", "source_list": "Slide 13"},
    {"symbol": "NTRK1", "name": "High affinity nerve growth factor receptor (TrkA)", "chembl_id": "CHEMBL2815", "uniprot": "P04629", "ensembl": "ENSG00000198400", "source_list": "Slide 13"},
    {"symbol": "NTRK2", "name": "BDNF/NT-3 growth factors receptor (TrkB)", "chembl_id": "CHEMBL4898", "uniprot": "Q16620", "ensembl": "ENSG00000148053", "source_list": "Slide 13"},
    {"symbol": "NTRK3", "name": "NT-3 growth factor receptor (TrkC)", "chembl_id": "CHEMBL5608", "uniprot": "Q16288", "ensembl": "ENSG00000140538", "source_list": "Slide 13"},
    {"symbol": "CDK2", "name": "Cyclin-dependent kinase 2", "chembl_id": "CHEMBL301", "uniprot": "P24941", "ensembl": "ENSG00000123374", "source_list": "Group initial"},
    {"symbol": "CDK9", "name": "Cyclin-dependent kinase 9", "chembl_id": "CHEMBL3116", "uniprot": "P50750", "ensembl": "ENSG00000136807", "source_list": "Group initial"},
]


def _get_json(url: str, params: dict[str, Any] | None = None, timeout: int = 45, retries: int = 4) -> dict[str, Any]:
    headers = {"Accept": "application/json", "User-Agent": "bioinfo-drug-repurposing/0.1"}
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, headers=headers, timeout=timeout)
            r.raise_for_status()
            return r.json()
        except Exception:
            if attempt == retries - 1:
                raise
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError("Unreachable")


def _post_json(url: str, payload: dict[str, Any], timeout: int = 45, retries: int = 4) -> dict[str, Any]:
    headers = {"Accept": "application/json", "Content-Type": "application/json", "User-Agent": "bioinfo-drug-repurposing/0.1"}
    for attempt in range(retries):
        try:
            r = requests.post(url, json=payload, headers=headers, timeout=timeout)
            r.raise_for_status()
            return r.json()
        except Exception:
            if attempt == retries - 1:
                raise
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError("Unreachable")


def get_chembl_status() -> dict[str, Any]:
    """Fetch ChEMBL release version and status metadata (with graceful fallback on transient 500)."""
    try:
        return _get_json(f"{CHEMBL_API}/status.json", timeout=15, retries=3)
    except Exception:
        return {"chembl_db_version": "ChEMBL_37", "status": "FALLBACK"}


def get_chembl_activity_count(target_chembl_id: str, standard_type: str = "IC50") -> int:
    """Fast total count of ChEMBL bioactivity records for a target and standard_type."""
    data = _get_json(
        f"{CHEMBL_API}/activity.json",
        params={"target_chembl_id": target_chembl_id, "standard_type": standard_type, "limit": 1},
    )
    return int(data.get("page_meta", {}).get("total_count", 0))


def fetch_chembl_mechanisms_approved(target_chembl_id: str) -> pd.DataFrame:
    """Fetch drugs with documented ChEMBL mechanism of action on `target_chembl_id` (deduplicated by parent molecule)."""
    rows: list[dict[str, Any]] = []
    offset = 0
    limit = 200
    while True:
        data = _get_json(
            f"{CHEMBL_API}/mechanism.json",
            params={"target_chembl_id": target_chembl_id, "limit": limit, "offset": offset},
            timeout=25,
        )
        mechs = data.get("mechanisms", [])
        rows.extend(mechs)
        if not data.get("page_meta", {}).get("next") or not mechs:
            break
        offset += limit

    if not rows:
        return pd.DataFrame(
            columns=["molecule_chembl_id", "pref_name", "max_phase", "molecule_type", "mechanism_of_action", "smiles"]
        )

    # Deduplicate by parent_molecule_chembl_id, keeping max phase
    by_parent: dict[str, dict[str, Any]] = {}
    for r in rows:
        pid = r.get("parent_molecule_chembl_id") or r.get("molecule_chembl_id")
        if not pid:
            continue
        mp = float(r["max_phase"]) if r.get("max_phase") is not None else -1.0
        if pid not in by_parent or mp > (by_parent[pid]["max_phase"] or -1.0):
            by_parent[pid] = {
                "molecule_chembl_id": pid,
                "max_phase": mp if mp >= 0 else None,
                "mechanism_of_action": r.get("mechanism_of_action"),
            }

    # Resolve pref_name, molecule_type, smiles for Phase 4 entries (fast: use local cache if available)
    mol_info: dict[str, dict[str, Any]] = {}
    from pathlib import Path

    cache_path = Path("data/raw/chembl_approved_small_molecules_raw.csv")
    if cache_path.exists():
        cdf = pd.read_csv(cache_path, low_memory=False)
        for _, cr in cdf.iterrows():
            mol_info[str(cr["molecule_chembl_id"])] = {
                "pref_name": cr.get("pref_name"),
                "molecule_type": cr.get("molecule_type"),
                "smiles": cr.get("smiles"),
            }

    missing_p4 = [pid for pid, v in by_parent.items() if v["max_phase"] == 4.0 and pid not in mol_info]
    for i in range(0, len(missing_p4), 25):
        chunk = missing_p4[i : i + 25]
        mdata = _get_json(
            f"{CHEMBL_API}/molecule.json",
            params={"molecule_chembl_id__in": ",".join(chunk), "limit": 50},
            timeout=25,
        )
        for m in mdata.get("molecules", []):
            structs = m.get("molecule_structures") or {}
            mol_info[m["molecule_chembl_id"]] = {
                "pref_name": m.get("pref_name"),
                "molecule_type": m.get("molecule_type"),
                "smiles": structs.get("canonical_smiles"),
            }

    out = []
    for pid, v in by_parent.items():
        info = mol_info.get(pid, {})
        out.append(
            {
                "molecule_chembl_id": pid,
                "pref_name": info.get("pref_name"),
                "max_phase": v["max_phase"],
                "molecule_type": info.get("molecule_type"),
                "mechanism_of_action": v["mechanism_of_action"],
                "smiles": info.get("smiles"),
            }
        )
    return pd.DataFrame(out)


def fetch_opentargets_approved_drugs(ensembl_id: str) -> pd.DataFrame:
    """Query Open Targets Platform GraphQL API for approved (APPROVAL / PHASE_4) drugs targeting `ensembl_id`."""
    query = """
    query TargetClinicalDrugs($ensemblId: String!) {
      target(ensemblId: $ensemblId) {
        id
        approvedSymbol
        drugAndClinicalCandidates {
          count
          rows {
            id
            maxClinicalStage
            drug {
              id
              name
              drugType
              maximumClinicalStage
              parentMolecule {
                id
                name
              }
            }
          }
        }
      }
    }
    """
    resp = _post_json(OPENTARGETS_API, {"query": query, "variables": {"ensemblId": ensembl_id}})
    target = (resp.get("data") or {}).get("target") or {}
    dcc = target.get("drugAndClinicalCandidates") or {}
    rows = dcc.get("rows") or []

    dedup: dict[str, dict[str, Any]] = {}
    for r in rows:
        stage = r.get("maxClinicalStage")
        drug = r.get("drug") or {}
        if stage not in {"APPROVAL", "PHASE_4"}:
            continue
        parent = drug.get("parentMolecule") or {}
        mol_id = parent.get("id") or drug.get("id")
        mol_name = parent.get("name") or drug.get("name")
        if not mol_id:
            continue
        if mol_id not in dedup:
            dedup[mol_id] = {
                "molecule_chembl_id": mol_id,
                "pref_name": mol_name,
                "drug_type": drug.get("drugType"),
                "stage_for_target": stage,
                "max_clinical_stage": drug.get("maximumClinicalStage"),
                "mechanism_of_action": None,
            }
    if not dedup:
        return pd.DataFrame(
            columns=[
                "molecule_chembl_id",
                "pref_name",
                "drug_type",
                "stage_for_target",
                "max_clinical_stage",
                "mechanism_of_action",
            ]
        )
    return pd.DataFrame(list(dedup.values()))


def fetch_pdb_holo_count(uniprot_id: str) -> tuple[int, list[str]]:
    """Query RCSB PDB for human crystal structures of `uniprot_id` with a non-polymer bound ligand."""
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                        "operator": "exact_match",
                        "value": uniprot_id,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.nonpolymer_entity_count",
                        "operator": "greater",
                        "value": 0,
                    },
                },
            ],
        },
        "request_options": {
            "paginate": {"start": 0, "rows": 5},
            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}],
        },
        "return_type": "entry",
    }
    try:
        r = requests.post(RCSB_SEARCH_API, json=query, timeout=(5, 10))
        if r.status_code == 204:
            return 0, []
        r.raise_for_status()
        data = r.json()
        total = int(data.get("total_count", 0))
        top_ids = [item["identifier"] for item in data.get("result_set", [])[:5]]
        return total, top_ids
    except Exception:
        return -1, []


def fetch_chembl_activities(
    target_chembl_id: str,
    standard_types: tuple[str, ...] = ("IC50",),
    assay_type: str | None = "B",
    page_size: int = 1000,
    max_workers: int = 4,
) -> pd.DataFrame:
    """Download ChEMBL bioactivity records for a target using lean fields (`only=...`) and parallel pagination."""
    from concurrent.futures import ThreadPoolExecutor, as_completed

    only_fields = ",".join(
        [
            "activity_id",
            "molecule_chembl_id",
            "canonical_smiles",
            "standard_type",
            "standard_relation",
            "standard_value",
            "standard_units",
            "pchembl_value",
            "assay_chembl_id",
            "assay_type",
            "target_organism",
            "data_validity_comment",
            "document_chembl_id",
        ]
    )

    base_params: dict[str, Any] = {
        "target_chembl_id": target_chembl_id,
        "standard_type__in": ",".join(standard_types),
        "only": only_fields,
        "limit": page_size,
    }
    if assay_type:
        base_params["assay_type"] = assay_type

    def _parse(acts: list[dict[str, Any]]) -> list[dict[str, Any]]:
        out = []
        for a in acts:
            out.append(
                {
                    "activity_id": a.get("activity_id"),
                    "molecule_chembl_id": a.get("molecule_chembl_id"),
                    "smiles": a.get("canonical_smiles"),
                    "standard_type": a.get("standard_type"),
                    "standard_relation": a.get("standard_relation"),
                    "standard_value": a.get("standard_value"),
                    "standard_units": a.get("standard_units"),
                    "pchembl_value": a.get("pchembl_value"),
                    "assay_chembl_id": a.get("assay_chembl_id"),
                    "assay_type": a.get("assay_type"),
                    "target_organism": a.get("target_organism"),
                    "data_validity_comment": a.get("data_validity_comment"),
                    "document_chembl_id": a.get("document_chembl_id"),
                }
            )
        return out

    first = _get_json(f"{CHEMBL_API}/activity.json", params={**base_params, "offset": 0}, timeout=45, retries=5)
    total = int(first.get("page_meta", {}).get("total_count", 0))
    records = _parse(first.get("activities", []))
    print(f"  page 1 fetched ({len(records)} / {total} records)...", flush=True)

    offsets = list(range(page_size, total, page_size))
    if offsets:
        def _fetch_page(off: int) -> list[dict[str, Any]]:
            d = _get_json(f"{CHEMBL_API}/activity.json", params={**base_params, "offset": off}, timeout=45, retries=5)
            return _parse(d.get("activities", []))

        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            futs = {ex.submit(_fetch_page, off): off for off in offsets}
            for fut in as_completed(futs):
                records.extend(fut.result())
                print(f"  fetched {len(records)} / {total} records...", flush=True)

    return pd.DataFrame(records)


def fetch_chembl_approved_drugs(page_size: int = 1000) -> pd.DataFrame:
    """Download all ChEMBL approved small-molecule drugs (`max_phase == 4`, `molecule_type == 'Small molecule'`)."""
    records: list[dict[str, Any]] = []
    offset = 0
    while True:
        data = _get_json(
            f"{CHEMBL_API}/molecule.json",
            params={
                "max_phase": 4,
                "molecule_type": "Small molecule",
                "limit": page_size,
                "offset": offset,
            },
        )
        mols = data.get("molecules", [])
        for m in mols:
            structs = m.get("molecule_structures") or {}
            records.append(
                {
                    "molecule_chembl_id": m.get("molecule_chembl_id"),
                    "pref_name": m.get("pref_name"),
                    "max_phase": m.get("max_phase"),
                    "first_approval": m.get("first_approval"),
                    "molecule_type": m.get("molecule_type"),
                    "smiles": structs.get("canonical_smiles"),
                    "standard_inchi_key": structs.get("standard_inchi_key"),
                }
            )
        if not data.get("page_meta", {}).get("next") or not mols:
            break
        offset += page_size
    return pd.DataFrame(records)
