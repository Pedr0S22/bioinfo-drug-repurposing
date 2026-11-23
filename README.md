# Bioinformatics Project 7 — Repurposing Approved Drugs Against One Target

**Research Question:** *For human acetylcholinesterase (**AChE**, `CHEMBL220`), which FDA/EMA-approved small-molecule drugs do ligand-based QSAR models rank as candidate binders—and does a leak-free ranking recover the approved AChE drugs near the top of the library?*

---

## Quick Setup (Managed with `uv`)

```bash
# 1. Create virtual environment (.venv) and install all locked dependencies + local package
uv sync

# 2. (Optional, Week 2) Install XGBoost / ChemBERTa extras
uv sync --extra week2

# 3. Verify Assignment 5 compatibility tests
uv run python tests/test_core.py

# 4. Reproduce live target screening & AChE downloads (cached in data/raw/)
uv run python scripts/screen_targets.py
uv run python scripts/fetch_target_and_approved.py --target AChE

# 5. Launch JupyterLab inside the uv environment
uv run jupyter lab
```

---

## Repository Architecture

```text
bioinfo-drug-repurposing/
├── pyproject.toml                       # uv project & dependency specification
├── uv.lock                              # Locked reproducible environment
├── src/bioinfo_repurposing/
│   ├── data_fetch.py                    # Open-access ChEMBL, Open Targets (GraphQL v4), RCSB PDB clients
│   ├── curation.py                      # Assignment 5 Task A curation + censored inactive rescue + InChIKey14 leakage control
│   ├── features.py                      # Assignment 5 Task B Bemis–Murcko scaffold_split, ECFP4, RDKit descriptors, Tanimoto AD
│   ├── metrics.py                       # Assignment 5 Task C enrichment_factor + BEDROC20, ROC-AUC, PR-AUC, rank stats
│   └── baselines.py                     # Classification & regression baselines + 3-tier leakage ablation
├── scripts/
│   ├── screen_targets.py                # Live 13-target screening (Slide 13 targets + CDK2 & CDK9 verification)
│   ├── fetch_target_and_approved.py     # Caches target bioactivities, approved small molecules, and positive controls
│   └── build_notebooks.py               # Reproducible notebook generator
├── notebooks/
│   ├── 01_target_screening_and_dataset_check.ipynb  # Executed Checkpoint 1 Dataset Check
│   └── 02_checkpoint1_baselines.ipynb               # Executed Checkpoint 1 Baselines & Leakage Ablation
├── docs/
│   ├── architecture_and_methodology_review.md       # Audit of Projeto 7 methodology, gaps fixed, CDK2/CDK9 proof
│   └── checkpoint1_summary.md                       # Complete Checkpoint 1 deliverable & presentation summary
├── data/
│   ├── raw/                             # Cached open-access ChEMBL 37 & Open Targets downloads
│   └── processed/                       # Curated leak-free datasets, screening table, and baseline metrics
├── results/figures/                     # Saved high-resolution Checkpoint 1 plots
└── tests/test_core.py                   # Unit tests verifying Assignment 5 function equivalence
```

---

## Checkpoint 1 Highlights (Real Data, `AChE` / `CHEMBL220`)

- **Why not CDK2 or CDK9?** Live queries against ChEMBL 37 and Open Targets (`ENSG00000123374`, `ENSG00000136807`) confirm **0 FDA/EMA-approved Phase 4 drugs** with primary MoA on CDK2 or CDK9 (see [docs/architecture_and_methodology_review.md](docs/architecture_and_methodology_review.md)).
- **Curated Data:** `9,882` raw ChEMBL binding `IC50` records → **`5,664`** unique exact-`pIC50` compounds (`5,613` leak-free) and **`4,338`** leak-free classification compounds (`2,095` Bemis–Murcko scaffolds, including `+585` rescued `> 10 µM` censored inactives).
- **Approved Library:** `2,027` curated organic small-molecule approved drugs from ChEMBL (`max_phase = 4`), containing **`18`** matched AChE positive controls (`14/18` of which would have leaked into naive regression training).
- **Leak-Free First Numbers (5 seeds):**
  - **Max Tanimoto (ECFP4):** `ROC-AUC = 0.847` · `PR-AUC = 0.467` · `EF@1% = 48.26` · `EF@5% = 13.25` · `BEDROC20 = 0.654` · Median Rank `21 / 2,027` (`66.7%` in top 5%)
  - **RF Classifier (ECFP4+Desc):** `ROC-AUC = 0.631` · `PR-AUC = 0.247` · `EF@1% = 27.88` · `EF@5% = 9.27` · `BEDROC20 = 0.426` (`46.7%` in top 5%)
  - **RF Regressor (ECFP4+Desc):** `ROC-AUC = 0.613` · `PR-AUC = 0.200` · `EF@1% = 19.30` · `EF@5% = 7.07` · `BEDROC20 = 0.322` (`35.6%` in top 5%)
  - **MolWt Control (size bias check):** `ROC-AUC = 0.378` · `EF@1% = 0.00` · `EF@5% = 0.00`

See [docs/checkpoint1_summary.md](docs/checkpoint1_summary.md) for full tables, error analysis, group role split (`A/B/C/D`), and Plan B.
