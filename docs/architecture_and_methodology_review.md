# Architecture & Methodology Review — Project 7 (Drug Repurposing)

This document audits the group's working guide (`Projeto 7 – Reposicionamento de fármacos aprovados contra um alvo.pdf`) and early prototype (`notebooks/codigo.ipynb`) against the official project rules (`Bioinformatics.Projects.pdf`) and the cheminformatics best practices taught in Assignment 5 (`assignment5_uc2023235452.ipynb`).

---

## 1. What Is Strong in the Original `Projeto 7` Guide

The methodology outlined in `Projeto 7` is well-aligned with the grading rubric:

- **Question-driven design**: Directly operationalizes Slide 13 (*does a ligand-based model recover the drugs already approved against our target?*).
- **Leakage awareness (#1 risk)**: Explicitly requires *leave-approved-out* validation and recommends reporting the **naive (leaky) vs. leak-free** comparison to quantify performance inflation.
- **Honest chemical splitting**: Uses Bemis–Murcko **scaffold splits** instead of random splits.
- **Virtual screening metrics**: Emphasizes `PR-AUC`, `EF@1%`, `EF@5%`, `BEDROC`, and positive-control ranks instead of accuracy alone.
- **Open-access Plan B**: Anticipates DrugBank licensing restrictions by proposing ChEMBL `max_phase = 4` + Open Targets Platform.
- **Clear 4-person role split**: Primary + backup responsibilities (`A/B/C/D`), satisfying the *Group organisation* rubric criterion.

---

## 2. Critical Bugs & Architectural Gaps Identified (and How We Fix Them)

### Gap 1 — Target Selection & Missing Positive Controls for CDK2 / CDK9
- **Issue**: Project 7's core validation rule is: *"if the drugs already approved against that target do not surface near the top, it is not evidence."* Neither **CDK2** (`CHEMBL301`) nor **CDK9** (`CHEMBL3116`) has **any FDA-approved (Phase 4) small-molecule drug** whose primary approved mechanism of action targets CDK2 or CDK9 (approved CDK drugs target **CDK4/CDK6**: palbociclib, ribociclib, abemaciclib, trilaciclib).
- **Bug in prototype `codigo.ipynb`**: 4 of the 6 ChEMBL target IDs in `cdk_targets` were incorrect targets (`CHEMBL302` = Rat ADORA2A, `CHEMBL359` = *Bacillus subtilis*, `CHEMBL303` = Rat GluN2D, `CHEMBL2972` = Human PFKP).
- **Fix**: Screen **all Slide 13 targets + CDK2 & CDK9** programmatically via live ChEMBL and Open Targets APIs, document the exact approved-drug counts (including live verification that CDK2 & CDK9 have 0 approved MoA drugs), and select a single Slide 13 target with both rich curated IC50 data and a strong positive-control set.

### Gap 2 — Overcomplicated Pre-Training Idea vs. Clean Single-Target Baselines
- **Issue**: `codigo.ipynb` attempted multi-target CDK-family pre-training + fine-tuning before even building a single-target baseline. Multi-target pre-training also risks **cross-target data leakage** (e.g., an approved multi-kinase drug appearing in a sister kinase's pre-training set).
- **Fix**: Keep Checkpoint 1 focused on clean, single-target, leak-free baselines (**Max Tanimoto similarity on ECFP4**, **Logistic Regression on ECFP4**, and **Random Forest regression on ECFP4 / RDKit descriptors**). Save richer models (XGBoost, ECFP4 + RDKit descriptors, ChemBERTa-77M-MLM embeddings, Tanimoto-kernel SVM) for Week 2 under the exact same splits and leakage controls.

### Gap 3 — Bioactivity Curation & Censored Inactives (`>` Records)
- **Issue in `codigo.ipynb`**:
  1. Mixed `IC50`, `Ki`, and `Kd` into a single `pIC50` column without justification.
  2. Used `drop_duplicates()` (keeping an arbitrary first record) instead of aggregating replicates by **median** and checking for >2 log-unit replicate discrepancies (Assignment 5 Task A & Q1b).
  3. Omitted **salt stripping & neutralization** (`largest_fragment_smiles` + `Uncharger`), allowing salt/free-base pairs to duplicate or leak across train and approved sets.
  4. Filtering on `pchembl_value is not None` silently discards all censored `>` assays (`IC50 > 10 µM`), creating severe active-class skew.
- **Fix**: Reuse Assignment 5's `curate_activity` (`IC50` only, `assay_type == 'B'`, `data_validity_comment is None`, molar unit conversion, salt stripping + neutralization, median aggregation, replicate-range diagnostics). Run **both**:
  - **Binary classification**: active `pIC50 >= 6.0`, inactive `pIC50 <= 5.0` (dropping the `5–6` grey zone), **rescuing unambiguous censored inactives** (`standard_relation` in `{'>', '>='}` with `IC50 >= 10,000 nM`).
  - **Continuous `pIC50` regression**: exact `=` `IC50` measurements only (retaining the `5–6` region and avoiding arbitrary binary cutoffs), ranking approved drugs by predicted `pIC50`.

### Gap 4 — Cross-Database Leakage Matching (Parent InChIKey14 + 3-Tier Ablation)
- **Issue**: Matching training compounds to approved drugs by raw SMILES or ChEMBL ID misses salt forms, protonation states, and stereoisomers. Also, removing only exact approved drugs leaves close medicinal-chemistry analogues (`Tanimoto >= 0.85`) in the training set.
- **Fix**: Standardize both ChEMBL bioactivities and the approved-drug library identically (largest fragment → neutralize → canonical SMILES → **14-character InChIKey connectivity layer**). Evaluate a **3-tier leakage ablation**:
  1. **Naive (leaky)**: approved drugs left inside training.
  2. **Exact leak-free**: all approved drugs (by InChIKey14 / canonical SMILES) removed from training.
  3. **Strict analogue-free**: exact approved drugs **plus** training molecules with ECFP4 `Tanimoto >= 0.85` to any approved drug removed from training.

### Gap 5 — Approved Drug Library Curation
- **Issue**: Raw approved drug tables contain biologics, peptides, inorganic salts (`Li+`, `NaCl`), organometallics, contrast agents, and tiny solvents that cannot be modeled meaningfully by small-molecule ECFP4 QSAR.
- **Fix**: Use 100% open-access **ChEMBL `max_phase == 4` & `molecule_type == 'Small molecule'`**, then filter to valid organic small molecules (`C` present, organic elements only, `100 <= MW <= 900 Da`, `heavy_atoms >= 6`), standardize, deduplicate by `InChIKey14`, and log every removed entry in the Dataset Check funnel.

### Gap 6 — Trivial Confounder Controls & Applicability Domain (Assignment 5 Q2c & Q3b)
- **Issue**: Without a trivial control, high global ROC-AUC can simply reflect molecular-weight bias; without an Applicability Domain check, it is unclear whether recovered positive controls or top repurposing hits are trivial analogues vs. true scaffold hops.
- **Fix**: Include **Molecular Weight (`MolWt`) alone** and **random ranking** in every evaluation table, and record each approved drug's **maximum Tanimoto similarity to training actives** (`AD_max_tanimoto`) alongside its rank.

### Gap 7 — Reproducible Engineering Architecture
- **Issue**: A single ad-hoc notebook without pinned dependencies, saved datasets, or reusable modules.
- **Fix**:
  - `pyproject.toml` + `uv.lock` managed with **`uv`** (reproducible virtual environment `.venv`).
  - Modular Python package under `src/bioinfo_repurposing/` (`data_fetch.py`, `curation.py`, `features.py`, `metrics.py`) cleanly wrapping Assignment 5's tested functions.
  - Numbered, executed notebooks under `notebooks/` (`01_target_screening_and_dataset_check.ipynb`, `02_checkpoint1_baselines.ipynb`) writing versioned CSVs and figures to `data/raw/`, `data/processed/`, and `results/figures/`.

---

## 3. How to Verify Directly That CDK2 and CDK9 Have 0 FDA-Approved Primary-MoA Drugs

You can double-check (and show your professor directly) using three independent open-access sources:

### A. Browser Verification (2-minute visual check)
1. **Open Targets Platform**:
   - **CDK2 (`ENSG00000123374`)**: `https://platform.opentargets.org/target/ENSG00000123374` → scroll to **Known Drugs** → filter **Phase = 4** → **0 drugs**. (All listed CDK2 inhibitors—dinaciclib, alvocidib/flavopiridol, seliciclib, milciclib, tagtociclib—are Phase I–III investigational.)
   - **CDK9 (`ENSG00000136807`)**: `https://platform.opentargets.org/target/ENSG00000136807` → scroll to **Known Drugs** → filter **Phase = 4** → **0 drugs**.
   - Contrast with **CDK4 (`ENSG00000135446`)** or **ABL1 (`ENSG00000097007`)** → immediately shows Phase 4 FDA-approved drugs (palbociclib/ribociclib/abemaciclib/trilaciclib for CDK4; imatinib/dasatinib/nilotinib/bosutinib/ponatinib/asciminib for ABL1).
2. **ChEMBL Web Interface**:
   - **CDK2 (`CHEMBL301`)**: `https://www.ebi.ac.uk/chembl/target_report_card/CHEMBL301/` → **Drug Mechanisms** tab → check `Max Phase` column (none at Phase 4 as approved CDK2 mechanism).
   - **CDK9 (`CHEMBL3116`)**: `https://www.ebi.ac.uk/chembl/target_report_card/CHEMBL3116/` → **Drug Mechanisms** tab → check `Max Phase` column.
3. **FDA Drugs@FDA / DrugCentral**:
   - Search `CDK2` or `CDK9` at `https://drugcentral.org/` → no FDA-approved drugs with CDK2 or CDK9 as primary approved pharmacological target (approved CDK drugs are strictly CDK4/6 inhibitors).

### B. Programmatic Verification (Reproducible API Query in This Repo)
Run our live screening module (also executed in `notebooks/01_target_screening_and_dataset_check.ipynb`):

```python
from bioinfo_repurposing.data_fetch import (
    fetch_chembl_mechanisms_approved,
    fetch_opentargets_approved_drugs,
)

for sym, cid, ens in [("CDK2", "CHEMBL301", "ENSG00000123374"), ("CDK9", "CHEMBL3116", "ENSG00000136807")]:
    mech = fetch_chembl_mechanisms_approved(cid)
    mech_p4 = mech[mech["max_phase"] == 4] if len(mech) else mech
    ot_p4 = fetch_opentargets_approved_drugs(ens)
    print(f"{sym}: ChEMBL mechanism Phase 4 = {len(mech_p4)}, Open Targets Phase 4 = {len(ot_p4)}")
```
