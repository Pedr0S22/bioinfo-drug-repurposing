"""Generate clean, reproducible Checkpoint 1 Jupyter notebooks using nbformat."""

from __future__ import annotations

import argparse
from pathlib import Path

import nbformat as nbf


def build_nb01(target_symbol: str, target_name: str, chembl_id: str, uniprot_id: str, ensembl_id: str) -> nbf.NotebookNode:
    nb = nbf.v4.new_notebook()
    nb.metadata["kernelspec"] = {"display_name": "Python 3 (.venv)", "language": "python", "name": "python3"}

    cells = [
        nbf.v4.new_markdown_cell(
            f"""# Checkpoint 1 · Notebook 01: Target Screening & Dataset Check

**Project 7 — Repurposing Approved Drugs Against One Target**

- **Selected Target:** **{target_symbol}** ({target_name})
- **Target IDs:** ChEMBL `{chembl_id}` · UniProt `{uniprot_id}` · Ensembl `{ensembl_id}`
- **Research Question (one sentence):** *For human **{target_symbol}** (`{chembl_id}`), which approved small-molecule drugs do ligand-based QSAR models rank as candidate binders—and does the leak-free ranking recover the drugs already approved against {target_symbol}?*

> **Reproducibility:** Managed with `uv` (`pyproject.toml` + `uv.lock`). All curation and splitting functions reuse our verified Assignment 5 (CAP5) implementations in `src/bioinfo_repurposing/`.
"""
        ),
        nbf.v4.new_code_cell(
            f"""import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from bioinfo_repurposing.curation import (
    curate_activity,
    curate_activity_classification,
    filter_approved_small_molecules,
    remove_leakage,
    replicate_disagreement_report,
)
from bioinfo_repurposing.features import ecfp4_bitvects, murcko_scaffold

sns.set_theme(style="whitegrid")
plt.rcParams.update({{"figure.dpi": 110, "axes.spines.top": False, "axes.spines.right": False}})

ROOT = Path("..").resolve()
RAW_DIR = ROOT / "data" / "raw"
PROC_DIR = ROOT / "data" / "processed"
FIG_DIR = ROOT / "results" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

TARGET = "{target_symbol}"
TARGET_LOW = TARGET.lower()
"""
        ),
        nbf.v4.new_markdown_cell(
            """## 1. Multi-Target Screening (All Slide 13 Targets + CDK2 & CDK9)

Before committing to one target, we queried **ChEMBL** (live REST API), **Open Targets Platform** (GraphQL `drugAndClinicalCandidates`), and **RCSB PDB** across all Slide 13 targets plus our group's initial candidates (**CDK2** and **CDK9**).

Key verification result:
- **CDK2 (`CHEMBL301`) and CDK9 (`CHEMBL3116`) have 0 FDA-approved (Phase 4 / APPROVAL) drugs** with primary approved mechanism of action in both ChEMBL and Open Targets (all CDK2/CDK9 agents such as dinaciclib, alvocidib, and seliciclib are Phase I–III investigational; approved CDK drugs target CDK4/6).
- Because Project 7's central validation test requires recovering drugs *already approved against the target*, we select the Slide 13 target offering the strongest balance of **rich curated IC50 bioactivity data**, **many approved positive controls**, and **abundant PDB co-crystal structures**.
"""
        ),
        nbf.v4.new_code_cell(
            """screen_df = pd.read_csv(PROC_DIR / "target_screening_summary.csv")
display(
    screen_df[
        [
            "symbol",
            "source_list",
            "chembl_id",
            "chembl_ic50_records",
            "chembl_ki_records",
            "chembl_mech_phase4_sm",
            "opentargets_phase4_sm",
            "approved_pos_controls_union",
            "pdb_holo_structures",
            "top_pdb_ids",
        ]
    ]
)

cdk_check = screen_df[screen_df["symbol"].isin(["CDK2", "CDK9"])]
print("\\n--- Explicit Verification: CDK2 & CDK9 Approved Drug Counts ---")
display(cdk_check[["symbol", "chembl_id", "ensembl", "chembl_mech_phase4_sm", "opentargets_phase4_sm", "approved_pos_controls_union"]])
"""
        ),
        nbf.v4.new_markdown_cell(
            f"""## 2. Raw Bioactivity Curation & Dataset Check for `{target_symbol}`

We apply Assignment 5 Task A's strict curation rules (`curate_activity`) for **exact-value `pIC50` regression** and extend them (`curate_activity_classification`) for **binary classification**:
1. Remove flagged dubious records (`data_validity_comment` not null) and restrict to biochemical/binding assays (`assay_type == 'B'`) and `IC50` only (avoiding unvalidated mixing of `IC50` and `Ki`).
2. Convert valid molar units (`nM`, `uM`, `M`) to `nM` and drop non-positive values.
3. Standardize molecules with RDKit (`largest_fragment_smiles` + charge neutralization) and compute the 14-character **InChIKey connectivity layer** (`inchikey14`) for leak-proof cross-database matching.
4. Aggregate replicate measurements by **median** `pIC50` and flag compounds with replicate range `> 2.0` log units (Assignment 5 Q1b bookkeeping/unit-error diagnostic).
5. For **classification**: active `pIC50 >= 6.0` (IC50 <= 1 uM), inactive `pIC50 <= 5.0` (IC50 >= 10 uM), dropping the `5–6` grey zone and **rescuing unambiguous censored inactives** (`standard_relation` in `{{'>', '>='}}` with `IC50 >= 10,000 nM`).
"""
        ),
        nbf.v4.new_code_cell(
            """raw_act = pd.read_csv(RAW_DIR / f"{TARGET_LOW}_chembl_activities_raw.csv", low_memory=False)
with open(RAW_DIR / f"{TARGET_LOW}_metadata.json", encoding="utf-8") as f:
    meta = json.load(f)

print(f"ChEMBL version: {meta['chembl_db_version']} | Raw IC50/Ki records: {len(raw_act):,}")
display(
    pd.concat(
        {col: raw_act[col].value_counts(dropna=False).head(6) for col in ["standard_type", "standard_relation", "standard_units", "assay_type"]},
        names=["field", "value"],
    ).to_frame("records")
)

df_reg, reg_funnel = curate_activity(raw_act, standard_type="IC50", return_funnel=True)
df_clf, clf_funnel = curate_activity_classification(
    raw_act, active_threshold=6.0, inactive_threshold=5.0, rescue_censored_inactives=True, standard_type="IC50"
)

funnel_df = pd.DataFrame(
    list(clf_funnel.items()), columns=["curation_step", "count"]
)
display(funnel_df)

disagree = replicate_disagreement_report(raw_act, threshold=2.0)
print(f"\\nExact IC50 compounds measured >1 time: {(df_reg['n_measurements'] > 1).sum():,}")
print(f"Compounds whose replicates disagree by > 2.0 log units: {len(disagree):,}")
if len(disagree):
    display(disagree[["pIC50", "n_measurements", "pIC50_min", "pIC50_max", "pIC50_range"]].head(8).round(2))
"""
        ),
        nbf.v4.new_markdown_cell(
            """## 3. Approved Small-Molecule Library & Positive-Control Matching (Leakage Check)

We curate the **100% open-access ChEMBL `max_phase == 4` small-molecule library** identically (salt stripping, neutralization, organic elements only, `100 <= MW <= 900 Da`, `heavy_atoms >= 6`, deduplicated by `inchikey14`).

Then we cross-match our positive controls (union of ChEMBL `mechanism` Phase 4 and Open Targets `APPROVAL / PHASE_4`) onto the curated approved library and measure **exact overlap** and **close-analogue overlap (`Tanimoto >= 0.85`)** with the ChEMBL training data.
"""
        ),
        nbf.v4.new_code_cell(
            """raw_drugs = pd.read_csv(RAW_DIR / "chembl_approved_small_molecules_raw.csv", low_memory=False)
raw_pos = pd.read_csv(RAW_DIR / f"{TARGET_LOW}_positive_controls_raw.csv")

approved_df, drug_funnel = filter_approved_small_molecules(raw_drugs)
display(pd.DataFrame(list(drug_funnel.items()), columns=["approved_library_step", "count"]))

# Match positive controls via both molecule_chembl_id and standardized InChIKey14
pos_ids = set(raw_pos["molecule_chembl_id"].dropna())
pos_raw_rows = raw_drugs[raw_drugs["molecule_chembl_id"].isin(pos_ids)].copy()
pos_std, _ = filter_approved_small_molecules(pos_raw_rows)
pos_iks = set(pos_std["inchikey14"].dropna())

approved_df["is_positive_control"] = (
    approved_df["molecule_chembl_id"].isin(pos_ids) | approved_df["inchikey14"].isin(pos_iks)
).astype(int)

print(f"Curated approved small-molecule library: {len(approved_df):,} unique drugs")
print(f"Positive controls matched in curated library: {approved_df['is_positive_control'].sum()} / {len(raw_pos)}")
display(approved_df[approved_df["is_positive_control"] == 1][["molecule_chembl_id", "pref_name", "first_approval", "MolWt", "inchikey14"]].head(25))

# Compute ECFP4 fingerprints for leakage audit
reg_fps = ecfp4_bitvects(df_reg.index)
clf_fps = ecfp4_bitvects(df_clf.index)
app_fps = ecfp4_bitvects(approved_df["canonical_smiles"])

df_reg_exact_clean, reg_leak_exact = remove_leakage(df_reg, approved_df, reg_fps, app_fps, tanimoto_cutoff=None)
df_reg_strict_clean, reg_leak_strict = remove_leakage(df_reg, approved_df, reg_fps, app_fps, tanimoto_cutoff=0.85)

df_clf_exact_clean, clf_leak_exact = remove_leakage(df_clf, approved_df, clf_fps, app_fps, tanimoto_cutoff=None)
df_clf_strict_clean, clf_leak_strict = remove_leakage(df_clf, approved_df, clf_fps, app_fps, tanimoto_cutoff=0.85)

# Check how many positive controls themselves appear in ChEMBL training!
pos_only_df = approved_df[approved_df["is_positive_control"] == 1]
pos_in_reg = int(df_reg["inchikey14"].isin(set(pos_only_df["inchikey14"])).sum())
pos_in_clf = int(df_clf["inchikey14"].isin(set(pos_only_df["inchikey14"])).sum())
print(f"\\nPositive controls present in raw ChEMBL regression set: {pos_in_reg} / {len(pos_only_df)}")
print(f"Positive controls present in raw ChEMBL classification set: {pos_in_clf} / {len(pos_only_df)}")

leak_summary = pd.DataFrame(
    [
        {"dataset": "Regression (exact IC50)", "tier": "Exact leak-free", **reg_leak_exact},
        {"dataset": "Regression (exact IC50)", "tier": "Strict (Tanimoto<0.85)", **reg_leak_strict},
        {"dataset": "Classification (<=5 vs >=6 + censored)", "tier": "Exact leak-free", **clf_leak_exact},
        {"dataset": "Classification (<=5 vs >=6 + censored)", "tier": "Strict (Tanimoto<0.85)", **clf_leak_strict},
    ]
)
display(leak_summary)
"""
        ),
        nbf.v4.new_markdown_cell(
            """## 4. Chemical Diversity (Bemis–Murcko Scaffolds), Plots & Saving Processed Datasets

We inspect the `pIC50` distribution, class balance, and Bemis–Murcko scaffold diversity (Assignment 5 Part 2), then save all versioned datasets to `data/processed/` for Notebook 02.
"""
        ),
        nbf.v4.new_code_cell(
            """scaffs_clf = pd.Series([murcko_scaffold(s) for s in df_clf_exact_clean.index], index=df_clf_exact_clean.index)
scaff_sizes = scaffs_clf.value_counts()
print(
    f"Leak-free classification set: {len(df_clf_exact_clean):,} compounds · "
    f"{scaffs_clf.nunique():,} distinct Bemis-Murcko scaffolds · "
    f"largest series: {scaff_sizes.iloc[0]} · singletons: {(scaff_sizes == 1).sum():,}"
)

fig, axes = plt.subplots(1, 3, figsize=(14, 4))
axes[0].hist(df_reg["pIC50"], bins=35, color="#2a78d6", alpha=0.85)
axes[0].axvline(5.0, color="#eb6834", ls="--", lw=1.2, label="Inactive <= 5.0")
axes[0].axvline(6.0, color="#1baf7a", ls="--", lw=1.2, label="Active >= 6.0")
axes[0].set(title=f"{TARGET} Curated Exact pIC50 (n={len(df_reg):,})", xlabel="pIC50", ylabel="Compounds")
axes[0].legend(frameon=False, fontsize=8)

cls_counts = df_clf_exact_clean.groupby(["y", "source_label"]).size().unstack(fill_value=0)
cls_counts.index = ["Inactive (0)", "Active (1)"]
cls_counts.plot(kind="bar", stacked=True, ax=axes[1], color=["#eb6834", "#2a78d6"])
axes[1].set(title=f"Leak-Free Classification Balance (n={len(df_clf_exact_clean):,})", ylabel="Compounds")
axes[1].tick_params(axis="x", rotation=0)

axes[2].bar(range(1, min(21, len(scaff_sizes) + 1)), scaff_sizes.iloc[:20].to_numpy(), color="#2a78d6")
axes[2].set(title="Top 20 Bemis-Murcko Scaffold Series", xlabel="Scaffold rank", ylabel="Compounds")

plt.tight_layout()
fig.savefig(FIG_DIR / f"{TARGET_LOW}_dataset_check.png", dpi=150, bbox_inches="tight")
plt.show()

# Save clean processed files
approved_df.to_csv(PROC_DIR / f"{TARGET_LOW}_approved_library.csv", index=False)
df_reg.to_csv(PROC_DIR / f"{TARGET_LOW}_reg_naive.csv")
df_reg_exact_clean.to_csv(PROC_DIR / f"{TARGET_LOW}_reg_leakfree.csv")
df_reg_strict_clean.to_csv(PROC_DIR / f"{TARGET_LOW}_reg_strict.csv")
df_clf.to_csv(PROC_DIR / f"{TARGET_LOW}_clf_naive.csv")
df_clf_exact_clean.to_csv(PROC_DIR / f"{TARGET_LOW}_clf_leakfree.csv")
df_clf_strict_clean.to_csv(PROC_DIR / f"{TARGET_LOW}_clf_strict.csv")
leak_summary.to_csv(PROC_DIR / f"{TARGET_LOW}_leakage_summary.csv", index=False)
funnel_df.to_csv(PROC_DIR / f"{TARGET_LOW}_curation_funnel.csv", index=False)

print("Saved all processed datasets and figures successfully.")
"""
        ),
    ]
    nb["cells"] = cells
    return nb


def build_nb02(target_symbol: str, target_name: str, chembl_id: str) -> nbf.NotebookNode:
    nb = nbf.v4.new_notebook()
    nb.metadata["kernelspec"] = {"display_name": "Python 3 (.venv)", "language": "python", "name": "python3"}

    cells = [
        nbf.v4.new_markdown_cell(
            f"""# Checkpoint 1 · Notebook 02: Classification & Regression Baselines + Leakage Ablation

**Project 7 — Repurposing Approved Drugs Against `{target_symbol}` ({target_name}, `{chembl_id}`)**

In this notebook we establish our **Checkpoint 1 baseline numbers** under both ligand-based formulations:
1. **Binary Classification** (`pIC50 >= 6.0` active vs `<= 5.0` inactive + rescued `> 10 uM` censored inactives):
   - **Max Tanimoto similarity (ECFP4)** to training actives (classic similarity baseline)
   - **Logistic Regression (ECFP4)** (linear baseline)
   - **Random Forest Classifier (ECFP4 + RDKit descriptors)**
2. **Continuous `pIC50` Regression** (exact `=` IC50 records, reusing Assignment 5 Part 2/3):
   - **Random Forest Regressor (ECFP4)**
   - **Random Forest Regressor (RDKit descriptors)**
   - **Random Forest Regressor (ECFP4 + RDKit descriptors)**
3. **Controls & Ablations:**
   - **Random split vs. Bemis–Murcko scaffold split** (`scaffold_split`, Assignment 5 Task B)
   - **Trivial size control:** `MolWt` alone (Assignment 5 Q3b)
   - **3-tier leakage ablation** on the approved-drug ranking:
     - *Naive (leaky)* vs *Exact leak-free* vs *Strict analogue-free (`Tanimoto < 0.85`)*
   - **Applicability Domain:** max ECFP4 Tanimoto similarity of each approved drug to training actives.
"""
        ),
        nbf.v4.new_code_cell(
            f"""from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from bioinfo_repurposing.baselines import (
    evaluate_internal_classification,
    evaluate_internal_regression,
    rank_approved_drugs_suite,
)
from bioinfo_repurposing.features import ecfp4_bitvects, max_tanimoto_similarity

sns.set_theme(style="whitegrid")
plt.rcParams.update({{"figure.dpi": 110, "axes.spines.top": False, "axes.spines.right": False}})

ROOT = Path("..").resolve()
PROC_DIR = ROOT / "data" / "processed"
FIG_DIR = ROOT / "results" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

TARGET = "{target_symbol}"
TARGET_LOW = TARGET.lower()
SEEDS = (0, 1, 2, 3, 4)

approved_df = pd.read_csv(PROC_DIR / f"{{TARGET_LOW}}_approved_library.csv")
df_clf_naive = pd.read_csv(PROC_DIR / f"{{TARGET_LOW}}_clf_naive.csv", index_col=0)
df_clf_clean = pd.read_csv(PROC_DIR / f"{{TARGET_LOW}}_clf_leakfree.csv", index_col=0)
df_clf_strict = pd.read_csv(PROC_DIR / f"{{TARGET_LOW}}_clf_strict.csv", index_col=0)

df_reg_naive = pd.read_csv(PROC_DIR / f"{{TARGET_LOW}}_reg_naive.csv", index_col=0)
df_reg_clean = pd.read_csv(PROC_DIR / f"{{TARGET_LOW}}_reg_leakfree.csv", index_col=0)
df_reg_strict = pd.read_csv(PROC_DIR / f"{{TARGET_LOW}}_reg_strict.csv", index_col=0)

print(f"Approved library: {{len(approved_df):,}} drugs ({{approved_df['is_positive_control'].sum()}} positive controls)")
print(f"Leak-free classification set: {{len(df_clf_clean):,}} | Leak-free regression set: {{len(df_reg_clean):,}}")
"""
        ),
        nbf.v4.new_markdown_cell(
            """## 1. Internal Held-Out Validation: Random Split vs. Bemis–Murcko Scaffold Split

Evaluated on the **exact leak-free** ChEMBL datasets across 5 random seeds (mean ± std).
"""
        ),
        nbf.v4.new_code_cell(
            """internal_clf = evaluate_internal_classification(df_clf_clean, seeds=SEEDS, test_fraction=0.25)
internal_reg = evaluate_internal_regression(df_reg_clean, seeds=SEEDS, test_fraction=0.25)

print("=== Internal Classification (Random vs Bemis-Murcko Scaffold Split) ===")
display(internal_clf.round(3))

print("\\n=== Internal pIC50 Regression (Random vs Bemis-Murcko Scaffold Split) ===")
display(internal_reg.round(3))

internal_clf.to_csv(PROC_DIR / f"{TARGET_LOW}_internal_clf_metrics.csv", index=False)
internal_reg.to_csv(PROC_DIR / f"{TARGET_LOW}_internal_reg_metrics.csv", index=False)
"""
        ),
        nbf.v4.new_markdown_cell(
            """## 2. Approved-Drug Ranking & 3-Tier Leakage Ablation

We rank the entire curated approved small-molecule library and measure recovery of the known `{target_symbol}` approved drugs (`ROC_AUC`, `PR_AUC`, `EF@1%`, `EF@5%`, `BEDROC20`, median rank, and % in top 1% / top 5%) across all three leakage tiers:
1. **Naive (leaky)**: approved drugs remain in ChEMBL training.
2. **Leak-free (exact)**: all approved drugs removed from training by parent `InChIKey14` / canonical SMILES.
3. **Strict (analogue-free)**: exact approved drugs **plus** any training molecule with ECFP4 `Tanimoto >= 0.85` to an approved drug removed.
"""
        ),
        nbf.v4.new_code_cell(
            """res_naive, scores_naive = rank_approved_drugs_suite(df_clf_naive, df_reg_naive, approved_df, "1_Naive_Leaky", seeds=SEEDS)
res_clean, scores_clean = rank_approved_drugs_suite(df_clf_clean, df_reg_clean, approved_df, "2_LeakFree_Exact", seeds=SEEDS)
res_strict, scores_strict = rank_approved_drugs_suite(df_clf_strict, df_reg_strict, approved_df, "3_Strict_Tanimoto<0.85", seeds=SEEDS)

all_rank_res = pd.concat([res_naive, res_clean, res_strict], ignore_index=True)
all_rank_res.to_csv(PROC_DIR / f"{TARGET_LOW}_approved_ranking_metrics.csv", index=False)

cols_show = [
    "leakage_tier",
    "formulation",
    "method",
    "ROC_AUC_mean",
    "PR_AUC_mean",
    "EF@1%_mean",
    "EF@5%_mean",
    "BEDROC20_mean",
    "median_rank_mean",
    "pct_in_top5_mean",
]
print("=== Main Checkpoint 1 Result: Leak-Free (Exact) Approved-Drug Ranking ===")
display(res_clean[cols_show].round(3))

print("\\n=== Full 3-Tier Leakage Ablation ===")
display(all_rank_res[cols_show].round(3))
"""
        ),
        nbf.v4.new_markdown_cell(
            """## 3. Enrichment Curves, Leakage Inflation Plot & Top Ranked Approved Drugs

We visualize (a) how much **data leakage inflates** `EF@5%` and `BEDROC20` across the three tiers, (b) cumulative enrichment curves for the leak-free baselines vs. `MolWt` control and random chance, and (c) inspect both the **recovered positive controls** and the **top non-approved repurposing candidates** alongside their **Applicability Domain** (`AD_max_tanimoto` to leak-free training actives).
"""
        ),
        nbf.v4.new_code_cell(
            """# Compute Applicability Domain (max Tanimoto to leak-free training actives)
app_bv = ecfp4_bitvects(approved_df["canonical_smiles"])
clean_act_bv = ecfp4_bitvects([s for s, y in zip(df_clf_clean.index, df_clf_clean["y"]) if y == 1])
approved_df["AD_max_tanimoto"] = max_tanimoto_similarity(app_bv, clean_act_bv)

# Attach leak-free average scores
for mname, sc in scores_clean.items():
    safe_col = "score_" + mname.lower().replace(" ", "_").replace("(", "").replace(")", "").replace(":", "").replace("+", "_")
    approved_df[safe_col] = sc

fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))

# Panel A: Leakage ablation on EF@5%
plot_methods = [
    "Max Tanimoto (ECFP4)",
    "Logistic Regression (ECFP4)",
    "RF Classifier (ECFP4+Desc)",
    "RF Regressor (ECFP4+Desc)",
    "Control: MolWt alone",
]
sub_ab = all_rank_res[all_rank_res["method"].isin(plot_methods)]
sns.barplot(data=sub_ab, y="method", x="EF@5%_mean", hue="leakage_tier", ax=axes[0], palette=["#eb6834", "#2a78d6", "#1baf7a"])
axes[0].axvline(1.0, color="black", ls="--", lw=1.0, label="Random (EF=1)")
axes[0].set(title=f"{TARGET}: Leakage Ablation on EF@5%", xlabel="Enrichment Factor @ 5%", ylabel="")
axes[0].legend(frameon=False, fontsize=8, loc="lower right")

# Panel B: Enrichment curves (Leak-Free Exact)
is_pos = approved_df["is_positive_control"].to_numpy(dtype=int)
n_tot = len(is_pos)
n_pos = is_pos.sum()
frac_x = np.arange(1, n_tot + 1) / n_tot * 100.0
axes[1].plot([0, 100], [0, 100], color="#8a8984", ls="--", lw=1.0, label="Random")

curve_styles = [
    ("Max Tanimoto (ECFP4)", "#eb6834"),
    ("Logistic Regression (ECFP4)", "#2a78d6"),
    ("RF Classifier (ECFP4+Desc)", "#1baf7a"),
    ("RF Regressor (ECFP4+Desc)", "#8e5bd0"),
    ("Control: MolWt alone", "#555555"),
]
for mname, col in curve_styles:
    order = np.argsort(-scores_clean[mname], kind="stable")
    cum_rec = np.cumsum(is_pos[order]) / n_pos * 100.0
    axes[1].plot(frac_x, cum_rec, color=col, lw=1.8, label=mname)

axes[1].set(
    title=f"{TARGET}: Leak-Free Cumulative Enrichment Curve (n={n_pos} pos controls)",
    xlabel="Top % of Approved Drug Library Screened",
    ylabel="% of Approved Target Drugs Recovered",
    xlim=(0, 30),
    ylim=(0, 105),
)
axes[1].legend(frameon=False, fontsize=8, loc="lower right")

plt.tight_layout()
fig.savefig(FIG_DIR / f"{TARGET_LOW}_checkpoint1_baselines.png", dpi=150, bbox_inches="tight")
plt.show()

# Show positive-control ranks & top 10 non-approved repurposing candidates under Logistic Regression (ECFP4) and RF Regressor (ECFP4+Desc)
approved_df["rank_lr_clf"] = approved_df["score_logistic_regression_ecfp4"].rank(ascending=False, method="min").astype(int)
approved_df["rank_rf_reg"] = approved_df["score_rf_regressor_ecfp4_desc"].rank(ascending=False, method="min").astype(int)
approved_df["rank_tanimoto"] = approved_df["score_max_tanimoto_ecfp4"].rank(ascending=False, method="min").astype(int)

print("\\n=== Top 15 Positive Controls (Sorted by Leak-Free Logistic Regression Rank) ===")
pos_view = approved_df[approved_df["is_positive_control"] == 1].sort_values("rank_lr_clf")
display(pos_view[["pref_name", "molecule_chembl_id", "rank_tanimoto", "rank_lr_clf", "rank_rf_reg", "AD_max_tanimoto"]].head(15).round(3))

print("\\n=== Top 10 Non-Approved Repurposing Hypotheses (Leak-Free Logistic Regression) ===")
cand_view = approved_df[approved_df["is_positive_control"] == 0].sort_values("rank_lr_clf")
display(cand_view[["pref_name", "molecule_chembl_id", "first_approval", "rank_tanimoto", "rank_lr_clf", "rank_rf_reg", "AD_max_tanimoto"]].head(10).round(3))

approved_df.to_csv(PROC_DIR / f"{TARGET_LOW}_approved_ranked_checkpoint1.csv", index=False)
print("Saved full ranked approved library to data/processed/.")
"""
        ),
    ]
    nb["cells"] = cells
    return nb


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--symbol", required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--chembl_id", required=True)
    parser.add_argument("--uniprot_id", required=True)
    parser.add_argument("--ensembl_id", required=True)
    args = parser.parse_args()

    nb_dir = Path("notebooks")
    nb_dir.mkdir(exist_ok=True)

    nb1 = build_nb01(args.symbol, args.name, args.chembl_id, args.uniprot_id, args.ensembl_id)
    nb2 = build_nb02(args.symbol, args.name, args.chembl_id)

    with open(nb_dir / "01_target_screening_and_dataset_check.ipynb", "w", encoding="utf-8") as f:
        nbf.write(nb1, f)
    with open(nb_dir / "02_checkpoint1_baselines.ipynb", "w", encoding="utf-8") as f:
        nbf.write(nb2, f)
    print("Generated notebooks/01_target_screening_and_dataset_check.ipynb and 02_checkpoint1_baselines.ipynb")


if __name__ == "__main__":
    main()
