"""Checkpoint 1 baseline models (classification, regression, and trivial controls)."""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, mean_absolute_error, r2_score, roc_auc_score
from sklearn.preprocessing import StandardScaler

from .features import descriptors, ecfp4, ecfp4_bitvects, max_tanimoto_similarity, random_split, scaffold_split
from .metrics import evaluate_ranking


def evaluate_internal_classification(
    df_clf: pd.DataFrame,
    seeds: tuple[int, ...] = (0, 1, 2, 3, 4),
    test_fraction: float = 0.25,
) -> pd.DataFrame:
    """Evaluate classification baselines on Random split vs Bemis-Murcko Scaffold split."""
    smiles = list(df_clf.index)
    y = df_clf["y"].to_numpy(dtype=int)
    X_fp = ecfp4(smiles)
    fps_bv = ecfp4_bitvects(smiles)
    X_desc = descriptors(smiles).to_numpy(dtype=float)
    X_comb = np.hstack([X_fp, X_desc])

    tr_scaf, te_scaf = scaffold_split(smiles, test_fraction=test_fraction)
    rows = []

    for split_name in ("random", "scaffold"):
        split_runs = [random_split(len(smiles), test_fraction=test_fraction, seed=s) for s in seeds] if split_name == "random" else [(tr_scaf, te_scaf)]
        for model_name in ("Max Tanimoto (ECFP4)", "Logistic Regression (ECFP4)", "RF Classifier (ECFP4+Desc)"):
            rocs, prs = [], []
            for s_idx, (tr, te) in enumerate(split_runs):
                if model_name == "Max Tanimoto (ECFP4)":
                    act_fps = [fps_bv[i] for i in tr if y[i] == 1]
                    te_fps = [fps_bv[i] for i in te]
                    scores = max_tanimoto_similarity(te_fps, act_fps)
                    rocs.append(roc_auc_score(y[te], scores))
                    prs.append(average_precision_score(y[te], scores))
                elif model_name == "Logistic Regression (ECFP4)":
                    seed_list = [seeds[s_idx]] if split_name == "random" else list(seeds)
                    for sd in seed_list:
                        clf = LogisticRegression(C=1.0, max_iter=1000, class_weight="balanced", random_state=sd)
                        clf.fit(X_fp[tr], y[tr])
                        scores = clf.predict_proba(X_fp[te])[:, 1]
                        rocs.append(roc_auc_score(y[te], scores))
                        prs.append(average_precision_score(y[te], scores))
                else:
                    seed_list = [seeds[s_idx]] if split_name == "random" else list(seeds)
                    for sd in seed_list:
                        clf = RandomForestClassifier(
                            n_estimators=150, max_features="sqrt", class_weight="balanced", random_state=sd, n_jobs=-1
                        )
                        clf.fit(X_comb[tr], y[tr])
                        scores = clf.predict_proba(X_comb[te])[:, 1]
                        rocs.append(roc_auc_score(y[te], scores))
                        prs.append(average_precision_score(y[te], scores))

            rows.append(
                {
                    "task": "classification",
                    "model": model_name,
                    "split": split_name,
                    "ROC_AUC_mean": float(np.mean(rocs)),
                    "ROC_AUC_std": float(np.std(rocs)),
                    "PR_AUC_mean": float(np.mean(prs)),
                    "PR_AUC_std": float(np.std(prs)),
                }
            )
    return pd.DataFrame(rows)


def evaluate_internal_regression(
    df_reg: pd.DataFrame,
    seeds: tuple[int, ...] = (0, 1, 2, 3, 4),
    test_fraction: float = 0.25,
) -> pd.DataFrame:
    """Evaluate continuous pIC50 regression baselines on Random vs Bemis-Murcko Scaffold split (Assignment 5 Part 2)."""
    smiles = list(df_reg.index)
    y = df_reg["pIC50"].to_numpy(dtype=float)
    X_fp = ecfp4(smiles)
    X_desc = descriptors(smiles).to_numpy(dtype=float)
    X_comb = np.hstack([X_fp, X_desc])

    tr_scaf, te_scaf = scaffold_split(smiles, test_fraction=test_fraction)
    reps = {
        "RF Regressor (ECFP4)": X_fp,
        "RF Regressor (RDKit Desc)": X_desc,
        "RF Regressor (ECFP4+Desc)": X_comb,
    }
    rows = []
    for split_name in ("random", "scaffold"):
        split_runs = [random_split(len(smiles), test_fraction=test_fraction, seed=s) for s in seeds] if split_name == "random" else [(tr_scaf, te_scaf)]
        for model_name, X in reps.items():
            r2s, maes = [], []
            for s_idx, (tr, te) in enumerate(split_runs):
                seed_list = [seeds[s_idx]] if split_name == "random" else list(seeds)
                for sd in seed_list:
                    mf = "sqrt" if X.shape[1] > 50 else 1.0
                    reg = RandomForestRegressor(n_estimators=150, max_features=mf, random_state=sd, n_jobs=-1)
                    reg.fit(X[tr], y[tr])
                    pred = reg.predict(X[te])
                    r2s.append(r2_score(y[te], pred))
                    maes.append(mean_absolute_error(y[te], pred))
            rows.append(
                {
                    "task": "regression",
                    "model": model_name,
                    "split": split_name,
                    "R2_mean": float(np.mean(r2s)),
                    "R2_std": float(np.std(r2s)),
                    "MAE_mean": float(np.mean(maes)),
                    "MAE_std": float(np.std(maes)),
                }
            )
    return pd.DataFrame(rows)


def rank_approved_drugs_suite(
    df_clf: pd.DataFrame,
    df_reg: pd.DataFrame,
    approved_df: pd.DataFrame,
    tier_label: str,
    seeds: tuple[int, ...] = (0, 1, 2, 3, 4),
) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    """Train all Checkpoint 1 classification & regression baselines and rank the approved-drug library."""
    app_smiles = list(approved_df["canonical_smiles"])
    is_pos = approved_df["is_positive_control"].to_numpy(dtype=int)

    X_app_fp = ecfp4(app_smiles)
    app_bv = ecfp4_bitvects(app_smiles)
    X_app_desc_df = descriptors(app_smiles)
    X_app_desc = X_app_desc_df.to_numpy(dtype=float)
    X_app_comb = np.hstack([X_app_fp, X_app_desc])

    # Classification training features
    clf_smiles = list(df_clf.index)
    y_clf = df_clf["y"].to_numpy(dtype=int)
    X_clf_fp = ecfp4(clf_smiles)
    clf_act_bv = ecfp4_bitvects([s for s, yi in zip(clf_smiles, y_clf) if yi == 1])
    X_clf_desc = descriptors(clf_smiles).to_numpy(dtype=float)
    X_clf_comb = np.hstack([X_clf_fp, X_clf_desc])

    # Regression training features
    reg_smiles = list(df_reg.index)
    y_reg = df_reg["pIC50"].to_numpy(dtype=float)
    X_reg_fp = ecfp4(reg_smiles)
    X_reg_desc = descriptors(reg_smiles).to_numpy(dtype=float)
    X_reg_comb = np.hstack([X_reg_fp, X_reg_desc])

    results = []
    score_dict: dict[str, np.ndarray] = {}

    def _record(method: str, formulation: str, seed_scores: list[np.ndarray]) -> None:
        metrics_per_seed = [evaluate_ranking(sc, is_pos, higher_is_better=True) for sc in seed_scores]
        avg_scores = np.mean(np.vstack(seed_scores), axis=0)
        score_dict[method] = avg_scores
        row = {"leakage_tier": tier_label, "formulation": formulation, "method": method}
        for k in ("ROC_AUC", "PR_AUC", "EF@1%", "EF@5%", "BEDROC20", "median_rank", "mean_rank", "pct_in_top1", "pct_in_top5"):
            vals = [m[k] for m in metrics_per_seed]
            row[f"{k}_mean"] = float(np.mean(vals))
            row[f"{k}_std"] = float(np.std(vals))
        results.append(row)

    # 1. Trivial control: Molecular Weight alone (Assignment 5 Q3b)
    _record("Control: MolWt alone", "control", [X_app_desc_df["MolWt"].to_numpy(dtype=float)])

    # 2. Classification baseline A: Max Tanimoto (ECFP4)
    tan_scores = max_tanimoto_similarity(app_bv, clf_act_bv)
    _record("Max Tanimoto (ECFP4)", "classification", [tan_scores])

    # 3. Classification baseline B: Logistic Regression (ECFP4)
    lr_scores = []
    for sd in seeds:
        lr = LogisticRegression(C=1.0, max_iter=1000, class_weight="balanced", random_state=sd)
        lr.fit(X_clf_fp, y_clf)
        lr_scores.append(lr.predict_proba(X_app_fp)[:, 1])
    _record("Logistic Regression (ECFP4)", "classification", lr_scores)

    # 4. Classification baseline C: RF Classifier (ECFP4 + RDKit Desc)
    rfc_scores = []
    for sd in seeds:
        rfc = RandomForestClassifier(
            n_estimators=150, max_features="sqrt", class_weight="balanced", random_state=sd, n_jobs=-1
        )
        rfc.fit(X_clf_comb, y_clf)
        rfc_scores.append(rfc.predict_proba(X_app_comb)[:, 1])
    _record("RF Classifier (ECFP4+Desc)", "classification", rfc_scores)

    # 5. Regression baseline A: RF Regressor (ECFP4)
    rfr_fp_scores = []
    for sd in seeds:
        reg = RandomForestRegressor(n_estimators=150, max_features="sqrt", random_state=sd, n_jobs=-1)
        reg.fit(X_reg_fp, y_reg)
        rfr_fp_scores.append(reg.predict(X_app_fp))
    _record("RF Regressor (ECFP4)", "regression", rfr_fp_scores)

    # 6. Regression baseline B: RF Regressor (RDKit Desc)
    rfr_desc_scores = []
    for sd in seeds:
        reg = RandomForestRegressor(n_estimators=150, max_features=1.0, random_state=sd, n_jobs=-1)
        reg.fit(X_reg_desc, y_reg)
        rfr_desc_scores.append(reg.predict(X_app_desc))
    _record("RF Regressor (RDKit Desc)", "regression", rfr_desc_scores)

    # 7. Regression baseline C: RF Regressor (ECFP4 + RDKit Desc) — Assignment 5 Part 3 pattern
    rfr_comb_scores = []
    for sd in seeds:
        reg = RandomForestRegressor(n_estimators=150, max_features="sqrt", random_state=sd, n_jobs=-1)
        reg.fit(X_reg_comb, y_reg)
        rfr_comb_scores.append(reg.predict(X_app_comb))
    _record("RF Regressor (ECFP4+Desc)", "regression", rfr_comb_scores)

    return pd.DataFrame(results), score_dict
