"""Ranking and virtual screening evaluation metrics.

Reuses Assignment 5 Task C (`enrichment_factor`) and adds BEDROC, PR-AUC, ROC-AUC, and rank summaries.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


def enrichment_factor(
    scores: np.ndarray | list[float],
    is_active: np.ndarray | list[int] | list[bool],
    fraction: float = 0.01,
    higher_is_better: bool = True,
) -> float:
    """Enrichment factor at top `fraction` of ranked library (reused from Assignment 5 Task C)."""
    s = np.asarray(scores, dtype=float)
    a = np.asarray(is_active, dtype=int)
    n = len(s)
    total_actives = int(a.sum())
    if n == 0 or total_actives == 0:
        return float("nan")
    k = max(1, int(np.ceil(fraction * n)))
    order = np.argsort(-s if higher_is_better else s, kind="stable")
    hits_in_top = int(a[order[:k]].sum())
    return float((hits_in_top / k) / (total_actives / n))


def bedroc_score(
    scores: np.ndarray | list[float],
    is_active: np.ndarray | list[int] | list[bool],
    alpha: float = 20.0,
    higher_is_better: bool = True,
) -> float:
    """Boltzmann-Enhanced Discrimination of ROC (Truchon & Bayly, 2007), bounded in [0, 1]."""
    s = np.asarray(scores, dtype=float)
    a = np.asarray(is_active, dtype=int)
    n = len(s)
    n_act = int(a.sum())
    if n == 0 or n_act == 0 or n_act == n:
        return float("nan")

    order = np.argsort(-s if higher_is_better else s, kind="stable")
    ranks = np.where(a[order] == 1)[0] + 1  # 1-indexed ranks
    ra = n_act / n
    rie_sum = np.sum(np.exp(-alpha * ranks / n))
    rie_min = (1.0 - np.exp(-alpha * ra)) / (ra * (np.exp(alpha / n) - 1.0)) * np.exp(-alpha * (1.0 - ra))
    rie_max = (1.0 - np.exp(-alpha * ra)) / (ra * (np.exp(alpha / n) - 1.0))
    # Standard exact formula via RIE:
    factor = (1.0 - np.exp(-alpha)) / (np.exp(alpha / n) - 1.0)
    rie = rie_sum / (ra * factor)
    rie_max_exact = (1.0 - np.exp(-alpha * ra)) / (ra * (1.0 - np.exp(-alpha)))
    rie_min_exact = (1.0 - np.exp(alpha * ra)) / (ra * (1.0 - np.exp(alpha)))
    denom = rie_max_exact - rie_min_exact
    if denom == 0:
        return float("nan")
    return float(np.clip((rie - rie_min_exact) / denom, 0.0, 1.0))


def positive_control_rank_summary(
    scores: np.ndarray | list[float],
    is_active: np.ndarray | list[int] | list[bool],
    higher_is_better: bool = True,
) -> dict[str, float]:
    """Compute rank statistics (1-indexed) and top-1%/5% recovery for positive controls."""
    s = np.asarray(scores, dtype=float)
    a = np.asarray(is_active, dtype=int)
    n = len(s)
    n_act = int(a.sum())
    if n == 0 or n_act == 0:
        return {
            "n_pos": 0,
            "median_rank": float("nan"),
            "mean_rank": float("nan"),
            "best_rank": float("nan"),
            "worst_rank": float("nan"),
            "pct_in_top1": float("nan"),
            "pct_in_top5": float("nan"),
        }
    order = np.argsort(-s if higher_is_better else s, kind="stable")
    pos_ranks = np.where(a[order] == 1)[0] + 1
    top1_cut = int(np.ceil(0.01 * n))
    top5_cut = int(np.ceil(0.05 * n))
    return {
        "n_pos": n_act,
        "median_rank": float(np.median(pos_ranks)),
        "mean_rank": float(np.mean(pos_ranks)),
        "best_rank": float(np.min(pos_ranks)),
        "worst_rank": float(np.max(pos_ranks)),
        "pct_in_top1": float(100.0 * np.mean(pos_ranks <= top1_cut)),
        "pct_in_top5": float(100.0 * np.mean(pos_ranks <= top5_cut)),
    }


def evaluate_ranking(
    scores: np.ndarray | list[float],
    is_active: np.ndarray | list[int] | list[bool],
    higher_is_better: bool = True,
) -> dict[str, float]:
    """Full virtual-screening evaluation suite for approved-drug ranking."""
    s = np.asarray(scores, dtype=float)
    a = np.asarray(is_active, dtype=int)
    s_dir = s if higher_is_better else -s
    res = {
        "ROC_AUC": float(roc_auc_score(a, s_dir)) if 0 < a.sum() < len(a) else float("nan"),
        "PR_AUC": float(average_precision_score(a, s_dir)) if 0 < a.sum() < len(a) else float("nan"),
        "EF@1%": enrichment_factor(s, a, fraction=0.01, higher_is_better=higher_is_better),
        "EF@5%": enrichment_factor(s, a, fraction=0.05, higher_is_better=higher_is_better),
        "BEDROC20": bedroc_score(s, a, alpha=20.0, higher_is_better=higher_is_better),
    }
    res.update(positive_control_rank_summary(s, a, higher_is_better=higher_is_better))
    return res
