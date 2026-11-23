"""Molecular featurization, scaffold splitting, and Tanimoto similarity.

Reuses Assignment 5 (CAP5) functions (`murcko_scaffold`, `scaffold_split`, `random_split`, `ecfp4`, `descriptors`).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Crippen, Descriptors, QED, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

_FPG = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def murcko_scaffold(smiles: str) -> str | None:
    """Bemis-Murcko scaffold of `smiles` as canonical SMILES (Assignment 5)."""
    m = Chem.MolFromSmiles(smiles)
    return None if m is None else MurckoScaffold.MurckoScaffoldSmiles(mol=m)


def scaffold_split(smiles: list[str] | pd.Index, test_fraction: float = 0.25) -> tuple[np.ndarray, np.ndarray]:
    """Deterministic Bemis-Murcko scaffold split (directly reused from Assignment 5 Task B).

    Groups molecules by Bemis-Murcko scaffold, orders groups largest-first (ties by scaffold SMILES),
    fills training until `round((1 - test_fraction) * n)` molecules, and puts remaining groups in test.
    """
    groups: dict[str, list[int]] = {}
    for i, s in enumerate(smiles):
        scaf = murcko_scaffold(s) or ""
        groups.setdefault(scaf, []).append(i)

    sorted_groups = sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    n_train_target = int(round((1.0 - test_fraction) * len(smiles)))

    train_idx: list[int] = []
    test_idx: list[int] = []
    for _, idx_list in sorted_groups:
        if len(train_idx) < n_train_target:
            train_idx.extend(idx_list)
        else:
            test_idx.extend(idx_list)

    return np.array(sorted(train_idx), dtype=int), np.array(sorted(test_idx), dtype=int)


def random_split(n: int, test_fraction: float = 0.25, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Random train/test split for comparison with scaffold split (Assignment 5)."""
    idx = np.random.default_rng(seed).permutation(n)
    cut = int(round((1.0 - test_fraction) * n))
    return np.sort(idx[:cut]), np.sort(idx[cut:])


def ecfp4(smiles: list[str] | pd.Index) -> np.ndarray:
    """ECFP4 (Morgan radius 2, 2048 bits) as a 0/1 float matrix (Assignment 5)."""
    return np.array([_FPG.GetFingerprintAsNumPy(Chem.MolFromSmiles(s)) for s in smiles], dtype=np.float32)


def ecfp4_bitvects(smiles: list[str] | pd.Index) -> list:
    """ECFP4 as RDKit ExplicitBitVect list for fast BulkTanimotoSimilarity."""
    return [_FPG.GetFingerprint(Chem.MolFromSmiles(s)) for s in smiles]


def descriptors(smiles: list[str] | pd.Index) -> pd.DataFrame:
    """Physicochemical RDKit descriptors matching Assignment 5."""
    rows = []
    for s in smiles:
        m = Chem.MolFromSmiles(s)
        rows.append(
            {
                "MolWt": Descriptors.MolWt(m),
                "MolLogP": Crippen.MolLogP(m),
                "TPSA": rdMolDescriptors.CalcTPSA(m),
                "NumHDonors": rdMolDescriptors.CalcNumHBD(m),
                "NumHAcceptors": rdMolDescriptors.CalcNumHBA(m),
                "NumRotatableBonds": rdMolDescriptors.CalcNumRotatableBonds(m),
                "HeavyAtoms": m.GetNumHeavyAtoms(),
                "AromaticRings": rdMolDescriptors.CalcNumAromaticRings(m),
                "FractionCSP3": rdMolDescriptors.CalcFractionCSP3(m),
                "RingCount": rdMolDescriptors.CalcNumRings(m),
                "QED": QED.qed(m),
                "LabuteASA": rdMolDescriptors.CalcLabuteASA(m),
            }
        )
    return pd.DataFrame(rows, index=list(smiles))


def max_tanimoto_similarity(query_fps: list, ref_fps: list) -> np.ndarray:
    """Maximum Tanimoto similarity of each query fingerprint to a reference set (baseline + Applicability Domain)."""
    if not ref_fps:
        return np.zeros(len(query_fps), dtype=float)
    return np.array([max(DataStructs.BulkTanimotoSimilarity(q, ref_fps)) for q in query_fps], dtype=float)
