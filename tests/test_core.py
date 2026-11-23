"""Unit tests verifying Assignment 5 compatibility for curation, scaffold splitting, and metrics."""

import numpy as np
import pandas as pd

from bioinfo_repurposing.curation import curate_activity, largest_fragment_smiles
from bioinfo_repurposing.features import murcko_scaffold, scaffold_split
from bioinfo_repurposing.metrics import bedroc_score, enrichment_factor


def test_curate_activity_matches_assignment5() -> None:
    raw = pd.DataFrame(
        {
            "compound_id": ["C1", "C1", "C2", "C3", "C4", "C5", "C6", "C1"],
            "smiles": ["c1ccccc1O.Cl", "c1ccccc1O", "CCN", "CCO", "CCC", "c1ccccc1N", "CCQ(bad", "c1ccccc1O"],
            "standard_type": ["IC50"] * 6 + ["IC50", "Ki"],
            "standard_relation": ["=", "=", ">", "=", "=", "=", "=", "="],
            "standard_value": [100.0, 0.3, 50.0, 1.0, 25.0, 2.0, 10.0, 1.0],
            "standard_units": ["nM", "uM", "nM", "uM", "ug.mL-1", "M", "nM", "nM"],
        }
    )
    got = curate_activity(raw)
    assert len(got) == 3
    assert got.loc[largest_fragment_smiles("c1ccccc1O"), "n_measurements"] == 2
    assert np.isclose(got.loc[largest_fragment_smiles("CCO"), "pIC50"], 6.0)


def test_scaffold_split_no_shared_scaffolds() -> None:
    smi = ["c1ccccc1C", "c1ccccc1CC", "c1ccccc1CCC", "c1ccncc1C", "c1ccncc1CC", "C1CCCCC1C", "CCCCO", "CCCO"]
    tr, te = scaffold_split(smi, 0.25)
    assert sorted(tr.tolist() + te.tolist()) == list(range(len(smi)))
    s_tr = {murcko_scaffold(smi[i]) for i in tr}
    s_te = {murcko_scaffold(smi[i]) for i in te}
    assert not (s_tr & s_te)


def test_enrichment_factor_and_bedroc() -> None:
    sc = np.array([-9.0, -8.5, -8.0, -7.5, -7.0, -6.5, -6.0, -5.5, -5.0, -4.5])
    act = np.array([1, 1, 0, 0, 1, 0, 0, 0, 0, 0])
    assert np.isclose(enrichment_factor(sc, act, 0.2, higher_is_better=False), 10.0 / 3.0)
    assert np.isclose(enrichment_factor(-sc, act, 0.2, higher_is_better=True), 10.0 / 3.0)
    b = bedroc_score(-sc, act, alpha=20.0, higher_is_better=True)
    assert 0.5 < b <= 1.0


if __name__ == "__main__":
    test_curate_activity_matches_assignment5()
    test_scaffold_split_no_shared_scaffolds()
    test_enrichment_factor_and_bedroc()
    print("[OK] All Assignment 5 compatibility tests passed!")
