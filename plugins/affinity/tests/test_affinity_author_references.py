"""The same PDB can encode different interfaces; never silently substitute one."""

import numpy as np
import pandas as pd
import pytest

from topoppi_affinity.author_references import (
    match_author_partners,
    read_author_partners,
    run_bindpred_features,
    run_proaffinity_graphs,
)


def test_author_graph_partner_assignment_preserves_multichain_groups(tmp_path):
    index = tmp_path / "index.txt"
    index.write_text("1abc\tL, H; I;\n")
    groups = read_author_partners(index)["1abc"]
    match_author_partners({"chains_a": "I", "chains_b": "H,L"}, groups)
    with pytest.raises(ValueError, match="partners differ"):
        match_author_partners({"chains_a": "H", "chains_b": "I"}, groups)


def test_reference_failure_remains_in_requested_cohort(tmp_path):
    pytest.importorskip("torch_geometric")
    repository = tmp_path / "author"
    (repository / "data").mkdir(parents=True)
    (repository / "data/index_example.txt").write_text("1abc\tA; B;\n")
    manifest = tmp_path / "manifest.csv"
    pd.DataFrame([
        {"sample_id": "wrong-interface", "pdb_id": "1abc", "chains_a": "A", "chains_b": "C"},
        {"sample_id": "missing-graph", "pdb_id": "1abc", "chains_a": "A", "chains_b": "B"},
    ]).to_csv(manifest, index=False)
    results = run_proaffinity_graphs(manifest, tmp_path / "result.csv", repository, device="cpu")
    assert results.sample_id.tolist() == ["wrong-interface", "missing-graph"]
    assert results.status.tolist() == ["error", "error"]
    assert results.predicted.isna().all()
    assert "partners differ" in results.iloc[0].error
    assert "FileNotFoundError" in results.iloc[1].error
    assert set(results.training_overlap) == {"possible"}


def test_bindpred_negates_log_kd_and_aligns_author_feature_rows(tmp_path):
    catboost = pytest.importorskip('catboost')
    rng = np.random.default_rng(4)
    matrix = np.zeros((12, 2560))
    matrix[:, :3] = rng.normal(size=(12, 3))
    model = catboost.CatBoostRegressor(iterations=3, depth=2, thread_count=1, verbose=False, allow_writing_files=False)
    model.fit(matrix, -8.+matrix[:, 0])
    model.save_model(str(tmp_path/'model.cbm'))
    np.savez(tmp_path/'features.npz', X=matrix, sample_ids=np.array([f's{i}' for i in range(12)]),
             feature_names=np.array([f'author.{i}' for i in range(2560)]))
    pd.DataFrame({'sample_id': ['missing', 's2', 's0'], 'pkd': 'unreadable_label'}).to_csv(tmp_path/'manifest.csv', index=False)
    result = run_bindpred_features(tmp_path/'manifest.csv', tmp_path/'features.npz',
                                   tmp_path/'model.cbm', tmp_path/'predictions.csv')
    assert result.status.tolist() == ['error', 'ok', 'ok']
    np.testing.assert_allclose(result.predicted.iloc[1:], -model.predict(matrix[[2, 0]]))
