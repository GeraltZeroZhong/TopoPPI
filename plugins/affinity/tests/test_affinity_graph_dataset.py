"""Input alignment, source identity, cropping and complete surface ablation."""
import json

import numpy as np
import pandas as pd
import pytest
from torch_geometric.loader import DataLoader

from topoppi_affinity.graph_dataset import GraphDataset
from topoppi_affinity.graph_features import FEATURE_NAMES, extract_complex_graph, save_graph


@pytest.fixture()
def inputs(tmp_path):
    from pathlib import Path

    graph = extract_complex_graph(Path(__file__).parent / "fixtures" / "tiny_complex.pdb", "A", "B", max_surface_points=16)
    # Isolate an actual model-input issue: a surface-only edge outside the CA
    # radius must disappear in the control, while cross-side edges remain.
    graph["edge_index"] = np.array([[0, 1, 0, 2], [1, 0, 2, 0]], dtype=np.int64)
    graph["edge_attr"] = np.array([[12, 12, 1, 0, 0, 0, 5, 3], [12, 12, 1, 0, 0, 0, 5, 3],
                                    [4, 4, 0, 1, 0, 3, 0, 0], [4, 4, 0, 1, 0, 3, 0, 0]], dtype=np.float32)
    path = tmp_path / "graph.npz"
    save_graph(graph, path)
    manifest = pd.DataFrame([{"source": source, "complex_id": "tiny:A:B", "graph_path": str(path), "status": "ok"}
                             for source in ["PDBbind", "SKEMPI"]])
    sequence_dir = tmp_path / "sequences"
    (sequence_dir / "embeddings").mkdir(parents=True)
    wt_rows = []
    for source, offset in [("PDBbind", 0), ("SKEMPI", 100)]:
        for chain in ["A", "B"]:
            indices = np.flatnonzero(graph["chain"] == chain)[::-1]
            identifier = f"{source}_{chain}"
            wt_rows.append({"source": source, "complex_id": "tiny:A:B", "chain": chain,
                            "sequence_id": identifier, "residue_keys": json.dumps(graph["residue_key"][indices].tolist())})
            values = np.repeat((indices + offset)[:, None], 1152, axis=1).astype(np.float16)
            np.savez_compressed(sequence_dir / "embeddings" / f"{identifier}.npz", embedding=values)
    pd.DataFrame(wt_rows).to_csv(sequence_dir / "wild_type_chains.csv", index=False)
    pd.DataFrame([{"sample_id": "mutant", "complex_id": "tiny:A:B", "chain": "A", "sequence_id": "mutant_A"}]).to_csv(
        sequence_dir / "mutant_chains.csv", index=False)
    mutant = np.repeat(np.array([11, 10])[:, None], 1152, axis=1).astype(np.float16)
    np.savez_compressed(sequence_dir / "embeddings" / "mutant_A.npz", embedding=mutant)
    frame = pd.DataFrame([{"source": "PDBbind", "complex_id": "tiny:A:B", "sample_id": "mutant", "mutations": "GA1W", "ddg": 1.2},
                          {"source": "SKEMPI", "complex_id": "tiny:A:B", "sample_id": "other", "mutations": np.nan, "ddg": 0.}])
    return frame, manifest, sequence_dir, graph


def test_exact_residue_mapping_mutant_chain_reuse_and_source_separation(inputs):
    frame, manifest, sequence_dir, graph = inputs
    dataset = GraphDataset(frame, manifest, sequence_dir)
    first, second = dataset[0], dataset[1]
    np.testing.assert_array_equal(first.esm[:, 0], [0, 1, 2, 3])
    np.testing.assert_array_equal(first.esm_mut[:, 0], [10, 11, 2, 3])
    np.testing.assert_array_equal(second.esm[:, 0], [100, 101, 102, 103])
    assert first.mutation_mask.tolist() == [True, False, False, False]
    assert first.x.shape == first.x_mut.shape == (4, len(FEATURE_NAMES))
    np.testing.assert_array_equal(first.node_weight, graph["x"][:, FEATURE_NAMES.index("surface.footprint_area")])
    batch = next(iter(DataLoader(dataset, batch_size=2)))
    assert batch.x.shape[0] == 8
    assert batch.edge_index.max().item() == 6
    assert batch.sample_id == ["mutant", "other"]
    assert batch.source == ["PDBbind", "SKEMPI"]
    assert batch.selected_nodes.tolist() == [0, 1, 2, 3, 0, 1, 2, 3]


def test_surface_and_plm_controls_cannot_retain_surface_edges_or_weights(inputs):
    frame, manifest, sequence_dir, graph = inputs
    control = GraphDataset(frame, manifest, sequence_dir, use_surface=False, use_plm=False)[0]
    assert control.x[:, graph["surface_feature_mask"]].count_nonzero() == 0
    assert control.x_mut[:, graph["surface_feature_mask"]].count_nonzero() == 0
    assert control.edge_index.shape[1] == 2
    assert control.edge_attr[:, -2:].count_nonzero() == 0
    assert control.node_weight.count_nonzero() == 0
    assert control.esm.shape == (4, 1152)
    assert control.esm.count_nonzero() == control.esm_mut.count_nonzero() == 0


def test_cropping_preserves_distant_mutations_and_reindexes_edges(inputs):
    frame, manifest, sequence_dir, _ = inputs
    frame.loc[0, "mutations"] = "GA1W,AB2D"
    dataset = GraphDataset(frame.iloc[:1], manifest, sequence_dir, use_plm=False, max_nodes=2)
    graph = dataset[0]
    assert graph.num_nodes == 2
    assert graph.selected_nodes.tolist() == [0, 3]
    assert graph.mutation_mask.all()
    assert graph.edge_index.shape == (2, 0)
    np.testing.assert_array_equal(graph.x[:, :21].sum(axis=1), [1, 1])
