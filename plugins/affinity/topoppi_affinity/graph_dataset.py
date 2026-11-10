"""Paired WT/MT residue graphs with frozen, residue-matched PLM embeddings."""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset
from torch_geometric.data import Data

from topoppi_affinity.graph_features import FEATURE_NAMES, load_graph, mutant_node_features

PLM_DIMENSION = 1152


@lru_cache(maxsize=64)
def _cached_graph(path):
    return load_graph(path)


@lru_cache(maxsize=64)
def _cached_embedding(path):
    with np.load(path, allow_pickle=False) as source:
        return source["embedding"]


def _read_frame(value):
    return value.copy() if isinstance(value, pd.DataFrame) else pd.read_csv(value, keep_default_na=False)


class GraphDataset(Dataset):
    """One PyG Data object contains both WT and MT states on identical nodes.

    Source is part of structural identity: PDBbind and SKEMPI can provide the
    same PDB chain groups with different cleaned residue numbering. This is
    independent of the broader complex/family groups used for data splitting.
    Feature scaling is deliberately left to the train-only fitting code.
    """

    def __init__(self, frame, graph_manifest, sequence_dir, use_surface=True,
                 use_plm=True, max_nodes=512, target="ddg"):
        self.frame = _read_frame(frame).fillna({"mutations": ""}).reset_index(drop=True)
        graphs = _read_frame(graph_manifest)
        self.graph_paths = {(str(row.source), str(row.complex_id)): str(row.graph_path)
                            for row in graphs.itertuples() if row.status == "ok"}
        self.sequence_dir = Path(sequence_dir) if sequence_dir is not None else None
        self.use_surface = use_surface
        self.use_plm = use_plm
        self.max_nodes = max_nodes
        self.target = target
        self.wt_chains = {}
        self.mt_chains = {}
        if use_plm:
            wild_type = pd.read_csv(self.sequence_dir / "wild_type_chains.csv", keep_default_na=False)
            mutants = pd.read_csv(self.sequence_dir / "mutant_chains.csv", keep_default_na=False)
            for row in wild_type.itertuples():
                key = (str(row.source), str(row.complex_id), str(row.chain))
                self.wt_chains[key] = (str(row.sequence_id), json.loads(row.residue_keys))
            self.mt_chains = {(str(row.sample_id), str(row.chain)): str(row.sequence_id)
                              for row in mutants.itertuples()}

    def __len__(self):
        return len(self.frame)

    def _nodes(self, graph, mutation_mask):
        count = len(graph["x"])
        if self.max_nodes is None or count <= self.max_nodes:
            return np.arange(count)
        distance = graph["x"][:, FEATURE_NAMES.index("contact.partner_distance")].copy()
        mutations = np.flatnonzero(mutation_mask)
        if len(mutations):
            mutation_distance = np.linalg.norm(graph["pos"][:, None] - graph["pos"][mutations][None], axis=2).min(axis=1)
            distance = np.minimum(distance, mutation_distance)
        remaining = np.flatnonzero(~mutation_mask)
        nearest = remaining[np.argsort(distance[remaining], kind="stable")[:max(self.max_nodes - len(mutations), 0)]]
        return np.sort(np.concatenate([mutations, nearest]))

    def _plm_features(self, row, graph, selected):
        wt = np.zeros((len(selected), PLM_DIMENSION), dtype=np.float32)
        mt = np.zeros_like(wt)
        if not self.use_plm:
            return wt, mt
        for chain in np.unique(graph["chain"][selected]):
            positions = np.flatnonzero(graph["chain"][selected] == chain)
            sequence_id, residue_keys = self.wt_chains[(str(row.source), str(row.complex_id), str(chain))]
            residue_lookup = {key: i for i, key in enumerate(residue_keys)}
            sequence_positions = np.asarray([residue_lookup[str(graph["residue_key"][selected[i]])] for i in positions])
            embedding = _cached_embedding(str(self.sequence_dir / "embeddings" / f"{sequence_id}.npz"))
            if len(embedding) != len(residue_keys):
                raise ValueError(f"{row.sample_id} chain {chain}: embedding length differs from structural residue map")
            wt[positions] = embedding[sequence_positions]
            mutant_id = self.mt_chains.get((str(row.sample_id), str(chain)), sequence_id)
            mutant_embedding = _cached_embedding(str(self.sequence_dir / "embeddings" / f"{mutant_id}.npz"))
            if len(mutant_embedding) != len(residue_keys):
                raise ValueError(f"{row.sample_id} chain {chain}: mutant embedding length differs from WT residue map")
            mt[positions] = mutant_embedding[sequence_positions]
        return wt, mt

    def __getitem__(self, index):
        row = self.frame.iloc[index]
        graph = _cached_graph(self.graph_paths[(str(row.source), str(row.complex_id))])
        x_mut, _, mutation_mask = mutant_node_features(graph, str(row.mutations))
        selected = self._nodes(graph, mutation_mask)
        x = graph["x"][selected].copy()
        x_mut = x_mut[selected]
        old_to_new = np.full(len(graph["x"]), -1, dtype=np.int64)
        old_to_new[selected] = np.arange(len(selected))
        edge_index = old_to_new[graph["edge_index"]]
        keep = (edge_index >= 0).all(axis=0)
        if not self.use_surface:
            attr = graph["edge_attr"]
            keep &= ((attr[:, 2] > 0) & (attr[:, 0] <= 10.)) | (attr[:, 3] > 0) | (attr[:, 4] > 0)
        edge_index = edge_index[:, keep]
        edge_attr = graph["edge_attr"][keep].copy()
        if not self.use_surface:
            x[:, graph["surface_feature_mask"]] = 0.
            x_mut[:, graph["surface_feature_mask"]] = 0.
            edge_attr[:, graph["edge_surface_feature_mask"]] = 0.
        esm, esm_mut = self._plm_features(row, graph, selected)
        partner_distance = graph["x"][selected, FEATURE_NAMES.index("contact.partner_distance")]
        surface_area = x[:, FEATURE_NAMES.index("surface.area")].copy()
        footprint_area = x[:, FEATURE_NAMES.index("surface.footprint_area")].copy()
        return Data(
            x=torch.from_numpy(x), x_mut=torch.from_numpy(x_mut),
            esm=torch.from_numpy(esm), esm_mut=torch.from_numpy(esm_mut),
            mutation_mask=torch.from_numpy(mutation_mask[selected]),
            edge_index=torch.from_numpy(edge_index), edge_attr=torch.from_numpy(edge_attr),
            side=torch.from_numpy(graph["side"][selected]), pos=torch.from_numpy(graph["pos"][selected]),
            contact_mask=torch.from_numpy(partner_distance <= 10.),
            surface_area=torch.from_numpy(surface_area), footprint_area=torch.from_numpy(footprint_area),
            node_weight=torch.from_numpy(footprint_area.copy()),
            y=torch.tensor([float(row[self.target])], dtype=torch.float32),
            sample_id=str(row.sample_id), complex_id=str(row.complex_id), source=str(row.source),
            selected_nodes=torch.from_numpy(selected),
        )
