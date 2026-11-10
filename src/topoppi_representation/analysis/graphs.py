"""Matched-node views, physical targets and grouped sample I/O."""

import copy
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from ..chemistry import HYDROPATHY, HYDROPHOBIC, chemical_indices
from ..source import face_areas, native_surface

FEATURE_NAMES = (
    "centered_x_10A",
    "centered_y_10A",
    "centered_z_10A",
    "qP",
    "qA",
    "qC",
    "qX",
    "hydrophobic_fraction",
    "hydropathy_KD",
    "log1p_mass_A2",
    "side",
)
TARGET_NAMES = (
    "O_side0",
    "O_side1",
    "X_PP",
    "X_AA",
    "X_CC",
    "X_PA",
    "X_PC",
    "X_AC",
    "X_PX",
    "X_AX",
    "X_CX",
    "X_XX",
    "X_none",
)
VIEWS = ("self", "native", "ordinary", "aware_uniform")
LEGACY_VIEWS = ("self", "native", "ordinary", "aware")
UNIFORM_PROTOCOL = "canonical_native_uniform_residue_optcuts20_v1"
VIEW_LABELS = {"self": "Self", "native": "N", "ordinary": "Uo", "aware_uniform": "Ua0", "aware": "Uac"}
X_STATES = ("PP", "AA", "CC", "PA", "PC", "AC", "PX", "AX", "CX", "XX", "NO_PROXIMITY")


@dataclass
class GraphProbeSample:
    sample_id: str
    family_id: str
    features: np.ndarray
    mass: np.ndarray
    side: np.ndarray
    edges: dict
    target: np.ndarray
    target_mask: np.ndarray
    metadata: dict = field(default_factory=dict)


@dataclass(frozen=True)
class ReaderConfig:
    width: int = 32
    depth: int = 1
    learning_rate: float = 0.001
    weight_decay: float = 0.0001
    batch_size: int = 8
    max_epochs: int = 200
    patience: int = 20
    seed: int = 17
    cpu_threads: int = 1
    device: str = "cpu"


def common_features(surface):
    """Label-free fixed-unit channels; no data-fitted normalization or IDs."""
    weights = surface.weights
    total = weights.sum()
    center = (surface.xyz * weights[:, None]).sum(axis=0) / total if total > 0 else surface.xyz.mean(axis=0)
    residue_chem = chemical_indices(surface.amino_acids)
    lookup = np.eye(4)[residue_chem]
    denom = np.where(weights > 0, weights, 1.0)
    q = np.asarray(surface.mass @ lookup) / denom[:, None]
    q[weights == 0, 3] = 1.0
    h = np.asarray(surface.mass @ np.array([a in HYDROPHOBIC for a in surface.amino_acids], float)) / denom
    hydro = np.asarray(surface.mass @ np.array([HYDROPATHY.get(a, 0.0) for a in surface.amino_acids])) / denom
    return np.column_stack([(surface.xyz - center) / 10.0, q, h, hydro, np.log1p(weights), surface.side]).astype(
        np.float32
    )


def _uv_patch_edges(patch, layout, k):
    """Union kNN at every cut appearance, using original common node indices.

    Charts are face components joined only across continuous UV native edges.
    A vertex cut more than once in a chart keeps all its UV appearances. Final
    common-node degrees can exceed k after this union; they are reported.
    """
    uv = np.asarray(patch.uv[layout], dtype=np.float64)
    if uv.shape != (len(patch.faces), 3, 2) or not np.all(np.isfinite(uv)):
        raise ValueError(f"Invalid {layout} face-corner UV on side {patch.side}, patch {patch.patch_id}")
    twice_area = (uv[:, 1, 0] - uv[:, 0, 0]) * (uv[:, 2, 1] - uv[:, 0, 1]) - (uv[:, 1, 1] - uv[:, 0, 1]) * (
        uv[:, 2, 0] - uv[:, 0, 0]
    )
    uv_area = np.abs(twice_area).sum() / 2
    if uv_area <= 0:
        raise ValueError("Cannot normalize zero-area UV layout")
    scale = np.sqrt(face_areas(patch).sum() / uv_area)
    uv = uv * scale
    # 1e-9 Å snapping is solely for equal shared endpoint identification.
    tags = np.rint(uv / 1e-9).astype(np.int64)
    parent = np.arange(len(patch.faces))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    edge_faces = {}
    for f, face in enumerate(patch.faces):
        for a, b in ((0, 1), (1, 2), (2, 0)):
            va, vb = int(face[a]), int(face[b])
            ta, tb = tuple(tags[f, a]), tuple(tags[f, b])
            key = (va, vb, ta, tb) if va < vb else (vb, va, tb, ta)
            if key in edge_faces:
                parent[find(f)] = find(edge_faces[key])
            else:
                edge_faces[key] = f
    charts = {}
    for f, face in enumerate(patch.faces):
        chart = charts.setdefault(find(f), {})
        for c, vertex in enumerate(face):
            key = (int(vertex), *tags[f, c].tolist())
            chart.setdefault(key, uv[f, c])
    result = set()
    appearance_count = 0
    for chart in charts.values():
        node_ids = np.array([key[0] for key in chart], dtype=np.int64)
        coords = np.array(list(chart.values()))
        appearance_count += len(coords)
        if len(coords) < 2:
            continue
        tree = cKDTree(coords)
        # Expand only when multiple cut appearances occupy the nearest positions.
        initial = min(len(coords), k + 1)
        _, neighbors = tree.query(coords, k=initial)
        if initial == 1:
            neighbors = neighbors[:, None]
        for i, candidates in enumerate(neighbors):

            def distinct_nodes(indices, node_ids=node_ids, index=i):
                chosen = []
                for j in np.atleast_1d(indices):
                    node = int(node_ids[j])
                    if node != node_ids[index] and node not in chosen:
                        chosen.append(node)
                    if len(chosen) == k:
                        break
                return chosen

            chosen = distinct_nodes(candidates)
            query_count = initial
            while len(chosen) < k and query_count < len(coords):
                query_count = min(len(coords), query_count * 2)
                chosen = distinct_nodes(tree.query(coords[i], k=query_count)[1])
            for node in chosen:
                result.add(tuple(sorted((int(node_ids[i]), node))))
    edges = np.asarray(sorted(result), dtype=np.int64).reshape(-1, 2)
    return edges, {
        "charts": len(charts),
        "appearances": appearance_count,
        "common_nodes": len(patch.xyz),
        "uv_to_A_scale": float(scale),
    }


def uniform_provenance(source_metadata):
    """Require the native-only producer record; contact-aware is never an alias."""
    detail = source_metadata.get("uniform_layout", {})
    expected = {
        "protocol": UNIFORM_PROTOCOL,
        "source_contract": "canonical_native_source_only",
        "residue_weighting": "uniform_one",
        "fragmentation_lambda": 20.0,
        "atomic_coordinates_read": False,
        "atomic_contacts_read": False,
        "previous_uv_used_for_initialization": False,
        "canonical_source_bytes_unchanged": True,
        "legacy_uv_arrays_unchanged": True,
    }
    if any(key not in detail or detail[key] != value for key, value in expected.items()):
        raise ValueError("Ua0 requires verified aware_uniform canonical-native provenance; contact-weighted aware is Uac")
    return copy.deepcopy(detail)


def prepare_matched_views(source, k=8, aware_layout="aware", include_ordinary=True):
    """Prepare common native nodes and four views; missing layouts stay missing."""
    if k < 1:
        raise ValueError("k must be positive")
    surface = native_surface(source)
    edges = {"self": np.empty((0, 2), dtype=np.int64), "native": surface.edges.copy()}
    if aware_layout not in ("aware", "aware_uniform"):
        raise ValueError("aware_layout must be aware (Uac) or aware_uniform (Ua0)")
    if aware_layout == "aware_uniform":
        uniform_provenance(source.metadata)
    metadata = {
        "uv_operator": "chart-local appearance kNN, symmetrized and deduplicated",
        "uv_k_per_appearance": k,
        "native_operator": "native mesh edges",
        "cross_partner_edges": False,
        "missing_views": [],
        "uv_patches": {},
        "aware_layout": aware_layout,
        "view_labels": VIEW_LABELS,
    }
    for layout in ("ordinary", aware_layout) if include_ordinary else (aware_layout,):
        if any(layout not in patch.uv for patch in source.patches):
            metadata["missing_views"].append(layout)
            continue
        parts, diagnostics, offset = [], [], 0
        for patch in source.patches:
            local, detail = _uv_patch_edges(patch, layout, k)
            parts.append(local + offset)
            diagnostics.append({"side": patch.side, "patch": patch.patch_id, **detail})
            offset += len(patch.xyz)
        edges[layout] = np.concatenate(parts) if parts else np.empty((0, 2), dtype=np.int64)
        metadata["uv_patches"][layout] = diagnostics
    metadata["degrees"] = {}
    for view, pairs in edges.items():
        degree = np.bincount(pairs.ravel(), minlength=len(surface.xyz))
        metadata["degrees"][view] = {
            "edges": len(pairs),
            "mean": float(degree.mean()),
            "max": int(degree.max(initial=0)),
        }
    return surface, edges, metadata


def make_probe_sample(source, target, target_mask=None, k=8, aware_layout="aware", include_ordinary=True):
    """Targets are externally supplied physics references, in TARGET_NAMES order."""
    surface, edges, metadata = prepare_matched_views(
        source, k=k, aware_layout=aware_layout, include_ordinary=include_ordinary
    )
    metadata["source_metadata"] = copy.deepcopy(source.metadata)
    target = np.asarray(target, dtype=np.float32)
    if target.shape != (13,):
        raise ValueError("Expected O_side0, O_side1 and eleven X probabilities")
    mask = np.isfinite(target) if target_mask is None else np.asarray(target_mask, bool) & np.isfinite(target)
    if mask.shape != (13,) or (np.any(mask[2:]) and not np.all(mask[2:])):
        raise ValueError("X distribution must be observed as a whole or wholly masked")
    target = np.where(mask, target, 0.0)
    return GraphProbeSample(
        source.sample_id,
        source.family_id,
        common_features(surface),
        surface.weights.astype(np.float32),
        surface.side.copy(),
        edges,
        target,
        mask,
        metadata,
    )


def targets_from_physics(result):
    """Adapt physics.physical_targets output without recomputing any labels."""
    if tuple(result["X_states"]) != X_STATES:
        raise ValueError("Physical X target state order differs from the reader contract")
    target = np.zeros(13, dtype=np.float32)
    mask = np.zeros(13, dtype=bool)
    for side in (0, 1):
        row = result["O_by_side"][str(side)]
        target[side] = row["L_H"]
        mask[side] = row["applicable"]
    target[2:] = result["X"]
    mask[2:] = result["X_valid"]
    return target, mask


def save_probe_sample(sample, path):
    """Cache one prepared case without recomputing reference targets at fitting."""
    path = Path(path).with_suffix(".npz")
    path.parent.mkdir(parents=True, exist_ok=True)
    arrays = {
        "features": sample.features,
        "mass": sample.mass,
        "side": sample.side,
        "target": sample.target,
        "target_mask": sample.target_mask,
    }
    for view, edges in sample.edges.items():
        arrays[f"edges_{view}"] = edges
    np.savez(path, **arrays)
    meta = {
        "sample_id": sample.sample_id,
        "family_id": sample.family_id,
        "views": list(sample.edges),
        "feature_names": FEATURE_NAMES,
        "target_names": TARGET_NAMES,
        "metadata": sample.metadata,
    }
    path.with_suffix(".json").write_text(json.dumps(meta, indent=2) + "\n")
    return path


def load_probe_sample(path):
    path = Path(path).with_suffix(".npz")
    meta = json.loads(path.with_suffix(".json").read_text())
    if tuple(meta["feature_names"]) != FEATURE_NAMES or tuple(meta["target_names"]) != TARGET_NAMES:
        raise ValueError("Prepared probe feature/target contract differs from current reader")
    with np.load(path, allow_pickle=False) as data:
        return GraphProbeSample(
            meta["sample_id"],
            meta["family_id"],
            data["features"],
            data["mass"],
            data["side"],
            {v: data[f"edges_{v}"] for v in meta["views"]},
            data["target"],
            data["target_mask"],
            meta["metadata"],
        )


def grouped_splits(samples, n_splits=5, excluded_test_groups=(), seed=2718, validation_fraction=0.2):
    """Fixed outer grouped folds and train-only inner validation groups.

    Explicitly excluded groups remain eligible for training and validation.
    No targets are accessed. Identical returned splits are reused by all views,
    depths and optimization seeds.
    """
    groups = np.array([s.family_id for s in samples])
    if any(not g for g in groups):
        raise ValueError("A known relationship group is required for every sample")
    excluded = set(excluded_test_groups)
    eligible = np.array(sorted(set(groups) - excluded))
    if len(eligible) < n_splits or n_splits < 2:
        raise ValueError("Need at least n_splits eligible relationship groups")
    rng = np.random.default_rng(seed)
    eligible = rng.permutation(eligible)
    result = []
    for fold, test_groups in enumerate(np.array_split(eligible, n_splits)):
        train_groups = np.array(sorted(set(groups) - set(test_groups)))
        if len(train_groups) < 2:
            raise ValueError("Outer training needs at least two groups for inner validation")
        train_groups = np.random.default_rng(seed + fold + 1).permutation(train_groups)
        n_validation = min(len(train_groups) - 1, max(1, int(np.ceil(len(train_groups) * validation_fraction))))
        validation_groups, fitting_groups = train_groups[:n_validation], train_groups[n_validation:]
        result.append(
            {
                "fold": fold,
                "train": np.flatnonzero(np.isin(groups, fitting_groups)).tolist(),
                "validation": np.flatnonzero(np.isin(groups, validation_groups)).tolist(),
                "test": np.flatnonzero(np.isin(groups, test_groups)).tolist(),
                "test_groups": test_groups.tolist(),
                "excluded_test_groups": sorted(excluded),
            }
        )
    return result
