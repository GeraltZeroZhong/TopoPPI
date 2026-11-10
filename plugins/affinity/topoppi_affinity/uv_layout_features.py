"""Matched spatial readouts on the same native interface, before/after UV cuts.

This representation experiment consumes distances *in* UV charts. It never convolves across
packed islands. All arms share 3-D sampling masses, residue embeddings, channel
projection, spatial scales and readout capacity. No affinity labels enter here.
"""
from __future__ import annotations

import json
import time
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components, dijkstra
from scipy.spatial.distance import cdist
from threadpoolctl import threadpool_limits
from topoppi.atlas.footprints import (
    geometric_contact_partner_map,
    mesh_vertex_residue_labels,
    residue_footprint_pieces,
    source_atom_residue_labels,
)
from topoppi.atlas.seams import uv_seam_topology
from topoppi.atlas.uv import as_corner_uv
from topoppi.config import OptCutsConfig, ParameterizationConfig, SurfaceConfig, TopologyConfig
from topoppi.io.io_loader import PDBLoader
from topoppi.mesh.parameterization import Parameterizer
from topoppi.mesh.provenance import OPTCUTS_GEOMETRY_VERTEX_IDS, SOURCE_FACE_IDS
from topoppi.mesh.surface import SurfaceGenerator
from topoppi.mesh.topology import TopologyManager
from topoppi.optimization.optcuts import OptCutsUVOptimizer

from topoppi_affinity.features import _chains, _footprint_areas, _group_data
from topoppi_affinity.sequence_baselines import SequenceFeatureStore

SPATIAL_VIEWS = ('xyz', 'geodesic', 'uv_harmonic', 'uv_geometry', 'uv_topoppi', 'uv_permuted')
LAYOUTS = ('uv_harmonic', 'uv_geometry', 'uv_topoppi')


def default_protocol(binary):
    return dict(version='uv_local_correlations_v1', optcuts_bin=str(Path(binary).resolve()),
                channels=32, projection_seed=20260925, sampling_seed=42,
                permutation_seed=17, points_per_side=256, minimum_points_per_patch=16,
                radii=[2., 4., 8., 16.], timeout_sec=120., distortion_bound=4.1,
                fragmentation_weights={'uv_geometry': 0., 'uv_topoppi': 20.},
                uv_scale='one isotropic factor per physical patch; total UV area equals 3D area',
                uv_neighbors='same uncut UV chart only; no atlas packing',
                geometry='same complete prepared 3D face domain in every arm',
                features='frozen random ESMC projection; four local symmetric channel correlations',
                targets_read=False)


@lru_cache(maxsize=4)
def projection(dimension, channels, seed):
    # Fixed, label-free random channels; no PCA fitted on a held-out structure.
    return np.random.default_rng(seed).normal(size=(dimension, channels)) / np.sqrt(dimension)


def sample_surface(mesh, count, seed):
    """Equal-mass stratified samples, reused at identical barycentric locations."""
    rng = np.random.default_rng(seed)
    areas = np.asarray(mesh.area_faces)
    total = float(areas.sum())
    if total <= 0 or count < 2:
        raise ValueError('Surface sampling needs positive area and at least two points')
    quantiles = (np.arange(count) + rng.random(count)) * total / count
    faces = np.searchsorted(np.cumsum(areas), quantiles, side='right')
    a, b = rng.random((2, count))
    root = np.sqrt(a)
    barycentric = np.column_stack((1-root, root*(1-b), root*b))
    xyz = np.einsum('ni,nij->nj', barycentric, np.asarray(mesh.triangles)[faces])
    return dict(face_indices=faces, barycentric=barycentric, xyz=xyz,
                mass=np.full(count, total/count))


def chart_ids(mesh, uv):
    topology = uv_seam_topology(mesh, uv)
    pairs = topology.incident_faces[topology.internal_mask & ~topology.seam_mask]
    rows, cols = (pairs[:, 0], pairs[:, 1]) if len(pairs) else ([], [])
    graph = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(mesh.faces), len(mesh.faces))).tocsr()
    return connected_components(graph, directed=False)[1]


def normalized_uv(mesh, uv):
    corners = as_corner_uv(mesh, uv).copy()
    left, right = corners[:, 1]-corners[:, 0], corners[:, 2]-corners[:, 0]
    area = np.abs(left[:, 0]*right[:, 1]-left[:, 1]*right[:, 0]).sum()/2.
    if not np.isfinite(area) or area <= 0:
        raise ValueError('UV area is not finite positive')
    return corners * np.sqrt(float(mesh.area)/area)


def uv_distances(mesh, uv, samples, *, permutation_seed=None):
    corners = normalized_uv(mesh, uv)
    locations = np.einsum('ni,nij->nj', samples['barycentric'], corners[samples['face_indices']])
    charts = chart_ids(mesh, corners)[samples['face_indices']]
    if permutation_seed is not None:
        # Move location AND chart ID together; channels/masses stay fixed.
        # Equal sample masses preserve the full channel distribution in a patch.
        order = np.random.default_rng(permutation_seed).permutation(len(locations))
        locations, charts = locations[order], charts[order]
    distances = cdist(locations, locations)
    distances[charts[:, None] != charts[None, :]] = np.inf
    return distances, locations, charts


def intrinsic_distances(mesh, samples):
    """Pre-cut face-dual shortest paths, with within-face Euclidean distances.

    Diskification copies are stitched by geometry IDs. This is a discretized
    intrinsic-distance control, not exact continuous surface geodesics.
    """
    identities = np.asarray(mesh.metadata.get(OPTCUTS_GEOMETRY_VERTEX_IDS, np.arange(len(mesh.vertices))))
    triangles = identities[np.asarray(mesh.faces)]
    edges = np.sort(triangles[:, ((0, 1), (1, 2), (2, 0))].reshape(-1, 2), axis=1)
    _, inverse, counts = np.unique(edges, axis=0, return_inverse=True, return_counts=True)
    if np.any(counts > 2):
        raise ValueError('Prepared geometric surface is not edge-manifold')
    order = np.argsort(inverse, kind='stable')
    starts = np.cumsum(np.r_[0, counts[:-1]])[counts == 2]
    pairs = np.column_stack((order[starts]//3, order[starts+1]//3))
    centers = np.asarray(mesh.triangles_center)
    lengths = np.linalg.norm(centers[pairs[:, 0]]-centers[pairs[:, 1]], axis=1)
    graph = coo_matrix((np.r_[lengths, lengths],
                        (np.r_[pairs[:, 0], pairs[:, 1]], np.r_[pairs[:, 1], pairs[:, 0]])),
                       shape=(len(centers), len(centers))).tocsr()
    faces, lookup = np.unique(samples['face_indices'], return_inverse=True)
    paths = dijkstra(graph, directed=False, indices=faces)
    distances = paths[:, faces][lookup][:, lookup]
    offset = np.linalg.norm(samples['xyz']-centers[samples['face_indices']], axis=1)
    distances += offset[:, None] + offset[None, :]
    same = samples['face_indices'][:, None] == samples['face_indices'][None, :]
    distances[same] = cdist(samples['xyz'], samples['xyz'])[same]
    return distances


def local_correlations(channels, distances, mass, radii):
    """Area-averaged symmetric channel products after radial local filtering."""
    z, mass = np.asarray(channels, dtype=float), np.asarray(mass, dtype=float)
    if distances.shape != (len(z), len(z)) or len(mass) != len(z) or np.any(mass <= 0):
        raise ValueError('Channels, distances and positive masses must describe the same samples')
    indices = np.triu_indices(z.shape[1])
    values = []
    for radius in radii:
        kernel = np.exp(-.5*(distances/float(radius))**2)
        kernel[distances > 3*float(radius)] = 0.
        weights = kernel * mass[None, :]
        # Self mass guarantees a defined filter even for an isolated sample.
        message = weights @ z / weights.sum(axis=1, keepdims=True)
        correlation = z.T @ (mass[:, None]*message) / mass.sum()
        values.append(((correlation+correlation.T)/2.)[indices])
    result = np.concatenate(values)
    if not np.isfinite(result).all():
        raise ValueError('Nonfinite local spatial readout')
    return result


def spatial_names(protocol):
    upper = list(zip(*np.triu_indices(protocol['channels']), strict=True))
    return np.array([f'spatial.{operation}.r{r:g}.{a}_{b}' for operation in ('sum', 'absdiff')
                     for r in protocol['radii'] for a, b in upper])


def _embedding_map(row, store, groups):
    result = {}
    for chain in sorted(set(groups[0]+groups[1])):
        sequence_id, keys = store.wild[(row['source'], row['complex_id'], chain)]
        _, embedding = store._embedding(sequence_id)
        if len(embedding) != len(keys):
            raise ValueError('Observed residue keys and embedding lengths differ')
        result.update(zip(keys, embedding, strict=True))
    return result


def _mass_check(mesh, uv, source_labels, expected):
    pieces = residue_footprint_pieces(mesh, uv, mesh_vertex_residue_labels(mesh, source_labels))
    actual = {key: sum(piece['footprint_mass'] for piece in value) for key, value in pieces.items()}
    if set(actual) != set(expected):
        raise ValueError('UV footprint reconstruction changed the residue domain')
    differences = np.array([actual[k]-expected[k] for k in expected])
    if not np.allclose([actual[k] for k in expected], list(expected.values()), rtol=1e-9, atol=1e-9):
        raise ValueError('Cut footprint pieces do not reconstruct the original 3D residue masses')
    return dict(max_residue_mass_error=float(np.abs(differences).max(initial=0.)),
                residue_count=len(pieces), piece_count=sum(len(p) for p in pieces.values()),
                fragmented_residue_count=sum(len(p) > 1 for p in pieces.values()))


def extract_case(row, sequence_dir, output, protocol):
    """Generate both surfaces and paired UV arms; no target field is accessed."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    groups = sorted((_chains(row['chains_a']), _chains(row['chains_b'])))
    store = SequenceFeatureStore(sequence_dir)
    embeddings = _embedding_map(row, store, groups)
    dimension = len(next(iter(embeddings.values())))
    transform = projection(dimension, protocol['channels'], protocol['projection_seed'])
    projected = {key: value @ transform for key, value in embeddings.items()}
    loader = PDBLoader(row['structure_path'])
    coordinates = [_group_data(loader, group) for group in groups]
    partner_map = geometric_contact_partner_map(*coordinates[0], *coordinates[1], distance_cutoff=6.)
    reverse = {}
    for left, values in partner_map.items():
        for right in values:
            reverse.setdefault(right, set()).add(left)
    degrees = [{key: len(value) for key, value in partner_map.items()},
               {key: len(value) for key, value in reverse.items()}]
    aggregates = {view: [] for view in SPATIAL_VIEWS}
    failures = {view: [] for view in SPATIAL_VIEWS}
    optimizers = {name: OptCutsUVOptimizer(OptCutsConfig(
        optcuts_bin=protocol['optcuts_bin'], optcuts_mode=100, timeout_sec=protocol['timeout_sec'],
        optcuts_distortion_bound=protocol['distortion_bound'], residue_fragmentation_weight=weight))
        for name, weight in protocol['fragmentation_weights'].items()}
    patch_reports, side_masses, weighted_embeddings = [], [], []
    for side, ((coords, atoms), (other, _)) in enumerate((coordinates, coordinates[::-1])):
        generator = SurfaceGenerator(coords, config=SurfaceConfig())
        mesh = generator.generate_mesh()
        if mesh is None:
            raise ValueError(f'Native surface failed: {generator.last_report}')
        labels = source_atom_residue_labels(atoms)
        manager = TopologyManager(mesh, other, config=TopologyConfig())
        patches = manager.get_interface_patches()
        if not patches:
            raise ValueError(f'No retained interface on partner {side}')
        side_area = sum(float(p.area) for p in patches)
        side_masses.append(side_area)
        weighted = np.zeros(dimension)
        side_vectors = {view: [] for view in SPATIAL_VIEWS}
        for index, patch in enumerate(patches):
            patch_started = time.monotonic()
            mass_before = _footprint_areas(patch, labels)
            for key, area in mass_before.items():
                weighted += area*embeddings[str(key)]
            face_ids = np.asarray(patch.metadata[SOURCE_FACE_IDS]).copy()
            triangles = np.asarray(patch.triangles).copy()
            parameterizer = Parameterizer(ParameterizationConfig())
            prepared = parameterizer.prepare_patch(patch)
            if prepared is None:
                raise ValueError('Common topology preparation failed')
            if (not np.array_equal(face_ids, prepared.metadata[SOURCE_FACE_IDS])
                    or not np.array_equal(triangles, prepared.triangles)):
                raise ValueError('Common topology preparation changed the original face domain')
            count = max(protocol['minimum_points_per_patch'],
                        int(round(protocol['points_per_side']*float(prepared.area)/side_area)))
            samples = sample_surface(prepared, count, protocol['sampling_seed']+index)
            vertex_keys = mesh_vertex_residue_labels(prepared, labels)
            vertex_channels = np.asarray([projected[str(key)] for key in vertex_keys])
            channels = np.einsum('ni,nij->nj', samples['barycentric'],
                                 vertex_channels[np.asarray(prepared.faces)[samples['face_indices']]])
            report = dict(side=side, patch=index, area_3d=float(prepared.area),
                          vertices=len(prepared.vertices), faces=len(prepared.faces), points=count,
                          domain_preserved=True, layouts={})
            arrays = dict(vertices=np.asarray(prepared.vertices), faces=np.asarray(prepared.faces),
                          residue_keys=np.asarray(vertex_keys, dtype=str), source_face_ids=face_ids,
                          geometry_vertex_ids=np.asarray(prepared.metadata[OPTCUTS_GEOMETRY_VERTEX_IDS]),
                          channels=channels, **samples)
            distance_inputs = {'xyz': cdist(samples['xyz'], samples['xyz']),
                               'geodesic': intrinsic_distances(prepared, samples)}
            for layout in LAYOUTS:
                layout_started = time.monotonic()
                try:
                    if layout == 'uv_harmonic':
                        uv = parameterizer.flatten_patch(prepared, method='harmonic')
                        if uv is None:
                            raise ValueError('Harmonic parameterization failed')
                        uv = as_corner_uv(prepared, uv)
                        execution = {}
                    else:
                        current = prepared.copy()
                        optimizer = optimizers[layout]
                        objective = {str(key): 1.+degrees[side].get(str(key), 0) for key in set(labels)}
                        optimizer.optimize_patches([current], pack=False, build_report=False,
                                                   source_residue_labels=labels, residue_weights=objective)
                        uv = as_corner_uv(current, key='uv_optcuts')
                        execution = current.metadata['optcuts_execution']
                    invariant = _mass_check(prepared, uv, labels, mass_before)
                    distances, locations, charts = uv_distances(prepared, uv, samples)
                    distance_inputs[layout] = distances
                    arrays.update({layout: uv, layout+'_points': locations, layout+'_charts': charts})
                    if layout == 'uv_topoppi':
                        distance_inputs['uv_permuted'] = uv_distances(
                            prepared, uv, samples, permutation_seed=protocol['permutation_seed']+index)[0]
                    report['layouts'][layout] = dict(status='ok', seconds=time.monotonic()-layout_started,
                                                      chart_count=int(len(np.unique(chart_ids(prepared, uv)))),
                                                      **invariant, execution=execution)
                except (OSError, RuntimeError, ValueError) as exc:
                    report['layouts'][layout] = dict(status='error', error=str(exc),
                                                      seconds=time.monotonic()-layout_started)
                    failures[layout].append(f'side {side} patch {index}: {exc}')
                    if layout == 'uv_topoppi':
                        failures['uv_permuted'].append(f'side {side} patch {index}: {exc}')
            for view, distances in distance_inputs.items():
                side_vectors[view].append(local_correlations(channels, distances, samples['mass'], protocol['radii'])
                                          * float(prepared.area)/side_area)
            for view in SPATIAL_VIEWS:
                if len(side_vectors[view]) != index+1 and not failures[view]:
                    raise ValueError(f'Unrecorded missing spatial view {view}')
            np.savez_compressed(output/f'side_{side}_patch_{index}.npz', **arrays)
            report['seconds'] = time.monotonic()-patch_started
            patch_reports.append(report)
        weighted_embeddings.append(weighted/max(side_area, 1.))
        for view in SPATIAL_VIEWS:
            aggregates[view].append(np.sum(side_vectors[view], axis=0) if len(side_vectors[view]) == len(patches) else None)
    size = len(spatial_names(protocol))
    vectors = {}
    for view in SPATIAL_VIEWS:
        if failures[view]:
            vectors[view] = np.full(size, np.nan)
        else:
            left, right = aggregates[view]
            vectors[view] = np.concatenate((left+right, np.abs(left-right)))
    w_left, w_right = weighted_embeddings
    native_w = np.concatenate((w_left+w_right, np.abs(w_left-w_right)))
    np.savez_compressed(output/'features.npz', **vectors, native_w=native_w,
                        feature_names=spatial_names(protocol))
    report = dict(sample_id=row['sample_id'], protocol=protocol, patches=patch_reports,
                  side_areas=side_masses, failures=failures,
                  statuses={key: 'error' if value else 'ok' for key, value in failures.items()},
                  seconds=time.monotonic()-started, affinity_labels_read=False)
    (output/'case.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


def extract_case_job(row, sequence_dir, output, protocol):
    """Resume complete cases; scientific input/layout failures keep their rows."""
    output = Path(output)
    if (output/'case.json').exists():
        report = json.loads((output/'case.json').read_text())
        if report['protocol'] != protocol:
            raise ValueError('Existing case uses another protocol; use a separate output root')
        return report
    started = time.monotonic()
    with threadpool_limits(limits=1):
        try:
            return extract_case(row, sequence_dir, output, protocol)
        except (OSError, RuntimeError, ValueError, KeyError) as exc:
            output.mkdir(parents=True, exist_ok=True)
            report = dict(sample_id=row['sample_id'], protocol=protocol, patches=[],
                          statuses={view: 'error' for view in SPATIAL_VIEWS},
                          failures={view: [f'{type(exc).__name__}: {exc}'] for view in SPATIAL_VIEWS},
                          seconds=time.monotonic()-started, affinity_labels_read=False)
            np.savez_compressed(output/'features.npz', **{v: np.full(len(spatial_names(protocol)), np.nan)
                                                         for v in SPATIAL_VIEWS},
                                feature_names=spatial_names(protocol))
            (output/'case.json').write_text(json.dumps(report, indent=2)+'\n')
            return report
