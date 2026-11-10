"""Prevent external-score sign mistakes and held-out selection mistakes."""
import math

import pandas as pd
import pytest

from topoppi_affinity.external_absolute import native_dg_to_pkd, read_prodigy_training_ids, select_records


def test_energy_sign_conventions_agree_for_same_nanomal_kd():
    magnitude = 0.001987 * 298.15 * math.log(1e9)
    assert native_dg_to_pkd(-magnitude) == pytest.approx(9.0)
    assert native_dg_to_pkd(magnitude, positive_affinity_score=True) == pytest.approx(9.0)
    with pytest.raises(ValueError):
        native_dg_to_pkd(float('nan'))


def test_original_publication_pdb_suffix_is_recognized_for_overlap(tmp_path):
    table = tmp_path / 'supp2.csv'
    table.write_text('PDB,Binding_affinity\n1A2K.pdb,9.3\n1ACB.pdb,13.1\n')
    assert read_prodigy_training_ids(table) == {'1a2k', '1acb'}


def test_external_heldout_union_contains_both_tests_without_validation_expansion():
    frame = pd.DataFrame({'sample_id': ['d', 'c', 'b', 'a'],
                          'split_complex': ['train', 'test', 'train', 'validation'],
                          'split_sequence': ['test', 'test', 'validation', 'train']})
    assert select_records(frame, 'heldout').sample_id.tolist() == ['c', 'd']


@pytest.mark.parametrize("subset,expected", [(None, ['a', 'b']), ('heldout', ['b'])])
def test_mirage_supplied_manifest_and_explicit_subset(tmp_path, monkeypatch, subset, expected):
    """Default user manifests need no benchmark split columns."""
    import sys
    from types import SimpleNamespace

    from topoppi_affinity import external_absolute as adapter

    frame = pd.DataFrame({'sample_id': ['b', 'a'], 'pdb_id': ['toy_b', 'toy_a'],
                          'structure_path': ['b.pdb', 'a.pdb'],
                          'chains_a': ['A', 'A'], 'chains_b': ['B', 'B']})
    if subset is not None:
        frame['split_complex'] = ['test', 'train']
        frame['split_sequence'] = ['train', 'train']
    manifest = tmp_path / 'manifest.csv'
    frame.to_csv(manifest, index=False)
    monkeypatch.setitem(sys.modules, 'torch', SimpleNamespace(set_num_threads=lambda n: None))
    monkeypatch.setattr(adapter.pd, 'read_excel', lambda _: pd.DataFrame({'PDB': ['other'], 'name': ['other']}))
    upstream = SimpleNamespace(ensure_feature_files=lambda *args: ('region', 'embedding'),
                               predict_dg=lambda *args, **kwargs: (8.0, {}))
    monkeypatch.setattr(adapter, '_load_mirage', lambda _: upstream)
    monkeypatch.setattr(adapter, 'write_partner_structure', lambda *args: None)
    monkeypatch.setattr(adapter, 'version', lambda _: 'test')
    selected = []
    monkeypatch.setattr(adapter, 'prepare_mirage_embeddings',
                        lambda records, *args, **kwargs: selected.extend(records.sample_id.tolist()))
    for name in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'TRANSFORMERS_OFFLINE'):
        monkeypatch.setenv(name, '1')
    kwargs = {} if subset is None else {'subset': subset}
    output = tmp_path / 'results/predictions.csv'
    result = adapter.run_mirage(manifest, output, tmp_path, device='cpu', **kwargs)
    assert result.sample_id.tolist() == selected == expected
    assert result.status.eq('ok').all()
    assert result.predicted.tolist() == pytest.approx([native_dg_to_pkd(8.0, positive_affinity_score=True)] * len(expected))
    import json
    metadata = json.loads(output.with_suffix('.json').read_text())
    assert metadata['subset'] == (subset or 'all')
    assert metadata['requested'] == metadata['successful'] == len(expected)


def test_mirage_region_preparation_accepts_unsplit_manifest(tmp_path, monkeypatch):
    from concurrent.futures import ThreadPoolExecutor

    from topoppi_affinity import external_absolute as adapter

    manifest = tmp_path / 'manifest.csv'
    pd.DataFrame({'sample_id': ['b', 'a']}).to_csv(manifest, index=False)
    monkeypatch.setattr(adapter, 'ProcessPoolExecutor', ThreadPoolExecutor)
    monkeypatch.setattr(adapter, '_mirage_region', lambda args: {'sample_id': args[0]['sample_id']})
    result = adapter.prepare_mirage_regions(manifest, tmp_path / 'output.csv', tmp_path, jobs=1)
    assert [row['sample_id'] for row in result] == ['a', 'b']


def test_region_remapping_drops_incomplete_backbone_without_shifting_identity():
    from topoppi_affinity.external_absolute import remap_region_positions
    original = pd.DataFrame({'pos': [1, 2, 3, 4], 'res': ['ALA', 'GLY', 'TYR', 'LYS'],
                             'region': ['SUR', 'COR', 'RIM', 'SUP'], 'asa': [9.0, 2.0, 7.0, 4.0]})
    aligned = remap_region_positions(original, [1, 3, 4])
    assert aligned.pos.tolist() == [1, 2, 3]
    assert aligned.res.tolist() == ['ALA', 'TYR', 'LYS']
    assert aligned.region.tolist() == ['SUR', 'RIM', 'SUP']
    assert aligned.asa.tolist() == [9.0, 7.0, 4.0]
    with pytest.raises(ValueError, match='exactly'):
        remap_region_positions(original.iloc[:2], [1, 3, 4])


def test_mirage_input_alignment_preserves_native_residue_identity_and_partner_chains(tmp_path, monkeypatch):
    """Exercise the PDB-to-FASTA/region adapter with the observed omission types."""
    from types import SimpleNamespace

    from Bio.PDB.Polypeptide import one_to_index

    from topoppi_affinity import external_absolute

    lines = []

    def add_residue(name, chain, number, atoms, *, insertion=" ", hetero=False):
        for atom in atoms:
            serial = len(lines) + 1
            record_type = "HETATM" if hetero else "ATOM"
            element = "CA" if name == "CA" else atom[0]
            lines.append(
                f"{record_type:<6}{serial:5d} {atom:^4} {name:>3} {chain}{number:4d}{insertion}   "
                f"{float(serial):8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{10.0:6.2f}          {element:>2}\n"
            )

    add_residue("ALA", "A", 10, ["N", "CA", "C"])
    add_residue("GLY", "A", 11, ["CA"])
    add_residue("CA", "A", 12, ["CA"], hetero=True)
    add_residue("MLY", "A", 13, ["N", "CA", "C"], hetero=True)
    add_residue("TYR", "B", 21, ["N", "CA", "C"], insertion="A")
    add_residue("LYS", "C", 3, ["N", "CA", "C"])
    path = tmp_path / "partners.pdb"
    path.write_text("".join(lines) + "END\n")
    original_region = tmp_path / "native.region"
    native_regions = pd.DataFrame({
        "pos": [1, 2, 3, 4], "res": ["ALA", "GLY", "TYR", "LYS"],
        "region": ["SUR", "COR", "RIM", "SUP"], "asa": [9., 2., 7., 4.],
    })
    native_regions.to_csv(original_region, index=False)

    def native_parser(parsed_path, *, chains):
        assert parsed_path == path
        assert chains == ["A", "B", "C"]
        # Fixture at the author parser boundary: unsupported/incomplete residues
        # are absent, while original chain order and amino-acid indices remain.
        return {"chain_id": ["A", "B", "C"], "aa": [one_to_index(aa) for aa in "AYK"]}

    upstream = SimpleNamespace(
        parse_pdb_atom14=native_parser,
        is_supported_amino_acid=lambda name: name in {"ALA", "GLY", "TYR", "LYS"},
    )
    monkeypatch.setattr(external_absolute, "_load_mirage", lambda root: upstream)
    sequences, audit = external_absolute.align_mirage_sequence_region(
        {"chains_a": "A", "chains_b": "B,C"}, path, original_region, tmp_path,
    )
    assert sequences == {"PDB|A": "A", "PDB|B": "Y", "PDB|C": "K"}
    assert audit["parsed_residues"] == 3
    assert audit["native_region_residues"] == 4
    assert [(row["resname"], row["reason"]) for row in audit["omitted_residues"]] == [
        ("GLY", "incomplete_N_CA_C"), ("CA", "unsupported_residue"), ("MLY", "unsupported_residue"),
    ]
    expected = native_regions.iloc[[0, 2, 3]].copy().reset_index(drop=True)
    expected["pos"] = [1, 2, 3]
    pd.testing.assert_frame_equal(pd.read_csv(path.with_suffix(".region")), expected)
    pd.testing.assert_frame_equal(pd.read_csv(original_region), native_regions)
