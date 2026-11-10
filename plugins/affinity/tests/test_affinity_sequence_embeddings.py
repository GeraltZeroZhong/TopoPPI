import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from topoppi_affinity import sequence_embeddings


def test_sequence_preparation_keeps_source_specific_numbering(tmp_path, monkeypatch):
    def records(path, *_, **kwargs):
        start = 15 if path == 'absolute.pdb' else 1
        return {'A': {f'A:{start}': ('A', start, ''), f'A:{start+1}': ('R', start + 1, '')},
                'B': {'B:3': ('G', 3, '')}}

    monkeypatch.setattr(sequence_embeddings, 'chain_records', records)
    rows = pd.DataFrame([
        {'sample_id': 'abs', 'source': 'PDBbind', 'complex_id': '1abc:A:B', 'structure_path': 'absolute.pdb',
         'chains_a': 'A', 'chains_b': 'B', 'mutations': ''},
        {'sample_id': 'mut', 'source': 'SKEMPI', 'complex_id': '1abc:A:B', 'structure_path': 'cleaned.pdb',
         'chains_a': 'A', 'chains_b': 'B', 'mutations': 'AA1V'},
    ])
    path = tmp_path / 'manifest.csv'
    rows.to_csv(path, index=False)
    output = tmp_path / 'sequences'
    sequence_embeddings.prepare_sequences([path], output)
    wild = pd.read_csv(output / 'wild_type_chains.csv')
    assert set(wild.source) == {'PDBbind', 'SKEMPI'}
    assert json.loads(wild[(wild.source == 'PDBbind') & (wild.chain == 'A')].iloc[0].residue_keys)[0] == 'A:15'
    mutants = pd.read_csv(output / 'mutant_chains.csv')
    sequences = pd.read_csv(output / 'sequences.csv').set_index('sequence_id')
    assert sequences.loc[mutants.iloc[0].sequence_id, 'sequence'] == 'VR'


def test_overlapping_plm_windows_cover_every_structural_residue():
    for length in (50, 1000, 1001, 1700, 2120):
        counts = np.zeros(length, dtype=int)
        for start in sequence_embeddings.window_starts(length):
            counts[start:start + 1000] += 1
        assert counts.min() >= 1
        assert len(counts) == length


def test_appending_structures_preserves_cached_sequence_ids(tmp_path, monkeypatch):
    def records(path, *_, **kwargs):
        aa = 'W' if path == 'original.pdb' else 'A'
        return {'A': {'A:1': (aa, 1, '')}, 'B': {'B:1': ('G', 1, '')}}

    monkeypatch.setattr(sequence_embeddings, 'chain_records', records)
    output = tmp_path / 'sequences'
    for number, name in enumerate(['original', 'auxiliary']):
        frame = pd.DataFrame([{'sample_id': name, 'source': 'PDBbind', 'complex_id': name,
            'structure_path': name + '.pdb', 'chains_a': 'A', 'chains_b': 'B', 'mutations': ''}])
        manifest = tmp_path / (name + '.csv')
        frame.to_csv(manifest, index=False)
        sequence_embeddings.prepare_sequences([manifest], output)
        table = pd.read_csv(output / 'sequences.csv').set_index('sequence').sequence_id.to_dict()
        if number == 0:
            original = table.copy()
    assert all(table[sequence] == identifier for sequence, identifier in original.items())
    assert table['A'] not in set(original.values())
    assert len(pd.read_csv(output / 'wild_type_chains.csv')) == 4
    assert pd.read_csv(output / 'mutant_chains.csv').empty


@pytest.mark.parametrize("old_mutations,new_mutations,expected_chains", [
    ("A:G1V", "B:G1V", {"B"}),
    ("A:G1V,B:G1V", "B:G1V", {"B"}),
    ("A:G1V,B:G1V", "", set()),
])
def test_reusing_variant_replaces_all_mutant_chains(tmp_path, old_mutations, new_mutations, expected_chains):
    structure = Path(__file__).parent / "fixtures/tiny_complex.pdb"
    manifest = tmp_path / "inputs.csv"
    output = tmp_path / "sequences"
    common = {"source": "user", "complex_id": "same", "structure_path": str(structure),
              "chains_a": "A", "chains_b": "B"}
    pd.DataFrame([
        dict(common, sample_id="variant_001", mutations=old_mutations),
        dict(common, sample_id="retained", mutations="B:G1F"),
    ]).to_csv(manifest, index=False)
    sequence_embeddings.prepare_sequences([manifest], output)
    before = pd.read_csv(output / "mutant_chains.csv", keep_default_na=False)
    retained = before[before.sample_id.eq("retained")].to_dict("records")
    cached_ids = pd.read_csv(output / "sequences.csv").set_index("sequence").sequence_id.to_dict()

    pd.DataFrame([dict(common, sample_id="variant_001", mutations=new_mutations)]).to_csv(manifest, index=False)
    sequence_embeddings.prepare_sequences([manifest], output)
    after = pd.read_csv(output / "mutant_chains.csv", keep_default_na=False)
    current = after[after.sample_id.eq("variant_001")]
    assert set(current.chain) == expected_chains
    assert after[after.sample_id.eq("retained")].to_dict("records") == retained
    sequence_table = pd.read_csv(output / "sequences.csv").set_index("sequence")
    assert all(sequence_table.loc[sequence, "sequence_id"] == identifier for sequence, identifier in cached_ids.items())

    # A reused registry must describe the same requested mutant sequences as a
    # fresh registry, while retaining unrelated samples and reusable embeddings.
    fresh = tmp_path / "fresh"
    sequence_embeddings.prepare_sequences([manifest], fresh)
    fresh_mutants = pd.read_csv(fresh / "mutant_chains.csv", keep_default_na=False)
    fresh_sequences = pd.read_csv(fresh / "sequences.csv").set_index("sequence_id").sequence
    reused_sequences = sequence_table.reset_index().set_index("sequence_id").sequence
    assert dict(zip(current.chain, current.sequence_id.map(reused_sequences), strict=True)) == dict(
        zip(fresh_mutants.chain, fresh_mutants.sequence_id.map(fresh_sequences), strict=True))
