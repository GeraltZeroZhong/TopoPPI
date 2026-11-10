# Synthetic fixtures

All fixtures in this directory are independently constructed software inputs.

- `tiny_complex.pdb` uses hand-set integer and half-integer coordinates for two
  two-residue chains and one water molecule. It tests geometry, chain selection,
  and heteroatom filtering. Its incomplete residues are intentional.
- `synthetic_complex.pdb` contains two artificial ALA-PHE peptides. The generator
  uses RDKit residue templates, distance geometry with seed 23, and UFF relaxation;
  chain B is a copy translated by 4 Å along the z axis. The arbitrary
  sequence and coordinates support chemical perception and GUI examples.
- `prolif_interactions.json` contains hand-written residue tokens and interaction
  labels to test normalization, insertion codes, invalid entries, and new labels.

Regenerate the peptide example with
`python tests/fixtures/generate_synthetic_example.py` in an environment containing
the project dependencies. Both the test fixture and bundled example are written.
RDKit releases may produce slightly different relaxed coordinates.

The separately documented experimental 1BVK structure is the package's default
interactive example. Synthetic tests use the fixtures above.
