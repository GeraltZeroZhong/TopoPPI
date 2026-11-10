"""Construct the artificial peptide assembly used by the public walkthrough.

The input is an arbitrary two-residue peptide, constructed from standard residue
templates. No experimental coordinates, sequence records, or research results
are used. Coordinates can vary slightly across RDKit releases.
"""

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem


def generate() -> str:
    molecule = Chem.AddHs(Chem.MolFromSequence("AF"))
    if AllChem.EmbedMolecule(molecule, randomSeed=23) != 0:
        raise RuntimeError("Could not embed the synthetic peptide")
    AllChem.UFFOptimizeMolecule(molecule)
    molecule = Chem.RemoveHs(molecule)
    conformer = molecule.GetConformer()
    lines = [
        "REMARK Independently generated synthetic ALA-PHE peptide pair.",
        "REMARK RDKit distance geometry seed 23; UFF; chain B offset (0, 0, 4).",
        "REMARK Artificial contact geometry for software examples and tests.",
    ]
    serial = 1
    for chain, offset in (("A", 0.0), ("B", 4.0)):
        for atom in molecule.GetAtoms():
            info = atom.GetPDBResidueInfo()
            point = conformer.GetAtomPosition(atom.GetIdx())
            lines.append(
                f"ATOM  {serial:5d} {info.GetName():4s} "
                f"{info.GetResidueName():3s} {chain}{info.GetResidueNumber():4d}    "
                f"{point.x:8.3f}{point.y:8.3f}{point.z + offset:8.3f}"
                f"{1.0:6.2f}{20.0:6.2f}          {atom.GetSymbol():>2s}  "
            )
            serial += 1
        lines.append("TER")
    lines.append("END")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[2]
    content = generate()
    for output in (
        root / "tests/fixtures/synthetic_complex.pdb",
        root / "src/topoppi/assets/examples/synthetic_complex.pdb",
    ):
        output.write_text(content, encoding="ascii")
