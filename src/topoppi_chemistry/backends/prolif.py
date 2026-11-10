"""ProLIF perception backend using copied Core atoms and unchanged defaults."""
from __future__ import annotations

import logging
import re
import warnings
from collections import defaultdict
from contextlib import contextmanager, nullcontext
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory

from topoppi.interactions.records import residue_sequence_token

PROLIF_TO_STANDARD = {
    "HydrogenBond": "HydrogenBond",
    "HBAcceptor": "HydrogenBond",
    "HBDonor": "HydrogenBond",
    "Hydrophobic": "Hydrophobic",
    "PiStacking": "PiStacking",
    "FaceToFace": "PiStacking",
    "EdgeToFace": "PiStacking",
    "PiCation": "PiCation",
    "CationPi": "PiCation",
    "Cationic": "Ionic",
    "Anionic": "Ionic",
    "Ionic": "Ionic",
    "HalogenBond": "HalogenBond",
    "XBAcceptor": "HalogenBond",
    "XBDonor": "HalogenBond",
    "MetalAcceptor": "MetalCoordination",
    "MetalDonor": "MetalCoordination",
    "MetalCoordination": "MetalCoordination",
    "VdWContact": "VdWContact",
}


_RESIDUE_TOKEN_RE = re.compile(r"(?<!\d)(-?\d+)([A-Za-z]?)")


_PROLIF_RESIDUE_RE = re.compile(r"^\s*([A-Za-z0-9 ]*?)(-?\d+)(?:\.[^.\s]+)?\s*$")


def normalize_interaction_name(name):
    if not name:
        return None
    compact = str(name).replace(" ", "").replace("-", "")
    return PROLIF_TO_STANDARD.get(compact, "Other")


def _prolif_residue_identity(value):
    match = _PROLIF_RESIDUE_RE.match(str(value))
    if match is None:
        return None
    return match.group(1).strip().upper(), int(match.group(2))


def _chain_residue_token_map(atom_group, residue_number_offset=0):
    mapping = defaultdict(set)
    for residue in atom_group.residues:
        sequence = int(residue.resid)
        insertion = str(getattr(residue, "icode", "") or "").strip().upper()
        token = f"{sequence}{insertion}"
        name = str(residue.resname).strip().upper()
        mapping[(name, sequence + residue_number_offset)].add(token)
        mapping[(None, sequence + residue_number_offset)].add(token)
    return mapping


def _resolve_generated_residue_token(value, token_map):
    fallback = residue_sequence_token(value)
    if fallback is None or token_map is None:
        return fallback
    identity = _prolif_residue_identity(value)
    if identity is not None:
        candidates = token_map.get(identity, set())
        if len(candidates) == 1:
            return next(iter(candidates))
        if len(candidates) > 1:
            return None
        sequence = identity[1]
    else:
        sequence = int(_RESIDUE_TOKEN_RE.search(fallback).group(1))
    candidates = token_map.get((None, sequence), set())
    if len(candidates) == 1:
        return next(iter(candidates))
    if len(candidates) > 1:
        return None
    return fallback


def _to_records(dataframe, residue_tokens_a=None, residue_tokens_b=None, log=None):
    records = []
    if dataframe is None or dataframe.empty:
        return records
    skipped_ambiguous = 0
    observed = set()
    for column, values in dataframe.items():
        if not isinstance(column, tuple) or len(column) < 3:
            continue
        if not values.fillna(False).astype(bool).any():
            continue
        ligand_res, protein_res, interaction_name = column[-3:]
        res_a = _resolve_generated_residue_token(protein_res, residue_tokens_a)
        res_b = _resolve_generated_residue_token(ligand_res, residue_tokens_b)
        if res_a is None or res_b is None:
            skipped_ambiguous += 1
            continue
        interaction = normalize_interaction_name(interaction_name)
        key = (res_a, res_b, interaction)
        if key in observed:
            continue
        observed.add(key)
        records.append({"res_a_seq": res_a, "res_b_seq": res_b, "interaction": interaction})
    if skipped_ambiguous and log is not None:
        log.warning(
            "Skipped %d ProLIF interaction columns whose insertion-coded residue could not be resolved uniquely.",
            skipped_ambiguous,
        )
    return records


def _select_chain_atoms(universe, chain):
    for query in (f"chainID {chain}", f"segid {chain}"):
        atoms = universe.select_atoms(query)
        if len(atoms):
            return atoms
    return None


@contextmanager
def _prolif_structure_input(structure_path, chain_a=None, chain_b=None):
    """Yield the same selected heavy atoms used for surface geometry.

    Resolve alternate conformers with the core loader, then clear their altloc
    markers in this temporary input so RDKit retains the selected coordinates.
    Surrogate chain IDs allow both PDB and mmCIF author chains to pass through
    the PDB-based chemical-perception step; output records keep the author IDs.
    """
    from Bio.PDB import PDBIO
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Model import Model
    from Bio.PDB.Structure import Structure

    from topoppi.io.atom_records import atom_records, restore_atoms
    from topoppi.io.io_loader import PDBLoader

    if hasattr(structure_path, "chain_data"):
        interface = _as_interface(structure_path)
        chain_a, chain_b = (partner.chains[0] for partner in interface.partners)
        data = {partner.chains[0]: (partner.coords, partner.atoms) for partner in interface.partners}
        get_chain_data = data.__getitem__
    else:
        loader = PDBLoader(Path(structure_path))
        get_chain_data = loader.get_chain_data
    converted = Structure("TopoPPI")
    model = Model(0)
    converted.add(model)
    surrogate_ids = ("A", "B")
    for source_id, surrogate_id in zip((str(chain_a), str(chain_b)), surrogate_ids, strict=True):
        coords, selected = get_chain_data(source_id)
        if not selected:
            raise ValueError(f"Chain {source_id!r} has no protein heavy atoms for interaction detection.")
        target_chain = Chain(surrogate_id)
        # Fresh Atom objects contain only the selected child of a disordered
        # atom; copying its Bio.PDB wrapper would retain every conformer.
        for atom in restore_atoms(atom_records(selected), coords):
            atom.set_altloc(" ")
            residue = atom.get_parent()
            if residue.id not in target_chain:
                residue.detach_parent()
                target_chain.add(residue)
        model.add(target_chain)

    with TemporaryDirectory(prefix="topoppi-prolif-") as tmpdir:
        converted_path = Path(tmpdir) / "selected-chains.pdb"
        writer = PDBIO()
        writer.set_structure(converted)
        writer.save(str(converted_path))
        yield str(converted_path), *surrogate_ids


def _atom_group_to_pdb_block(atom_group):
    """Serialize an AtomGroup in memory for RDKit's PDB chemical perception."""

    from MDAnalysis.coordinates.PDB import PDBWriter
    from MDAnalysis.lib.util import NamedStream

    buffer = StringIO()
    stream = NamedStream(buffer, "topoppi-prolif-chain.pdb", close=False)
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=r"Found no information for attr: 'formalcharges'.*",
            )
            warnings.filterwarnings(
                "ignore",
                message=r"Unit cell dimensions not found.*",
            )
            with PDBWriter(
                stream,
                n_atoms=atom_group.n_atoms,
                bonds=None,
                multiframe=False,
            ) as writer:
                writer.write(atom_group)
        return buffer.getvalue()
    finally:
        buffer.close()


def _mda_to_prolif_with_explicit_hydrogen(atom_group, chain_label, plf, log, *, residue_number_offset=0):
    """Build an isolated, hydrogenated ProLIF molecule for fingerprinting.

    The in-memory PDB round-trip lets RDKit perceive standard protein bond
    orders before hydrogens are added.  Both operations create fingerprint-only
    molecule objects; the source AtomGroup and its coordinates are read but
    never modified.
    """

    from rdkit import Chem, rdBase

    pdb_block = _atom_group_to_pdb_block(atom_group)
    quiet_rdkit = logging.getLogger().getEffectiveLevel() > logging.DEBUG
    with rdBase.BlockLogs() if quiet_rdkit else nullcontext():
        mol = Chem.MolFromPDBBlock(
            pdb_block,
            sanitize=True,
            removeHs=False,
            proximityBonding=True,
        )
    if mol is None:
        raise ValueError(
            f"Could not prepare Chain {str(chain_label).strip()!r} for interaction detection. "
            "Check that the chain contains complete protein residues, or provide a ProLIF JSON file."
        )

    input_hydrogen_count = sum(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms())
    with rdBase.BlockLogs() if quiet_rdkit else nullcontext():
        prepared = Chem.AddHs(
            mol,
            addCoords=True,
            addResidueInfo=True,
        )
    explicit_hydrogen_count = sum(atom.GetAtomicNum() == 1 for atom in prepared.GetAtoms())
    if not explicit_hydrogen_count:
        raise ValueError(
            f"Could not add hydrogens to Chain {str(chain_label).strip()!r}. "
            "Check the selected structure, or provide a ProLIF JSON file."
        )

    log.debug(
        "Prepared an isolated ProLIF fingerprint molecule for Chain %s with %d explicit hydrogens (%d added).",
        chain_label,
        explicit_hydrogen_count,
        explicit_hydrogen_count - input_hydrogen_count,
    )
    if residue_number_offset:
        # ProLIF stores residue numbers as unsigned integers. Shift only its
        # private RDKit copy; the token map restores author residue numbers.
        for atom in prepared.GetAtoms():
            info = atom.GetPDBResidueInfo()
            if info is not None:
                info.SetResidueNumber(info.GetResidueNumber() + residue_number_offset)
    return plf.Molecule(prepared)


def _as_interface(value):
    from topoppi.interface import Interface, PreparedInterface

    interface = value.interface if isinstance(value, PreparedInterface) else value
    if not isinstance(interface, Interface):
        raise TypeError("ProLIF detection requires a Core Interface or GeometryRepresentation.")
    if any(len(partner.chains) != 1 for partner in interface.partners):
        raise ValueError("ProLIF detection currently requires one author chain per partner.")
    return interface


def detect(interface, *, context=None):
    """Detect contacts from copied Core atoms without reopening their source."""
    from topoppi.execution import ExecutionContext
    from topoppi.interactions.records import InteractionResult

    context = context or ExecutionContext()
    context.check_cancelled()
    interface = _as_interface(interface)
    chain_a, chain_b = (partner.chains[0] for partner in interface.partners)
    log = context.logger
    quiet_third_party = logging.getLogger().getEffectiveLevel() > logging.DEBUG
    try:
        with warnings.catch_warnings():
            if quiet_third_party:
                warnings.simplefilter("ignore", DeprecationWarning)
            import MDAnalysis as mda

            if quiet_third_party:
                warnings.simplefilter("ignore", DeprecationWarning)
            import prolif as plf
    except ImportError as exc:
        raise RuntimeError(
            "ProLIF auto-generation requires MDAnalysis, ProLIF, and RDKit. "
            "Install the chemistry dependencies: pip install 'topoppi[chemistry]'."
        ) from exc

    context.emit("chemistry.preparing", backend="prolif")
    with _prolif_structure_input(interface) as (mda_path, mda_chain_a, mda_chain_b):
        universe = mda.Universe(mda_path)
        chain_a_atoms = _select_chain_atoms(universe, mda_chain_a)
        chain_b_atoms = _select_chain_atoms(universe, mda_chain_b)
        if chain_a_atoms is None or chain_b_atoms is None:
            raise ValueError(
                f"ProLIF atom selection failed for chain {str(chain_a).strip()!r} or {str(chain_b).strip()!r}."
            )
        offset_a = max(0, -min(int(residue.resid) for residue in chain_a_atoms.residues))
        offset_b = max(0, -min(int(residue.resid) for residue in chain_b_atoms.residues))
        mol_a = _mda_to_prolif_with_explicit_hydrogen(
            chain_a_atoms, chain_a, plf, log, residue_number_offset=offset_a)
        mol_b = _mda_to_prolif_with_explicit_hydrogen(
            chain_b_atoms, chain_b, plf, log, residue_number_offset=offset_b)
        context.check_cancelled()
        context.emit("chemistry.detecting", backend="prolif")
        fp = plf.Fingerprint()
        fp.run_from_iterable([mol_b], mol_a, progress=False, n_jobs=1)
        records = _to_records(
            fp.to_dataframe(),
            residue_tokens_a=_chain_residue_token_map(chain_a_atoms, offset_a),
            residue_tokens_b=_chain_residue_token_map(chain_b_atoms, offset_b),
            log=log,
        )
    context.check_cancelled()
    return InteractionResult("prolif", chain_a, chain_b, records,
                             {"analysis_id": interface.analysis_id, "source": dict(interface.source)})
