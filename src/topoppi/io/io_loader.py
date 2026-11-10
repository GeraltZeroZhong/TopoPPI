"""Protein structure loading and deterministic chain-pair selection."""

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.Polypeptide import is_aa
from scipy.spatial import cKDTree


class PDBLoader:
    def __init__(self, file_path):
        self.file_path = file_path
        if str(file_path).lower().endswith((".cif", ".mmcif")):
            self.parser = MMCIFParser(QUIET=True)
        else:
            self.parser = PDBParser(QUIET=True)

        try:
            self.structure = self.parser.get_structure("P", file_path)
        except KeyError as exc:
            raise ValueError(f"Structure is missing required coordinate field {exc}: {file_path}") from exc
        self.model = next(self.structure.get_models(), None)
        if self.model is None:
            raise ValueError(f"Structure contains no coordinate model: {file_path}. "
                             "Provide a PDB or mmCIF file containing atomic coordinates.")

    def get_chain_data(self, chain_id):
        """
        Extract coordinates and atom objects for a specific chain.
        Returns:
            coords (np.ndarray): (N, 3) coordinates
            atoms (list): List of Bio.PDB.Atom objects (aligned with coords)
        """
        if chain_id not in self.model:
            raise ValueError(f"Chain {chain_id} not found in PDB.")

        chain = self.model[chain_id]

        # Filter atoms: remove HOH, use only N, CA, C, O, CB etc.
        # We generally keep all heavy atoms for surface generation
        atoms = []
        coords = []

        for residue in chain:
            # Keep canonical and recognized modified amino-acid residues, while
            # excluding waters, ions, cofactors, and unrelated hetero groups.
            if not is_aa(residue, standard=False):
                continue

            for atom in residue:
                if self._is_hydrogen(atom):
                    continue
                atoms.append(atom)
                coords.append(atom.get_coord())

        if len(coords) == 0:
            return np.empty((0, 3), dtype=float), atoms
        return np.asarray(coords, dtype=float), atoms

    @staticmethod
    def _is_hydrogen(atom) -> bool:
        """Identify hydrogen/deuterium atoms even when the element is absent."""

        element = str(getattr(atom, "element", "") or "").strip().upper()
        if element:
            return element in {"H", "D"}
        name = str(atom.get_name()).strip().upper().lstrip("0123456789")
        return name.startswith(("H", "D"))

    @staticmethod
    def _residue_identity(atom) -> tuple[str, int, str, str]:
        residue = atom.get_parent()
        chain = residue.get_parent()
        return (
            str(chain.id),
            int(residue.id[1]),
            str(residue.id[2]).strip(),
            str(residue.get_resname()),
        )

    def get_protein_chain_ids(self):
        """
        Return chain IDs that contain at least one recognized amino-acid residue.
        """
        chain_ids = []
        for chain in self.model:
            has_protein_residue = any(is_aa(residue, standard=False) for residue in chain)
            if has_protein_residue:
                chain_ids.append(chain.id)
        return chain_ids

    def get_chain_residue_count(self, chain_id):
        """
        Return the number of recognized amino-acid residues in a chain.
        """
        if chain_id not in self.model:
            raise ValueError(f"Chain {chain_id} not found in PDB.")
        chain = self.model[chain_id]
        return sum(1 for residue in chain if is_aa(residue, standard=False))

    def _chain_contact_data(self, chain_ids, min_chain_residues):
        data = {}
        for chain_id in chain_ids:
            residue_count = self.get_chain_residue_count(chain_id)
            if residue_count < int(min_chain_residues):
                continue
            coords, atoms = self.get_chain_data(chain_id)
            if not len(coords):
                continue
            data[chain_id] = {
                "coords": coords,
                "atoms": atoms,
                "residue_count": residue_count,
                "atom_count": int(len(coords)),
            }
        return data

    def _chain_pair_contact_record(self, left, right, chain_data, trees, distance_cutoff):
        left_coords = chain_data[left]["coords"]
        tree = trees[right]
        neighborhoods = tree.query_ball_point(left_coords, r=float(distance_cutoff))
        residue_pairs = {
            (
                self._residue_identity(chain_data[left]["atoms"][left_atom_index]),
                self._residue_identity(chain_data[right]["atoms"][right_atom_index]),
            )
            for left_atom_index, right_indices in enumerate(neighborhoods)
            for right_atom_index in right_indices
        }
        return {
            "left": left,
            "right": right,
            "contact_atom_pair_count": int(sum(len(indices) for indices in neighborhoods)),
            "contact_residue_pair_count": int(len(residue_pairs)),
            "minimum_atom_distance": float(np.min(tree.query(left_coords, k=1)[0])),
            "combined_residue_count": int(chain_data[left]["residue_count"] + chain_data[right]["residue_count"]),
        }

    @staticmethod
    def _best_contact_record(candidates):
        return min(
            candidates,
            key=lambda item: (
                -item["contact_residue_pair_count"],
                -item["contact_atom_pair_count"],
                item["minimum_atom_distance"],
                -item["combined_residue_count"],
                str(item["left"]),
                str(item["right"]),
            ),
        )

    def select_contact_chain_pair(self, min_chain_residues=1, distance_cutoff=9.0):
        """Select by contact-residue pairs, then atom contacts and stable ties."""

        chain_data = self._chain_contact_data(self.get_protein_chain_ids(), min_chain_residues)
        eligible = list(chain_data)
        if len(eligible) < 2:
            raise ValueError(f"Need at least two protein chains with >= {min_chain_residues} residues and heavy atoms.")

        trees = {chain_id: cKDTree(data["coords"]) for chain_id, data in chain_data.items()}
        candidates = []
        for left_index, left in enumerate(eligible):
            for right in eligible[left_index + 1 :]:
                candidates.append(self._chain_pair_contact_record(left, right, chain_data, trees, distance_cutoff))

        selected = self._best_contact_record(candidates)
        left, right = str(selected["left"]), str(selected["right"])
        # Use the larger chain as the surface chain; ties are lexical.
        left_size = int(chain_data[left]["atom_count"])
        right_size = int(chain_data[right]["atom_count"])
        if right_size > left_size or (right_size == left_size and right < left):
            left, right = right, left
        return (
            left,
            right,
            {
                "selection_mode": "auto_contact",
                "distance_cutoff_angstrom": float(distance_cutoff),
                "contact_atom_pair_count": int(selected["contact_atom_pair_count"]),
                "contact_residue_pair_count": int(selected["contact_residue_pair_count"]),
                "minimum_atom_distance_angstrom": float(selected["minimum_atom_distance"]),
                "eligible_chain_count": int(len(eligible)),
            },
        )

    def select_contact_chain_pair_between_groups(
        self,
        receptor_chain_ids,
        ligand_chain_ids,
        *,
        min_chain_residues=1,
        distance_cutoff=6.0,
    ):
        """Select the dominant contacting pair across declared partner groups."""

        available = set(self.get_protein_chain_ids())
        receptor = tuple(dict.fromkeys(str(chain_id) for chain_id in receptor_chain_ids))
        ligand = tuple(dict.fromkeys(str(chain_id) for chain_id in ligand_chain_ids))
        if not receptor or not ligand:
            raise ValueError("Both receptor and ligand chain groups must be non-empty.")
        if set(receptor) & set(ligand):
            raise ValueError("Receptor and ligand chain groups must be disjoint.")
        missing = sorted((set(receptor) | set(ligand)) - available)
        if missing:
            raise ValueError(f"Declared protein chains are missing: {missing}; available={sorted(available)}")

        requested = (*receptor, *ligand)
        chain_data = self._chain_contact_data(requested, min_chain_residues)
        ineligible = [chain_id for chain_id in requested if chain_id not in chain_data]
        if ineligible:
            raise ValueError(
                f"Declared chains have fewer than {min_chain_residues} recognized amino-acid residues "
                f"or no heavy atoms: {ineligible}"
            )
        trees = {chain_id: cKDTree(data["coords"]) for chain_id, data in chain_data.items()}
        candidates = [
            self._chain_pair_contact_record(left, right, chain_data, trees, distance_cutoff)
            for left in receptor
            for right in ligand
        ]
        selected = self._best_contact_record(candidates)
        total_atom_contacts = int(sum(item["contact_atom_pair_count"] for item in candidates))
        total_residue_contacts = int(sum(item["contact_residue_pair_count"] for item in candidates))
        return (
            str(selected["left"]),
            str(selected["right"]),
            {
                "selection_mode": "declared_partner_groups_dominant_contact_pair",
                "distance_cutoff_angstrom": float(distance_cutoff),
                "receptor_chain_ids": list(receptor),
                "ligand_chain_ids": list(ligand),
                "candidate_chain_pair_count": int(len(candidates)),
                "contact_atom_pair_count": int(selected["contact_atom_pair_count"]),
                "contact_residue_pair_count": int(selected["contact_residue_pair_count"]),
                "minimum_atom_distance_angstrom": float(selected["minimum_atom_distance"]),
                "total_cross_group_atom_pair_count": total_atom_contacts,
                "total_cross_group_residue_pair_count": total_residue_contacts,
                "selected_atom_contact_fraction": (
                    float(selected["contact_atom_pair_count"] / total_atom_contacts) if total_atom_contacts else 0.0
                ),
                "selected_residue_contact_fraction": (
                    float(selected["contact_residue_pair_count"] / total_residue_contacts)
                    if total_residue_contacts
                    else 0.0
                ),
            },
        )
