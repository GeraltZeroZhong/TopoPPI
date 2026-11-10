"""Legacy geometric labels, distinct from formal chemical perception.

This utility supplies display annotations for explicitly selected geometric
contacts. It does not claim to identify hydrogen bonds or other chemical
interactions and imports no ProLIF dependencies.
"""

_CHARGED_POS = frozenset({"ARG", "LYS", "HIS"})
_CHARGED_NEG = frozenset({"ASP", "GLU"})
_AROMATIC = frozenset({"PHE", "TYR", "TRP", "HIS"})
_CATION_ATOMS = frozenset({"NZ", "NH1", "NH2", "ND1", "NE2"})
_ANION_ATOMS = frozenset({"OD1", "OD2", "OE1", "OE2", "OXT"})
_POLAR_ATOMS = frozenset({"N", "O", "S", "F"})
_BACKBONE_ATOMS = frozenset({"CA", "C", "O", "N"})


def classify_contact(atom_a, atom_b, distance, config, *, context=None):
    """Classify a geometric contact using the established annotation heuristic."""
    if context is not None:
        context.check_cancelled()
    if distance > max(
        config.vdw_distance,
        config.ionic_distance,
        config.polar_contact_distance,
    ):
        return None
    res_a = atom_a.get_parent().get_resname()
    res_b = atom_b.get_parent().get_resname()
    name_a = atom_a.get_name()
    name_b = atom_b.get_name()
    elem_a = atom_a.element.upper()
    elem_b = atom_b.element.upper()
    is_ani_a = (res_a in _CHARGED_NEG and name_a in _ANION_ATOMS) or name_a == "OXT"
    is_ani_b = (res_b in _CHARGED_NEG and name_b in _ANION_ATOMS) or name_b == "OXT"
    is_cat_a = res_a in _CHARGED_POS and name_a in _CATION_ATOMS
    is_cat_b = res_b in _CHARGED_POS and name_b in _CATION_ATOMS
    is_aro_a = res_a in _AROMATIC and name_a not in _BACKBONE_ATOMS
    is_aro_b = res_b in _AROMATIC and name_b not in _BACKBONE_ATOMS
    opposite_charge = (is_cat_a and is_ani_b) or (is_ani_a and is_cat_b)
    if distance < config.ionic_distance and opposite_charge:
        return "Ionic"
    if distance < config.ionic_distance and ((is_aro_a and is_cat_b) or (is_cat_a and is_aro_b)):
        return "PiCation"
    polar_pair = elem_a in _POLAR_ATOMS and elem_b in _POLAR_ATOMS
    if distance < config.polar_contact_distance and polar_pair:
        return "PolarContact"
    if distance < config.vdw_distance:
        nonpolar_a = elem_a in {"C", "S"} and name_a not in _BACKBONE_ATOMS
        nonpolar_b = elem_b in {"C", "S"} and name_b not in _BACKBONE_ATOMS
        if nonpolar_a and nonpolar_b:
            return "Hydrophobic"
        return "VdWContact"
    return None

