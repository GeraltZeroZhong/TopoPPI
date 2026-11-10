"""Residue dictionaries for surface organization and proximity descriptors."""

import numpy as np

CHEMICAL_CLASSES = ("P", "A", "C", "X")
POLAR = frozenset("CHNQSTYW")
APOLAR = frozenset("AFGIVLMP")
CHARGED = frozenset("DEKR")
HYDROPHOBIC = frozenset("ACFILMVWY")
AA3_TO_1 = dict(
    zip(
        "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
        "ARNDCQEGHILKMFPSTWYV", strict=False,
    )
)
HYDROPATHY = dict(
    zip(
        "ACDEFGHIKLMNPQRSTVWY",
        [
            1.8,
            2.5,
            -3.5,
            -3.5,
            2.8,
            -0.4,
            -3.2,
            4.5,
            -3.9,
            3.8,
            1.9,
            -3.5,
            -1.6,
            -3.5,
            -4.5,
            -0.8,
            -0.7,
            4.2,
            -0.9,
            -1.3,
        ], strict=False,
    )
)


def aa_from_residue_key(key):
    """Read CHAIN:RESNAME:POSITION, optionally prefixed by side0/side1."""
    name = str(key).split(":")[-2].upper()
    return name if name in HYDROPATHY else AA3_TO_1.get(name, "X")


def classify_aa(aa):
    aa = str(aa).upper()
    if aa in POLAR:
        return "P"
    if aa in APOLAR:
        return "A"
    if aa in CHARGED:
        return "C"
    return "X"


def is_hydrophobic(aa):
    return str(aa).upper() in HYDROPHOBIC


def chemical_indices(amino_acids):
    return np.asarray([CHEMICAL_CLASSES.index(classify_aa(a)) for a in amino_acids], dtype=np.int8)
