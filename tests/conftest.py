import os.path

from pydesc.config import ConfigManager
from pydesc.mers.full_atom import Ion
from pydesc.mers.full_atom import Nucleotide
from pydesc.mers.full_atom import Residue

ConfigManager.warnings.class_filters.UnknownParticleName = "ignore"

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_STRUCTURES_DIR = os.path.join(BASE_DIR, "data", "test_structures")
TEST_CMAPS_DIR = os.path.join(BASE_DIR, "data", "validated_cmaps")

CIF_FILES_WITH_TYPE = []


FILE_TYPES = []
PDB_FILES_WITH_TYPE = []
PDB_FILES_WITH_TYPE_SHORT = []
PDB_FILES_DICT = {}

for type_ in os.listdir(TEST_STRUCTURES_DIR):
    FILE_TYPES.append(type_)

    files_list = os.listdir(os.path.join(TEST_STRUCTURES_DIR, type_))

    PDB_FILES_WITH_TYPE_SHORT.append((type_, files_list[0]))
    PDB_FILES_DICT[type_] = []

    for structure_file in files_list:
        if structure_file.endswith(".pdb"):
            PDB_FILES_WITH_TYPE.append((type_, structure_file))
            PDB_FILES_DICT[type_].append(structure_file)
        elif structure_file.endswith(".cif"):
            CIF_FILES_WITH_TYPE.append((type_, structure_file))

PDB_FILES_WITH_PURE_TYPE = [
    (type_, file_name) for type_, file_name in PDB_FILES_WITH_TYPE if "only" in type_
]
PDB_FILES_WITH_PURE_TYPE_SHORT = [
    (type_, file_name)
    for type_, file_name in PDB_FILES_WITH_TYPE_SHORT
    if "only" in type_
]

PURE_TYPES_2_MERS_DICT = {
    "dna_only": Nucleotide,
    "dna_only_nmr": Nucleotide,
    "rna_only": Nucleotide,
    "rna_only_nmr": Nucleotide,
    "prots_only": Residue,
    "prots_only_nmr": Residue,
    "PorCA_only": Ion,
}

DIR_DICT = {v: k for k, v in list(PURE_TYPES_2_MERS_DICT.items())}
