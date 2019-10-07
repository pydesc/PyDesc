import os.path

from pydesc.config import ConfigManager
from pydesc.mers import Ion
from pydesc.mers import Nucleotide
from pydesc.mers import Residue

ConfigManager.warnings.class_filters.UnknownParticleName = 'ignore'

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_STRUCTURES_DIR = os.path.join(BASE_DIR, 'data', 'test_structures')
TEST_CMAPS_DIR = os.path.join(BASE_DIR, 'data', 'validated_cmaps')

PDB_FILE_TYPES = []
PDB_FILES_WITH_TYPE = []
PDB_FILES_WITH_TYPE_SHORT = []
PDB_FILES_DICT = {}

for type_ in os.listdir(TEST_STRUCTURES_DIR):
    PDB_FILE_TYPES.append(type_)
    files_list = os.listdir(os.path.join(TEST_STRUCTURES_DIR, type_))
    PDB_FILES_WITH_TYPE_SHORT.append((type_, files_list[0]))
    PDB_FILES_DICT[type_] = []
    for struc_file in files_list:
        PDB_FILES_WITH_TYPE.append((type_, struc_file))
        PDB_FILES_DICT[type_].append(struc_file)

TYPE_DICT = {
    'dna_only': Nucleotide,
    'dna_only_nmr': Nucleotide,
    'rna_only': Nucleotide,
    'rna_only_nmr': Nucleotide,
    'prots_only': Residue,
    'prots_only_nmr': Residue,
    'PorCA_only': Ion,
}

DIR_DICT = {v: k for k, v in TYPE_DICT.items()}
