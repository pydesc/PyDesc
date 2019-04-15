import os.path

from pydesc.monomer import Nucleotide
from pydesc.monomer import Residue
from pydesc.monomer import Ion

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_STRUCTURES_DIR = os.path.join(BASE_DIR, 'data', 'test_structures')

PDB_FILES_WITH_TYPE = []
PDB_FILES_WITH_TYPE_SHORT = []
for type_ in os.listdir(TEST_STRUCTURES_DIR):
    files_list = os.listdir(os.path.join(TEST_STRUCTURES_DIR, type_))
    PDB_FILES_WITH_TYPE_SHORT.append((type_, files_list[0]))
    for struc_file in files_list:
        PDB_FILES_WITH_TYPE.append((type_, struc_file))

TYPE_DICT = {'dna_only': Nucleotide,
             'rna_only': Nucleotide,
             'prots_only': Residue,
             'PorCA_only': Ion}
