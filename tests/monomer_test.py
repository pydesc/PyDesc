import pytest
import os.path
import Bio.PDB

from pydesc.monomer import MonomerFactory
from pydesc.monomer import Nucleotide
from pydesc.monomer import Residue


BASE_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_STRUCTURES = os.path.join(BASE_DIR, 'data', 'test_structures')

PDB_FILES_WITH_TYPE = [(type_, struc_file) for type_ in os.listdir(TEST_STRUCTURES) for struc_file in os.listdir(os.path.join(TEST_STRUCTURES, type_))]


@pytest.mark.parametrize('type_, struc_file', PDB_FILES_WITH_TYPE)
def test_default_monomer_factory(type_, struc_file):

    factory = MonomerFactory()

    types_dict = {'dna_only': Nucleotide,
                  'rna_only': Nucleotide,
                  'prots_only': Residue}

    pdb_structure = Bio.PDB.PDBParser(QUIET=True).get_structure(struc_file,
                                                                os.path.join(TEST_STRUCTURES,
                                                                             type_,
                                                                             struc_file)
                                                                )

    expected_type = types_dict[type_]

    points = {True: 0, False: 1}
    length = 0

    for model in pdb_structure:
        for chain in model:
            for residue in chain:
                result, warns = factory.create_from_biopdb(residue, warn_in_place=False)

                if result is None:
                    assert warns is None
                    assert residue.get_id()[0] == 'W'
                else:
                    length += 1
                    points[expected_type in result] += 1

    assert length - points[True] < 10, '%i out of %i residues are of expected type' % (points[True], length)
