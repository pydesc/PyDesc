import pytest
import os.path
import Bio.PDB

from pydesc.monomer import MonomerFactory

from conftest import PDB_FILES_WITH_TYPE
from conftest import TEST_STRUCTURES_DIR
from conftest import TYPE_DICT


@pytest.mark.parametrize('type_, struc_file', PDB_FILES_WITH_TYPE)
def test_default_monomer_factory_create_from_pdbres(type_, struc_file):

    factory = MonomerFactory()

    pdb_structure = Bio.PDB.PDBParser(QUIET=True).get_structure(struc_file,
                                                                os.path.join(TEST_STRUCTURES_DIR,
                                                                             type_,
                                                                             struc_file)
                                                                )

    expected_type = TYPE_DICT[type_]

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

    assert length - points[True] < 5, '%i out of %i residues are of expected type' % (points[True], length)
