import pytest
import os.path

from conftest import TEST_STRUCTURES_DIR
from conftest import PDB_FILES_WITH_TYPE

from pydesc.structure import StructureLoader


@pytest.mark.parametrize('type_, struc_file', PDB_FILES_WITH_TYPE)
def test_default_structure_loader_load_local(type_, struc_file):

    sl = StructureLoader()

    structures = sl.load_structures(path=os.path.join(TEST_STRUCTURES_DIR, type_, struc_file))

    for structure in structures:
        structure.set_contact_map()

        assert structure.contact_map