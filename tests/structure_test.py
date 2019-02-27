import pytest
import os.path

from conftest import PDB_FILES_WITH_TYPE
from conftest import TEST_STRUCTURES_DIR
from conftest import TYPE_DICT

from pydesc.structure import StructureLoader


@pytest.mark.parametrize('type_, struc_file', PDB_FILES_WITH_TYPE)
def test_default_structure_loader(type_, struc_file):

    sl = StructureLoader()

    struc = sl.load_structure(path=os.path.join(TEST_STRUCTURES_DIR, struc_file))

    expected_main_mer_type = TYPE_DICT[type_]

    types = map(type, struc)
    counts = {i: 0 for i in types}
    for i in types:
        counts[i] += 1

    assert max(counts, key=lambda x: counts[x]) is expected_main_mer_type
