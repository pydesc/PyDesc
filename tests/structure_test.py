import pytest
import os.path

from tests.conftest import PDB_FILES_WITH_TYPE
from tests.conftest import TEST_STRUCTURES_DIR
from tests.conftest import TYPE_DICT

from pydesc.structure import StructureLoader

from pydesc.config import ConfigManager

ConfigManager.warnings.set("quiet", True)


def assert_structure(structure):
    assert len(structure) > 0
    assert len(structure.chains) > 0
    for chain in structure.chains:
        assert len(chain) > 0


@pytest.mark.parametrize('type_, struc_file', PDB_FILES_WITH_TYPE)
def test_default_structure_loader_load_local(type_, struc_file):

    sl = StructureLoader()
    structures = sl.load_structures(path=os.path.join(TEST_STRUCTURES_DIR, type_, struc_file))
    expected_main_mer_type = TYPE_DICT[type_]

    for structure in structures:
        types = map(type, structure)
        counts = {i: 0 for i in types}

        for res in structure:
            counts[type(res)] += 1

        assert max(counts, key=lambda x: counts[x]) is expected_main_mer_type
        assert_structure(structure)


@pytest.mark.parametrize('type_, struc_file', PDB_FILES_WITH_TYPE)
def test_default_structure_loader_load_from_pdb(type_, struc_file):
    sl = StructureLoader()
    code = os.path.splitext(struc_file)[0]
    structures = sl.load_structures(code=code)
    expected_main_mer_type = TYPE_DICT[type_]

    for structure in structures:
        types = map(type, structure)
        counts = {i: 0 for i in types}

        for res in structure:
            counts[type(res)] += 1

        assert max(counts, key=lambda x: counts[x]) is expected_main_mer_type
        assert_structure(structure)
