import os.path

from pydesc.config import ConfigManager
from pydesc.structure import StructureLoader
from tests.conftest import PURE_TYPES_2_MERS_DICT
from tests.conftest import TEST_STRUCTURES_DIR

ConfigManager.warnings.set("quiet", True)


def assert_structure(structure):
    assert len(structure) > 0
    assert len(structure.chains) > 0
    for chain in structure.chains:
        assert len(chain) > 0


def test_default_structure_loader_load_local(structure_file_w_pure_type):
    sl = StructureLoader()
    path_str = os.path.join(TEST_STRUCTURES_DIR, structure_file_w_pure_type)
    structures = sl.load_structures(path=path_str)
    type_ = os.path.split(structure_file_w_pure_type)[0]
    expected_main_mer_type = PURE_TYPES_2_MERS_DICT[type_]

    if "nmr" in structure_file_w_pure_type:
        assert len(structures) > 1

    for structure in structures:
        types = list(map(type, structure))
        counts = {i: 0 for i in types}

        for res in structure:
            counts[type(res)] += 1

        assert max(counts, key=lambda x: counts[x]) is expected_main_mer_type
        assert_structure(structure)


def test_default_structure_loader_load_from_pdb(structure_file_w_pure_type):
    sl = StructureLoader()
    type_, file_ = os.path.split(structure_file_w_pure_type)
    code = os.path.splitext(file_)[0]
    structures = sl.load_structures(code=code)
    expected_main_mer_type = PURE_TYPES_2_MERS_DICT[type_]

    if "nmr" in structure_file_w_pure_type:
        assert len(structures) > 1

    for structure in structures:
        types = list(map(type, structure))
        counts = {i: 0 for i in types}

        for res in structure:
            counts[type(res)] += 1

        most_occurring_mer_type = max(counts, key=lambda x: counts[x])
        assert most_occurring_mer_type is expected_main_mer_type
        assert_structure(structure)
