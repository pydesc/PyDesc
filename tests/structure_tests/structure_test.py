from pydesc.config import ConfigManager
from pydesc.dbhandler import MetaHandler
from pydesc.structure import StructureLoader

ConfigManager.warnings.set("quiet", True)


def assert_structure(structure):
    assert len(structure) > 0
    assert len(structure.chains) > 0
    for chain in structure.chains:
        assert len(chain) > 0


def test_default_structure_loader_load_local(pure_file, pure_types_2_mers):
    sl = StructureLoader()
    with open(pure_file) as fh:
        structures = sl.load_structures([fh])
    type_ = pure_file.parent.stem
    expected_main_mer_type = pure_types_2_mers[type_]

    if "nmr" in str(pure_file):
        assert len(structures) > 1

    for structure in structures:
        types = list(map(type, structure))
        counts = {i: 0 for i in types}

        for res in structure:
            counts[type(res)] += 1

        assert max(counts, key=lambda x: counts[x]) is expected_main_mer_type
        assert_structure(structure)


def test_default_structure_loader_load_from_pdb(pure_file, pure_types_2_mers):
    type_ = pure_file.parent.stem
    code = pure_file.stem
    sl = StructureLoader()

    with MetaHandler().open(code) as files:
        structures = sl.load_structures(files)
    expected_main_mer_type = pure_types_2_mers[type_]

    if "nmr" in str(pure_file):
        assert len(structures) > 1

    for structure in structures:
        types = list(map(type, structure))
        counts = {i: 0 for i in types}

        for res in structure:
            counts[type(res)] += 1

        most_occurring_mer_type = max(counts, key=lambda x: counts[x])
        assert most_occurring_mer_type is expected_main_mer_type
        assert_structure(structure)
