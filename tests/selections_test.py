import pytest
import os.path

from conftest import TEST_STRUCTURES_DIR
from conftest import PDB_FILES_WITH_TYPE_SHORT

from pydesc.structure import StructureLoader, AbstractStructure
from pydesc import selection


@pytest.mark.parametrize('type_, struc_file', PDB_FILES_WITH_TYPE_SHORT)
def test_everything_sele(type_, struc_file):

    sl = StructureLoader()

    structures = sl.load_structures(path=os.path.join(TEST_STRUCTURES_DIR, type_, struc_file))

    for structure in structures:

        sel = selection.Everything().specify(structure)
        assert type(sel) is selection.Set
        assert len(tuple(sel)) == len(structure)

        stc = selection.Everything().create_structure(structure)
        assert isinstance(stc, AbstractStructure)
        assert len(stc) == len(structure)
        assert stc[0] is structure[0]

        # new_stc = selection.Everything().create_new_structure(structure)
        # assert isinstance(new_stc, AbstractStructure)
        # assert len(new_stc) == len(stc)
        # assert new_stc[0] is not structure[0]


@pytest.mark.parametrize('type_, struc_file', PDB_FILES_WITH_TYPE_SHORT)
def test_set_sele(type_, struc_file):

    sl = StructureLoader()

    structures = sl.load_structures(path=os.path.join(TEST_STRUCTURES_DIR, type_, struc_file))

    for structure in structures:

        sel = selection.Set([i.get_pdb_id() for i in tuple(structure)[:6]])
        assert type(sel) is selection.Set
        assert len(tuple(sel)) == 6

        new_sel = sel.specify(structure)
        assert type(new_sel) is selection.Set
        assert len(tuple(new_sel)) == 6

        stc = sel.create_structure(structure)
        assert isinstance(stc, AbstractStructure)
        assert len(stc) == 6
        assert stc[0] is structure[0]


