import os.path
import time
import pytest

from tests.conftest import TEST_STRUCTURES_DIR
from tests.conftest import PDB_FILES_WITH_TYPE_SHORT

from pydesc.structure import StructureLoader, AbstractStructure
from pydesc import selection
from pydesc.config import ConfigManager

ConfigManager.warnings.quiet = True


@pytest.mark.parametrize('type_, struc_file', PDB_FILES_WITH_TYPE_SHORT)
class TestEverythingSelection:
    # TODO: test on substructures

    @staticmethod
    def get_structure(type_, structure_file):
        sl = StructureLoader()
        structures = sl.load_structures(path=os.path.join(TEST_STRUCTURES_DIR, type_, structure_file))
        return structures

    def test_everything_create_structure(self, type_, struc_file):
        structures = self.get_structure(type_, struc_file)
        for structure in structures:
            stc = selection.Everything().create_structure(structure)
            assert isinstance(stc, AbstractStructure)
            assert len(stc) == len(structure)
            assert stc[0] is structure[0]
            assert sorted(stc[0].atoms.keys()) == sorted(structure[0].atoms.keys())

    @pytest.mark.long
    def test_everything_create_new_structure(self, type_, struc_file):
        load_start = time.time()
        structures = self.get_structure(type_, struc_file)
        load_end = time.time()
        load_time = load_end - load_start
        for structure in structures:
            select_start = time.time()
            new_stc = selection.Everything().create_new_structure(structure)
            select_end = time.time()
            select_time = select_end - select_start
            assert isinstance(new_stc, AbstractStructure)
            assert len(new_stc) == len(structure)
            assert new_stc[0] is not structure[0]
            assert sorted(new_stc[0].atoms.keys()) == sorted(structure[0].atoms.keys())
            assert sorted(new_stc[0].pseudoatoms.keys()) == sorted(structure[0].pseudoatoms.keys())
            assert select_time <= load_time

    def test_everything_specify(self, type_, struc_file):
        structures = self.get_structure(type_, struc_file)
        for structure in structures:
            sel = selection.Everything().specify(structure)
            assert type(sel) is selection.Set
            assert len(tuple(sel)) == len(structure)


@pytest.mark.parametrize('type_, struc_file', PDB_FILES_WITH_TYPE_SHORT)
class TestSetSelection:

    @staticmethod
    def get_structure(type_, structure_file):
        sl = StructureLoader()
        structures = sl.load_structures(path=os.path.join(TEST_STRUCTURES_DIR, type_, structure_file))
        return structures

    def test_set_specify(self, type_, struc_file):
        structures = self.get_structure(type_, struc_file)
        for structure in structures:
            sel = selection.Set([i.get_pdb_id() for i in tuple(structure)[:6]])
            assert type(sel) is selection.Set
            assert len(tuple(sel)) == 6

            new_sel = sel.specify(structure)
            assert type(new_sel) is selection.Set
            assert len(tuple(new_sel)) == 6

    def test_set_create_structure(self, type_, struc_file):
        structures = self.get_structure(type_, struc_file)
        for structure in structures:
            sel = selection.Set([i.get_pdb_id() for i in tuple(structure)[:6]])
            stc = sel.create_structure(structure)
            assert isinstance(stc, AbstractStructure)
            assert len(stc) == 6
            assert stc[0] is structure[0]
            assert sorted(stc[0].atoms.keys()) == sorted(structure[0].atoms.keys())

    def test_set_create_new_structure(self, type_, struc_file):
        structures = self.get_structure(type_, struc_file)
        for structure in structures:
            sel = selection.Set([i.get_pdb_id() for i in tuple(structure)[:6]])
            stc = sel.create_new_structure(structure)
            assert isinstance(stc, AbstractStructure)
            assert len(stc) == 6
            assert stc[0] is not structure[0]
            assert sorted(stc[0].atoms.keys()) == sorted(structure[0].atoms.keys())
