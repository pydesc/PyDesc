import itertools
import os.path

import pytest

from pydesc import selection
from pydesc.config import ConfigManager
from pydesc.mers import Ion
from pydesc.mers import Ligand
from pydesc.mers import MerChainable
from pydesc.mers import Nucleotide
from pydesc.mers import Residue
from pydesc.mers import WrongMerType
from pydesc.structure import AbstractStructure
from pydesc.structure import PartialStructure
from pydesc.structure import StructureLoader
from pydesc.warnexcept import DiscontinuityError
from tests.conftest import PDB_FILES_WITH_TYPE_SHORT
from tests.conftest import TEST_STRUCTURES_DIR

ConfigManager.warnings.quiet = True


@pytest.fixture(scope="module", params=PDB_FILES_WITH_TYPE_SHORT)
def structure(request):
    sl = StructureLoader()
    type_, structure_file = request.param
    path_str = os.path.join(TEST_STRUCTURES_DIR, type_, structure_file)
    structure = sl.load_structures(path=path_str)[0]
    return structure


@pytest.fixture(scope="module")
def stc_2dlc():
    sl = StructureLoader()
    pth = os.path.join(TEST_STRUCTURES_DIR, "mixed", "2DLC.cif")
    return sl.load_structures(path=pth)[0]


class TestEverythingSelection:
    def test_create_structure(self, structure):
        stc = selection.Everything().create_structure(structure)
        assert isinstance(stc, AbstractStructure)
        assert len(stc) == len(structure)
        assert stc[0] is structure[0]
        assert sorted(stc[0].atoms.keys()) == sorted(structure[0].atoms.keys())

    @pytest.mark.long
    def test_create_new_structure(self, structure):
        new_stc = selection.Everything().create_new_structure(structure)
        assert isinstance(new_stc, AbstractStructure)
        assert len(new_stc) == len(structure)
        assert new_stc[0] is not structure[0]
        assert sorted(new_stc[0].atoms.keys()) == sorted(structure[0].atoms.keys())
        assert sorted(new_stc[0].pseudoatoms.keys()) == sorted(
            structure[0].pseudoatoms.keys()
        )

    def test_specify(self, structure):
        sel = selection.Everything().specify(structure)
        assert type(sel) is selection.Set
        assert len(tuple(sel)) == len(structure)

        segment_6 = PartialStructure(tuple(structure)[:6])
        sel = selection.Everything().specify(segment_6)
        assert type(sel) is selection.Set
        assert len(tuple(sel)) == 6


class TestSetSelection:
    def test_specify(self, structure):
        sel = selection.Set([i.get_pdb_id() for i in tuple(structure)[:6]])
        assert type(sel) is selection.Set
        assert len(tuple(sel)) == 6

        new_sel = sel.specify(structure)
        assert type(new_sel) is selection.Set
        assert len(tuple(new_sel)) == 6

    def test_create_structure(self, structure):
        sel = selection.Set([i.get_pdb_id() for i in tuple(structure)[:6]])
        stc = sel.create_structure(structure)
        assert isinstance(stc, AbstractStructure)
        assert len(stc) == 6
        assert stc[0] is structure[0]
        assert sorted(stc[0].atoms.keys()) == sorted(structure[0].atoms.keys())

    def test_create_new_structure(self, structure):
        sel = selection.Set([i.get_pdb_id() for i in tuple(structure)[:6]])

        stc = sel.create_new_structure(structure)
        assert isinstance(stc, AbstractStructure)
        assert len(stc) == 6
        assert stc[0] is not structure[0]
        assert sorted(stc[0].atoms.keys()) == sorted(structure[0].atoms.keys())


class TestRangeSelection:
    @staticmethod
    def create_6_mer_range(structure):
        start = structure[0].get_pdb_id()
        end = structure[5].get_pdb_id()
        return selection.Range(start, end), start, end

    def test_range_specify(self, structure):
        range_selection, start, end = self.create_6_mer_range(structure)
        new_sel = range_selection.specify(structure)
        assert type(new_sel) is selection.Set
        assert len(tuple(new_sel)) == 6
        assert start in new_sel.ids
        assert end in new_sel.ids

    def test_create_segment(self, structure):
        chainable = True
        for mer in tuple(structure)[0:6]:
            if not isinstance(mer, MerChainable):
                chainable = False
        range_selection, start, end = self.create_6_mer_range(structure)
        if chainable:
            segment = range_selection.create_segment(structure)
            assert type(segment).__name__ == "Segment"
            assert len(segment) == 6
        else:
            with pytest.raises(DiscontinuityError):
                range_selection.create_segment(structure)
            assert True

    @pytest.mark.system
    def test_range_on_discontinuity_chain(self):
        sl = StructureLoader()
        pth = os.path.join(TEST_STRUCTURES_DIR, "prots_only", "3NPU.pdb")
        stc = sl.load_structures(path=pth)[0]
        # discontinuity occurs between A19 and A24 (1! res missing)
        msg = "Something is wrong with structure that supposed have broken " "backbone."
        assert str(stc[17].get_pdb_id()) == "A19", msg
        assert str(stc[18].get_pdb_id()) == "A24", msg

        range_selection = selection.Range(stc[16].get_pdb_id(), stc[19].get_pdb_id())
        new_sel = range_selection.specify(stc)

        assert len(tuple(new_sel)) == 4
        with pytest.raises(DiscontinuityError):
            range_selection.create_segment(stc)


class TestChainSelection:
    def test_specify(self, structure):
        for chain in structure.chains:
            chain_selection = selection.ChainSelection(chain.chain_name)
            new_selection = chain_selection.specify(structure)
            assert len(tuple(new_selection)) == len(chain)
            for mer in chain:
                assert mer.get_pdb_id() in new_selection.ids

    def test_not_existing_chain(self, structure):
        chain_selection = selection.ChainSelection("fake_name")
        new_selection = chain_selection.specify(structure)
        assert len(tuple(new_selection)) == 0

    def test_chain_from_mixed_structure(self, structure):
        if len(structure.chains) > 2:
            pytest.skip("Not enough chains to perform test.")
        chain_samples = {
            chain.chain_name: tuple(chain)[-5:] for chain in structure.chains
        }
        mixed_mers = tuple(itertools.chain(*list(chain_samples.values())))
        mixed_partial = PartialStructure(mixed_mers, structure.converter)

        test_chain_name = max(chain_samples)

        chain_selection = selection.ChainSelection(test_chain_name)

        new_sel = chain_selection.specify(mixed_partial)

        assert len(tuple(new_sel)) == len(chain_samples[test_chain_name])
        for mer in chain_samples[test_chain_name]:
            assert mer.get_pdb_id() in new_sel.ids


class TestMerNameSelection:
    def test_specify(self, structure):
        mers = tuple(structure)[:6]
        for mer in mers:
            mer_name_selection = selection.MerName(mer.name)
            new_selection = mer_name_selection.specify(structure)
            assert mer.get_pdb_id() in new_selection
            assert len(tuple(new_selection)) >= 1
            new_stc = mer_name_selection.create_structure(structure)
            for selected_mer in new_stc:
                assert selected_mer.name == mer.name


class TestMerExactTypeSelection:
    def test_wrong_class(self):
        with pytest.raises(WrongMerType):
            selection.MerExactType(type(None))

    def test_specify_residue(self, stc_2dlc):
        residue_selection = selection.MerExactType(Residue)
        residues_set = residue_selection.specify(stc_2dlc)
        assert len(tuple(residues_set)) > 0
        for mer in stc_2dlc.get_chain("X"):
            if not mer.is_chainable():
                continue
            assert mer.get_pdb_id() in residues_set.ids

    def test_specify_nucleotide(self, stc_2dlc):
        nucleotide_selection = selection.MerExactType(Nucleotide)
        nucleotides_set = nucleotide_selection.specify(stc_2dlc)
        assert len(tuple(nucleotides_set)) > 0
        for mer in stc_2dlc.get_chain("Y"):
            if not mer.is_chainable():
                continue
            assert mer.get_pdb_id() in nucleotides_set.ids

    def test_specify_ligand(self, stc_2dlc):
        ligand_selection = selection.MerExactType(Ligand)
        ligands_set = ligand_selection.specify(stc_2dlc)
        assert len(tuple(ligands_set)) > 0

    def test_specify_ion(self, stc_2dlc):
        ion_selection = selection.MerExactType(Ion)
        ions_set = ion_selection.specify(stc_2dlc)
        assert len(tuple(ions_set)) > 0


class TestNothing:
    """It is testing nothing, what did you expect?"""

    def test_specify(self, structure):
        nothing_sel = selection.Nothing()
        empty_set = nothing_sel.specify(structure)
        assert len(tuple(empty_set)) == 0


@pytest.mark.system
class TestComplexSelections:
    def test_name_and_type_intersection(self, structure):
        the_mer = max(structure, key=lambda mer: mer.is_chainable())
        test_name = the_mer.name
        test_type = type(the_mer)
        sel1 = selection.MerName(test_name)
        sel2 = selection.MerExactType(test_type)
        intersection = sel1 * sel2
        assert isinstance(intersection, selection.SelectionsIntersection)

        set_sel = intersection.specify(structure)
        assert the_mer.get_pdb_id() in set_sel.ids

    def test_chains_union(self, structure):
        if len(structure.chains) < 2:
            pytest.skip("Not enough chains to perform test.")
        chain_sels = []
        for chain in structure.chains:
            chain_sels.append(selection.ChainSelection(chain.chain_name))

        union = selection.SelectionsUnion(chain_sels)
        new_sel = union.specify(structure)

        for mer in structure:
            assert mer.get_pdb_id() in new_sel.ids

        union = chain_sels[0] + chain_sels[1]
        assert isinstance(union, selection.SelectionsUnion)

    def test_complement(self, stc_2dlc):
        res_sel = selection.MerExactType(Residue)
        prot_chain = selection.ChainSelection("X")
        diff = prot_chain - res_sel
        assert isinstance(diff, selection.SelectionsComplement)
        new_sel = diff.specify(stc_2dlc)
        assert len(new_sel.ids) == 1  # single Mg ion

        nuc_sel = selection.MerExactType(Nucleotide)
        nuc_chain = selection.ChainSelection("Y")
        diff2 = nuc_chain - nuc_sel
        new_stc = diff2.create_structure(stc_2dlc)
        assert len(new_stc) > 0
        for mer in new_stc:
            assert not isinstance(mer, Nucleotide)
