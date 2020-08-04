import itertools
import os.path

import pytest

from pydesc import selection
from pydesc.config import ConfigManager
from pydesc.chemistry.base import Mer
from pydesc.chemistry.factories import CopyingFactor
from pydesc.chemistry.factories import WrongAtomSetType
from pydesc.chemistry.full_atom import MonoatomicIon
from pydesc.chemistry.full_atom import Compound
from pydesc.chemistry.full_atom import Nucleotide
from pydesc.chemistry.full_atom import Residue
from pydesc.structure import StructureLoader
from pydesc.structure.topology import AbstractStructure
from pydesc.structure.topology import PartialStructure
from pydesc.warnexcept import DiscontinuityError

ConfigManager.warnings.quiet = True


@pytest.fixture(scope="module")
def structure(structure_file_w_type_short):
    sl = StructureLoader()
    structure = sl.load_structures(path=structure_file_w_type_short)[0]
    return structure


@pytest.fixture(scope="module")
def stc_2dlc(structures_dir):
    sl = StructureLoader()
    pth = os.path.join(structures_dir, "mixed", "2DLC.cif")
    return sl.load_structures(path=pth)[0]


class SelectionTestBase:
    @staticmethod
    def assert_atoms(new_mer, old_mer):
        new_sorted_atoms = sorted(new_mer.atoms.keys())
        old_sorted_atoms = sorted(old_mer.atoms.keys())
        assert new_sorted_atoms == old_sorted_atoms


class TestSelectorCreateNewStructure(SelectionTestBase):
    picker = selection.Selector(CopyingFactor())

    def test_everything(self, structure):
        sel = selection.Everything()
        new_structure = self.picker.create_new_structure(sel, structure)
        assert new_structure != structure
        for mer1, mer2 in zip(structure, new_structure):
            assert type(mer1) == type(mer2)
            for at1, at2 in zip(mer1, mer2):
                assert tuple(at1.vector) == tuple(at2.vector)
        assert isinstance(new_structure, AbstractStructure)
        assert len(new_structure) == len(structure)
        assert new_structure[0] is not structure[0]
        self.assert_atoms(new_structure[0], structure[0])
        new_sorted_pseudoatoms = sorted(new_structure[0].pseudoatoms.keys())
        old_sorted_pseudoatoms = sorted(structure[0].pseudoatoms.keys())
        assert new_sorted_pseudoatoms == old_sorted_pseudoatoms
        assert new_structure

    def test_set(self, structure):
        get_id = structure.converter.get_pdb_id
        ids = [get_id(mer.ind) for mer in tuple(structure)[:6]]
        sel = selection.Set(ids)
        stc = self.picker.create_new_structure(sel, structure)
        assert isinstance(stc, AbstractStructure)
        assert len(stc) == 6
        assert stc[0] is not structure[0]
        self.assert_atoms(stc[0], structure[0])

    def test_range(self, structure):
        get_id = structure.converter.get_pdb_id
        chain0 = structure.chains[0]
        start_mer, *dummy, end_mer = chain0
        start_pdb = get_id(start_mer.ind)
        end_pdb = get_id(end_mer.ind)
        sel = selection.Range(start_pdb, end_pdb)
        new_structure = self.picker.create_new_structure(sel, structure)
        assert len(new_structure) == len(chain0)
        self.assert_atoms(new_structure[0], chain0[0])

    def test_chain(self, structure):
        chain0 = structure.chains[0]
        chain_name = chain0.chain_name
        sel = selection.ChainSelection(chain_name)
        new_structure = self.picker.create_new_structure(sel, structure)
        assert len(new_structure) == len(chain0)
        self.assert_atoms(new_structure[0], chain0[0])

    def test_mer_name(self, structure):
        names = {mer.name: mer for mer in structure}
        for name, mer in names.items():
            sele = selection.AtomSetName(name)
            new_structure = self.picker.create_new_structure(sele, structure)
            self.assert_atoms(new_structure[-1], mer)

    def test_mer_type(self, structure):
        types = {type(mer): mer for mer in structure}
        for type_, mer in types.items():
            sele = selection.AtomSetExactType(type_)
            new_structure = self.picker.create_new_structure(sele, structure)
            self.assert_atoms(new_structure[-1], mer)

    def test_nothing(self, structure):
        sele = selection.Nothing()
        new_structure = self.picker.create_new_structure(sele, structure)
        assert len(new_structure) == 0


class TestEverythingSelection(SelectionTestBase):
    """Test methods from Everything selection."""

    def test_specify(self, structure):
        sel = selection.Everything().specify(structure)
        assert type(sel) is selection.Set
        assert len(tuple(sel)) == len(structure)

        segment_6 = PartialStructure(structure, tuple(structure)[:6])
        sel = selection.Everything().specify(segment_6)
        assert type(sel) is selection.Set
        assert len(tuple(sel)) == 6

    def test_create_structure(self, structure):
        sel = selection.Everything()
        new_structure = sel.create_structure(structure)
        assert new_structure == structure
        assert new_structure[0] is structure[0]
        self.assert_atoms(new_structure[0], structure[0])


class TestSetSelection(SelectionTestBase):
    """Test methods from Set type of selection."""

    def test_specify(self, structure):
        get_id = structure.converter.get_pdb_id
        sel = selection.Set([get_id(mer.ind) for mer in tuple(structure)[:6]])
        assert type(sel) is selection.Set
        assert len(tuple(sel)) == 6

        new_sel = sel.specify(structure)
        assert type(new_sel) is selection.Set
        assert len(tuple(new_sel)) == 6

    def test_create_structure(self, structure):
        get_id = structure.converter.get_pdb_id
        sel = selection.Set([get_id(mer.ind) for mer in tuple(structure)[:6]])
        new_structure = sel.create_structure(structure)
        assert isinstance(new_structure, AbstractStructure)
        assert len(new_structure) == 6
        assert new_structure[0] is structure[0]
        self.assert_atoms(new_structure[0], structure[0])


class TestRangeSelection:
    """Test methods from Range selection."""

    @staticmethod
    def create_6_mer_range(structure):
        get_id = structure.converter.get_pdb_id
        start = get_id(structure[0].ind)
        end = get_id(structure[5].ind)
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
            if not isinstance(mer, Mer):
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
    def test_range_on_discontinuity_chain(self, structures_dir):
        sl = StructureLoader()
        pth = os.path.join(structures_dir, "prots_only", "3NPU.pdb")
        stc = sl.load_structures(path=pth)[0]
        get_id = stc.converter.get_pdb_id
        # discontinuity occurs between A19 and A24 (1! res missing)
        msg = "Something is wrong with structure that supposed have broken " "backbone."
        assert str(get_id(17)) == "A19", msg
        assert str(get_id(18)) == "A24", msg

        range_selection = selection.Range(get_id(16), get_id(19))
        new_sel = range_selection.specify(stc)

        assert len(tuple(new_sel)) == 4
        with pytest.raises(DiscontinuityError):
            range_selection.create_segment(stc)


class TestChainSelection:
    def test_specify(self, structure):
        get_id = structure.converter.get_pdb_id
        for chain in structure.chains:
            chain_selection = selection.ChainSelection(chain.chain_name)
            new_selection = chain_selection.specify(structure)
            assert len(tuple(new_selection)) == len(chain)
            for mer in chain:
                assert get_id(mer.ind) in new_selection.ids

    def test_not_existing_chain(self, structure):
        chain_selection = selection.ChainSelection("fake_name")
        new_selection = chain_selection.specify(structure)
        assert len(tuple(new_selection)) == 0

    def test_chain_from_mixed_structure(self, structure):
        get_id = structure.converter.get_pdb_id
        if len(structure.chains) > 2:
            pytest.skip("Not enough chains to perform test.")
        chain_samples = {
            chain.chain_name: tuple(chain)[-5:] for chain in structure.chains
        }
        mixed_mers = tuple(itertools.chain(*list(chain_samples.values())))
        mixed_partial = PartialStructure(structure, mixed_mers)

        test_chain_name = max(chain_samples)

        chain_selection = selection.ChainSelection(test_chain_name)

        new_sel = chain_selection.specify(mixed_partial)

        assert len(tuple(new_sel)) == len(chain_samples[test_chain_name])
        for mer in chain_samples[test_chain_name]:
            assert get_id(mer.ind) in new_sel.ids


class TestAtomSetNameSelection:
    def test_specify(self, structure):
        get_id = structure.converter.get_pdb_id
        mers = tuple(structure)[:6]
        for mer in mers:
            mer_name_selection = selection.AtomSetName(mer.name)
            new_selection = mer_name_selection.specify(structure)
            assert get_id(mer.ind) in new_selection
            assert len(tuple(new_selection)) >= 1
            new_stc = mer_name_selection.create_structure(structure)
            for selected_mer in new_stc:
                assert selected_mer.name == mer.name


class TestAtomSetExactTypeSelection:
    def test_wrong_class(self):
        with pytest.raises(WrongAtomSetType):
            selection.AtomSetExactType(type(None))

    def test_specify_residue(self, stc_2dlc):
        get_id = stc_2dlc.converter.get_pdb_id
        residue_selection = selection.AtomSetExactType(Residue)
        residues_set = residue_selection.specify(stc_2dlc)
        assert len(tuple(residues_set)) > 0
        for mer in stc_2dlc.get_chain("X"):
            if not mer.is_chainable():
                continue
            assert get_id(mer.ind) in residues_set.ids

    def test_specify_nucleotide(self, stc_2dlc):
        get_id = stc_2dlc.converter.get_pdb_id
        nucleotide_selection = selection.AtomSetExactType(Nucleotide)
        nucleotides_set = nucleotide_selection.specify(stc_2dlc)
        assert len(tuple(nucleotides_set)) > 0
        for mer in stc_2dlc.get_chain("Y"):
            if not mer.is_chainable():
                continue
            assert get_id(mer.ind) in nucleotides_set.ids

    def test_specify_ligand(self, stc_2dlc):
        ligand_selection = selection.AtomSetExactType(Compound)
        ligands_set = ligand_selection.specify(stc_2dlc)
        assert len(tuple(ligands_set)) > 0

    def test_specify_ion(self, stc_2dlc):
        ion_selection = selection.AtomSetExactType(MonoatomicIon)
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
    """Test Union, Complement and Intersection selections."""

    def test_name_and_type_intersection(self, structure):
        get_id = structure.converter.get_pdb_id
        the_mer = max(structure, key=lambda mer: mer.is_chainable())
        test_name = the_mer.name
        test_type = type(the_mer)
        sel1 = selection.AtomSetName(test_name)
        sel2 = selection.AtomSetExactType(test_type)
        intersection = sel1 * sel2
        assert isinstance(intersection, selection.SelectionsIntersection)

        set_sel = intersection.specify(structure)
        assert get_id(the_mer.ind) in set_sel.ids

    def test_chains_union(self, structure):
        get_id = structure.converter.get_pdb_id
        if len(structure.chains) < 2:
            pytest.skip("Not enough chains to perform test.")
        chain_sels = []
        for chain in structure.chains:
            chain_sels.append(selection.ChainSelection(chain.chain_name))

        union = selection.SelectionsUnion(chain_sels)
        new_sel = union.specify(structure)

        for mer in structure:
            assert get_id(mer.ind) in new_sel.ids

        union = chain_sels[0] + chain_sels[1]
        assert isinstance(union, selection.SelectionsUnion)

    def test_complement(self, stc_2dlc):
        res_sel = selection.AtomSetExactType(Residue)
        prot_chain = selection.ChainSelection("X")
        diff = prot_chain - res_sel
        assert isinstance(diff, selection.SelectionsComplement)
        new_sel = diff.specify(stc_2dlc)
        assert len(new_sel.ids) == 1  # single Mg ion

        nuc_sel = selection.AtomSetExactType(Nucleotide)
        nuc_chain = selection.ChainSelection("Y")
        diff2 = nuc_chain - nuc_sel
        new_stc = diff2.create_structure(stc_2dlc)
        assert len(new_stc) > 0
        for mer in new_stc:
            assert not isinstance(mer, Nucleotide)
