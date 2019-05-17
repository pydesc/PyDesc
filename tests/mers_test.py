import pytest
import os.path
import Bio.PDB

from pydesc.mers import MonomerFactory
from pydesc.mers import Nucleotide
from pydesc.mers import Residue
from pydesc.mers import Ion

from tests.conftest import PDB_FILES_WITH_TYPE
from tests.conftest import TEST_STRUCTURES_DIR
from tests.conftest import TYPE_DICT
from tests.conftest import DIR_DICT

TYPE_THRESHOLDS = {Nucleotide: .25,
                   Residue: .01,
                   Ion: .0,
                   }


class TestMonomerFactory(object):

    @pytest.mark.parametrize('type_, struc_file', PDB_FILES_WITH_TYPE)
    def test_default_monomer_factory_create_from_pdb_res(self, type_, struc_file):

        factory = MonomerFactory()
        pdb_structure = Bio.PDB.PDBParser(QUIET=True).get_structure(
            struc_file,
            os.path.join(
                TEST_STRUCTURES_DIR,
                type_,
                struc_file)
        )

        expected_type = TYPE_DICT[type_]
        type_threshold = TYPE_THRESHOLDS[expected_type]

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

        assert (length - points[True]) / float(length) <= type_threshold, \
            '%i out of %i residues are of expected type' % (points[True], length)

    def test_create_residue_from_pdb_res(self):
        factory = MonomerFactory()
        pdb_structure = Bio.PDB.PDBParser(QUIET=True).get_structure(
            '5MPV.pdb',
            os.path.join(
                TEST_STRUCTURES_DIR,
                'prots_only',
                '5MPV.pdb')
        )

        r1 = pdb_structure[0]['D'][15]
        r2 = pdb_structure[0]['D'][16]
        r3 = pdb_structure[0]['D'][17]
        prev = factory.create_from_biopdb(r1, warn_in_place=False)[0][Residue]
        res = factory.create_from_biopdb(r2, warn_in_place=False)[0][Residue]
        next = factory.create_from_biopdb(r3, warn_in_place=False)[0][Residue]

        res.previous_mer = prev
        res.next_mer = next
        res.finalize()

        assert len(tuple(res.iter_bb_atoms())) == 3
        assert len(tuple(res.iter_nbb_atoms())) > 0
        assert 'cbx' in res.pseudoatoms
        assert 'rc' in res.pseudoatoms
        assert 'backbone_average' in res.pseudoatoms
        assert 'CA' in res.atoms

        res.calculate_angles()
        res.calculate_cbx()
        res.calculate_rc()

        assert [round(i, 2) for i in res.angles] == [0.62, 1.30]
        assert res.next_mer is next
        assert res.previous_mer is prev


class MerTest(object):

    @staticmethod
    def iter_structure(dir_, file_, class_):
        factory = MonomerFactory()
        pdb_structure = Bio.PDB.PDBParser(QUIET=True).get_structure(
            file_,
            os.path.join(
                TEST_STRUCTURES_DIR,
                dir_,
                file_)
        )
        to_return = []
        for model in pdb_structure:
            for chain in model:
                for residue in chain:
                    result, warns = factory.create_from_biopdb(residue, warn_in_place=False)
                    if result is None or class_ not in result:
                        continue
                    result = result[class_]
                    to_return.append(result)
        return to_return


class TestResidue(MerTest):

    @pytest.mark.parametrize('type_, structure_file', [i for i in PDB_FILES_WITH_TYPE if i[0] == DIR_DICT[Residue]])
    def test_calculate_cbx(self, type_, structure_file):
        checked = 0
        for result in self.iter_structure(type_, structure_file, Residue):
            if result.name == 'GLY':
                continue
            assert 'cbx' in result.pseudoatoms
            del result.pseudoatoms['cbx']
            Residue.calculate_cbx(result)
            assert round((result.CB - result.cbx).calculate_length(), 2) == 1.
            unit_a_bx = (result.cbx - result.CA).get_unit_vector()
            unit_a_b = (result.CB - result.CA).get_unit_vector()
            dot_prod = unit_a_b.dot(unit_a_bx)
            assert pytest.approx(dot_prod) == 1.
            checked += 1
        assert checked != 0


class TestNucleotide(MerTest):

    @pytest.mark.parametrize('type_, structure_file', [i for i in PDB_FILES_WITH_TYPE if i[0] == DIR_DICT[Nucleotide]])
    def test_calculate_features(self, type_, structure_file):
        checked = 0
        for result in self.iter_structure(type_, structure_file, Nucleotide):
            result.finalize()
            for feat, mth in (
                    ('nx', Nucleotide.calculate_nx),
                    ('prc', Nucleotide.calculate_proximate_ring_center),
                    ('ring_center', Nucleotide.calculate_ring_center),
            ):
                assert feat in result.pseudoatoms
                del result.pseudoatoms[feat]
                mth(result)
                if result.pseudoatoms[feat] is None:  # like 'prc'
                    continue
                assert (result.rc - getattr(result, feat)).calculate_length() < 10
            checked += 1
        assert checked != 0
