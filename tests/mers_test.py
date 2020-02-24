import os.path

import Bio.PDB
import pytest

from pydesc.mers.factories import BioPythonMerFactory
from pydesc.mers.full_atom import Ion
from pydesc.mers.full_atom import Nucleotide
from pydesc.mers.full_atom import Residue
from tests.conftest import PURE_TYPES_2_MERS_DICT
from tests.conftest import TEST_STRUCTURES_DIR

TYPE_THRESHOLDS = {Nucleotide: 0.25, Residue: 0.01, Ion: 0.0}


class TestMonomerFactory(object):
    def test_default_mer_factory_create_from_pdb_res(self, structure_file_w_pure_type):
        type_, fname = os.path.split(structure_file_w_pure_type)
        factory = BioPythonMerFactory()
        pth = os.path.join(TEST_STRUCTURES_DIR, structure_file_w_pure_type)
        pdb_structure = Bio.PDB.PDBParser(QUIET=True).get_structure(fname, pth)

        expected_type = PURE_TYPES_2_MERS_DICT[type_]
        type_threshold = TYPE_THRESHOLDS[expected_type]

        points = {True: 0, False: 1}
        length = 0

        for model in pdb_structure:
            for chain in model:
                for residue in chain:
                    result, warns = factory.create(residue, warn_in_place=False)

                    if result is None:
                        assert warns is None
                        assert residue.get_id()[0] == "W"
                    else:
                        length += 1
                        points[expected_type in result] += 1

        assert (length - points[True]) / float(
            length
        ) <= type_threshold, "%i out of %i mers are of expected type" % (
            points[True],
            length,
        )

    def test_create_residue_from_pdb_res(self):
        factory = BioPythonMerFactory()
        pdb_structure = Bio.PDB.PDBParser(QUIET=True).get_structure(
            "5MPV.pdb", os.path.join(TEST_STRUCTURES_DIR, "prots_only", "5MPV.pdb")
        )

        r1 = pdb_structure[0]["D"][15]
        r2 = pdb_structure[0]["D"][16]
        r3 = pdb_structure[0]["D"][17]
        prev = factory.create(r1, warn_in_place=False)[0][Residue]
        res = factory.create(r2, warn_in_place=False)[0][Residue]
        next = factory.create(r3, warn_in_place=False)[0][Residue]

        res.previous_mer = prev
        res.next_mer = next

        assert len(tuple(res.iter_bb_atoms())) == 3
        assert len(tuple(res.iter_nbb_atoms())) > 0

        assert res.angles
        assert res.cbx
        assert res.rc
        assert res.backbone_average
        assert res.ca

        assert [round(i, 2) for i in res.angles] == [0.62, 1.30]
        assert res.next_mer is next
        assert res.previous_mer is prev


class MerTest(object):
    @staticmethod
    def iter_structure(file_, class_):
        factory = BioPythonMerFactory()
        pdb_structure = Bio.PDB.PDBParser(QUIET=True).get_structure(
            "test", os.path.join(TEST_STRUCTURES_DIR, file_)
        )
        to_return = []
        for model in pdb_structure:
            for chain in model:
                for residue in chain:
                    result, warns = factory.create(residue, warn_in_place=False)
                    if result is None or class_ not in result:
                        continue
                    result = result[class_]
                    to_return.append(result)
        return to_return


class TestResidue(MerTest):
    def test_calculate_cbx(self, protein_file):
        checked = 0
        for result in self.iter_structure(protein_file, Residue):
            if result.name == "GLY":
                continue
            assert result.cbx
            del result.pseudoatoms["cbx"]
            assert result.cbx
            length = (result.CB - result.cbx).calculate_length()
            assert round(length, 2) == 1.0
            unit_a_bx = (result.cbx - result.CA).get_unit_vector()
            unit_a_b = (result.CB - result.CA).get_unit_vector()
            dot_prod = unit_a_b.dot(unit_a_bx)
            assert pytest.approx(dot_prod) == 1.0
            checked += 1
        assert checked != 0


class TestNucleotide(MerTest):
    def test_calculate_features(self, nuclei_file):
        checked = 0
        for result in self.iter_structure(nuclei_file, Nucleotide):
            assert result.nx
            assert result.ring_center
            assert result.prc
            for pseudoatom in result.nx, result.ring_center, result.prc:
                dist_from_rc = (result.rc - pseudoatom).calculate_length()
                assert dist_from_rc < 10
            checked += 1
        assert checked != 0
